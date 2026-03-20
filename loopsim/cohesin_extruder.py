"""
Cohesin Loop Extrusion Engine — Stochastic Simulation of DSB-Anchored Extrusion
=================================================================================

This is the computational heart of LoopSim. It models the **loop extrusion**
process that occurs after a DNA double-strand break, as described in the
landmark study:

    Marin-Gonzalez A, et al. "Cohesin drives chromatin scanning during the
    RAD51-mediated homology search." Science, 2025.

What Happens Biologically After a DSB
--------------------------------------
When Cas9 (or any nuclease) creates a DSB, a carefully orchestrated sequence
of events unfolds at the break site:

1. **DNA damage sensing (seconds).**
   The MRN complex (MRE11-RAD50-NBS1) recognizes the broken DNA ends within
   seconds. MRN recruits and activates the ATM kinase.

2. **gamma-H2AX spreading (minutes).**
   ATM phosphorylates the histone variant H2AX at serine 139 (creating
   gamma-H2AX) across a domain of ~1-2 megabases flanking the break. This
   spreading is mediated by a positive feedback loop: ATM -> gamma-H2AX ->
   MDC1 -> more ATM. The gamma-H2AX domain is visible as a "focus" by
   immunofluorescence microscopy.

3. **Cohesin loading at the DSB (minutes).**
   The gamma-H2AX/MDC1 platform recruits NIPBL (the cohesin loader), which
   deposits the **cohesin complex** at the break site. Cohesin is a large
   ring-shaped protein complex consisting of:
   - **SMC1** and **SMC3**: two ATPase subunits that form a V-shape
   - **RAD21** (also called SCC1/Mcd1): the kleisin subunit that bridges
     the SMC heads, closing the ring
   - **SA1 or SA2** (STAG1/STAG2): regulatory subunits

   The cohesin ring topologically encircles DNA — it physically embraces
   the chromatin fiber inside its lumen (~40 nm inner diameter).

   Source: Arnould C, et al. "Loop extrusion as a mechanism for formation
   of DNA damage repair foci." Nature, 2021.

4. **Loop extrusion begins.**
   Once loaded, the cohesin complex acts as a **molecular motor** that
   performs bidirectional loop extrusion:

   - The cohesin ring stays anchored at the DSB (the "anchor point").
   - Chromatin from BOTH sides of the DSB is simultaneously reeled through
     the ring, creating an expanding chromatin loop.
   - The motor is powered by ATP hydrolysis at the SMC1/SMC3 heads.

   Think of it like this: imagine holding a loop of rope and pulling more
   rope through from both ends simultaneously. The point where you hold
   the rope = the DSB anchor. The rope being pulled through = chromatin.

   Single-molecule experiments (Davidson IF, et al., Science, 2019; Kim Y,
   et al., Science, 2019) measured the extrusion rate at ~0.5-2 kb/sec
   in vitro. We use 1 kb/sec as the default.

5. **CTCF stalling.**
   As the extrusion arms advance along the fiber, they may encounter bound
   CTCF proteins. CTCF acts as a **directional roadblock**:

   - CTCF's N-terminus interacts with cohesin's SA2-RAD21 interface.
   - This interaction ONLY occurs when cohesin approaches from the correct
     direction (N-terminal side of CTCF).
   - In our convention:
     * Forward CTCF (+ strand) stalls the RIGHT-ward extrusion arm
     * Reverse CTCF (- strand) stalls the LEFT-ward extrusion arm

   Stalling is **probabilistic**: in a population of cells, a CTCF site
   stalls cohesin in a fraction of cells proportional to its occupancy/
   strength. Weak sites allow bypass more often.

   Source: Li Y, et al. "The structural basis for cohesin-CTCF-anchored
   loops." Nature, 2020.

6. **Search domain defined.**
   The genomic region that is reeled through the DSB anchor during
   extrusion constitutes the **homology search domain**. The RAD51
   filament (formed on the resected 3' ssDNA overhang at the DSB) can
   sample the passing chromatin for sequence complementarity. This is how
   the cell finds homologous sequences for HDR.

   The search domain is bounded by:
   - CTCF stall points (the "walls" of the loop)
   - The time available for extrusion before repair completes or the cell
     progresses through the cell cycle
   - The physical boundaries of the chromatin fiber (centromere, telomere)

Simulation Algorithm (per run)
------------------------------
Each call to ``extrude()`` performs one stochastic simulation:

    1. Place DSB at genomic position X -> convert to bead index.
    2. Initialize LEFT arm at bead X, RIGHT arm at bead X.
    3. For each time step (dt = 1 second):
       a. LEFT ARM (extending leftward, toward lower genomic coords):
          - Compute advance probability:
            * Open chromatin: advance with prob = 1.0
            * Closed chromatin: advance with prob = heterochromatin_slowdown
          - If advancing: move left arm 1 bead to the left.
          - Check: does the new bead have a REVERSE-oriented CTCF site?
            * If yes: draw random number. If < CTCF_strength: STALL
              (this arm stops permanently for the rest of the simulation).
            * If no, or if bypass: continue.
          - Check: has the arm reached the left edge of the fiber?
            * If yes: stop this arm.
       b. RIGHT ARM (extending rightward, same logic but mirrored):
          - Advance rightward.
          - Check for FORWARD-oriented CTCF sites.
       c. Record current domain: [left_arm_position, right_arm_position].
       d. TERMINATION: stop when both arms are stalled/stopped, OR
          when max_time_sec is reached.
    4. Return ExtrusionResult with final boundaries and metadata.

Monte Carlo Ensemble
--------------------
The ``simulate()`` method repeats the above N times (default 1000) with
different random seeds. Because CTCF stalling is probabilistic and
heterochromatin slowdown is stochastic, each run gives a slightly different
domain. The ensemble captures the **cell-to-cell variability** in the
search domain, which is biologically real — different cells in a population
will extrude to different extents.

Performance
-----------
Each simulation is O(max_time) per run, and the time steps are simple
integer operations with a few random draws. For 1000 runs of 3600 time
steps each, this completes in ~1-2 seconds on a modern laptop.

For the full Monte Carlo ensemble, we also provide a vectorized path
(``_simulate_vectorized``) that runs all simulations simultaneously using
NumPy broadcasting, achieving ~10-50x speedup.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from .chromatin_fiber import ChromatinFiber


# ---------------------------------------------------------------------------
# Result data classes
# ---------------------------------------------------------------------------

@dataclass
class ExtrusionResult:
    """Result of a single cohesin loop extrusion simulation.

    This dataclass stores the outcome of one stochastic extrusion run:
    where the left and right arms ended up, whether they stalled at CTCF,
    and how long the extrusion took.

    The **search domain** is defined as the genomic interval
    [left_boundary_bp, right_boundary_bp]. Any endogenous sequence
    (e.g., sister chromatid homology) falling within this interval had the
    opportunity to be scanned by the RAD51 filament during extrusion.

    Attributes
    ----------
    dsb_position_bp : int
        Genomic coordinate of the DSB (the cohesin loading / anchor site).
    left_boundary_bp : int
        Leftmost genomic coordinate reached by the left extrusion arm.
    right_boundary_bp : int
        Rightmost genomic coordinate reached by the right extrusion arm.
    domain_size_bp : int
        Total size of the search domain (right - left boundary).
    left_stalled_at_ctcf : bool
        True if the left arm was stopped by a CTCF site (rather than by
        reaching the fiber edge or running out of time).
    right_stalled_at_ctcf : bool
        True if the right arm was stopped by a CTCF site.
    extrusion_time_sec : float
        Total simulated time in seconds until both arms stopped (or
        max_time was reached).
    n_ctcf_bypassed : int
        Number of CTCF sites that cohesin encountered but successfully
        bypassed (the random draw exceeded the site's stall strength).
        Higher values indicate a more permeable chromatin landscape.
    time_trace_left : np.ndarray or None
        If recorded, the per-time-step position of the left arm (bead index).
        Used for kymograph visualization. None if not recorded.
    time_trace_right : np.ndarray or None
        If recorded, the per-time-step position of the right arm.
    """

    dsb_position_bp: int
    left_boundary_bp: int
    right_boundary_bp: int
    domain_size_bp: int
    left_stalled_at_ctcf: bool
    right_stalled_at_ctcf: bool
    extrusion_time_sec: float
    n_ctcf_bypassed: int
    time_trace_left: Optional[np.ndarray] = field(default=None, repr=False)
    time_trace_right: Optional[np.ndarray] = field(default=None, repr=False)


@dataclass
class SimulationResults:
    """Aggregated results from a Monte Carlo ensemble of extrusion simulations.

    Contains the full list of individual ``ExtrusionResult`` objects plus
    pre-computed summary statistics. This is the primary data structure
    returned by ``CohesinExtruder.simulate()``.

    Attributes
    ----------
    individual_results : list of ExtrusionResult
        All individual simulation results.
    n_simulations : int
        Number of simulations in the ensemble.
    domain_sizes_bp : np.ndarray
        Array of domain sizes from each simulation (for histograms).
    left_boundaries_bp : np.ndarray
        Array of left boundary positions (for probability calculations).
    right_boundaries_bp : np.ndarray
        Array of right boundary positions.
    mean_domain_size_bp : float
        Mean search domain size across all simulations.
    median_domain_size_bp : float
        Median search domain size.
    std_domain_size_bp : float
        Standard deviation of domain sizes.
    fraction_left_stalled_ctcf : float
        Fraction of simulations where the left arm stalled at CTCF.
    fraction_right_stalled_ctcf : float
        Fraction where the right arm stalled at CTCF.
    mean_ctcf_bypassed : float
        Average number of CTCF sites bypassed per simulation.
    """

    individual_results: List[ExtrusionResult]
    n_simulations: int
    domain_sizes_bp: np.ndarray
    left_boundaries_bp: np.ndarray
    right_boundaries_bp: np.ndarray
    mean_domain_size_bp: float
    median_domain_size_bp: float
    std_domain_size_bp: float
    fraction_left_stalled_ctcf: float
    fraction_right_stalled_ctcf: float
    mean_ctcf_bypassed: float


# ---------------------------------------------------------------------------
# CohesinExtruder
# ---------------------------------------------------------------------------

class CohesinExtruder:
    """Stochastic simulator for cohesin-mediated loop extrusion at a DSB.

    This class implements the core simulation engine. Given a ChromatinFiber
    (the "substrate") and a DSB position (where cohesin is loaded), it
    simulates bidirectional loop extrusion with stochastic CTCF stalling
    and heterochromatin slowdown.

    The simulation captures the key finding of Marin-Gonzalez et al. (2025):
    cohesin loaded at a DSB creates an expanding search domain whose
    boundaries are determined by the nearest CTCF barriers. This domain
    is the region of chromatin scanned by RAD51 during the homology search.

    Parameters
    ----------
    extrusion_rate_kb_per_sec : float, optional
        Speed of cohesin loop extrusion in kilobases per second. Default 1.0.

        Source: Single-molecule measurements by Davidson et al. (Science,
        2019) and Kim et al. (Science, 2019) measured rates of 0.5-2 kb/s
        in vitro. In vivo rates may be somewhat lower due to the crowded
        nuclear environment and nucleosome obstacles. 1 kb/s is a consensus
        estimate.

    ctcf_stall_probability : float, optional
        Global multiplier applied to CTCF site strength when determining
        stall probability. Default 0.8. The effective stall probability at
        a given site is: min(1.0, site.strength * ctcf_stall_probability).

        Set to 1.0 to use each site's intrinsic strength as-is.
        Set to 0.0 to ignore all CTCF sites (no stalling).

    max_time_sec : int, optional
        Maximum simulation time in seconds. Default 3600 (1 hour).
        After this time, extrusion is terminated regardless of stalling.

        In the cell, extrusion is terminated by WAPL-mediated cohesin
        release (residence time ~20-30 minutes for most cohesin) or by
        cell cycle progression. 1 hour is a generous upper bound; most
        biologically relevant extrusion completes within 30 minutes.

        At 1 kb/sec for 3600 sec, maximum possible domain = 7.2 Mb
        (3.6 Mb per arm). This is larger than most TADs (0.5-2 Mb),
        so CTCF stalling usually terminates extrusion well before the
        time limit.

    heterochromatin_slowdown : float, optional
        Probability multiplier for advancing through closed chromatin.
        Default 0.5 (cohesin advances through heterochromatin at half the
        normal rate).

        When cohesin encounters a bead marked as "closed" (heterochromatin),
        it advances with probability = heterochromatin_slowdown per time step,
        rather than advancing deterministically. This models the physical
        obstruction caused by densely packed H3K9me3-marked nucleosomes and
        HP1-compacted chromatin.

    record_traces : bool, optional
        Whether to record per-time-step arm positions (for kymograph plots).
        Default False (saves memory for large ensembles). Set to True when
        you want to visualize individual extrusion trajectories.

    Examples
    --------
    >>> fiber = ChromatinFiber(0, 5_000_000)
    >>> fiber.add_ctcf_sites_uniform()
    >>> extruder = CohesinExtruder(extrusion_rate_kb_per_sec=1.0)
    >>> result = extruder.extrude(fiber, dsb_position_bp=2_500_000)
    >>> print(f"Search domain: {result.domain_size_bp / 1e6:.2f} Mb")
    """

    def __init__(
        self,
        extrusion_rate_kb_per_sec: float = 1.0,
        ctcf_stall_probability: float = 0.8,
        max_time_sec: int = 3600,
        heterochromatin_slowdown: float = 0.5,
        record_traces: bool = False,
    ) -> None:
        if extrusion_rate_kb_per_sec <= 0:
            raise ValueError("Extrusion rate must be positive.")
        if not 0.0 <= ctcf_stall_probability <= 1.0:
            raise ValueError("CTCF stall probability must be in [0, 1].")
        if max_time_sec <= 0:
            raise ValueError("max_time_sec must be positive.")
        if not 0.0 <= heterochromatin_slowdown <= 1.0:
            raise ValueError("heterochromatin_slowdown must be in [0, 1].")

        self.extrusion_rate_kb_per_sec: float = extrusion_rate_kb_per_sec
        self.ctcf_stall_probability: float = ctcf_stall_probability
        self.max_time_sec: int = max_time_sec
        self.heterochromatin_slowdown: float = heterochromatin_slowdown
        self.record_traces: bool = record_traces

    # -----------------------------------------------------------------
    # Single extrusion simulation
    # -----------------------------------------------------------------

    def extrude(
        self,
        fiber: ChromatinFiber,
        dsb_position_bp: int,
        rng: Optional[np.random.Generator] = None,
    ) -> ExtrusionResult:
        """Run a single stochastic extrusion simulation.

        Cohesin is loaded at the DSB position and begins extruding
        bidirectionally along the chromatin fiber. Each time step (1 second),
        each arm attempts to advance by one bead. Stochastic events:
        - CTCF stalling (probabilistic, based on site strength)
        - Heterochromatin slowdown (probabilistic, slower advance)

        Parameters
        ----------
        fiber : ChromatinFiber
            The chromatin substrate on which extrusion occurs.
        dsb_position_bp : int
            Genomic coordinate of the DSB (cohesin loading site).
        rng : numpy.random.Generator, optional
            Random number generator. If None, creates a new default RNG.

        Returns
        -------
        ExtrusionResult
            Result of this single simulation, including domain boundaries,
            stall status, and extrusion time.
        """
        if rng is None:
            rng = np.random.default_rng()

        # Convert DSB to bead index
        dsb_bead = fiber.get_bead_index(dsb_position_bp)

        # Compute beads per time step from extrusion rate
        # If resolution = 1000 bp/bead and rate = 1 kb/sec, then 1 bead/sec
        beads_per_sec = self.extrusion_rate_kb_per_sec * 1000.0 / fiber.resolution_bp

        # Pre-fetch fiber arrays for fast access (avoid property copies in loop)
        ctcf_fwd = fiber._ctcf_forward
        ctcf_rev = fiber._ctcf_reverse
        chrom_state = fiber._chromatin_state
        n_beads = fiber.n_beads

        # Arm positions (bead indices)
        left_pos = dsb_bead
        right_pos = dsb_bead

        # Stall flags
        left_stalled = False
        right_stalled = False
        left_stalled_ctcf = False
        right_stalled_ctcf = False

        # Counters
        n_ctcf_bypassed = 0
        extrusion_time = 0.0

        # Optional time traces for kymograph
        if self.record_traces:
            trace_left = np.zeros(self.max_time_sec + 1, dtype=np.int32)
            trace_right = np.zeros(self.max_time_sec + 1, dtype=np.int32)
            trace_left[0] = left_pos
            trace_right[0] = right_pos
        else:
            trace_left = None
            trace_right = None

        # Extrusion accumulator for sub-bead steps
        # When rate != 1 bead/sec, we accumulate fractional steps
        left_accumulator = 0.0
        right_accumulator = 0.0

        for t in range(1, self.max_time_sec + 1):
            if left_stalled and right_stalled:
                extrusion_time = float(t - 1)
                # Fill remaining trace if recording
                if self.record_traces:
                    trace_left[t:] = left_pos
                    trace_right[t:] = right_pos
                break

            # ------- LEFT ARM -------
            if not left_stalled:
                # Determine advance probability
                # Check chromatin state at the next bead (one to the left)
                next_left = left_pos - 1
                if next_left < 0:
                    # Reached left edge of fiber
                    left_stalled = True
                else:
                    advance_prob = 1.0
                    if chrom_state[next_left] == 1:  # heterochromatin
                        advance_prob = self.heterochromatin_slowdown

                    # Stochastic advance
                    if rng.random() < advance_prob:
                        left_accumulator += beads_per_sec
                    # Advance by integer beads accumulated
                    while left_accumulator >= 1.0 and not left_stalled:
                        left_accumulator -= 1.0
                        new_left = left_pos - 1
                        if new_left < 0:
                            left_stalled = True
                            break

                        left_pos = new_left

                        # Check for reverse CTCF at new position
                        # (reverse CTCF blocks leftward extrusion)
                        if ctcf_rev[left_pos] > 0:
                            effective_stall = min(
                                1.0,
                                ctcf_rev[left_pos] * self.ctcf_stall_probability,
                            )
                            if rng.random() < effective_stall:
                                left_stalled = True
                                left_stalled_ctcf = True
                                break
                            else:
                                n_ctcf_bypassed += 1

            # ------- RIGHT ARM -------
            if not right_stalled:
                next_right = right_pos + 1
                if next_right >= n_beads:
                    right_stalled = True
                else:
                    advance_prob = 1.0
                    if chrom_state[next_right] == 1:
                        advance_prob = self.heterochromatin_slowdown

                    if rng.random() < advance_prob:
                        right_accumulator += beads_per_sec
                    while right_accumulator >= 1.0 and not right_stalled:
                        right_accumulator -= 1.0
                        new_right = right_pos + 1
                        if new_right >= n_beads:
                            right_stalled = True
                            break

                        right_pos = new_right

                        # Check for forward CTCF at new position
                        # (forward CTCF blocks rightward extrusion)
                        if ctcf_fwd[right_pos] > 0:
                            effective_stall = min(
                                1.0,
                                ctcf_fwd[right_pos] * self.ctcf_stall_probability,
                            )
                            if rng.random() < effective_stall:
                                right_stalled = True
                                right_stalled_ctcf = True
                                break
                            else:
                                n_ctcf_bypassed += 1

            # Record traces
            if self.record_traces:
                trace_left[t] = left_pos
                trace_right[t] = right_pos

            extrusion_time = float(t)

        else:
            # max_time reached without both arms stalling
            extrusion_time = float(self.max_time_sec)

        # Convert bead indices back to genomic coordinates
        left_bp = fiber.get_bead_position_bp(left_pos)
        right_bp = fiber.get_bead_position_bp(right_pos)
        domain_size = right_bp - left_bp

        # Trim traces to actual duration
        actual_steps = int(extrusion_time) + 1
        if self.record_traces and trace_left is not None:
            trace_left = trace_left[:actual_steps]
            trace_right = trace_right[:actual_steps]

        return ExtrusionResult(
            dsb_position_bp=dsb_position_bp,
            left_boundary_bp=left_bp,
            right_boundary_bp=right_bp,
            domain_size_bp=domain_size,
            left_stalled_at_ctcf=left_stalled_ctcf,
            right_stalled_at_ctcf=right_stalled_ctcf,
            extrusion_time_sec=extrusion_time,
            n_ctcf_bypassed=n_ctcf_bypassed,
            time_trace_left=trace_left,
            time_trace_right=trace_right,
        )

    # -----------------------------------------------------------------
    # Monte Carlo ensemble
    # -----------------------------------------------------------------

    def simulate(
        self,
        fiber: ChromatinFiber,
        dsb_position_bp: int,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> SimulationResults:
        """Run a Monte Carlo ensemble of extrusion simulations.

        Each simulation uses an independent random seed, capturing the
        stochastic variability in CTCF occupancy and heterochromatin
        effects that exists across a **population of cells**.

        In a real experiment (e.g., CRISPR editing of 10^6 cells), each cell
        will have slightly different CTCF binding patterns (CTCF binding is
        ~60-80% occupied at strong sites — meaning 20-40% of cells lack CTCF
        at any given site). This Monte Carlo ensemble approximates that
        cell-to-cell heterogeneity.

        Parameters
        ----------
        fiber : ChromatinFiber
            The chromatin substrate.
        dsb_position_bp : int
            Genomic coordinate of the DSB.
        n_simulations : int, optional
            Number of independent simulations to run. Default 1000.
            For quick exploratory runs, 100-500 is sufficient.
            For publication-quality statistics, use 5000-10000.
        seed : int, optional
            Base random seed for reproducibility. Each simulation uses
            seed + i as its individual seed. If None, non-reproducible.

        Returns
        -------
        SimulationResults
            Aggregated results including distributions and summary statistics.

        Examples
        --------
        >>> fiber = ChromatinFiber(0, 5_000_000)
        >>> fiber.add_ctcf_sites_uniform()
        >>> extruder = CohesinExtruder()
        >>> results = extruder.simulate(fiber, 2_500_000, n_simulations=1000)
        >>> print(f"Mean domain: {results.mean_domain_size_bp / 1e6:.2f} Mb")
        >>> print(f"Left stalled at CTCF: {results.fraction_left_stalled_ctcf:.0%}")
        """
        results: List[ExtrusionResult] = []

        for i in range(n_simulations):
            if seed is not None:
                rng = np.random.default_rng(seed + i)
            else:
                rng = np.random.default_rng()

            result = self.extrude(fiber, dsb_position_bp, rng=rng)
            results.append(result)

        return self._aggregate_results(results)

    # -----------------------------------------------------------------
    # Vectorized ensemble (for speed)
    # -----------------------------------------------------------------

    def simulate_vectorized(
        self,
        fiber: ChromatinFiber,
        dsb_position_bp: int,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> SimulationResults:
        """Run the Monte Carlo ensemble using vectorized NumPy operations.

        This is a performance-optimized version of ``simulate()`` that runs
        all N simulations simultaneously by operating on arrays of shape
        (n_simulations,) at each time step. This achieves ~10-50x speedup
        over the sequential loop version for large ensembles.

        The trade-off is that per-step time traces are NOT recorded (the
        ``record_traces`` setting is ignored). Use the sequential
        ``simulate()`` with ``record_traces=True`` if you need kymograph data.

        Parameters
        ----------
        fiber : ChromatinFiber
            The chromatin substrate.
        dsb_position_bp : int
            Genomic coordinate of the DSB.
        n_simulations : int, optional
            Number of simulations. Default 1000.
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        SimulationResults
            Aggregated results (same structure as ``simulate()``).
        """
        rng = np.random.default_rng(seed)

        dsb_bead = fiber.get_bead_index(dsb_position_bp)
        beads_per_sec = self.extrusion_rate_kb_per_sec * 1000.0 / fiber.resolution_bp
        n_beads = fiber.n_beads

        # Pre-fetch fiber arrays
        ctcf_fwd = fiber._ctcf_forward
        ctcf_rev = fiber._ctcf_reverse
        chrom_state = fiber._chromatin_state

        N = n_simulations

        # Current arm positions (bead indices) for all simulations
        left_pos = np.full(N, dsb_bead, dtype=np.int32)
        right_pos = np.full(N, dsb_bead, dtype=np.int32)

        # Stall flags
        left_stalled = np.zeros(N, dtype=bool)
        right_stalled = np.zeros(N, dtype=bool)
        left_stalled_ctcf = np.zeros(N, dtype=bool)
        right_stalled_ctcf = np.zeros(N, dtype=bool)

        # Counters
        n_ctcf_bypassed = np.zeros(N, dtype=np.int32)
        extrusion_time = np.full(N, float(self.max_time_sec))

        # Accumulators for sub-bead stepping
        left_acc = np.zeros(N, dtype=np.float64)
        right_acc = np.zeros(N, dtype=np.float64)

        for t in range(1, self.max_time_sec + 1):
            both_stalled = left_stalled & right_stalled
            if np.all(both_stalled):
                break

            # Mark time for newly-both-stalled simulations
            newly_done = both_stalled & (extrusion_time == float(self.max_time_sec))
            # Actually we track extrusion_time at the end

            # ---- LEFT ARM ----
            active_left = ~left_stalled
            if np.any(active_left):
                next_left = left_pos - 1
                # Boundary check
                at_edge = next_left < 0
                left_stalled[active_left & at_edge] = True
                active_left = active_left & ~at_edge

                if np.any(active_left):
                    # Chromatin state at next bead
                    nl = np.clip(next_left, 0, n_beads - 1)
                    is_closed = chrom_state[nl] == 1

                    # Advance probability
                    advance_prob = np.where(
                        is_closed, self.heterochromatin_slowdown, 1.0
                    )
                    rand_advance = rng.random(N)
                    do_advance = active_left & (rand_advance < advance_prob)

                    # Add to accumulator
                    left_acc[do_advance] += beads_per_sec

                    # Process integer steps
                    can_step = (left_acc >= 1.0) & active_left
                    while np.any(can_step & ~left_stalled):
                        stepping = can_step & ~left_stalled
                        left_acc[stepping] -= 1.0
                        new_pos = left_pos[stepping] - 1

                        # Boundary
                        hit_edge = new_pos < 0
                        left_stalled[np.where(stepping)[0][hit_edge]] = True

                        # Update positions for non-edge
                        valid = ~hit_edge
                        valid_idx = np.where(stepping)[0][valid]
                        left_pos[valid_idx] = new_pos[valid]

                        # CTCF check at new positions
                        for idx in valid_idx:
                            bp = left_pos[idx]
                            if ctcf_rev[bp] > 0:
                                eff = min(
                                    1.0,
                                    ctcf_rev[bp] * self.ctcf_stall_probability,
                                )
                                if rng.random() < eff:
                                    left_stalled[idx] = True
                                    left_stalled_ctcf[idx] = True
                                else:
                                    n_ctcf_bypassed[idx] += 1

                        can_step = (left_acc >= 1.0) & ~left_stalled

            # ---- RIGHT ARM ----
            active_right = ~right_stalled
            if np.any(active_right):
                next_right = right_pos + 1
                at_edge = next_right >= n_beads
                right_stalled[active_right & at_edge] = True
                active_right = active_right & ~at_edge

                if np.any(active_right):
                    nr = np.clip(next_right, 0, n_beads - 1)
                    is_closed = chrom_state[nr] == 1
                    advance_prob = np.where(
                        is_closed, self.heterochromatin_slowdown, 1.0
                    )
                    rand_advance = rng.random(N)
                    do_advance = active_right & (rand_advance < advance_prob)
                    right_acc[do_advance] += beads_per_sec

                    can_step = (right_acc >= 1.0) & active_right
                    while np.any(can_step & ~right_stalled):
                        stepping = can_step & ~right_stalled
                        right_acc[stepping] -= 1.0
                        new_pos = right_pos[stepping] + 1

                        hit_edge = new_pos >= n_beads
                        right_stalled[np.where(stepping)[0][hit_edge]] = True

                        valid = ~hit_edge
                        valid_idx = np.where(stepping)[0][valid]
                        right_pos[valid_idx] = new_pos[valid]

                        for idx in valid_idx:
                            bp = right_pos[idx]
                            if ctcf_fwd[bp] > 0:
                                eff = min(
                                    1.0,
                                    ctcf_fwd[bp] * self.ctcf_stall_probability,
                                )
                                if rng.random() < eff:
                                    right_stalled[idx] = True
                                    right_stalled_ctcf[idx] = True
                                else:
                                    n_ctcf_bypassed[idx] += 1

                        can_step = (right_acc >= 1.0) & ~right_stalled

            # Update extrusion time for newly stopped simulations
            just_done = (left_stalled & right_stalled) & (
                extrusion_time == float(self.max_time_sec)
            )
            extrusion_time[just_done] = float(t)

        # Build individual ExtrusionResult objects
        results: List[ExtrusionResult] = []
        for i in range(N):
            l_bp = fiber.get_bead_position_bp(int(left_pos[i]))
            r_bp = fiber.get_bead_position_bp(int(right_pos[i]))
            results.append(
                ExtrusionResult(
                    dsb_position_bp=dsb_position_bp,
                    left_boundary_bp=l_bp,
                    right_boundary_bp=r_bp,
                    domain_size_bp=r_bp - l_bp,
                    left_stalled_at_ctcf=bool(left_stalled_ctcf[i]),
                    right_stalled_at_ctcf=bool(right_stalled_ctcf[i]),
                    extrusion_time_sec=float(extrusion_time[i]),
                    n_ctcf_bypassed=int(n_ctcf_bypassed[i]),
                    time_trace_left=None,
                    time_trace_right=None,
                )
            )

        return self._aggregate_results(results)

    # -----------------------------------------------------------------
    # Aggregation helper
    # -----------------------------------------------------------------

    @staticmethod
    def _aggregate_results(results: List[ExtrusionResult]) -> SimulationResults:
        """Compute summary statistics from a list of individual results.

        Parameters
        ----------
        results : list of ExtrusionResult
            Individual simulation outcomes.

        Returns
        -------
        SimulationResults
            Aggregated statistics.
        """
        n = len(results)
        domains = np.array([r.domain_size_bp for r in results], dtype=np.float64)
        lefts = np.array([r.left_boundary_bp for r in results], dtype=np.float64)
        rights = np.array([r.right_boundary_bp for r in results], dtype=np.float64)

        return SimulationResults(
            individual_results=results,
            n_simulations=n,
            domain_sizes_bp=domains,
            left_boundaries_bp=lefts,
            right_boundaries_bp=rights,
            mean_domain_size_bp=float(np.mean(domains)),
            median_domain_size_bp=float(np.median(domains)),
            std_domain_size_bp=float(np.std(domains)),
            fraction_left_stalled_ctcf=sum(
                1 for r in results if r.left_stalled_at_ctcf
            )
            / n,
            fraction_right_stalled_ctcf=sum(
                1 for r in results if r.right_stalled_at_ctcf
            )
            / n,
            mean_ctcf_bypassed=float(
                np.mean([r.n_ctcf_bypassed for r in results])
            ),
        )
