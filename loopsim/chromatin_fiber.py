"""
Chromatin Fiber Model — 1D Lattice Representation of Interphase Chromatin
==========================================================================

This module models interphase chromatin as a **one-dimensional bead-on-string
lattice**, where each bead represents a fixed length of genomic DNA (default:
1 kb per bead). This coarse-grained representation captures the essential
features needed to simulate cohesin loop extrusion:

    - **Genomic position**: each bead maps to a contiguous stretch of the
      genome (e.g., bead 0 = chr1:0-999, bead 1 = chr1:1000-1999, ...).

    - **CTCF occupancy**: certain beads are occupied by the CCCTC-binding
      factor (CTCF). CTCF is a zinc-finger protein that binds a 19-bp core
      motif in an **oriented** manner (the motif is asymmetric). CTCF acts
      as a directional "roadblock" for cohesin: the extrusion motor stalls
      when it encounters a CTCF site whose N-terminus points toward the
      approaching cohesin. In the genome, convergent CTCF pairs (pointing
      inward toward each other) demarcate the boundaries of Topologically
      Associating Domains (TADs).

    - **CTCF orientation**: Each CTCF site is annotated as "forward" or
      "reverse." By convention:
        * "forward" = the CTCF motif is on the + strand, N-terminus pointing
          rightward. This blocks cohesin approaching from the LEFT (i.e.,
          the rightward-extruding arm).
        * "reverse" = motif on the - strand, N-terminus pointing leftward.
          This blocks cohesin approaching from the RIGHT (i.e., the
          leftward-extruding arm).

    - **CTCF strength**: Not all CTCF sites are equally effective barriers.
      In vivo, CTCF occupancy is stochastic (~60-80% at strong sites,
      lower at weaker ones). We model this as a probability (0-1) that the
      site will actually stall cohesin on any given encounter. A strength
      of 1.0 = always stalls; 0.5 = stalls half the time.

    - **Chromatin state**: Beads are labeled "open" (euchromatin) or
      "closed" (heterochromatin). Open chromatin permits normal extrusion
      speed. Closed/compacted chromatin slows the extrusion motor, because
      the densely packed nucleosomes present a physical obstacle.

Why a 1D lattice?
-----------------
The 3D folding of chromatin is important for many processes, but cohesin
loop extrusion is fundamentally a **1D process** — the motor walks along
the chromatin fiber in cis. The 3D spatial position of the fiber matters
for *inter-chromosomal* contacts (e.g., an exogenous cssDNA donor finding
the DSB by diffusion), but the extrusion itself is accurately captured by
a 1D model. This keeps the simulation fast enough for Monte Carlo
ensembles of thousands of runs.

Resolution
----------
The default resolution of 1 kb/bead is chosen because:
- Cohesin extrudes at ~1 kb/sec, so 1 bead = 1 second of extrusion.
- CTCF sites are spaced ~50-200 kb apart, so we need hundreds of beads
  to resolve a typical TAD (~500 kb - 2 Mb).
- Finer resolution (e.g., 200 bp = 1 nucleosome) would increase compute
  cost 5-fold with minimal gain for the loop extrusion question.

References
----------
- Rao SSP, et al. "A 3D map of the human genome at kilobase resolution..."
  Cell, 2014. (CTCF/cohesin define TAD boundaries.)
- Nora EP, et al. "Spatial partitioning of the regulatory landscape of the
  X-inactivation centre." Nature, 2012. (TAD discovery.)
- Fudenberg G, et al. "Formation of chromosomal domains by loop extrusion."
  Cell Reports, 2016. (Foundational loop extrusion model.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class CTCFSite:
    """A single CTCF binding site on the chromatin fiber.

    CTCF (CCCTC-binding factor) is an 11-zinc-finger protein that binds a
    ~20 bp DNA motif. The motif is **asymmetric**, so CTCF binding has an
    inherent orientation. This orientation is critical for loop extrusion:

    - Cohesin approaching a CTCF site from the **N-terminal side** of the
      bound CTCF protein is stalled. The N-terminal domain of CTCF directly
      interacts with the SA2-SCC1 subunits of the cohesin ring, creating
      a physical barrier.

    - Cohesin approaching from the **C-terminal side** can bypass the site
      (CTCF is "transparent" from this direction).

    In the standard genomic convention:
        - "forward" (+ strand motif) => stalls RIGHT-ward extrusion arm
        - "reverse" (- strand motif) => stalls LEFT-ward extrusion arm

    A pair of convergent CTCF sites (one forward, one reverse) creates a
    stable TAD boundary when they face each other:

        ← reverse ... forward →       (convergent pair = TAD boundary)

    Attributes
    ----------
    position_bp : int
        Genomic coordinate (in base pairs) of the CTCF binding site.
    orientation : str
        Either "forward" or "reverse". See above for the biological meaning.
    strength : float
        Probability (0.0 to 1.0) that this CTCF site will successfully stall
        cohesin on any given encounter. Reflects in vivo occupancy and the
        quality of the motif match.

        Typical values:
        - Strong constitutive sites (e.g., insulator elements): 0.8-1.0
        - Medium sites at TAD sub-boundaries: 0.4-0.7
        - Weak / cell-type-specific sites: 0.1-0.3

        Source: Li et al., Genome Biology, 2020 — CTCF occupancy measured
        by ChIP-seq correlates with loop extrusion barrier strength.
    """

    position_bp: int
    orientation: str  # "forward" or "reverse"
    strength: float = 0.8

    def __post_init__(self) -> None:
        if self.orientation not in ("forward", "reverse"):
            raise ValueError(
                f"CTCF orientation must be 'forward' or 'reverse', "
                f"got '{self.orientation}'. "
                f"Forward = + strand motif (blocks rightward extrusion). "
                f"Reverse = - strand motif (blocks leftward extrusion)."
            )
        if not 0.0 <= self.strength <= 1.0:
            raise ValueError(
                f"CTCF strength must be between 0 and 1, got {self.strength}. "
                f"This represents the probability of stalling cohesin."
            )


# ---------------------------------------------------------------------------
# ChromatinFiber class
# ---------------------------------------------------------------------------

class ChromatinFiber:
    """A 1D lattice model of a chromatin fiber for loop extrusion simulation.

    The fiber is represented as a NumPy array of ``n_beads`` elements, where
    each bead spans ``resolution_bp`` base pairs of genomic DNA. Overlaid on
    this array are annotations for CTCF binding sites and chromatin states.

    This is the "substrate" that the cohesin extruder walks along. Think of
    it as laying out the genome on a ruler and marking the important features:

        |------|------|------|------|------|------|------|------|------|
        0kb    1kb    2kb    3kb    4kb    5kb    6kb    7kb    8kb

                  ^CTCF(fwd)       ^CTCF(rev)         ^CTCF(fwd)
        [  open  ][  open  ][ closed ][ closed ][  open  ][  open  ]

    When cohesin is loaded at a DSB and begins extruding, it walks left and
    right along this fiber, checking each bead for:
    - CTCF sites that might stall the motor
    - Chromatin state that might slow the motor

    Parameters
    ----------
    start_bp : int
        Genomic coordinate of the left end of the modeled region (bp).
    end_bp : int
        Genomic coordinate of the right end of the modeled region (bp).
    resolution_bp : int, optional
        Number of base pairs per bead. Default is 1000 (1 kb). At this
        resolution, each simulation time step (1 second) corresponds to
        one bead of extrusion, because cohesin extrudes at ~1 kb/s.

    Attributes
    ----------
    n_beads : int
        Total number of beads in the fiber lattice.
    resolution : int
        Resolution in base pairs per bead.
    ctcf_sites : list of CTCFSite
        All CTCF sites placed on this fiber.

    Examples
    --------
    >>> fiber = ChromatinFiber(start_bp=0, end_bp=5_000_000, resolution_bp=1000)
    >>> fiber.n_beads
    5000
    >>> fiber.add_ctcf_site(position_bp=1_000_000, orientation="forward", strength=0.9)
    >>> fiber.add_ctcf_site(position_bp=2_000_000, orientation="reverse", strength=0.7)
    >>> fiber.set_chromatin_state(start_bp=1_500_000, end_bp=1_800_000, state="closed")
    """

    def __init__(
        self,
        start_bp: int,
        end_bp: int,
        resolution_bp: int = 1000,
    ) -> None:
        if end_bp <= start_bp:
            raise ValueError(
                f"end_bp ({end_bp}) must be greater than start_bp ({start_bp})."
            )
        if resolution_bp <= 0:
            raise ValueError(
                f"resolution_bp must be positive, got {resolution_bp}."
            )

        self.start_bp: int = start_bp
        self.end_bp: int = end_bp
        self.resolution_bp: int = resolution_bp

        # Number of beads in the lattice
        self._n_beads: int = (end_bp - start_bp) // resolution_bp

        # -----------------------------------------------------------------
        # Internal arrays (one element per bead)
        # -----------------------------------------------------------------

        # Chromatin state: 0 = open (euchromatin), 1 = closed (heterochromatin)
        # Default: all open.
        self._chromatin_state: np.ndarray = np.zeros(self._n_beads, dtype=np.int8)

        # CTCF binding: per-bead arrays for fast lookup during simulation.
        # _ctcf_forward[i] = stall strength if bead i has a forward CTCF site, else 0.
        # _ctcf_reverse[i] = stall strength if bead i has a reverse CTCF site, else 0.
        self._ctcf_forward: np.ndarray = np.zeros(self._n_beads, dtype=np.float64)
        self._ctcf_reverse: np.ndarray = np.zeros(self._n_beads, dtype=np.float64)

        # Keep a human-readable list of all CTCF site objects
        self._ctcf_sites: List[CTCFSite] = []

    # -----------------------------------------------------------------
    # Properties
    # -----------------------------------------------------------------

    @property
    def n_beads(self) -> int:
        """Total number of beads in the fiber lattice."""
        return self._n_beads

    @property
    def resolution(self) -> int:
        """Resolution in base pairs per bead."""
        return self.resolution_bp

    @property
    def ctcf_sites(self) -> List[CTCFSite]:
        """List of all CTCF sites placed on this fiber (read-only copy)."""
        return list(self._ctcf_sites)

    @property
    def chromatin_state(self) -> np.ndarray:
        """Per-bead chromatin state array (0 = open, 1 = closed). Read-only copy."""
        return self._chromatin_state.copy()

    @property
    def ctcf_forward(self) -> np.ndarray:
        """Per-bead forward-CTCF stall strength. Read-only copy."""
        return self._ctcf_forward.copy()

    @property
    def ctcf_reverse(self) -> np.ndarray:
        """Per-bead reverse-CTCF stall strength. Read-only copy."""
        return self._ctcf_reverse.copy()

    # -----------------------------------------------------------------
    # Coordinate conversion
    # -----------------------------------------------------------------

    def get_bead_index(self, position_bp: int) -> int:
        """Convert a genomic coordinate (bp) to a bead index.

        The bead index is a zero-based integer indicating which bead in the
        lattice contains the given genomic position.

        Parameters
        ----------
        position_bp : int
            Genomic coordinate in base pairs.

        Returns
        -------
        int
            Bead index (0-based). Clipped to [0, n_beads - 1].

        Raises
        ------
        ValueError
            If position_bp is outside the modeled region.

        Examples
        --------
        >>> fiber = ChromatinFiber(start_bp=100_000, end_bp=200_000, resolution_bp=1000)
        >>> fiber.get_bead_index(150_000)
        50
        """
        if position_bp < self.start_bp or position_bp >= self.end_bp:
            raise ValueError(
                f"Position {position_bp} bp is outside the modeled region "
                f"[{self.start_bp}, {self.end_bp})."
            )
        idx = (position_bp - self.start_bp) // self.resolution_bp
        return int(np.clip(idx, 0, self._n_beads - 1))

    def get_bead_position_bp(self, bead_index: int) -> int:
        """Convert a bead index back to a genomic coordinate (bp).

        Returns the start coordinate of the genomic interval represented
        by this bead.

        Parameters
        ----------
        bead_index : int
            Zero-based bead index.

        Returns
        -------
        int
            Genomic coordinate (bp) at the start of this bead's interval.
        """
        if bead_index < 0 or bead_index >= self._n_beads:
            raise ValueError(
                f"Bead index {bead_index} is out of range [0, {self._n_beads})."
            )
        return self.start_bp + bead_index * self.resolution_bp

    # -----------------------------------------------------------------
    # CTCF site placement
    # -----------------------------------------------------------------

    def add_ctcf_site(
        self,
        position_bp: int,
        orientation: str,
        strength: float = 0.8,
    ) -> None:
        """Place a single CTCF binding site on the fiber.

        CTCF sites are the "stop signs" for cohesin loop extrusion. In the
        mammalian genome, there are ~50,000-80,000 CTCF sites, spaced every
        ~50-200 kb on average, but with highly non-uniform spacing. Strong
        CTCF sites at TAD boundaries have high occupancy (strength ~ 0.8-1.0),
        while weaker sites within TADs have lower occupancy.

        The **orientation** determines which direction of extrusion is blocked:
        - "forward" (+ strand motif): blocks the RIGHT-ward extrusion arm.
          Cohesin extruding to the right encounters this CTCF and stalls.
        - "reverse" (- strand motif): blocks the LEFT-ward extrusion arm.
          Cohesin extruding to the left encounters this CTCF and stalls.

        A convergent pair (one reverse on the left, one forward on the right)
        traps cohesin between them, creating a stable loop / TAD:

            ←[reverse CTCF]........loop........[forward CTCF]→

        Parameters
        ----------
        position_bp : int
            Genomic coordinate of the CTCF site (bp).
        orientation : str
            "forward" or "reverse".
        strength : float, optional
            Probability [0, 1] that this site stalls cohesin. Default 0.8.
        """
        site = CTCFSite(
            position_bp=position_bp,
            orientation=orientation,
            strength=strength,
        )
        self._ctcf_sites.append(site)

        # Update the per-bead lookup arrays
        idx = self.get_bead_index(position_bp)
        if orientation == "forward":
            self._ctcf_forward[idx] = max(self._ctcf_forward[idx], strength)
        else:
            self._ctcf_reverse[idx] = max(self._ctcf_reverse[idx], strength)

    def add_ctcf_sites_uniform(
        self,
        spacing_bp: int = 100_000,
        probability: float = 0.7,
        strength_mean: float = 0.7,
        strength_std: float = 0.15,
        rng: Optional[np.random.Generator] = None,
    ) -> int:
        """Generate CTCF sites at approximately uniform spacing with random orientation.

        In the real genome, CTCF sites are not uniformly spaced — they cluster
        at TAD boundaries and are sparser within TAD interiors. However, for a
        first-approximation simulation (or when the user does not have ChIP-seq
        data for their locus), uniform placement with stochastic occupancy is a
        reasonable starting model.

        At each potential site (placed every ``spacing_bp``), the site is
        occupied with ``probability`` (mimicking the stochastic nature of CTCF
        binding in a cell population). Orientation is assigned randomly (50/50
        forward/reverse), reflecting the roughly equal prevalence of each
        orientation genome-wide.

        Parameters
        ----------
        spacing_bp : int, optional
            Average distance between potential CTCF sites (bp). Default
            100,000 (100 kb). In the human genome, CTCF sites are spaced
            ~50-200 kb on average. TAD-boundary CTCF sites are typically
            ~200-500 kb apart.
        probability : float, optional
            Probability (0-1) that each potential site is actually occupied.
            Default 0.7. Reflects the fact that CTCF binding is stochastic
            in individual cells — a given site might be bound in 70% of cells
            in a population.
        strength_mean : float, optional
            Mean stall strength for occupied sites. Default 0.7.
        strength_std : float, optional
            Std dev of stall strength. Default 0.15.
        rng : numpy.random.Generator, optional
            Random number generator for reproducibility. If None, uses the
            default numpy RNG.

        Returns
        -------
        int
            Number of CTCF sites actually placed (after stochastic filtering).

        Examples
        --------
        >>> fiber = ChromatinFiber(0, 5_000_000)
        >>> n_sites = fiber.add_ctcf_sites_uniform(spacing_bp=120_000, probability=0.7)
        >>> print(f"Placed {n_sites} CTCF sites")
        """
        if rng is None:
            rng = np.random.default_rng()

        n_placed = 0
        # Generate candidate positions
        positions = np.arange(
            self.start_bp + spacing_bp,
            self.end_bp,
            spacing_bp,
        )

        for pos in positions:
            # Stochastic occupancy
            if rng.random() > probability:
                continue

            # Random orientation (50/50)
            orientation = "forward" if rng.random() < 0.5 else "reverse"

            # Random strength drawn from a clipped normal distribution
            strength = float(
                np.clip(rng.normal(strength_mean, strength_std), 0.05, 1.0)
            )

            self.add_ctcf_site(
                position_bp=int(pos),
                orientation=orientation,
                strength=strength,
            )
            n_placed += 1

        return n_placed

    # -----------------------------------------------------------------
    # Chromatin state
    # -----------------------------------------------------------------

    def set_chromatin_state(
        self,
        start_bp: int,
        end_bp: int,
        state: str,
    ) -> None:
        """Mark a genomic region as "open" (euchromatin) or "closed" (heterochromatin).

        Chromatin state affects the speed of cohesin extrusion:

        - **Open chromatin** (euchromatin): Loosely packed nucleosomes with
          histone modifications like H3K4me3, H3K27ac. Cohesin extrudes at
          its normal rate (~1 kb/sec) through open chromatin.

        - **Closed chromatin** (heterochromatin): Densely packed nucleosomes
          with H3K9me3, H3K27me3 (Polycomb). The compact fiber presents a
          physical obstacle, slowing the cohesin motor. In the simulation,
          heterochromatin reduces the extrusion advance probability per time
          step by a configurable factor (default 50%).

        This effect is biologically motivated by single-molecule studies
        showing that cohesin extrusion is slower on reconstituted chromatin
        with higher nucleosome density (Kim Y, et al., Science, 2019).

        By default, the entire fiber is initialized as "open."

        Parameters
        ----------
        start_bp : int
            Start of the region to annotate (bp).
        end_bp : int
            End of the region to annotate (bp).
        state : str
            "open" or "closed".

        Raises
        ------
        ValueError
            If state is not "open" or "closed".

        Examples
        --------
        >>> fiber = ChromatinFiber(0, 5_000_000)
        >>> # Mark a 500 kb heterochromatin block
        >>> fiber.set_chromatin_state(2_000_000, 2_500_000, "closed")
        """
        if state not in ("open", "closed"):
            raise ValueError(
                f"Chromatin state must be 'open' or 'closed', got '{state}'. "
                f"Open = euchromatin (normal extrusion speed). "
                f"Closed = heterochromatin (slower extrusion)."
            )

        # Clamp to fiber boundaries
        start_bp = max(start_bp, self.start_bp)
        end_bp = min(end_bp, self.end_bp)

        start_idx = self.get_bead_index(start_bp)
        end_idx = self.get_bead_index(end_bp - 1) + 1  # inclusive end
        end_idx = min(end_idx, self._n_beads)

        state_val = 0 if state == "open" else 1
        self._chromatin_state[start_idx:end_idx] = state_val

    # -----------------------------------------------------------------
    # Utility / display
    # -----------------------------------------------------------------

    def __repr__(self) -> str:
        region_mb = (self.end_bp - self.start_bp) / 1e6
        n_fwd = int(np.sum(self._ctcf_forward > 0))
        n_rev = int(np.sum(self._ctcf_reverse > 0))
        pct_closed = 100.0 * np.mean(self._chromatin_state)
        return (
            f"ChromatinFiber("
            f"region={region_mb:.1f} Mb, "
            f"beads={self._n_beads}, "
            f"resolution={self.resolution_bp} bp/bead, "
            f"CTCF_fwd={n_fwd}, CTCF_rev={n_rev}, "
            f"heterochromatin={pct_closed:.1f}%)"
        )

    def describe(self) -> str:
        """Return a human-readable summary of the fiber and its annotations.

        Useful for quick sanity checks before running a simulation.

        Returns
        -------
        str
            Multi-line description of the fiber.
        """
        region_size = self.end_bp - self.start_bp
        n_fwd = int(np.sum(self._ctcf_forward > 0))
        n_rev = int(np.sum(self._ctcf_reverse > 0))
        n_closed = int(np.sum(self._chromatin_state == 1))
        pct_closed = 100.0 * n_closed / self._n_beads if self._n_beads > 0 else 0.0

        lines = [
            "=== Chromatin Fiber Summary ===",
            f"Genomic region:     {self.start_bp:,} - {self.end_bp:,} bp "
            f"({region_size / 1e6:.2f} Mb)",
            f"Resolution:         {self.resolution_bp:,} bp per bead",
            f"Number of beads:    {self._n_beads:,}",
            f"CTCF sites:         {len(self._ctcf_sites)} total "
            f"({n_fwd} forward, {n_rev} reverse)",
            f"Heterochromatin:    {n_closed:,} beads ({pct_closed:.1f}% of fiber)",
        ]
        if self._ctcf_sites:
            strengths = [s.strength for s in self._ctcf_sites]
            lines.append(
                f"CTCF strength:      mean={np.mean(strengths):.2f}, "
                f"range=[{np.min(strengths):.2f}, {np.max(strengths):.2f}]"
            )
        return "\n".join(lines)
