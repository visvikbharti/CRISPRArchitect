"""
LoopSimulator — Main User Interface for the Cohesin Extrusion Simulation
=========================================================================

This module provides ``LoopSimulator``, the primary entry point for
biologists who want to run cohesin loop extrusion simulations without
worrying about the internal machinery. It chains together:

1. **Fiber setup** — create a 1D chromatin lattice with CTCF sites
2. **DSB placement** — specify one or two DSB positions
3. **Simulation** — run Monte Carlo extrusion ensembles
4. **Analysis** — compute search domains, overlap probabilities
5. **Visualization** — plot results as genomic tracks and histograms

Typical Workflow
----------------
A molecular biologist planning a CRISPR experiment would use this module
to answer questions like:

- "What genomic region will be scanned for homology after my Cas9 cut?"
- "If I make two cuts 1 Mb apart, do their search domains overlap?"
  (Usually NO — confirming the risk of large deletions.)
- "Does it matter for my cssDNA donor whether the DSB is near a CTCF site?"
  (For exogenous donors: NO, because the donor finds the DSB by diffusion.)

Example
-------
>>> from crisprarchitect.loopsim import LoopSimulator
>>>
>>> # Simulate a 5 Mb region with a DSB at the center
>>> sim = LoopSimulator(region_start_bp=0, region_end_bp=5_000_000)
>>> sim.setup_fiber(ctcf_spacing=120_000, ctcf_strength=0.8)
>>> sim.add_dsb(2_500_000)
>>> results = sim.run(n_simulations=500)
>>> sim.summary()
>>> sim.plot_search_domain()
>>>
>>> # Dual-DSB scenario: two cuts 1 Mb apart
>>> sim2 = LoopSimulator(0, 5_000_000)
>>> sim2.setup_fiber()
>>> sim2.add_dsb(2_000_000)
>>> sim2.add_second_dsb(3_000_000)
>>> results2 = sim2.run(n_simulations=500)
>>> sim2.summary()  # Shows non-overlapping domains
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .chromatin_fiber import ChromatinFiber
from .cohesin_extruder import CohesinExtruder
from .homology_search import HomologySearchPredictor, SearchDomainResult
from .visualize import (
    plot_extrusion_kymograph,
    plot_scan_probability_track,
    plot_dual_dsb_comparison,
)


# ---------------------------------------------------------------------------
# Result data class
# ---------------------------------------------------------------------------

@dataclass
class LoopSimResults:
    """Complete results from a LoopSimulator run.

    Attributes
    ----------
    dsb_positions : list of int
        Genomic coordinates of all DSBs simulated.
    per_dsb_results : dict
        Mapping from DSB position (int) to its SearchDomainResult.
    dual_dsb_overlap : float or None
        If two DSBs were simulated, the fraction of simulations where their
        search domains overlap. None if only one DSB.
    dual_dsb_analysis : dict or None
        Full dual-DSB analysis dictionary (from HomologySearchPredictor).
        None if only one DSB.
    interpretation : str
        Human-readable biological interpretation of the results.
    """

    dsb_positions: List[int]
    per_dsb_results: Dict[int, SearchDomainResult]
    dual_dsb_overlap: Optional[float] = None
    dual_dsb_analysis: Optional[Dict] = field(default=None, repr=False)
    interpretation: str = ""


# ---------------------------------------------------------------------------
# LoopSimulator
# ---------------------------------------------------------------------------

class LoopSimulator:
    """Main interface for cohesin loop extrusion simulation.

    This class orchestrates the full simulation pipeline: building the
    chromatin fiber, running extrusion simulations, analyzing results,
    and generating visualizations. It is designed for molecular biologists
    who want predictions about the homology search domain at their DSB
    of interest, without needing to interact with the lower-level simulation
    engine directly.

    The simulation is based on the model from Marin-Gonzalez et al.
    (Science, 2025), which demonstrated that cohesin loaded at a DSB
    performs bidirectional loop extrusion to create a chromatin scanning
    domain for RAD51-mediated homology search. CTCF sites act as directional
    barriers that define the boundaries of this domain.

    Parameters
    ----------
    region_start_bp : int
        Left boundary of the genomic region to simulate (bp).
    region_end_bp : int
        Right boundary of the genomic region (bp).
    resolution_bp : int, optional
        Resolution of the chromatin lattice (bp per bead). Default 1000
        (1 kb per bead). This matches the extrusion rate of ~1 kb/sec,
        so each time step = one bead.

    Examples
    --------
    >>> sim = LoopSimulator(0, 5_000_000)
    >>> sim.setup_fiber(ctcf_spacing=100_000, ctcf_strength=0.8)
    >>> sim.add_dsb(2_500_000)
    >>> results = sim.run(n_simulations=1000)
    >>> sim.plot_search_domain()
    """

    def __init__(
        self,
        region_start_bp: int,
        region_end_bp: int,
        resolution_bp: int = 1000,
    ) -> None:
        self.region_start_bp: int = region_start_bp
        self.region_end_bp: int = region_end_bp
        self.resolution_bp: int = resolution_bp

        # Components (created lazily or by setup methods)
        self.fiber: Optional[ChromatinFiber] = None
        self.extruder: Optional[CohesinExtruder] = None
        self.predictor: Optional[HomologySearchPredictor] = None

        # DSB positions
        self._dsb_positions: List[int] = []

        # Results (populated by run())
        self.results: Optional[LoopSimResults] = None

        # Simulation parameters (set by run())
        self._n_simulations: int = 1000
        self._seed: Optional[int] = None

    # -----------------------------------------------------------------
    # Fiber setup
    # -----------------------------------------------------------------

    def setup_fiber(
        self,
        ctcf_spacing: int = 100_000,
        ctcf_strength: float = 0.8,
        ctcf_occupancy: float = 0.7,
        heterochromatin_regions: Optional[List[tuple]] = None,
        seed: Optional[int] = None,
    ) -> ChromatinFiber:
        """Initialize the chromatin fiber with CTCF sites and chromatin states.

        Creates a ChromatinFiber spanning the simulation region and populates
        it with CTCF sites at approximately uniform spacing (with stochastic
        occupancy). Optionally marks regions as heterochromatin.

        In the real genome, CTCF sites are NOT uniformly distributed — they
        cluster at TAD boundaries and are sparser within TAD interiors. For
        a locus-specific simulation, use the returned fiber object to add
        custom CTCF sites from ChIP-seq data. For exploratory simulations
        (e.g., "what is the typical search domain size?"), uniform placement
        is a reasonable approximation.

        Parameters
        ----------
        ctcf_spacing : int, optional
            Average spacing between potential CTCF sites (bp). Default
            100,000 (100 kb). Genomic average is ~50-200 kb.
        ctcf_strength : float, optional
            Mean barrier strength of CTCF sites. Default 0.8.
        ctcf_occupancy : float, optional
            Probability that each potential site is occupied. Default 0.7.
        heterochromatin_regions : list of (start, end) tuples, optional
            Regions to mark as heterochromatin. Each tuple is (start_bp,
            end_bp). Default None (all euchromatin).
        seed : int, optional
            Random seed for CTCF placement.

        Returns
        -------
        ChromatinFiber
            The configured fiber (also stored as self.fiber).
        """
        self.fiber = ChromatinFiber(
            start_bp=self.region_start_bp,
            end_bp=self.region_end_bp,
            resolution_bp=self.resolution_bp,
        )

        rng = np.random.default_rng(seed)
        n_placed = self.fiber.add_ctcf_sites_uniform(
            spacing_bp=ctcf_spacing,
            probability=ctcf_occupancy,
            strength_mean=ctcf_strength,
            strength_std=0.15,
            rng=rng,
        )

        if heterochromatin_regions:
            for start, end in heterochromatin_regions:
                self.fiber.set_chromatin_state(start, end, "closed")

        return self.fiber

    # -----------------------------------------------------------------
    # DSB placement
    # -----------------------------------------------------------------

    def add_dsb(self, position_bp: int) -> None:
        """Place a DSB at the specified genomic position.

        Cohesin will be loaded at this position during the simulation.
        In the Marin-Gonzalez et al. (2025) model, the DSB is the anchor
        point from which bidirectional loop extrusion initiates.

        Parameters
        ----------
        position_bp : int
            Genomic coordinate of the DSB (bp). Must be within the
            simulated region.

        Raises
        ------
        ValueError
            If position is outside the simulated region or a DSB is
            already placed at this position.
        """
        if position_bp < self.region_start_bp or position_bp >= self.region_end_bp:
            raise ValueError(
                f"DSB position {position_bp} is outside the simulated region "
                f"[{self.region_start_bp}, {self.region_end_bp})."
            )
        if position_bp in self._dsb_positions:
            raise ValueError(
                f"A DSB is already placed at position {position_bp}."
            )
        if len(self._dsb_positions) >= 2:
            raise ValueError(
                "Maximum 2 DSBs supported. Use add_dsb() for the first and "
                "add_second_dsb() for the second."
            )
        self._dsb_positions.append(position_bp)

    def add_second_dsb(self, position_bp: int) -> None:
        """Place a second DSB for dual-cut scenario analysis.

        In CRISPR experiments, dual-DSB strategies are sometimes used for:
        - Large segment deletion (excision of a gene or regulatory element)
        - Paired nickase strategies (two Cas9-D10A nicks -> effective DSB)
        - Replacement of a large segment with a donor

        The dual-DSB simulation will assess whether the search domains from
        the two DSBs overlap. For DSBs separated by more than ~500 kb (a
        typical TAD), the domains almost never overlap, meaning each DSB
        is repaired independently. This has important implications:

        - The intervening segment is at high risk of deletion (if both DSBs
          are repaired by NHEJ with deletion of the fragment between them).
        - Each DSB needs its own donor if HDR is desired at both sites.

        Parameters
        ----------
        position_bp : int
            Genomic coordinate of the second DSB.
        """
        if len(self._dsb_positions) == 0:
            raise ValueError("Place the first DSB with add_dsb() before adding a second.")
        if len(self._dsb_positions) >= 2:
            raise ValueError("A second DSB is already placed.")
        self.add_dsb(position_bp)

    # -----------------------------------------------------------------
    # Run simulation
    # -----------------------------------------------------------------

    def run(
        self,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
        extrusion_rate: float = 1.0,
        ctcf_stall_prob: float = 0.8,
        max_time_sec: int = 3600,
        heterochromatin_slowdown: float = 0.5,
    ) -> LoopSimResults:
        """Run the full simulation for all placed DSBs.

        Creates the CohesinExtruder with the specified parameters and runs
        a Monte Carlo ensemble for each DSB. If two DSBs are placed, also
        computes the overlap analysis.

        Parameters
        ----------
        n_simulations : int, optional
            Number of Monte Carlo simulations per DSB. Default 1000.
        seed : int, optional
            Random seed for reproducibility.
        extrusion_rate : float, optional
            Cohesin extrusion rate in kb/sec. Default 1.0.
        ctcf_stall_prob : float, optional
            Global CTCF stall probability multiplier. Default 0.8.
        max_time_sec : int, optional
            Maximum extrusion time. Default 3600 (1 hour).
        heterochromatin_slowdown : float, optional
            Slowdown factor in heterochromatin. Default 0.5.

        Returns
        -------
        LoopSimResults
            Complete results including per-DSB search domains and
            optional dual-DSB overlap analysis.

        Raises
        ------
        RuntimeError
            If fiber or DSBs have not been set up.
        """
        if self.fiber is None:
            raise RuntimeError(
                "Chromatin fiber not set up. Call setup_fiber() first."
            )
        if not self._dsb_positions:
            raise RuntimeError(
                "No DSBs placed. Call add_dsb() first."
            )

        self._n_simulations = n_simulations
        self._seed = seed

        # Create extruder and predictor
        self.extruder = CohesinExtruder(
            extrusion_rate_kb_per_sec=extrusion_rate,
            ctcf_stall_probability=ctcf_stall_prob,
            max_time_sec=max_time_sec,
            heterochromatin_slowdown=heterochromatin_slowdown,
            record_traces=False,
        )
        self.predictor = HomologySearchPredictor(extruder=self.extruder)

        # Run per-DSB simulations
        per_dsb: Dict[int, SearchDomainResult] = {}
        for i, dsb_pos in enumerate(self._dsb_positions):
            dsb_seed = (seed + i * n_simulations) if seed is not None else None
            domain = self.predictor.predict_search_domain(
                self.fiber,
                dsb_pos,
                n_simulations=n_simulations,
                seed=dsb_seed,
            )
            per_dsb[dsb_pos] = domain

        # Dual-DSB analysis
        dual_overlap: Optional[float] = None
        dual_analysis: Optional[Dict] = None
        if len(self._dsb_positions) == 2:
            dual_analysis = self.predictor.dual_dsb_overlap_probability(
                self._dsb_positions[0],
                self._dsb_positions[1],
                self.fiber,
                n_simulations=n_simulations,
                seed=seed,
            )
            dual_overlap = dual_analysis["any_overlap_prob"]

        # Build interpretation text
        interpretation = self._build_interpretation(per_dsb, dual_overlap)

        self.results = LoopSimResults(
            dsb_positions=list(self._dsb_positions),
            per_dsb_results=per_dsb,
            dual_dsb_overlap=dual_overlap,
            dual_dsb_analysis=dual_analysis,
            interpretation=interpretation,
        )

        return self.results

    # -----------------------------------------------------------------
    # Analysis methods
    # -----------------------------------------------------------------

    def can_single_search_reach_both_dsbs(
        self,
        dsb1_bp: Optional[int] = None,
        dsb2_bp: Optional[int] = None,
    ) -> float:
        """Check if a single search domain from DSB1 can reach DSB2.

        This answers the question: "If cohesin is loaded at DSB1, what
        fraction of the time does the extruded domain extend all the way
        to DSB2?"

        For DSBs separated by > 500 kb, this fraction is almost always
        very low (< 5%), confirming that distant DSBs have independent
        repair dynamics. This is a key finding of the Marin-Gonzalez
        et al. (2025) work applied to genome editing strategy.

        Parameters
        ----------
        dsb1_bp : int, optional
            First DSB position. If None, uses the first placed DSB.
        dsb2_bp : int, optional
            Second DSB position. If None, uses the second placed DSB.

        Returns
        -------
        float
            Fraction of simulations where DSB1's search domain includes
            DSB2's position.

        Raises
        ------
        RuntimeError
            If results have not been computed yet.
        """
        if self.results is None:
            raise RuntimeError("Run the simulation first with run().")

        if dsb1_bp is None:
            dsb1_bp = self._dsb_positions[0]
        if dsb2_bp is None:
            if len(self._dsb_positions) < 2:
                raise ValueError("Need two DSBs for this analysis.")
            dsb2_bp = self._dsb_positions[1]

        # Look up DSB2 bead in DSB1's scan probability
        if dsb1_bp in self.results.per_dsb_results:
            domain = self.results.per_dsb_results[dsb1_bp]
            dsb2_bead = self.fiber.get_bead_index(dsb2_bp)
            return float(domain.position_scan_probability[dsb2_bead])

        # If not precomputed, run a new prediction
        if self.predictor is None or self.fiber is None:
            raise RuntimeError("Set up and run the simulation first.")
        return self.predictor.probability_donor_in_domain(
            dsb1_bp, dsb2_bp, self.fiber, n_simulations=self._n_simulations
        )

    # -----------------------------------------------------------------
    # Visualization
    # -----------------------------------------------------------------

    def plot_search_domain(self, figsize: tuple = (14, 6), save_path: Optional[str] = None):
        """Plot the search domain as a genomic track.

        Creates a publication-quality figure showing:
        - Top panel: per-position scan probability (blue gradient)
        - CTCF sites marked as triangles (orientation shown by direction)
        - DSB position(s) marked with red vertical lines
        - For dual DSBs: both search domains overlaid

        Parameters
        ----------
        figsize : tuple, optional
            Figure size (width, height) in inches. Default (14, 6).
        save_path : str, optional
            If provided, save figure to this path.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        if self.results is None:
            raise RuntimeError("Run the simulation first with run().")

        return plot_scan_probability_track(
            results=self.results,
            fiber=self.fiber,
            figsize=figsize,
            save_path=save_path,
        )

    def plot_domain_size_distribution(
        self,
        figsize: tuple = (10, 5),
        save_path: Optional[str] = None,
    ):
        """Plot histograms of search domain sizes.

        Generates one histogram per DSB showing the distribution of
        extruded domain sizes across the Monte Carlo ensemble. The
        distribution width reflects cell-to-cell variability in CTCF
        occupancy and extrusion dynamics.

        Parameters
        ----------
        figsize : tuple, optional
            Figure size. Default (10, 5).
        save_path : str, optional
            If provided, save figure to this path.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt

        if self.results is None:
            raise RuntimeError("Run the simulation first with run().")

        n_dsbs = len(self.results.dsb_positions)
        fig, axes = plt.subplots(1, n_dsbs, figsize=figsize, squeeze=False)

        for i, dsb_pos in enumerate(self.results.dsb_positions):
            ax = axes[0, i]
            domain = self.results.per_dsb_results[dsb_pos]

            sizes_mb = domain.domain_sizes / 1e6

            ax.hist(
                sizes_mb,
                bins=50,
                color="#4A90D9",
                edgecolor="white",
                alpha=0.85,
                density=True,
            )
            ax.axvline(
                domain.mean_domain_size_bp / 1e6,
                color="#C0392B",
                linestyle="--",
                linewidth=1.5,
                label=f"Mean = {domain.mean_domain_size_bp / 1e6:.2f} Mb",
            )
            ax.axvline(
                domain.median_domain_size_bp / 1e6,
                color="#27AE60",
                linestyle=":",
                linewidth=1.5,
                label=f"Median = {domain.median_domain_size_bp / 1e6:.2f} Mb",
            )

            ax.set_xlabel("Search Domain Size (Mb)", fontsize=11)
            ax.set_ylabel("Density", fontsize=11)
            ax.set_title(
                f"DSB at {dsb_pos / 1e6:.2f} Mb\n"
                f"({domain.n_simulations} simulations)",
                fontsize=12,
            )
            ax.legend(fontsize=9, framealpha=0.9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.suptitle(
            "Cohesin Extrusion Search Domain Size Distribution",
            fontsize=14,
            fontweight="bold",
            y=1.02,
        )
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")

        return fig

    def plot_kymograph(
        self,
        dsb_position_bp: Optional[int] = None,
        n_traces: int = 5,
        figsize: tuple = (12, 6),
        save_path: Optional[str] = None,
    ):
        """Plot a kymograph (space-time diagram) of loop extrusion.

        Requires running a small number of simulations with trace recording
        enabled. This method creates a new extruder with ``record_traces=True``
        and runs ``n_traces`` simulations to generate the kymograph data.

        The kymograph shows genomic position (x-axis) vs. time (y-axis), with
        the extruded domain expanding outward from the DSB like a growing
        "V" shape. CTCF stall points appear as flat horizontal lines where
        the domain boundary stops advancing.

        Parameters
        ----------
        dsb_position_bp : int, optional
            DSB position. If None, uses the first placed DSB.
        n_traces : int, optional
            Number of extrusion traces to overlay. Default 5.
        figsize : tuple, optional
            Figure size. Default (12, 6).
        save_path : str, optional
            If provided, save figure to this path.

        Returns
        -------
        matplotlib.figure.Figure
            The generated kymograph figure.
        """
        if self.fiber is None:
            raise RuntimeError("Set up the fiber first with setup_fiber().")
        if dsb_position_bp is None:
            if not self._dsb_positions:
                raise RuntimeError("Place a DSB first with add_dsb().")
            dsb_position_bp = self._dsb_positions[0]

        # Create a trace-recording extruder
        trace_extruder = CohesinExtruder(
            extrusion_rate_kb_per_sec=(
                self.extruder.extrusion_rate_kb_per_sec
                if self.extruder
                else 1.0
            ),
            ctcf_stall_probability=(
                self.extruder.ctcf_stall_probability
                if self.extruder
                else 0.8
            ),
            max_time_sec=(
                self.extruder.max_time_sec if self.extruder else 3600
            ),
            heterochromatin_slowdown=(
                self.extruder.heterochromatin_slowdown
                if self.extruder
                else 0.5
            ),
            record_traces=True,
        )

        trace_results = []
        for i in range(n_traces):
            seed = (self._seed + 99999 + i) if self._seed is not None else None
            rng = np.random.default_rng(seed)
            result = trace_extruder.extrude(self.fiber, dsb_position_bp, rng=rng)
            trace_results.append(result)

        return plot_extrusion_kymograph(
            trace_results,
            fiber=self.fiber,
            figsize=figsize,
            save_path=save_path,
        )

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------

    def summary(self) -> str:
        """Print and return a comprehensive summary with biological interpretation.

        Generates a detailed text report covering:
        - Fiber configuration (region size, CTCF sites, chromatin states)
        - Per-DSB search domain statistics
        - Dual-DSB overlap analysis (if applicable)
        - Biological interpretation and recommendations

        Returns
        -------
        str
            The summary text (also printed to stdout).
        """
        if self.results is None:
            raise RuntimeError("Run the simulation first with run().")

        lines = []
        lines.append("=" * 72)
        lines.append("  LoopSim: Cohesin Loop Extrusion Simulation Results")
        lines.append("  Based on Marin-Gonzalez et al., Science, 2025")
        lines.append("=" * 72)
        lines.append("")

        # Fiber info
        if self.fiber is not None:
            lines.append(self.fiber.describe())
            lines.append("")

        # Per-DSB results
        for dsb_pos in self.results.dsb_positions:
            domain = self.results.per_dsb_results[dsb_pos]
            lines.append(f"--- DSB at {dsb_pos:,} bp ({dsb_pos / 1e6:.2f} Mb) ---")
            lines.append(
                f"  Search domain size:  "
                f"mean = {domain.mean_domain_size_bp / 1e6:.2f} Mb, "
                f"median = {domain.median_domain_size_bp / 1e6:.2f} Mb, "
                f"std = {domain.std_domain_size_bp / 1e6:.2f} Mb"
            )

            # Boundary statistics
            mean_left = np.mean(domain.left_boundaries)
            mean_right = np.mean(domain.right_boundaries)
            lines.append(
                f"  Mean boundaries:     "
                f"left = {mean_left / 1e6:.2f} Mb, "
                f"right = {mean_right / 1e6:.2f} Mb"
            )

            # CTCF stalling (from simulation_results)
            sr = domain.simulation_results
            lines.append(
                f"  Left stalled at CTCF:  {sr.fraction_left_stalled_ctcf:.0%}"
            )
            lines.append(
                f"  Right stalled at CTCF: {sr.fraction_right_stalled_ctcf:.0%}"
            )
            lines.append(
                f"  Mean CTCF bypassed:    {sr.mean_ctcf_bypassed:.1f} sites"
            )
            lines.append(
                f"  Simulations:           {domain.n_simulations}"
            )
            lines.append("")

        # Dual-DSB analysis
        if self.results.dual_dsb_overlap is not None:
            lines.append("--- Dual-DSB Overlap Analysis ---")
            dist = abs(
                self.results.dsb_positions[1] - self.results.dsb_positions[0]
            )
            lines.append(f"  Distance between DSBs: {dist:,} bp ({dist / 1e6:.2f} Mb)")
            lines.append(
                f"  Overlap probability:   {self.results.dual_dsb_overlap:.1%}"
            )
            if self.results.dual_dsb_analysis:
                lines.append(
                    f"  DSB1 reaches DSB2:     "
                    f"{self.results.dual_dsb_analysis['dsb1_reaches_dsb2_prob']:.1%}"
                )
                lines.append(
                    f"  DSB2 reaches DSB1:     "
                    f"{self.results.dual_dsb_analysis['dsb2_reaches_dsb1_prob']:.1%}"
                )
            lines.append("")

        # Interpretation
        lines.append("--- Biological Interpretation ---")
        lines.append(self.results.interpretation)
        lines.append("")
        lines.append("=" * 72)

        text = "\n".join(lines)
        print(text)
        return text

    # -----------------------------------------------------------------
    # Private helpers
    # -----------------------------------------------------------------

    def _build_interpretation(
        self,
        per_dsb: Dict[int, SearchDomainResult],
        dual_overlap: Optional[float],
    ) -> str:
        """Build a biological interpretation of the simulation results.

        Parameters
        ----------
        per_dsb : dict
            Per-DSB SearchDomainResult objects.
        dual_overlap : float or None
            Dual-DSB overlap probability, if applicable.

        Returns
        -------
        str
            Multi-line interpretation text.
        """
        parts = []

        # Single-DSB interpretation
        for dsb_pos, domain in per_dsb.items():
            mean_mb = domain.mean_domain_size_bp / 1e6
            parts.append(
                f"The DSB at {dsb_pos / 1e6:.2f} Mb generates a cohesin-extruded "
                f"search domain of ~{mean_mb:.1f} Mb (mean). This domain is bounded "
                f"by the nearest CTCF barriers. Within this domain, the RAD51 "
                f"filament can scan for homologous sequences along the chromatin "
                f"fiber as it is reeled through the DSB anchor point."
            )

        # Dual-DSB interpretation
        if dual_overlap is not None and len(self._dsb_positions) == 2:
            dist = abs(self._dsb_positions[1] - self._dsb_positions[0])
            dist_mb = dist / 1e6
            if dual_overlap < 0.05:
                parts.append(
                    f"\nDual-DSB Analysis: The two DSBs are {dist_mb:.2f} Mb apart "
                    f"with search domain overlap in only {dual_overlap:.1%} of cells. "
                    f"This means the DSBs will undergo INDEPENDENT repair in virtually "
                    f"all cells. The intervening {dist_mb:.1f} Mb segment is at high "
                    f"risk of deletion if both DSBs are repaired by end-joining. "
                    f"This confirms the literature finding that widely-spaced dual-DSB "
                    f"strategies carry significant risk of large chromosomal deletions."
                )
            elif dual_overlap < 0.5:
                parts.append(
                    f"\nDual-DSB Analysis: The two DSBs are {dist_mb:.2f} Mb apart. "
                    f"Their search domains overlap in {dual_overlap:.0%} of cells. "
                    f"This partial overlap suggests they are near the boundary of the "
                    f"same TAD. In cells where overlap occurs, coordinated repair may "
                    f"be possible, but in most cells the DSBs will be independent."
                )
            else:
                parts.append(
                    f"\nDual-DSB Analysis: The two DSBs are {dist_mb:.2f} Mb apart "
                    f"and their search domains overlap in {dual_overlap:.0%} of cells. "
                    f"They appear to be within the same TAD, so cohesin extruded from "
                    f"one DSB can frequently reach the other."
                )

        # General note about exogenous donors
        parts.append(
            "\nNote for exogenous donor (cssDNA/ssODN/AAV) experiments: "
            "The cohesin search domain does NOT limit exogenous donor access to "
            "the DSB. Exogenous donors reach the DSB by 3D diffusion through the "
            "nucleoplasm, independent of loop extrusion. The search domain is "
            "relevant primarily for sister chromatid recombination and loss-of-"
            "heterozygosity risk assessment."
        )

        return "\n".join(parts)
