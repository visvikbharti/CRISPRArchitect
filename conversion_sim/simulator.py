"""
simulator.py — Main Monte Carlo Engine for HDR Gene Conversion Tract Simulation
================================================================================

This module chains together the three mechanistic sub-models (resection,
filament formation, synthesis) into a complete end-to-end simulation of
HDR-mediated gene conversion after a CRISPR-induced DSB.

Overview of the biological pathway being simulated
---------------------------------------------------

When a CRISPR nuclease creates a double-strand break (DSB), the cell
activates the DNA Damage Response (DDR).  Cells in S or G2 phase of the
cell cycle preferentially repair DSBs by **Homology-Directed Repair (HDR)**
— the only pathway that can incorporate user-provided donor-template
information into the genome.  The HDR sub-pathway we model here is
**Synthesis-Dependent Strand Annealing (SDSA)**, the dominant form of HDR
in mitotic mammalian cells.

SDSA proceeds through the following steps, each of which is simulated by
a dedicated sub-module:

1. **End resection** (``resection.py``)
   5'→3' nucleases (MRN/CtIP, then EXO1 or BLM-DNA2) chew back the
   5'-terminated strands on both sides of the DSB, exposing 3' single-
   stranded DNA (ssDNA) tails.

2. **RAD51 filament formation** (``filament.py``)
   BRCA2 loads RAD51 monomers onto the 3' ssDNA tails, forming a helical
   nucleoprotein filament that is the engine of homology search.

3. **Homology search and strand invasion**
   The RAD51 filament samples dsDNA sequences throughout the nucleus.
   When it encounters the homologous donor template, it catalyses strand
   invasion: the 3' ssDNA tail invades the donor duplex (or displaces the
   complementary strand in ssDNA donors), forming a D-loop.

   We model invasion success as a probability that depends on:
   - Whether the filament is long enough (≥ 15 nt minimum)
   - Whether there is sufficient homology between the filament and the
     donor arm (filament length vs. homology arm length)
   - The cell type's baseline HDR efficiency

4. **Template-directed DNA synthesis** (``synthesis.py``)
   DNA Pol δ extends the invading 3' end by copying the donor template.
   The length of synthesis = the gene conversion tract length.

5. **D-loop collapse and strand annealing**
   Helicases (RTEL1, BLM) dismantle the D-loop.  The newly synthesised
   strand, now carrying donor-derived sequence, anneals back to the
   second resected end.  Gap-filling and ligation complete the repair.

The output of the simulation is a distribution of gene conversion tract
lengths, representing the range of outcomes across a population of cells.

Architecture notes
------------------
All per-simulation computations are vectorised with NumPy.  The ``run()``
method processes all *n* simulations in parallel arrays — no Python-level
for-loops over individual cells.  This makes 100 000-simulation runs
complete in well under 1 second.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Literal, Optional, Union

import numpy as np

# ---------------------------------------------------------------------------
# Import sub-models
# ---------------------------------------------------------------------------
from .resection import ResectionSimulator
from .filament import FilamentModel
from .synthesis import SynthesisSimulator

# ---------------------------------------------------------------------------
# Import biological constants
# ---------------------------------------------------------------------------
import sys, os
_PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

from utils.constants import (
    MIN_HOMOLOGY_FOR_INVASION_BP,
    CELL_TYPE_PARAMS,
    DONOR_TOPOLOGY_MULTIPLIER,
    HDR_ENHANCEMENT_PER_BP_OVERHANG,
)


@dataclass
class SimulationResults:
    """Container for the full output of a ConversionSimulator run.

    This dataclass holds *per-simulation* arrays (one element per Monte
    Carlo trial) plus summary metadata.

    Attributes
    ----------
    hdr_success : np.ndarray, shape (n,), dtype bool
        Whether HDR succeeded in each simulation.  A simulation is counted
        as "HDR success" if the RAD51 filament was long enough, there was
        sufficient homology overlap, strand invasion succeeded, and at
        least ``CONVERSION_TRACT_MIN_BP`` of donor sequence was synthesised.
    tract_lengths_bp : np.ndarray, shape (n,), dtype float64
        Gene conversion tract length (bp) for each simulation.  Set to 0
        for simulations where HDR failed.
    resection_left_bp : np.ndarray, shape (n,), dtype float64
        Total resection on the left side of the DSB (bp).
    resection_right_bp : np.ndarray, shape (n,), dtype float64
        Total resection on the right side of the DSB (bp).
    filament_length_nt : np.ndarray, shape (n,), dtype float64
        RAD51 filament length on the invading strand (nt).
    invasion_success : np.ndarray, shape (n,), dtype bool
        Whether strand invasion succeeded (prerequisite for synthesis).
    n_simulations : int
        Total number of Monte Carlo trials.
    parameters : dict
        Dictionary recording all input parameters for reproducibility.
    """

    hdr_success: np.ndarray
    tract_lengths_bp: np.ndarray
    resection_left_bp: np.ndarray
    resection_right_bp: np.ndarray
    filament_length_nt: np.ndarray
    invasion_success: np.ndarray
    n_simulations: int
    parameters: Dict


class ConversionSimulator:
    """End-to-end Monte Carlo simulator of HDR gene conversion tracts.

    This is the main entry point for the ConversionSim package.  It
    orchestrates the three mechanistic sub-models (resection, filament,
    synthesis) to produce a population-level distribution of gene
    conversion tract lengths for a given experimental configuration.

    Parameters
    ----------
    cut_type : {"blunt", "staggered_5prime"}
        Geometry of the DSB produced by the nuclease.

        * ``"blunt"`` — flush-ended DSB, e.g. SpCas9.  Both strands are
          cut at (approximately) the same position.

        * ``"staggered_5prime"`` — the non-target strand is cut further
          from the PAM, producing a 5' overhang.  Examples: Cas12a (~5 bp),
          enFnCas9 (~3 bp), vCas9 (~6 bp).

    overhang_length : int
        Length of the 5' overhang in bp.  Must be 0 for blunt cuts and > 0
        for staggered cuts.

    donor_topology : {"circular_ssDNA", "linear_ssDNA", "linear_dsDNA", "plasmid"}
        Physical form of the donor template.  Affects D-loop stability
        (and hence tract length) and the baseline HDR efficiency multiplier.

    homology_arm_length : int
        Length of each homology arm on the donor template, in bp.
        Default 300 bp, which is optimal for cssDNA donors
        (Iyer et al., CRISPR J. 2022).

    cell_type : str
        Cell type being edited.  Must be a key in ``CELL_TYPE_PARAMS``
        from constants.py (e.g., "iPSC", "HEK293T", "K562", "T_cell",
        "HSC").  Determines baseline HDR efficiency and cell-cycle
        parameters.

    n_simulations : int
        Number of Monte Carlo trials.  Default 10 000.  For publication-
        quality statistics, use 50 000–100 000.

    seed : int or None
        Master random seed.  Sub-models derive their seeds from this.

    Examples
    --------
    >>> sim = ConversionSimulator(
    ...     cut_type="staggered_5prime",
    ...     overhang_length=3,
    ...     donor_topology="circular_ssDNA",
    ...     homology_arm_length=300,
    ...     cell_type="iPSC",
    ...     n_simulations=50_000,
    ...     seed=42,
    ... )
    >>> results = sim.run()
    >>> sim.summary()
    """

    def __init__(
        self,
        cut_type: Literal["blunt", "staggered_5prime"] = "blunt",
        overhang_length: int = 0,
        donor_topology: Literal[
            "circular_ssDNA", "linear_ssDNA", "linear_dsDNA", "plasmid"
        ] = "linear_dsDNA",
        homology_arm_length: int = 300,
        cell_type: str = "HEK293T",
        n_simulations: int = 10_000,
        seed=None,
    ) -> None:
        # ---- Validate inputs ----
        if cut_type not in ("blunt", "staggered_5prime"):
            raise ValueError(f"Invalid cut_type: {cut_type!r}")
        if cut_type == "blunt" and overhang_length != 0:
            raise ValueError("overhang_length must be 0 for blunt cuts")
        if cut_type == "staggered_5prime" and overhang_length <= 0:
            raise ValueError("overhang_length must be > 0 for staggered cuts")
        if donor_topology not in (
            "circular_ssDNA", "linear_ssDNA", "linear_dsDNA", "plasmid"
        ):
            raise ValueError(f"Invalid donor_topology: {donor_topology!r}")
        if cell_type not in CELL_TYPE_PARAMS:
            raise ValueError(
                f"Unknown cell_type {cell_type!r}. "
                f"Available: {list(CELL_TYPE_PARAMS.keys())}"
            )
        if homology_arm_length < 20:
            raise ValueError("homology_arm_length must be >= 20 bp")

        self.cut_type = cut_type
        self.overhang_length = overhang_length
        self.donor_topology = donor_topology
        self.homology_arm_length = homology_arm_length
        self.cell_type = cell_type
        self.n_simulations = n_simulations
        self.seed = seed

        # ---- Retrieve cell-type parameters ----
        self._cell_params = CELL_TYPE_PARAMS[cell_type]

        # ---- Create master RNG and derive child seeds ----
        master_rng = np.random.default_rng(seed)
        child_seeds = master_rng.integers(0, 2**31, size=3)

        # ---- Instantiate sub-models ----
        self._resection_sim = ResectionSimulator(seed=int(child_seeds[0]))
        self._filament_model = FilamentModel(seed=int(child_seeds[1]))
        self._synthesis_sim = SynthesisSimulator(seed=int(child_seeds[2]))

        # ---- Internal RNG for this class's own stochastic decisions ----
        self._rng = np.random.default_rng(
            master_rng.integers(0, 2**31)
        )

        # ---- Results placeholder ----
        self._results: Optional[SimulationResults] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> SimulationResults:
        """Execute the full Monte Carlo simulation.

        For each of *n_simulations* independent trials, the following
        biological steps are simulated:

        1. **End resection**: Draw resection lengths for both sides of the
           DSB.  The invading strand uses the *longer* resected side (the
           cell uses whichever 3' tail finds the donor first; biased toward
           the longer tail because it has a larger RAD51 filament and
           greater homology-search reach).

        2. **RAD51 filament formation**: Determine how much of the resected
           ssDNA is coated by RAD51.

        3. **Homology check**: Verify that:
           (a) The filament is at least ``MIN_HOMOLOGY_FOR_INVASION_BP``
               long (otherwise the filament is too short for D-loop
               formation).
           (b) The filament does not extend far beyond the homology arm on
               the donor.  If the filament is longer than the homology arm,
               only the portion overlapping the arm can base-pair with the
               donor; the excess filament hangs off the end.  This reduces
               invasion efficiency.  We model this as: if filament >
               homology arm, invasion probability drops linearly from 1.0
               to 0.3 as the excess grows.

        4. **Strand invasion**: A probabilistic event.  Base probability
           comes from the cell type's HDR efficiency, modified by donor
           topology and cut-stagger effects.  The homology check from
           step 3 further modulates this probability.

        5. **DNA synthesis (tract length)**: For cells where invasion
           succeeds, draw the gene conversion tract length from the SDSA
           model.

        Returns
        -------
        SimulationResults
            Dataclass with per-simulation arrays and metadata.
        """
        n = self.n_simulations

        # ==================================================================
        # STEP 1: End resection
        # ==================================================================
        # Simulate the 5'->3' degradation of the broken DNA ends.
        # This produces 3' ssDNA overhangs that are the substrate for
        # RAD51 filament formation.
        resection_result = self._resection_sim.simulate(
            cut_type=self.cut_type,
            stagger_bp=self.overhang_length if self.cut_type == "staggered_5prime" else 0,
            n_simulations=n,
        )

        # For the purpose of strand invasion, we use the LONGER resected
        # side.  Biologically, either end can invade the donor template,
        # and the longer filament has a kinetic advantage (larger search
        # volume, more stable D-loop).  In reality the cell doesn't
        # "choose" — both ends search simultaneously and the first
        # productive invasion wins.  Using the max is a simplification
        # that captures this bias.
        invading_resection = np.maximum(
            resection_result.left_bp,
            resection_result.right_bp,
        )

        # ==================================================================
        # STEP 2: RAD51 filament formation
        # ==================================================================
        # The RAD51 filament assembles on the invading 3' ssDNA tail.
        # Not all of the ssDNA is coated (stochastic coverage 70-95 %).
        filament_result = self._filament_model.form_filament(invading_resection)

        # ==================================================================
        # STEP 3: Homology check
        # ==================================================================
        # Two conditions must be met for strand invasion:
        #
        # (a) MINIMUM FILAMENT LENGTH
        #     The RAD51 filament must be at least MIN_HOMOLOGY_FOR_INVASION_BP
        #     (~15 nt) long.  Shorter filaments cannot form a stable D-loop.
        filament_long_enough = filament_result.filament_complete  # bool array

        # (b) FILAMENT vs. HOMOLOGY ARM OVERLAP
        #     The filament can only base-pair with the portion of the donor
        #     that is homologous to the break-flanking sequence.  This is
        #     the homology arm.  If the filament is shorter than the arm,
        #     no problem — there is plenty of donor sequence to pair with.
        #     If the filament is LONGER than the arm, the excess filament
        #     "hangs off" the end of the donor arm and cannot form
        #     heteroduplex.  This excess destabilises the D-loop and
        #     reduces invasion efficiency.
        #
        #     Model: invasion efficiency scales from 1.0 (filament <= arm)
        #     down to 0.3 (filament >> arm) using a sigmoid-like function.
        arm_len = float(self.homology_arm_length)
        filament_nt = filament_result.filament_length_nt

        # Compute the ratio of excess filament to arm length.
        # excess_ratio = 0 when filament <= arm, increases as filament
        # extends beyond the arm.
        excess = np.maximum(filament_nt - arm_len, 0.0)
        excess_ratio = excess / arm_len  # normalised to arm length

        # Homology overlap factor: 1.0 when ratio=0, decays to 0.3 as
        # excess grows.  We use an exponential decay: factor = 0.3 + 0.7 * exp(-2 * excess_ratio).
        # At excess_ratio = 0: factor = 1.0
        # At excess_ratio = 1 (filament is 2x arm): factor ~ 0.39
        # At excess_ratio = 2 (filament is 3x arm): factor ~ 0.31
        homology_overlap_factor = 0.3 + 0.7 * np.exp(-2.0 * excess_ratio)

        # ==================================================================
        # STEP 4: Strand invasion probability
        # ==================================================================
        # The probability that the RAD51 filament successfully invades the
        # donor template depends on multiple factors:
        #
        # (i)   Cell type baseline HDR efficiency
        #       Different cell types have very different HDR rates, due to
        #       differences in cell-cycle distribution (only S/G2 cells do
        #       HDR), BRCA2/RAD51 expression levels, and p53 status.
        #       From constants.py CELL_TYPE_PARAMS.
        base_hdr = self._cell_params["hdr_base_efficiency"]

        # (ii)  Donor topology multiplier
        #       Some donor forms are more effective templates.  cssDNA is
        #       ~3x better than linear dsDNA because:
        #       - Exonuclease resistance (longer intracellular half-life)
        #       - Single-stranded: no need for donor strand separation
        #       - May be preferentially channelled into the HDR pathway
        #       From constants.py DONOR_TOPOLOGY_MULTIPLIER.
        topology_key = self.donor_topology
        # Map "plasmid" to "plasmid_dsDNA" for the lookup
        if topology_key == "plasmid":
            topology_key = "plasmid_dsDNA"
        donor_mult = DONOR_TOPOLOGY_MULTIPLIER.get(topology_key, 1.0)

        # (iii) Stagger enhancement
        #       5' overhangs from staggered cuts enhance HDR.  The
        #       enhancement scales with overhang length.
        #       From constants.py HDR_ENHANCEMENT_PER_BP_OVERHANG.
        stagger_mult = 1.0
        if self.cut_type == "staggered_5prime":
            stagger_mult = 1.0 + HDR_ENHANCEMENT_PER_BP_OVERHANG * self.overhang_length
            # e.g., for 5 bp overhang: 1.0 + 0.15 * 5 = 1.75

        # (iv)  Cell-cycle gating
        #       Only cells in S or G2 phase have the sister chromatid
        #       available as endogenous template and express high levels of
        #       HR factors.  The fraction of cells in S/G2 acts as an
        #       upper bound on HDR-competent cells.
        sg2_fraction = self._cell_params["cell_cycle_s_g2_fraction"]

        # Combined invasion probability (per cell):
        #   P(invasion) = base_hdr * donor_mult * stagger_mult * homology_factor
        #
        # But we also gate by cell cycle: a cell not in S/G2 cannot do HDR.
        # We first draw which cells are in S/G2, then apply the invasion
        # probability only to those cells.

        # Draw cell-cycle phase: True = in S/G2 (HDR-competent)
        in_sg2 = self._rng.random(n) < sg2_fraction  # bool array

        # Compute invasion probability for each cell (can exceed 1.0 from
        # multiplication, so we clip).
        p_invasion = np.clip(
            base_hdr * donor_mult * stagger_mult * homology_overlap_factor,
            0.0,
            1.0,
        )

        # Draw invasion success.  A cell succeeds if it is in S/G2 AND
        # the filament is long enough AND the stochastic invasion roll
        # succeeds.
        invasion_roll = self._rng.random(n) < p_invasion
        invasion_success = filament_long_enough & in_sg2 & invasion_roll

        # ==================================================================
        # STEP 5: DNA synthesis (gene conversion tract length)
        # ==================================================================
        # For cells where strand invasion succeeded, simulate the SDSA
        # synthesis tract.  For failed cells, tract length = 0.
        synthesis_result = self._synthesis_sim.simulate_synthesis(
            donor_topology=self.donor_topology,
            is_staggered=(self.cut_type == "staggered_5prime"),
            n_simulations=n,
        )

        # Mask: only cells with successful invasion get a nonzero tract.
        tract_lengths = np.where(
            invasion_success,
            synthesis_result.tract_lengths_bp,
            0.0,
        )

        # Additionally, the tract length cannot exceed the homology arm
        # length on the donor.  If the polymerase synthesises past the end
        # of the homology arm, it runs into non-homologous donor sequence
        # (e.g., the plasmid backbone or the circularisation junction).
        # For cssDNA donors, the arm defines the boundary of homology.
        # Tract lengths beyond the arm can still result in partial
        # incorporation (the portion within the arm is converted), so we
        # do NOT hard-clip here.  Instead, we note that the effective
        # conversion at a given distance from the cut is limited by the
        # arm length.  For the purpose of this simulation, we allow tract
        # lengths to exceed the arm — the ``probability_at_distance``
        # method will handle the interpretation.

        # ==================================================================
        # HDR success flag
        # ==================================================================
        hdr_success = invasion_success & (tract_lengths > 0)

        # ==================================================================
        # Package results
        # ==================================================================
        self._results = SimulationResults(
            hdr_success=hdr_success,
            tract_lengths_bp=tract_lengths,
            resection_left_bp=resection_result.left_bp,
            resection_right_bp=resection_result.right_bp,
            filament_length_nt=filament_result.filament_length_nt,
            invasion_success=invasion_success,
            n_simulations=n,
            parameters={
                "cut_type": self.cut_type,
                "overhang_length": self.overhang_length,
                "donor_topology": self.donor_topology,
                "homology_arm_length": self.homology_arm_length,
                "cell_type": self.cell_type,
                "n_simulations": n,
                "seed": self.seed,
            },
        )

        return self._results

    def probability_at_distance(self, distance_bp: float) -> float:
        """Probability that the conversion tract extends at least *distance_bp*.

        This answers the question: "If I place a donor-encoded edit
        (e.g., a silent SNP, a reporter gene, a loxP site) at a distance
        of *distance_bp* from the cut site, what fraction of HDR events
        will incorporate that edit?"

        The answer is computed empirically from the Monte Carlo results:
        among all simulations where HDR succeeded, what fraction had a
        tract length ≥ *distance_bp*?

        Parameters
        ----------
        distance_bp : float
            Distance from the cut site to the edit of interest, in bp.

        Returns
        -------
        float
            Fraction of HDR-successful simulations with tract ≥ distance_bp.
            Returns 0.0 if no HDR events occurred.

        Raises
        ------
        RuntimeError
            If ``run()`` has not been called yet.
        """
        if self._results is None:
            raise RuntimeError("Call run() before querying results.")

        # Consider only HDR-successful simulations.
        successful_tracts = self._results.tract_lengths_bp[
            self._results.hdr_success
        ]
        if len(successful_tracts) == 0:
            return 0.0

        # Fraction with tract >= distance_bp
        return float(np.mean(successful_tracts >= distance_bp))

    def plot_tract_distribution(self, ax=None, show: bool = True):
        """Plot a histogram of gene conversion tract lengths.

        Produces an annotated histogram showing:
        - The distribution of tract lengths among HDR-successful cells
        - Vertical lines at the median and 95th percentile
        - An annotation with the HDR success rate

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None
            Axes to plot on.  If None, a new figure is created.
        show : bool
            If True (default), call ``plt.show()``.

        Returns
        -------
        matplotlib.axes.Axes
            The axes object with the plot.

        Raises
        ------
        RuntimeError
            If ``run()`` has not been called yet.
        """
        if self._results is None:
            raise RuntimeError("Call run() before plotting.")

        import matplotlib.pyplot as plt

        successful_tracts = self._results.tract_lengths_bp[
            self._results.hdr_success
        ]

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        if len(successful_tracts) == 0:
            ax.text(
                0.5, 0.5, "No HDR events in simulation",
                ha="center", va="center", transform=ax.transAxes, fontsize=14,
            )
            ax.set_title("Gene Conversion Tract Length Distribution")
            if show:
                plt.show()
            return ax

        # ---- Histogram ----
        # Use 80 bins spanning the range of observed tract lengths.
        # The distribution is typically right-skewed (geometric/exponential).
        bins = np.linspace(0, np.percentile(successful_tracts, 99.5), 80)
        ax.hist(
            successful_tracts,
            bins=bins,
            color="#4C72B0",
            edgecolor="white",
            linewidth=0.5,
            alpha=0.85,
            density=True,
            label="Tract length density",
        )

        # ---- Vertical annotations ----
        median_val = float(np.median(successful_tracts))
        p95_val = float(np.percentile(successful_tracts, 95))
        mean_val = float(np.mean(successful_tracts))

        ax.axvline(
            median_val, color="#DD8452", linestyle="--", linewidth=2,
            label=f"Median = {median_val:.0f} bp",
        )
        ax.axvline(
            p95_val, color="#C44E52", linestyle=":", linewidth=2,
            label=f"95th pctl = {p95_val:.0f} bp",
        )
        ax.axvline(
            mean_val, color="#55A868", linestyle="-.", linewidth=2,
            label=f"Mean = {mean_val:.0f} bp",
        )

        # ---- Homology arm boundary ----
        ax.axvline(
            self.homology_arm_length, color="gray", linestyle="-",
            linewidth=1.5, alpha=0.7,
            label=f"Homology arm = {self.homology_arm_length} bp",
        )

        # ---- Labels and annotations ----
        hdr_rate = float(np.mean(self._results.hdr_success)) * 100
        ax.set_xlabel("Gene Conversion Tract Length (bp)", fontsize=12)
        ax.set_ylabel("Probability Density", fontsize=12)
        ax.set_title(
            f"Gene Conversion Tract Distribution\n"
            f"{self.cut_type} cut | {self.donor_topology} donor | "
            f"{self.cell_type} cells | HDR rate = {hdr_rate:.1f}%",
            fontsize=13,
        )
        ax.legend(fontsize=10, loc="upper right")

        # Text box with summary stats
        textstr = (
            f"n = {len(successful_tracts):,} HDR events\n"
            f"Mean = {mean_val:.0f} bp\n"
            f"Median = {median_val:.0f} bp\n"
            f"P(tract >= 500 bp) = "
            f"{np.mean(successful_tracts >= 500) * 100:.1f}%"
        )
        props = dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.7)
        ax.text(
            0.97, 0.65, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment="top", horizontalalignment="right",
            bbox=props,
        )

        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        plt.tight_layout()

        if show:
            plt.show()

        return ax

    def summary(self) -> Dict[str, float]:
        """Print and return key statistics from the simulation.

        This method provides a comprehensive summary including:
        - HDR success rate (fraction of all simulations)
        - Mean, median, and 95th percentile tract lengths
        - Probability of donor incorporation at key distances from the cut
        - Resection and filament statistics

        Returns
        -------
        dict
            Dictionary of summary statistics.

        Raises
        ------
        RuntimeError
            If ``run()`` has not been called yet.
        """
        if self._results is None:
            raise RuntimeError("Call run() before requesting summary.")

        res = self._results
        successful_tracts = res.tract_lengths_bp[res.hdr_success]

        # ---- Compute statistics ----
        stats: Dict[str, float] = {}

        stats["n_simulations"] = res.n_simulations
        stats["n_hdr_success"] = int(np.sum(res.hdr_success))
        stats["hdr_success_rate"] = float(np.mean(res.hdr_success))

        if len(successful_tracts) > 0:
            stats["tract_mean_bp"] = float(np.mean(successful_tracts))
            stats["tract_median_bp"] = float(np.median(successful_tracts))
            stats["tract_std_bp"] = float(np.std(successful_tracts))
            stats["tract_p5_bp"] = float(np.percentile(successful_tracts, 5))
            stats["tract_p95_bp"] = float(np.percentile(successful_tracts, 95))

            # Probability at key distances
            for dist in [100, 200, 300, 500, 800, 1000, 1500, 2000]:
                key = f"p_conversion_at_{dist}bp"
                stats[key] = float(np.mean(successful_tracts >= dist))
        else:
            stats["tract_mean_bp"] = 0.0
            stats["tract_median_bp"] = 0.0
            stats["tract_std_bp"] = 0.0
            stats["tract_p5_bp"] = 0.0
            stats["tract_p95_bp"] = 0.0

        # Resection statistics
        stats["resection_left_mean_bp"] = float(np.mean(res.resection_left_bp))
        stats["resection_right_mean_bp"] = float(np.mean(res.resection_right_bp))
        stats["filament_mean_nt"] = float(np.mean(res.filament_length_nt))
        stats["invasion_success_rate"] = float(np.mean(res.invasion_success))

        # ---- Print human-readable summary ----
        print("=" * 65)
        print("  ConversionSim — Monte Carlo HDR Simulation Summary")
        print("=" * 65)
        print(f"  Configuration:")
        print(f"    Cut type:          {self.cut_type}")
        if self.cut_type == "staggered_5prime":
            print(f"    Overhang:          {self.overhang_length} bp (5' overhang)")
        print(f"    Donor topology:    {self.donor_topology}")
        print(f"    Homology arms:     {self.homology_arm_length} bp each")
        print(f"    Cell type:         {self.cell_type}")
        print(f"    Simulations:       {res.n_simulations:,}")
        print("-" * 65)
        print(f"  HDR Outcomes:")
        print(f"    HDR success rate:  {stats['hdr_success_rate'] * 100:.1f}%"
              f"  ({stats['n_hdr_success']:,} / {res.n_simulations:,})")
        print(f"    Invasion success:  {stats['invasion_success_rate'] * 100:.1f}%")
        print("-" * 65)

        if len(successful_tracts) > 0:
            print(f"  Tract Length Statistics (among HDR-successful cells):")
            print(f"    Mean:              {stats['tract_mean_bp']:.0f} bp")
            print(f"    Median:            {stats['tract_median_bp']:.0f} bp")
            print(f"    Std dev:           {stats['tract_std_bp']:.0f} bp")
            print(f"    5th percentile:    {stats['tract_p5_bp']:.0f} bp")
            print(f"    95th percentile:   {stats['tract_p95_bp']:.0f} bp")
            print("-" * 65)
            print(f"  Conversion Probability at Distance from Cut:")
            for dist in [100, 200, 300, 500, 800, 1000, 1500, 2000]:
                key = f"p_conversion_at_{dist}bp"
                if key in stats:
                    bar_len = int(stats[key] * 30)
                    bar = "#" * bar_len + "." * (30 - bar_len)
                    print(f"    {dist:>5d} bp:  [{bar}]  {stats[key] * 100:5.1f}%")
        else:
            print("  No HDR events occurred in this simulation.")

        print("-" * 65)
        print(f"  Resection (mean):    L={stats['resection_left_mean_bp']:.0f} bp"
              f"  R={stats['resection_right_mean_bp']:.0f} bp")
        print(f"  Filament (mean):     {stats['filament_mean_nt']:.0f} nt")
        print("=" * 65)

        return stats
