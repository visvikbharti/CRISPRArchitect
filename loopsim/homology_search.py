"""
Homology Search Predictor — Mapping Loop Extrusion to HDR Donor Accessibility
===============================================================================

This module translates the raw output of cohesin loop extrusion simulations
into biologically meaningful predictions about the **homology search** step
of homology-directed repair (HDR).

The Central Question
--------------------
After a CRISPR-induced DSB, the cell must find a homologous template to
perform HDR. The homology search is mediated by the **RAD51 filament** —
a helical nucleoprotein filament formed on the 3' single-stranded DNA
(ssDNA) overhang created by end resection at the break.

The RAD51 filament searches for complementary sequences by:
1. **Sampling**: The filament transiently contacts duplex DNA and tests
   for sequence complementarity using an 8-nucleotide "microhomology
   sampling" mechanism (Qi Z, et al., Cell, 2015).
2. **Stable pairing**: If a region of sufficient homology is found
   (>=15-20 bp of continuous match), the filament catalyzes strand
   invasion, forming a displacement loop (D-loop).

But the nucleus is ENORMOUS relative to the search target. The diploid
human genome is ~6.4 billion bp, and the target homology might be only
300-1000 bp. How does RAD51 search efficiently?

The answer, revealed by Marin-Gonzalez et al. (2025), is that **cohesin
loop extrusion dramatically narrows the search space**:

- Cohesin loaded at the DSB reels chromatin through the break site.
- As chromatin passes through, the DSB-anchored RAD51 filament can scan
  the passing duplex DNA for complementarity.
- This converts a 3D search problem (finding a needle in a nuclear haystack)
  into a 1D scanning problem (reading a tape as it passes through a head).

Two Categories of Donor Template
---------------------------------
For genome editing, we must distinguish:

1. **Endogenous / cis donors** (on the same chromosome):
   - Sister chromatid (available only in S/G2 phase)
   - Ectopic homology elsewhere on the same chromosome
   - For these donors, loop extrusion is CRITICAL: the donor must be
     within the extruded search domain to be found efficiently.
   - If the donor is outside the domain (beyond the nearest CTCF barrier),
     the probability of homologous pairing drops dramatically.

2. **Exogenous donors** (delivered to the cell):
   - ssODN, lssDNA, cssDNA, plasmid dsDNA, AAV
   - These donors are FREE IN THE NUCLEOPLASM — they are not part of the
     chromatin fiber and are not subject to loop extrusion.
   - The exogenous donor finds the DSB by **3D diffusion** through the
     nuclear volume. Once it arrives at the DSB, it can pair with the
     resected ssDNA directly, without requiring cohesin-mediated scanning.
   - Therefore, loop extrusion is NOT the rate-limiting step for exogenous
     donor-mediated HDR. The bottleneck is: (a) donor delivery/stability,
     (b) donor arrival by diffusion, (c) competition with NHEJ.

   This is a crucial distinction for the CRISPRArchitect toolkit: when
   designing a cssDNA donor for knock-in, the cohesin search domain is
   not directly relevant — the cssDNA donor will find the DSB regardless
   of CTCF landscape. However, loop extrusion IS relevant for:
   - Understanding why paired-DSB strategies might fail (if two DSBs are
     far apart, their search domains do not overlap, so one DSB cannot
     "help" the other find homology)
   - Predicting sister chromatid recombination rates (important for LOH
     risk assessment)

References
----------
- Marin-Gonzalez A, et al. "Cohesin drives chromatin scanning during
  the RAD51-mediated homology search." Science, 2025.
- Qi Z, et al. "DNA sequence alignment by microhomology sampling during
  homologous recombination." Cell, 2015.
- Arnould C, et al. "Loop extrusion as a mechanism for formation of DNA
  damage repair foci." Nature, 2021.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from .chromatin_fiber import ChromatinFiber
from .cohesin_extruder import CohesinExtruder, SimulationResults


# ---------------------------------------------------------------------------
# Result data class
# ---------------------------------------------------------------------------

@dataclass
class SearchDomainResult:
    """Predicted homology search domain from a Monte Carlo extrusion ensemble.

    This is the biologically interpretable output: for a given DSB on a given
    chromatin fiber, what region of the genome is scanned by RAD51 during
    cohesin-mediated loop extrusion?

    The ``position_scan_probability`` array is the key output. It is a
    per-bead probability (values 0-1) indicating the fraction of Monte Carlo
    simulations in which that bead fell within the extruded search domain.

    Interpretation:
    - Beads near the DSB will have probability ~1.0 (always scanned).
    - Beads at or beyond strong CTCF barriers will have low probability
      (only scanned when CTCF is bypassed).
    - The probability gradient at domain edges reflects CTCF barrier
      permeability across the cell population.

    Attributes
    ----------
    mean_domain_size_bp : float
        Average search domain size across simulations (bp).
    median_domain_size_bp : float
        Median search domain size (bp). More robust to outliers than mean.
    std_domain_size_bp : float
        Standard deviation of domain sizes (bp).
    left_boundaries : np.ndarray
        Left boundary position (bp) from each simulation. Shape: (N,).
    right_boundaries : np.ndarray
        Right boundary position (bp) from each simulation. Shape: (N,).
    position_scan_probability : np.ndarray
        Per-bead probability of being within the search domain. Shape:
        (n_beads,). This is the primary output for visualization.
    domain_sizes : np.ndarray
        Domain size (bp) from each simulation. Shape: (N,).
    dsb_position_bp : int
        The DSB position used for this prediction.
    n_simulations : int
        Number of Monte Carlo simulations in the ensemble.
    simulation_results : SimulationResults
        The raw simulation results (for advanced analysis).
    """

    mean_domain_size_bp: float
    median_domain_size_bp: float
    std_domain_size_bp: float
    left_boundaries: np.ndarray
    right_boundaries: np.ndarray
    position_scan_probability: np.ndarray
    domain_sizes: np.ndarray
    dsb_position_bp: int
    n_simulations: int
    simulation_results: SimulationResults = field(repr=False)


# ---------------------------------------------------------------------------
# HomologySearchPredictor
# ---------------------------------------------------------------------------

class HomologySearchPredictor:
    """Predicts the RAD51 homology search domain based on cohesin loop extrusion.

    This class wraps the CohesinExtruder and translates its raw output into
    predictions about which genomic regions are accessible to the RAD51
    filament. It answers questions like:

    - "If I make a DSB at position X, what region of the chromosome will
      be scanned for homology?"
    - "If I have a cis donor at position Y, what is the probability that
      loop extrusion will bring it into the search domain?"
    - "How do exogenous (cssDNA) and endogenous (sister chromatid) donors
      differ in their reliance on loop extrusion?"

    The predictor uses a Monte Carlo ensemble of extrusion simulations to
    capture cell-to-cell variability. Each simulation represents one cell
    in the population, with stochastic CTCF occupancy and heterochromatin
    effects.

    Parameters
    ----------
    extruder : CohesinExtruder, optional
        A configured CohesinExtruder instance. If None, a default extruder
        is created with standard parameters (1 kb/s extrusion rate, 0.8
        CTCF stall probability, 3600 s max time).

    Examples
    --------
    >>> fiber = ChromatinFiber(0, 5_000_000)
    >>> fiber.add_ctcf_sites_uniform()
    >>> predictor = HomologySearchPredictor()
    >>> result = predictor.predict_search_domain(fiber, 2_500_000)
    >>> print(f"Mean domain: {result.mean_domain_size_bp / 1e6:.2f} Mb")
    """

    def __init__(self, extruder: Optional[CohesinExtruder] = None) -> None:
        if extruder is None:
            self.extruder = CohesinExtruder()
        else:
            self.extruder = extruder

    # -----------------------------------------------------------------
    # Core prediction
    # -----------------------------------------------------------------

    def predict_search_domain(
        self,
        fiber: ChromatinFiber,
        dsb_position_bp: int,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> SearchDomainResult:
        """Predict the homology search domain for a DSB on a chromatin fiber.

        Runs a Monte Carlo ensemble of cohesin loop extrusion simulations and
        computes the per-bead probability of being scanned by the RAD51
        filament. This is the main entry point for search domain prediction.

        The output ``position_scan_probability`` array can be plotted as a
        genomic track to visualize the search domain, similar to a ChIP-seq
        enrichment profile but for "extrusion coverage."

        Parameters
        ----------
        fiber : ChromatinFiber
            The chromatin fiber model (substrate for extrusion).
        dsb_position_bp : int
            Genomic coordinate of the DSB where cohesin will be loaded.
        n_simulations : int, optional
            Number of Monte Carlo simulations. Default 1000.
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        SearchDomainResult
            Comprehensive prediction including per-bead scan probabilities,
            domain size distribution, and boundary distributions.

        Notes
        -----
        The scan probability at position X is computed as:

            P(X scanned) = (number of simulations where left_boundary <= X
                            AND X <= right_boundary) / n_simulations

        Positions very close to the DSB will have P ~ 1.0. Positions beyond
        the typical CTCF boundary will have P << 1.0 (only scanned when
        CTCF is bypassed).
        """
        # Run the ensemble
        sim_results = self.extruder.simulate(
            fiber, dsb_position_bp, n_simulations=n_simulations, seed=seed
        )

        # Compute per-bead scan probability
        # For each simulation, the domain spans [left_bead, right_bead].
        # Count how many simulations include each bead.
        scan_counts = np.zeros(fiber.n_beads, dtype=np.float64)

        for result in sim_results.individual_results:
            left_bead = fiber.get_bead_index(result.left_boundary_bp)
            right_bead = fiber.get_bead_index(result.right_boundary_bp)
            scan_counts[left_bead : right_bead + 1] += 1.0

        scan_probability = scan_counts / n_simulations

        return SearchDomainResult(
            mean_domain_size_bp=sim_results.mean_domain_size_bp,
            median_domain_size_bp=sim_results.median_domain_size_bp,
            std_domain_size_bp=sim_results.std_domain_size_bp,
            left_boundaries=sim_results.left_boundaries_bp,
            right_boundaries=sim_results.right_boundaries_bp,
            position_scan_probability=scan_probability,
            domain_sizes=sim_results.domain_sizes_bp,
            dsb_position_bp=dsb_position_bp,
            n_simulations=n_simulations,
            simulation_results=sim_results,
        )

    # -----------------------------------------------------------------
    # Cis donor accessibility
    # -----------------------------------------------------------------

    def probability_donor_in_domain(
        self,
        dsb_position_bp: int,
        donor_position_bp: int,
        fiber: ChromatinFiber,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> float:
        """Compute the probability that a cis donor falls within the search domain.

        For a donor template located on the same chromosome as the DSB (e.g.,
        sister chromatid locus, ectopic homology), this method asks: "In what
        fraction of cells does the cohesin-extruded domain extend far enough
        to include the donor position?"

        This is only meaningful for **endogenous / cis** donors. For exogenous
        donors (cssDNA, ssODN, AAV), the donor finds the DSB by 3D diffusion,
        not by loop extrusion — so this probability is irrelevant.

        Parameters
        ----------
        dsb_position_bp : int
            Genomic coordinate of the DSB.
        donor_position_bp : int
            Genomic coordinate of the cis donor template.
        fiber : ChromatinFiber
            The chromatin fiber model.
        n_simulations : int, optional
            Number of Monte Carlo simulations. Default 1000.
        seed : int, optional
            Random seed.

        Returns
        -------
        float
            Fraction of simulations (0.0 to 1.0) in which the donor position
            falls within the extruded search domain.

        Notes
        -----
        Expected values:
        - Donor very close to DSB (< 100 kb): P ~ 1.0 (almost always scanned)
        - Donor within same TAD (100-500 kb): P ~ 0.3-0.9 (depends on CTCF)
        - Donor in adjacent TAD (500 kb - 2 Mb): P ~ 0.01-0.2 (low, needs
          CTCF bypass)
        - Donor very far (> 2 Mb): P ~ 0.0 (essentially never reached)

        Examples
        --------
        >>> # Sister chromatid homology is at the same position -> always in domain
        >>> p = predictor.probability_donor_in_domain(2_500_000, 2_500_000, fiber)
        >>> # p ~ 1.0

        >>> # Ectopic donor 500 kb away
        >>> p = predictor.probability_donor_in_domain(2_500_000, 3_000_000, fiber)
        >>> # p depends on CTCF landscape
        """
        sim_results = self.extruder.simulate(
            fiber, dsb_position_bp, n_simulations=n_simulations, seed=seed
        )

        count_in_domain = 0
        for result in sim_results.individual_results:
            if result.left_boundary_bp <= donor_position_bp <= result.right_boundary_bp:
                count_in_domain += 1

        return count_in_domain / n_simulations

    # -----------------------------------------------------------------
    # Exogenous vs. endogenous comparison
    # -----------------------------------------------------------------

    def compare_exogenous_vs_endogenous(
        self,
        dsb_position_bp: int,
        fiber: ChromatinFiber,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> Dict[str, object]:
        """Compare the role of loop extrusion for exogenous vs. endogenous donors.

        This method generates a structured comparison explaining how the
        cohesin search domain affects HDR differently depending on donor type.
        It is meant to be educational — the output includes explanatory text
        alongside quantitative predictions.

        Key biological insight:

        **Exogenous donor (e.g., cssDNA):**
        The donor is a short (~1-3 kb) single-stranded DNA molecule
        delivered by electroporation. It exists free in the nucleoplasm,
        not incorporated into chromatin. It finds the DSB by 3D diffusion
        through the nuclear volume. Loop extrusion does NOT help or hinder
        this process — the cssDNA simply diffuses until it encounters the
        DSB-associated RAD51 filament. The rate-limiting steps are:
        (1) donor delivery and nuclear entry, (2) donor stability (cssDNA
        resists exonuclease degradation, giving it a longer half-life than
        linear ssDNA), (3) competition with NHEJ for the DSB.

        **Endogenous donor (sister chromatid):**
        The sister chromatid is the identical copy of the broken chromosome,
        available in S/G2 phase. It is part of the chromatin fiber and can
        only be found by the RAD51 filament if loop extrusion brings it
        into proximity with the DSB. Since the sister chromatid homology is
        at the exact same genomic position, it is ALWAYS within the search
        domain (probability = 1.0). However, ectopic homologies (e.g.,
        duplicated genes) may or may not be within the domain.

        Parameters
        ----------
        dsb_position_bp : int
            Genomic coordinate of the DSB.
        fiber : ChromatinFiber
            The chromatin fiber model.
        n_simulations : int, optional
            Number of Monte Carlo simulations. Default 1000.
        seed : int, optional
            Random seed.

        Returns
        -------
        dict
            Dictionary with keys:
            - "search_domain": SearchDomainResult
            - "exogenous_donor_note": str (explanation)
            - "endogenous_donor_note": str (explanation)
            - "sister_chromatid_in_domain_prob": float
            - "mean_domain_size_mb": float
            - "recommendation": str
        """
        domain_result = self.predict_search_domain(
            fiber, dsb_position_bp, n_simulations=n_simulations, seed=seed
        )

        mean_mb = domain_result.mean_domain_size_bp / 1e6
        median_mb = domain_result.median_domain_size_bp / 1e6

        return {
            "search_domain": domain_result,

            "exogenous_donor_note": (
                "EXOGENOUS DONORS (cssDNA, ssODN, plasmid, AAV):\n"
                "Loop extrusion is NOT relevant for exogenous donor-mediated "
                "HDR. The donor molecule is free in the nucleoplasm and reaches "
                "the DSB by 3D diffusion, independent of the cohesin-extruded "
                "chromatin loop. HDR efficiency with exogenous donors depends on: "
                "(1) donor intracellular concentration, (2) donor stability "
                "(cssDNA > lssDNA due to exonuclease resistance), (3) homology "
                "arm design, and (4) competition with NHEJ. The cohesin search "
                "domain does not constrain exogenous donor access."
            ),

            "endogenous_donor_note": (
                f"ENDOGENOUS DONORS (sister chromatid, ectopic homology):\n"
                f"Loop extrusion creates a search domain of {mean_mb:.2f} Mb "
                f"(mean) / {median_mb:.2f} Mb (median) around the DSB. "
                f"The sister chromatid (at the identical genomic position) is "
                f"ALWAYS within this domain. However, ectopic homologies must "
                f"fall within the domain boundaries (set by CTCF barriers) to "
                f"be efficiently found. Sequences beyond the nearest strong "
                f"CTCF sites have low scan probability."
            ),

            "sister_chromatid_in_domain_prob": 1.0,  # Always at DSB position

            "mean_domain_size_mb": mean_mb,

            "recommendation": (
                "For cssDNA-based knock-in experiments, do NOT rely on loop "
                "extrusion for donor access. Focus on optimizing donor design "
                "(homology arm length, strand selection, circularization) and "
                "delivery (electroporation parameters, donor concentration). "
                "Loop extrusion is relevant only when assessing risk of "
                "unintended sister chromatid recombination (loss of "
                "heterozygosity) or when evaluating dual-DSB strategies where "
                "one break's search domain might (or might not) overlap the other."
            ),
        }

    # -----------------------------------------------------------------
    # Dual-DSB overlap analysis
    # -----------------------------------------------------------------

    def dual_dsb_overlap_probability(
        self,
        dsb1_position_bp: int,
        dsb2_position_bp: int,
        fiber: ChromatinFiber,
        n_simulations: int = 1000,
        seed: Optional[int] = None,
    ) -> Dict[str, object]:
        """Analyze whether search domains from two DSBs overlap.

        For dual-cut CRISPR strategies (e.g., paired Cas9 nickases, or two
        separate guides cutting ~1 Mb apart), this method assesses whether
        the cohesin-extruded search domain from one DSB can reach the other.

        **Key prediction:** For DSBs separated by more than the typical TAD
        size (~500 kb - 2 Mb), the search domains will almost NEVER overlap.
        This means:
        - Each DSB has its own independent repair dynamics.
        - The segment between two distant DSBs is at high risk of being
          lost (deletion) or rearranged (inversion, translocation).
        - Exogenous donors at each DSB site must be provided independently.

        This is a critical finding from the Marin-Gonzalez et al. (2025) work:
        the cohesin search domain is bounded by CTCF sites and typically spans
        only one TAD (~500 kb - 2 Mb). Two DSBs in different TADs will have
        non-overlapping search domains.

        Parameters
        ----------
        dsb1_position_bp : int
            Genomic coordinate of the first DSB.
        dsb2_position_bp : int
            Genomic coordinate of the second DSB.
        fiber : ChromatinFiber
            The chromatin fiber model.
        n_simulations : int, optional
            Number of Monte Carlo simulations. Default 1000.
        seed : int, optional
            Random seed.

        Returns
        -------
        dict
            Dictionary with keys:
            - "dsb1_domain": SearchDomainResult
            - "dsb2_domain": SearchDomainResult
            - "distance_bp": int (distance between the two DSBs)
            - "dsb1_reaches_dsb2_prob": float (fraction of DSB1 domains that
              include DSB2 position)
            - "dsb2_reaches_dsb1_prob": float (vice versa)
            - "any_overlap_prob": float (fraction of simulations where the
              domains overlap at all)
            - "interpretation": str (biological interpretation)
        """
        seed1 = seed if seed is not None else None
        seed2 = (seed + n_simulations) if seed is not None else None

        domain1 = self.predict_search_domain(
            fiber, dsb1_position_bp, n_simulations=n_simulations, seed=seed1
        )
        domain2 = self.predict_search_domain(
            fiber, dsb2_position_bp, n_simulations=n_simulations, seed=seed2
        )

        distance = abs(dsb2_position_bp - dsb1_position_bp)

        # Check: does DSB1's domain include DSB2's position?
        dsb2_bead = fiber.get_bead_index(dsb2_position_bp)
        dsb1_reaches_dsb2 = float(domain1.position_scan_probability[dsb2_bead])

        dsb1_bead = fiber.get_bead_index(dsb1_position_bp)
        dsb2_reaches_dsb1 = float(domain2.position_scan_probability[dsb1_bead])

        # Check for any overlap between the two domains (pairwise)
        # For each pair of simulations (one from each DSB), check if
        # the domains overlap. Since simulations are independent, we check
        # whether the range from DSB1's right boundary can reach past DSB2's
        # left boundary (or vice versa), for PAIRED simulations.
        n = min(len(domain1.left_boundaries), len(domain2.left_boundaries))
        overlap_count = 0
        for i in range(n):
            # Domains overlap if one's right >= other's left AND vice versa
            d1_left = domain1.left_boundaries[i]
            d1_right = domain1.right_boundaries[i]
            d2_left = domain2.left_boundaries[i]
            d2_right = domain2.right_boundaries[i]
            if d1_right >= d2_left and d2_right >= d1_left:
                overlap_count += 1
        any_overlap_prob = overlap_count / n if n > 0 else 0.0

        # Biological interpretation
        dist_mb = distance / 1e6
        if any_overlap_prob > 0.5:
            interp = (
                f"The two DSBs are {dist_mb:.2f} Mb apart. Their search domains "
                f"overlap in {any_overlap_prob:.0%} of simulations, indicating "
                f"they are likely within the same TAD or CTCF-bounded domain. "
                f"Cohesin from one DSB can potentially scan past the other DSB "
                f"position, which may facilitate coordinated repair."
            )
        elif any_overlap_prob > 0.05:
            interp = (
                f"The two DSBs are {dist_mb:.2f} Mb apart. Their search domains "
                f"overlap in only {any_overlap_prob:.0%} of simulations — they are "
                f"separated by CTCF barriers that usually prevent the extrusion "
                f"domain of one DSB from reaching the other. These DSBs will "
                f"predominantly undergo INDEPENDENT repair, increasing the risk "
                f"of intervening segment deletion or rearrangement."
            )
        else:
            interp = (
                f"The two DSBs are {dist_mb:.2f} Mb apart. Their search domains "
                f"essentially NEVER overlap ({any_overlap_prob:.1%} of simulations). "
                f"Multiple CTCF barriers separate these loci, placing them in "
                f"different TADs. Each DSB will be repaired independently. "
                f"If both DSBs are on the same chromosome, the intervening "
                f"segment (~{dist_mb:.1f} Mb) has a high probability of being "
                f"deleted. This confirms that widely-spaced dual-DSB strategies "
                f"carry significant risk of large deletions and should be avoided "
                f"unless segment excision is the intended outcome."
            )

        return {
            "dsb1_domain": domain1,
            "dsb2_domain": domain2,
            "distance_bp": distance,
            "dsb1_reaches_dsb2_prob": dsb1_reaches_dsb2,
            "dsb2_reaches_dsb1_prob": dsb2_reaches_dsb1,
            "any_overlap_prob": any_overlap_prob,
            "interpretation": interp,
        }
