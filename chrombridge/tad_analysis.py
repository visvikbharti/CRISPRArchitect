"""
TAD (Topologically Associating Domain) Boundary Analysis
========================================================

What are TADs?
--------------
If you look at a Hi-C contact map (a matrix showing how often each pair of
genomic loci physically touch each other), you will see a striking pattern:
the map is not uniform. Instead, it is organized into **blocks** along the
diagonal — square-ish patches where contacts are frequent WITHIN the block
but rare BETWEEN adjacent blocks. These blocks are called **Topologically
Associating Domains (TADs)**.

Think of TADs like neighborhoods in a city. People within one neighborhood
interact frequently (they share the same grocery store, park, and school),
but interactions with people in the next neighborhood over are much less
common, even though they might be geographically close. The "boundary" between
neighborhoods is like a major road that people rarely cross for daily
errands.

Key properties of TADs:

    - **Size**: Typically 200 kb to 2 Mb, with a mean of ~800 kb in human
      cells. (Source: Dixon et al., Nature, 2012; Rao et al., Cell, 2014.)

    - **Boundaries**: TAD boundaries are enriched for CTCF binding sites
      (in convergent orientation) and cohesin. CTCF acts like a "stop sign"
      for the loop extrusion machinery (cohesin). When cohesin slides along
      the chromatin fiber, it stops at CTCF sites, creating a loop. The
      base of the loop defines the TAD boundary.

    - **Conservation**: TAD boundaries are relatively conserved across cell
      types and even across species, though there is variation (especially
      in stem cells vs. differentiated cells).

    - **Functional significance**: Genes within the same TAD tend to be
      co-regulated. Enhancers primarily activate promoters within their
      own TAD. Disrupting a TAD boundary (e.g., by deleting the CTCF site)
      can cause enhancer "adoption" — an enhancer from one TAD now activates
      genes in the adjacent TAD. This has been linked to limb malformations
      (Lupianez et al., Cell, 2015) and cancer (Hnisz et al., Science, 2016).

Why do TADs matter for CRISPR editing?
--------------------------------------

1. **Homology search is TAD-bounded**: Recent work (Marin-Gonzalez et al.,
   Science, 2025) shows that cohesin-mediated homology search during HDR
   occurs primarily within the TAD containing the DSB. The broken end
   "scans" for homologous sequences by being extruded through the cohesin
   ring, but this extrusion stops at TAD boundaries. This means a donor
   template placed (or homologous to a region) in a different TAD may be
   less efficiently used for HDR.

2. **Translocation risk within TADs**: Two loci within the same TAD have
   much higher contact frequency than two loci in adjacent TADs. This
   means that if you create two DSBs within the same TAD, the risk of
   deletion/inversion is HIGHER than predicted by genomic distance alone.
   The TAD "forces" the broken ends to stay close together.

3. **Translocation risk across TAD boundaries**: Conversely, if two DSBs
   are in different TADs, the contact frequency drops by ~2-5 fold compared
   to same-TAD loci at the same genomic distance. This somewhat reduces
   rearrangement risk.

Limitations of this module
--------------------------
Without actual Hi-C data for the specific cell type and locus, we cannot
determine exact TAD boundaries. Instead, this module provides **estimates**
based on the statistical properties of TADs (mean size ~800 kb). For
precise TAD boundary information, the user should consult:

    - Published Hi-C datasets (4D Nucleome, ENCODE, Rao et al. 2014)
    - Online browsers: 3D Genome Browser (http://3dgenome.fsm.northwestern.edu)
    - HiGlass (https://higlass.io)

This module flags when loci are near estimated TAD boundaries and
recommends consulting Hi-C data for precise analysis.

References
----------
- Dixon et al., Nature, 2012 (discovery of TADs).
- Rao et al., Cell, 2014 (high-resolution Hi-C, sub-TAD structure).
- Lupianez et al., Cell, 2015 (TAD disruption causes disease).
- Hnisz et al., Science, 2016 (TAD boundaries and oncogene activation).
- Fudenberg et al., Cell Rep, 2016 (loop extrusion model for TADs).
- Marin-Gonzalez et al., Science, 2025 (cohesin-mediated homology search).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple


# =============================================================================
# TAD size parameters
# =============================================================================

# Mean TAD size in human cells (base pairs).
# Source: Dixon et al., 2012; Rao et al., 2014.
# The distribution of TAD sizes is approximately log-normal with
# mean ~800 kb and standard deviation ~400 kb.
MEAN_TAD_SIZE_BP = 800_000
TAD_SIZE_STD_BP = 400_000

# Minimum and maximum plausible TAD sizes
MIN_TAD_SIZE_BP = 200_000    # 200 kb — below this, it's a "sub-TAD" or loop
MAX_TAD_SIZE_BP = 2_500_000  # 2.5 Mb — above this, rare but observed

# Contact frequency fold-change at TAD boundaries
# Within a TAD, contact frequency is ~2-5x higher than across a TAD boundary
# at the same genomic distance. (Source: Dixon et al., 2012.)
TAD_BOUNDARY_CONTACT_FOLD_CHANGE = 3.0  # Central estimate


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class TADInfo:
    """Estimated TAD membership for a genomic locus.

    Since we do not have actual Hi-C data, TAD boundaries are estimated
    using the known statistical properties of TADs. The estimated boundaries
    should be treated as approximate — within ~200-400 kb of the true
    boundaries in most cases.

    Attributes
    ----------
    coordinate_bp : int
        The input genomic coordinate.
    chromosome : str
        Chromosome identifier.
    estimated_tad_start : int
        Estimated start coordinate of the containing TAD.
    estimated_tad_end : int
        Estimated end coordinate of the containing TAD.
    estimated_tad_size : int
        Size of the estimated TAD in bp.
    estimated_tad_index : int
        Index of this TAD along the chromosome (0-based), based on
        uniform-size tiling.
    distance_to_nearest_boundary : int
        Distance from the locus to the nearest estimated TAD boundary.
        Loci near boundaries (<100 kb) may be at the edges of real TADs.
    near_boundary : bool
        True if the locus is within 100 kb of an estimated TAD boundary.
    note : str
        Important caveats about the estimation method.
    """

    coordinate_bp: int
    chromosome: str
    estimated_tad_start: int
    estimated_tad_end: int
    estimated_tad_size: int
    estimated_tad_index: int
    distance_to_nearest_boundary: int
    near_boundary: bool
    note: str


# =============================================================================
# TAD Analyzer
# =============================================================================

class TADAnalyzer:
    """Estimate TAD membership and analyze boundary effects.

    This class provides approximate TAD analysis without requiring Hi-C data.
    It uses a simple tiling model: the chromosome is divided into uniform
    TADs of a given size (default: 800 kb). While real TADs vary in size
    and position, this approximation is useful for:

    1. Determining whether two loci are likely in the same or different TADs.
    2. Estimating how many TAD boundaries lie between two loci.
    3. Flagging loci that are near estimated boundaries (where the estimate
       is least reliable).

    For precise TAD boundary information, consult published Hi-C data.

    Analogy
    -------
    Imagine you know that a city has neighborhoods averaging 1 km in diameter,
    but you don't have a map showing exact neighborhood boundaries. If someone
    asks "Are addresses 100 Main St and 200 Main St in the same neighborhood?",
    you can give a pretty good answer if they are only 100 m apart (yes) or
    5 km apart (no). But if they are 800 m apart, you'd say "probably, but
    I'd need to check a map". That's essentially what this module does.

    Usage
    -----
    >>> analyzer = TADAnalyzer()
    >>> tad_info = analyzer.estimate_tad_membership(50_000_000, "chr7")
    >>> print(f"Estimated TAD: {tad_info.estimated_tad_start:,} - "
    ...       f"{tad_info.estimated_tad_end:,}")

    >>> same = analyzer.same_tad(50_000_000, 50_500_000)
    >>> print(f"Same TAD? {same}")  # Likely True (500 kb apart, within one TAD)

    >>> same = analyzer.same_tad(50_000_000, 52_000_000)
    >>> print(f"Same TAD? {same}")  # Likely False (2 Mb apart, ~2-3 TADs)
    """

    def __init__(self, tad_size_bp: int = MEAN_TAD_SIZE_BP) -> None:
        """
        Parameters
        ----------
        tad_size_bp : int
            Assumed uniform TAD size in bp. Default: 800,000 (800 kb).
            This can be adjusted if the user has prior knowledge about
            TAD sizes in their region of interest (e.g., from published
            Hi-C data showing larger or smaller TADs near their locus).
        """
        self.tad_size_bp = tad_size_bp

    def estimate_tad_membership(
        self,
        coordinate_bp: int,
        chromosome: str = "unknown",
    ) -> TADInfo:
        """Estimate which TAD a genomic locus belongs to.

        Uses a simple uniform-tiling model: the chromosome is divided into
        non-overlapping TADs of size ``self.tad_size_bp``, starting from
        position 0.

        Parameters
        ----------
        coordinate_bp : int
            Genomic coordinate of the locus in bp (absolute position on
            the chromosome, e.g., 50000000 for a locus at 50 Mb).
        chromosome : str
            Chromosome identifier (e.g., "chr7", "chrX"). Used for
            labeling only — does not affect the calculation.

        Returns
        -------
        TADInfo
            Estimated TAD boundaries and metadata.

        Notes
        -----
        The note field in the returned TADInfo explicitly states that these
        are estimates. If the locus is near an estimated boundary (within
        100 kb), the ``near_boundary`` flag is set, indicating that the
        assignment is particularly uncertain and Hi-C data should be
        consulted.
        """
        tad_index = coordinate_bp // self.tad_size_bp
        tad_start = tad_index * self.tad_size_bp
        tad_end = tad_start + self.tad_size_bp

        # Distance to nearest boundary
        dist_to_start = coordinate_bp - tad_start
        dist_to_end = tad_end - coordinate_bp
        dist_to_boundary = min(dist_to_start, dist_to_end)
        near_boundary = dist_to_boundary < 100_000  # 100 kb threshold

        note = (
            "TAD boundaries are ESTIMATED using a uniform-tiling model with "
            f"TAD size = {self.tad_size_bp / 1000:.0f} kb. Real TADs vary in "
            "size (200 kb to 2.5 Mb) and their boundaries depend on CTCF/cohesin "
            "binding, which is cell-type-specific. For precise boundaries, "
            "consult Hi-C data for your cell type (e.g., from 4D Nucleome, "
            "ENCODE, or Rao et al., Cell, 2014). "
        )
        if near_boundary:
            note += (
                "WARNING: This locus is within 100 kb of an estimated TAD "
                "boundary. The TAD assignment is particularly uncertain. "
                "Consulting Hi-C data is strongly recommended."
            )

        return TADInfo(
            coordinate_bp=coordinate_bp,
            chromosome=chromosome,
            estimated_tad_start=tad_start,
            estimated_tad_end=tad_end,
            estimated_tad_size=self.tad_size_bp,
            estimated_tad_index=tad_index,
            distance_to_nearest_boundary=dist_to_boundary,
            near_boundary=near_boundary,
            note=note,
        )

    def same_tad(
        self,
        coord1: int,
        coord2: int,
        chromosome: Optional[str] = None,
    ) -> bool:
        """Estimate whether two loci are within the same TAD.

        Two loci are assigned to the same TAD if they fall within the same
        uniform-tiling interval. This is an approximation — the true answer
        depends on actual CTCF/cohesin binding and loop extrusion dynamics.

        Parameters
        ----------
        coord1 : int
            Genomic coordinate of the first locus (bp).
        coord2 : int
            Genomic coordinate of the second locus (bp).
        chromosome : str, optional
            Chromosome (for annotation only).

        Returns
        -------
        bool
            True if both loci are estimated to be in the same TAD.

        Notes
        -----
        If the two loci are in different TADs, the contact frequency between
        them is expected to be ~2-5x lower than if they were in the same TAD
        at the same genomic distance. This affects both homology search
        efficiency (Marin-Gonzalez et al., Science, 2025) and translocation
        risk.
        """
        tad_index1 = coord1 // self.tad_size_bp
        tad_index2 = coord2 // self.tad_size_bp
        return tad_index1 == tad_index2

    def inter_tad_distance(
        self,
        coord1: int,
        coord2: int,
    ) -> int:
        """Count the number of TAD boundaries crossed between two loci.

        Parameters
        ----------
        coord1 : int
            Genomic coordinate of the first locus (bp).
        coord2 : int
            Genomic coordinate of the second locus (bp).

        Returns
        -------
        int
            Number of estimated TAD boundaries between the two loci.
            0 means they are in the same TAD.
            1 means they are in adjacent TADs (one boundary crossed).
            etc.

        Notes
        -----
        Each TAD boundary crossed reduces the expected contact frequency
        by a factor of ~2-5x (relative to same-TAD contacts at that genomic
        distance). Multiple boundary crossings compound this reduction.

        This has practical implications:
            - Homology search during HDR is less efficient across TAD
              boundaries (Marin-Gonzalez et al., Science, 2025).
            - Translocation risk is reduced when breaks are in different
              TADs (lower contact frequency = lower probability of end-joining).
        """
        tad_index1 = coord1 // self.tad_size_bp
        tad_index2 = coord2 // self.tad_size_bp
        return abs(tad_index2 - tad_index1)

    def contact_frequency_modifier(
        self,
        coord1: int,
        coord2: int,
    ) -> float:
        """Estimate how TAD boundaries modify contact frequency.

        Returns a multiplicative factor that adjusts the "bare" contact
        frequency (from polymer models) to account for TAD boundaries.

        - Same TAD: factor = 1.0 (no modification; polymer model is ~correct).
        - Across 1 boundary: factor ~ 1/3 (contact frequency reduced).
        - Across N boundaries: factor ~ (1/3)^N.

        This is a rough estimate. Real TAD boundaries vary in "strength"
        (how much they reduce contact frequency), and some boundaries are
        leaky while others are strong insulators.

        Parameters
        ----------
        coord1 : int
            First locus coordinate (bp).
        coord2 : int
            Second locus coordinate (bp).

        Returns
        -------
        float
            Multiplicative modifier for contact frequency.
            1.0 = same TAD; < 1.0 = different TADs.
        """
        n_boundaries = self.inter_tad_distance(coord1, coord2)
        if n_boundaries == 0:
            return 1.0
        return (1.0 / TAD_BOUNDARY_CONTACT_FOLD_CHANGE) ** n_boundaries

    def tad_aware_analysis(
        self,
        coord1: int,
        coord2: int,
        chromosome: str = "unknown",
    ) -> str:
        """Generate a formatted TAD-aware analysis of two loci.

        This produces a human-readable report about the TAD context of
        two genomic loci, with implications for CRISPR editing.

        Parameters
        ----------
        coord1 : int
            First locus coordinate (bp).
        coord2 : int
            Second locus coordinate (bp).
        chromosome : str
            Chromosome identifier.

        Returns
        -------
        str
            Formatted analysis report.
        """
        tad1 = self.estimate_tad_membership(coord1, chromosome)
        tad2 = self.estimate_tad_membership(coord2, chromosome)
        same = self.same_tad(coord1, coord2)
        n_boundaries = self.inter_tad_distance(coord1, coord2)
        contact_mod = self.contact_frequency_modifier(coord1, coord2)
        genomic_dist = abs(coord2 - coord1)

        lines = []
        lines.append("TAD-AWARE LOCUS ANALYSIS")
        lines.append("=" * 60)
        lines.append("")

        # Format coordinates
        def _fmt(bp: int) -> str:
            if bp >= 1_000_000:
                return f"{bp / 1_000_000:.2f} Mb"
            elif bp >= 1_000:
                return f"{bp / 1_000:.1f} kb"
            return f"{bp} bp"

        lines.append(f"Chromosome: {chromosome}")
        lines.append(f"Locus 1: {coord1:,} bp ({_fmt(coord1)})")
        lines.append(f"Locus 2: {coord2:,} bp ({_fmt(coord2)})")
        lines.append(f"Genomic distance: {_fmt(genomic_dist)}")
        lines.append("")

        lines.append("--- Estimated TAD Membership ---")
        lines.append(
            f"  Locus 1: TAD #{tad1.estimated_tad_index} "
            f"({_fmt(tad1.estimated_tad_start)} - {_fmt(tad1.estimated_tad_end)})"
        )
        if tad1.near_boundary:
            lines.append(
                f"    ** Near boundary ({_fmt(tad1.distance_to_nearest_boundary)} "
                f"from edge) — assignment uncertain **"
            )

        lines.append(
            f"  Locus 2: TAD #{tad2.estimated_tad_index} "
            f"({_fmt(tad2.estimated_tad_start)} - {_fmt(tad2.estimated_tad_end)})"
        )
        if tad2.near_boundary:
            lines.append(
                f"    ** Near boundary ({_fmt(tad2.distance_to_nearest_boundary)} "
                f"from edge) — assignment uncertain **"
            )

        lines.append("")
        lines.append("--- TAD Relationship ---")
        if same:
            lines.append(
                "  Both loci are in the SAME estimated TAD."
            )
            lines.append(
                "  Implications:"
            )
            lines.append(
                "    - High contact frequency between the loci."
            )
            lines.append(
                "    - Cohesin-mediated homology search can span both loci"
            )
            lines.append(
                "      (Marin-Gonzalez et al., Science, 2025)."
            )
            lines.append(
                "    - If two DSBs are made: ELEVATED deletion/inversion risk"
            )
            lines.append(
                "      due to high intra-TAD contact frequency."
            )
        else:
            lines.append(
                f"  Loci are in DIFFERENT TADs, separated by {n_boundaries} "
                f"estimated TAD boundary/boundaries."
            )
            lines.append(
                f"  Contact frequency modifier: {contact_mod:.4f} "
                f"({contact_mod * 100:.2f}% of same-TAD frequency)"
            )
            lines.append(
                "  Implications:"
            )
            lines.append(
                "    - Lower contact frequency between the loci than expected"
            )
            lines.append(
                "      from genomic distance alone."
            )
            lines.append(
                "    - Cohesin-mediated homology search is LESS efficient across"
            )
            lines.append(
                "      TAD boundaries (extrusion stops at CTCF sites)."
            )
            lines.append(
                "    - If two DSBs are made: rearrangement risk is LOWER than"
            )
            lines.append(
                "      same-TAD breaks at the same genomic distance."
            )
            lines.append(
                "    - However, the absolute risk still depends on genomic"
            )
            lines.append(
                "      distance (see translocation module)."
            )

        lines.append("")
        lines.append("--- Caveats ---")
        lines.append(
            f"  TAD boundaries are estimated using a {self.tad_size_bp / 1000:.0f} kb"
            " uniform-tiling model."
        )
        lines.append(
            "  Real TADs vary in size (200 kb - 2.5 Mb) and boundary positions"
        )
        lines.append(
            "  depend on CTCF/cohesin binding patterns specific to each cell type."
        )
        lines.append(
            "  For definitive TAD boundary information, consult Hi-C data:"
        )
        lines.append(
            "    - 3D Genome Browser: http://3dgenome.fsm.northwestern.edu"
        )
        lines.append(
            "    - 4D Nucleome: https://data.4dnucleome.org"
        )
        lines.append(
            "    - ENCODE: https://www.encodeproject.org"
        )

        return "\n".join(lines)
