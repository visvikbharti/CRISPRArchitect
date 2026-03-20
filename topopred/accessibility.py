"""
Homology Arm Accessibility Scorer for cssDNA Donor Templates
=============================================================

BIOLOGICAL BACKGROUND — WHY ACCESSIBILITY MATTERS FOR HDR
----------------------------------------------------------

Homology-directed repair (HDR) using a cssDNA donor template follows these
key molecular steps:

  1. DSB RESECTION: The broken chromosome ends are resected by nucleases
     (MRN/CtIP, then EXO1 or BLM-DNA2) to expose 3' single-stranded
     overhangs of ~200-2000 nt.

  2. RAD51 FILAMENT FORMATION: The recombinase RAD51 (with mediators BRCA2
     and PALB2) polymerizes along the 3' ssDNA overhang to form a
     nucleoprotein filament.

  3. HOMOLOGY SEARCH: The RAD51 filament searches for a complementary
     sequence. When it encounters the cssDNA donor's homology arm, it
     initiates strand invasion.

  4. STRAND INVASION: The 3' end of the broken chromosome invades the donor,
     base-pairing with the homology arm. This creates a D-loop structure.

  5. DNA SYNTHESIS & RESOLUTION: A DNA polymerase extends from the invading
     3' end, copying the donor's insert region. The D-loop is then resolved
     by SDSA (synthesis-dependent strand annealing) or other mechanisms.

THE CRITICAL REQUIREMENT: HOMOLOGY ARM ACCESSIBILITY
------------------------------------------------------

Step 3 (homology search) REQUIRES that the donor's homology arms be
SINGLE-STRANDED and UNSTRUCTURED. Here is why:

  - RAD51 performs homology search by sampling short segments (8 nt) of the
    donor for complementarity (Qi et al., Cell, 2015). If those 8-nt windows
    are base-paired in a hairpin or sequestered in a G-quadruplex, the RAD51
    filament cannot "read" them, and the homology is invisible.

  - Even if homology is found, strand invasion requires the donor strand to
    be displaced from its intramolecular structure and paired with the
    invading strand. This is thermodynamically unfavorable if the donor's
    secondary structure is stable — the strand invasion must COMPETE with
    the pre-existing hairpin or G4.

QUANTIFYING ACCESSIBILITY:
  We define "accessibility" as the fraction of nucleotides in a homology arm
  that are NOT involved in any secondary structure (hairpin stem or G4 motif).

    Accessibility = 1.0 - (structured_nucleotides / total_arm_nucleotides)

  Interpretation:
    1.0 = perfectly accessible (no secondary structure)
    0.0 = completely inaccessible (entire arm is structured)

  For efficient HDR, both arms should have accessibility > 0.6 (less than
  40% structured). Arms with accessibility < 0.6 are flagged for redesign.

STRUCTURED vs. UNSTRUCTURED NUCLEOTIDES:
  "Structured" nucleotides are those involved in:
    - Hairpin STEMS (base-paired, behaving as dsDNA — invisible to RAD51)
    - G-quadruplex motifs (folded into G4 tetrads — inaccessible for invasion)

  "Unstructured" nucleotides are:
    - Hairpin LOOPS (single-stranded, still accessible despite being in a loop)
    - Regions with no predicted secondary structure
    - Mismatches and bulges in imperfect hairpins

References:
  - Qi et al., Cell, 2015 (8-nt homology sampling by RAD51)
  - Richardson et al., Nature Biotechnology, 2016 (donor design for HDR)
  - Iyer et al., CRISPR Journal, 2022 (cssDNA donor optimization)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict

import numpy as np

from .g_quadruplex import GQuadruplexScanner, G4Motif
from .hairpin import HairpinPredictor, Hairpin


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class AccessibilityReport:
    """
    Comprehensive report on the structural accessibility of a cssDNA donor.

    This report tells the molecular biologist whether the donor template's
    homology arms are sufficiently accessible for productive HDR, and if not,
    what structures are causing problems.

    Attributes
    ----------
    left_arm_accessibility : float
        Fraction of the LEFT homology arm nucleotides that are accessible
        (not involved in secondary structures). Range: 0.0 to 1.0.
        Values > 0.6 are acceptable; < 0.6 means >40% of the arm is
        structured and HDR efficiency is likely to be impaired.

    right_arm_accessibility : float
        Same as above for the RIGHT homology arm.

    overall_score : float
        Weighted average of left and right arm accessibility.
        The minimum of the two arms is weighted more heavily because HDR
        requires BOTH arms to function — the weaker arm is the bottleneck.

    structured_regions : list of (start, end, structure_type)
        All structured regions in the full donor sequence.
        structure_type is one of: "hairpin_stem", "g_quadruplex"
        These regions represent nucleotides that are base-paired or in G4
        structures and therefore inaccessible to RAD51.

    g4_motifs : list of G4Motif
        All G-quadruplex motifs found in the sequence.

    hairpins : list of Hairpin
        All stable hairpins predicted in the sequence.

    warnings : list of str
        Human-readable warnings about structural problems that may impair HDR.
        Empty list = no problems detected.

    recommendation : str
        Overall recommendation for the donor design:
        - "Good donor design" — both arms are highly accessible
        - "Acceptable with caveats" — minor structures present, HDR may be reduced
        - "Consider redesigning" — significant structures in one arm
        - "Redesign strongly recommended" — severe structures in both arms
    """
    left_arm_accessibility: float = 1.0
    right_arm_accessibility: float = 1.0
    overall_score: float = 1.0
    structured_regions: List[Tuple[int, int, str]] = field(default_factory=list)
    g4_motifs: List[G4Motif] = field(default_factory=list)
    hairpins: List[Hairpin] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    recommendation: str = "Good donor design"

    def __repr__(self) -> str:
        return (
            f"AccessibilityReport(\n"
            f"  left_arm={self.left_arm_accessibility:.2f}, "
            f"right_arm={self.right_arm_accessibility:.2f}, "
            f"overall={self.overall_score:.2f},\n"
            f"  structures={len(self.structured_regions)}, "
            f"warnings={len(self.warnings)},\n"
            f"  recommendation='{self.recommendation}'\n"
            f")"
        )


# ---------------------------------------------------------------------------
# Accessibility Scorer
# ---------------------------------------------------------------------------

class AccessibilityScorer:
    """
    Scores the structural accessibility of cssDNA donor homology arms.

    This class integrates G-quadruplex scanning and hairpin prediction to
    determine what fraction of each homology arm is available for RAD51-
    mediated strand invasion during HDR.

    The fundamental question this class answers:
      "Can the RAD51 filament at the DSB find and invade this donor's
       homology arms, or are the arms folded into secondary structures
       that block strand invasion?"

    WORKFLOW:
    ---------
    1. Scan the full donor sequence for G4 motifs (GQuadruplexScanner)
    2. Predict stable hairpins in the full donor sequence (HairpinPredictor)
    3. Merge all "structured" positions (G4 motifs + hairpin stems)
    4. For each homology arm, calculate what fraction of its nucleotides
       fall in structured regions
    5. Generate warnings and recommendations

    Parameters
    ----------
    g4_scanner : GQuadruplexScanner, optional
        Custom G4 scanner. Default scanner is used if not provided.
    hairpin_predictor : HairpinPredictor, optional
        Custom hairpin predictor. Default predictor is used if not provided.
    accessibility_threshold : float
        Minimum acceptable accessibility for a homology arm (default: 0.6).
        Arms below this threshold are flagged as "potentially impaired."
        The value of 0.6 means that up to 40% of the arm can be structured
        before we raise a concern. This threshold is based on the observation
        that RAD51 can tolerate some gaps in its filament as long as
        contiguous ssDNA stretches of >= 15 bp remain for nucleation
        (RAD51_NUCLEATION_MIN_MONOMERS * RAD51_FOOTPRINT_NT from constants.py).
    """

    def __init__(
        self,
        g4_scanner: Optional[GQuadruplexScanner] = None,
        hairpin_predictor: Optional[HairpinPredictor] = None,
        accessibility_threshold: float = 0.6,
    ):
        self.g4_scanner = g4_scanner or GQuadruplexScanner()
        self.hairpin_predictor = hairpin_predictor or HairpinPredictor()
        self.accessibility_threshold = accessibility_threshold

    # -----------------------------------------------------------------
    # Main scoring method
    # -----------------------------------------------------------------

    def score_accessibility(
        self,
        sequence: str,
        left_arm_region: Tuple[int, int],
        right_arm_region: Tuple[int, int],
    ) -> AccessibilityReport:
        """
        Score the structural accessibility of a cssDNA donor's homology arms.

        This is the primary method for evaluating whether a donor template
        is suitable for HDR. It combines G4 and hairpin analysis to provide
        a comprehensive assessment.

        The method follows these steps:

        STEP 1: G-QUADRUPLEX SCANNING
          Find all G4 motifs in the sequence. G4 structures are the most
          thermodynamically stable secondary structures in DNA (Tm > 90 C)
          and pose the greatest threat to HDR efficiency.

        STEP 2: HAIRPIN PREDICTION
          Find all stable hairpins (dG < -3.0 kcal/mol). Hairpins are less
          stable than G4s but much more common, especially in GC-rich regions.

        STEP 3: STRUCTURED POSITION MAPPING
          Merge all "structured" nucleotide positions from both G4 motifs
          and hairpin stems. A nucleotide is "structured" if it is:
            - Within a G4 motif (participating in G-tetrad or G4 loop)
            - Within a hairpin STEM (base-paired, acting as dsDNA)
          Note: hairpin LOOP nucleotides are NOT counted as structured
          because they remain single-stranded and accessible.

        STEP 4: ARM ACCESSIBILITY CALCULATION
          For each homology arm, count how many of its nucleotides are
          structured vs. accessible.

        STEP 5: WARNING GENERATION
          Flag arms that fall below the accessibility threshold (default 0.6).
          Provide specific information about which structures are problematic.

        Parameters
        ----------
        sequence : str
            The full cssDNA donor sequence (A, T, C, G).
        left_arm_region : tuple of (int, int)
            (start, end) indices of the left homology arm in the sequence.
            The left arm is typically the 5' end of the donor, homologous
            to the genomic region upstream of the DSB.
        right_arm_region : tuple of (int, int)
            (start, end) indices of the right homology arm in the sequence.
            The right arm is typically the 3' end of the donor, homologous
            to the genomic region downstream of the DSB.

        Returns
        -------
        AccessibilityReport
            Comprehensive report including per-arm accessibility scores,
            structured regions, warnings, and redesign recommendations.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)

        # STEP 1: Find all G4 motifs.
        g4_motifs = self.g4_scanner.find_g4_motifs(sequence)

        # STEP 2: Find all stable hairpins (non-overlapping, most stable set).
        hairpins = self.hairpin_predictor.predict_structure(sequence)

        # STEP 3: Map structured positions.
        # We track the structure type for each position so we can report
        # what kind of structure is causing the accessibility problem.
        structured_positions: Dict[int, str] = {}

        # 3a. G4 motif positions.
        # ALL nucleotides within a G4 motif are considered structured,
        # including the loop nucleotides, because the G4 loop is constrained
        # within the tetraplex structure and is not freely accessible for
        # strand invasion.
        for motif in g4_motifs:
            # Only count sense-strand G4 motifs for accessibility purposes,
            # because the cssDNA donor IS the sense strand. Antisense motifs
            # would only form on the complementary strand during synthesis.
            if motif.strand == "sense":
                for pos in range(motif.start, motif.end):
                    if 0 <= pos < seq_len:
                        structured_positions[pos] = "g_quadruplex"

        # 3b. Hairpin stem positions.
        # Only the STEM nucleotides are structured (base-paired).
        # Loop nucleotides remain single-stranded and accessible.
        for hp in hairpins:
            # 5' half of stem
            for pos in range(hp.stem_start, hp.loop_start):
                if 0 <= pos < seq_len and pos not in structured_positions:
                    structured_positions[pos] = "hairpin_stem"
            # 3' half of stem
            for pos in range(hp.loop_end, hp.stem_end):
                if 0 <= pos < seq_len and pos not in structured_positions:
                    structured_positions[pos] = "hairpin_stem"

        # Build structured_regions list (contiguous regions of the same type).
        structured_regions = self._build_structured_regions(
            structured_positions, seq_len
        )

        # STEP 4: Calculate per-arm accessibility.
        left_access = self._arm_accessibility(
            structured_positions, left_arm_region
        )
        right_access = self._arm_accessibility(
            structured_positions, right_arm_region
        )

        # Overall score: weighted toward the WEAKER arm, because HDR
        # requires both arms to function. The weaker arm is the bottleneck.
        # Formula: 0.6 * min(left, right) + 0.4 * max(left, right)
        # This penalizes designs where one arm is good but the other is poor.
        overall = (
            0.6 * min(left_access, right_access)
            + 0.4 * max(left_access, right_access)
        )

        # STEP 5: Generate warnings and recommendation.
        warnings = self._generate_warnings(
            left_access, right_access,
            g4_motifs, hairpins,
            left_arm_region, right_arm_region,
        )
        recommendation = self._generate_recommendation(left_access, right_access)

        return AccessibilityReport(
            left_arm_accessibility=left_access,
            right_arm_accessibility=right_access,
            overall_score=overall,
            structured_regions=structured_regions,
            g4_motifs=g4_motifs,
            hairpins=hairpins,
            warnings=warnings,
            recommendation=recommendation,
        )

    # -----------------------------------------------------------------
    # Per-nucleotide accessibility map
    # -----------------------------------------------------------------

    def generate_accessibility_map(self, sequence: str) -> np.ndarray:
        """
        Generate a per-nucleotide accessibility map for the donor sequence.

        This produces a numpy array where each element represents the
        structural accessibility of the corresponding nucleotide:
          1.0 = fully accessible (single-stranded, available for RAD51)
          0.0 = fully structured (in a hairpin stem or G4 motif)

        This array can be visualized as a HEATMAP along the donor sequence,
        giving the molecular biologist an intuitive view of which regions
        are accessible (green/high) and which are structured (red/low).

        The heatmap is particularly useful for:
          - Identifying the exact positions causing accessibility problems
          - Guiding codon optimization (target the red regions)
          - Comparing different donor designs side by side

        Parameters
        ----------
        sequence : str
            The full donor sequence.

        Returns
        -------
        numpy.ndarray
            1D array of float, length = len(sequence), values in [0.0, 1.0].
            Index i corresponds to nucleotide i in the sequence.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)

        # Initialize all positions as fully accessible.
        accessibility = np.ones(seq_len, dtype=float)

        # Find G4 motifs and mark their positions as structured.
        g4_motifs = self.g4_scanner.find_g4_motifs(sequence)
        for motif in g4_motifs:
            if motif.strand == "sense":
                # Scale the accessibility reduction by the motif's stability.
                # More stable G4 → lower accessibility.
                for pos in range(motif.start, min(motif.end, seq_len)):
                    # Use the motif's stability score to modulate accessibility.
                    # A stability of 1.0 → accessibility = 0.0 (fully structured)
                    # A stability of 0.3 → accessibility = 0.7 (mostly accessible)
                    accessibility[pos] = min(
                        accessibility[pos],
                        1.0 - motif.stability_score
                    )

        # Find hairpins and mark stem positions as structured.
        hairpins = self.hairpin_predictor.predict_structure(sequence)
        for hp in hairpins:
            # The degree of inaccessibility depends on the hairpin's stability.
            # We normalize the free energy to a 0-1 scale:
            #   dG = -3.0 kcal/mol → barely stable → accessibility ~ 0.5
            #   dG = -10.0 kcal/mol → very stable → accessibility ~ 0.0
            # Clamp the range to [0.0, 0.7] (never claim a position is
            # completely inaccessible from hairpins alone — RAD51 with
            # BRCA2 can displace some hairpins).
            stability_factor = min(
                0.9,
                abs(hp.free_energy_kcal_mol) / 15.0
            )

            # Mark 5' stem positions.
            for pos in range(hp.stem_start, min(hp.loop_start, seq_len)):
                accessibility[pos] = min(
                    accessibility[pos],
                    1.0 - stability_factor
                )
            # Mark 3' stem positions.
            for pos in range(hp.loop_end, min(hp.stem_end, seq_len)):
                accessibility[pos] = min(
                    accessibility[pos],
                    1.0 - stability_factor
                )
            # Loop positions remain accessible (they are single-stranded).

        return accessibility

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    def _arm_accessibility(
        self,
        structured_positions: Dict[int, str],
        arm_region: Tuple[int, int],
    ) -> float:
        """
        Calculate the accessibility of a single homology arm.

        Accessibility = fraction of arm nucleotides NOT in structured_positions.

        Parameters
        ----------
        structured_positions : dict
            Mapping of position -> structure_type for all structured nucleotides.
        arm_region : tuple of (int, int)
            (start, end) indices of the homology arm.

        Returns
        -------
        float
            Accessibility score from 0.0 (entirely structured) to 1.0 (fully accessible).
        """
        arm_start, arm_end = arm_region
        arm_length = arm_end - arm_start
        if arm_length <= 0:
            return 1.0  # Empty arm is trivially "accessible."

        # Count structured positions within the arm.
        structured_count = sum(
            1 for pos in range(arm_start, arm_end)
            if pos in structured_positions
        )

        return 1.0 - (structured_count / arm_length)

    def _build_structured_regions(
        self,
        structured_positions: Dict[int, str],
        seq_len: int,
    ) -> List[Tuple[int, int, str]]:
        """
        Merge adjacent structured positions into contiguous regions.

        This converts a position-level dictionary into a list of (start, end, type)
        tuples representing contiguous stretches of the same structure type.

        Parameters
        ----------
        structured_positions : dict
            Position -> structure_type mapping.
        seq_len : int
            Total sequence length.

        Returns
        -------
        list of (int, int, str)
            Contiguous structured regions as (start, end, type).
        """
        if not structured_positions:
            return []

        regions: List[Tuple[int, int, str]] = []
        sorted_positions = sorted(structured_positions.keys())

        # Walk through sorted positions, merging adjacent ones of the same type.
        current_start = sorted_positions[0]
        current_end = sorted_positions[0] + 1
        current_type = structured_positions[sorted_positions[0]]

        for pos in sorted_positions[1:]:
            stype = structured_positions[pos]
            if pos == current_end and stype == current_type:
                # Extend the current region.
                current_end = pos + 1
            else:
                # Save the current region and start a new one.
                regions.append((current_start, current_end, current_type))
                current_start = pos
                current_end = pos + 1
                current_type = stype

        # Don't forget the last region.
        regions.append((current_start, current_end, current_type))
        return regions

    def _generate_warnings(
        self,
        left_access: float,
        right_access: float,
        g4_motifs: List[G4Motif],
        hairpins: List[Hairpin],
        left_arm_region: Tuple[int, int],
        right_arm_region: Tuple[int, int],
    ) -> List[str]:
        """
        Generate human-readable warnings about structural problems.

        These warnings are designed for molecular biologists and explain
        the biological consequence of each structural problem.
        """
        warnings: List[str] = []

        # Check left arm accessibility.
        if left_access < self.accessibility_threshold:
            pct_structured = (1.0 - left_access) * 100
            warnings.append(
                f"LEFT homology arm: {pct_structured:.0f}% of nucleotides are "
                f"in secondary structures. This exceeds the {(1 - self.accessibility_threshold) * 100:.0f}% "
                f"threshold and may impair RAD51 filament formation on this arm, "
                f"reducing strand invasion efficiency."
            )

        if right_access < self.accessibility_threshold:
            pct_structured = (1.0 - right_access) * 100
            warnings.append(
                f"RIGHT homology arm: {pct_structured:.0f}% of nucleotides are "
                f"in secondary structures. This exceeds the {(1 - self.accessibility_threshold) * 100:.0f}% "
                f"threshold and may impair RAD51 filament formation on this arm, "
                f"reducing strand invasion efficiency."
            )

        # Check for G4 motifs in homology arms (especially dangerous).
        for motif in g4_motifs:
            if motif.strand != "sense":
                continue
            in_left = (motif.start < left_arm_region[1] and
                       left_arm_region[0] < motif.end)
            in_right = (motif.start < right_arm_region[1] and
                        right_arm_region[0] < motif.end)
            if in_left or in_right:
                arm_name = "LEFT" if in_left else "RIGHT"
                if in_left and in_right:
                    arm_name = "BOTH"
                warnings.append(
                    f"G-QUADRUPLEX in {arm_name} homology arm at position "
                    f"{motif.start}-{motif.end} (stability={motif.stability_score:.2f}). "
                    f"G4 structures are extremely stable (Tm > 90C) and will "
                    f"block RAD51 loading at this position. RAD51 cannot displace "
                    f"G4 structures, making this region permanently inaccessible."
                )

        # Check for very stable hairpins.
        for hp in hairpins:
            if hp.free_energy_kcal_mol < -8.0:
                warnings.append(
                    f"Very stable hairpin at position {hp.stem_start}-{hp.stem_end} "
                    f"(dG = {hp.free_energy_kcal_mol:.1f} kcal/mol, "
                    f"stem = {hp.stem_length}bp). This hairpin has a Tm well above "
                    f"37C and will persist under physiological conditions, sequestering "
                    f"{hp.stem_length * 2} nucleotides from RAD51 binding."
                )

        return warnings

    def _generate_recommendation(
        self,
        left_access: float,
        right_access: float,
    ) -> str:
        """
        Generate an overall recommendation based on arm accessibility scores.

        The recommendation tiers are designed to be actionable:
          - "Good donor design": Proceed with this donor as-is.
          - "Acceptable with caveats": Minor issues; donor may work but consider
            alternatives if HDR efficiency is lower than expected.
          - "Consider redesigning": Significant structural issues likely to
            reduce HDR by >50%. Codon optimization recommended.
          - "Redesign strongly recommended": Both arms severely compromised.
            This donor is unlikely to support efficient HDR.
        """
        min_access = min(left_access, right_access)
        avg_access = (left_access + right_access) / 2.0

        if min_access >= 0.8:
            return "Good donor design"
        elif min_access >= 0.6:
            return "Acceptable with caveats"
        elif min_access >= 0.4 or avg_access >= 0.6:
            return "Consider redesigning"
        else:
            return "Redesign strongly recommended"
