"""
cssDNA-TopoPred: Secondary Structure Analysis for cssDNA Donor Templates
=========================================================================

This module is part of the CRISPRArchitect toolkit. It provides computational
tools for predicting and mitigating secondary structures in circular single-
stranded DNA (cssDNA) donor templates used for homology-directed repair (HDR).

WHY THIS MODULE EXISTS — THE BIOLOGICAL PROBLEM
-------------------------------------------------

cssDNA donors are among the most effective templates for CRISPR-mediated HDR
(Iyer et al., CRISPR Journal, 2022), offering ~3x higher editing efficiency
than linear dsDNA donors. However, cssDNA is ENTIRELY single-stranded, which
means it is free to fold into complex secondary structures:

  - G-QUADRUPLEXES: Four-stranded structures formed by guanine-rich sequences.
    Extremely stable (Tm > 90 C). Block RAD51 loading completely.

  - HAIRPINS (STEM-LOOPS): The ssDNA folds back on itself when complementary
    regions base-pair. Sequesters nucleotides in dsDNA-like stems.

  - HIGHER-ORDER STRUCTURES: Multiple hairpins and G4s can interact to fold
    the donor into compact, non-functional topologies.

These structures are problematic because HDR requires the donor's homology
arms to be SINGLE-STRANDED and ACCESSIBLE:

  1. RAD51 filament formation requires unstructured ssDNA substrate
  2. Homology search requires bases to be "readable" by the RAD51 filament
  3. Strand invasion requires the donor strand to be displaceable

If homology arms are locked in secondary structures, HDR efficiency drops.

MODULE COMPONENTS
------------------

  GQuadruplexScanner   — G4 motif detection and risk scoring
  HairpinPredictor     — Hairpin/stem-loop prediction with free energy
  AccessibilityScorer  — Per-nucleotide and per-arm accessibility scoring
  DonorOptimizer       — Synonymous codon optimization to fix structures
  DonorAnalyzer        — One-stop analysis pipeline (chains everything)

DATA CLASSES
------------

  G4Motif              — A single G-quadruplex motif with stability score
  Hairpin              — A hairpin structure with thermodynamic parameters
  AccessibilityReport  — Comprehensive donor accessibility assessment
  Optimization         — A suggested codon change
  OptimizedResult      — Full optimization result with before/after comparison

QUICK START
-----------

    from crisprarchitect.topopred import DonorAnalyzer

    analyzer = DonorAnalyzer()
    report = analyzer.analyze(
        sequence="ATGCCC...your_donor_sequence...GGGCAT",
        left_arm=(0, 300),
        right_arm=(700, 1000),
        coding_regions=[(300, 700)],
    )
    print(report)

References:
  - Iyer et al., CRISPR Journal, 2022 (cssDNA donor design)
  - Burge et al., Nucleic Acids Res, 2006 (G4 motifs)
  - SantaLucia, PNAS, 1998 (DNA thermodynamics)
  - Richardson et al., Nature Biotechnology, 2016 (ssDNA donor optimization)
"""

from __future__ import annotations

from typing import List, Optional, Tuple, Dict

# --- Import all public classes and data classes ---
from .g_quadruplex import GQuadruplexScanner, G4Motif
from .hairpin import HairpinPredictor, Hairpin
from .accessibility import AccessibilityScorer, AccessibilityReport
from .optimizer import DonorOptimizer, Optimization, OptimizedResult

# --- Public API ---
__all__ = [
    # Core analysis classes
    "GQuadruplexScanner",
    "HairpinPredictor",
    "AccessibilityScorer",
    "DonorOptimizer",
    "DonorAnalyzer",
    # Data classes
    "G4Motif",
    "Hairpin",
    "AccessibilityReport",
    "Optimization",
    "OptimizedResult",
]


# ===========================================================================
# DonorAnalyzer — One-Stop Analysis Pipeline
# ===========================================================================

class DonorAnalyzer:
    """
    One-stop analysis of a cssDNA donor template for secondary structures.

    This class chains together all the components of the TopoPred module to
    provide a comprehensive assessment of a donor template in a single call.
    It is designed for molecular biologists who want a quick answer to:

        "Is my cssDNA donor going to work for HDR, or will secondary
         structures get in the way?"

    THE ANALYSIS PIPELINE:
    ----------------------

    Step 1: G-QUADRUPLEX SCANNING
      The GQuadruplexScanner identifies all G4 motifs in the donor using
      regex-based pattern matching (G{3+}N{1-7}G{3+}N{1-7}G{3+}N{1-7}G{3+}).
      Both the sense and antisense strands are scanned. G4 motifs are scored
      for stability based on the number of G-tetrads and loop lengths.

    Step 2: HAIRPIN PREDICTION
      The HairpinPredictor uses a sliding-window complementarity scan to find
      all possible hairpin structures, then selects the most stable non-
      overlapping set (greedy algorithm). Free energies are calculated using
      a simplified nearest-neighbor model.

    Step 3: ACCESSIBILITY SCORING
      The AccessibilityScorer integrates G4 and hairpin data to compute the
      fraction of each homology arm that is accessible (single-stranded and
      unstructured) vs. structured (base-paired or in G4). Arms with >40%
      structured nucleotides are flagged.

    Step 4: OPTIMIZATION (optional)
      If coding_regions are provided, the DonorOptimizer suggests synonymous
      codon substitutions to disrupt problematic structures. Only coding
      regions are modified (non-coding regions are flagged but not changed).

    Step 5: REPORT GENERATION
      A comprehensive human-readable report is generated with:
        - Per-arm accessibility scores
        - List of all secondary structures found
        - Warnings about critical structures
        - Optimization suggestions (if applicable)
        - Overall recommendation (proceed / redesign / etc.)

    Parameters
    ----------
    g4_scanner : GQuadruplexScanner, optional
        Custom G4 scanner with user-specified parameters.
    hairpin_predictor : HairpinPredictor, optional
        Custom hairpin predictor with user-specified parameters.
    accessibility_scorer : AccessibilityScorer, optional
        Custom accessibility scorer.
    optimizer : DonorOptimizer, optional
        Custom donor optimizer.

    Example
    -------
    >>> analyzer = DonorAnalyzer()
    >>> report = analyzer.analyze(
    ...     sequence="ATGC..." * 250,  # 1000 nt donor
    ...     left_arm=(0, 300),
    ...     right_arm=(700, 1000),
    ...     coding_regions=[(300, 700)],
    ... )
    >>> print(report["recommendation"])
    """

    def __init__(
        self,
        g4_scanner: Optional[GQuadruplexScanner] = None,
        hairpin_predictor: Optional[HairpinPredictor] = None,
        accessibility_scorer: Optional[AccessibilityScorer] = None,
        optimizer: Optional[DonorOptimizer] = None,
    ):
        self.g4_scanner = g4_scanner or GQuadruplexScanner()
        self.hairpin_predictor = hairpin_predictor or HairpinPredictor()
        self.accessibility_scorer = accessibility_scorer or AccessibilityScorer(
            g4_scanner=self.g4_scanner,
            hairpin_predictor=self.hairpin_predictor,
        )
        self.optimizer = optimizer or DonorOptimizer(
            g4_scanner=self.g4_scanner,
            hairpin_predictor=self.hairpin_predictor,
            accessibility_scorer=self.accessibility_scorer,
        )

    def analyze(
        self,
        sequence: str,
        left_arm: Tuple[int, int],
        right_arm: Tuple[int, int],
        coding_regions: Optional[List[Tuple[int, int]]] = None,
    ) -> Dict:
        """
        Perform comprehensive secondary structure analysis of a cssDNA donor.

        This is the main entry point for donor analysis. It runs all four
        analysis steps and returns a structured report.

        Parameters
        ----------
        sequence : str
            The full cssDNA donor sequence (A, T, C, G characters only).
            This should include the left homology arm, insert region, and
            right homology arm.

        left_arm : tuple of (int, int)
            (start, end) indices of the LEFT homology arm.
            The left arm is homologous to the genomic region UPSTREAM
            (5') of the DSB cut site. For cssDNA donors, the optimal
            arm length is ~300 nt (Iyer et al., 2022).

        right_arm : tuple of (int, int)
            (start, end) indices of the RIGHT homology arm.
            The right arm is homologous to the genomic region DOWNSTREAM
            (3') of the DSB cut site.

        coding_regions : list of (int, int), optional
            Coding regions within the donor where synonymous codon
            substitutions are safe. If provided, the optimizer will
            suggest changes to disrupt problematic structures.
            If None, optimization is skipped.

            NOTE: The start position of each coding region must be
            codon-aligned (i.e., the first nucleotide of a codon).

        Returns
        -------
        dict
            Comprehensive analysis report with the following keys:

            "sequence_length" : int
                Length of the input sequence.

            "g4_motifs" : list of G4Motif
                All G-quadruplex motifs found in the sequence.

            "g4_risk_score" : float
                Overall G4 risk score (0-1).

            "g4_arm_map" : dict
                G4 motifs mapped to homology arm regions.

            "hairpins" : list of Hairpin
                All stable hairpins predicted (non-overlapping set).

            "accessibility_report" : AccessibilityReport
                Full accessibility assessment with per-arm scores.

            "accessibility_map" : numpy.ndarray
                Per-nucleotide accessibility values (for heatmap plotting).

            "optimization" : OptimizedResult or None
                Optimization results (if coding_regions were provided).

            "summary" : str
                Human-readable summary of all findings.

            "recommendation" : str
                Overall recommendation for the donor design.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)

        # ---------------------------------------------------------------
        # STEP 1: G-Quadruplex Scanning
        # ---------------------------------------------------------------
        # Scan for G4 motifs that could block RAD51 filament formation.
        # G4s are the most stable secondary structures in DNA and pose
        # the greatest threat to HDR efficiency.
        g4_motifs = self.g4_scanner.find_g4_motifs(sequence)
        g4_risk = self.g4_scanner.score_g4_risk(sequence)
        g4_arm_map = self.g4_scanner.map_g4_positions(
            sequence, (left_arm, right_arm)
        )

        # ---------------------------------------------------------------
        # STEP 2: Hairpin Prediction
        # ---------------------------------------------------------------
        # Find stable hairpins that sequester homology arm nucleotides
        # in base-paired stems, making them unavailable for strand invasion.
        hairpins = self.hairpin_predictor.predict_structure(sequence)

        # ---------------------------------------------------------------
        # STEP 3: Accessibility Scoring
        # ---------------------------------------------------------------
        # Integrate G4 and hairpin data to determine what fraction of
        # each homology arm is accessible for RAD51 binding.
        accessibility_report = self.accessibility_scorer.score_accessibility(
            sequence, left_arm, right_arm
        )
        accessibility_map = self.accessibility_scorer.generate_accessibility_map(
            sequence
        )

        # ---------------------------------------------------------------
        # STEP 4: Optimization (if coding regions provided)
        # ---------------------------------------------------------------
        optimization_result = None
        if coding_regions:
            optimization_result = self.optimizer.optimize_sequence(
                sequence,
                coding_regions,
                left_arm_region=left_arm,
                right_arm_region=right_arm,
            )

        # ---------------------------------------------------------------
        # STEP 5: Summary Generation
        # ---------------------------------------------------------------
        summary = self._generate_summary(
            seq_len=seq_len,
            g4_motifs=g4_motifs,
            g4_risk=g4_risk,
            hairpins=hairpins,
            accessibility_report=accessibility_report,
            optimization_result=optimization_result,
            left_arm=left_arm,
            right_arm=right_arm,
        )

        recommendation = accessibility_report.recommendation

        return {
            "sequence_length": seq_len,
            "g4_motifs": g4_motifs,
            "g4_risk_score": g4_risk,
            "g4_arm_map": g4_arm_map,
            "hairpins": hairpins,
            "accessibility_report": accessibility_report,
            "accessibility_map": accessibility_map,
            "optimization": optimization_result,
            "summary": summary,
            "recommendation": recommendation,
        }

    def _generate_summary(
        self,
        seq_len: int,
        g4_motifs: List[G4Motif],
        g4_risk: float,
        hairpins: List[Hairpin],
        accessibility_report: AccessibilityReport,
        optimization_result: Optional[OptimizedResult],
        left_arm: Tuple[int, int],
        right_arm: Tuple[int, int],
    ) -> str:
        """
        Generate a comprehensive human-readable summary of the analysis.

        This summary is designed for molecular biologists who need to quickly
        understand whether their cssDNA donor template is suitable for HDR.
        """
        lines: List[str] = []

        # Header
        lines.append("=" * 70)
        lines.append("  cssDNA-TopoPred: Donor Template Secondary Structure Analysis")
        lines.append("=" * 70)
        lines.append("")

        # Sequence overview
        lines.append(f"Sequence length: {seq_len} nt")
        lines.append(
            f"Left homology arm:  positions {left_arm[0]}-{left_arm[1]} "
            f"({left_arm[1] - left_arm[0]} nt)"
        )
        lines.append(
            f"Right homology arm: positions {right_arm[0]}-{right_arm[1]} "
            f"({right_arm[1] - right_arm[0]} nt)"
        )
        insert_start = left_arm[1]
        insert_end = right_arm[0]
        lines.append(
            f"Insert region:      positions {insert_start}-{insert_end} "
            f"({insert_end - insert_start} nt)"
        )
        lines.append("")

        # G4 results
        lines.append("-" * 40)
        lines.append("G-QUADRUPLEX ANALYSIS")
        lines.append("-" * 40)
        sense_g4 = [m for m in g4_motifs if m.strand == "sense"]
        antisense_g4 = [m for m in g4_motifs if m.strand == "antisense"]
        lines.append(f"G4 motifs on sense strand:     {len(sense_g4)}")
        lines.append(f"G4 motifs on antisense strand: {len(antisense_g4)}")
        lines.append(f"Overall G4 risk score:         {g4_risk:.2f}")

        if sense_g4:
            lines.append("")
            lines.append("Sense-strand G4 motifs (directly affect donor folding):")
            for m in sense_g4:
                lines.append(
                    f"  Position {m.start}-{m.end}: "
                    f"{m.num_tetrads} tetrads, "
                    f"loops={m.loop_lengths}, "
                    f"stability={m.stability_score:.2f}"
                )
        lines.append("")

        # Hairpin results
        lines.append("-" * 40)
        lines.append("HAIRPIN ANALYSIS")
        lines.append("-" * 40)
        lines.append(f"Stable hairpins found: {len(hairpins)}")

        if hairpins:
            total_structured = sum(hp.stem_length * 2 for hp in hairpins)
            lines.append(
                f"Total nucleotides in stems: {total_structured} "
                f"({total_structured / max(seq_len, 1) * 100:.1f}% of sequence)"
            )
            lines.append("")
            lines.append("Hairpin details:")
            for hp in hairpins:
                lines.append(
                    f"  Position {hp.stem_start}-{hp.stem_end}: "
                    f"stem={hp.stem_length}bp, loop={hp.loop_length}nt, "
                    f"dG={hp.free_energy_kcal_mol:.1f} kcal/mol, "
                    f"GC={hp.gc_content:.0%}"
                )
        lines.append("")

        # Accessibility results
        lines.append("-" * 40)
        lines.append("HOMOLOGY ARM ACCESSIBILITY")
        lines.append("-" * 40)
        lines.append(
            f"Left arm accessibility:  {accessibility_report.left_arm_accessibility:.1%}"
        )
        lines.append(
            f"Right arm accessibility: {accessibility_report.right_arm_accessibility:.1%}"
        )
        lines.append(
            f"Overall accessibility:   {accessibility_report.overall_score:.1%}"
        )

        if accessibility_report.warnings:
            lines.append("")
            lines.append("WARNINGS:")
            for w in accessibility_report.warnings:
                lines.append(f"  ! {w}")
        lines.append("")

        # Optimization results
        if optimization_result is not None:
            lines.append("-" * 40)
            lines.append("OPTIMIZATION RESULTS")
            lines.append("-" * 40)
            n_changes = len(optimization_result.changes_made)
            lines.append(f"Synonymous codon changes applied: {n_changes}")

            if n_changes > 0:
                lines.append(
                    f"Left arm improvement:  "
                    f"{optimization_result.improvement_left:+.1%}"
                )
                lines.append(
                    f"Right arm improvement: "
                    f"{optimization_result.improvement_right:+.1%}"
                )
                lines.append("")
                lines.append("Changes made:")
                for opt in optimization_result.changes_made:
                    lines.append(
                        f"  Position {opt.position}: "
                        f"{opt.original_codon} -> {opt.suggested_codon} "
                        f"({opt.amino_acid}, {opt.target_structure})"
                    )

            if optimization_result.non_coding_warnings:
                lines.append("")
                lines.append("Non-coding region warnings:")
                for w in optimization_result.non_coding_warnings:
                    lines.append(f"  ! {w}")
            lines.append("")

        # Recommendation
        lines.append("=" * 40)
        lines.append(f"RECOMMENDATION: {accessibility_report.recommendation}")
        lines.append("=" * 40)

        return "\n".join(lines)


# ===========================================================================
# Example: Analyze a Realistic cssDNA Donor Template
# ===========================================================================

if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # EXAMPLE: Analyzing a cssDNA donor for a GFP knock-in experiment
    # -----------------------------------------------------------------------
    #
    # This example models a ~960 nt cssDNA donor with:
    #   - Left homology arm (300 nt): genomic sequence upstream of DSB
    #   - Insert region (~360 nt): partial GFP coding sequence
    #   - Right homology arm (300 nt): genomic sequence downstream of DSB
    #
    # The donor contains several deliberate features for demonstration:
    #   1. A G4 motif in the left homology arm (simulates a G-rich genomic locus)
    #   2. A stable hairpin in the right homology arm (complementary sequences)
    #   3. Normal sequence in the insert region
    #
    # In practice, the user would supply their actual donor sequence.
    # -----------------------------------------------------------------------

    # --- Left homology arm: 300 nt ---
    # Contains a G4-forming sequence around position 100-130.
    # This simulates a donor for a G-rich genomic locus (common in promoter regions).
    left_arm = (
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 0-59
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"                       # 60-99
        "GGGTAGGGTTAGGGTTAGGGTTTTTTTTTTTTTTTTTTTT"                       # 100-139 G4 motif!
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 140-199
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 200-259
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"                       # 260-299
    )

    # --- Insert region: ~360 nt (partial GFP CDS) ---
    # This is a coding region — synonymous codon optimization is safe here.
    insert = (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"   # 300-359
        "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"   # 360-419
        "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"   # 420-479
        "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"   # 480-539
        "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"   # 540-599
        "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"   # 600-659
    )

    # --- Right homology arm: 300 nt ---
    # Contains a region that forms a stable hairpin around position 720-780.
    # The first 30 nt are complementary to nt at positions 740-770 (reversed).
    right_arm = (
        "GCTAGCTAGCTAGCTAGCAATTGCCGATCGA"                                 # 660-689
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 690-749
        "TCGATCGGCAATTGCTAGCTAGCTAGCTAGC"                                 # 750-779 (complement of 660-689!)
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 780-839
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 840-899
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"   # 900-959
    )

    # Assemble the full donor.
    donor_sequence = left_arm + insert + right_arm
    print(f"Donor sequence length: {len(donor_sequence)} nt")
    print()

    # Define regions.
    left_arm_region = (0, 300)
    right_arm_region = (660, 960)
    coding_region = [(300, 660)]  # The GFP insert is a coding sequence.

    # Run the analysis.
    analyzer = DonorAnalyzer()
    report = analyzer.analyze(
        sequence=donor_sequence,
        left_arm=left_arm_region,
        right_arm=right_arm_region,
        coding_regions=coding_region,
    )

    # Print the full summary.
    print(report["summary"])
    print()

    # Print additional details.
    print(f"Number of G4 motifs found: {len(report['g4_motifs'])}")
    print(f"Number of stable hairpins: {len(report['hairpins'])}")
    print(f"G4 risk score: {report['g4_risk_score']:.2f}")
    print()

    # Show the accessibility map (first 50 positions).
    acc_map = report["accessibility_map"]
    print("Accessibility map (first 50 nt, 1.0=accessible, 0.0=structured):")
    for i in range(0, min(50, len(acc_map)), 10):
        chunk = acc_map[i:i+10]
        values = " ".join(f"{v:.1f}" for v in chunk)
        print(f"  nt {i:3d}-{i+len(chunk)-1:3d}: {values}")
    print()

    # Show optimization results if any changes were made.
    if report["optimization"] and report["optimization"].changes_made:
        print("Optimization suggestions applied:")
        for opt in report["optimization"].changes_made:
            print(f"  {opt.rationale}")
        print()
        print(
            f"Optimized left arm accessibility: "
            f"{report['optimization'].optimized_report.left_arm_accessibility:.1%}"
        )
        print(
            f"Optimized right arm accessibility: "
            f"{report['optimization'].optimized_report.right_arm_accessibility:.1%}"
        )
    else:
        print("No optimization changes were needed or possible.")
