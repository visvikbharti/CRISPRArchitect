"""
Donor Sequence Optimizer for Reducing Problematic Secondary Structures
========================================================================

BIOLOGICAL BACKGROUND — SYNONYMOUS CODON OPTIMIZATION
------------------------------------------------------

When a cssDNA donor template contains problematic secondary structures
(G-quadruplexes or stable hairpins) in its homology arms or insert region,
we can sometimes fix these structures by making SYNONYMOUS CODON CHANGES.

WHAT ARE SYNONYMOUS CODONS?
  The genetic code is DEGENERATE — most amino acids are encoded by multiple
  codons. For example, glycine (Gly) is encoded by GGU, GGC, GGA, and GGG.
  Changing GGG to GGA in a coding sequence still produces glycine, but it
  replaces one G with an A, which can disrupt a G-run that forms part of a
  G-quadruplex motif.

  This is the key insight: we can modify the DNA sequence WITHOUT changing
  the encoded protein, specifically to disrupt problematic secondary
  structures.

WHY ONLY CODING REGIONS?
  Synonymous changes are "safe" only in protein-coding regions, where the
  genetic code provides a mapping from codon to amino acid. In non-coding
  regions (introns, UTRs, regulatory elements):
    - Introns contain splice signals (GT...AG, branch point sequences) that
      are essential for correct mRNA processing. Mutations here can cause
      exon skipping, intron retention, or cryptic splice site activation.
    - 5' UTRs contain regulatory elements (Kozak sequence, upstream ORFs,
      IRES elements) that control translation initiation.
    - 3' UTRs contain mRNA stability signals, miRNA binding sites, and
      polyadenylation signals.
    - Enhancers and silencers contain transcription factor binding sites.

  Therefore, for non-coding regions, we FLAG the structural problem but
  DO NOT automatically suggest nucleotide changes. The molecular biologist
  must evaluate whether changes are safe on a case-by-case basis.

OPTIMIZATION STRATEGY FOR G-QUADRUPLEXES:
  A G4 motif requires four G-runs of 3+ guanines: G{3+}N{1-7}G{3+}N{1-7}G{3+}N{1-7}G{3+}
  To disrupt a G4, we need to break at least one G-run below 3 consecutive Gs.
  Strategy: In each G-run that overlaps a coding region, find a codon containing
  GG or GGG and substitute a synonymous codon with fewer Gs.

  Example:
    Original:  ...GGG GGC GGA...  (Gly-Gly-Gly)
    Optimized: ...GGA GGC GGA...  (Gly-Gly-Gly)  — same protein, fewer G-runs

OPTIMIZATION STRATEGY FOR HAIRPINS:
  A hairpin requires complementary sequences. To disrupt a stem:
  Strategy: In one half of the stem that overlaps a coding region, substitute
  codons to break complementarity.

  Example:
    Original:  ...ATG GCC TAA... paired with ...TTA GGC CAT...
    Optimized: ...ATG GCT TAA... (Ala→Ala, but C→T breaks one bp in the stem)

SAFETY CONSIDERATIONS:
  Even in coding regions, not all synonymous changes are truly neutral:
    - Codon usage bias: Some synonymous codons are translated more efficiently
      than others. We preferentially choose common codons.
    - Exonic splicing enhancers/silencers (ESE/ESS): Some exonic sequences
      regulate splicing. Disrupting an ESE can cause exon skipping.
    - mRNA secondary structure: Changes can affect mRNA folding and stability.

  These risks are MUCH lower than the HDR impairment from a stable G4 or
  hairpin, but users should be aware of them. The optimizer flags these
  caveats in its output.

References:
  - Plotkin & Kudla, Nature Reviews Genetics, 2011 (codon usage bias)
  - Cartegni et al., Nature Reviews Genetics, 2002 (exonic splicing signals)
  - Kudla et al., Science, 2009 (codon optimization and mRNA structure)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Set

from .g_quadruplex import GQuadruplexScanner, G4Motif
from .hairpin import HairpinPredictor, Hairpin
from .accessibility import AccessibilityScorer, AccessibilityReport


# ===========================================================================
# STANDARD GENETIC CODE
# ===========================================================================
#
# The universal genetic code maps 64 codons to 20 amino acids + 1 stop signal.
# Each entry is codon -> amino acid (single letter code).
# Stop codons are represented as '*'.
#
# This dictionary is the foundation of synonymous codon substitution:
# for any codon, we can look up its amino acid and find all other codons
# that encode the same amino acid (synonymous codons).

GENETIC_CODE: Dict[str, str] = {
    # Phenylalanine (F) — 2 codons
    "TTT": "F", "TTC": "F",
    # Leucine (L) — 6 codons
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine (I) — 3 codons
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine (M) — 1 codon (no synonymous substitution possible)
    "ATG": "M",
    # Valine (V) — 4 codons
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine (S) — 6 codons
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline (P) — 4 codons
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine (T) — 4 codons
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine (A) — 4 codons
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine (Y) — 2 codons
    "TAT": "Y", "TAC": "Y",
    # Stop codons (*)
    "TAA": "*", "TAG": "*", "TGA": "*",
    # Histidine (H) — 2 codons
    "CAT": "H", "CAC": "H",
    # Glutamine (Q) — 2 codons
    "CAA": "Q", "CAG": "Q",
    # Asparagine (N) — 2 codons
    "AAT": "N", "AAC": "N",
    # Lysine (K) — 2 codons
    "AAA": "K", "AAG": "K",
    # Aspartic acid (D) — 2 codons
    "GAT": "D", "GAC": "D",
    # Glutamic acid (E) — 2 codons
    "GAA": "E", "GAG": "E",
    # Cysteine (C) — 2 codons
    "TGT": "C", "TGC": "C",
    # Tryptophan (W) — 1 codon (no synonymous substitution possible)
    "TGG": "W",
    # Arginine (R) — 6 codons
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine (G) — 4 codons
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Build a reverse map: amino acid -> list of codons
# This is used to find all synonymous alternatives for a given codon.
AMINO_ACID_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in GENETIC_CODE.items():
    AMINO_ACID_TO_CODONS.setdefault(_aa, []).append(_codon)

# Human codon usage frequency (approximate, codons per 1000).
# Source: Kazusa Codon Usage Database for Homo sapiens.
# Higher-frequency codons are preferred because they are translated more
# efficiently by the ribosome (more abundant tRNAs).
# We use these to choose the "best" synonymous alternative.
HUMAN_CODON_FREQUENCY: Dict[str, float] = {
    "TTT": 17.6, "TTC": 20.3,
    "TTA": 7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5,
    "ATG": 22.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4, "AGT": 12.1, "AGC": 19.5,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 12.2, "TAC": 15.3,
    "TAA": 1.0, "TAG": 0.8, "TGA": 1.6,
    "CAT": 10.9, "CAC": 15.1,
    "CAA": 12.3, "CAG": 34.2,
    "AAT": 17.0, "AAC": 19.1,
    "AAA": 24.4, "AAG": 31.9,
    "GAT": 21.8, "GAC": 25.1,
    "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6,
    "TGG": 13.2,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4, "AGA": 12.2, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class Optimization:
    """
    Represents a single suggested codon optimization.

    Each optimization describes a specific codon change that would help
    disrupt a problematic secondary structure (G4 or hairpin) while
    preserving the encoded amino acid.

    Attributes
    ----------
    position : int
        0-based position of the FIRST nucleotide of the codon in the
        donor sequence.
    original_codon : str
        The original codon (3 nucleotides).
    suggested_codon : str
        The suggested replacement codon (same amino acid, different sequence).
    amino_acid : str
        The amino acid encoded by both codons (single letter code).
    target_structure : str
        Which structural problem this optimization addresses:
        "g_quadruplex" or "hairpin"
    structure_position : Tuple[int, int]
        (start, end) of the problematic structure in the donor sequence.
    rationale : str
        Human-readable explanation of why this change helps.
    safety_notes : List[str]
        Any caveats about this change (codon usage, ESE disruption, etc.).
    """
    position: int = 0
    original_codon: str = ""
    suggested_codon: str = ""
    amino_acid: str = ""
    target_structure: str = ""
    structure_position: Tuple[int, int] = (0, 0)
    rationale: str = ""
    safety_notes: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"Optimization(pos={self.position}, "
            f"{self.original_codon}->{self.suggested_codon} [{self.amino_acid}], "
            f"target='{self.target_structure}')"
        )


@dataclass
class OptimizedResult:
    """
    Result of donor sequence optimization.

    Attributes
    ----------
    original_sequence : str
        The original donor sequence before optimization.
    optimized_sequence : str
        The modified donor sequence after applying safe optimizations.
    changes_made : List[Optimization]
        All codon changes that were applied.
    original_report : AccessibilityReport
        Accessibility report for the original sequence.
    optimized_report : AccessibilityReport
        Accessibility report for the optimized sequence.
    improvement_left : float
        Change in left arm accessibility (positive = improvement).
    improvement_right : float
        Change in right arm accessibility (positive = improvement).
    non_coding_warnings : List[str]
        Warnings about problematic structures in non-coding regions
        that could not be automatically optimized.
    """
    original_sequence: str = ""
    optimized_sequence: str = ""
    changes_made: List[Optimization] = field(default_factory=list)
    original_report: Optional[AccessibilityReport] = None
    optimized_report: Optional[AccessibilityReport] = None
    improvement_left: float = 0.0
    improvement_right: float = 0.0
    non_coding_warnings: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"OptimizedResult(\n"
            f"  changes={len(self.changes_made)},\n"
            f"  left_arm: {self.original_report.left_arm_accessibility:.2f} -> "
            f"{self.optimized_report.left_arm_accessibility:.2f} "
            f"(+{self.improvement_left:.2f}),\n"
            f"  right_arm: {self.original_report.right_arm_accessibility:.2f} -> "
            f"{self.optimized_report.right_arm_accessibility:.2f} "
            f"(+{self.improvement_right:.2f}),\n"
            f"  non_coding_warnings={len(self.non_coding_warnings)}\n"
            f")"
        )


# ---------------------------------------------------------------------------
# Donor Optimizer
# ---------------------------------------------------------------------------

class DonorOptimizer:
    """
    Optimizes cssDNA donor sequences to reduce problematic secondary structures.

    This optimizer makes SYNONYMOUS codon substitutions in coding regions to
    disrupt G-quadruplex motifs and hairpin stems, while preserving the encoded
    protein sequence. For non-coding regions, it flags problems but does not
    make automatic changes (too risky without domain knowledge).

    HOW THE OPTIMIZER WORKS:
    ------------------------

    1. IDENTIFY PROBLEMS:
       Use GQuadruplexScanner and HairpinPredictor to find all secondary
       structures in the donor sequence.

    2. LOCATE CODING OVERLAPS:
       For each problematic structure, check if it overlaps with any
       user-specified coding region. Only coding-region overlaps are
       candidates for automatic optimization.

    3. FIND SYNONYMOUS ALTERNATIVES:
       For each codon in the problematic region, look up its amino acid
       and find all synonymous codons. Score alternatives by:
         a. How many problematic bases they change (more disruption = better)
         b. Human codon usage frequency (prefer common codons)
         c. Reduction in G-content (for G4 disruption)
         d. Reduction in complementarity (for hairpin disruption)

    4. APPLY AND VERIFY:
       Apply the best changes and re-score accessibility to confirm improvement.

    Parameters
    ----------
    g4_scanner : GQuadruplexScanner, optional
        Custom G4 scanner.
    hairpin_predictor : HairpinPredictor, optional
        Custom hairpin predictor.
    accessibility_scorer : AccessibilityScorer, optional
        Custom accessibility scorer.
    """

    def __init__(
        self,
        g4_scanner: Optional[GQuadruplexScanner] = None,
        hairpin_predictor: Optional[HairpinPredictor] = None,
        accessibility_scorer: Optional[AccessibilityScorer] = None,
    ):
        self.g4_scanner = g4_scanner or GQuadruplexScanner()
        self.hairpin_predictor = hairpin_predictor or HairpinPredictor()
        self.accessibility_scorer = accessibility_scorer or AccessibilityScorer(
            g4_scanner=self.g4_scanner,
            hairpin_predictor=self.hairpin_predictor,
        )

    # -----------------------------------------------------------------
    # Suggest optimizations
    # -----------------------------------------------------------------

    def suggest_optimizations(
        self,
        sequence: str,
        coding_regions: List[Tuple[int, int]],
        problems: Optional[List] = None,
    ) -> List[Optimization]:
        """
        Suggest synonymous codon changes to disrupt problematic structures.

        This method analyzes each problematic secondary structure and, for
        those overlapping coding regions, suggests specific codon substitutions
        that would disrupt the structure while preserving the protein.

        The suggestions are CONSERVATIVE — they prioritize safety:
          - Only synonymous codons are suggested (same amino acid)
          - Common codons (high human usage frequency) are preferred
          - Stop codons are never introduced
          - Changes near splice sites are flagged with caution
          - The minimum number of changes needed is suggested

        Parameters
        ----------
        sequence : str
            The donor DNA sequence.
        coding_regions : list of (int, int)
            List of (start, end) tuples defining coding regions in the donor.
            These are the ONLY regions where automatic synonymous substitutions
            are safe. The start position must be codon-aligned (i.e., the first
            position of a codon in the reading frame).
        problems : list, optional
            Pre-computed list of problematic structures (G4Motif and Hairpin objects).
            If None, the method will find them automatically.

        Returns
        -------
        list of Optimization
            Suggested codon changes, sorted by position. Each optimization
            includes the rationale and safety notes.
        """
        sequence = sequence.upper().strip()
        suggestions: List[Optimization] = []

        # Find problems if not provided.
        if problems is None:
            g4_motifs = self.g4_scanner.find_g4_motifs(sequence)
            # Only sense-strand G4 motifs are relevant for the donor.
            g4_motifs = [m for m in g4_motifs if m.strand == "sense"]
            hairpins = self.hairpin_predictor.predict_structure(sequence)
        else:
            g4_motifs = [p for p in problems if isinstance(p, G4Motif) and p.strand == "sense"]
            hairpins = [p for p in problems if isinstance(p, Hairpin)]

        # --- Process G4 motifs ---
        for motif in g4_motifs:
            motif_suggestions = self._suggest_g4_disruption(
                sequence, motif, coding_regions
            )
            suggestions.extend(motif_suggestions)

        # --- Process hairpins ---
        for hp in hairpins:
            if hp.is_stable:
                hp_suggestions = self._suggest_hairpin_disruption(
                    sequence, hp, coding_regions
                )
                suggestions.extend(hp_suggestions)

        # Remove duplicate suggestions (same position).
        seen_positions: Set[int] = set()
        unique_suggestions: List[Optimization] = []
        for opt in suggestions:
            if opt.position not in seen_positions:
                seen_positions.add(opt.position)
                unique_suggestions.append(opt)

        # Sort by position.
        unique_suggestions.sort(key=lambda o: o.position)
        return unique_suggestions

    # -----------------------------------------------------------------
    # Full optimization pipeline
    # -----------------------------------------------------------------

    def optimize_sequence(
        self,
        sequence: str,
        coding_regions: List[Tuple[int, int]],
        left_arm_region: Optional[Tuple[int, int]] = None,
        right_arm_region: Optional[Tuple[int, int]] = None,
    ) -> OptimizedResult:
        """
        Apply safe optimizations to reduce secondary structures in a donor.

        This method runs the full optimization pipeline:
          1. Score the original sequence's accessibility
          2. Find and suggest optimizations
          3. Apply the safe codon changes
          4. Re-score the optimized sequence
          5. Report the improvement

        Parameters
        ----------
        sequence : str
            The original donor sequence.
        coding_regions : list of (int, int)
            Coding regions where synonymous changes are safe.
        left_arm_region : tuple of (int, int), optional
            Left homology arm boundaries. If None, defaults to first 300 nt.
        right_arm_region : tuple of (int, int), optional
            Right homology arm boundaries. If None, defaults to last 300 nt.

        Returns
        -------
        OptimizedResult
            Complete result with original/optimized sequences, changes made,
            before/after accessibility scores, and non-coding warnings.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)

        # Default arm regions if not specified.
        # Convention: left arm = first 300 nt, right arm = last 300 nt.
        # (300 nt is the optimal HA length for cssDNA; see constants.py)
        if left_arm_region is None:
            left_arm_region = (0, min(300, seq_len))
        if right_arm_region is None:
            right_arm_region = (max(0, seq_len - 300), seq_len)

        # Score original accessibility.
        original_report = self.accessibility_scorer.score_accessibility(
            sequence, left_arm_region, right_arm_region
        )

        # Get optimization suggestions.
        all_problems = original_report.g4_motifs + original_report.hairpins
        suggestions = self.suggest_optimizations(
            sequence, coding_regions, problems=all_problems
        )

        # Apply changes to the sequence.
        # We apply changes from right to left (highest position first)
        # to avoid position shifts when we modify the string.
        optimized_seq = list(sequence)
        applied_changes: List[Optimization] = []

        for opt in sorted(suggestions, key=lambda o: -o.position):
            pos = opt.position
            original_codon = opt.original_codon
            new_codon = opt.suggested_codon

            # Verify the original codon is still at the expected position.
            # (It might have been shifted by a previous change, but since
            # we work right-to-left, this shouldn't happen.)
            current = "".join(optimized_seq[pos:pos + 3])
            if current == original_codon:
                optimized_seq[pos:pos + 3] = list(new_codon)
                applied_changes.append(opt)

        optimized_sequence = "".join(optimized_seq)

        # Re-score optimized accessibility.
        optimized_report = self.accessibility_scorer.score_accessibility(
            optimized_sequence, left_arm_region, right_arm_region
        )

        # Calculate improvement.
        improvement_left = (
            optimized_report.left_arm_accessibility
            - original_report.left_arm_accessibility
        )
        improvement_right = (
            optimized_report.right_arm_accessibility
            - original_report.right_arm_accessibility
        )

        # Flag non-coding structural problems.
        non_coding_warnings = self._flag_noncoding_problems(
            sequence, coding_regions, all_problems
        )

        return OptimizedResult(
            original_sequence=sequence,
            optimized_sequence=optimized_sequence,
            changes_made=list(reversed(applied_changes)),  # Restore left-to-right order
            original_report=original_report,
            optimized_report=optimized_report,
            improvement_left=improvement_left,
            improvement_right=improvement_right,
            non_coding_warnings=non_coding_warnings,
        )

    # -----------------------------------------------------------------
    # G4 disruption suggestions
    # -----------------------------------------------------------------

    def _suggest_g4_disruption(
        self,
        sequence: str,
        motif: G4Motif,
        coding_regions: List[Tuple[int, int]],
    ) -> List[Optimization]:
        """
        Suggest synonymous codon changes to disrupt a G-quadruplex motif.

        STRATEGY:
        To destroy a G4, we need to break at least ONE G-run below 3
        consecutive guanines. We scan the G4 region for codons in coding
        regions that contain G-rich codons (like GGG, GGC, GGA, GGT) and
        suggest alternatives with fewer Gs.

        The ideal substitution:
          GGG (Gly) -> GGA (Gly): replaces one G with A, potentially
          breaking a G-run from GGGGG to GGG-A-G, which disrupts the G4.

        Parameters
        ----------
        sequence : str
            The donor sequence.
        motif : G4Motif
            The G4 motif to disrupt.
        coding_regions : list of (int, int)
            Coding regions where changes are safe.

        Returns
        -------
        list of Optimization
            Suggested changes to disrupt this G4.
        """
        suggestions: List[Optimization] = []

        # Find all codons within the G4 motif that are in coding regions.
        for cds_start, cds_end in coding_regions:
            # Iterate over codons in this coding region.
            # Codons are 3-nt blocks starting from cds_start.
            codon_start = cds_start
            while codon_start + 3 <= cds_end:
                codon_end = codon_start + 3

                # Check if this codon overlaps with the G4 motif.
                if codon_start < motif.end and motif.start < codon_end:
                    codon_seq = sequence[codon_start:codon_end]

                    if codon_seq not in GENETIC_CODE:
                        codon_start += 3
                        continue

                    amino_acid = GENETIC_CODE[codon_seq]
                    if amino_acid == "*":
                        # Don't modify stop codons.
                        codon_start += 3
                        continue

                    # Count Gs in this codon — we want to reduce G content.
                    g_count = codon_seq.count("G")
                    if g_count == 0:
                        codon_start += 3
                        continue

                    # Find synonymous alternatives with fewer Gs.
                    alternatives = AMINO_ACID_TO_CODONS.get(amino_acid, [])
                    best_alt = None
                    best_score = -1

                    for alt_codon in alternatives:
                        if alt_codon == codon_seq:
                            continue  # Skip the original codon.

                        alt_g_count = alt_codon.count("G")
                        if alt_g_count >= g_count:
                            continue  # Not an improvement.

                        # Score: prefer fewer Gs and higher usage frequency.
                        g_reduction = g_count - alt_g_count
                        usage = HUMAN_CODON_FREQUENCY.get(alt_codon, 5.0)
                        # Normalize usage to [0, 1] range (max ~40).
                        usage_norm = usage / 40.0
                        score = 0.6 * g_reduction + 0.4 * usage_norm

                        if score > best_score:
                            best_score = score
                            best_alt = alt_codon

                    if best_alt is not None:
                        suggestions.append(Optimization(
                            position=codon_start,
                            original_codon=codon_seq,
                            suggested_codon=best_alt,
                            amino_acid=amino_acid,
                            target_structure="g_quadruplex",
                            structure_position=(motif.start, motif.end),
                            rationale=(
                                f"Replace {codon_seq} ({amino_acid}) with "
                                f"{best_alt} ({amino_acid}) to reduce G-content "
                                f"and disrupt G-quadruplex at position "
                                f"{motif.start}-{motif.end}. This change removes "
                                f"{codon_seq.count('G') - best_alt.count('G')} "
                                f"guanine(s) from the G-run."
                            ),
                            safety_notes=self._codon_safety_notes(
                                codon_seq, best_alt, amino_acid
                            ),
                        ))

                codon_start += 3

        return suggestions

    # -----------------------------------------------------------------
    # Hairpin disruption suggestions
    # -----------------------------------------------------------------

    def _suggest_hairpin_disruption(
        self,
        sequence: str,
        hairpin: Hairpin,
        coding_regions: List[Tuple[int, int]],
    ) -> List[Optimization]:
        """
        Suggest synonymous codon changes to disrupt a hairpin stem.

        STRATEGY:
        To weaken a hairpin, we need to introduce mismatches in the stem.
        We focus on the 5' half of the stem (positions stem_start to loop_start)
        and look for codons in coding regions where a synonymous change would
        break a base pair with the corresponding position in the 3' half.

        Even a single mismatch in a short stem (4-6 bp) can destabilize the
        hairpin enough to prevent it from forming at 37 C.

        Parameters
        ----------
        sequence : str
            The donor sequence.
        hairpin : Hairpin
            The hairpin to disrupt.
        coding_regions : list of (int, int)
            Coding regions where changes are safe.

        Returns
        -------
        list of Optimization
            Suggested changes to disrupt this hairpin.
        """
        suggestions: List[Optimization] = []

        # We target the 5' half of the stem for modifications.
        # Each position i in the 5' stem pairs with position j in the 3' stem,
        # where j = loop_end + (loop_start - 1 - i).
        # Changing position i to break complementarity with j weakens the stem.

        for cds_start, cds_end in coding_regions:
            codon_start = cds_start
            while codon_start + 3 <= cds_end:
                codon_end = codon_start + 3

                # Check if this codon overlaps with the 5' stem half.
                overlaps_5prime = (
                    codon_start < hairpin.loop_start
                    and hairpin.stem_start < codon_end
                )
                # Also check overlap with 3' stem half.
                overlaps_3prime = (
                    codon_start < hairpin.stem_end
                    and hairpin.loop_end < codon_end
                )

                if not overlaps_5prime and not overlaps_3prime:
                    codon_start += 3
                    continue

                codon_seq = sequence[codon_start:codon_end]
                if codon_seq not in GENETIC_CODE:
                    codon_start += 3
                    continue

                amino_acid = GENETIC_CODE[codon_seq]
                if amino_acid == "*":
                    codon_start += 3
                    continue

                # Find synonymous alternatives that break complementarity
                # at the stem positions this codon covers.
                alternatives = AMINO_ACID_TO_CODONS.get(amino_acid, [])
                best_alt = None
                best_score = -1

                for alt_codon in alternatives:
                    if alt_codon == codon_seq:
                        continue

                    # Count how many stem base pairs this change would break.
                    mismatches_introduced = self._count_stem_disruption(
                        sequence, codon_start, codon_seq, alt_codon, hairpin
                    )

                    if mismatches_introduced == 0:
                        continue  # This alternative doesn't help.

                    usage = HUMAN_CODON_FREQUENCY.get(alt_codon, 5.0)
                    usage_norm = usage / 40.0
                    score = 0.6 * mismatches_introduced + 0.4 * usage_norm

                    if score > best_score:
                        best_score = score
                        best_alt = alt_codon

                if best_alt is not None:
                    mismatches = self._count_stem_disruption(
                        sequence, codon_start, codon_seq, best_alt, hairpin
                    )
                    suggestions.append(Optimization(
                        position=codon_start,
                        original_codon=codon_seq,
                        suggested_codon=best_alt,
                        amino_acid=amino_acid,
                        target_structure="hairpin",
                        structure_position=(hairpin.stem_start, hairpin.stem_end),
                        rationale=(
                            f"Replace {codon_seq} ({amino_acid}) with "
                            f"{best_alt} ({amino_acid}) to introduce "
                            f"{mismatches} mismatch(es) in the hairpin stem at "
                            f"position {hairpin.stem_start}-{hairpin.stem_end} "
                            f"(dG = {hairpin.free_energy_kcal_mol:.1f} kcal/mol). "
                            f"Breaking stem base pairs destabilizes the hairpin, "
                            f"freeing the nucleotides for RAD51 binding."
                        ),
                        safety_notes=self._codon_safety_notes(
                            codon_seq, best_alt, amino_acid
                        ),
                    ))

                codon_start += 3

        return suggestions

    # -----------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------

    def _count_stem_disruption(
        self,
        sequence: str,
        codon_start: int,
        original_codon: str,
        new_codon: str,
        hairpin: Hairpin,
    ) -> int:
        """
        Count how many hairpin stem base pairs a codon change would break.

        For each position in the codon that overlaps the hairpin stem:
          - Find the paired position in the opposite stem half
          - Check if the original base pairs with the partner
          - Check if the new base STILL pairs with the partner
          - If the change breaks the pairing, count it as a disruption

        Parameters
        ----------
        sequence : str
            The full donor sequence.
        codon_start : int
            Position of the first nucleotide of the codon.
        original_codon : str
            The original 3-nt codon.
        new_codon : str
            The proposed alternative codon.
        hairpin : Hairpin
            The hairpin being targeted.

        Returns
        -------
        int
            Number of stem base pairs disrupted by this codon change.
        """
        _WC = {"A": "T", "T": "A", "G": "C", "C": "G"}
        disruptions = 0
        seq_len = len(sequence)

        for offset in range(3):
            pos = codon_start + offset
            orig_base = original_codon[offset]
            new_base = new_codon[offset]

            if orig_base == new_base:
                continue  # No change at this position.

            # Determine if this position is in the 5' or 3' stem half,
            # and find its paired position.
            paired_pos = None

            if hairpin.stem_start <= pos < hairpin.loop_start:
                # Position is in the 5' half of the stem.
                # It pairs with the 3' half in antiparallel orientation.
                # The i-th base from the 5' stem start pairs with the i-th
                # base from the 3' stem END (counting inward).
                offset_from_5prime = pos - hairpin.stem_start
                paired_pos = hairpin.stem_end - 1 - offset_from_5prime

            elif hairpin.loop_end <= pos < hairpin.stem_end:
                # Position is in the 3' half of the stem.
                offset_from_3prime_end = hairpin.stem_end - 1 - pos
                paired_pos = hairpin.stem_start + offset_from_3prime_end

            if paired_pos is None or paired_pos < 0 or paired_pos >= seq_len:
                continue

            partner_base = sequence[paired_pos]

            # Check if the original base paired with the partner.
            orig_pairs = (_WC.get(orig_base) == partner_base)
            # Check if the new base still pairs with the partner.
            new_pairs = (_WC.get(new_base) == partner_base)

            if orig_pairs and not new_pairs:
                disruptions += 1

        return disruptions

    def _codon_safety_notes(
        self,
        original: str,
        replacement: str,
        amino_acid: str,
    ) -> List[str]:
        """
        Generate safety notes about a codon substitution.

        These notes alert the user to potential non-obvious effects of
        the synonymous change.
        """
        notes: List[str] = []

        # Check usage frequency difference.
        orig_freq = HUMAN_CODON_FREQUENCY.get(original, 10.0)
        new_freq = HUMAN_CODON_FREQUENCY.get(replacement, 10.0)

        if new_freq < orig_freq * 0.5:
            notes.append(
                f"Codon usage: {replacement} ({new_freq:.1f}/1000) is less "
                f"common than {original} ({orig_freq:.1f}/1000) in human cells. "
                f"This may slightly reduce translation efficiency."
            )

        # Check for CpG creation (can be methylated and silenced).
        # CpG dinucleotides are targets for DNA methyltransferases, which
        # can lead to transcriptional silencing if the donor integrates
        # near a CpG island.
        if "CG" in replacement and "CG" not in original:
            notes.append(
                "Creates a CpG dinucleotide. CpG sites can be methylated "
                "in mammalian cells, potentially affecting expression if the "
                "donor integrates into the genome."
            )

        return notes

    def _flag_noncoding_problems(
        self,
        sequence: str,
        coding_regions: List[Tuple[int, int]],
        problems: List,
    ) -> List[str]:
        """
        Flag structural problems in non-coding regions that cannot be
        automatically optimized.

        Non-coding regions include introns, UTRs, and regulatory elements.
        Changes to these regions risk disrupting splicing, translation
        regulation, or gene expression. We report the problems but leave
        the decision to the molecular biologist.
        """
        warnings: List[str] = []

        for problem in problems:
            if isinstance(problem, G4Motif):
                if problem.strand != "sense":
                    continue
                p_start, p_end = problem.start, problem.end
                structure_type = "G-quadruplex"
            elif isinstance(problem, Hairpin):
                if not problem.is_stable:
                    continue
                p_start, p_end = problem.stem_start, problem.stem_end
                structure_type = "hairpin"
            else:
                continue

            # Check if ANY part of this problem is outside all coding regions.
            problem_positions = set(range(p_start, p_end))
            coding_positions = set()
            for cds_start, cds_end in coding_regions:
                coding_positions.update(range(cds_start, cds_end))

            noncoding_overlap = problem_positions - coding_positions
            if noncoding_overlap:
                nc_start = min(noncoding_overlap)
                nc_end = max(noncoding_overlap) + 1
                warnings.append(
                    f"{structure_type} at position {p_start}-{p_end} extends "
                    f"into non-coding region ({nc_start}-{nc_end}). "
                    f"Automatic optimization is not safe here because changes "
                    f"to non-coding regions may affect splicing signals, "
                    f"regulatory elements, or mRNA stability. Manual review "
                    f"is recommended."
                )

        return warnings
