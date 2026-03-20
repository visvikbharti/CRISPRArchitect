"""
G-Quadruplex (G4) Motif Scanner for cssDNA Donor Templates
===========================================================

BIOLOGICAL BACKGROUND — WHY G-QUADRUPLEXES MATTER FOR HDR
----------------------------------------------------------

G-quadruplexes (G4s) are non-canonical four-stranded DNA secondary structures
that form in guanine-rich sequences. They arise when four runs of two or more
consecutive guanines (G-tracts) fold into stacked planar tetrads of four
guanines, each stabilized by Hoogsteen hydrogen bonding and coordinated by a
monovalent cation (typically K+ or Na+) in the central channel.

    G --- G
    |     |       A single G-tetrad: four guanines arranged in a planar
    G --- G       ring via Hoogsteen base pairing, with a central K+ ion.

A G-quadruplex typically consists of 3 or more stacked G-tetrads, giving rise
to EXTRAORDINARY thermodynamic stability. The melting temperature (Tm) of a G4
can exceed 90 degrees C — far above physiological temperature — meaning that
once a G4 forms, it is essentially permanent under cellular conditions.

WHY THIS IS CRITICAL FOR cssDNA-MEDIATED HDR:
----------------------------------------------

1. RAD51 LOADING BLOCKADE:
   Homology-directed repair (HDR) requires the RAD51 recombinase to polymerize
   along single-stranded DNA (ssDNA) to form a nucleoprotein filament. This
   filament then performs homology search and strand invasion. However, RAD51
   can only bind to unstructured ssDNA. If a G4 structure forms within the
   donor's homology arm, RAD51 cannot load past it, creating a gap in the
   filament that prevents productive strand invasion.

   Experimental evidence: Schiavone et al. (Nature Chemical Biology, 2014)
   showed that G4 structures block replication fork progression, and by
   analogy, they would sterically block RAD51 filament growth.

2. DONOR TRAPPING IN NON-PRODUCTIVE CONFORMATIONS:
   G4 structures in the BACKBONE (insert region) of a cssDNA donor can fold
   the circular donor into a compact topology that sequesters the homology
   arms, making them geometrically inaccessible for strand invasion even if
   they themselves are unstructured.

3. INTRAMOLECULAR vs. INTERMOLECULAR G4:
   In cssDNA donors, we are primarily concerned with INTRAMOLECULAR G4
   formation — the single strand folding back on itself. This is more likely
   than intermolecular G4 because the effective local concentration of G-tracts
   within the same molecule is very high.

CANONICAL G4 MOTIF:
   G{3+} - N{1-7} - G{3+} - N{1-7} - G{3+} - N{1-7} - G{3+}

   Where each G{3+} is a run of 3 or more consecutive guanines (forming
   the G-tetrad stacks), and N{1-7} are the loop regions of 1-7 nucleotides
   connecting the G-tracts. Shorter loops generally yield MORE stable G4s
   because they reduce the entropic cost of folding.

STABILITY SCORING RATIONALE:
   - More G-runs → more G-tetrads → more stacking → higher stability
   - Shorter loops → lower entropic penalty → higher stability
   - G4s with loops of 1 nt are the most stable (loop entropy is minimal)
   - G4s with loops of 7 nt are at the detection threshold

References:
   - Burge et al., Nucleic Acids Res, 2006 (G4 motif definition)
   - Huppert & Balasubramanian, Nucleic Acids Res, 2005 (G4 prediction)
   - Schiavone et al., Nature Chemical Biology, 2014 (G4 and recombination)
   - Biffi et al., Nature Chemistry, 2013 (G4 in living cells)
"""

from __future__ import annotations

import re
import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class G4Motif:
    """
    Represents a single predicted G-quadruplex motif found in a DNA sequence.

    The stability_score is a heuristic that combines information about the
    number and length of G-runs with the loop lengths. Higher scores indicate
    G4s that are more likely to form spontaneously and resist unfolding under
    physiological conditions.

    Attributes
    ----------
    start : int
        0-based index of the first nucleotide in the G4 motif.
    end : int
        0-based index one past the last nucleotide (Python slice convention).
    sequence : str
        The actual nucleotide sequence of the motif.
    g_run_lengths : list of int
        The length of each G-run (e.g., [3, 4, 3, 3] for a 4-tract G4).
        The number of stacked G-tetrads equals min(g_run_lengths), because
        the shortest G-run limits how many complete tetrads can form.
    loop_lengths : list of int
        The length of each loop connecting successive G-runs.
        For a canonical G4, this list has exactly 3 elements.
    stability_score : float
        Heuristic score from 0.0 (barely a G4) to 1.0 (maximally stable).
        See _compute_stability_score() for the scoring formula.
    strand : str
        'sense' if found on the input strand, 'antisense' if found on the
        reverse complement. For cssDNA donors (which are single-stranded),
        the 'sense' hits are the primary concern; 'antisense' hits indicate
        potential G4 formation if the complementary strand were generated
        during synthesis or if the donor is used as a template.
    """
    start: int
    end: int
    sequence: str
    g_run_lengths: List[int] = field(default_factory=list)
    loop_lengths: List[int] = field(default_factory=list)
    stability_score: float = 0.0
    strand: str = "sense"

    @property
    def num_tetrads(self) -> int:
        """
        The number of stacked G-tetrads this G4 can form.

        This equals the length of the shortest G-run, because each tetrad
        requires one guanine from each of the four G-runs. If the shortest
        run has 3 Gs, only 3 complete tetrads can form (the remaining Gs
        in longer runs may participate in alternative conformations or
        remain unpaired).
        """
        if not self.g_run_lengths:
            return 0
        return min(self.g_run_lengths)

    @property
    def length(self) -> int:
        """Total length of the G4 motif in nucleotides."""
        return self.end - self.start

    def __repr__(self) -> str:
        return (
            f"G4Motif(start={self.start}, end={self.end}, "
            f"tetrads={self.num_tetrads}, stability={self.stability_score:.2f}, "
            f"strand='{self.strand}')"
        )


# ---------------------------------------------------------------------------
# G-Quadruplex Scanner
# ---------------------------------------------------------------------------

class GQuadruplexScanner:
    """
    Scans a DNA sequence for potential G-quadruplex (G4) forming motifs.

    This scanner uses a regex-based approach to identify sequences matching
    the canonical G4 motif pattern:

        G{3,} N{1,7} G{3,} N{1,7} G{3,} N{1,7} G{3,}

    The scanner also examines the reverse complement to identify C-rich
    sequences on the input strand that could form G4s on the complementary
    strand (relevant when the cssDNA donor serves as a synthesis template).

    The biological significance of each hit is assessed via a stability score
    that considers:
      - Number of G-tetrads (more tetrads = more stacking = more stable)
      - Loop lengths (shorter loops = lower entropy cost = more stable)
      - Total G-run length (longer G-runs allow more tetrad conformations)

    Usage Example
    -------------
    >>> scanner = GQuadruplexScanner()
    >>> motifs = scanner.find_g4_motifs("ATGGGTAGGGTTAGGGTTAGGGTTT")
    >>> for m in motifs:
    ...     print(m)

    Parameters
    ----------
    min_g_run : int
        Minimum number of consecutive Gs in each G-tract (default: 3).
        The canonical G4 requires at least 3 Gs per tract to form a
        three-tetrad stack, which is the minimum for appreciable stability.
    max_loop : int
        Maximum loop length between G-runs (default: 7).
        Loops longer than 7 nt impose a large entropic penalty, making G4
        formation thermodynamically unfavorable. Some studies use 12, but
        7 is the widely accepted cutoff for high-confidence predictions.
    min_loop : int
        Minimum loop length (default: 1).
        A loop of 0 nt would mean adjacent G-runs merge into one long run,
        which is handled separately. The minimum physical loop is 1 nt.
    """

    # -------------------------------------------------------------------
    # The canonical G4 regex pattern.
    #
    # Explanation of the pattern:
    #   (G{3,})     -> First G-run: 3 or more consecutive guanines
    #   ([ATCG]{1,7}) -> First loop: 1-7 nucleotides of any base
    #   (G{3,})     -> Second G-run
    #   ([ATCG]{1,7}) -> Second loop
    #   (G{3,})     -> Third G-run
    #   ([ATCG]{1,7}) -> Third loop
    #   (G{3,})     -> Fourth G-run
    #
    # Capturing groups allow us to extract G-run and loop lengths.
    # -------------------------------------------------------------------

    def __init__(
        self,
        min_g_run: int = 3,
        max_loop: int = 7,
        min_loop: int = 1,
    ):
        self.min_g_run = min_g_run
        self.max_loop = max_loop
        self.min_loop = min_loop

        # Build the regex pattern dynamically based on parameters.
        # Each G-run is captured, and each loop is captured, so we can
        # extract their lengths for stability scoring.
        g_run = f"(G{{{self.min_g_run},}})"
        loop = f"([ATCG]{{{self.min_loop},{self.max_loop}}})"
        self._pattern_str = f"{g_run}{loop}{g_run}{loop}{g_run}{loop}{g_run}"
        self._pattern = re.compile(self._pattern_str, re.IGNORECASE)

    # -----------------------------------------------------------------
    # Reverse complement utility
    # -----------------------------------------------------------------

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """
        Return the reverse complement of a DNA sequence.

        This is needed because a C-rich sequence on the input strand
        corresponds to a G-rich sequence on the complementary strand.
        For cssDNA donors, the input strand IS the donor; but during HDR
        synthesis, a complementary strand is generated, and G4 motifs on
        THAT strand could stall the extending polymerase.

        Watson-Crick complement:  A <-> T,  G <-> C
        """
        complement_map = str.maketrans("ATCGatcg", "TAGCtagc")
        return sequence.translate(complement_map)[::-1]

    # -----------------------------------------------------------------
    # Core motif finding
    # -----------------------------------------------------------------

    def find_g4_motifs(self, sequence: str) -> List[G4Motif]:
        """
        Scan a DNA sequence for all G-quadruplex forming motifs.

        This method searches both the sense strand (the input sequence)
        and the antisense strand (reverse complement). For cssDNA donors:

          - SENSE hits are the most critical: they represent G4 structures
            that can form directly on the donor molecule as delivered.
          - ANTISENSE hits are secondary: they indicate G4 potential on the
            complementary strand that would be synthesized during HDR
            (relevant if the newly synthesized strand is displaced as ssDNA
            during SDSA, where it could form G4 and impede second-end capture).

        The method uses overlapping scanning to catch nested or overlapping
        G4 motifs. In G-rich regions, multiple G4 conformations may compete,
        and we want to report all of them.

        Parameters
        ----------
        sequence : str
            The DNA sequence to scan. Should contain only A, T, C, G characters.
            Case-insensitive.

        Returns
        -------
        list of G4Motif
            All G4 motifs found, sorted by start position. Each motif includes
            stability scoring and strand annotation.
        """
        sequence = sequence.upper().strip()
        motifs: List[G4Motif] = []

        # --- Scan the sense strand ---
        # We use finditer with overlapping by advancing the search start
        # position after each match. re.finditer normally returns non-
        # overlapping matches, so we use a manual loop.
        motifs.extend(self._scan_strand(sequence, strand="sense"))

        # --- Scan the antisense strand ---
        # The reverse complement of the input sequence. Hits here correspond
        # to C-rich regions on the input strand. We convert positions back
        # to the sense strand coordinate system for consistent reporting.
        rc_sequence = self._reverse_complement(sequence)
        rc_motifs = self._scan_strand(rc_sequence, strand="antisense")

        # Convert antisense positions to sense strand coordinates.
        # If the motif spans [start, end) on the reverse complement,
        # it corresponds to [len - end, len - start) on the sense strand.
        seq_len = len(sequence)
        for motif in rc_motifs:
            original_start = seq_len - motif.end
            original_end = seq_len - motif.start
            motif.start = original_start
            motif.end = original_end
        motifs.extend(rc_motifs)

        # Sort all motifs by start position for easy visualization.
        motifs.sort(key=lambda m: (m.start, m.strand))
        return motifs

    def _scan_strand(self, sequence: str, strand: str = "sense") -> List[G4Motif]:
        """
        Internal method: scan a single strand for G4 motifs with overlapping matches.

        Standard regex finditer returns non-overlapping matches, which can miss
        G4 motifs that share G-runs with adjacent motifs. To catch these, we
        advance the search position by 1 after each match rather than jumping
        to the end of the match.

        Parameters
        ----------
        sequence : str
            The strand sequence to scan (already uppercased).
        strand : str
            Label for the motif ('sense' or 'antisense').

        Returns
        -------
        list of G4Motif
            G4 motifs found on this strand.
        """
        motifs: List[G4Motif] = []
        # Track already-found motifs by (start, end) to avoid duplicates
        # when the overlapping search finds the same motif again.
        seen: set = set()
        pos = 0

        while pos < len(sequence):
            match = self._pattern.search(sequence, pos)
            if match is None:
                break  # No more motifs downstream

            start = match.start()
            end = match.end()
            key = (start, end, strand)

            if key not in seen:
                seen.add(key)

                # Extract the G-run lengths and loop lengths from
                # the capturing groups.
                # Groups: (G-run1)(Loop1)(G-run2)(Loop2)(G-run3)(Loop3)(G-run4)
                # Indices:    1       2       3       4       5       6       7
                g_run_lengths = [
                    len(match.group(1)),
                    len(match.group(3)),
                    len(match.group(5)),
                    len(match.group(7)),
                ]
                loop_lengths = [
                    len(match.group(2)),
                    len(match.group(4)),
                    len(match.group(6)),
                ]

                motif_seq = match.group(0)

                # Compute stability score based on biophysical heuristics.
                stability = self._compute_stability_score(
                    g_run_lengths, loop_lengths
                )

                motifs.append(G4Motif(
                    start=start,
                    end=end,
                    sequence=motif_seq,
                    g_run_lengths=g_run_lengths,
                    loop_lengths=loop_lengths,
                    stability_score=stability,
                    strand=strand,
                ))

            # Advance by 1 to allow overlapping matches.
            # This is necessary because in G-rich regions like
            # GGGAGGGAGGGAGGGAGGGAGGG, multiple valid G4 conformations
            # can start at different positions.
            pos = start + 1

        return motifs

    def _compute_stability_score(
        self,
        g_run_lengths: List[int],
        loop_lengths: List[int],
    ) -> float:
        """
        Compute a heuristic stability score for a G4 motif (0.0 to 1.0).

        The score approximates the thermodynamic stability of the G4 based on
        two biophysical principles:

        1. G-TETRAD STACKING ENERGY:
           Each stacked G-tetrad contributes approximately -20 to -25 kcal/mol
           of stacking free energy (Hazel et al., JACS, 2004). More tetrads =
           more stable. A 3-tetrad G4 (the minimum) has dG ~ -60 kcal/mol,
           while a 4-tetrad G4 has dG ~ -80 kcal/mol.

           Score component: we use min(g_run_lengths) as the number of tetrads,
           and scale from 3 (baseline) to 6+ (maximum commonly observed).

        2. LOOP ENTROPY PENALTY:
           Longer loops impose a larger entropic cost for G4 formation, following
           approximately dS ~ -R * ln(loop_length). Short loops (1-2 nt) have
           minimal entropy penalty and allow very stable G4s. Loops of 7 nt are
           at the upper limit of favorable G4 formation.

           Score component: average loop length, scaled so that loop=1 gives
           maximum score and loop=7 gives minimum score.

        The two components are combined with weights reflecting their relative
        contribution to overall G4 stability.

        Parameters
        ----------
        g_run_lengths : list of int
            Length of each G-run.
        loop_lengths : list of int
            Length of each connecting loop.

        Returns
        -------
        float
            Stability score between 0.0 and 1.0.
        """
        # --- Tetrad component ---
        # min(g_run_lengths) = number of G-tetrads.
        # 3 tetrads → score ~0.5 (baseline, moderately stable)
        # 4 tetrads → score ~0.75 (very stable)
        # 5+ tetrads → score ~1.0 (extremely stable)
        num_tetrads = min(g_run_lengths)
        # Normalize: 3 tetrads = 0.4, 4 = 0.7, 5 = 0.9, 6+ = 1.0
        tetrad_score = min(1.0, (num_tetrads - 2) / 4.0)

        # --- Loop penalty component ---
        # Average loop length: shorter is better.
        # avg_loop = 1.0 → loop_score = 1.0 (best)
        # avg_loop = 7.0 → loop_score = 0.0 (worst, barely a G4)
        avg_loop = sum(loop_lengths) / max(len(loop_lengths), 1)
        # Linear scaling from 1 (score=1) to 7 (score=0)
        loop_score = max(0.0, 1.0 - (avg_loop - 1.0) / 6.0)

        # --- Combined score ---
        # Weight tetrad component more heavily because the number of tetrads
        # is the dominant determinant of G4 thermal stability.
        # Tetrad weight: 0.6, Loop weight: 0.4
        combined = 0.6 * tetrad_score + 0.4 * loop_score

        # Clamp to [0, 1]
        return max(0.0, min(1.0, combined))

    # -----------------------------------------------------------------
    # Risk scoring
    # -----------------------------------------------------------------

    def score_g4_risk(self, sequence: str) -> float:
        """
        Compute an overall G4 risk score for a cssDNA donor sequence.

        This single number summarizes the likelihood that G4 structures in the
        donor will impair HDR. The score ranges from 0.0 (no G4 risk — the
        sequence is devoid of G-rich motifs) to 1.0 (extremely high risk —
        the sequence contains multiple highly stable G4 motifs likely to form
        spontaneously and block RAD51 filament assembly).

        The risk score accounts for:
          - The number of G4 motifs (more motifs = more chance of at least one forming)
          - The stability of the most stable motif (the worst-case scenario)
          - The fraction of the sequence covered by G4 motifs (extensive G4-prone
            regions are worse than isolated motifs)

        Interpretation guide for the end user:
          0.0 - 0.2  : LOW RISK    — No significant G4 formation expected.
          0.2 - 0.5  : MODERATE    — G4 possible; monitor if in homology arms.
          0.5 - 0.8  : HIGH RISK   — G4 likely; consider sequence optimization.
          0.8 - 1.0  : VERY HIGH   — Strong G4 predicted; redesign donor.

        Parameters
        ----------
        sequence : str
            The donor DNA sequence.

        Returns
        -------
        float
            G4 risk score between 0.0 and 1.0.
        """
        motifs = self.find_g4_motifs(sequence)

        if not motifs:
            return 0.0

        # --- Component 1: Maximum stability of any single motif ---
        # The most stable G4 is the one most likely to form and persist.
        # A single very stable G4 in a homology arm is enough to impair HDR.
        max_stability = max(m.stability_score for m in motifs)

        # --- Component 2: Coverage fraction ---
        # What fraction of the sequence is covered by G4 motifs?
        # Higher coverage means more of the donor is at risk.
        # We merge overlapping intervals to avoid double-counting.
        covered_positions = set()
        for m in motifs:
            covered_positions.update(range(m.start, m.end))
        coverage_fraction = len(covered_positions) / max(len(sequence), 1)

        # --- Component 3: Motif count penalty ---
        # Multiple G4 motifs compound the risk because even if one doesn't
        # form, another might. Saturates at ~5 motifs.
        count_factor = min(1.0, len(motifs) / 5.0)

        # --- Combine ---
        # Max stability is the strongest predictor; coverage and count
        # are secondary modifiers.
        risk = (
            0.50 * max_stability
            + 0.30 * coverage_fraction
            + 0.20 * count_factor
        )
        return max(0.0, min(1.0, risk))

    # -----------------------------------------------------------------
    # Homology arm mapping
    # -----------------------------------------------------------------

    def map_g4_positions(
        self,
        sequence: str,
        homology_arms: Tuple[Tuple[int, int], Tuple[int, int]],
    ) -> Dict:
        """
        Map G4 motifs to homology arm regions of a cssDNA donor.

        This is the most important method for donor design assessment. A G4
        anywhere in the donor is suboptimal, but a G4 WITHIN a homology arm
        is CRITICAL because:

          1. The homology arm must be single-stranded and accessible for RAD51
             filament formation. A G4 within the arm directly blocks this.
          2. The insert region (between the arms) does not participate in
             strand invasion — G4s there are less harmful (though they can
             still affect donor topology).

        This method returns a structured report showing:
          - Which G4 motifs overlap with the LEFT homology arm
          - Which G4 motifs overlap with the RIGHT homology arm
          - Which G4 motifs are in the insert region (lower priority)
          - A criticality flag for motifs in the arms

        Parameters
        ----------
        sequence : str
            The full donor DNA sequence.
        homology_arms : tuple of two (start, end) tuples
            ((left_arm_start, left_arm_end), (right_arm_start, right_arm_end))
            where start/end are 0-based indices in the sequence.
            The LEFT arm is the homology region 5' of the insert.
            The RIGHT arm is the homology region 3' of the insert.

        Returns
        -------
        dict
            {
                "left_arm_motifs": [G4Motif, ...],
                "right_arm_motifs": [G4Motif, ...],
                "insert_motifs": [G4Motif, ...],
                "left_arm_risk": float,  # 0-1, G4 risk in left arm
                "right_arm_risk": float, # 0-1, G4 risk in right arm
                "critical_motifs": [G4Motif, ...],  # Motifs in arms only
                "summary": str,
            }
        """
        left_arm = homology_arms[0]
        right_arm = homology_arms[1]

        # Find all G4 motifs in the full sequence.
        all_motifs = self.find_g4_motifs(sequence)

        left_arm_motifs: List[G4Motif] = []
        right_arm_motifs: List[G4Motif] = []
        insert_motifs: List[G4Motif] = []

        for motif in all_motifs:
            # Check overlap with left arm.
            # Two intervals [a, b) and [c, d) overlap if a < d and c < b.
            overlaps_left = motif.start < left_arm[1] and left_arm[0] < motif.end
            overlaps_right = motif.start < right_arm[1] and right_arm[0] < motif.end

            if overlaps_left:
                left_arm_motifs.append(motif)
            if overlaps_right:
                right_arm_motifs.append(motif)
            if not overlaps_left and not overlaps_right:
                insert_motifs.append(motif)

        # Compute per-arm risk scores.
        # These are based on the most stable G4 in each arm, adjusted
        # for coverage of the arm region.
        left_arm_risk = self._arm_risk(left_arm_motifs, left_arm)
        right_arm_risk = self._arm_risk(right_arm_motifs, right_arm)

        # Critical motifs are those that overlap with either homology arm.
        critical_motifs = left_arm_motifs + right_arm_motifs

        # Generate human-readable summary.
        summary = self._generate_summary(
            left_arm_motifs, right_arm_motifs, insert_motifs,
            left_arm_risk, right_arm_risk,
        )

        return {
            "left_arm_motifs": left_arm_motifs,
            "right_arm_motifs": right_arm_motifs,
            "insert_motifs": insert_motifs,
            "left_arm_risk": left_arm_risk,
            "right_arm_risk": right_arm_risk,
            "critical_motifs": critical_motifs,
            "summary": summary,
        }

    def _arm_risk(
        self, motifs: List[G4Motif], arm_region: Tuple[int, int]
    ) -> float:
        """
        Compute G4 risk score for a specific homology arm.

        The risk depends on:
          - Maximum stability of any G4 motif overlapping this arm
          - The fraction of the arm covered by G4 motifs

        A single high-stability G4 in the arm is enough to block RAD51
        loading in that region. The coverage fraction indicates whether
        the blockade is localized (a few bases) or extensive.

        Parameters
        ----------
        motifs : list of G4Motif
            G4 motifs overlapping with this arm.
        arm_region : tuple of (int, int)
            (start, end) of the arm in the sequence.

        Returns
        -------
        float
            Risk score from 0.0 (no G4 in arm) to 1.0 (severe blockade).
        """
        if not motifs:
            return 0.0

        arm_start, arm_end = arm_region
        arm_length = arm_end - arm_start

        # Maximum stability of G4 motifs in this arm.
        max_stab = max(m.stability_score for m in motifs)

        # Coverage: fraction of arm positions covered by G4 motifs.
        covered = set()
        for m in motifs:
            # Only count positions within the arm bounds.
            overlap_start = max(m.start, arm_start)
            overlap_end = min(m.end, arm_end)
            if overlap_start < overlap_end:
                covered.update(range(overlap_start, overlap_end))
        coverage = len(covered) / max(arm_length, 1)

        # Combine: max stability is dominant, coverage is a modifier.
        return min(1.0, 0.7 * max_stab + 0.3 * coverage)

    def _generate_summary(
        self,
        left_motifs: List[G4Motif],
        right_motifs: List[G4Motif],
        insert_motifs: List[G4Motif],
        left_risk: float,
        right_risk: float,
    ) -> str:
        """Generate a human-readable summary of G4 findings."""
        lines: List[str] = ["G-Quadruplex Analysis Summary"]
        lines.append("=" * 40)

        total = len(left_motifs) + len(right_motifs) + len(insert_motifs)
        lines.append(f"Total G4 motifs found: {total}")

        if total == 0:
            lines.append("No G-quadruplex motifs detected. This is favorable for HDR.")
            return "\n".join(lines)

        lines.append(f"  Left homology arm:  {len(left_motifs)} motif(s), risk={left_risk:.2f}")
        lines.append(f"  Right homology arm: {len(right_motifs)} motif(s), risk={right_risk:.2f}")
        lines.append(f"  Insert region:      {len(insert_motifs)} motif(s)")

        # Warnings for critical motifs.
        if left_risk > 0.5 or right_risk > 0.5:
            lines.append("")
            lines.append("WARNING: High-risk G4 motifs detected in homology arm(s)!")
            lines.append("These G4 structures are predicted to block RAD51 loading,")
            lines.append("preventing productive strand invasion during HDR.")
            lines.append("Consider synonymous codon optimization to disrupt G-runs.")

        return "\n".join(lines)
