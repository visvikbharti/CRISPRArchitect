"""
Prime Editing Feasibility Engine for CRISPRArchitect v2
========================================================

Assesses whether a variant is correctable by prime editing (PE2/PE3),
designs the pegRNA components (spacer, PBS, RT template), and searches
for a PE3 nicking guide.

Biological background
---------------------
Prime editing uses a Cas9 H840A nickase fused to an engineered M-MLV
reverse transcriptase (RT).  The prime editing guide RNA (pegRNA)
contains three functional segments:

  1. **Spacer** (20 nt): directs Cas9 to nick the PAM-containing strand
     (the non-edited strand).
  2. **Primer Binding Site (PBS)** (10-17 nt): hybridises to the nicked
     strand upstream (3') of the nick, priming reverse transcription.
  3. **RT template** (10-30+ nt): encodes the desired edit plus a short
     homology region downstream.

PE3 strategy adds a second nicking guide 40-100 bp away on the opposite
strand to stimulate mismatch repair in favour of the edited strand,
boosting efficiency 2-5x (Anzalone et al., Nature, 2019).

Size limits (empirically determined):
  - Deletions up to ~80 bp (Anzalone et al., 2019)
  - Insertions up to ~40 bp (Anzalone et al., 2019; larger inserts have
    sharply reduced efficiency)

PBS design rules (Anzalone et al., 2019; Nelson et al., Nat Biotechnol, 2022):
  - Default 13 nt; range 10-17 nt
  - GC content 40-60% preferred for PBS stability

RT template design:
  - Must encode the edit + downstream homology (flap resolution)
  - Shorter templates for substitutions (~10-15 nt total)
  - Longer templates for insertions (insert + 10 nt downstream)

Python 3.9 compatible -- no slots=True, no match statements, no X | Y unions.
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional, Tuple

from core.models import (
    FeasibilityLabel,
    GuideCandidate,
    NormalizedVariant,
    PrimeEditingFeasibility,
)
from core.feasibility.pam_scan import EnhancedPAMScanner

try:
    from utils.sequence import gc_content, reverse_complement
    from utils.constants import NUCLEASE_PARAMS
except ImportError:
    from crisprarchitect.utils.sequence import gc_content, reverse_complement
    from crisprarchitect.utils.constants import NUCLEASE_PARAMS


# ─── Size limits (Anzalone et al., Nature, 2019) ─────────────────────
MAX_DELETION_BP = 80   # PE efficiency drops sharply beyond this
MAX_INSERTION_BP = 40   # PE efficiency drops sharply beyond this

# ─── PBS design parameters ───────────────────────────────────────────
PBS_DEFAULT_LEN = 13
PBS_MIN_LEN = 10
PBS_MAX_LEN = 17

# ─── RT template parameters ─────────────────────────────────────────
RT_DOWNSTREAM_HOMOLOGY = 10  # nt of homology past the edit
RT_MIN_HOMOLOGY = 10
RT_MAX_HOMOLOGY = 15

# ─── PE3 nicking guide search window ────────────────────────────────
PE3_MIN_DISTANCE = 40   # bp from primary nick
PE3_MAX_DISTANCE = 100  # bp from primary nick

# ─── Scoring weights ────────────────────────────────────────────────
_W_DISTANCE = 0.30   # primary guide distance to edit
_W_PBS_GC = 0.20     # PBS GC quality
_W_RT_LEN = 0.15     # shorter RT is better
_W_PE3 = 0.35        # PE3 availability bonus


class PrimeEditingEngine:
    """Assess prime editing feasibility for a variant.

    Parameters
    ----------
    pam_scanner : EnhancedPAMScanner
        Pre-configured scanner (carries nuclease info).

    References
    ----------
    Anzalone et al., Nature, 2019 (prime editing, PE2/PE3 design)
    Nelson et al., Nat Biotechnol, 2022 (PrimeDesign rules)
    Chen et al., Cell, 2021 (PE4/PE5, engineered pegRNAs)
    """

    def __init__(self, pam_scanner: EnhancedPAMScanner) -> None:
        self.scanner = pam_scanner

    def check_feasibility(
        self,
        variant: NormalizedVariant,
        local_sequence: str,
        edit_pos_in_sequence: int,
        nuclease: str = "SpCas9",
    ) -> PrimeEditingFeasibility:
        """Evaluate whether *variant* can be corrected by prime editing.

        Parameters
        ----------
        variant : NormalizedVariant
            Fully annotated variant.
        local_sequence : str
            Genomic sequence window around the variant.
        edit_pos_in_sequence : int
            0-based index of the edit start within *local_sequence*.
        nuclease : str
            Nuclease to use (default ``"SpCas9"``).

        Returns
        -------
        PrimeEditingFeasibility
        """
        ref = variant.input.ref_allele.upper()
        alt = variant.input.alt_allele.upper()

        # ── Classify edit type and size ───────────────────────────────
        edit_type, edit_size = self._classify_edit(ref, alt)

        # ── Size-limit checks (Anzalone et al., 2019) ────────────────
        if edit_type == "deletion" and edit_size > MAX_DELETION_BP:
            return PrimeEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                edit_type=edit_type,
                rejection_reason=(
                    f"Deletion of {edit_size} bp exceeds the PE limit of "
                    f"{MAX_DELETION_BP} bp (Anzalone et al., 2019)."
                ),
            )
        if edit_type == "insertion" and edit_size > MAX_INSERTION_BP:
            return PrimeEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                edit_type=edit_type,
                rejection_reason=(
                    f"Insertion of {edit_size} bp exceeds the PE limit of "
                    f"{MAX_INSERTION_BP} bp (Anzalone et al., 2019)."
                ),
            )

        # ── Scan for nicking guides (PE uses Cas9 nickase) ───────────
        scanner = (
            self.scanner
            if self.scanner.nuclease == nuclease
            else EnhancedPAMScanner(nuclease)
        )
        # Search in a moderate window -- PE nick should be near the edit
        guides = scanner.scan(local_sequence, edit_pos_in_sequence, window_bp=60)

        if not guides:
            return PrimeEditingFeasibility(
                label=FeasibilityLabel.NOT_FEASIBLE,
                edit_type=edit_type,
                rejection_reason="No suitable nicking guide found near the edit.",
            )

        # ── Design pegRNA for the best primary guide ──────────────────
        best_guide = guides[0]
        pbs_len, pbs_seq = self._design_pbs(
            local_sequence, best_guide, edit_pos_in_sequence
        )
        rt_len, rt_seq = self._design_rt_template(
            local_sequence, best_guide, edit_pos_in_sequence,
            ref, alt, edit_type,
        )

        # ── Search for PE3 nicking guide ──────────────────────────────
        pe3_guide, pe3_distance = self._find_pe3_nick(
            local_sequence, edit_pos_in_sequence, best_guide, scanner,
        )

        # ── Score ─────────────────────────────────────────────────────
        score = self._score(
            best_guide, pbs_seq, rt_len, pe3_guide is not None, edit_size,
        )

        # ── Determine label ───────────────────────────────────────────
        warnings: List[str] = []
        label = FeasibilityLabel.FEASIBLE

        if pe3_guide is None:
            warnings.append(
                "No PE3 nicking guide found 40-100 bp away; "
                "PE2-only efficiency may be lower."
            )
            label = FeasibilityLabel.MARGINAL

        pbs_gc = gc_content(pbs_seq) if pbs_seq else 0.0
        if pbs_gc < 0.40 or pbs_gc > 0.60:
            warnings.append(
                f"PBS GC content ({pbs_gc:.0%}) is outside the 40-60% optimal "
                f"range (Anzalone et al., Nature, 2019)."
            )

        if edit_size > 20:
            warnings.append(
                f"Large edit ({edit_size} bp) -- PE efficiency may be reduced."
            )
            if label == FeasibilityLabel.FEASIBLE:
                label = FeasibilityLabel.MARGINAL

        return PrimeEditingFeasibility(
            label=label,
            best_guide=best_guide,
            pbs_length=pbs_len,
            rt_template_length=rt_len,
            rt_template_sequence=rt_seq,
            pe3_nick_guide=pe3_guide,
            pe3_nick_distance=pe3_distance,
            edit_type=edit_type,
            score=round(score, 4),
            warnings=warnings,
        )

    # ── internals ─────────────────────────────────────────────────────

    @staticmethod
    def _classify_edit(ref: str, alt: str) -> Tuple[str, int]:
        """Return (edit_type, size_bp)."""
        if len(ref) == 1 and len(alt) == 1:
            return ("substitution", 1)
        elif len(ref) > len(alt):
            # deletion (patient has shorter allele -- but we're correcting
            # alt->ref, so the correction is an insertion if ref is longer)
            # Actually: ref is longer means patient deleted bases.
            # To correct: we need to *insert* the missing bases back.
            # But from PE perspective the edit_type should reflect what PE does.
            # Correction direction: alt -> ref.  ref longer = insertion by PE.
            return ("insertion", len(ref) - len(alt))
        elif len(alt) > len(ref):
            # alt is longer, correction removes extra bases = deletion by PE
            return ("deletion", len(alt) - len(ref))
        else:
            # complex: same length > 1 (MNV)
            return ("substitution", len(ref))

    @staticmethod
    def _design_pbs(
        sequence: str,
        guide: GuideCandidate,
        edit_pos: int,
    ) -> Tuple[int, str]:
        """Design the Primer Binding Site (PBS).

        The PBS is complementary to the target strand immediately upstream
        (3' direction on the nicked strand) of the nick site.

        Returns (pbs_length, pbs_sequence).
        """
        nick_pos = guide.cut_position  # 0-based in sequence
        seq_upper = sequence.upper()

        # The PBS binds upstream of the nick on the PAM-containing strand.
        # For a + strand guide: nick is on the non-target (coding) strand,
        # PBS hybridises to the target strand going 5' from the nick.
        # We take bases immediately 5' (upstream) of the nick on the target strand.
        pbs_len = PBS_DEFAULT_LEN

        if guide.strand == "+":
            pbs_start = max(0, nick_pos - pbs_len)
            pbs_template = seq_upper[pbs_start:nick_pos]
            # PBS on pegRNA is the reverse complement (binds to target strand)
            pbs_seq = reverse_complement(pbs_template)
        else:
            pbs_start = nick_pos + 1
            pbs_end = min(len(seq_upper), pbs_start + pbs_len)
            pbs_template = seq_upper[pbs_start:pbs_end]
            pbs_seq = pbs_template  # already in the correct orientation

        actual_len = len(pbs_seq)
        # Clamp to allowed range
        if actual_len < PBS_MIN_LEN:
            actual_len = len(pbs_seq)  # can't extend beyond sequence
        elif actual_len > PBS_MAX_LEN:
            pbs_seq = pbs_seq[:PBS_MAX_LEN]
            actual_len = PBS_MAX_LEN

        return (actual_len, pbs_seq)

    @staticmethod
    def _design_rt_template(
        sequence: str,
        guide: GuideCandidate,
        edit_pos: int,
        ref: str,
        alt: str,
        edit_type: str,
    ) -> Tuple[int, str]:
        """Design the RT template portion of the pegRNA.

        The RT template encodes:
          1. Any sequence between the nick and the edit
          2. The desired correction (ref allele replacing alt)
          3. Downstream homology (10-15 nt)

        Returns (rt_length, rt_sequence).
        """
        seq_upper = sequence.upper()
        nick_pos = guide.cut_position

        # For + strand guide the RT template is read from the nick going 3'
        # on the target strand (which is 5'->3' in genomic coordinates from
        # the nick to downstream).
        # The template must encode: nick->edit region + corrected bases + downstream

        if guide.strand == "+":
            # Region from nick to just before edit
            pre_edit = seq_upper[nick_pos:edit_pos]
            # The corrected sequence (ref replacing alt)
            correction = ref
            # Downstream homology
            after_edit = edit_pos + len(alt)  # skip the alt allele
            downstream = seq_upper[after_edit:after_edit + RT_DOWNSTREAM_HOMOLOGY]
            rt_seq = reverse_complement(pre_edit + correction + downstream)
        else:
            # For - strand guide, directions are swapped
            after_nick = edit_pos
            pre_edit = seq_upper[after_nick:nick_pos]
            correction = ref
            before_edit = max(0, edit_pos - RT_DOWNSTREAM_HOMOLOGY)
            downstream = seq_upper[before_edit:edit_pos]
            rt_seq = downstream + correction + pre_edit

        rt_len = len(rt_seq)
        return (rt_len, rt_seq)

    @staticmethod
    def _find_pe3_nick(
        sequence: str,
        edit_pos: int,
        primary_guide: GuideCandidate,
        scanner: EnhancedPAMScanner,
    ) -> Tuple[Optional[GuideCandidate], int]:
        """Find a PE3 nicking guide 40-100 bp away on the opposite strand.

        Returns (pe3_guide, distance) or (None, 0).
        """
        # Scan the full sequence
        all_guides = scanner.scan(sequence, edit_pos, window_bp=200)

        opposite_strand = "-" if primary_guide.strand == "+" else "+"
        best_pe3 = None  # type: Optional[GuideCandidate]
        best_distance = 0

        for g in all_guides:
            if g.strand != opposite_strand:
                continue
            dist = abs(g.cut_position - primary_guide.cut_position)
            if PE3_MIN_DISTANCE <= dist <= PE3_MAX_DISTANCE:
                if best_pe3 is None or g.score > best_pe3.score:
                    best_pe3 = g
                    best_distance = dist

        return (best_pe3, best_distance)

    @staticmethod
    def _score(
        guide: GuideCandidate,
        pbs_seq: str,
        rt_len: int,
        has_pe3: bool,
        edit_size: int,
    ) -> float:
        """Composite PE feasibility score in [0, 1].

        Components:
        - Distance of primary nick to edit (closer = better)
        - PBS GC quality (40-60% optimal)
        - RT template length (shorter = better for PE efficiency)
        - PE3 availability bonus
        """
        # Distance component (reuse guide.score which already encodes distance)
        dist_score = guide.score  # already in [0, 1]

        # PBS GC component
        pbs_gc = gc_content(pbs_seq) if pbs_seq else 0.5
        gc_dev = abs(pbs_gc - 0.50)
        pbs_gc_score = max(0.0, 1.0 - (gc_dev / 0.20) ** 2)

        # RT length component (shorter is better; penalty for very long templates)
        if rt_len <= 15:
            rt_score = 1.0
        elif rt_len <= 30:
            rt_score = 1.0 - (rt_len - 15) / 30.0
        else:
            rt_score = 0.5 - (rt_len - 30) / 60.0
        rt_score = max(0.0, rt_score)

        # PE3 bonus
        pe3_score = 1.0 if has_pe3 else 0.0

        score = (
            _W_DISTANCE * dist_score
            + _W_PBS_GC * pbs_gc_score
            + _W_RT_LEN * rt_score
            + _W_PE3 * pe3_score
        )
        return score


# ═══════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    from core.models import (
        GenomicVariantInput,
        TranscriptInfo,
        TranscriptCoordinate,
        CodingAnnotation,
        ReferenceValidation,
        ExonRecord,
    )

    print("=" * 60)
    print("PrimeEditingEngine  --  self-test")
    print("=" * 60)

    def _make_variant(ref: str, alt: str) -> NormalizedVariant:
        return NormalizedVariant(
            input=GenomicVariantInput(
                chromosome="17", position=100, ref_allele=ref, alt_allele=alt,
                gene_symbol="TEST",
            ),
            transcript=TranscriptInfo(
                transcript_id="ENST00000000001", gene_symbol="TEST",
                gene_id="ENSG00000000001", chromosome="17",
                start=1, end=1000, strand=1, biotype="protein_coding",
                is_canonical=True, exons=[
                    ExonRecord("ENSE001", 1, 1, 1000, 1, "17"),
                ],
            ),
            transcript_coord=TranscriptCoordinate(
                genomic_position=100, exon_number=1, transcript_position=100,
                cds_position=100, codon_index=34, codon_position=1,
                reference_codon="ATG", reference_aa="M",
                distance_to_exon_start=99, distance_to_exon_end=900,
            ),
            coding=CodingAnnotation(consequence=ConsequenceType.MISSENSE),
            ref_validation=ReferenceValidation(
                is_valid=True, expected_ref_genomic=ref,
                provided_ref_genomic=ref, expected_ref_transcript=ref,
                provided_ref_transcript=ref,
            ),
        )

    # Build a synthetic 400-bp sequence with PAM sites
    import random
    random.seed(42)
    seq_parts = []
    for i in range(400):
        seq_parts.append(random.choice("ATCG"))
    synth = list("".join(seq_parts))
    # Break poly-T
    s = "".join(synth)
    while "TTTT" in s:
        s = s.replace("TTTT", "ATCG", 1)
    synth = list(s)

    # Place a good NGG PAM near the edit (pos 200)
    # protospacer at 185-205, PAM at 205-208
    synth[205] = "A"
    synth[206] = "G"
    synth[207] = "G"
    # Ensure decent GC in protospacer region
    for i in range(185, 205):
        synth[i] = ["A", "G", "C", "T"][(i - 185) % 4]
    test_seq = "".join(synth)

    edit_pos = 200
    scanner = EnhancedPAMScanner("SpCas9")
    engine = PrimeEditingEngine(scanner)

    from core.models import ConsequenceType

    # Test 1: SNV (substitution)
    print("\n--- Test 1: SNV substitution ---")
    v_snv = _make_variant(ref="G", alt="A")
    r = engine.check_feasibility(v_snv, test_seq, edit_pos)
    print(f"  Label: {r.label.value}")
    print(f"  Edit type: {r.edit_type}")
    print(f"  PBS length: {r.pbs_length}")
    print(f"  RT length: {r.rt_template_length}")
    print(f"  PE3 found: {r.pe3_nick_guide is not None}")
    print(f"  Score: {r.score:.4f}")
    print(f"  Warnings: {r.warnings}")

    # Test 2: Large deletion (>80 bp) -> NOT_FEASIBLE
    print("\n--- Test 2: Deletion > 80 bp ---")
    v_del = _make_variant(ref="A", alt="A" + "G" * 85)
    r_del = engine.check_feasibility(v_del, test_seq, edit_pos)
    print(f"  Label: {r_del.label.value}")
    print(f"  Rejection: {r_del.rejection_reason}")
    assert r_del.label == FeasibilityLabel.NOT_FEASIBLE, \
        "Deletion >80 bp must be NOT_FEASIBLE"

    # Test 3: Large insertion (>40 bp) -> NOT_FEASIBLE
    print("\n--- Test 3: Insertion > 40 bp ---")
    v_ins = _make_variant(ref="A" + "G" * 45, alt="A")
    r_ins = engine.check_feasibility(v_ins, test_seq, edit_pos)
    print(f"  Label: {r_ins.label.value}")
    print(f"  Rejection: {r_ins.rejection_reason}")
    assert r_ins.label == FeasibilityLabel.NOT_FEASIBLE, \
        "Insertion >40 bp must be NOT_FEASIBLE"

    # Test 4: Small indel (within limits)
    print("\n--- Test 4: Small 5-bp deletion (PE-amenable) ---")
    v_small_del = _make_variant(ref="A", alt="AGTCGT")
    r_small = engine.check_feasibility(v_small_del, test_seq, edit_pos)
    print(f"  Label: {r_small.label.value}")
    print(f"  Edit type: {r_small.edit_type}")
    print(f"  Score: {r_small.score:.4f}")

    # Test 5: Boundary exactly at 80 bp deletion -> NOT rejected by size
    # len(alt) - len(ref) = 81 - 1 = 80.  80 > 80 is False -> size check passes.
    print("\n--- Test 5: Boundary 80-bp deletion ---")
    v_boundary = _make_variant(ref="A", alt="A" + "G" * 80)
    r_boundary = engine.check_feasibility(v_boundary, test_seq, edit_pos)
    print(f"  Label: {r_boundary.label.value}")
    # Size check should NOT reject (80 is not > 80)
    assert "exceeds the PE limit" not in r_boundary.rejection_reason, \
        "80-bp deletion must NOT be rejected by size limit (80 > 80 is False)"
    print("  80-bp deletion not rejected by size limit. Correct.")

    # Test 5b: 81 bp deletion -> MUST be rejected by size
    v_over = _make_variant(ref="A", alt="A" + "G" * 81)
    r_over = engine.check_feasibility(v_over, test_seq, edit_pos)
    assert r_over.label == FeasibilityLabel.NOT_FEASIBLE, \
        "81-bp deletion must be NOT_FEASIBLE"
    print("  81-bp deletion correctly rejected.")

    print("\nAll assertions passed.  PASS")
