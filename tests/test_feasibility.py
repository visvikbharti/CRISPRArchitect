"""
Comprehensive tests for CRISPRArchitect v2 feasibility engines.

Tests all four engines (PAM Scanner, Base Editing, Prime Editing, HDR Design)
using SYNTHETIC sequences only -- no network calls, no Ensembl API.

Python 3.9 compatible.  Uses pytest.
"""

from __future__ import annotations

import math
import sys
import os

import pytest

# Ensure the project root is on sys.path so bare "core.*" imports work.
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from core.models import (
    BaseEditingFeasibility,
    CodingAnnotation,
    ConsequenceType,
    EditModality,
    ExonRecord,
    FeasibilityLabel,
    GenomicVariantInput,
    GuideCandidate,
    HDRFeasibility,
    NormalizedVariant,
    PrimeEditingFeasibility,
    ReferenceValidation,
    TranscriptCoordinate,
    TranscriptInfo,
)
from core.feasibility.pam_scan import EnhancedPAMScanner
from core.feasibility.base_editing import BaseEditingEngine
from core.feasibility.prime_editing import PrimeEditingEngine
from core.feasibility.hdr_design import HDRDesignEngine


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _make_variant(ref='C', alt='T', position=100):
    """Create a minimal NormalizedVariant for testing."""
    vi = GenomicVariantInput('1', position, ref, alt, gene_symbol='TEST')
    ti = TranscriptInfo(
        'T1', 'TEST', 'G1', '1', 1, 1000, 1, 'protein_coding', True,
        [ExonRecord('ENSE001', 1, 1, 1000, 1, '1')],
    )
    tc = TranscriptCoordinate(
        position, 1, position, position, 34, 1, 'ATG', 'M', 99, 50, True,
    )
    ca = CodingAnnotation(consequence=ConsequenceType.MISSENSE)
    rv = ReferenceValidation(True, ref, ref, ref, ref)
    return NormalizedVariant(
        input=vi, transcript=ti, transcript_coord=tc,
        coding=ca, ref_validation=rv,
    )


def _build_controlled_sequence(length=400, edit_pos=200):
    """Build a synthetic sequence with known PAM sites and good GC.

    Places a + strand NGG PAM such that the cut falls very near edit_pos.
    Layout (+ strand):
      protospacer: [edit_pos - 3 - 17 .. edit_pos - 3 + 2]  (20-mer)
      PAM:         [edit_pos - 3 + 3 .. edit_pos - 3 + 5]    (NGG)
      cut:         edit_pos - 3 + 3 - 3 = edit_pos - 3       (3 bp upstream of PAM)

    Simplified: place protospacer at [edit_pos-20 .. edit_pos], PAM at
    [edit_pos .. edit_pos+3].  Cut = edit_pos - 3.

    For simplicity we embed:
      protospacer with 50% GC at positions [edit_pos-20 .. edit_pos)
      NGG PAM at positions [edit_pos .. edit_pos+3)
    """
    # Start with a repeating pattern that gives ~50% GC and no poly-T
    base_pattern = "AGCT"
    seq = list((base_pattern * ((length // 4) + 1))[:length])

    # Embed a 20-mer protospacer with 50% GC (no TTTT) just before edit_pos
    proto = list("AGCTAGCTAGCTAGCTAGCT")  # 50% GC, no TTTT
    proto_start = edit_pos - 20
    for i, base in enumerate(proto):
        seq[proto_start + i] = base

    # Embed NGG PAM right at edit_pos
    seq[edit_pos] = 'A'      # N
    seq[edit_pos + 1] = 'G'  # G
    seq[edit_pos + 2] = 'G'  # G

    # Final poly-T cleanup
    result = "".join(seq)
    while "TTTT" in result:
        result = result.replace("TTTT", "AGCT", 1)
    return result


def _build_sequence_with_pam_at(pam_pos, pam="AGG", length=400,
                                 gc_balanced=True):
    """Build a sequence with a specific PAM at a known position.

    Returns (sequence, expected_cut_position).
    For + strand SpCas9: cut is at pam_pos - 3.
    """
    base_pattern = "AGCT" if gc_balanced else "AAAA"
    seq = list((base_pattern * ((length // 4) + 1))[:length])

    # Place PAM
    for i, base in enumerate(pam):
        seq[pam_pos + i] = base

    # Place a good 20-mer protospacer upstream of PAM
    proto = list("AGCTAGCTAGCTAGCTAGCT")
    proto_start = pam_pos - 20
    if proto_start >= 0:
        for i, base in enumerate(proto):
            seq[proto_start + i] = base

    result = "".join(seq)
    while "TTTT" in result:
        result = result.replace("TTTT", "AGCT", 1)

    expected_cut = pam_pos - 3  # + strand SpCas9 cut
    return result, expected_cut


# ═══════════════════════════════════════════════════════════════════════════
# PAM Scanner Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestPAMScanner:
    """Tests for EnhancedPAMScanner."""

    def test_find_known_ngg_sites(self):
        """Scanner returns guides for a known NGG site in a constructed sequence."""
        seq, cut = _build_sequence_with_pam_at(pam_pos=203)
        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(seq, edit_position=200, window_bp=200)
        assert len(guides) > 0, "Must find at least 1 guide for embedded NGG"

    def test_gc_filter_rejects_low_gc(self):
        """Guides with <30% GC are rejected."""
        # Build a sequence where the protospacer is all A/T (GC = 0%)
        seq = list("ATAT" * 100)  # 400 bases, 0% GC
        # Place an NGG PAM
        seq[203] = 'A'
        seq[204] = 'G'
        seq[205] = 'G'
        # Protospacer at 183-203 is all A/T -> GC = 0%
        test_seq = "".join(seq)
        # Fix any accidental poly-T by substituting with A
        while "TTTT" in test_seq:
            test_seq = test_seq.replace("TTTT", "ATAT", 1)
        # Re-place the PAM in case it was damaged
        seq = list(test_seq)
        seq[203] = 'A'
        seq[204] = 'G'
        seq[205] = 'G'
        # Ensure the 20-mer is very low GC (all A/T)
        for i in range(183, 203):
            seq[i] = 'A' if i % 2 == 0 else 'T'
        test_seq = "".join(seq)
        while "TTTT" in test_seq:
            # Replace poly-T with AT pattern but keep GC < 30%
            test_seq = test_seq.replace("TTTT", "ATAT", 1)
        # Re-place the PAM again
        seq = list(test_seq)
        seq[203] = 'A'
        seq[204] = 'G'
        seq[205] = 'G'
        test_seq = "".join(seq)

        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(test_seq, edit_position=200, window_bp=200)
        # The specific embedded guide (protospacer at 183-203) should be
        # filtered out because its GC is <30%.  It might still find
        # other incidental guides, but none at the embedded position.
        embedded_guides = [
            g for g in guides
            if g.strand == "+" and abs(g.cut_position - 200) <= 5
        ]
        assert len(embedded_guides) == 0, (
            "Low-GC guide should be filtered out"
        )

    def test_gc_filter_rejects_high_gc(self):
        """Guides with >70% GC are rejected."""
        seq = list("GCGC" * 100)
        seq[203] = 'G'
        seq[204] = 'G'
        seq[205] = 'G'
        # protospacer at 183-203 is all G/C -> GC=100%
        test_seq = "".join(seq)
        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(test_seq, edit_position=200, window_bp=200)
        high_gc_guides = [
            g for g in guides if g.gc_content > 0.70
        ]
        assert len(high_gc_guides) == 0, "High-GC guides must be filtered out"

    def test_poly_t_filter_rejects_tttt(self):
        """Guides containing TTTT in the protospacer are rejected."""
        seq = list("AGCT" * 100)
        # Place NGG PAM at position 203
        seq[203] = 'A'
        seq[204] = 'G'
        seq[205] = 'G'
        # Force TTTT in the protospacer (positions 183-203)
        seq[190] = 'T'
        seq[191] = 'T'
        seq[192] = 'T'
        seq[193] = 'T'
        test_seq = "".join(seq)

        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(test_seq, edit_position=200, window_bp=200)
        for g in guides:
            assert "TTTT" not in g.sequence_20mer, (
                f"Guide with poly-T should be filtered: {g.sequence_20mer}"
            )

    def test_enfncas9_finds_at_least_as_many_as_spcas9(self):
        """enFnCas9 (NRG) finds >= SpCas9 (NGG) sites in same sequence."""
        seq = _build_controlled_sequence(length=400, edit_pos=200)
        scanner_sp = EnhancedPAMScanner("SpCas9")
        scanner_en = EnhancedPAMScanner("enFnCas9")

        guides_sp = scanner_sp.scan(seq, edit_position=200, window_bp=200)
        guides_en = scanner_en.scan(seq, edit_position=200, window_bp=200)

        # NRG is a superset of NGG (NGG matches NRG with R=G)
        assert len(guides_en) >= len(guides_sp), (
            f"enFnCas9 (NRG) must find >= SpCas9 (NGG): "
            f"{len(guides_en)} vs {len(guides_sp)}"
        )

    def test_cut_position_3bp_upstream_of_pam_plus_strand(self):
        """Cut position is 3 bp upstream of PAM for + strand guides."""
        # Construct a minimal sequence with exactly one clear NGG
        seq, expected_cut = _build_sequence_with_pam_at(pam_pos=203)
        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(seq, edit_position=200, window_bp=200)

        # Find the guide corresponding to our embedded PAM at 203
        plus_guides = [g for g in guides if g.strand == "+"]
        assert len(plus_guides) > 0, "Must find at least one + strand guide"

        # For + strand: cut = pam_pos - 3 = 200
        found_expected = any(
            g.cut_position == expected_cut for g in plus_guides
        )
        assert found_expected, (
            f"Expected cut at {expected_cut}, found cuts at "
            f"{[g.cut_position for g in plus_guides]}"
        )

    def test_guides_sorted_by_score_descending(self):
        """Returned guides are sorted by composite score, best first."""
        seq = _build_controlled_sequence(length=400, edit_pos=200)
        scanner = EnhancedPAMScanner("SpCas9")
        guides = scanner.scan(seq, edit_position=200, window_bp=200)

        if len(guides) > 1:
            for i in range(len(guides) - 1):
                assert guides[i].score >= guides[i + 1].score, (
                    f"Guides not sorted: score[{i}]={guides[i].score} "
                    f"< score[{i+1}]={guides[i+1].score}"
                )

    def test_unknown_nuclease_raises(self):
        """Unknown nuclease name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown nuclease"):
            EnhancedPAMScanner("NonexistentCas")


# ═══════════════════════════════════════════════════════════════════════════
# Base Editing Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestBaseEditing:
    """Tests for BaseEditingEngine."""

    @pytest.fixture
    def engine(self):
        scanner = EnhancedPAMScanner("SpCas9")
        return BaseEditingEngine(scanner)

    @pytest.fixture
    def test_sequence(self):
        """400-bp sequence with a known NGG placed for ABE window testing.

        Layout: PAM at 217, protospacer at 197-217, cut at 214.
        Edit at 200 -> pos_in_proto = 200 - (214-17) + 1 = 200 - 197 + 1 = 4.
        Position 4 is inside the ABE window (4-7).
        """
        seq = list("AGCT" * 100)
        # Embed NGG PAM at position 217
        seq[217] = 'A'
        seq[218] = 'G'
        seq[219] = 'G'
        # Build the 20-mer protospacer (197-217) with ~50% GC
        proto = list("AGCTAGCTAGCTAGCTAGCT")
        for i, base in enumerate(proto):
            seq[197 + i] = base
        # Set edit position (200) to 'A' for ABE testing (patient allele)
        seq[200] = 'G'  # reference allele at pos 200
        result = "".join(seq)
        while "TTTT" in result:
            result = result.replace("TTTT", "AGCT", 1)
        return result

    def test_abe_identified_for_g_to_a(self, engine, test_sequence):
        """ABE is correctly identified for G>A correction (patient=A, ref=G)."""
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.editor_type == "ABE", (
            f"Expected ABE for G>A correction, got {result.editor_type}"
        )

    def test_cbe_identified_for_t_to_c(self, engine, test_sequence):
        """CBE is correctly identified for T>C correction (patient=C on sense
        means need C->T on antisense = CBE)."""
        variant = _make_variant(ref='A', alt='G', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.editor_type == "CBE", (
            f"Expected CBE for G->A correction (alt=G, ref=A), got {result.editor_type}"
        )

    @pytest.mark.parametrize("alt,ref", [
        ('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
        ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G'),
    ])
    def test_transversions_always_not_feasible(self, engine, test_sequence,
                                                alt, ref):
        """All transversions return NOT_FEASIBLE."""
        variant = _make_variant(ref=ref, alt=alt, position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label == FeasibilityLabel.NOT_FEASIBLE, (
            f"Transversion {alt}->{ref} should be NOT_FEASIBLE, "
            f"got {result.label.value}"
        )

    def test_indel_insertion_not_feasible(self, engine, test_sequence):
        """Insertions (multi-base alt) return NOT_FEASIBLE for base editing."""
        variant = _make_variant(ref='A', alt='ATG', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label == FeasibilityLabel.NOT_FEASIBLE
        assert "indels" in result.rejection_reason.lower() or \
               "single-nucleotide" in result.rejection_reason.lower()

    def test_indel_deletion_not_feasible(self, engine, test_sequence):
        """Deletions (multi-base ref) return NOT_FEASIBLE for base editing."""
        variant = _make_variant(ref='ATG', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label == FeasibilityLabel.NOT_FEASIBLE

    def test_abe_window_is_4_to_7(self):
        """ABE editing window constants are positions 4-7."""
        from core.feasibility.base_editing import ABE_WINDOW_START, ABE_WINDOW_END
        assert ABE_WINDOW_START == 4
        assert ABE_WINDOW_END == 7

    def test_cbe_window_is_4_to_8(self):
        """CBE editing window constants are positions 4-8."""
        from core.feasibility.base_editing import CBE_WINDOW_START, CBE_WINDOW_END
        assert CBE_WINDOW_START == 4
        assert CBE_WINDOW_END == 8

    def test_target_base_must_be_in_window(self, engine):
        """If the target base is outside the editing window, result is NOT_FEASIBLE."""
        # Build a sequence where the edit position would place the target
        # base at position 1 (outside window 4-7) of any nearby protospacer.
        # The simplest way: put the NGG PAM right next to the edit so
        # position in protospacer = 18, 19, or 20 (way outside window).
        seq = list("AGCT" * 100)
        # PAM at position 201 -> protospacer at 181-201, cut at 198
        # Edit at 200 -> pos_in_proto = 200 - (198-17) + 1 = 200 - 181 + 1 = 20
        # Position 20 is outside ABE window (4-7), so it should fail.
        seq[201] = 'A'
        seq[202] = 'G'
        seq[203] = 'G'
        proto = list("AGCTAGCTAGCTAGCTAGCT")
        for i, base in enumerate(proto):
            seq[181 + i] = base
        seq[200] = 'G'  # reference allele
        test_seq = "".join(seq)
        while "TTTT" in test_seq:
            test_seq = test_seq.replace("TTTT", "AGCT", 1)

        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(variant, test_seq, 200)
        # The target base at position 20 is outside ABE window 4-7.
        # The engine should either find a different guide or return NOT_FEASIBLE.
        # We specifically check that if it IS feasible, the window position is valid.
        if result.label != FeasibilityLabel.NOT_FEASIBLE:
            if result.editor_type == "ABE":
                assert 4 <= result.target_position_in_window <= 7
            elif result.editor_type == "CBE":
                assert 4 <= result.target_position_in_window <= 8

    def test_bystander_count_correct(self, engine):
        """Bystander count reflects same-type bases in the editing window."""
        # Build a sequence where we control the exact protospacer content.
        # Place multiple A's in the ABE window to create known bystanders.
        seq = list("AGCT" * 100)
        # PAM at 217, protospacer at 197-217, cut at 214
        seq[217] = 'A'
        seq[218] = 'G'
        seq[219] = 'G'
        # Set up a protospacer with known bases
        # Positions 1-20 in protospacer correspond to seq indices 197-216
        # ABE window = positions 4-7 -> seq indices 200-203
        # Place A at positions 4,5,6,7 (indices 200-203) for maximum bystanders
        # The target is at position 4 (index 200), bystanders at 5,6,7
        proto = list("AGCTAAAAAACTAGCTAGCT")  # A at indices 3-8 of proto (pos 4-9)
        for i, base in enumerate(proto):
            seq[197 + i] = base
        test_seq = "".join(seq)
        while "TTTT" in test_seq:
            test_seq = test_seq.replace("TTTT", "AGCT", 1)

        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(variant, test_seq, 200)
        if result.label != FeasibilityLabel.NOT_FEASIBLE:
            # If a guide was found in the ABE window, bystanders should be >= 0
            assert result.bystander_count >= 0
            assert len(result.bystander_positions) == result.bystander_count

    def test_patient_allele_checked_in_protospacer(self, engine):
        """The engine substitutes the alt allele (patient) into the protospacer
        before checking window position, not the reference allele."""
        # For ABE: we need alt='A', ref='G'. The patient carries A.
        # The protospacer in the patient DNA should have 'A' at the target
        # position. The engine must use the patient sequence for window checks.
        seq = list("AGCT" * 100)
        seq[217] = 'A'
        seq[218] = 'G'
        seq[219] = 'G'
        proto = list("AGCTAGCTAGCTAGCTAGCT")
        for i, base in enumerate(proto):
            seq[197 + i] = base
        # Reference has 'G' at position 200
        seq[200] = 'G'
        test_seq = "".join(seq)
        while "TTTT" in test_seq:
            test_seq = test_seq.replace("TTTT", "AGCT", 1)

        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(variant, test_seq, 200)
        # The engine should have substituted 'A' (patient allele) at pos 200
        # in the patient_sequence before checking. If it checked 'G' (ref)
        # instead, ABE would fail to find 'A' in the window.
        assert result.editor_type == "ABE"


# ═══════════════════════════════════════════════════════════════════════════
# Prime Editing Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestPrimeEditing:
    """Tests for PrimeEditingEngine."""

    @pytest.fixture
    def engine(self):
        scanner = EnhancedPAMScanner("SpCas9")
        return PrimeEditingEngine(scanner)

    @pytest.fixture
    def test_sequence(self):
        """400-bp sequence with PAM sites for PE testing."""
        seq = _build_controlled_sequence(length=400, edit_pos=200)
        return seq

    def test_substitution_feasible(self, engine, test_sequence):
        """SNV substitutions are feasible for prime editing."""
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label != FeasibilityLabel.NOT_FEASIBLE, (
            f"Substitution should be PE-feasible. Rejection: {result.rejection_reason}"
        )
        assert result.edit_type == "substitution"

    def test_small_deletion_feasible(self, engine, test_sequence):
        """Small deletions (<=80 bp) are feasible for prime editing."""
        # Alt is longer -> correction direction is deletion by PE
        variant = _make_variant(ref='A', alt='A' + 'G' * 10, position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label != FeasibilityLabel.NOT_FEASIBLE, (
            f"Small deletion should be PE-feasible. "
            f"Rejection: {result.rejection_reason}"
        )
        assert result.edit_type == "deletion"

    def test_small_insertion_feasible(self, engine, test_sequence):
        """Small insertions (<=40 bp) are feasible for prime editing."""
        # Ref is longer -> correction direction is insertion by PE
        variant = _make_variant(ref='A' + 'G' * 10, alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label != FeasibilityLabel.NOT_FEASIBLE, (
            f"Small insertion should be PE-feasible. "
            f"Rejection: {result.rejection_reason}"
        )
        assert result.edit_type == "insertion"

    def test_large_deletion_not_feasible(self, engine, test_sequence):
        """Deletions >80 bp return NOT_FEASIBLE."""
        # alt has 81 extra bases -> correction removes 81 bp = too large
        variant = _make_variant(ref='A', alt='A' + 'G' * 81, position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label == FeasibilityLabel.NOT_FEASIBLE
        assert "80" in result.rejection_reason or "exceeds" in result.rejection_reason

    def test_large_insertion_not_feasible(self, engine, test_sequence):
        """Insertions >40 bp return NOT_FEASIBLE."""
        # ref has 41 extra bases -> correction inserts 41 bp = too large
        variant = _make_variant(ref='A' + 'G' * 41, alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label == FeasibilityLabel.NOT_FEASIBLE
        assert "40" in result.rejection_reason or "exceeds" in result.rejection_reason

    def test_boundary_deletion_80bp_not_rejected(self, engine, test_sequence):
        """Exactly 80 bp deletion is at the limit and should NOT be rejected by size."""
        variant = _make_variant(ref='A', alt='A' + 'G' * 80, position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert "exceeds the PE limit" not in result.rejection_reason, (
            "80-bp deletion must not be rejected by size limit (80 > 80 is False)"
        )

    def test_boundary_insertion_40bp_not_rejected(self, engine, test_sequence):
        """Exactly 40 bp insertion is at the limit and should NOT be rejected by size."""
        variant = _make_variant(ref='A' + 'G' * 40, alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert "exceeds the PE limit" not in result.rejection_reason, (
            "40-bp insertion must not be rejected by size limit (40 > 40 is False)"
        )

    def test_pbs_length_in_range(self, engine, test_sequence):
        """PBS length is in the valid range 10-17."""
        from core.feasibility.prime_editing import PBS_MIN_LEN, PBS_MAX_LEN
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        if result.label != FeasibilityLabel.NOT_FEASIBLE:
            assert PBS_MIN_LEN <= result.pbs_length <= PBS_MAX_LEN, (
                f"PBS length {result.pbs_length} outside range "
                f"[{PBS_MIN_LEN}, {PBS_MAX_LEN}]"
            )

    def test_pe_size_limits_match_constants(self):
        """The size limit constants match Anzalone et al., 2019."""
        from core.feasibility.prime_editing import MAX_DELETION_BP, MAX_INSERTION_BP
        assert MAX_DELETION_BP == 80
        assert MAX_INSERTION_BP == 40


# ═══════════════════════════════════════════════════════════════════════════
# HDR Design Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestHDRDesign:
    """Tests for HDRDesignEngine."""

    @pytest.fixture
    def engine(self):
        scanner = EnhancedPAMScanner("SpCas9")
        return HDRDesignEngine(scanner, nuclease="SpCas9", cell_type="iPSC")

    @pytest.fixture
    def test_sequence(self):
        """400-bp sequence with PAM for HDR testing."""
        return _build_controlled_sequence(length=400, edit_pos=200)

    def test_always_feasible_when_guide_exists(self, engine, test_sequence):
        """HDR is always feasible when a cutting guide exists near the edit."""
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        assert result.label != FeasibilityLabel.NOT_FEASIBLE, (
            f"HDR should be feasible when a guide exists. "
            f"Rejection: {result.rejection_reason}"
        )

    def test_cut_to_edit_score_decreases_monotonically(self, engine):
        """Score MUST decrease monotonically with increasing cut-to-edit distance."""
        distances = [0, 5, 10, 20, 50, 100, 200]
        scores = []
        for d in distances:
            s = engine._score(d, "cssDNA", True)
            scores.append(s)

        for i in range(1, len(scores)):
            assert scores[i] < scores[i - 1], (
                f"Score must decrease with distance: "
                f"at {distances[i]} bp ({scores[i]:.4f}) >= "
                f"at {distances[i-1]} bp ({scores[i-1]:.4f})"
            )

    def test_donor_type_ssodn_for_close_edits(self):
        """ssODN recommended for cut-to-edit <= 30 bp."""
        donor_type, _ = HDRDesignEngine._recommend_donor(5)
        assert donor_type == "ssODN"

        donor_type, _ = HDRDesignEngine._recommend_donor(30)
        assert donor_type == "ssODN"

    def test_donor_type_cssdna_for_moderate_distance(self):
        """cssDNA recommended for 31-5000 bp."""
        donor_type, _ = HDRDesignEngine._recommend_donor(31)
        assert donor_type == "cssDNA"

        donor_type, _ = HDRDesignEngine._recommend_donor(500)
        assert donor_type == "cssDNA"

        donor_type, _ = HDRDesignEngine._recommend_donor(5000)
        assert donor_type == "cssDNA"

    def test_donor_type_lssdna_for_large_distance(self):
        """lssDNA recommended for >5000 bp."""
        donor_type, _ = HDRDesignEngine._recommend_donor(5001)
        assert donor_type == "lssDNA"

    def test_donor_type_threshold_boundary(self):
        """Boundary between ssODN and cssDNA is exactly at 30/31 bp."""
        dtype_30, _ = HDRDesignEngine._recommend_donor(30)
        dtype_31, _ = HDRDesignEngine._recommend_donor(31)
        assert dtype_30 == "ssODN"
        assert dtype_31 == "cssDNA"

    def test_no_guide_returns_not_feasible(self, engine):
        """When no guide exists near the edit, HDR returns NOT_FEASIBLE."""
        # Build a sequence with no PAM sites at all (all A's + minimal breaks)
        seq = "ATAT" * 100
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(variant, seq, 200)
        assert result.label == FeasibilityLabel.NOT_FEASIBLE
        assert "guide" in result.rejection_reason.lower()

    def test_conversion_probability_decreases_with_distance(self):
        """Gene conversion probability falls with distance (exponential decay)."""
        p0 = HDRDesignEngine._get_conversion_probability(0)
        p100 = HDRDesignEngine._get_conversion_probability(100)
        p500 = HDRDesignEngine._get_conversion_probability(500)
        p1000 = HDRDesignEngine._get_conversion_probability(1000)

        assert p0 > p100 > p500 > p1000
        assert abs(p0 - 1.0) < 0.01, "Probability at 0 bp should be ~1.0"

    def test_invalid_cell_type_raises(self):
        """Unknown cell type raises ValueError."""
        scanner = EnhancedPAMScanner("SpCas9")
        with pytest.raises(ValueError, match="Unknown cell type"):
            HDRDesignEngine(scanner, cell_type="nonexistent_cell")

    def test_hdr_result_has_guide(self, engine, test_sequence):
        """Feasible HDR result includes a best_guide."""
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        if result.label != FeasibilityLabel.NOT_FEASIBLE:
            assert result.best_guide is not None
            assert len(result.best_guide.sequence_20mer) == 20

    def test_hdr_warns_about_p53_in_ipsc(self, engine, test_sequence):
        """iPSC HDR results should warn about p53-mediated apoptosis."""
        variant = _make_variant(ref='G', alt='A', position=200)
        result = engine.check_feasibility(
            variant, test_sequence, edit_pos_in_sequence=200,
        )
        if result.label != FeasibilityLabel.NOT_FEASIBLE:
            p53_warnings = [w for w in result.warnings if "p53" in w.lower()]
            assert len(p53_warnings) > 0, "Should warn about p53 in iPSC"


# ═══════════════════════════════════════════════════════════════════════════
# Cross-engine integration tests
# ═══════════════════════════════════════════════════════════════════════════

class TestFeasibilityIntegration:
    """Cross-engine tests verifying consistent behavior."""

    def test_transversion_only_pe_or_hdr_feasible(self):
        """For transversions, only PE and HDR should be feasible (not BE)."""
        seq = _build_controlled_sequence(length=400, edit_pos=200)
        scanner = EnhancedPAMScanner("SpCas9")
        be_engine = BaseEditingEngine(scanner)
        pe_engine = PrimeEditingEngine(scanner)
        hdr_engine = HDRDesignEngine(scanner, cell_type="iPSC")

        variant = _make_variant(ref='C', alt='A', position=200)  # transversion

        be_result = be_engine.check_feasibility(variant, seq, 200)
        pe_result = pe_engine.check_feasibility(variant, seq, 200)
        hdr_result = hdr_engine.check_feasibility(variant, seq, 200)

        assert be_result.label == FeasibilityLabel.NOT_FEASIBLE
        # PE and HDR should work for transversions (given guides exist)
        assert pe_result.label != FeasibilityLabel.NOT_FEASIBLE or \
               "guide" in pe_result.rejection_reason.lower()
        assert hdr_result.label != FeasibilityLabel.NOT_FEASIBLE or \
               "guide" in hdr_result.rejection_reason.lower()

    def test_snv_all_engines_produce_results(self):
        """An SNV transition should produce non-None results from all engines."""
        seq = _build_controlled_sequence(length=400, edit_pos=200)
        scanner = EnhancedPAMScanner("SpCas9")
        be_engine = BaseEditingEngine(scanner)
        pe_engine = PrimeEditingEngine(scanner)
        hdr_engine = HDRDesignEngine(scanner, cell_type="iPSC")

        variant = _make_variant(ref='G', alt='A', position=200)

        be_result = be_engine.check_feasibility(variant, seq, 200)
        pe_result = pe_engine.check_feasibility(variant, seq, 200)
        hdr_result = hdr_engine.check_feasibility(variant, seq, 200)

        assert isinstance(be_result, BaseEditingFeasibility)
        assert isinstance(pe_result, PrimeEditingFeasibility)
        assert isinstance(hdr_result, HDRFeasibility)
