"""
Tests for CRISPRArchitect v3 Multi-Nuclease Engine
====================================================

Tests the multi-editor evaluation system that systematically checks
ABE8e/BE4max with SpCas9, enFnCas9, SpCas9-NG, and SpRY.

References
----------
Richter et al., Nat Biotechnol, 2020 (ABE8e window: 3-9)
Nishimasu et al., Science, 2018 (SpCas9-NG: NG PAM)
Walton et al., Science, 2020 (SpRY: near-PAMless)
Acharya et al., Nat Commun 15:5471, 2024 (enFnCas9: NRG PAM)
"""

from __future__ import annotations

import pytest

from core.models import (
    BaseEditingFeasibility,
    CodingAnnotation,
    ConsequenceType,
    ExonRecord,
    FeasibilityLabel,
    GenomicVariantInput,
    NormalizedVariant,
    ReferenceValidation,
    TranscriptCoordinate,
    TranscriptInfo,
)
from core.feasibility.pam_scan import EnhancedPAMScanner
from core.feasibility.base_editing import (
    BaseEditingEngine,
    _get_editor_window,
    _get_nucleases_for_editor,
    ABE_WINDOW_START,
    ABE_WINDOW_END,
    CBE_WINDOW_START,
    CBE_WINDOW_END,
)


# ── Helpers ──────────────────────────────────────────────────────────

def _make_variant(ref, alt):
    """Build a minimal NormalizedVariant for testing."""
    return NormalizedVariant(
        input=GenomicVariantInput(
            chromosome="17", position=100, ref_allele=ref, alt_allele=alt,
            gene_symbol="TEST",
        ),
        transcript=TranscriptInfo(
            transcript_id="ENST0001", gene_symbol="TEST",
            gene_id="ENSG0001", chromosome="17",
            start=1, end=1000, strand=1, biotype="protein_coding",
            is_canonical=True,
            exons=[ExonRecord("ENSE001", 1, 1, 1000, 1, "17")],
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


def _make_test_sequence_with_ngg():
    """Build a sequence with a known NGG PAM placing edit at ABE window pos 4."""
    seq_list = list("ATCGATCG" * 50)  # 400 bases
    seq_list[200] = "A"  # edit position
    seq_list[217] = "A"
    seq_list[218] = "G"
    seq_list[219] = "G"  # NGG PAM at 217
    for i in range(197, 217):
        if i != 200:
            seq_list[i] = ["A", "G", "C", "T"][(i - 197) % 4]
    test_seq = "".join(seq_list)
    while "TTTT" in test_seq:
        test_seq = test_seq.replace("TTTT", "ATCG", 1)
    return test_seq


# ── Tests: Editor Profile Configuration ──────────────────────────────

class TestEditorProfiles:
    """Tests for BASE_EDITOR_PROFILES configuration."""

    def test_abe8e_window_is_broader_than_abe7(self):
        """ABE8e has window 3-9 vs ABE7.10 window 4-7 (Richter et al., 2020)."""
        abe7_win = _get_editor_window("ABE7.10")
        abe8e_win = _get_editor_window("ABE8e")
        assert abe8e_win[0] <= abe7_win[0], "ABE8e start should be <= ABE7.10 start"
        assert abe8e_win[1] >= abe7_win[1], "ABE8e end should be >= ABE7.10 end"

    def test_abe8e_window_exact_values(self):
        """ABE8e window is 3-9 per Richter et al., Nat Biotechnol, 2020."""
        start, end = _get_editor_window("ABE8e")
        assert start == 3
        assert end == 9

    def test_be4max_window_matches_cbe(self):
        """BE4max window should match standard CBE window."""
        start, end = _get_editor_window("BE4max")
        assert start == CBE_WINDOW_START
        assert end == CBE_WINDOW_END

    def test_legacy_abe_fallback(self):
        """Legacy 'ABE' name falls back to ABE7.10 window."""
        start, end = _get_editor_window("ABE")
        assert start == ABE_WINDOW_START
        assert end == ABE_WINDOW_END

    def test_abe8e_supports_multiple_nucleases(self):
        """ABE8e should be compatible with SpCas9, enFnCas9, SpCas9-NG, SpRY."""
        nucleases = _get_nucleases_for_editor("ABE8e")
        assert "SpCas9" in nucleases
        assert "enFnCas9" in nucleases
        assert "SpCas9-NG" in nucleases
        assert "SpRY" in nucleases

    def test_abe7_only_spcas9(self):
        """ABE7.10 is only validated with SpCas9."""
        nucleases = _get_nucleases_for_editor("ABE7.10")
        assert nucleases == ["SpCas9"]

    def test_unknown_editor_fallback(self):
        """Unknown editor name falls back to SpCas9."""
        nucleases = _get_nucleases_for_editor("NonexistentEditor")
        assert nucleases == ["SpCas9"]


# ── Tests: PAM Scanner with New Nucleases ────────────────────────────

class TestNewNucleases:
    """Tests for SpCas9-NG and SpRY nuclease support."""

    def test_spcas9_ng_scanner_initializes(self):
        """SpCas9-NG should be a valid nuclease."""
        scanner = EnhancedPAMScanner("SpCas9-NG")
        assert scanner.pam == "NG"

    def test_spy_scanner_initializes(self):
        """SpRY should be a valid nuclease."""
        scanner = EnhancedPAMScanner("SpRY")
        assert scanner.pam == "NNN"

    def test_spcas9_ng_finds_more_pams_than_spcas9(self):
        """SpCas9-NG (NG PAM) should find >= SpCas9 (NGG PAM) guides."""
        seq = _make_test_sequence_with_ngg()
        sp_guides = EnhancedPAMScanner("SpCas9").scan(seq, 200, window_bp=50)
        ng_guides = EnhancedPAMScanner("SpCas9-NG").scan(seq, 200, window_bp=50)
        assert len(ng_guides) >= len(sp_guides)

    def test_spy_finds_most_pams(self):
        """SpRY (NNN PAM) should find the most guides of any nuclease."""
        seq = _make_test_sequence_with_ngg()
        sp_guides = EnhancedPAMScanner("SpCas9").scan(seq, 200, window_bp=50)
        ry_guides = EnhancedPAMScanner("SpRY").scan(seq, 200, window_bp=50)
        assert len(ry_guides) >= len(sp_guides)

    def test_enfncas9_finds_more_pams_than_spcas9(self):
        """enFnCas9 (NRG) should find >= SpCas9 (NGG) sites."""
        seq = _make_test_sequence_with_ngg()
        sp_guides = EnhancedPAMScanner("SpCas9").scan(seq, 200, window_bp=50)
        en_guides = EnhancedPAMScanner("enFnCas9").scan(seq, 200, window_bp=50)
        assert len(en_guides) >= len(sp_guides)


# ── Tests: Multi-Editor Evaluation ───────────────────────────────────

class TestMultiEditorEvaluation:
    """Tests for BaseEditingEngine.check_all_editors()."""

    def setup_method(self):
        self.engine = BaseEditingEngine(EnhancedPAMScanner("SpCas9"))
        self.seq = _make_test_sequence_with_ngg()
        self.edit_pos = 200

    def test_abe_variant_returns_multiple_results(self):
        """ABE-amenable variant should produce results for multiple editors."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        assert len(results) > 1, "Should have results from multiple editor-nuclease combos"

    def test_abe_variant_has_feasible_results(self):
        """At least one editor-nuclease combo should be feasible for ABE variant."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        feasible = [r for r in results if r.label != FeasibilityLabel.NOT_FEASIBLE]
        assert len(feasible) >= 1

    def test_cbe_variant_returns_cbe_editors(self):
        """CBE-amenable variant should return CBE editors, not ABE."""
        v = _make_variant("T", "C")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        for r in results:
            if r.editor_type:
                assert r.editor_type == "CBE"

    def test_transversion_rejected_by_all_editors(self):
        """Transversion should be rejected universally."""
        v = _make_variant("C", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        assert len(results) == 1
        assert results[0].label == FeasibilityLabel.NOT_FEASIBLE

    def test_indel_rejected(self):
        """Indels should be rejected by all base editors."""
        v = _make_variant("A", "ATG")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        assert len(results) == 1
        assert results[0].label == FeasibilityLabel.NOT_FEASIBLE

    def test_results_sorted_by_score(self):
        """Results should be sorted: feasible first, then by score descending."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        feasible = [r for r in results if r.label != FeasibilityLabel.NOT_FEASIBLE]
        if len(feasible) > 1:
            for i in range(len(feasible) - 1):
                assert feasible[i].score >= feasible[i + 1].score

    def test_metadata_contains_editor_and_nuclease(self):
        """Each result should have editor and nuclease in metadata."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        for r in results:
            assert "editor" in r.metadata
            assert "nuclease" in r.metadata

    def test_evidence_tier_tagged(self):
        """Each result should have evidence_tier in metadata."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        for r in results:
            assert "evidence_tier" in r.metadata
            assert r.metadata["evidence_tier"] in ("A", "B")

    def test_abe8e_broader_window_finds_more_positions(self):
        """ABE8e window (3-9) should potentially find guides that ABE7.10 (4-7) misses."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(v, self.seq, self.edit_pos)
        abe7_results = [
            r for r in results
            if r.metadata.get("editor") == "ABE7.10"
            and r.label != FeasibilityLabel.NOT_FEASIBLE
        ]
        abe8e_spcas9_results = [
            r for r in results
            if r.metadata.get("editor") == "ABE8e"
            and r.metadata.get("nuclease") == "SpCas9"
            and r.label != FeasibilityLabel.NOT_FEASIBLE
        ]
        # ABE8e should find at least as many positions as ABE7.10 with same nuclease
        assert len(abe8e_spcas9_results) >= len(abe7_results)

    def test_custom_editor_list(self):
        """check_all_editors should accept a custom editor list."""
        v = _make_variant("G", "A")
        results = self.engine.check_all_editors(
            v, self.seq, self.edit_pos,
            editor_names=["ABE8e"],
        )
        editors_used = set(r.metadata.get("editor") for r in results)
        # Should only have ABE8e variants
        for ed in editors_used:
            assert "ABE8e" in ed or ed == "ABE8e"


# ── Tests: Backward Compatibility ────────────────────────────────────

class TestBackwardCompatibility:
    """Ensure legacy check_feasibility still works alongside check_all_editors."""

    def setup_method(self):
        self.engine = BaseEditingEngine(EnhancedPAMScanner("SpCas9"))
        self.seq = _make_test_sequence_with_ngg()
        self.edit_pos = 200

    def test_legacy_check_feasibility_still_works(self):
        """The original check_feasibility method should still function."""
        v = _make_variant("G", "A")
        result = self.engine.check_feasibility(v, self.seq, self.edit_pos)
        assert isinstance(result, BaseEditingFeasibility)
        assert result.editor_type == "ABE"

    def test_legacy_and_multi_agree_on_feasibility(self):
        """Legacy and multi-editor should agree on basic feasibility."""
        v = _make_variant("G", "A")
        legacy = self.engine.check_feasibility(v, self.seq, self.edit_pos)
        multi = self.engine.check_all_editors(v, self.seq, self.edit_pos)

        # If legacy says FEASIBLE, at least one multi result should too
        if legacy.label != FeasibilityLabel.NOT_FEASIBLE:
            feasible_multi = [r for r in multi if r.label != FeasibilityLabel.NOT_FEASIBLE]
            assert len(feasible_multi) >= 1
