"""
Tests for CRISPRArchitect v2 core models.

Verifies that all dataclasses and enums instantiate correctly,
that default values are biologically sensible, and that the
Strategy/Scoring model hierarchy works.
"""

from __future__ import annotations

import pytest

from core.models import (
    BenchmarkCase,
    BenchmarkResult,
    BenchmarkSummary,
    BaseEditingFeasibility,
    CodingAnnotation,
    ConsequenceType,
    EditModality,
    EvidenceTier,
    ExonRecord,
    FeasibilityBundle,
    FeasibilityLabel,
    GenomicVariantInput,
    GuideCandidate,
    HDRFeasibility,
    NormalizedVariant,
    PipelineResult,
    PrimeEditingFeasibility,
    ReferenceValidation,
    RiskLevel,
    ScoredStrategy,
    Strategy,
    StrategyStep,
    TranscriptCoordinate,
    TranscriptInfo,
)


class TestEnums:
    def test_consequence_types_complete(self):
        """All standard consequence types are defined."""
        expected = {
            "synonymous", "missense", "nonsense",
            "splice_donor", "splice_acceptor", "splice_region",
            "frameshift", "inframe_insertion", "inframe_deletion",
            "intronic", "5_prime_UTR", "3_prime_UTR",
            "non_coding", "unknown",
        }
        actual = {c.value for c in ConsequenceType}
        assert expected == actual

    def test_edit_modalities_complete(self):
        expected = {"ABE", "CBE", "PE", "HDR_ssODN", "HDR_cssDNA",
                    "HDR_lssDNA", "HDR_dsDNA", "exon_deletion"}
        actual = {m.value for m in EditModality}
        assert expected == actual

    def test_feasibility_labels(self):
        assert FeasibilityLabel.FEASIBLE.value == "feasible"
        assert FeasibilityLabel.NOT_FEASIBLE.value == "not_feasible"

    def test_evidence_tiers(self):
        assert EvidenceTier.A.value == "A"
        assert EvidenceTier.C.value == "C"


class TestTranscriptModels:
    def test_exon_record_length(self):
        exon = ExonRecord("E1", 1, 100, 200, 1, "1")
        assert exon.length == 101  # inclusive

    def test_transcript_info_properties(self):
        t = TranscriptInfo(
            transcript_id="ENST00000358273",
            gene_symbol="NF1", gene_id="ENSG00000196712",
            chromosome="17", start=31094927, end=31377677,
            strand=-1, biotype="protein_coding", is_canonical=True,
            exons=[
                ExonRecord("E1", 1, 31094927, 31095068, -1, "17"),
                ExonRecord("E2", 2, 31100000, 31100200, -1, "17"),
            ],
        )
        assert t.n_exons == 2
        assert t.span_bp == 282751

    def test_transcript_coordinate_defaults(self):
        tc = TranscriptCoordinate(
            genomic_position=100, exon_number=1,
            transcript_position=50, cds_position=50,
            codon_index=17, codon_position=2,
            reference_codon="ATG", reference_aa="M",
            distance_to_exon_start=49, distance_to_exon_end=100,
        )
        assert tc.in_cds is True


class TestVariantModels:
    def test_genomic_variant_input_defaults(self):
        v = GenomicVariantInput(
            chromosome="17", position=31232193,
            ref_allele="C", alt_allele="T",
        )
        assert v.gene_symbol is None
        assert v.species == "homo_sapiens"

    def test_coding_annotation_defaults(self):
        ca = CodingAnnotation(consequence=ConsequenceType.MISSENSE)
        assert ca.hgvs_c == ""
        assert ca.splice_proximity is None


class TestFeasibilityModels:
    def test_guide_candidate_defaults(self):
        gc = GuideCandidate(sequence_20mer="ATCGATCGATCGATCGATCG")
        assert len(gc.sequence_20mer) == 20
        assert gc.gc_content == 0.0  # not computed by default

    def test_base_editing_default_not_feasible(self):
        be = BaseEditingFeasibility()
        assert be.label == FeasibilityLabel.NOT_FEASIBLE
        assert be.editor_type is None

    def test_feasibility_bundle_best_modality_empty(self):
        v = _make_minimal_variant()
        bundle = FeasibilityBundle(variant=v, mutation_index=0)
        assert bundle.best_modality() is None

    def test_feasibility_bundle_best_modality_with_be(self):
        v = _make_minimal_variant()
        bundle = FeasibilityBundle(
            variant=v, mutation_index=0,
            base_editing_results=[
                BaseEditingFeasibility(
                    label=FeasibilityLabel.FEASIBLE,
                    editor_type="ABE",
                    score=0.85,
                ),
            ],
        )
        assert bundle.best_modality() == EditModality.ABE


class TestStrategyModels:
    def test_strategy_not_rejected_by_default(self):
        s = Strategy(name="Test")
        assert not s.is_rejected

    def test_strategy_rejected_with_reasons(self):
        s = Strategy(
            name="Bad Strategy",
            rejection_reasons=["No PAM available"],
        )
        assert s.is_rejected

    def test_scored_strategy_properties(self):
        s = Strategy(name="Test BE", evidence_tier=EvidenceTier.A)
        scored = ScoredStrategy(
            strategy=s, overall_score=0.85, rank=1,
        )
        assert scored.strategy_name == "Test BE"
        assert scored.confidence == "A"


class TestBenchmarkModels:
    def test_benchmark_summary_to_dict(self):
        summary = BenchmarkSummary(
            n_cases=10, top1_accuracy=0.8,
            top3_accuracy=0.9, rejection_accuracy=1.0,
        )
        d = summary.to_dict()
        assert d["n_cases"] == 10
        assert d["top1_accuracy"] == 0.8

    def test_pipeline_result_top_strategy(self):
        s = Strategy(name="Top")
        scored = ScoredStrategy(strategy=s, overall_score=0.9, rank=1)
        result = PipelineResult(
            transcript=TranscriptInfo(
                "T1", "G", "G1", "1", 1, 100, 1, "pc", True, [],
            ),
            strategies=[scored],
        )
        assert result.top_strategy.strategy_name == "Top"

    def test_pipeline_result_no_strategies(self):
        result = PipelineResult(
            transcript=TranscriptInfo(
                "T1", "G", "G1", "1", 1, 100, 1, "pc", True, [],
            ),
        )
        assert result.top_strategy is None


class TestScoringIntegration:
    """Test the v2 scoring engine produces biologically sensible rankings."""

    def test_be_ranks_above_hdr_in_ipsc(self):
        """DSB-free BE should always rank above HDR in iPSC context."""
        from core.pipeline.strategy_stage import StrategyScorer

        scorer = StrategyScorer()
        be = Strategy(
            name="Base Editing", num_dsbs=0, num_rounds=1,
            evidence_tier=EvidenceTier.A,
            modality_prior_score=0.90,
            rearrangement_risk=RiskLevel.LOW,
        )
        hdr = Strategy(
            name="HDR", num_dsbs=1, num_rounds=1, num_donors=1,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.72,
            rearrangement_risk=RiskLevel.LOW,
            p53_active=True,
        )
        ranked = scorer.rank([hdr, be], [])
        assert ranked[0].strategy_name == "Base Editing"

    def test_sequential_safer_than_simultaneous(self):
        """Sequential editing should be safer than simultaneous dual DSB."""
        from core.pipeline.strategy_stage import StrategyScorer

        scorer = StrategyScorer()
        seq = Strategy(
            name="Sequential HDR", num_dsbs=1, num_rounds=2,
            simultaneous_dsbs=False,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.65,
            p53_active=True,
        )
        sim = Strategy(
            name="Simultaneous HDR", num_dsbs=2, num_rounds=1,
            simultaneous_dsbs=True,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.65,
            p53_active=True,
        )
        scored_seq = scorer.score_strategy(seq, [])
        scored_sim = scorer.score_strategy(sim, [])
        assert scored_seq.safety_score > scored_sim.safety_score


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _make_minimal_variant() -> NormalizedVariant:
    return NormalizedVariant(
        input=GenomicVariantInput("1", 100, "C", "T"),
        transcript=TranscriptInfo("T1", "G", "G1", "1", 1, 100, 1, "pc", True, []),
        transcript_coord=TranscriptCoordinate(
            100, 1, 100, 100, 34, 1, "ATG", "M", 99, 50,
        ),
        coding=CodingAnnotation(consequence=ConsequenceType.MISSENSE),
        ref_validation=ReferenceValidation(
            True, "C", "C", "C", "C",
        ),
    )
