"""
Comprehensive tests for CRISPRArchitect v2 strategy generator and scorer.

Tests the StrategyGenerator (bundle -> strategy enumeration) and
StrategyScorer (multi-objective scoring + ranking) using synthetic
feasibility bundles.  No network calls, no Ensembl API.

Python 3.9 compatible.  Uses pytest.
"""

from __future__ import annotations

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
from core.mosaic.generator import StrategyGenerator
from core.pipeline.strategy_stage import StrategyScorer


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


def _make_guide():
    """Create a minimal GuideCandidate."""
    return GuideCandidate(
        sequence_20mer="AGCTAGCTAGCTAGCTAGCT",
        pam_sequence="AGG",
        strand="+",
        cut_position=200,
        distance_to_edit=3,
        gc_content=0.50,
        score=0.80,
    )


def _make_be_feasible(editor_type="ABE", score=0.85, bystanders=0):
    """Create a FEASIBLE BaseEditingFeasibility result."""
    return BaseEditingFeasibility(
        label=FeasibilityLabel.FEASIBLE,
        editor_type=editor_type,
        best_guide=_make_guide(),
        target_position_in_window=5,
        bystander_count=bystanders,
        bystander_positions=list(range(6, 6 + bystanders)),
        bystander_consequences=[ConsequenceType.UNKNOWN] * bystanders,
        compatible_nucleases=["SpCas9"],
        score=score,
        metadata={"editor": editor_type, "editable_bystanders": bystanders},
    )


def _make_be_not_feasible():
    """Create a NOT_FEASIBLE BaseEditingFeasibility result."""
    return BaseEditingFeasibility(
        label=FeasibilityLabel.NOT_FEASIBLE,
        rejection_reason="Transversion cannot be corrected by ABE or CBE.",
    )


def _make_pe_feasible(score=0.75):
    """Create a FEASIBLE PrimeEditingFeasibility result."""
    return PrimeEditingFeasibility(
        label=FeasibilityLabel.FEASIBLE,
        best_guide=_make_guide(),
        pbs_length=13,
        rt_template_length=15,
        rt_template_sequence="AGCTAGCTAGCTAGC",
        edit_type="substitution",
        score=score,
        metadata={"editor": "Prime Editor"},
    )


def _make_pe_not_feasible():
    """Create a NOT_FEASIBLE PrimeEditingFeasibility result."""
    return PrimeEditingFeasibility(
        label=FeasibilityLabel.NOT_FEASIBLE,
        rejection_reason="Deletion exceeds PE limit.",
    )


def _make_hdr_feasible(score=0.70, cut_to_edit=5, donor_type="cssDNA"):
    """Create a FEASIBLE HDRFeasibility result."""
    return HDRFeasibility(
        label=FeasibilityLabel.FEASIBLE,
        best_guide=_make_guide(),
        cut_to_edit_distance=cut_to_edit,
        recommended_donor_type=donor_type,
        donor_length_estimate=601,
        homology_arm_length=300,
        conversion_probability=0.90,
        pam_disruption_possible=True,
        score=score,
        metadata={"nuclease": "SpCas9", "cell_type": "iPSC",
                  "hdr_base_efficiency": 0.08, "nuclease_hdr_multiplier": 1.0},
    )


def _make_hdr_not_feasible():
    """Create a NOT_FEASIBLE HDRFeasibility result."""
    return HDRFeasibility(
        label=FeasibilityLabel.NOT_FEASIBLE,
        rejection_reason="No cutting guide found.",
    )


def _make_bundle(
    variant=None,
    mutation_index=0,
    be_results=None,
    pe_result=None,
    hdr_result=None,
):
    """Assemble a FeasibilityBundle with specified results."""
    if variant is None:
        variant = _make_variant()
    bundle = FeasibilityBundle(
        variant=variant,
        mutation_index=mutation_index,
    )
    if be_results is not None:
        bundle.base_editing_results = be_results
    if pe_result is not None:
        bundle.prime_editing_result = pe_result
    if hdr_result is not None:
        bundle.hdr_result = hdr_result
    return bundle


# ═══════════════════════════════════════════════════════════════════════════
# Generator Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestStrategyGenerator:
    """Tests for StrategyGenerator.generate()."""

    @pytest.fixture
    def generator(self):
        return StrategyGenerator(p53_active=True)

    def test_empty_bundles_produce_no_strategies(self, generator):
        """Empty input produces zero strategies."""
        strategies = generator.generate([])
        assert len(strategies) == 0

    def test_single_be_feasible_produces_be_strategy(self, generator):
        """A single BE-feasible bundle produces a Single-step Base Editing strategy."""
        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE")],
        )
        strategies = generator.generate([bundle])
        be_strategies = [s for s in strategies if "Base Editing" in s.name]
        assert len(be_strategies) >= 1, (
            f"Expected at least 1 BE strategy, got {len(be_strategies)}. "
            f"All strategies: {[s.name for s in strategies]}"
        )
        # Verify the strategy properties
        be_strat = be_strategies[0]
        assert be_strat.num_dsbs == 0
        assert be_strat.num_rounds == 1
        assert not be_strat.is_rejected

    def test_single_pe_feasible_produces_pe_strategy(self, generator):
        """A single PE-feasible bundle produces a Single-step Prime Editing strategy."""
        bundle = _make_bundle(
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert len(pe_strategies) >= 1, (
            f"Expected PE strategy. All: {[s.name for s in strategies]}"
        )
        pe_strat = pe_strategies[0]
        assert pe_strat.num_dsbs == 0
        assert pe_strat.num_rounds == 1

    def test_single_hdr_feasible_produces_hdr_strategy(self, generator):
        """A single HDR-feasible bundle produces a Single-step HDR strategy."""
        bundle = _make_bundle(
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([bundle])
        hdr_strategies = [s for s in strategies if "HDR" in s.name]
        assert len(hdr_strategies) >= 1, (
            f"Expected HDR strategy. All: {[s.name for s in strategies]}"
        )
        hdr_strat = hdr_strategies[0]
        assert hdr_strat.num_dsbs == 1
        assert hdr_strat.num_donors == 1

    def test_two_be_feasible_produces_dual_be(self, generator):
        """Two BE-feasible bundles produce a Dual Base Editing strategy."""
        b1 = _make_bundle(
            mutation_index=0,
            be_results=[_make_be_feasible("ABE", score=0.85)],
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='A', alt='G', position=200),
            mutation_index=1,
            be_results=[_make_be_feasible("CBE", score=0.80)],
        )
        strategies = generator.generate([b1, b2])
        dual_be = [s for s in strategies if "Dual Base Editing" in s.name]
        assert len(dual_be) >= 1, (
            f"Expected Dual BE strategy. All: {[s.name for s in strategies]}"
        )
        dual = dual_be[0]
        assert dual.num_dsbs == 0
        assert dual.num_distinct_guides == 2
        assert len(dual.steps) == 2

    def test_infeasible_bundles_produce_no_strategies(self, generator):
        """Bundles where all modalities are NOT_FEASIBLE produce no strategies."""
        bundle = _make_bundle(
            be_results=[_make_be_not_feasible()],
            pe_result=_make_pe_not_feasible(),
            hdr_result=_make_hdr_not_feasible(),
        )
        strategies = generator.generate([bundle])
        # All results are NOT_FEASIBLE, so no strategy should be generated
        assert len(strategies) == 0, (
            f"Expected 0 strategies for fully infeasible bundle, "
            f"got {len(strategies)}: {[s.name for s in strategies]}"
        )

    def test_all_modalities_feasible_produces_multiple(self, generator):
        """A bundle with all modalities feasible produces BE + PE + HDR strategies."""
        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE")],
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([bundle])
        names = [s.name for s in strategies]
        assert any("Base Editing" in n for n in names)
        assert any("Prime Editing" in n for n in names)
        assert any("HDR" in n for n in names)

    def test_two_bundles_generate_hybrid_strategies(self, generator):
        """Two bundles with mixed modalities produce hybrid BE+HDR or PE+HDR."""
        b1 = _make_bundle(
            mutation_index=0,
            be_results=[_make_be_feasible("ABE")],
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='C', alt='A', position=300),
            mutation_index=1,
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([b1, b2])
        hybrid = [s for s in strategies if "Hybrid" in s.name]
        assert len(hybrid) >= 1, (
            f"Expected hybrid strategy. All: {[s.name for s in strategies]}"
        )

    def test_sequential_hdr_for_two_hdr_bundles(self, generator):
        """Two HDR-feasible bundles produce a Sequential HDR strategy."""
        b1 = _make_bundle(
            mutation_index=0,
            hdr_result=_make_hdr_feasible(score=0.70),
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='A', alt='T', position=300),
            mutation_index=1,
            hdr_result=_make_hdr_feasible(score=0.65),
        )
        strategies = generator.generate([b1, b2])
        seq_hdr = [s for s in strategies if "Sequential HDR" in s.name]
        assert len(seq_hdr) >= 1
        assert seq_hdr[0].num_rounds == 2
        assert seq_hdr[0].num_donors == 2


# ═══════════════════════════════════════════════════════════════════════════
# Scorer Tests
# ═══════════════════════════════════════════════════════════════════════════

class TestStrategyScorer:
    """Tests for StrategyScorer."""

    @pytest.fixture
    def scorer(self):
        return StrategyScorer()

    def test_be_scores_higher_than_hdr_in_ipsc(self, scorer):
        """BE (0 DSBs) scores higher than HDR (1 DSB) in iPSC context."""
        be = Strategy(
            name="Single-step Base Editing",
            num_dsbs=0,
            num_rounds=1,
            evidence_tier=EvidenceTier.A,
            modality_prior_score=0.90,
            rearrangement_risk=RiskLevel.LOW,
        )
        hdr = Strategy(
            name="Single-step HDR",
            num_dsbs=1,
            num_rounds=1,
            num_donors=1,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.72,
            rearrangement_risk=RiskLevel.LOW,
            p53_active=True,
        )
        scored_be = scorer.score_strategy(be, [])
        scored_hdr = scorer.score_strategy(hdr, [])
        assert scored_be.overall_score > scored_hdr.overall_score, (
            f"BE ({scored_be.overall_score:.3f}) must score higher than "
            f"HDR ({scored_hdr.overall_score:.3f}) in iPSC"
        )

    def test_sequential_lower_complexity_than_single_round(self, scorer):
        """Sequential (2 rounds) has higher complexity score (worse) than single-round."""
        single = Strategy(
            name="Single HDR",
            num_dsbs=1,
            num_rounds=1,
            num_donors=1,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.72,
        )
        sequential = Strategy(
            name="Sequential HDR",
            num_dsbs=1,
            num_rounds=2,
            num_donors=2,
            num_distinct_guides=2,
            evidence_tier=EvidenceTier.B,
            modality_prior_score=0.68,
        )
        scored_single = scorer.score_strategy(single, [])
        scored_seq = scorer.score_strategy(sequential, [])
        # Higher complexity_score = worse (more complex)
        assert scored_seq.complexity_score > scored_single.complexity_score, (
            f"Sequential complexity ({scored_seq.complexity_score:.3f}) "
            f"must be higher than single ({scored_single.complexity_score:.3f})"
        )

    def test_safety_score_1_for_dsb_free(self, scorer):
        """Safety score = 1.0 for DSB-free strategies (BE, PE)."""
        be = Strategy(
            name="BE", num_dsbs=0, num_rounds=1,
            modality_prior_score=0.90,
        )
        scored = scorer.score_strategy(be, [])
        assert scored.safety_score == 1.0, (
            f"DSB-free safety should be 1.0, got {scored.safety_score}"
        )

    def test_safety_score_less_than_1_for_dsb(self, scorer):
        """Safety score < 1.0 for DSB-requiring strategies."""
        hdr = Strategy(
            name="HDR", num_dsbs=1, num_rounds=1,
            num_donors=1, modality_prior_score=0.72,
            p53_active=True,
        )
        scored = scorer.score_strategy(hdr, [])
        assert scored.safety_score < 1.0, (
            f"DSB strategy safety should be < 1.0, got {scored.safety_score}"
        )

    def test_dual_simultaneous_dsb_lower_safety_than_single(self, scorer):
        """Dual simultaneous DSB has lower safety than single DSB."""
        single_dsb = Strategy(
            name="Single HDR", num_dsbs=1, simultaneous_dsbs=False,
            num_rounds=1, modality_prior_score=0.72,
            p53_active=True,
        )
        dual_dsb = Strategy(
            name="Dual HDR", num_dsbs=2, simultaneous_dsbs=True,
            num_rounds=1, modality_prior_score=0.65,
            p53_active=True,
        )
        scored_single = scorer.score_strategy(single_dsb, [])
        scored_dual = scorer.score_strategy(dual_dsb, [])
        assert scored_dual.safety_score < scored_single.safety_score, (
            f"Dual sim DSB safety ({scored_dual.safety_score:.3f}) "
            f"must be < single DSB ({scored_single.safety_score:.3f})"
        )

    def test_pe_scores_between_be_and_hdr(self, scorer):
        """PE (nick, no DSB) should score between BE and HDR."""
        be = Strategy(
            name="BE", num_dsbs=0, num_rounds=1,
            evidence_tier=EvidenceTier.A,
            modality_prior_score=0.90,
        )
        pe = Strategy(
            name="PE", num_dsbs=0, num_rounds=1,
            evidence_tier=EvidenceTier.A,
            modality_prior_score=0.82,
        )
        hdr = Strategy(
            name="HDR", num_dsbs=1, num_rounds=1,
            num_donors=1, evidence_tier=EvidenceTier.B,
            modality_prior_score=0.72,
            p53_active=True,
        )
        s_be = scorer.score_strategy(be, [])
        s_pe = scorer.score_strategy(pe, [])
        s_hdr = scorer.score_strategy(hdr, [])
        assert s_be.overall_score > s_pe.overall_score, "BE > PE"
        assert s_pe.overall_score > s_hdr.overall_score, "PE > HDR"

    def test_rejected_strategies_excluded_from_ranking(self, scorer):
        """Rejected strategies (with rejection_reasons) are excluded from ranking."""
        good = Strategy(
            name="Good BE", num_dsbs=0, modality_prior_score=0.90,
        )
        rejected = Strategy(
            name="Rejected",
            rejection_reasons=["No PAM available"],
            modality_prior_score=0.50,
        )
        ranked = scorer.rank([good, rejected], [])
        names = [r.strategy_name for r in ranked]
        assert "Rejected" not in names
        assert "Good BE" in names

    def test_evidence_tier_a_scores_higher_confidence(self, scorer):
        """Tier A evidence yields higher confidence score than Tier C."""
        tier_a = Strategy(
            name="Tier A", evidence_tier=EvidenceTier.A,
            modality_prior_score=0.80,
        )
        tier_c = Strategy(
            name="Tier C", evidence_tier=EvidenceTier.C,
            modality_prior_score=0.80,
        )
        s_a = scorer.score_strategy(tier_a, [])
        s_c = scorer.score_strategy(tier_c, [])
        assert s_a.confidence_score > s_c.confidence_score

    def test_high_rearrangement_risk_penalized(self, scorer):
        """High rearrangement risk increases the risk score (worse)."""
        low_risk = Strategy(
            name="Low Risk", rearrangement_risk=RiskLevel.LOW,
            modality_prior_score=0.80,
        )
        high_risk = Strategy(
            name="High Risk", rearrangement_risk=RiskLevel.HIGH,
            modality_prior_score=0.80,
        )
        s_low = scorer.score_strategy(low_risk, [])
        s_high = scorer.score_strategy(high_risk, [])
        assert s_high.risk_score > s_low.risk_score


# ═══════════════════════════════════════════════════════════════════════════
# Integration Tests: Generator + Scorer
# ═══════════════════════════════════════════════════════════════════════════

class TestGeneratorScorerIntegration:
    """Tests that generator + scorer together produce correctly ranked output."""

    def test_ranked_list_is_correctly_ordered(self):
        """Generator + Scorer produces a list where rank 1 has the highest score."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        # Bundle with all modalities feasible
        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE", score=0.85)],
            pe_result=_make_pe_feasible(score=0.75),
            hdr_result=_make_hdr_feasible(score=0.65),
        )
        strategies = gen.generate([bundle])
        assert len(strategies) > 0

        ranked = scorer.rank(strategies, [bundle])
        assert len(ranked) > 0
        assert ranked[0].rank == 1

        # Verify monotonically non-increasing overall_score
        for i in range(1, len(ranked)):
            assert ranked[i - 1].overall_score >= ranked[i].overall_score, (
                f"Rank {i} score ({ranked[i-1].overall_score:.3f}) < "
                f"rank {i+1} ({ranked[i].overall_score:.3f})"
            )

    def test_rank_1_has_highest_overall_score(self):
        """Rank 1 always has the highest overall_score."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        b1 = _make_bundle(
            mutation_index=0,
            be_results=[_make_be_feasible("ABE")],
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='A', alt='G', position=300),
            mutation_index=1,
            be_results=[_make_be_feasible("CBE", score=0.80)],
            pe_result=_make_pe_feasible(score=0.72),
            hdr_result=_make_hdr_feasible(score=0.60),
        )
        strategies = gen.generate([b1, b2])
        ranked = scorer.rank(strategies, [b1, b2])

        if len(ranked) >= 2:
            assert ranked[0].overall_score >= ranked[1].overall_score
            assert ranked[0].rank == 1

    def test_rejected_excluded_from_ranking(self):
        """Rejected strategies are separated and excluded from the ranked list."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        # One feasible + one infeasible bundle
        b1 = _make_bundle(
            mutation_index=0,
            be_results=[_make_be_feasible("ABE")],
        )
        strategies = gen.generate([b1])

        # Manually add a rejected strategy
        rejected_strat = Strategy(
            name="Manually Rejected",
            rejection_reasons=["Test rejection"],
            modality_prior_score=0.50,
        )
        all_strategies = strategies + [rejected_strat]

        ranked = scorer.rank(all_strategies, [b1])
        ranked_names = [r.strategy_name for r in ranked]
        assert "Manually Rejected" not in ranked_names

    def test_be_ranked_first_over_hdr_with_ipsc_weights(self):
        """In iPSC context, BE strategy should rank above HDR strategy."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE", score=0.85)],
            hdr_result=_make_hdr_feasible(score=0.70),
        )
        strategies = gen.generate([bundle])
        ranked = scorer.rank(strategies, [bundle])

        assert len(ranked) >= 2
        # The top strategy should be BE (DSB-free, higher safety)
        assert ranked[0].strategy_name == "Single-step Base Editing", (
            f"Expected BE as rank 1, got {ranked[0].strategy_name}"
        )

    def test_dual_be_ranked_above_sequential_hdr(self):
        """Dual BE (DSB-free) should rank above Sequential HDR (2 rounds, DSBs)."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        b1 = _make_bundle(
            mutation_index=0,
            be_results=[_make_be_feasible("ABE", score=0.85)],
            hdr_result=_make_hdr_feasible(score=0.70),
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='A', alt='G', position=300),
            mutation_index=1,
            be_results=[_make_be_feasible("CBE", score=0.80)],
            hdr_result=_make_hdr_feasible(score=0.65),
        )
        strategies = gen.generate([b1, b2])
        ranked = scorer.rank(strategies, [b1, b2])

        # Find the Dual BE and Sequential HDR in the ranked list
        dual_be_rank = None
        seq_hdr_rank = None
        for r in ranked:
            if "Dual Base Editing" in r.strategy_name:
                dual_be_rank = r.rank
            if "Sequential HDR" in r.strategy_name:
                seq_hdr_rank = r.rank

        if dual_be_rank is not None and seq_hdr_rank is not None:
            assert dual_be_rank < seq_hdr_rank, (
                f"Dual BE (rank {dual_be_rank}) should rank above "
                f"Sequential HDR (rank {seq_hdr_rank})"
            )

    def test_all_ranked_strategies_have_positive_score(self):
        """All ranked strategies should have overall_score > 0."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE")],
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = gen.generate([bundle])
        ranked = scorer.rank(strategies, [bundle])

        for r in ranked:
            assert r.overall_score > 0, (
                f"Strategy '{r.strategy_name}' has non-positive score: "
                f"{r.overall_score}"
            )

    def test_ranks_are_contiguous_from_1(self):
        """Ranks should be 1, 2, 3, ... with no gaps."""
        gen = StrategyGenerator(p53_active=True)
        scorer = StrategyScorer()

        bundle = _make_bundle(
            be_results=[_make_be_feasible("ABE")],
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = gen.generate([bundle])
        ranked = scorer.rank(strategies, [bundle])

        expected_ranks = list(range(1, len(ranked) + 1))
        actual_ranks = [r.rank for r in ranked]
        assert actual_ranks == expected_ranks, (
            f"Ranks should be contiguous: expected {expected_ranks}, "
            f"got {actual_ranks}"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Prime editing capability gate tests
# ═══════════════════════════════════════════════════════════════════════════


class TestPrimeEditingCapabilityGate:
    """Regression tests for the hard PE capability gate.

    Motivated by the v3 benchmark: HDR-required large-deletion cases
    (DMD, NF1, FBN1) had PE ranked first despite PE being incapable of
    executing edits beyond its empirical ~50 bp single-pegRNA limit.
    The gate rejects PE for such variants instead of silently ranking it.
    """

    @pytest.fixture
    def generator(self):
        return StrategyGenerator(p53_active=True)

    def _variant_with_span(self, span_bp):
        """Build a NormalizedVariant whose structural span is explicit."""
        vi = GenomicVariantInput(
            '1', 100, 'N', '-',
            gene_symbol='TEST',
            structural_span_bp=span_bp,
        )
        ti = TranscriptInfo(
            'T1', 'TEST', 'G1', '1', 1, 1000, 1, 'protein_coding', True,
            [ExonRecord('ENSE001', 1, 1, 1000, 1, '1')],
        )
        tc = TranscriptCoordinate(
            100, 1, 100, 100, 34, 1, 'ATG', 'M', 99, 50, True,
        )
        ca = CodingAnnotation(consequence=ConsequenceType.INFRAME_DELETION)
        rv = ReferenceValidation(True, 'N', 'N', 'N', 'N')
        return NormalizedVariant(
            input=vi, transcript=ti, transcript_coord=tc,
            coding=ca, ref_validation=rv,
        )

    def test_pe_rejected_for_large_structural_deletion(self, generator):
        """Variant with structural_span_bp > 50 → PE strategy is rejected."""
        bundle = _make_bundle(
            variant=self._variant_with_span(70),
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert len(pe_strategies) == 1, (
            "PE strategy should still be emitted (as rejected), "
            "so users see it was considered. "
            f"Got {[s.name for s in strategies]}"
        )
        assert pe_strategies[0].is_rejected, (
            "PE strategy on a large deletion must be marked rejected."
        )
        reason = " ".join(pe_strategies[0].rejection_reasons).lower()
        assert "exceed" in reason or "capability" in reason or "pegrna" in reason, (
            f"Rejection reason should cite the capability limit. "
            f"Got: {pe_strategies[0].rejection_reasons}"
        )

    def test_pe_kept_for_substitution(self, generator):
        """Point substitution (1 bp) → PE strategy is kept (not rejected)."""
        # Use a normal 1-bp substitution variant (no structural_span_bp)
        bundle = _make_bundle(
            variant=_make_variant(ref='C', alt='T'),
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert len(pe_strategies) == 1
        assert not pe_strategies[0].is_rejected, (
            "PE on a 1-bp substitution must not be gated. "
            f"Rejection reasons: {pe_strategies[0].rejection_reasons}"
        )

    def test_pe_kept_at_threshold_boundary(self, generator):
        """structural_span_bp == 50 (== threshold) → PE still kept."""
        bundle = _make_bundle(
            variant=self._variant_with_span(50),
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert len(pe_strategies) == 1
        assert not pe_strategies[0].is_rejected, (
            "PE at the 50-bp boundary must not be gated "
            "(gate fires strictly above 50)."
        )

    def test_pe_rejected_just_above_threshold(self, generator):
        """structural_span_bp == 51 (1 bp above threshold) → PE is rejected."""
        bundle = _make_bundle(
            variant=self._variant_with_span(51),
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert len(pe_strategies) == 1
        assert pe_strategies[0].is_rejected

    def test_pe_rejected_for_long_ref_allele(self, generator):
        """Fallback path: long ref_allele (no structural_span_bp) → gated."""
        vi = GenomicVariantInput(
            '1', 100,
            ref_allele='A' * 100, alt_allele='-',
            gene_symbol='TEST',
        )
        ti = TranscriptInfo(
            'T1', 'TEST', 'G1', '1', 1, 1000, 1, 'protein_coding', True,
            [ExonRecord('ENSE001', 1, 1, 1000, 1, '1')],
        )
        tc = TranscriptCoordinate(
            100, 1, 100, 100, 34, 1, 'ATG', 'M', 99, 50, True,
        )
        ca = CodingAnnotation(consequence=ConsequenceType.INFRAME_DELETION)
        rv = ReferenceValidation(True, 'A' * 100, 'A' * 100, 'A' * 100, 'A' * 100)
        variant = NormalizedVariant(
            input=vi, transcript=ti, transcript_coord=tc,
            coding=ca, ref_validation=rv,
        )
        bundle = _make_bundle(
            variant=variant,
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([bundle])
        pe_strategies = [s for s in strategies if "Prime Editing" in s.name]
        assert pe_strategies and pe_strategies[0].is_rejected, (
            "Long explicit ref_allele should trigger the gate "
            "even without structural_span_bp."
        )

    def test_rejected_pe_is_not_scored_or_ranked(self, generator):
        """Rejected PE strategy should be excluded from ranker output."""
        scorer = StrategyScorer()
        bundle = _make_bundle(
            variant=self._variant_with_span(200),
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([bundle])
        ranked = scorer.rank(strategies, [bundle])
        ranked_names = [r.strategy.name for r in ranked]
        # Rejected PE must not appear in the ranked (scored) list
        assert "Single-step Prime Editing" not in ranked_names, (
            f"Rejected PE must not be scored/ranked. "
            f"Ranked strategies: {ranked_names}"
        )
        # HDR should now be the top-ranked (only viable candidate here)
        assert len(ranked) >= 1
        assert "HDR" in ranked[0].strategy.name, (
            f"With PE gated and no BE, HDR should be top-1. "
            f"Got: {ranked[0].strategy.name}"
        )

    def test_dual_pe_skipped_when_either_variant_gated(self, generator):
        """Two-mutation case: Dual PE is skipped if either variant exceeds gate."""
        b1 = _make_bundle(
            variant=self._variant_with_span(200),  # gated
            mutation_index=0,
            pe_result=_make_pe_feasible(),
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='G', alt='A', position=200),
            mutation_index=1,
            pe_result=_make_pe_feasible(),
        )
        strategies = generator.generate([b1, b2])
        dual_pe = [s for s in strategies if s.name == "Dual Prime Editing"]
        assert len(dual_pe) == 0, (
            "Dual PE must not be generated when one variant is PE-gated."
        )

    def test_hybrid_pe_hdr_skipped_when_pe_target_is_gated(self, generator):
        """Hybrid PE+HDR: skipped when the PE-target variant exceeds gate."""
        b1 = _make_bundle(
            variant=self._variant_with_span(200),  # gated — PE incapable
            mutation_index=0,
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        b2 = _make_bundle(
            variant=_make_variant(ref='G', alt='A', position=200),
            mutation_index=1,
            pe_result=_make_pe_feasible(),
            hdr_result=_make_hdr_feasible(),
        )
        strategies = generator.generate([b1, b2])
        # The hybrid with b1 as PE-target must not appear; the reverse
        # direction (b2 as PE-target, b1 as HDR-target) is still allowed.
        hybrid_pe_hdr = [
            s for s in strategies if s.name == "Hybrid Prime Editing + HDR"
        ]
        for s in hybrid_pe_hdr:
            # All remaining hybrids must have their PE step on b2 (index 1),
            # not on the gated b1 (index 0).
            pe_step = next(
                st for st in s.steps if st.modality == EditModality.PE
            )
            assert pe_step.target_mutation_index == 1, (
                "Hybrid PE+HDR must not assign PE to the gated variant."
            )
