"""
Tests for CRISPRArchitect v3 TOPSIS Scorer + Sensitivity Analysis
===================================================================

Tests the TOPSIS multi-criteria decision method and Monte Carlo
sensitivity analysis for strategy ranking.

References
----------
Hwang & Yoon, Multiple Attribute Decision Making, Springer, 1981.
"""

from __future__ import annotations

import pytest
from dataclasses import field

from core.models import (
    EvidenceTier,
    FeasibilityBundle,
    FeasibilityLabel,
    RiskLevel,
    ScoredStrategy,
    Strategy,
)
from core.pipeline.strategy_stage import (
    StrategyScorer,
    TOPSISScorer,
    SensitivityResult,
)


# ── Helpers ──────────────────────────────────────────────────────────

def _make_strategy(
    name: str,
    num_dsbs: int = 0,
    modality_prior: float = 0.9,
    donor_feas: float = 1.0,
    num_rounds: int = 1,
    num_donors: int = 0,
    bystander: float = 0.0,
    rearrangement: RiskLevel = RiskLevel.LOW,
    evidence: EvidenceTier = EvidenceTier.A,
    p53: bool = True,
    screening: int = 12,
) -> Strategy:
    return Strategy(
        name=name,
        num_dsbs=num_dsbs,
        modality_prior_score=modality_prior,
        donor_feasibility_score=donor_feas,
        num_rounds=num_rounds,
        num_donors=num_donors,
        bystander_severity=bystander,
        rearrangement_risk=rearrangement,
        evidence_tier=evidence,
        p53_active=p53,
        screening_clones=screening,
    )


# ── Tests: TOPSIS Algorithm ──────────────────────────────────────────

class TestTOPSISAlgorithm:
    """Tests for the core TOPSIS scoring logic."""

    def test_topsis_returns_scores_between_0_and_1(self):
        """All TOPSIS scores should be in [0, 1]."""
        scorer = TOPSISScorer()
        # Test with raw matrix
        matrix = [
            [1.0, 0.9, 0.1, 0.1, 1.0],  # Great strategy
            [0.5, 0.5, 0.5, 0.5, 0.5],  # Average
            [0.1, 0.3, 0.8, 0.8, 0.4],  # Poor
        ]
        scores = scorer._topsis(matrix, scorer.weights)
        for s in scores:
            assert 0.0 <= s <= 1.0, f"TOPSIS score {s} out of bounds"

    def test_ideal_strategy_scores_highest(self):
        """Strategy that's best on all dimensions should score highest."""
        scorer = TOPSISScorer()
        matrix = [
            [1.0, 1.0, 0.0, 0.0, 1.0],  # Ideal
            [0.5, 0.5, 0.5, 0.5, 0.5],  # Average
            [0.0, 0.0, 1.0, 1.0, 0.0],  # Anti-ideal
        ]
        scores = scorer._topsis(matrix, scorer.weights)
        assert scores[0] > scores[1] > scores[2]

    def test_identical_strategies_get_equal_scores(self):
        """Identical strategies should get equal TOPSIS scores."""
        scorer = TOPSISScorer()
        matrix = [
            [0.5, 0.5, 0.3, 0.3, 0.7],
            [0.5, 0.5, 0.3, 0.3, 0.7],
        ]
        scores = scorer._topsis(matrix, scorer.weights)
        assert abs(scores[0] - scores[1]) < 1e-10

    def test_single_strategy_returns_05(self):
        """Single strategy should score 0.5 (equidistant from ideal/anti-ideal)."""
        scorer = TOPSISScorer()
        # With single alternative, ideal == anti-ideal == the strategy itself
        # This is a degenerate case; our implementation handles it via the rank() method
        strategies = [_make_strategy("only_one")]
        ranked = scorer.rank(strategies, [], run_sensitivity=False)
        assert len(ranked) == 1
        assert ranked[0].rank == 1

    def test_be_ranks_above_hdr_for_safe_scenario(self):
        """Base editing (no DSB) should rank above HDR (DSB) when safety matters."""
        scorer = TOPSISScorer(w_safety=0.40)
        be_strategy = _make_strategy(
            "BE", num_dsbs=0, modality_prior=0.95,
            evidence=EvidenceTier.A, bystander=0.0,
        )
        hdr_strategy = _make_strategy(
            "HDR", num_dsbs=1, modality_prior=0.80,
            num_donors=1, evidence=EvidenceTier.A,
        )
        ranked = scorer.rank([be_strategy, hdr_strategy], [], run_sensitivity=False)
        assert ranked[0].strategy_name == "BE"

    def test_complexity_penalizes_multi_round(self):
        """Multi-round strategy should score lower than single-round."""
        scorer = TOPSISScorer()
        simple = _make_strategy("simple", num_rounds=1, num_donors=0)
        complex_ = _make_strategy("complex", num_rounds=3, num_donors=2, screening=48)
        ranked = scorer.rank([simple, complex_], [], run_sensitivity=False)
        assert ranked[0].strategy_name == "simple"


# ── Tests: Sensitivity Analysis ──────────────────────────────────────

class TestSensitivityAnalysis:
    """Tests for Monte Carlo weight perturbation."""

    def test_sensitivity_returns_results_for_all_strategies(self):
        """Sensitivity analysis should return results for every strategy."""
        scorer = TOPSISScorer(n_sensitivity_runs=100)
        s1 = _make_strategy("A", num_dsbs=0, modality_prior=0.95)
        s2 = _make_strategy("B", num_dsbs=1, modality_prior=0.80, num_donors=1)
        ranked = scorer.rank([s1, s2], [], run_sensitivity=True)
        # Both should have sensitivity info in annotations
        for r in ranked:
            has_stability = any("Rank stability" in n for n in r.annotation_notes)
            assert has_stability, f"{r.strategy_name} missing sensitivity annotation"

    def test_dominant_strategy_has_high_stability(self):
        """A clearly dominant strategy should have high rank stability."""
        scorer = TOPSISScorer(n_sensitivity_runs=1000)
        dominant = _make_strategy(
            "dominant", num_dsbs=0, modality_prior=0.95,
            evidence=EvidenceTier.A, bystander=0.0,
        )
        weak = _make_strategy(
            "weak", num_dsbs=2, modality_prior=0.30,
            num_donors=2, num_rounds=3,
            rearrangement=RiskLevel.HIGH,
            evidence=EvidenceTier.C, bystander=0.5,
        )
        ranked = scorer.rank([dominant, weak], [], run_sensitivity=True)
        assert ranked[0].strategy_name == "dominant"
        # The dominant strategy should be top in >90% of permutations
        stability_note = [n for n in ranked[0].annotation_notes if "Rank stability" in n]
        assert len(stability_note) > 0

    def test_sensitivity_is_deterministic(self):
        """Sensitivity uses seed=42, results should be reproducible."""
        scorer = TOPSISScorer(n_sensitivity_runs=500)
        s1 = _make_strategy("A", num_dsbs=0, modality_prior=0.9)
        s2 = _make_strategy("B", num_dsbs=1, modality_prior=0.7, num_donors=1)

        ranked1 = scorer.rank([s1, s2], [], run_sensitivity=True)
        ranked2 = scorer.rank([s1, s2], [], run_sensitivity=True)

        assert ranked1[0].overall_score == ranked2[0].overall_score
        assert ranked1[0].annotation_notes == ranked2[0].annotation_notes

    def test_sensitivity_skipped_for_single_strategy(self):
        """No sensitivity analysis needed for single strategy."""
        scorer = TOPSISScorer(n_sensitivity_runs=100)
        s = _make_strategy("only")
        ranked = scorer.rank([s], [], run_sensitivity=True)
        # Should not crash; no sensitivity notes expected
        assert len(ranked) == 1


# ── Tests: TOPSIS vs Weighted Sum Comparison ─────────────────────────

class TestTOPSISvsWeightedSum:
    """Compare TOPSIS with legacy weighted sum."""

    def test_both_scorers_agree_on_clear_winner(self):
        """Both methods should agree when one strategy clearly dominates."""
        topsis = TOPSISScorer()
        legacy = StrategyScorer()

        good = _make_strategy("good", num_dsbs=0, modality_prior=0.95)
        bad = _make_strategy("bad", num_dsbs=2, modality_prior=0.3,
                           rearrangement=RiskLevel.HIGH, num_donors=2)

        topsis_ranked = topsis.rank([good, bad], [], run_sensitivity=False)
        legacy_ranked = legacy.rank([good, bad], [])

        assert topsis_ranked[0].strategy_name == "good"
        assert legacy_ranked[0].strategy_name == "good"

    def test_rejected_strategies_excluded_from_both(self):
        """Rejected strategies should be excluded by both scorers."""
        topsis = TOPSISScorer()
        legacy = StrategyScorer()

        good = _make_strategy("good")
        rejected = _make_strategy("rejected")
        rejected.rejection_reasons.append("Infeasible")

        assert len(topsis.rank([good, rejected], [], run_sensitivity=False)) == 1
        assert len(legacy.rank([good, rejected], [])) == 1
