"""Regression tests for Fix #4 (2026-04-20) — compound-het completeness.

Verifies:
    1. ``Strategy.completeness_ratio`` correctly reports the fraction of
       addressed variant indices for single-variant and compound-het cases.
    2. ``StrategyScorer.score_strategy`` applies the safety penalty only
       when ``len(bundles) >= 2`` AND ``completeness_ratio < 1.0``.
    3. ``COMPLETENESS_PENALTY_COEF`` is pinned at 0.30.
    4. ``benchmarks.evaluator._strategy_matches`` rejects ``Single-step X``
       matches against ``hybrid``/``sequential`` truth labels when the case
       has ``n_variants >= 2``.

Reference: COL7A1_AUDIT_2026-04-19.md §7 option B + C. Design rationale
and benchmark implications documented in FIX_NOTES_2026-04-20_fix4.md.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import List

import pytest

from benchmarks.evaluator import _strategy_matches
from core.models import (
    EvidenceTier,
    FeasibilityBundle,
    RiskLevel,
    Strategy,
    StrategyStep,
    EditModality,
)
from core.pipeline.strategy_stage import (
    COMPLETENESS_PENALTY_COEF,
    StrategyScorer,
)


# ── Helpers ────────────────────────────────────────────────────────


def _strategy(
    name: str,
    target_indices: List[int],
    num_dsbs: int = 0,
    bystander: float = 0.0,
) -> Strategy:
    """Minimal Strategy with ``steps`` populated for completeness checks."""
    steps = [
        StrategyStep(
            modality=EditModality.ABE,
            target_mutation_index=i,
        )
        for i in target_indices
    ]
    return Strategy(
        name=name,
        steps=steps,
        num_dsbs=num_dsbs,
        modality_prior_score=0.9,
        donor_feasibility_score=1.0,
        num_rounds=1,
        num_donors=0,
        bystander_severity=bystander,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.A,
        p53_active=False,
        screening_clones=12,
    )


def _bundles(n: int) -> List:
    """Produce n minimal objects that pass ``isinstance(FeasibilityBundle)``
    context-free; scorer paths we exercise read only ``len(bundles)``
    plus per-element attributes on subsets we skip here (splice_proximity)."""
    return [
        SimpleNamespace(
            variant=SimpleNamespace(
                coding=SimpleNamespace(splice_proximity=None)
            )
        )
        for _ in range(n)
    ]


# ── Constant guard ─────────────────────────────────────────────────


def test_completeness_coefficient_is_0_30():
    assert COMPLETENESS_PENALTY_COEF == 0.30


# ── Strategy.completeness_ratio ────────────────────────────────────


class TestCompletenessRatio:
    def test_single_variant_case_returns_1(self):
        s = _strategy("Single", target_indices=[0])
        assert s.completeness_ratio(total_variants=1) == 1.0

    def test_zero_total_variants_returns_1(self):
        s = _strategy("Empty-like", target_indices=[])
        assert s.completeness_ratio(total_variants=0) == 1.0

    def test_empty_steps_in_compound_case_returns_0(self):
        s = _strategy("Empty", target_indices=[])
        assert s.completeness_ratio(total_variants=2) == 0.0

    def test_single_step_in_compound_het_returns_half(self):
        s = _strategy("Single-step BE", target_indices=[0])
        assert s.completeness_ratio(total_variants=2) == 0.5

    def test_two_step_hybrid_in_compound_het_returns_1(self):
        s = _strategy("Hybrid BE+HDR", target_indices=[0, 1])
        assert s.completeness_ratio(total_variants=2) == 1.0

    def test_duplicate_indices_count_as_one(self):
        # A strategy with two steps both targeting index 0 still only
        # addresses one variant in a compound-het case.
        s = _strategy("Double BE", target_indices=[0, 0])
        assert s.completeness_ratio(total_variants=2) == 0.5

    def test_invalid_indices_are_dropped(self):
        s = _strategy("OOB", target_indices=[0, 5, -1])
        assert s.completeness_ratio(total_variants=2) == 0.5


# ── Safety penalty via score_strategy ──────────────────────────────


class TestScoreStrategyCompletenessPenalty:
    scorer = StrategyScorer()

    def test_single_variant_case_no_penalty(self):
        s = _strategy("Single-step BE", target_indices=[0])
        scored = self.scorer.score_strategy(s, _bundles(1))
        # Safety should equal the 0-DSB base score (1.0), not penalised
        assert scored.safety_score == pytest.approx(1.0)

    def test_compound_het_complete_strategy_no_penalty(self):
        s = _strategy("Hybrid BE+HDR", target_indices=[0, 1])
        scored = self.scorer.score_strategy(s, _bundles(2))
        assert scored.safety_score == pytest.approx(1.0)

    def test_compound_het_incomplete_strategy_penalised(self):
        s = _strategy("Single-step BE", target_indices=[0])
        scored = self.scorer.score_strategy(s, _bundles(2))
        # completeness = 0.5 → penalty = (1 - 0.5) * 0.30 = 0.15
        assert scored.safety_score == pytest.approx(1.0 - 0.15)

    def test_compound_het_zero_coverage_maximally_penalised(self):
        s = _strategy("Mystery", target_indices=[])
        scored = self.scorer.score_strategy(s, _bundles(2))
        # completeness = 0 → penalty = 0.30
        assert scored.safety_score == pytest.approx(1.0 - 0.30)

    def test_triple_het_one_of_three_coverage(self):
        s = _strategy("Single-step BE", target_indices=[0])
        scored = self.scorer.score_strategy(s, _bundles(3))
        # completeness = 1/3 → penalty = (2/3) * 0.30 = 0.20
        assert scored.safety_score == pytest.approx(1.0 - 0.20)

    def test_penalty_floors_safety_at_zero(self):
        # Start from a very low base (2-DSB + p53 + simultaneous -> 0.1)
        # then compound-het incomplete (-0.15) would push below 0; clamp.
        from core.models import RiskLevel as R
        steps = [
            StrategyStep(
                modality=EditModality.ABE, target_mutation_index=0
            )
        ]
        s = Strategy(
            name="ugly", steps=steps,
            num_dsbs=2, simultaneous_dsbs=True, p53_active=True,
            num_rounds=1, num_donors=0, bystander_severity=0.0,
            modality_prior_score=0.5, donor_feasibility_score=1.0,
            rearrangement_risk=R.LOW, evidence_tier=EvidenceTier.A,
            screening_clones=12,
        )
        scored = self.scorer.score_strategy(s, _bundles(2))
        assert scored.safety_score >= 0.0


# ── Evaluator tightening (Option B) ────────────────────────────────


class TestEvaluatorCompoundHetAwareMatch:
    def test_single_step_does_not_match_hybrid_truth_in_compound_case(self):
        assert not _strategy_matches(
            "Single-step Base Editing",
            ["hybrid base editing + hdr"],
            n_variants=2,
        )

    def test_single_step_does_not_match_sequential_truth_in_compound_case(self):
        assert not _strategy_matches(
            "Single-step HDR",
            ["sequential hdr"],
            n_variants=2,
        )

    def test_single_step_still_matches_single_step_truth_in_compound_case(self):
        # Compound-het but truth label happens to be "base editing" — match
        assert _strategy_matches(
            "Single-step Base Editing",
            ["base editing"],
            n_variants=2,
        )

    def test_hybrid_strategy_matches_hybrid_truth(self):
        assert _strategy_matches(
            "Hybrid Base Editing + HDR",
            ["hybrid base editing + hdr"],
            n_variants=2,
        )

    def test_sequential_strategy_matches_sequential_truth(self):
        assert _strategy_matches(
            "Sequential HDR",
            ["sequential hdr"],
            n_variants=2,
        )

    def test_single_step_matches_hybrid_in_single_variant_case(self):
        # Edge case: n_variants=1 means no compound-het guard applies
        assert _strategy_matches(
            "Single-step Base Editing",
            ["hybrid base editing + hdr"],
            n_variants=1,
        )

    def test_default_n_variants_preserves_legacy_matching(self):
        # Callers who don't pass n_variants should see unchanged behaviour
        assert _strategy_matches(
            "Single-step Base Editing",
            ["base editing"],
        )
        assert _strategy_matches(
            "Sequential HDR",
            ["sequential hdr"],
        )
