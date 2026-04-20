"""Regression tests for Fix #2 (2026-04-20).

Verifies that:
    1. ``bystander_severity`` penalises the safety score for every DSB
       count, including the 0-DSB case (where it previously returned a
       ceiling of 1.0).
    2. ``_compute_consequence_penalty`` no longer includes the bystander
       term (it now only captures splice proximity).
    3. Clamping at 0 still holds for pathological bystander values.
    4. The published coefficient ``BYSTANDER_SAFETY_COEF`` is 0.08 —
       a pre-registered constant that must not drift post-refactor.

Reference: REVIEW_NOTES_2026-04-18.md §4.2 and FIX_NOTES_2026-04-20.md.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from core.models import (
    EvidenceTier,
    RiskLevel,
    Strategy,
)
from core.pipeline.strategy_stage import (
    BYSTANDER_SAFETY_COEF,
    StrategyScorer,
)


def _strategy(num_dsbs: int = 0, bystander: float = 0.0, p53: bool = False,
              simultaneous: bool = False) -> Strategy:
    """Minimal Strategy builder covering the fields _score_safety reads."""
    return Strategy(
        name="test",
        num_dsbs=num_dsbs,
        modality_prior_score=0.9,
        donor_feasibility_score=1.0,
        num_rounds=1,
        num_donors=0,
        bystander_severity=bystander,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.A,
        p53_active=p53,
        simultaneous_dsbs=simultaneous,
        screening_clones=12,
    )


# ── Constant guard ──────────────────────────────────────────────────


def test_bystander_coefficient_is_0_08():
    """If this fails, someone has tampered with a published constant.

    The 0.08 figure is cited in the review notes and in any paper
    describing v3+ scoring. Changes post-publication must be
    deliberate, justified against literature, and reflected in the
    paper's methods section.
    """
    assert BYSTANDER_SAFETY_COEF == 0.08


# ── Safety score now depends on bystander at every DSB count ────────


class TestSafetyBystanderContribution:
    scorer = StrategyScorer()

    def test_zero_dsb_no_bystander_is_1(self):
        s = _strategy(num_dsbs=0, bystander=0.0)
        assert self.scorer._score_safety(s) == pytest.approx(1.0)

    def test_zero_dsb_with_bystander_drops_below_1(self):
        s = _strategy(num_dsbs=0, bystander=0.5)
        assert self.scorer._score_safety(s) == pytest.approx(1.0 - 0.5 * 0.08)

    def test_zero_dsb_full_bystander(self):
        # bystander_severity is bounded to [0, 1] by convention
        s = _strategy(num_dsbs=0, bystander=1.0)
        assert self.scorer._score_safety(s) == pytest.approx(1.0 - 0.08)

    def test_one_dsb_bystander_stacks_with_p53(self):
        s = _strategy(num_dsbs=1, bystander=0.4, p53=True)
        # base = 0.6; p53 -= 0.1 -> 0.5; bystander 0.4*0.08 = 0.032 -> 0.468
        assert self.scorer._score_safety(s) == pytest.approx(0.5 - 0.4 * 0.08)

    def test_two_dsb_bystander_stacks_with_simultaneous_and_p53(self):
        s = _strategy(num_dsbs=2, bystander=0.2, p53=True, simultaneous=True)
        # base = 0.3; simultaneous -= 0.1 -> 0.2; p53 -= 0.1 -> 0.1; bystander 0.016 -> 0.084
        assert self.scorer._score_safety(s) == pytest.approx(0.1 - 0.2 * 0.08)

    def test_safety_clamped_non_negative(self):
        # Pathological bystander > 1 pushed into a low-base-score branch
        # should still clamp at 0 rather than go negative.
        s = _strategy(num_dsbs=2, bystander=1.0, p53=True, simultaneous=True)
        # base: 0.3 - 0.1 - 0.1 = 0.1; bystander 1.0*0.08 = 0.08; 0.02 >= 0
        assert self.scorer._score_safety(s) == pytest.approx(0.1 - 0.08)

    def test_safety_deep_negative_floor(self):
        # Force deep negative base; must clamp at 0, never negative.
        s = _strategy(num_dsbs=2, bystander=10.0, p53=True, simultaneous=True)
        assert self.scorer._score_safety(s) == 0.0


# ── Consequence penalty no longer includes bystander ────────────────


class TestConsequencePenaltyOnlySpliceProximity:
    scorer = StrategyScorer()

    def test_consequence_penalty_zero_when_no_splice_context(self):
        s = _strategy(num_dsbs=0, bystander=0.5)
        # No feasibility bundles with splice_proximity issues; old code
        # would return 0.5 * 0.08 = 0.04. Post-Fix-#2 returns 0.
        assert self.scorer._compute_consequence_penalty(s, bundles=[]) == 0.0

    def test_consequence_penalty_unaffected_by_bystander_value(self):
        low = _strategy(num_dsbs=0, bystander=0.0)
        high = _strategy(num_dsbs=0, bystander=1.0)
        assert (
            self.scorer._compute_consequence_penalty(low, bundles=[])
            == self.scorer._compute_consequence_penalty(high, bundles=[])
        )

    def test_splice_proximity_still_penalised(self):
        # Splice-proximity penalties stay in _compute_consequence_penalty.
        # Build a minimal bundle-like object with the fields accessed.
        bundle = SimpleNamespace(
            variant=SimpleNamespace(
                coding=SimpleNamespace(splice_proximity="near_exon_start")
            )
        )
        s = _strategy(num_dsbs=0, bystander=0.0)
        assert self.scorer._compute_consequence_penalty(s, bundles=[bundle]) == 0.05


# ── Net architectural property ─────────────────────────────────────


class TestPEvsBEDistinguishableOnSafety:
    """The core value of Fix #2: two 0-DSB modalities can no longer be
    safety-tied when one has bystanders and the other does not."""

    scorer = StrategyScorer()

    def test_pe_beats_be_on_safety_when_be_has_bystanders(self):
        pe = _strategy(num_dsbs=0, bystander=0.0)
        be = _strategy(num_dsbs=0, bystander=0.4)
        assert self.scorer._score_safety(pe) > self.scorer._score_safety(be)

    def test_two_0dsb_modalities_with_equal_bystander_are_safety_tied(self):
        pe = _strategy(num_dsbs=0, bystander=0.3)
        be = _strategy(num_dsbs=0, bystander=0.3)
        assert self.scorer._score_safety(pe) == pytest.approx(
            self.scorer._score_safety(be)
        )
