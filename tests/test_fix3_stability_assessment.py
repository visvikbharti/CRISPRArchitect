"""Regression tests for Fix #3 (2026-04-20) — rank stability surfacing.

Verifies:
    1. ``assess_stability`` classifies into
       {"robust", "stable", "flip_sensitive", "unknown"} per the
       pre-committed thresholds at REVIEW_NOTES_2026-04-18.md §4.3.
    2. Flip-sensitive assessments include the runner-up, a score gap,
       per-dimension deltas, and human-readable preferential reasoning.
    3. Dimension-delta direction handling: benefit dimensions prefer the
       higher value; cost dimensions prefer the lower value.
    4. Threshold boundaries behave as specified (0.80 is robust; 0.70 is
       stable; 0.699 is flip-sensitive).
    5. Constants are frozen at their published values.
"""

from __future__ import annotations

import pytest

from core.models import EvidenceTier, RiskLevel, ScoredStrategy, Strategy
from core.pipeline.strategy_stage import (
    DIMENSION_DELTA_MATERIAL,
    RANK_STABILITY_ROBUST,
    RANK_STABILITY_STABLE,
    DimensionDelta,
    StabilityAssessment,
    _classify_stability_level,
    _dimension_comparison,
    assess_stability,
)


def _scored(
    name: str,
    overall: float = 0.7,
    stability: float = None,
    safety: float = 0.9,
    feasibility: float = 0.85,
    complexity: float = 0.2,
    risk: float = 0.1,
    confidence: float = 1.0,
) -> ScoredStrategy:
    strat = Strategy(
        name=name,
        num_dsbs=0,
        modality_prior_score=0.9,
        donor_feasibility_score=1.0,
        num_rounds=1,
        num_donors=0,
        bystander_severity=0.0,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.A,
        p53_active=False,
        screening_clones=12,
    )
    return ScoredStrategy(
        strategy=strat,
        safety_score=safety,
        feasibility_score=feasibility,
        complexity_score=complexity,
        risk_score=risk,
        confidence_score=confidence,
        overall_score=overall,
        rank_stability=stability,
    )


# ── Constant guard ──────────────────────────────────────────────────


def test_published_thresholds_unchanged():
    """Tampering with these without a re-review breaks the paper's claim."""
    assert RANK_STABILITY_ROBUST == 0.80
    assert RANK_STABILITY_STABLE == 0.70
    assert DIMENSION_DELTA_MATERIAL == 0.05


# ── Stability-level classification ──────────────────────────────────


@pytest.mark.parametrize(
    "stability, expected_level",
    [
        (1.00, "robust"),
        (0.85, "robust"),
        (0.80, "robust"),       # boundary: inclusive
        (0.799, "stable"),
        (0.75, "stable"),
        (0.70, "stable"),        # boundary: inclusive
        (0.699, "flip_sensitive"),
        (0.50, "flip_sensitive"),
        (0.00, "flip_sensitive"),
        (None, "unknown"),
    ],
)
def test_stability_level_classification(stability, expected_level):
    assert _classify_stability_level(stability) == expected_level


# ── Dimension deltas ───────────────────────────────────────────────


def test_dimension_comparison_benefit_vs_cost_directions():
    top = _scored("Top", safety=1.0, feasibility=0.9, complexity=0.3, risk=0.1, confidence=1.0)
    ru = _scored("Alt", safety=0.85, feasibility=0.95, complexity=0.1, risk=0.2, confidence=0.7)
    deltas = {d.name: d for d in _dimension_comparison(top, ru)}

    # safety: benefit, top wins (1.0 > 0.85, delta 0.15 ≥ 0.05)
    assert deltas["safety"].prefers == "top"
    # feasibility: benefit, alt wins (0.95 > 0.9, delta -0.05 not material — tie!)
    # Delta is exactly -0.05, which is NOT material (strict < 0.05 test on abs)
    assert deltas["feasibility"].prefers is None
    # complexity: cost (lower better), alt wins (0.1 < 0.3, delta 0.2 material)
    assert deltas["complexity"].prefers == "runner_up"
    # risk: cost, top wins (0.1 < 0.2, delta -0.1)
    assert deltas["risk"].prefers == "top"
    # confidence: benefit, top wins (1.0 vs 0.7)
    assert deltas["confidence"].prefers == "top"


def test_dimension_comparison_ties_surface_as_none():
    top = _scored("Top", safety=0.90, feasibility=0.90, complexity=0.20, risk=0.10, confidence=1.0)
    ru = _scored("Alt", safety=0.92, feasibility=0.88, complexity=0.22, risk=0.12, confidence=1.0)
    # All differences < 0.05 → all ties
    deltas = _dimension_comparison(top, ru)
    for d in deltas:
        assert d.prefers is None, f"{d.name}: expected tie, got prefers={d.prefers}"


# ── Full assessment: robust / stable paths ─────────────────────────


def test_robust_assessment_has_no_runner_up():
    strategies = [
        _scored("PE", overall=0.80, stability=0.92),
        _scored("BE", overall=0.55, stability=0.05),
    ]
    a = assess_stability(strategies)
    assert a is not None
    assert a.level == "robust"
    assert a.runner_up is None
    assert a.score_gap is None
    assert a.dimension_deltas == []
    assert a.preferential_reasoning == []
    assert a.human_label == "ROBUST"


def test_stable_assessment_has_no_runner_up():
    strategies = [
        _scored("PE", overall=0.72, stability=0.75),
        _scored("BE", overall=0.65, stability=0.22),
    ]
    a = assess_stability(strategies)
    assert a.level == "stable"
    assert a.runner_up is None
    assert not a.is_flip_sensitive


def test_unknown_assessment_when_stability_not_run():
    strategies = [_scored("only", overall=0.8, stability=None)]
    a = assess_stability(strategies)
    assert a.level == "unknown"
    assert a.runner_up is None


# ── Full assessment: flip-sensitive path ───────────────────────────


def test_flip_sensitive_assessment_populates_runner_up():
    strategies = [
        _scored("PE", overall=0.71, stability=0.58,
                safety=1.00, feasibility=0.82, complexity=0.25, risk=0.05, confidence=1.0),
        _scored("BE", overall=0.69, stability=0.30,
                safety=0.90, feasibility=0.95, complexity=0.15, risk=0.10, confidence=0.7),
    ]
    a = assess_stability(strategies)
    assert a.is_flip_sensitive
    assert a.runner_up is not None
    assert a.runner_up.strategy_name == "BE"
    assert a.score_gap == pytest.approx(0.71 - 0.69, abs=1e-9)
    # 5 dimensions compared
    assert len(a.dimension_deltas) == 5
    # Reasoning mentions both
    reasoning_text = "\n".join(a.preferential_reasoning)
    assert "PE" in reasoning_text
    assert "BE" in reasoning_text


def test_flip_sensitive_preferential_reasoning_surfaces_safety_win_for_top():
    strategies = [
        _scored("PE", overall=0.71, stability=0.58,
                safety=1.00, feasibility=0.82, complexity=0.25, risk=0.10, confidence=1.0),
        _scored("BE", overall=0.69, stability=0.30,
                safety=0.85, feasibility=0.95, complexity=0.20, risk=0.10, confidence=0.7),
    ]
    a = assess_stability(strategies)
    # Top (PE) should win on safety + confidence; runner_up (BE) on feasibility
    lines = "\n".join(a.preferential_reasoning)
    assert "safety" in lines
    # At least one line should attribute wins to each side
    assert "PE wins" in lines
    assert "BE wins" in lines


def test_flip_sensitive_with_no_material_differences_says_either_defensible():
    strategies = [
        _scored("A", overall=0.71, stability=0.55,
                safety=0.90, feasibility=0.90, complexity=0.20, risk=0.10, confidence=1.0),
        _scored("B", overall=0.70, stability=0.45,
                safety=0.92, feasibility=0.88, complexity=0.22, risk=0.12, confidence=1.0),
    ]
    a = assess_stability(strategies)
    assert a.is_flip_sensitive
    lines = "\n".join(a.preferential_reasoning)
    assert "No material dimension differences" in lines


def test_flip_sensitive_single_strategy_degrades_gracefully():
    # Only one strategy available → can't surface a runner-up
    strategies = [_scored("solo", overall=0.5, stability=0.40)]
    a = assess_stability(strategies)
    # level is still computed, but no runner-up fields
    assert a.level == "flip_sensitive"
    assert a.runner_up is None


# ── Edge cases ─────────────────────────────────────────────────────


def test_empty_strategy_list_returns_none():
    assert assess_stability([]) is None


def test_boundary_exact_0_80_is_robust():
    s = [_scored("t", stability=RANK_STABILITY_ROBUST)]
    a = assess_stability(s)
    assert a.level == "robust"


def test_boundary_exact_0_70_is_stable():
    s = [_scored("t", stability=RANK_STABILITY_STABLE)]
    a = assess_stability(s)
    assert a.level == "stable"


def test_boundary_just_below_0_70_is_flip_sensitive():
    s = [_scored("t", stability=RANK_STABILITY_STABLE - 0.001)]
    a = assess_stability(s)
    assert a.level == "flip_sensitive"
