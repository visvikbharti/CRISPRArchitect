# CRISPRArchitect — Fix Notes #3 (rank stability surfacing)

**Date.** 2026-04-20
**Scope.** Surface Monte-Carlo rank-stability at inference time. Output-only — no scoring changes, no behavioural regression on the benchmark. Complements FIX_NOTES_2026-04-18.md (Fix #1) and FIX_NOTES_2026-04-20.md (Fix #2).
**Reference.** REVIEW_NOTES_2026-04-18.md §4.3.

---

## Motivation (from the external review)

> "Ship rank stability at inference time. Breaks the 'feels redundant' feeling — on flip-sensitive cases the tool is answering a question the user couldn't answer without it."

Before Fix #3, the CLI reported a single ``rank stability`` percentage for the top strategy (e.g. ``"58.0% (top-ranked in 58.0% of 10,000 weight permutations)"``). The number was present, but the user had no way to know whether 58% was "good" or "bad," what the runner-up modality would be, or what specific tradeoff the top vs. runner-up choice represented.

Fix #3 interprets the number for the user and, when the top is flip-sensitive, surfaces:
1. The runner-up modality and its own rank-stability.
2. The TOPSIS score gap between top and runner-up.
3. A per-dimension comparison (which dimensions each strategy wins on, ignoring dimensions where the delta is sub-threshold).
4. Human-readable reasoning lines for each strategy.

## Pre-registered thresholds

Three levels, mapping ``rank_stability`` to a categorical assessment:

| Level | Range | Action |
|---|---|---|
| ``robust`` | ≥ 0.80 | Unconditional recommendation |
| ``stable`` | 0.70 ≤ … < 0.80 | Single recommendation, note uncertainty |
| ``flip_sensitive`` | < 0.70 | Surface runner-up + dimension-level tradeoff |
| ``unknown`` | ``None`` | Sensitivity analysis disabled |

The thresholds are pre-committed module-level constants:

```python
RANK_STABILITY_ROBUST = 0.80
RANK_STABILITY_STABLE = 0.70
DIMENSION_DELTA_MATERIAL = 0.05  # below this, a dimension delta is a tie
```

A regression test (``test_published_thresholds_unchanged``) guards against post-publication tampering.

## What changed

### Files modified

| File | Change |
|---|---|
| ``core/pipeline/strategy_stage.py`` | Added module constants ``RANK_STABILITY_ROBUST``, ``RANK_STABILITY_STABLE``, ``DIMENSION_DELTA_MATERIAL``. Added ``StabilityAssessment`` dataclass, ``DimensionDelta`` dataclass, ``_classify_stability_level``, ``_dimension_comparison``, ``_preferential_reasoning``, and ``assess_stability`` public function. |
| ``cli.py`` | ``cmd_analyze_v3`` (``_analyze_v3`` reporting block) now calls ``assess_stability`` and prints the richer recommendation block instead of a bare percentage. |
| ``tests/test_fix3_stability_assessment.py`` | New file, 24 regression tests: threshold boundaries, benefit-vs-cost dimension direction handling, tie handling (< DIMENSION_DELTA_MATERIAL), robust/stable/flip-sensitive assessment construction, runner-up surfacing, and empty/single-strategy edge cases. |

### Example output — robust case

```
Recommendation: Single-step Prime Editing
  Evidence tier: A
  Rank stability: 92.0% — ROBUST (top-ranked in 92.0% of 10,000 weight permutations)
    Above 80% threshold: unconditional recommendation.
```

### Example output — stable case

```
Recommendation: Single-step Prime Editing
  Evidence tier: A
  Rank stability: 75.0% — STABLE (top-ranked in 75.0% of 10,000 weight permutations)
    In [70%, 80%): single recommendation but weight perturbations sometimes promote alternatives.
```

### Example output — flip-sensitive case (the value-add)

```
Recommendation: Single-step Prime Editing
  Evidence tier: A
  Rank stability: 58.0% — FLIP-SENSITIVE (top-ranked in 58.0% of 10,000 weight permutations)
    Below 70%: the choice between top and alternative depends on context. Consider the tradeoff below.

  Alternative considered: Single-step Base Editing
    Alternative stability: 30.0%
    TOPSIS score gap: +0.020 (top = 0.710, alt = 0.690)

  Why each might be preferred:
    - Single-step Prime Editing wins on: safety (1.00 vs 0.85), confidence (1.00 vs 0.70)
    - Single-step Base Editing wins on: feasibility (0.95 vs 0.82)
```

This is the case the tool is *for*. A user looking at this output can see concretely what tradeoff they are accepting by taking the top recommendation, and can adjust their weights or choose the alternative intentionally.

## Benchmark delta

Output-only refactor. The scoring pipeline is untouched: ``assess_stability`` reads a ranked ``List[ScoredStrategy]`` and returns an interpretation object. It does not modify any score, rank, or flag.

- 269/269 tests pass (245 pre-Fix-#3 + 24 new).
- No benchmark rerun needed — the ranker, sensitivity analysis, and per-case top-1 recommendations are all identical by construction.

## Paragraph for the paper

Proposed methods-section text:

> **Rank stability.** Every recommendation is accompanied by a rank-stability score: the fraction of 10,000 Monte-Carlo weight perturbations under which the top-ranked strategy retains rank 1. When stability is ≥80%, the recommendation is reported as unconditional. When stability is in [70%, 80%), the single recommendation stands but we note that weight perturbations occasionally promote alternatives. When stability is <70%, we classify the case as flip-sensitive: the CLI surfaces the runner-up modality, the TOPSIS-score gap, and a per-dimension comparison of where each candidate wins, so that the user can choose intentionally on their context. This behavior is implemented by ``assess_stability`` in ``core/pipeline/strategy_stage.py`` with thresholds ``RANK_STABILITY_ROBUST = 0.80`` and ``RANK_STABILITY_STABLE = 0.70``.

## Still not fixed

- **Compound-het completeness penalty.** Per COL7A1_AUDIT_2026-04-19.md §7 option C: a strategy that addresses only 1 of 2 pathogenic variants in a compound-heterozygous case should be penalised for incompleteness. Distinct architectural gap; not addressed by any of Fix #1-#3.
- **Webapp surfacing.** Fix #3 updates the CLI. The Streamlit v3 page (``webapp/app_v3_page.py``) still shows the older rank-stability block. Porting ``assess_stability`` into the webapp is a ~half-day follow-up.

## Commit

- ``core/pipeline/strategy_stage.py``: constants + assess_stability infrastructure.
- ``cli.py``: use assess_stability in the v3 analyze reporting block.
- ``tests/test_fix3_stability_assessment.py``: 24 new regression tests.
- ``FIX_NOTES_2026-04-20_fix3.md``: this document.

---

*Fix #3 closes Review §4.3. Compound-het completeness and webapp parity remain open.*
