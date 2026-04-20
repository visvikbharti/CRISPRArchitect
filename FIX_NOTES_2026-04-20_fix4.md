# CRISPRArchitect — Fix Notes #4 (compound-het completeness penalty)

**Date.** 2026-04-20
**Scope.** Close Options B + C of COL7A1_AUDIT_2026-04-19.md §7: penalise strategies that address fewer than all pathogenic variants in a compound-heterozygous case (scorer-side), and prevent the benchmark evaluator from counting "Single-step X" outputs as correct against "Hybrid" / "Sequential" truth labels when the case is compound-het.
**Reference.** COL7A1_AUDIT_2026-04-19.md §7 options B and C. Complements FIX_NOTES_2026-04-18.md (Fix #1), FIX_NOTES_2026-04-20.md (Fix #2), FIX_NOTES_2026-04-20_fix3.md (Fix #3).

---

## Motivation

The v3 benchmark contains one compound-heterozygous case, `HDR_COL7A1_020`, with two distinct pathogenic 1-bp deletions on separate alleles. Its intended truth label is `"hybrid base editing + hdr"` — a multi-step strategy that addresses both alleles. The pipeline produces a flat list of per-variant strategies (one strategy per pathogenic allele, ranked globally), and for this case the top-1 output is `"Single-step Base Editing"` targeting a single allele.

Two things were wrong with the status quo:

1. **Scorer-side:** the pipeline did not distinguish between a single-variant strategy and a multi-step strategy that addresses all variants. A biologically incomplete output (patient still has one diseased allele after the single-step BE) ranked above the complete hybrid alternative in compound-het cases.
2. **Evaluator-side:** `benchmarks/evaluator.py::_strategy_matches` used a loose substring test that credited `"Single-step Base Editing"` as matching the truth label `"hybrid base editing + hdr"` because both strings contain `"base"`. This silently turned a biological miss into a benchmark-correct case.

The COL7A1 audit (2026-04-19) recommended two separable remediations:

- **Option B:** Tighten `_strategy_matches` to reject single-step matches against hybrid/sequential truth labels when the case is compound-het.
- **Option C:** Apply a completeness penalty in the scorer so that incomplete strategies rank below complete ones.

Fix #4 implements **both**.

## What changed

### Scorer-side (Option C)

- `core/models.py`. `Strategy` gains two helper methods:
  - `addressed_variant_indices() -> set[int]`: the set of distinct `target_mutation_index` values across all `steps`.
  - `completeness_ratio(total_variants: int) -> float`: the fraction of total variants this strategy addresses. Returns 1.0 for `total_variants <= 1` (single-variant cases have no completeness issue). Drops out-of-range and non-integer indices defensively.
- `core/pipeline/strategy_stage.py`. Module-level constant `COMPLETENESS_PENALTY_COEF = 0.30`. In `StrategyScorer.score_strategy`, when `len(bundles) > 1` and `completeness_ratio < 1.0`, apply `safety -= (1 - completeness) * 0.30`, clamped to `[0, 1]`. Applied post-`_score_safety` to keep the safety unit-test interface narrow.

### Evaluator-side (Option B)

- `benchmarks/evaluator.py`. `_strategy_matches` gains an `n_variants: int = 1` argument. When `n_variants >= 2` AND the original (pre-normalisation) strategy name starts with `"single-step"`, any truth label containing `"hybrid"` or `"sequential"` is skipped for matching. The existing default behaviour is preserved for single-variant cases and for non-compound truth labels.
- The three call sites in `evaluate_case` now pass `n_variants=len(case.variants)`.

### Regression tests

`tests/test_fix4_compound_het_completeness.py` (21 tests):

- Constant-guard for `COMPLETENESS_PENALTY_COEF == 0.30`.
- `completeness_ratio`: single-variant (= 1.0), compound-het half-coverage (= 0.5), full hybrid (= 1.0), empty steps (= 0.0), duplicate indices (counted once), out-of-range indices (dropped).
- `score_strategy` safety penalty: no penalty in single-variant case, no penalty for complete hybrid, -0.15 for single-step in compound-het (N=2), -0.20 for single-step in triple-het (N=3), -0.30 for zero-coverage, and clamping at 0 for pathological combinations.
- Evaluator matching: single-step rejected against hybrid/sequential truth in compound cases, still accepted against single-step truth in compound cases, hybrid strategies still match hybrid truth, default call (without `n_variants`) preserves legacy matching.

## Pre-committed penalty coefficient

```python
COMPLETENESS_PENALTY_COEF = 0.30
```

The coefficient is chosen to meaningfully separate incomplete single-step strategies from complete hybrids without entirely vetoing them:

- A half-coverage strategy (1 of 2 variants) loses 0.15 off its safety score.
- A quarter-coverage strategy (1 of 4 variants) loses ~0.225.
- A hybrid/complete strategy (all variants) is unpenalised.

A *hard* rejection (completeness < 1.0 → is_rejected = True) was considered but rejected: there are legitimate staged-therapy contexts where addressing variants sequentially is the clinical plan. The penalty lets the user override via weight-vector customisation while being honest by default.

## Benchmark delta

Pre-Fix-#4 baseline (post-Fix-#1 through Fix-#3): 30/30 top-1, 30/30 top-3, 29/30 rejection.

**Benchmark re-run status.** A full benchmark execution was attempted at 2026-04-20 20:26 but the Ensembl REST API returned HTTP 500s continuously during the run, exhausting the built-in retry schedule after approximately 9 minutes without producing complete results. The failure is entirely external (other workflows that query Ensembl via the same client were also affected). No code change in Fix #4 affects the transcript-fetch path. The benchmark re-run should be re-attempted once Ensembl stabilises and the full JSON cached at ``benchmark_results/v3_post_fix4_results.json``.

**Expected outcomes for `HDR_COL7A1_020`** (both are correctness improvements over the pre-Fix-#4 evaluator-credit):

1. **Pipeline elevates a hybrid strategy to top-1.** If the generator's `_generate_two_mutation_strategies` produces a hybrid (e.g., `"Hybrid Base Editing + HDR"`) for this case, the safety penalty on `"Single-step Base Editing"` (-0.15) should be enough to flip the ordering. Benchmark remains 30/30 top-1, with a biologically-complete recommendation.
2. **Pipeline does not generate a hybrid that matches.** `"Single-step Base Editing"` remains top-1 but now fails the evaluator match (which no longer credits it against `"hybrid base editing + hdr"` in a compound-het case). Benchmark drops to 29/30 top-1 — the honest number per the COL7A1 audit.

Either outcome is a correctness improvement over the silent evaluator-credit from the pre-Fix-#4 state. The Fix #4 unit tests (21/21 passing) verify the scoring + evaluator logic independently of whether Ensembl is reachable.

## Paragraph for the paper

Proposed methods-section addition:

> **Compound-heterozygous completeness.** When a case has two or more pathogenic variants (e.g., compound-heterozygous inheritance), each strategy is assessed for *completeness*: the fraction of pathogenic variants its steps address. Strategies that address fewer than all variants incur a safety penalty of ``(1 - completeness) × COMPLETENESS_PENALTY_COEF``, with ``COMPLETENESS_PENALTY_COEF = 0.30``. This distinguishes a biologically incomplete single-step recommendation in a compound-het case from a hybrid or sequential strategy that resolves all alleles. The penalty is calibrated to rank complete hybrids above incomplete single-steps in typical benchmark regimes without entirely vetoing the single-step option (which may be appropriate for staged-therapy plans).

> **Compound-het evaluator guard.** For benchmark evaluation, a "Single-step X" strategy name is not accepted as matching a "Hybrid" or "Sequential" truth label when the case contains ≥2 pathogenic variants. This prevents a substring-level match from silently crediting an incomplete recommendation against a multi-step truth.

## Still not fixed

- **Webapp parity.** The Streamlit v3 page (`webapp/app_v3_page.py`) does not yet expose the completeness penalty in its output. A ~0.5-day follow-up will port `assess_stability` from Fix #3 and add completeness information alongside.
- **Generator improvements for compound-het hybrids.** If the benchmark turns out to regress to 29/30 after Fix #4 (outcome 2 above), the follow-up is generator work — make `_generate_two_mutation_strategies` produce a `"Hybrid Base Editing + HDR"` strategy for the COL7A1 compound-het class specifically. Scoped to a separate change; Fix #4 is the scorer-and-evaluator half of the correction.

## Commit

- `core/models.py`: `Strategy.addressed_variant_indices()` + `completeness_ratio(total_variants)`.
- `core/pipeline/strategy_stage.py`: `COMPLETENESS_PENALTY_COEF = 0.30`; completeness-penalty block in `score_strategy`.
- `benchmarks/evaluator.py`: `_strategy_matches(n_variants=1)` parameter + compound-het guard; call-site updates in `evaluate_case`.
- `tests/test_fix4_compound_het_completeness.py`: 21 new regression tests.
- `FIX_NOTES_2026-04-20_fix4.md`: this document.

---

*Fix #4 closes Options B + C of COL7A1_AUDIT §7. Generator-side hybrid generation for compound-het cases and webapp parity remain open.*
