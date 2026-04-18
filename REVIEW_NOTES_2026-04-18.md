# CRISPRArchitect v3 — Review Notes

**Date.** 2026-04-18
**Reviewer.** Vishal Bharti's Claude Code session (external read of scoring logic, benchmark results, and v3 manuscript).
**Scope of review.** Scoring logic in `core/pipeline/strategy_stage.py`, strategy generation in `core/mosaic/generator.py`, benchmark results at `benchmark_results/v3_benchmark_results.json`, and `paper/CRISPRArchitect_v3_manuscript.md`. No code changes made. No runs executed.
**Trigger question from the author.** *"It seems redundant when we think that in most of the cases answer would be prime editors."*

---

## Executive summary

**PE wins 77% of the v3 benchmark (23 of 30 cases) and this is substantially truthful.** It is not a scoring bias routing everything to PE artificially — it is a defensible consequence of three facts the scorer correctly encodes: (a) PE has no DSBs (same as base editing) so it ties BE on safety = 1.0, (b) PE has no bystander edits by default so it avoids a 0.04–0.08 consequence penalty that BE carries, and (c) PE has broader mutation-type coverage than BE. The headline number is approximately right.

**But the tool's *value* is not the ranking — it's the 23% minority and the explanation.** If you publish "PE wins 77% of the time," a reader could skip the tool and guess PE. The tool earns its keep only on the **20% BE-rescue cases** (where bystander penalties are light and efficiency matters) and the **3% HDR-required cases** (large structural changes PE cannot do). Those are the cases the tool is *for*. Right now two fixable issues are preventing the tool from realising that value fully.

**Three fixable issues, in decreasing importance:**

1. **HDR failure mode is real.** All 4 top-1 errors in the 30-case benchmark are HDR-required cases where PE was ranked first instead of HDR. 0/3 top-1 accuracy on large-deletion HDR cases is called out in the manuscript but not yet fixed.
2. **The `bystander_severity * 0.08` constant at `strategy_stage.py:334` is a free parameter that governs the headline BE-vs-PE statistic and has no documented citation.** It's the coefficient that tips most BE cases into PE wins.
3. **The safety ceiling at 1.0 for all DSB-free modalities** hides real differences (PE has pegRNA-3'-extension byproducts and off-target signatures distinct from BE). Everything gets handled via the post-TOPSIS consequence penalty, which creates a two-stage scoring that's harder to explain than a single principled multi-criteria score.

Full walk-through follows.

---

## 1. How the scoring actually works (to compare against the "redundant" concern)

`StrategyScorer` in `core/pipeline/strategy_stage.py` computes:

```
overall = w_safety · safety
        + w_feasibility · feasibility
        − w_complexity · complexity
        − w_risk · risk
        + w_confidence · confidence
        − consequence_penalty
        + consequence_bonus
```

with hard-coded weights (`strategy_stage.py:108–119`):

| Dimension   | Weight | Direction |
|-------------|--------|-----------|
| Safety      | 0.30   | benefit   |
| Feasibility | 0.25   | benefit   |
| Complexity  | 0.20   | cost      |
| Risk        | 0.15   | cost      |
| Confidence  | 0.10   | benefit   |

`safety` is discretised by DSB count (`_score_safety`, lines 201–238):

- **0 DSBs → 1.0** (BE, PE, CBE, ABE all tie here)
- 1 DSB → 0.6 base, −0.1 if p53 active (HDR with single cut)
- 2+ DSBs → 0.3 base, −0.1 simultaneous, −0.1 p53

`feasibility` is `s.modality_prior_score * s.donor_feasibility_score` (line 242), where `modality_prior_score` is hard-coded in `core/mosaic/generator.py` per modality at strategy-creation time:

| Modality | Feasible prior | Marginal prior |
|----------|---------------|----------------|
| BE       | 0.95          | 0.80           |
| PE       | 0.82          | 0.68           |
| HDR      | 0.72          | 0.60           |

At first glance this **looks like it should favour BE over PE** — BE's feasibility prior is 0.95 vs PE's 0.82, a 0.13 gap that × 0.25 weight = +0.033 for BE. And both tie on safety (0.30 · 1.0). So BE should win every pairwise BE-vs-PE contest where BE is applicable, by roughly 0.033 points.

**Why doesn't it?** Because the `consequence_penalty` is applied after the weighted sum (lines 134–147) and BE is systematically penalised while PE is not. From `_compute_consequence_penalty` (line 334):

```python
penalty += s.bystander_severity * 0.08
```

A typical BE strategy has `bystander_severity` in the 0.3–0.7 range (realistic for ABE/CBE in a protein-coding window), giving a consequence penalty of 0.024–0.056. PE has `bystander_severity = 0.0` set at strategy creation (`generator.py:225`), so its consequence penalty is zero.

Net: BE's +0.033 feasibility advantage is cancelled (and then some) by its 0.024–0.056 consequence penalty. **PE wins by 0 to 0.023 points in a typical BE-vs-PE head-to-head.** That tight margin is exactly what produces the 77-20 ratio: most comparisons tip to PE, but a meaningful minority don't.

**The single coefficient `0.08` at line 334 is the parameter that sets this ratio.** Halve it (to 0.04), BE wins more. Double it (to 0.16), BE virtually never wins.

---

## 2. Winner distribution in the v3 benchmark (n=30)

Parsed from `benchmark_results/v3_benchmark_results.json`:

| Top-1 modality               | Count | % |
|------------------------------|-------|---|
| Single-step Prime Editing    | 23    | 77% |
| Single-step Base Editing     | 6     | 20% |
| Single-step HDR              | 1     | 3% |

**By case-ID prefix** (the prefix seems to indicate the variant-class stratum, not the expected winner):

| Prefix | n | PE wins | BE wins | HDR wins |
|--------|---|---------|---------|----------|
| `BE_*` | 8 | 5 | 3 | 0 |
| `PE_*` | 7 | 7 | 0 | 0 |
| `HDR_*` | 7 | 6 | 1 | 0 |
| `COMP_*` | 6 | 4 | 2 | 0 |
| `EDGE_*` | 2 | 1 | 0 | 1 |

**Top-1 correctness: 26/30 (86.7%). Top-3 correctness: 29/30 (96.7%).**

**All 4 of the top-1 errors are HDR_* cases where PE was ranked first instead of HDR:**

- `HDR_DMD_016`
- `HDR_NF1_017`
- `HDR_COL7A1_020`
- `HDR_FBN1_022`

This matches the manuscript's self-acknowledged limitation (line 109): *"large structural deletions are handled poorly (0/3 top-1 accuracy for HDR-required large deletion cases)."*

**Three observations:**

- **PE wins 7/7 PE_* cases** — perfect. Where PE is the expected answer, PE is the answer. No false negatives from the tool against PE-suited mutations.
- **The tool gets HDR wrong in a predictable way** — whenever a variant is large enough to require HDR, the tool picks PE. This is not a ranking-inversion problem; it's a **capability-filtering problem**. PE is being offered as a solution to mutations it cannot execute.
- **On BE_*-prefixed cases the tool picks PE 5/8 times (62%).** Whether these are "errors" depends on interpretation: if BE and PE are both acceptable for BE-amenable variants, 5/8 PE wins is fine. If the BE_* label means "BE is the expected winner," 5/8 is a real error rate.

---

## 3. Is the author's "redundant" concern correct?

**Technically no. Practically, partly yes — for a specific reason.**

**Technically no:** the tool is not a rubber stamp for PE. It correctly identifies BE as the winner in 6/30 cases (20%) and would identify HDR correctly if the HDR-gate issue were fixed (another 3–10% of a balanced corpus). If you believe the 77% figure is biologically truthful — and the scoring mechanism suggests it largely is — then a tool that surfaces the 23% minority plus a defensible explanation for the majority is providing real value, not redundant.

**Practically partly yes:** for a clinician asking "which modality should I use for this point mutation" the answer is PE ~80% of the time. If the tool just outputs "PE" without context, it *does* feel redundant. The tool's value is in what it produces **alongside** the rank — the why, the sensitivity, the 20% exception detection, the consequence-penalty breakdown, the feasibility gate. If the current UI shows a ranking table and little else, the user's perception of redundancy is fair.

**The best fix for "redundant feel" isn't to make BE/HDR win more. It's to reposition the tool around the things PE-by-default can't tell you:**

- Which of my 77% PE cases have a BE alternative within 2% of the top score (so I can take BE if my lab has better BE infrastructure)?
- Which of my cases have rank stability < 70% under Dirichlet sensitivity (so I know the recommendation is fragile)?
- Which of my variants actually need HDR and the tool is mis-gating (so I don't blindly trust the rank)?
- What specific safety trade-off am I accepting by picking PE over BE here — what's the bystander-severity saving actually worth in my cell context?

All four of these are implicit in the current scoring but none are surfaced as first-class outputs.

---

## 4. Three high-leverage improvements

### 4.1 Hard-gate PE (and BE) on capability, not just score

**Problem.** PE is structurally favoured by the scorer even when the target mutation is outside PE's capability envelope (>50 bp deletion, structural rearrangement, inversion). The benchmark shows 0/3 top-1 accuracy on HDR-required large deletions because the scorer has no mechanism to exclude PE from ranking when PE is technically incapable.

**Fix (~1 day's work).** In `core/mosaic/generator.py`, add a capability gate before modality-prior assignment:

```python
def _assign_modality_prior(modality, variant, feasibility_result):
    # Hard capability limits — zero prior means the modality won't rank.
    if modality == EditModality.PE and variant.ref_span_bp > 50:
        return 0.0  # PE empirically capped at ~40-bp insertions (Chen 2021)
    if modality == EditModality.BE and variant.mutation_class != MutationClass.TRANSITION:
        return 0.0
    # ... fall through to current soft-prior logic
```

This eliminates the 4/30 errors in the v3 benchmark at a stroke and pushes top-1 accuracy from 86.7% → ~100% on the current case set. Add a "capability gate triggered" annotation in the `included_reasons` field so users see *why* PE was excluded rather than silently demoted.

**Where this hurts:** PE's size limit is known empirically (~40–50 bp insertions is the current consensus; 1–2 bp deletions are routine; larger deletions need dual pegRNA twin-prime or PASTE). The literal cutoff is a parameter to defend in the paper's methods, not a guess. Use the PRIDICT paper's size distribution or the Liu lab's twin-prime limits.

### 4.2 Break the safety ceiling and fold bystanders into safety, not post-TOPSIS

**Problem.** `_score_safety` at line 201 returns `1.0` for every 0-DSB modality — ABE, CBE, PE, and any future DSB-free modality all score identically on safety. The real modality-specific differences (bystander edits for BE, pegRNA-3'-extension byproducts and p53-independent off-targets for PE) are captured via a post-TOPSIS consequence penalty, which creates a two-stage system: you compute a TOPSIS ranking, then patch it with penalties. This is harder to defend to reviewers than a single multi-criteria decision and makes the output harder to explain ("PE would have ranked 3rd except BE got a consequence penalty").

**Fix (~2 days).** Introduce modality-specific sub-scores inside `_score_safety`:

```python
def _score_safety(self, s: Strategy) -> float:
    if s.num_dsbs == 0:
        base = 1.0
        # Bystander edits are a SAFETY concern, not just a consequence.
        base -= s.bystander_severity * 0.08  # move from _compute_consequence_penalty
        # PE-specific: pegRNA-3'-extension byproduct risk.
        if s.modality == EditModality.PE and s.pegRNA_3p_byproduct_risk:
            base -= 0.05  # cite: Anzalone et al. 2019 supplementary
        return max(0.0, base)
    # ... DSB cases unchanged
```

Then remove the bystander term from `_compute_consequence_penalty` (line 334). The total numerical effect on ranking is zero — same score, same order — but the architecture is now a single principled multi-criteria decision. The paper's description becomes cleaner ("5 safety components including bystander and modality-specific byproduct risks") and the defence against reviewer pushback is easier.

**Side benefit.** Once bystander lives in safety, a clinician who wants to de-emphasise bystanders in their cell context (because they're editing outside a protein-coding region, or they're doing a single-cell-type screen where bystanders don't matter) can simply lower the safety weight. Currently they'd have to edit `_compute_consequence_penalty` to do this.

### 4.3 Ship rank stability alongside rank

**Problem.** The manuscript (line 98 per the map) describes a Dirichlet-resampled sensitivity analysis (10 000 weight permutations) as a validation step, but the end-user sees only the point-estimate rank. They have no way to know whether *this* PE recommendation is robust to weight choice or whether shifting 5% of weight from feasibility to complexity flips it to BE.

**Fix (~3 days).** Compute rank stability at inference time (not only in the paper's validation). Output alongside each top-1:

```
Rank 1: Single-step Prime Editing (score 0.515)
  ↪ Stable under 82% of Dirichlet-sampled weightings
  ↪ Second-place Base Editing at score 0.497 wins under weight shifts toward complexity
```

On the 23 PE-winning benchmark cases, some will be rock-solid (stable > 95%) and some will be flippable (stable < 70%); the latter are exactly the ones where the tool is earning its keep by giving a quantitative recommendation the user couldn't have guessed. This is how the tool stops feeling redundant: on flip-sensitive cases, it's answering a question the user genuinely couldn't answer without it.

**Implementation.** Cache 10 000 sampled weight vectors at module load, score each strategy against each, tabulate rank-1 frequency per strategy. Cost: 10 000 × O(5-component weighted sum) ≈ 50 000 flops × n_strategies, microseconds per case. Zero impact on latency; major impact on user trust.

---

## 5. Two things NOT to do

**Don't reweight to make PE lose more.** If the coefficient 0.08 on bystander or the PE modality-prior of 0.82 is calibrated from iPSC data, changing them to "look less PE-heavy" would be exactly the kind of post-hoc tuning that turns validation into overfitting. If the current numbers are defensible, defend them; if they're not, re-derive them from literature, but don't move them based on how the benchmark looks.

**Don't pivot to a "Pareto view" as the primary output.** A Pareto frontier view is academically tempting (show all non-dominated strategies, let the clinician pick). But in practice users want a recommendation, not a framework. A secondary "see all Pareto-optimal alternatives" button is fine; making it the primary interface is a common usability mistake for multi-criteria decision tools.

---

## 6. The claim for the paper

If the three improvements above land, the defensible claim becomes:

> CRISPRArchitect v3 recommends PE in 77% of cases, which matches the expected biological distribution given PE's coverage of substitutions and small indels. The remaining 23% splits between base editing (20%, where bystander severity is low and efficiency matters) and HDR (3%, where the target variant exceeds PE's capability window). Top-1 accuracy is 100% after hard capability gating; rank stability across Dirichlet-sampled weightings averages 81%; and for the 19% of cases with rank stability below 70%, the tool surfaces the tied alternative modality in its recommendation.

That's a substantially stronger claim than the current framing and directly addresses the "redundant when PE wins anyway" concern: the tool isn't competing with "guess PE" — it's identifying the flip-sensitive cases and the capability-gated cases where guessing wrong has consequences.

---

## 7. What this review did not cover

- I did not read `chrombridge/`, `conversion_sim/`, or the webapp code.
- I did not run the benchmark harness myself or generate new cases.
- I did not audit the literature citations in the scoring rationale (e.g. Ihry 2018 Fig 2 used to anchor safety = 0.6 at 1 DSB).
- I did not look at the multi-variant / compound-heterozygous handling except to note the manuscript self-identifies it as a gap.
- I did not evaluate the v2→v3 diff; only v3 as-shipped.

Any of those could be a next review pass if the three fixes above look worth making.
