"""
CRISPRArchitect v2 Pipeline Orchestrator
==========================================

This is the single entry point for the v2 pipeline. It takes raw
genomic variant inputs and produces a ranked list of editing strategies
with consequence-aware scoring.

Pipeline stages
----------------
1. Fetch transcript (from Ensembl) for the gene
2. Build TranscriptMapper for coordinate translation
3. Normalize each variant (validate ref, map to transcript, annotate)
4. Fetch local genomic sequences around each variant
5. Run feasibility checks (BE, PE, HDR) for each variant
6. Generate strategies from feasibility bundles
7. Score strategies with consequence-aware adjustments
8. Rank and return

Usage
------
    from core.pipeline.strategy_stage import StrategyPipeline
    from core.models import GenomicVariantInput

    pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")
    result = pipeline.run([
        GenomicVariantInput(
            chromosome="17", position=31232193,
            ref_allele="C", alt_allele="T",
            gene_symbol="NF1", name="c.910C>T",
        ),
    ])

    for s in result.strategies:
        print(s.strategy_name, s.overall_score, s.confidence)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from core.models import (
    EditModality,
    EvidenceTier,
    FeasibilityBundle,
    FeasibilityLabel,
    GenomicVariantInput,
    NormalizedVariant,
    PipelineResult,
    RiskLevel,
    ScoredStrategy,
    Strategy,
    TranscriptInfo,
)

logger = logging.getLogger(__name__)


# ─── Safety-scoring constants ────────────────────────────────────────────
# BYSTANDER_SAFETY_COEF: coefficient by which ``bystander_severity``
# penalises the safety score. As of Fix #2 (2026-04-20), this coefficient
# lives in the safety dimension; prior to Fix #2 the same coefficient was
# applied inside ``_compute_consequence_penalty`` (the consequence TOPSIS
# dimension), which created a safety ceiling at 1.0 for every 0-DSB
# modality and made the bystander contribution difficult to explain in
# the paper. See REVIEW_NOTES_2026-04-18.md §4.2 and
# FIX_NOTES_2026-04-20.md for the refactor rationale.
BYSTANDER_SAFETY_COEF = 0.08


# ─── Rank-stability thresholds (Fix #3, 2026-04-20) ──────────────────────
# Pre-committed per REVIEW_NOTES_2026-04-18.md §4.3: a TOPSIS recommendation
# whose rank is sensitive to weight perturbation is exactly the case where
# the tool earns its keep (vs a user who would have picked the obvious
# answer by inspection). These thresholds define when to print the
# runner-up and the dimension-level tradeoff.
RANK_STABILITY_ROBUST = 0.80       # >= : unconditional recommendation
RANK_STABILITY_STABLE = 0.70       # >= : single recommendation, note uncertainty
# < RANK_STABILITY_STABLE : "flip-sensitive" — surface alternatives.

# Dimension-delta threshold for surfacing in flip-sensitive explanations.
# A delta smaller than this between top and runner-up on a given dimension
# is considered a tie (not a meaningful preference signal).
DIMENSION_DELTA_MATERIAL = 0.05


# ─── Compound-het completeness penalty (Fix #4, 2026-04-20) ──────────────
# COMPLETENESS_PENALTY_COEF: coefficient applied to the safety score for a
# strategy that addresses fewer than all pathogenic variants in a compound-
# heterozygous case. Per COL7A1_AUDIT_2026-04-19.md §7 option C: a strategy
# that resolves only 1 of N pathogenic variants is biologically incomplete
# (the patient still has one diseased allele) and should rank below hybrid
# or sequential strategies that address all N. The penalty enters safety
# as ``base -= (1 - completeness_ratio) * COMPLETENESS_PENALTY_COEF``.
#
# The coefficient 0.30 is chosen to meaningfully separate incomplete
# single-step strategies from complete hybrids without entirely vetoing
# them (that would be a hard rejection, the wrong response when the user
# may have a reason — e.g., staged therapy — to address one variant at
# a time).
COMPLETENESS_PENALTY_COEF = 0.30


# ═══════════════════════════════════════════════════════════════════════════
# Scoring Engine
# ═══════════════════════════════════════════════════════════════════════════

class StrategyScorer:
    """Multi-objective scoring engine for v2/v3 strategies.

    Scoring function:
        Score = w1*Safety + w2*Feasibility - w3*Complexity - w4*Risk + w5*Confidence

    Default weights — rationale for iPSC context
    -----------------------------------------------
    These weights are NOT arbitrary round numbers. They encode a deliberate
    priority ordering for iPSC editing, justified as follows:

      Safety (0.30):  Highest weight because DSBs in iPSCs trigger p53-mediated
        selection that can enrich TP53-mutant clones (Ihry et al., Nat Med, 2018;
        Haapaniemi et al., Nat Med, 2018). This is a patient-safety concern for
        therapeutic applications. Safety is the most consequential axis.

      Feasibility (0.25): Second-highest because a strategy that cannot be
        executed (no PAM, target outside editing window) is useless regardless
        of other scores. PAM-window constraints are the binding bottleneck
        (demonstrated by 0/30 BE recommendations with SpCas9 alone in v2).

      Complexity (0.20): Third because iPSC editing is resource-intensive
        (each round requires electroporation, clonal expansion, genotyping,
        karyotyping). Multiple rounds multiply time, cost, and passage-related
        karyotype risk (Baker et al., Nat Biotechnol, 2007).

      Risk (0.15): Fourth because rearrangement and bystander risks are
        conditional on DSB-based strategies, which are already penalized by
        the safety dimension. Risk captures residual concerns.

      Confidence (0.10): Lowest weight because all Tier A modalities (ABE, CBE,
        PE with SpCas9) are well-established. Confidence primarily distinguishes
        Tier B combinations (e.g., ABE8e-enFnCas9) from Tier A.

    IMPORTANT: These defaults are explored via TOPSIS sensitivity analysis
    (10,000 Dirichlet-sampled weight permutations). The rank stability
    metric indicates how robust the recommendation is to weight variation.
    Users should consult rank stability rather than relying solely on the
    point estimate from default weights.
    """

    def __init__(
        self,
        w_safety: float = 0.30,
        w_feasibility: float = 0.25,
        w_complexity: float = 0.20,
        w_risk: float = 0.15,
        w_confidence: float = 0.10,
    ):
        total = w_safety + w_feasibility + w_complexity + w_risk + w_confidence
        self.w_safety = w_safety / total
        self.w_feasibility = w_feasibility / total
        self.w_complexity = w_complexity / total
        self.w_risk = w_risk / total
        self.w_confidence = w_confidence / total

    def score_strategy(
        self,
        strategy: Strategy,
        bundles: List[FeasibilityBundle],
    ) -> ScoredStrategy:
        """Score a single strategy across all dimensions."""
        safety = self._score_safety(strategy)

        # Fix #4 (2026-04-20): compound-het completeness penalty.
        # A strategy addressing fewer than all pathogenic variants in a
        # multi-variant case is biologically incomplete; it should rank
        # below hybrid / sequential strategies that address all variants.
        total_variants = len(bundles)
        if total_variants > 1:
            completeness = strategy.completeness_ratio(total_variants)
            if completeness < 1.0:
                safety = max(
                    0.0,
                    safety - (1.0 - completeness) * COMPLETENESS_PENALTY_COEF,
                )

        feasibility = self._score_feasibility(strategy)
        complexity = self._score_complexity(strategy)
        risk = self._score_risk(strategy, bundles)
        confidence = self._score_confidence(strategy)

        # Consequence adjustments
        consequence_penalty = self._compute_consequence_penalty(
            strategy, bundles
        )
        consequence_bonus = self._compute_consequence_bonus(strategy, bundles)

        overall = (
            self.w_safety * safety
            + self.w_feasibility * feasibility
            - self.w_complexity * complexity
            - self.w_risk * risk
            + self.w_confidence * confidence
            - consequence_penalty
            + consequence_bonus
        )
        overall = max(0.0, min(1.0, overall))

        notes: List[str] = []
        if consequence_penalty > 0:
            notes.append(
                f"Consequence penalty: -{consequence_penalty:.3f}"
            )
        if consequence_bonus > 0:
            notes.append(
                f"Consequence bonus: +{consequence_bonus:.3f}"
            )

        return ScoredStrategy(
            strategy=strategy,
            safety_score=safety,
            feasibility_score=feasibility,
            complexity_score=complexity,
            risk_score=risk,
            confidence_score=confidence,
            consequence_penalty=consequence_penalty,
            consequence_bonus=consequence_bonus,
            overall_score=overall,
            annotation_notes=notes,
        )

    def rank(
        self,
        strategies: List[Strategy],
        bundles: List[FeasibilityBundle],
    ) -> List[ScoredStrategy]:
        """Score and rank all strategies.

        Rejected strategies are excluded from ranking but preserved
        in the PipelineResult.rejected_strategies list.
        """
        scored = []
        for s in strategies:
            if not s.is_rejected:
                scored.append(self.score_strategy(s, bundles))

        # Sort by overall_score descending, then safety descending
        scored.sort(
            key=lambda x: (x.overall_score, x.safety_score),
            reverse=True,
        )

        for i, s in enumerate(scored, 1):
            s.rank = i

        return scored

    # ─── Dimension scoring ───

    def _score_safety(self, s: Strategy) -> float:
        """0-1 where 1 = safest (no DSBs, no rearrangement risk, no bystanders).

        Safety scoring rationale
        -------------------------
        Base score by DSB count:

        0 DSBs = 1.0 : Base editing and prime editing introduce no DSBs.
            No p53 activation, no translocation risk, no large deletions.

        1 DSB = 0.6 : Single DSB triggers ATM/p53 signaling. In iPSCs,
            ~45% cell death (Ihry et al., Nat Med, 2018 Fig 2).
            p53 penalty: -0.1 reflects the additional risk of TP53-mutant
            clone enrichment specific to p53-active cell types.
            Value rationale: 0.6 is the midpoint of viable HDR outcomes
            (editing succeeds in ~8% of surviving cells, vs. 30-70% for BE).

        2+ DSBs = 0.3 : Two simultaneous DSBs create translocation risk
            (~1-10% for same-chromosome loci; Kosicki et al., 2018).
            Simultaneous penalty: -0.1 (translocations only form between
            co-temporal DSBs; sequential rounds eliminate this risk).
            p53 penalty: -0.1 (dual DSBs cause stronger p53 activation).
            Value rationale: 0.3 reflects that dual-DSB strategies produce
            viable, correctly-edited clones in <5% of attempts in iPSCs.

        Bystander contribution (added 2026-04-20, Fix #2):

        Bystander off-target edits reduce the safety score by
        ``bystander_severity * BYSTANDER_SAFETY_COEF``. Before Fix #2, this
        term lived as a post-hoc penalty in the separate ``consequence``
        TOPSIS dimension (weight 0.08). It now lives in ``safety`` (weight
        0.30) alongside DSB-derived risk, making safety a single principled
        multi-criteria descriptor instead of a ceiling.

        Side-effect: two DSB-free modalities with different bystander
        profiles (e.g., ABE at a C-rich locus vs PE at the same locus)
        now differ on the safety axis, which previously could not
        distinguish them. This was the entire point of the refactor
        per REVIEW_NOTES_2026-04-18.md §4.2.
        """
        if s.num_dsbs == 0:
            base = 1.0
        elif s.num_dsbs == 1 and not s.simultaneous_dsbs:
            base = 0.6
            if s.p53_active:
                base -= 0.1
        elif s.num_dsbs >= 2:
            base = 0.3
            if s.simultaneous_dsbs:
                base -= 0.1
            if s.p53_active:
                base -= 0.1
        else:
            base = 0.5

        base -= s.bystander_severity * BYSTANDER_SAFETY_COEF
        return max(0.0, base)

    def _score_feasibility(self, s: Strategy) -> float:
        """0-1 based on modality prior score and donor quality."""
        return s.modality_prior_score * s.donor_feasibility_score

    def _score_complexity(self, s: Strategy) -> float:
        """0-1 where 1 = most complex (worst). Lower is better.

        Complexity sub-weight rationale (iPSC context)
        ------------------------------------------------
        Editing rounds (weight 0.35, penalty 0.3/round):
            Each round requires electroporation -> 48h recovery -> clonal
            expansion (2-3 weeks) -> genotyping -> karyotyping. The 0.3
            penalty per additional round reflects ~30% of maximum complexity
            per round, so a 2-round strategy scores 0.3 and a 4-round
            strategy saturates at 1.0. This is the dominant complexity
            factor because rounds are the primary time bottleneck.

        Donors (weight 0.25, penalty 0.15/donor):
            Each donor template requires separate design, synthesis, and QC.
            cssDNA production takes ~30 hours (Iyer et al., 2022). Multiple
            donors also complicate screening (must verify each independently).

        Guides (weight 0.20, penalty 0.1/guide):
            Additional guides require separate synthesis and validation but
            are less costly than additional donors or rounds.

        Screening (weight 0.20, normalized to 100 clones):
            iPSC clonal screening requires single-cell sorting, expansion,
            genotyping, and karyotyping. 100 clones is the practical upper
            limit for most labs; beyond this, the experiment is impractical.
        """
        rounds_penalty = min(1.0, (s.num_rounds - 1) * 0.3)
        donor_penalty = min(1.0, s.num_donors * 0.15)
        guide_penalty = min(1.0, (s.num_distinct_guides - 1) * 0.1)
        screening_penalty = min(1.0, s.screening_clones / 100.0)
        return (
            0.35 * rounds_penalty
            + 0.25 * donor_penalty
            + 0.20 * guide_penalty
            + 0.20 * screening_penalty
        )

    def _score_risk(self, s: Strategy, bundles: List[FeasibilityBundle]) -> float:
        """0-1 where 1 = highest risk (worst). Lower is better.

        Risk now captures ONLY rearrangement/translocation risk (structural
        genotoxicity). Bystander edit risk has been moved to the dedicated
        consequence dimension (6th TOPSIS dimension) to avoid double-counting.

        Previous versions counted bystander severity in BOTH the risk
        dimension AND the consequence penalty, creating a triple-counting
        effect that artificially penalized base editing relative to prime
        editing. This was the root cause of the degenerate "always PE"
        pattern in the v2 benchmark (29/30 PE top-ranked).
        """
        risk = 0.0

        # Rearrangement risk (structural: deletions, inversions, translocations)
        risk_map = {
            RiskLevel.LOW: 0.0,
            RiskLevel.MODERATE: 0.3,
            RiskLevel.HIGH: 0.6,
            RiskLevel.VERY_HIGH: 0.9,
        }
        risk += risk_map.get(s.rearrangement_risk, 0.0)

        # NOTE: bystander_severity is NO LONGER included here.
        # It is captured in the consequence dimension (6th TOPSIS dim).

        return min(1.0, risk)

    def _score_confidence(self, s: Strategy) -> float:
        """0-1 based on evidence tier."""
        tier_map = {
            EvidenceTier.A: 1.0,
            EvidenceTier.B: 0.7,
            EvidenceTier.C: 0.4,
        }
        return tier_map.get(s.evidence_tier, 0.5)

    def _compute_consequence_penalty(
        self, s: Strategy, bundles: List[FeasibilityBundle]
    ) -> float:
        """Penalty for splice-proximal edits.

        As of 2026-04-20 (Fix #2 per REVIEW_NOTES_2026-04-18.md §4.2),
        bystander severity has been moved from this penalty into the
        ``safety`` dimension (see ``_score_safety``). This function now
        captures only splice-proximity consequences — it remains the
        source of the ``consequence`` TOPSIS dimension, which is still
        meaningful because splice-site effects are an editing-outcome
        consequence distinct from the safety / rearrangement axis.
        """
        penalty = 0.0

        for bundle in bundles:
            v = bundle.variant
            if v.coding.splice_proximity == "near_exon_start":
                penalty += 0.05
            elif v.coding.splice_proximity == "near_exon_end":
                penalty += 0.05

        return min(0.3, penalty)  # cap total penalty

    def _compute_consequence_bonus(
        self, s: Strategy, bundles: List[FeasibilityBundle]
    ) -> float:
        """Bonus for clean designs.

        Note: The clean design bonus is deliberately small (0.03) and now
        serves only as a tie-breaker within the legacy weighted-sum scorer.
        In the TOPSIS pipeline, clean designs are captured by the consequence
        dimension score being 1.0 (no bystander penalties), which is the
        mathematically principled representation.
        """
        bonus = 0.0

        # No DSBs + no bystanders = inherently clean
        if s.num_dsbs == 0 and s.bystander_severity == 0.0:
            bonus += 0.03

        return min(0.1, bonus)


# ═══════════════════════════════════════════════════════════════════════════
# TOPSIS Scorer (v3)
# ═══════════════════════════════════════════════════════════════════════════

import math
import numpy as np


@dataclass
class SensitivityResult:
    """Result of Monte Carlo sensitivity analysis on strategy rankings.

    Attributes
    ----------
    strategy_name : str
        Name of the strategy.
    rank_stability : float
        Fraction of weight permutations where this strategy is top-ranked.
        Range [0, 1]. Higher = more robust ranking.
    mean_rank : float
        Mean rank across all permutations (1 = best).
    rank_distribution : Dict[int, float]
        Distribution of ranks: {rank: fraction_of_time}.
    """
    strategy_name: str
    rank_stability: float = 0.0
    mean_rank: float = 1.0
    rank_distribution: Dict[int, float] = field(default_factory=dict)


class TOPSISScorer:
    """TOPSIS multi-criteria decision scorer with sensitivity analysis.

    TOPSIS (Technique for Order Preference by Similarity to Ideal Solution)
    ranks alternatives based on their geometric distance to the ideal and
    anti-ideal solutions in normalized decision space.

    This is more principled than a simple weighted sum because:
    1. It penalizes strategies that are terrible on ANY single axis
       (even if great on others)
    2. It handles dimensions that shouldn't trade off linearly
    3. It's a well-established MCDM method (Hwang & Yoon, 1981)

    Sensitivity analysis (Monte Carlo weight perturbation) reports how
    stable each ranking is across 10,000 random weight vectors. This is
    the genuinely novel contribution for CRISPR tool scoring.

    Parameters
    ----------
    w_safety : float
        Default weight for safety dimension.
    w_feasibility : float
        Default weight for feasibility dimension.
    w_complexity : float
        Default weight for complexity dimension.
    w_risk : float
        Default weight for risk dimension.
    w_confidence : float
        Default weight for confidence dimension.
    n_sensitivity_runs : int
        Number of Monte Carlo permutations for sensitivity analysis.

    References
    ----------
    Hwang & Yoon, Multiple Attribute Decision Making, Springer, 1981.
    """

    # 6 dimensions: safety, feasibility, complexity, risk, confidence, consequence
    # Benefit dimensions (higher = better): safety, feasibility, confidence, consequence
    # Cost dimensions (lower = better): complexity, risk
    _BENEFIT_DIMS = {0, 1, 4, 5}  # safety, feasibility, confidence, consequence
    _COST_DIMS = {2, 3}            # complexity, risk
    _DIM_NAMES = [
        "safety", "feasibility", "complexity", "risk", "confidence", "consequence"
    ]

    def __init__(
        self,
        w_safety: float = 0.30,
        w_feasibility: float = 0.25,
        w_complexity: float = 0.20,
        w_risk: float = 0.15,
        w_confidence: float = 0.10,
        w_consequence: float = 0.08,
        n_sensitivity_runs: int = 10000,
    ):
        total = (w_safety + w_feasibility + w_complexity + w_risk
                 + w_confidence + w_consequence)
        self.weights = [
            w_safety / total,
            w_feasibility / total,
            w_complexity / total,
            w_risk / total,
            w_confidence / total,
            w_consequence / total,
        ]
        self.n_sensitivity_runs = n_sensitivity_runs
        # Legacy scorer for computing per-dimension scores
        self._legacy_scorer = StrategyScorer(
            w_safety=w_safety,
            w_feasibility=w_feasibility,
            w_complexity=w_complexity,
            w_risk=w_risk,
            w_confidence=w_confidence,
        )

    def rank(
        self,
        strategies: List[Strategy],
        bundles: List[FeasibilityBundle],
        run_sensitivity: bool = True,
    ) -> List[ScoredStrategy]:
        """Score and rank strategies using TOPSIS.

        Parameters
        ----------
        strategies : list of Strategy
            Candidate strategies (rejected ones are excluded).
        bundles : list of FeasibilityBundle
            Feasibility results for consequence adjustments.
        run_sensitivity : bool
            Whether to run Monte Carlo sensitivity analysis.

        Returns
        -------
        List[ScoredStrategy]
            Ranked strategies with TOPSIS scores, sensitivity in metadata.
        """
        active = [s for s in strategies if not s.is_rejected]
        if not active:
            return []

        # Step 1: Compute raw dimension scores for each strategy
        # 6 dimensions: safety, feasibility, complexity, risk, confidence, consequence
        # The consequence dimension replaces the post-hoc additive penalty.
        # It is a benefit dimension: higher = fewer bystander/splice issues.
        dim_matrix = []
        scored_list = []
        for s in active:
            scored = self._legacy_scorer.score_strategy(s, bundles)
            scored_list.append(scored)
            # Consequence score: 1.0 = no issues, decreases with penalties
            consequence_score = max(0.0, min(1.0,
                1.0 - scored.consequence_penalty + scored.consequence_bonus
            ))
            dim_matrix.append([
                scored.safety_score,
                scored.feasibility_score,
                scored.complexity_score,
                scored.risk_score,
                scored.confidence_score,
                consequence_score,
            ])

        if len(active) == 1:
            scored_list[0].rank = 1
            return scored_list

        # Step 2: TOPSIS ranking (consequences now integrated as 6th dimension)
        topsis_scores = self._topsis(dim_matrix, self.weights)

        # Step 3: Assign TOPSIS scores directly (no post-hoc adjustment needed)
        for i, scored in enumerate(scored_list):
            scored.overall_score = round(topsis_scores[i], 4)

        # Step 4: Sort by TOPSIS score
        scored_list.sort(
            key=lambda x: (x.overall_score, x.safety_score),
            reverse=True,
        )
        for i, s in enumerate(scored_list, 1):
            s.rank = i

        # Step 5: Pareto front analysis
        # Identify non-dominated strategies across all 6 dimensions.
        # A strategy is Pareto-dominated if another strategy is at least
        # as good on ALL dimensions and strictly better on at least one.
        pareto_flags = self._pareto_front(dim_matrix)
        for i, scored in enumerate(scored_list):
            # Find which index in dim_matrix corresponds to this scored strategy
            orig_idx = active.index(scored.strategy)
            is_pareto = pareto_flags[orig_idx]
            if is_pareto:
                scored.annotation_notes.append("Pareto non-dominated")
            else:
                scored.annotation_notes.append("Pareto dominated")

        # Step 6: Sensitivity analysis
        if run_sensitivity and len(active) >= 2:
            sensitivity = self._sensitivity_analysis(
                dim_matrix, scored_list
            )
            for scored in scored_list:
                sr = sensitivity.get(scored.strategy_name)
                if sr:
                    scored.rank_stability = sr.rank_stability
                    scored.annotation_notes.append(
                        f"Rank stability: {sr.rank_stability:.1%} "
                        f"(top-ranked in {sr.rank_stability:.1%} of "
                        f"{self.n_sensitivity_runs} weight permutations)"
                    )
                    scored.annotation_notes.append(
                        f"Mean rank: {sr.mean_rank:.1f}"
                    )

        return scored_list

    def _topsis(
        self,
        matrix: List[List[float]],
        weights: List[float],
    ) -> List[float]:
        """Core TOPSIS algorithm.

        Steps:
        1. Normalize the decision matrix (vector normalization)
        2. Apply weights
        3. Determine ideal (A+) and anti-ideal (A-) solutions
        4. Compute Euclidean distances to A+ and A-
        5. Calculate relative closeness: C = D- / (D+ + D-)

        Returns list of TOPSIS scores (0-1, higher = better).
        """
        n_alts = len(matrix)
        n_dims = len(matrix[0])

        # Step 1: Vector normalization
        # Each value normalized by sqrt(sum of squares in that column)
        norm = [[0.0] * n_dims for _ in range(n_alts)]
        for j in range(n_dims):
            col_sum_sq = sum(matrix[i][j] ** 2 for i in range(n_alts))
            col_norm = math.sqrt(col_sum_sq) if col_sum_sq > 0 else 1.0
            for i in range(n_alts):
                norm[i][j] = matrix[i][j] / col_norm

        # Step 2: Weighted normalized matrix
        weighted = [[0.0] * n_dims for _ in range(n_alts)]
        for i in range(n_alts):
            for j in range(n_dims):
                weighted[i][j] = norm[i][j] * weights[j]

        # Step 3: Ideal and anti-ideal solutions
        ideal = [0.0] * n_dims
        anti_ideal = [0.0] * n_dims
        for j in range(n_dims):
            col_vals = [weighted[i][j] for i in range(n_alts)]
            if j in self._BENEFIT_DIMS:
                # Benefit: higher is better
                ideal[j] = max(col_vals)
                anti_ideal[j] = min(col_vals)
            else:
                # Cost: lower is better
                ideal[j] = min(col_vals)
                anti_ideal[j] = max(col_vals)

        # Step 4: Euclidean distances
        dist_to_ideal = []
        dist_to_anti = []
        for i in range(n_alts):
            d_plus = math.sqrt(sum(
                (weighted[i][j] - ideal[j]) ** 2 for j in range(n_dims)
            ))
            d_minus = math.sqrt(sum(
                (weighted[i][j] - anti_ideal[j]) ** 2 for j in range(n_dims)
            ))
            dist_to_ideal.append(d_plus)
            dist_to_anti.append(d_minus)

        # Step 5: Relative closeness
        scores = []
        for i in range(n_alts):
            denom = dist_to_ideal[i] + dist_to_anti[i]
            if denom == 0:
                scores.append(0.5)
            else:
                scores.append(dist_to_anti[i] / denom)

        return scores

    def _pareto_front(
        self,
        dim_matrix: List[List[float]],
    ) -> List[bool]:
        """Identify Pareto-non-dominated strategies.

        A strategy A dominates strategy B if A is at least as good as B on
        ALL dimensions and strictly better on at least one. Pareto-non-
        dominated strategies form the Pareto front — the set of strategies
        where no other strategy is universally better.

        This analysis is weight-independent: it does not depend on the
        relative importance of dimensions. A strategy on the Pareto front
        is defensible under SOME weighting scheme.

        For benefit dimensions (safety, feasibility, confidence, consequence):
        higher is better. For cost dimensions (complexity, risk): lower is better.

        Parameters
        ----------
        dim_matrix : List of List[float]
            Decision matrix, shape (n_alternatives, n_dimensions).

        Returns
        -------
        List[bool]
            True for each strategy that is Pareto-non-dominated.
        """
        n = len(dim_matrix)
        if n <= 1:
            return [True] * n

        # Convert cost dimensions to benefit (negate them) for uniform comparison
        converted = []
        for row in dim_matrix:
            new_row = []
            for j, val in enumerate(row):
                if j in self._COST_DIMS:
                    new_row.append(-val)  # lower cost -> higher converted value
                else:
                    new_row.append(val)
            converted.append(new_row)

        is_pareto = [True] * n
        for i in range(n):
            if not is_pareto[i]:
                continue
            for k in range(n):
                if i == k:
                    continue
                # Check if k dominates i:
                # k >= i on all dims AND k > i on at least one
                at_least_as_good = all(
                    converted[k][j] >= converted[i][j]
                    for j in range(len(converted[0]))
                )
                strictly_better = any(
                    converted[k][j] > converted[i][j]
                    for j in range(len(converted[0]))
                )
                if at_least_as_good and strictly_better:
                    is_pareto[i] = False
                    break

        return is_pareto

    def _sensitivity_analysis(
        self,
        dim_matrix: List[List[float]],
        scored_list: List[ScoredStrategy],
    ) -> Dict[str, SensitivityResult]:
        """Monte Carlo sensitivity analysis: perturb weights and re-rank.

        Generates n_sensitivity_runs random weight vectors from a Dirichlet
        distribution centered on the default weights, runs TOPSIS with each,
        and counts how often each strategy is top-ranked.

        Uses Dirichlet(alpha) where alpha_j = max(w_j * concentration, min_alpha).
        The min_alpha floor prevents small-weight dimensions (e.g., confidence
        at w=0.10) from having near-Uniform marginals, which would create
        asymmetric perturbation where small weights dominate some draws.

        Concentration=20 was chosen so that 95% of sampled weights for the
        largest dimension (safety, w=0.30) fall within [0.15, 0.50], which
        represents a reasonable range of "safety-prioritized" to "balanced"
        weighting schemes. This provides meaningful sensitivity exploration
        without pathological weight vectors.
        """
        n_alts = len(dim_matrix)
        if n_alts < 2:
            return {}

        strategy_names = [s.strategy_name for s in scored_list]
        rank_counts = {name: {} for name in strategy_names}  # type: Dict[str, Dict[int, int]]
        top1_counts = {name: 0 for name in strategy_names}

        # Dirichlet concentration: higher = closer to default weights
        # min_alpha=2.0 ensures no dimension has a near-Uniform marginal
        concentration = 20.0
        min_alpha = 2.0
        alphas = [max(w * concentration, min_alpha) for w in self.weights]

        # Use NumPy's PCG64 generator for better quality randomness and
        # vectorized Dirichlet sampling. Fixed seed for reproducibility.
        np_rng = np.random.default_rng(42)

        # Pre-generate all weight vectors at once (much faster than per-iteration)
        all_weights = np_rng.dirichlet(alphas, size=self.n_sensitivity_runs)

        for run_idx in range(self.n_sensitivity_runs):
            perturbed_weights = all_weights[run_idx].tolist()

            # Run TOPSIS with perturbed weights
            topsis_scores = self._topsis(dim_matrix, perturbed_weights)

            # Rank
            indexed = list(enumerate(topsis_scores))
            indexed.sort(key=lambda x: x[1], reverse=True)

            for rank, (idx, _score) in enumerate(indexed, 1):
                name = strategy_names[idx]
                rank_counts[name][rank] = rank_counts[name].get(rank, 0) + 1
                if rank == 1:
                    top1_counts[name] += 1

        # Build results with standard errors on rank stability
        results = {}
        n = self.n_sensitivity_runs
        for name in strategy_names:
            counts = rank_counts[name]
            mean_rank = sum(r * c for r, c in counts.items()) / n
            rank_dist = {r: c / n for r, c in sorted(counts.items())}
            p_top1 = top1_counts[name] / n
            # SE of a binomial proportion
            se_top1 = float(np.sqrt(p_top1 * (1.0 - p_top1) / n))
            results[name] = SensitivityResult(
                strategy_name=name,
                rank_stability=p_top1,
                mean_rank=round(mean_rank, 2),
                rank_distribution=rank_dist,
            )

        return results


# ═══════════════════════════════════════════════════════════════════════════
# Alternative MCDM Methods (for method-robustness comparison)
# ═══════════════════════════════════════════════════════════════════════════


def vikor_rank(
    matrix: List[List[float]],
    weights: List[float],
    benefit_dims: set,
    cost_dims: set,
    v: float = 0.5,
) -> List[float]:
    """VIKOR (VlseKriterijumska Optimizacija I Kompromisno Resenje) ranking.

    VIKOR identifies the compromise solution closest to the ideal by
    minimizing the group utility (S) and individual regret (R).

    Unlike TOPSIS (which uses Euclidean distance), VIKOR uses L1 (Manhattan)
    distance for group utility and Linfinity (Chebyshev) distance for
    individual regret. The parameter v controls the trade-off:
    - v=1.0: pure group utility (best average across all dimensions)
    - v=0.0: pure individual regret (minimize worst-case dimension)
    - v=0.5: balanced compromise (default, recommended by Opricovic & Tzeng, 2004)

    Parameters
    ----------
    matrix : List of List[float]
        Decision matrix (n_alternatives x n_dimensions).
    weights : List[float]
        Dimension weights (must sum to 1).
    benefit_dims : set
        Indices of benefit dimensions (higher = better).
    cost_dims : set
        Indices of cost dimensions (lower = better).
    v : float
        Weight of group utility vs individual regret (default 0.5).

    Returns
    -------
    List[float]
        VIKOR Q scores (lower = better, range [0, 1]).

    References
    ----------
    Opricovic S, Tzeng GH. Compromise solution by MCDM methods: A
    comparative analysis of VIKOR and TOPSIS. Eur J Oper Res. 2004;156:445-455.
    """
    n_alts = len(matrix)
    n_dims = len(matrix[0])
    if n_alts < 2:
        return [0.5] * n_alts

    # Step 1: Determine best (f*) and worst (f-) values per dimension
    f_best = [0.0] * n_dims
    f_worst = [0.0] * n_dims
    for j in range(n_dims):
        col = [matrix[i][j] for i in range(n_alts)]
        if j in benefit_dims:
            f_best[j] = max(col)
            f_worst[j] = min(col)
        else:
            f_best[j] = min(col)
            f_worst[j] = max(col)

    # Step 2: Compute S (group utility) and R (individual regret)
    S = []  # Manhattan weighted distance
    R = []  # Chebyshev weighted distance
    for i in range(n_alts):
        s_i = 0.0
        r_i = 0.0
        for j in range(n_dims):
            denom = f_best[j] - f_worst[j]
            if abs(denom) < 1e-12:
                continue
            ratio = weights[j] * abs(f_best[j] - matrix[i][j]) / denom
            s_i += ratio
            r_i = max(r_i, ratio)
        S.append(s_i)
        R.append(r_i)

    # Step 3: Compute Q (compromise score)
    S_best, S_worst = min(S), max(S)
    R_best, R_worst = min(R), max(R)

    Q = []
    for i in range(n_alts):
        s_norm = (S[i] - S_best) / (S_worst - S_best) if S_worst != S_best else 0.0
        r_norm = (R[i] - R_best) / (R_worst - R_best) if R_worst != R_best else 0.0
        Q.append(v * s_norm + (1.0 - v) * r_norm)

    return Q


def wpm_rank(
    matrix: List[List[float]],
    weights: List[float],
    benefit_dims: set,
    cost_dims: set,
) -> List[float]:
    """Weighted Product Model (WPM) ranking.

    WPM is a multiplicative MCDM method where the score of each alternative
    is the weighted geometric mean of its dimension values. Unlike additive
    methods (weighted sum) and distance methods (TOPSIS), WPM is
    non-compensatory: a zero on any dimension zeros the total score.

    This matches the biological reality that a strategy with zero safety
    (lethal DSBs) should never be recommended regardless of efficiency.

    Parameters
    ----------
    matrix : List of List[float]
        Decision matrix (n_alternatives x n_dimensions).
    weights : List[float]
        Dimension weights (must sum to 1).
    benefit_dims : set
        Indices of benefit dimensions.
    cost_dims : set
        Indices of cost dimensions.

    Returns
    -------
    List[float]
        WPM scores (higher = better). Unnormalized.

    References
    ----------
    Bridgman PW. Dimensional Analysis. Yale University Press, 1922.
    Triantaphyllou E. Multi-Criteria Decision Making Methods. Springer, 2000.
    """
    n_alts = len(matrix)
    n_dims = len(matrix[0])

    scores = []
    for i in range(n_alts):
        product = 1.0
        for j in range(n_dims):
            val = max(matrix[i][j], 1e-10)  # avoid log(0)
            if j in benefit_dims:
                product *= val ** weights[j]
            else:
                # Cost dimension: invert so lower is better
                product *= (1.0 / val) ** weights[j]
        scores.append(product)

    return scores


def cross_method_comparison(
    matrix: List[List[float]],
    weights: List[float],
    strategy_names: List[str],
    benefit_dims: set = None,
    cost_dims: set = None,
) -> Dict[str, Dict[str, Any]]:
    """Run TOPSIS, VIKOR, and WPM on the same decision matrix.

    Returns a comparison dictionary showing the rank assigned by each
    method for each strategy, plus rank concordance metrics.

    Parameters
    ----------
    matrix : decision matrix (n_alts x n_dims)
    weights : dimension weights
    strategy_names : names for each alternative
    benefit_dims : set of benefit dimension indices
    cost_dims : set of cost dimension indices

    Returns
    -------
    Dict with keys:
        "per_strategy": {name: {"topsis_rank": int, "vikor_rank": int, "wpm_rank": int}}
        "concordance": {"topsis_vikor": float, "topsis_wpm": float, "vikor_wpm": float}
        "all_agree_on_top1": bool
    """
    if benefit_dims is None:
        benefit_dims = {0, 1, 4, 5}  # safety, feasibility, confidence, consequence
    if cost_dims is None:
        cost_dims = {2, 3}  # complexity, risk

    n = len(matrix)
    if n < 2:
        return {"per_strategy": {strategy_names[0]: {"topsis_rank": 1, "vikor_rank": 1, "wpm_rank": 1}},
                "concordance": {"topsis_vikor": 1.0, "topsis_wpm": 1.0, "vikor_wpm": 1.0},
                "all_agree_on_top1": True}

    # Run all three methods
    topsis_scorer = TOPSISScorer()
    topsis_scores = topsis_scorer._topsis(matrix, weights)
    vikor_scores = vikor_rank(matrix, weights, benefit_dims, cost_dims)
    wpm_scores = wpm_rank(matrix, weights, benefit_dims, cost_dims)

    # Rank each (TOPSIS: higher=better, VIKOR: lower=better, WPM: higher=better)
    def rank_scores(scores, higher_is_better=True):
        indexed = list(enumerate(scores))
        indexed.sort(key=lambda x: x[1], reverse=higher_is_better)
        ranks = [0] * len(scores)
        for rank, (idx, _) in enumerate(indexed, 1):
            ranks[idx] = rank
        return ranks

    t_ranks = rank_scores(topsis_scores, higher_is_better=True)
    v_ranks = rank_scores(vikor_scores, higher_is_better=False)  # VIKOR: lower Q = better
    w_ranks = rank_scores(wpm_scores, higher_is_better=True)

    # Per-strategy results
    per_strategy = {}
    for i, name in enumerate(strategy_names):
        per_strategy[name] = {
            "topsis_rank": t_ranks[i],
            "vikor_rank": v_ranks[i],
            "wpm_rank": w_ranks[i],
            "topsis_score": round(topsis_scores[i], 4),
            "vikor_score": round(vikor_scores[i], 4),
            "wpm_score": round(wpm_scores[i], 4),
        }

    # Rank concordance (Spearman-like: fraction of pairwise agreements)
    def rank_concordance(r1, r2):
        pairs = 0
        agree = 0
        for i in range(n):
            for j in range(i + 1, n):
                pairs += 1
                if (r1[i] < r1[j]) == (r2[i] < r2[j]):
                    agree += 1
                elif r1[i] == r1[j] and r2[i] == r2[j]:
                    agree += 1
        return agree / pairs if pairs > 0 else 1.0

    concordance = {
        "topsis_vikor": round(rank_concordance(t_ranks, v_ranks), 3),
        "topsis_wpm": round(rank_concordance(t_ranks, w_ranks), 3),
        "vikor_wpm": round(rank_concordance(v_ranks, w_ranks), 3),
    }

    # Do all methods agree on top-1?
    t_top1 = t_ranks.index(1) if 1 in t_ranks else -1
    v_top1 = v_ranks.index(1) if 1 in v_ranks else -1
    w_top1 = w_ranks.index(1) if 1 in w_ranks else -1

    return {
        "per_strategy": per_strategy,
        "concordance": concordance,
        "all_agree_on_top1": (t_top1 == v_top1 == w_top1),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Pipeline Orchestrator
# ═══════════════════════════════════════════════════════════════════════════

class StrategyPipeline:
    """End-to-end pipeline: genomic variants → ranked strategies.

    This orchestrator connects all v2 modules in sequence:
    1. TranscriptFetcher → TranscriptInfo
    2. TranscriptMapper → TranscriptCoordinate
    3. ReferenceValidator → ReferenceValidation
    4. CodingAnnotator → CodingAnnotation
    5. VariantNormalizer → NormalizedVariant
    6. PAM Scanner + BE/PE/HDR engines → FeasibilityBundle
    7. StrategyGenerator → List[Strategy]
    8. StrategyScorer → List[ScoredStrategy]

    Parameters
    ----------
    cell_type : str
        Cell type for scoring context (default "iPSC").
    nuclease : str
        Nuclease for PAM scanning (default "SpCas9").
    species : str
        Species for Ensembl API (default "homo_sapiens").
    flank_size : int
        Flanking sequence to fetch around each variant (default 200 bp).
    """

    def __init__(
        self,
        cell_type: str = "iPSC",
        nuclease: str = "SpCas9",
        species: str = "homo_sapiens",
        flank_size: int = 200,
    ):
        self.cell_type = cell_type
        self.nuclease = nuclease
        self.species = species
        self.flank_size = flank_size

        # Lazy imports — these modules are created by other phases
        self._fetcher = None
        self._scorer = StrategyScorer()
        self._topsis_scorer = TOPSISScorer()

    def run(
        self,
        variants: List[GenomicVariantInput],
    ) -> PipelineResult:
        """Execute full pipeline.

        Parameters
        ----------
        variants : List[GenomicVariantInput]
            One or more variants to design strategies for.
            All variants must be in the same gene.

        Returns
        -------
        PipelineResult
            Complete output with ranked strategies.
        """
        if not variants:
            return PipelineResult(
                transcript=_empty_transcript(),
                warnings=["No variants provided."],
            )

        warnings: List[str] = []
        gene_symbol = variants[0].gene_symbol

        # Stage 1: Fetch transcript
        try:
            from core.sequence.fetcher import TranscriptFetcher
            fetcher = TranscriptFetcher(species=self.species)
            transcript_id = variants[0].transcript_id
            if transcript_id:
                transcript = fetcher.fetch_by_transcript_id(transcript_id)
            else:
                transcript = fetcher.fetch_by_gene(gene_symbol)
        except Exception as e:
            logger.error(f"Failed to fetch transcript: {e}")
            return PipelineResult(
                transcript=_empty_transcript(),
                warnings=[f"Transcript fetch failed: {e}"],
            )

        # Stage 2: Normalize variants
        normalized: List[NormalizedVariant] = []
        try:
            from core.sequence.variant_normalizer import VariantNormalizer
            normalizer = VariantNormalizer(species=self.species)
            for v in variants:
                try:
                    nv = normalizer.normalize(v, flank_size=self.flank_size)
                    normalized.append(nv)
                except Exception as e:
                    logger.warning(f"Failed to normalize variant {v.name}: {e}")
                    warnings.append(f"Variant {v.name} normalization failed: {e}")
        except ImportError:
            logger.warning(
                "VariantNormalizer not available. "
                "Using minimal normalization."
            )
            warnings.append("Full normalization not available; using minimal mode.")
            for v in variants:
                normalized.append(_minimal_normalize(v, transcript))

        if not normalized:
            return PipelineResult(
                transcript=transcript,
                warnings=warnings + ["No variants could be normalized."],
            )

        # Stage 3: Run feasibility checks
        bundles: List[FeasibilityBundle] = []
        try:
            from core.feasibility.pam_scan import EnhancedPAMScanner
            from core.feasibility.base_editing import BaseEditingEngine
            from core.feasibility.prime_editing import PrimeEditingEngine
            from core.feasibility.hdr_design import HDRDesignEngine

            scanner = EnhancedPAMScanner(nuclease=self.nuclease)
            be_engine = BaseEditingEngine(pam_scanner=scanner)
            pe_engine = PrimeEditingEngine(pam_scanner=scanner)
            hdr_engine = HDRDesignEngine(
                pam_scanner=scanner,
                nuclease=self.nuclease,
                cell_type=self.cell_type,
            )

            for idx, nv in enumerate(normalized):
                bundle = FeasibilityBundle(
                    variant=nv,
                    mutation_index=idx,
                )

                local_seq = nv.local_sequence
                edit_idx = nv.local_seq_edit_index

                if local_seq and edit_idx >= 0:
                    # Base editing — multi-nuclease evaluation (v3)
                    # Systematically tests all editor-nuclease combinations:
                    # ABE8e/BE4max with SpCas9, enFnCas9, SpCas9-NG, SpRY
                    try:
                        multi_be_results = be_engine.check_all_editors(
                            nv, local_seq, edit_idx
                        )
                        bundle.base_editing_results.extend(multi_be_results)
                    except Exception as e:
                        logger.warning(f"BE multi-editor check failed for variant {idx}: {e}")
                        # Fallback to legacy single-nuclease check
                        try:
                            be_result = be_engine.check_feasibility(
                                nv, local_seq, edit_idx, nuclease=self.nuclease
                            )
                            bundle.base_editing_results.append(be_result)
                        except Exception as e2:
                            logger.warning(f"BE fallback also failed: {e2}")

                    # Prime editing
                    try:
                        pe_result = pe_engine.check_feasibility(
                            nv, local_seq, edit_idx, nuclease=self.nuclease
                        )
                        bundle.prime_editing_result = pe_result
                    except Exception as e:
                        logger.warning(f"PE check failed for variant {idx}: {e}")

                    # HDR
                    try:
                        hdr_result = hdr_engine.check_feasibility(
                            nv, local_seq, edit_idx
                        )
                        bundle.hdr_result = hdr_result
                    except Exception as e:
                        logger.warning(f"HDR check failed for variant {idx}: {e}")
                else:
                    warnings.append(
                        f"No local sequence for variant {idx}; "
                        "feasibility checks skipped."
                    )

                bundles.append(bundle)

        except ImportError as e:
            logger.warning(f"Feasibility modules not available: {e}")
            warnings.append(f"Feasibility modules not available: {e}")
            # Create minimal bundles
            for idx, nv in enumerate(normalized):
                bundles.append(
                    FeasibilityBundle(variant=nv, mutation_index=idx)
                )

        # Stage 4: Generate strategies
        from core.mosaic.generator import StrategyGenerator

        p53_active = True
        try:
            from utils.constants import CELL_TYPE_PARAMS
            ct_params = CELL_TYPE_PARAMS.get(self.cell_type, {})
            p53_active = ct_params.get("p53_active", True)
        except ImportError:
            pass

        generator = StrategyGenerator(p53_active=p53_active)
        all_strategies = generator.generate(bundles)

        # Separate rejected strategies
        rejected = [s for s in all_strategies if s.is_rejected]
        viable = [s for s in all_strategies if not s.is_rejected]

        # Stage 5: Score and rank using TOPSIS (v3) with sensitivity analysis
        ranked = self._topsis_scorer.rank(
            viable, bundles, run_sensitivity=True
        )

        # Also compute legacy weighted-sum scores for comparison
        legacy_ranked = self._scorer.rank(viable, bundles)
        legacy_scores = {
            s.strategy_name: s.overall_score for s in legacy_ranked
        }

        result = PipelineResult(
            transcript=transcript,
            variants=normalized,
            bundles=bundles,
            strategies=ranked,
            rejected_strategies=rejected,
            metadata={
                "cell_type": self.cell_type,
                "nuclease": self.nuclease,
                "scoring_method": "TOPSIS",
                "sensitivity_runs": self._topsis_scorer.n_sensitivity_runs,
                "n_variants": len(normalized),
                "n_strategies_generated": len(all_strategies),
                "n_strategies_rejected": len(rejected),
                "n_strategies_ranked": len(ranked),
                "legacy_weighted_sum_scores": legacy_scores,
            },
            warnings=warnings,
        )

        # Stage 6: Post-ranking delivery annotations
        # Annotates each ranked strategy with delivery feasibility,
        # donor format recommendations, and cell-type-specific warnings.
        # Does NOT change TOPSIS scores or rankings.
        try:
            from core.feasibility.delivery_advisor import DeliveryAdvisor
            delivery_advisor = DeliveryAdvisor(cell_type=self.cell_type)
            delivery_result = delivery_advisor.advise(result)
            result.metadata["delivery_advisory"] = {
                "n_annotated": len(delivery_result.annotations),
                "global_warnings": delivery_result.global_warnings,
                "annotations": [
                    {
                        "strategy": ann.strategy_name,
                        "deliverable": ann.is_deliverable,
                        "method": ann.delivery_method,
                        "complexity": ann.delivery_complexity,
                        "donor_format": (
                            ann.donor_recommendation.format
                            if ann.donor_recommendation else None
                        ),
                        "warnings": ann.warnings,
                        "enhancers": ann.viability_enhancers,
                        "violations": ann.hard_constraint_violations,
                    }
                    for ann in delivery_result.annotations
                ],
            }
            # Add delivery warnings to pipeline warnings
            for ann in delivery_result.annotations:
                for violation in ann.hard_constraint_violations:
                    warnings.append(
                        f"[Delivery] {ann.strategy_name}: {violation}"
                    )
        except Exception as e:
            logger.warning(f"Delivery advisory failed: {e}")
            result.metadata["delivery_advisory"] = {
                "error": str(e),
            }

        return result


# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════

def _empty_transcript() -> TranscriptInfo:
    """Create a placeholder TranscriptInfo for error cases."""
    return TranscriptInfo(
        transcript_id="",
        gene_symbol="",
        gene_id="",
        chromosome="",
        start=0,
        end=0,
        strand=1,
        biotype="",
        is_canonical=False,
        exons=[],
    )


def _minimal_normalize(
    variant: GenomicVariantInput,
    transcript: TranscriptInfo,
) -> NormalizedVariant:
    """Create a minimal NormalizedVariant without full annotation.

    Used when the VariantNormalizer is not available (e.g., during
    early development or when Ensembl API is unreachable).
    """
    from core.models import (
        CodingAnnotation,
        ConsequenceType,
        ReferenceValidation,
        TranscriptCoordinate,
    )

    coord = TranscriptCoordinate(
        genomic_position=variant.position,
        exon_number=-1,
        transcript_position=-1,
        cds_position=-1,
        codon_index=-1,
        codon_position=-1,
        reference_codon="",
        reference_aa="",
        distance_to_exon_start=999,
        distance_to_exon_end=999,
        in_cds=False,
    )

    coding = CodingAnnotation(
        consequence=ConsequenceType.UNKNOWN,
        message="Minimal normalization — full annotation not available.",
    )

    ref_val = ReferenceValidation(
        is_valid=True,
        expected_ref_genomic=variant.ref_allele,
        provided_ref_genomic=variant.ref_allele,
        expected_ref_transcript=variant.ref_allele,
        provided_ref_transcript=variant.ref_allele,
        message="Not validated (minimal mode).",
    )

    return NormalizedVariant(
        input=variant,
        transcript=transcript,
        transcript_coord=coord,
        coding=coding,
        ref_validation=ref_val,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Rank-stability assessment (Fix #3, 2026-04-20)
# ═══════════════════════════════════════════════════════════════════════════


@dataclass
class DimensionDelta:
    """One dimension's top-vs-runner-up comparison."""

    name: str
    top_value: float
    runner_up_value: float
    delta: float           # top - runner_up (signed)
    direction: str         # "benefit" (higher better) or "cost" (lower better)
    prefers: Optional[str] # "top" | "runner_up" | None (tie)


@dataclass
class StabilityAssessment:
    """Interpretation of a top-ranked strategy's rank_stability.

    ``level`` is one of:
      - ``"robust"``      — rank_stability >= RANK_STABILITY_ROBUST
      - ``"stable"``      — in [RANK_STABILITY_STABLE, RANK_STABILITY_ROBUST)
      - ``"flip_sensitive"`` — < RANK_STABILITY_STABLE
      - ``"unknown"``     — sensitivity analysis not run (None)

    ``runner_up``, ``score_gap``, ``dimension_deltas`` and
    ``preferential_reasoning`` are only populated when ``level ==
    "flip_sensitive"`` and at least one alternative strategy exists.
    """

    level: str
    top: ScoredStrategy
    runner_up: Optional[ScoredStrategy] = None
    score_gap: Optional[float] = None
    dimension_deltas: List[DimensionDelta] = field(default_factory=list)
    preferential_reasoning: List[str] = field(default_factory=list)

    @property
    def is_flip_sensitive(self) -> bool:
        return self.level == "flip_sensitive"

    @property
    def human_label(self) -> str:
        return {
            "robust": "ROBUST",
            "stable": "STABLE",
            "flip_sensitive": "FLIP-SENSITIVE",
            "unknown": "UNKNOWN (sensitivity analysis not run)",
        }.get(self.level, self.level.upper())


def _classify_stability_level(rank_stability: Optional[float]) -> str:
    """Map a rank_stability probability to a categorical level."""
    if rank_stability is None:
        return "unknown"
    if rank_stability >= RANK_STABILITY_ROBUST:
        return "robust"
    if rank_stability >= RANK_STABILITY_STABLE:
        return "stable"
    return "flip_sensitive"


def _dimension_comparison(
    top: ScoredStrategy, runner_up: ScoredStrategy
) -> List[DimensionDelta]:
    """Compute per-dimension deltas between top and runner-up scores."""
    # Dimension direction: whether higher = better.
    specs = [
        ("safety",      top.safety_score,      runner_up.safety_score,      "benefit"),
        ("feasibility", top.feasibility_score, runner_up.feasibility_score, "benefit"),
        ("complexity",  top.complexity_score,  runner_up.complexity_score,  "cost"),
        ("risk",        top.risk_score,        runner_up.risk_score,        "cost"),
        ("confidence",  top.confidence_score,  runner_up.confidence_score,  "benefit"),
    ]
    out: List[DimensionDelta] = []
    for name, tv, rv, direction in specs:
        delta = tv - rv
        if abs(delta) < DIMENSION_DELTA_MATERIAL:
            prefers: Optional[str] = None
        elif direction == "benefit":
            prefers = "top" if delta > 0 else "runner_up"
        else:  # cost: lower is better
            prefers = "top" if delta < 0 else "runner_up"
        out.append(
            DimensionDelta(
                name=name,
                top_value=tv,
                runner_up_value=rv,
                delta=delta,
                direction=direction,
                prefers=prefers,
            )
        )
    return out


def _preferential_reasoning(
    top: ScoredStrategy,
    runner_up: ScoredStrategy,
    deltas: List[DimensionDelta],
) -> List[str]:
    """Human-readable lines explaining why each strategy might be preferred."""
    top_wins = [d for d in deltas if d.prefers == "top"]
    runner_up_wins = [d for d in deltas if d.prefers == "runner_up"]

    lines: List[str] = []
    if top_wins:
        dims = ", ".join(
            f"{d.name} ({d.top_value:.2f} vs {d.runner_up_value:.2f})"
            for d in top_wins
        )
        lines.append(f"{top.strategy_name} wins on: {dims}")
    if runner_up_wins:
        dims = ", ".join(
            f"{d.name} ({d.top_value:.2f} vs {d.runner_up_value:.2f})"
            for d in runner_up_wins
        )
        lines.append(f"{runner_up.strategy_name} wins on: {dims}")
    if not top_wins and not runner_up_wins:
        lines.append(
            "No material dimension differences — the top ordering is driven "
            "by fine-grained score differences only. Either choice is "
            "defensible in this context."
        )
    return lines


def assess_stability(strategies: List[ScoredStrategy]) -> Optional[StabilityAssessment]:
    """Interpret a ranked strategy list's top-strategy rank stability.

    Returns ``None`` if the list is empty. If ``rank_stability`` is ``None``
    on the top strategy, returns an assessment with ``level="unknown"``.

    Flip-sensitive assessments include the runner-up strategy, its score
    gap from the top, a per-dimension delta list, and human-readable
    ``preferential_reasoning`` lines suitable for printing in the CLI or
    rendering in the webapp.
    """
    if not strategies:
        return None
    top = strategies[0]
    level = _classify_stability_level(top.rank_stability)

    if level != "flip_sensitive" or len(strategies) < 2:
        return StabilityAssessment(level=level, top=top)

    runner_up = strategies[1]
    deltas = _dimension_comparison(top, runner_up)
    reasoning = _preferential_reasoning(top, runner_up, deltas)
    return StabilityAssessment(
        level=level,
        top=top,
        runner_up=runner_up,
        score_gap=top.overall_score - runner_up.overall_score,
        dimension_deltas=deltas,
        preferential_reasoning=reasoning,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Testing strategy_stage.py...")

    # Test scorer
    scorer = StrategyScorer()
    s = Strategy(
        name="Test BE",
        num_dsbs=0,
        num_rounds=1,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.A,
        modality_prior_score=0.90,
    )
    scored = scorer.score_strategy(s, [])
    assert scored.safety_score == 1.0, "DSB-free should have safety=1.0"
    assert scored.overall_score > 0.5, f"Score too low: {scored.overall_score}"
    print(f"  BE strategy score: {scored.overall_score:.3f}")

    # Test HDR strategy (lower expected score due to DSB)
    s2 = Strategy(
        name="Test HDR",
        num_dsbs=1,
        num_rounds=1,
        num_donors=1,
        rearrangement_risk=RiskLevel.LOW,
        evidence_tier=EvidenceTier.B,
        modality_prior_score=0.72,
        p53_active=True,
    )
    scored2 = scorer.score_strategy(s2, [])
    assert scored2.safety_score < scored.safety_score, "HDR should be less safe than BE"
    print(f"  HDR strategy score: {scored2.overall_score:.3f}")
    print(f"  BE > HDR: {scored.overall_score > scored2.overall_score}")

    # Test ranking
    ranked = scorer.rank([s, s2], [])
    assert ranked[0].rank == 1
    assert ranked[0].strategy_name == "Test BE"
    print(f"  Ranking: {[r.strategy_name for r in ranked]}")

    print("PASS")
