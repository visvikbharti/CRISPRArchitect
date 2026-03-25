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


# ═══════════════════════════════════════════════════════════════════════════
# Scoring Engine
# ═══════════════════════════════════════════════════════════════════════════

class StrategyScorer:
    """Multi-objective scoring engine for v2 strategies.

    Scoring function:
        Score = w1*Safety + w2*Feasibility - w3*Complexity - w4*Risk + w5*Confidence

    Default weights (iPSC-optimized):
        Safety:     0.30  (DSB-free is safer)
        Feasibility: 0.25 (PAM availability, window compatibility)
        Complexity:  0.20 (donors, rounds, screening)
        Risk:        0.15 (rearrangement, bystander, splice)
        Confidence:  0.10 (evidence tier, clean design)
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
        """0-1 where 1 = safest (no DSBs, no rearrangement risk)."""
        if s.num_dsbs == 0:
            return 1.0
        if s.num_dsbs == 1 and not s.simultaneous_dsbs:
            base = 0.6
            if s.p53_active:
                base -= 0.1  # p53 penalty
            return max(0.0, base)
        if s.num_dsbs >= 2:
            base = 0.3
            if s.simultaneous_dsbs:
                base -= 0.1  # simultaneous is riskier
            if s.p53_active:
                base -= 0.1
            return max(0.0, base)
        return 0.5

    def _score_feasibility(self, s: Strategy) -> float:
        """0-1 based on modality prior score and donor quality."""
        return s.modality_prior_score * s.donor_feasibility_score

    def _score_complexity(self, s: Strategy) -> float:
        """0-1 where 1 = most complex (worst). Lower is better."""
        # Factors: rounds, donors, distinct guides, screening
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
        """0-1 where 1 = highest risk (worst). Lower is better."""
        risk = 0.0

        # Rearrangement risk
        risk_map = {
            RiskLevel.LOW: 0.0,
            RiskLevel.MODERATE: 0.3,
            RiskLevel.HIGH: 0.6,
            RiskLevel.VERY_HIGH: 0.9,
        }
        risk += risk_map.get(s.rearrangement_risk, 0.0)

        # Bystander risk
        risk += s.bystander_severity * 0.3

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
        """Penalty for splice-proximal edits or damaging bystanders."""
        penalty = 0.0

        for bundle in bundles:
            v = bundle.variant
            if v.coding.splice_proximity == "near_exon_start":
                penalty += 0.05
            elif v.coding.splice_proximity == "near_exon_end":
                penalty += 0.05

        # Bystander severity is already encoded in the strategy
        penalty += s.bystander_severity * 0.08

        return min(0.3, penalty)  # cap total penalty

    def _compute_consequence_bonus(
        self, s: Strategy, bundles: List[FeasibilityBundle]
    ) -> float:
        """Bonus for clean designs."""
        bonus = 0.0

        # No DSBs = inherently clean
        if s.num_dsbs == 0 and s.bystander_severity == 0.0:
            bonus += 0.03

        return min(0.1, bonus)


# ═══════════════════════════════════════════════════════════════════════════
# TOPSIS Scorer (v3)
# ═══════════════════════════════════════════════════════════════════════════

import math
import random as _random


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

    # Benefit dimensions (higher = better): safety, feasibility, confidence
    # Cost dimensions (lower = better): complexity, risk
    _BENEFIT_DIMS = {0, 1, 4}  # safety, feasibility, confidence
    _COST_DIMS = {2, 3}        # complexity, risk
    _DIM_NAMES = ["safety", "feasibility", "complexity", "risk", "confidence"]

    def __init__(
        self,
        w_safety: float = 0.30,
        w_feasibility: float = 0.25,
        w_complexity: float = 0.20,
        w_risk: float = 0.15,
        w_confidence: float = 0.10,
        n_sensitivity_runs: int = 10000,
    ):
        total = w_safety + w_feasibility + w_complexity + w_risk + w_confidence
        self.weights = [
            w_safety / total,
            w_feasibility / total,
            w_complexity / total,
            w_risk / total,
            w_confidence / total,
        ]
        self.n_sensitivity_runs = n_sensitivity_runs
        # Also keep a weighted-sum scorer for consequence adjustments
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
        dim_matrix = []  # List of [safety, feas, complex, risk, conf]
        scored_list = []
        for s in active:
            scored = self._legacy_scorer.score_strategy(s, bundles)
            scored_list.append(scored)
            dim_matrix.append([
                scored.safety_score,
                scored.feasibility_score,
                scored.complexity_score,
                scored.risk_score,
                scored.confidence_score,
            ])

        if len(active) == 1:
            # Single strategy: skip TOPSIS, just use legacy score
            scored_list[0].rank = 1
            return scored_list

        # Step 2: TOPSIS ranking
        topsis_scores = self._topsis(dim_matrix, self.weights)

        # Step 3: Apply consequence adjustments to TOPSIS score
        for i, (s, scored) in enumerate(zip(active, scored_list)):
            consequence_adj = scored.consequence_bonus - scored.consequence_penalty
            # TOPSIS score is in [0, 1]; consequence adjustment is small (<0.1)
            adjusted = max(0.0, min(1.0, topsis_scores[i] + consequence_adj))
            scored.overall_score = round(adjusted, 4)

        # Step 4: Sort by TOPSIS score
        scored_list.sort(
            key=lambda x: (x.overall_score, x.safety_score),
            reverse=True,
        )
        for i, s in enumerate(scored_list, 1):
            s.rank = i

        # Step 5: Sensitivity analysis
        if run_sensitivity and len(active) >= 2:
            sensitivity = self._sensitivity_analysis(
                dim_matrix, scored_list
            )
            for scored in scored_list:
                sr = sensitivity.get(scored.strategy_name)
                if sr:
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

    def _sensitivity_analysis(
        self,
        dim_matrix: List[List[float]],
        scored_list: List[ScoredStrategy],
    ) -> Dict[str, SensitivityResult]:
        """Monte Carlo sensitivity analysis: perturb weights and re-rank.

        Generates n_sensitivity_runs random weight vectors from a Dirichlet
        distribution centered on the default weights, runs TOPSIS with each,
        and counts how often each strategy is top-ranked.

        Uses Dirichlet(alpha) where alpha = default_weight * concentration.
        Concentration=10 means weights vary moderately around defaults.
        """
        n_alts = len(dim_matrix)
        if n_alts < 2:
            return {}

        strategy_names = [s.strategy_name for s in scored_list]
        rank_counts = {name: {} for name in strategy_names}  # type: Dict[str, Dict[int, int]]
        top1_counts = {name: 0 for name in strategy_names}

        # Dirichlet concentration: higher = closer to default weights
        concentration = 10.0
        alphas = [w * concentration for w in self.weights]

        rng = _random.Random(42)  # deterministic for reproducibility

        for _ in range(self.n_sensitivity_runs):
            # Sample from Dirichlet by sampling Gamma and normalizing
            raw = []
            for a in alphas:
                # Gamma sampling via Marsaglia's method (stdlib random.gammavariate)
                raw.append(rng.gammavariate(a, 1.0))
            total = sum(raw)
            perturbed_weights = [r / total for r in raw]

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

        # Build results
        results = {}
        n = self.n_sensitivity_runs
        for name in strategy_names:
            counts = rank_counts[name]
            mean_rank = sum(r * c for r, c in counts.items()) / n
            rank_dist = {r: c / n for r, c in sorted(counts.items())}
            results[name] = SensitivityResult(
                strategy_name=name,
                rank_stability=top1_counts[name] / n,
                mean_rank=round(mean_rank, 2),
                rank_distribution=rank_dist,
            )

        return results


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

        return PipelineResult(
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
