"""
Strategy Scoring and Ranking for MOSAIC
========================================

This module assigns multi-dimensional scores to each editing strategy
and produces a final ranked list.  The scoring system is designed to
reflect the priorities of a molecular biologist working with iPSCs:

  1. **Safety** is paramount.  iPSCs with p53 mutations or chromosomal
     rearrangements are unsuitable for downstream applications
     (differentiation, disease modelling, cell therapy).
  2. **Efficiency** determines screening burden — how many clones must
     be picked and genotyped to find one perfect clone.
  3. **Time** is important because iPSC experiments are slow (each
     editing round takes ~4 weeks including clonal expansion).
  4. **Cost** includes reagent synthesis (donors, guides), screening
     (genotyping), and cell culture.

Scoring philosophy
------------------
Each dimension is normalised to [0, 1] where 1 is best.  The overall
score is a weighted average:

    overall = w_eff * efficiency_score
            + w_safe * safety_score
            + w_time * time_score
            + w_cost * cost_score

Default weights for iPSC work:
    w_safe = 0.40   (safety is the top priority)
    w_eff  = 0.30   (efficiency determines practical feasibility)
    w_time = 0.15   (time matters but is secondary)
    w_cost = 0.15   (cost is important but usually not the bottleneck)

For HEK293T or other immortalised lines, safety can be down-weighted
because p53 selection and karyotype stability are less critical.

Key safety considerations for iPSCs
-------------------------------------
1. **p53 selection**: Every DSB activates the p53 pathway.  In iPSCs,
   this causes massive apoptosis (up to 50-70% cell death per DSB).
   Surviving clones may have inactivating TP53 mutations acquired
   under this selective pressure, rendering them unsafe for clinical
   or disease-modelling use.
   Source: Ihry et al., "p53 inhibits CRISPR-Cas9 engineering in
   human pluripotent stem cells," Nature Medicine, 2018.

2. **Translocation risk**: Two simultaneous DSBs can generate balanced
   or unbalanced translocations.  On the same chromosome, the
   intervening segment can be deleted (interstitial deletion) or
   inverted.  These events are undetectable by standard PCR genotyping
   and require karyotyping, FISH, or whole-genome sequencing.
   Source: Kosicki et al., Nature Biotechnology, 2018;
   Leibowitz et al., Nature Genetics, 2021.

3. **Large deletions**: Even a single DSB can cause kilobase-scale
   deletions centered on the cut site, detectable by long-read
   sequencing but often missed by short-range PCR genotyping.
   Source: Kosicki et al., Nature Biotechnology, 2018.

4. **Chromothripsis**: Rarely, a DSB can trigger chromothripsis-like
   rearrangements, especially in p53-deficient backgrounds.
   Source: Leibowitz et al., Nature Genetics, 2021.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .strategy_enumerator import EditingStrategy

# Import cell-type parameters for context-dependent scoring.
try:
    from ..utils.constants import CELL_TYPE_PARAMS
except ImportError:
    CELL_TYPE_PARAMS = {
        "iPSC": {
            "hdr_base_efficiency": 0.08,
            "p53_active": True,
            "viability_single_dsb": 0.55,
            "viability_dual_dsb": 0.30,
        },
        "HEK293T": {
            "hdr_base_efficiency": 0.25,
            "p53_active": False,
            "viability_single_dsb": 0.85,
            "viability_dual_dsb": 0.70,
        },
    }


# ===========================================================================
# ScoredStrategy: a strategy with scores attached
# ===========================================================================

@dataclass
class ScoredStrategy:
    """An ``EditingStrategy`` with numerical scores attached.

    Attributes
    ----------
    strategy : EditingStrategy
        The underlying strategy.
    efficiency_score : float
        0-1 score for predicted editing efficiency.
    safety_score : float
        0-1 score for biological safety.
    time_score : float
        0-1 score inversely proportional to weeks needed.
    cost_score : float
        0-1 score based on reagent and screening costs.
    overall_score : float
        Weighted combination of all four scores.
    rank : int
        1-based rank (1 = best).
    """

    strategy: EditingStrategy
    efficiency_score: float = 0.0
    safety_score: float = 0.0
    time_score: float = 0.0
    cost_score: float = 0.0
    overall_score: float = 0.0
    rank: int = 0


# ===========================================================================
# StrategyScorer
# ===========================================================================

class StrategyScorer:
    """Score and rank editing strategies.

    The scorer evaluates each strategy across four dimensions and
    computes a weighted overall score.  Strategies are then ranked from
    best (rank 1) to worst.

    Parameters
    ----------
    weight_efficiency : float
        Weight for efficiency score in the overall calculation.
        Default 0.30.
    weight_safety : float
        Weight for safety score.  Default 0.40 (highest for iPSC work).
    weight_time : float
        Weight for time score.  Default 0.15.
    weight_cost : float
        Weight for cost score.  Default 0.15.

    Notes
    -----
    Weights must sum to 1.0.  If they do not, they are normalised
    automatically.
    """

    def __init__(
        self,
        weight_efficiency: float = 0.30,
        weight_safety: float = 0.40,
        weight_time: float = 0.15,
        weight_cost: float = 0.15,
    ) -> None:
        total = weight_efficiency + weight_safety + weight_time + weight_cost
        # Normalise weights so they sum to 1.0.
        self.w_eff: float = weight_efficiency / total
        self.w_safe: float = weight_safety / total
        self.w_time: float = weight_time / total
        self.w_cost: float = weight_cost / total

    # -----------------------------------------------------------------
    # Scoring individual strategies
    # -----------------------------------------------------------------

    def score(
        self,
        strategy: EditingStrategy,
        cell_type: str = "iPSC",
    ) -> ScoredStrategy:
        """Score a single strategy across all dimensions.

        Parameters
        ----------
        strategy : EditingStrategy
            The strategy to score.
        cell_type : str
            Cell type context (affects safety scoring).

        Returns
        -------
        ScoredStrategy
            The strategy with all scores populated.
        """
        ct = CELL_TYPE_PARAMS.get(cell_type, CELL_TYPE_PARAMS.get("iPSC"))

        eff_score = self._score_efficiency(strategy)
        safe_score = self._score_safety(strategy, ct)
        time_score = self._score_time(strategy)
        cost_score = self._score_cost(strategy)

        overall = (
            self.w_eff * eff_score
            + self.w_safe * safe_score
            + self.w_time * time_score
            + self.w_cost * cost_score
        )

        return ScoredStrategy(
            strategy=strategy,
            efficiency_score=round(eff_score, 4),
            safety_score=round(safe_score, 4),
            time_score=round(time_score, 4),
            cost_score=round(cost_score, 4),
            overall_score=round(overall, 4),
        )

    # -----------------------------------------------------------------
    # Dimension-specific scoring functions
    # -----------------------------------------------------------------

    def _score_efficiency(self, strategy: EditingStrategy) -> float:
        """Compute efficiency score (0-1).

        The efficiency score directly reflects the predicted fraction of
        transfected cells that will yield a correctly edited clone at ALL
        target sites.  Higher is better.

        Mapping:
            estimated_efficiency = 0.0  ->  score = 0.0
            estimated_efficiency = 1.0  ->  score = 1.0

        We apply a square-root transform to compress the dynamic range,
        because the difference between 1% and 5% efficiency is far more
        meaningful (in terms of screening burden) than the difference
        between 50% and 54%.

        The square root maps:
            0.01 -> 0.10
            0.05 -> 0.22
            0.10 -> 0.32
            0.25 -> 0.50
            0.50 -> 0.71
        """
        eff = np.clip(strategy.estimated_efficiency, 0.0, 1.0)
        return float(np.sqrt(eff))

    def _score_safety(
        self,
        strategy: EditingStrategy,
        cell_type_params: Dict,
    ) -> float:
        """Compute safety score (0-1).

        Safety scoring accounts for multiple risk factors:

        1. **Number of DSBs** (most important factor):
           - 0 DSBs: score starts at 1.0 (DSB-free editing is safest)
           - 1 DSB: score starts at 0.6 (single DSB carries p53 risk
             in iPSCs and large-deletion risk)
           - 2+ DSBs: score starts at 0.3 (compounded risk)
           Source: General principle; Ihry et al., 2018 for iPSC-specific
           p53 toxicity data.

        2. **p53 activation penalty** (iPSC-specific):
           - If the cell type has active p53 AND the strategy introduces
             DSBs, apply an additional 30% reduction.
           - Rationale: p53 activation causes massive cell death AND
             selects for p53-inactivated clones, which are genomically
             unstable and unsafe for clinical use.
           Source: Ihry et al., Nature Medicine, 2018.

        3. **Translocation risk penalty** (for simultaneous 2+ DSBs):
           - If 2+ DSBs are introduced in ONE round, apply an additional
             penalty proportional to the translocation probability.
           - Even if the estimated translocation rate is low (~0.1%),
             the *consequence* is catastrophic (chromosomal rearrangement
             undetectable by standard PCR).
           Source: Kosicki et al., Nature Biotechnology, 2018;
           Leibowitz et al., Nature Genetics, 2021.

        4. **Simultaneous DSB bonus/penalty**:
           - Sequential HDR (one DSB per round, never simultaneous)
             gets a small bonus relative to simultaneous dual HDR.

        Parameters
        ----------
        strategy : EditingStrategy
            The strategy to evaluate.
        cell_type_params : dict
            Cell-type-specific parameters from ``CELL_TYPE_PARAMS``.

        Returns
        -------
        float
            Safety score in [0, 1].
        """
        # Start with a base score determined by DSB count.
        if strategy.num_dsbs == 0:
            # DSB-free approaches (base editing, prime editing).
            # These do not activate the DNA-damage checkpoint and do not
            # risk large deletions or translocations.
            base_safety = 1.0
        elif strategy.num_dsbs == 1:
            # Single DSB: moderate risk.
            # Even a single DSB can cause:
            #   - Large deletions at the cut site (~10-20% of alleles
            #     may have >100 bp deletions; Kosicki et al., 2018)
            #   - p53 activation in iPSCs
            #   - On-target mutagenesis via NHEJ on the unedited allele
            base_safety = 0.6
        else:
            # Two or more DSBs.
            base_safety = 0.3

        # --- p53 penalty ---
        # In cells with active p53 (iPSCs, primary cells), DSBs trigger
        # apoptosis and select for p53-null clones.
        p53_penalty = 0.0
        if cell_type_params.get("p53_active", False) and strategy.num_dsbs > 0:
            # The penalty increases with DSB count.
            # 1 DSB: 20% penalty (moderate p53 activation)
            # 2 DSBs: 35% penalty (strong p53 activation)
            # Source: Ihry et al., 2018 — ~50% cell death per DSB in iPSCs
            p53_penalty = 0.15 + 0.10 * (strategy.num_dsbs - 1)

        # --- Translocation risk penalty ---
        # Only applies when 2+ DSBs are introduced SIMULTANEOUSLY
        # (in the same editing round).
        translocation_penalty = 0.0
        if strategy.num_dsbs >= 2 and strategy.num_editing_rounds == 1:
            # Simultaneous DSBs on the same chromosome: high risk of
            # interstitial deletion.
            # We apply a flat penalty because the exact translocation
            # probability depends on genomic distance (which is encoded
            # in the strategy's feasibility_warnings but not directly
            # accessible here as a number).  The enumerator already set
            # the strategy.safety_score with distance-aware calculations;
            # here we add an additional penalty for the simultaneous
            # nature of the DSBs.
            translocation_penalty = 0.15

        # --- Sequential DSB bonus ---
        # If DSBs occur in separate rounds (sequential editing), there
        # is no translocation risk (the first DSB is fully repaired
        # before the second is introduced).
        sequential_bonus = 0.0
        if (strategy.num_dsbs >= 2
                and strategy.num_editing_rounds >= 2
                and strategy.requires_dsb):
            # Sequential HDR: safer than simultaneous
            sequential_bonus = 0.10

        safety = (
            base_safety
            - p53_penalty
            - translocation_penalty
            + sequential_bonus
        )

        # Clamp to [0.05, 1.0] — nothing gets a perfect 0 safety score
        # because even the riskiest strategy has *some* probability of
        # producing a karyotypically normal clone.
        return float(np.clip(safety, 0.05, 1.0))

    def _score_time(self, strategy: EditingStrategy) -> float:
        """Compute time score (0-1).  Shorter is better.

        Mapping:
            4 weeks (1 round)  -> 1.0
            8 weeks (2 rounds) -> 0.5
            12 weeks           -> 0.33
            16 weeks           -> 0.25

        We use an inverse mapping: score = 4 / time_weeks.
        The reference (best-case) is 4 weeks = 1 editing round.

        Biology
        -------
        Each editing round in iPSCs takes approximately 4 weeks:
          Week 1: electroporation + recovery (48 h) + initial growth
          Week 2: single-cell plating + clonal expansion starts
          Week 3-4: colony picking, genomic DNA extraction, genotyping
                    (PCR + Sanger sequencing or amplicon-seq)

        Additional time for validation (karyotyping, off-target analysis)
        is not included here but should be budgeted separately (~2-4
        weeks per validated clone).

        Source: Standard iPSC editing workflows (Ran et al., Nat Protocols,
        2013; Bruntraeger et al., Methods, 2019).
        """
        if strategy.time_weeks <= 0:
            return 1.0
        # Reference time: 4 weeks (one round)
        reference_weeks = 4.0
        score = reference_weeks / strategy.time_weeks
        return float(np.clip(score, 0.0, 1.0))

    def _score_cost(self, strategy: EditingStrategy) -> float:
        """Compute cost score (0-1).  Lower cost is better.

        Cost drivers considered:
          1. **Donor templates**: cssDNA synthesis costs ~$500-2000 per
             donor depending on length.  ssODNs are ~$50.  No cost for
             base/prime editing (no donor).
          2. **Number of editing rounds**: each round costs reagents
             (electroporation supplies, culture media, genotyping) —
             approximately $500-1000 per round.
          3. **Screening burden**: each clone screened costs ~$10-20 for
             PCR genotyping.  High screening burdens (>100 clones)
             become expensive.
          4. **Selection markers**: if required, add ~$500 for selection
             reagents + time.

        We use a simple cost model:
            cost_units = (donor_templates * 2)
                       + (editing_rounds * 1)
                       + (screening_burden / 50)  # normalised
                       + (1 if selection else 0)

        Then: score = 1 / (1 + cost_units / 5)

        This gives:
            cost_units = 0 -> score ~ 1.0
            cost_units = 5 -> score ~ 0.5
            cost_units = 15 -> score ~ 0.25

        Source: Approximate costs from commercial CRISPR reagent
        suppliers (IDT, Genscript, Addgene) as of 2024-2025.
        """
        cost_units = 0.0

        # Donor templates: ~2 cost units each (cssDNA is expensive)
        cost_units += strategy.donor_templates_needed * 2.0

        # Editing rounds: ~1 cost unit each (reagents + culture)
        cost_units += strategy.num_editing_rounds * 1.0

        # Screening burden: normalised by 50 clones
        # (screening 50 clones is "normal"; >100 is expensive)
        cost_units += strategy.screening_burden / 50.0

        # Selection marker
        if strategy.requires_selection:
            cost_units += 1.0

        # Convert to 0-1 score via inverse mapping
        score = 1.0 / (1.0 + cost_units / 5.0)

        return float(np.clip(score, 0.0, 1.0))

    # -----------------------------------------------------------------
    # Ranking
    # -----------------------------------------------------------------

    def rank_strategies(
        self,
        strategies: List[EditingStrategy],
        cell_type: str = "iPSC",
    ) -> List[ScoredStrategy]:
        """Score and rank a list of strategies.

        Parameters
        ----------
        strategies : list of EditingStrategy
            Strategies to evaluate (from ``StrategyEnumerator``).
        cell_type : str
            Cell type context for scoring.

        Returns
        -------
        list of ScoredStrategy
            Strategies sorted by overall_score (descending).
            ``rank`` is set to 1 for the best strategy.

        Notes
        -----
        Ties are broken by safety_score (prefer safer strategies when
        overall scores are equal).
        """
        scored: List[ScoredStrategy] = []
        for s in strategies:
            scored.append(self.score(s, cell_type=cell_type))

        # Sort by overall_score descending, then safety_score descending
        scored.sort(
            key=lambda x: (x.overall_score, x.safety_score),
            reverse=True,
        )

        # Assign ranks
        for i, ss in enumerate(scored):
            ss.rank = i + 1

        return scored

    def format_ranking(
        self,
        ranked: List[ScoredStrategy],
    ) -> str:
        """Format the ranking as a human-readable table.

        Parameters
        ----------
        ranked : list of ScoredStrategy
            Ranked strategies (from ``rank_strategies``).

        Returns
        -------
        str
            Formatted table string.
        """
        header = (
            f"{'Rank':<5} {'Strategy':<40} {'Overall':>8} "
            f"{'Safety':>7} {'Effic.':>7} {'Time':>5} {'Cost':>5} "
            f"{'DSBs':>5} {'Weeks':>5} {'Screen':>7}"
        )
        sep = "-" * len(header)
        lines = [header, sep]

        for ss in ranked:
            s = ss.strategy
            lines.append(
                f"{ss.rank:<5} {s.name:<40} "
                f"{ss.overall_score:>8.3f} "
                f"{ss.safety_score:>7.3f} "
                f"{ss.efficiency_score:>7.3f} "
                f"{ss.time_score:>5.2f} "
                f"{ss.cost_score:>5.2f} "
                f"{s.num_dsbs:>5} "
                f"{s.time_weeks:>5} "
                f"{s.screening_burden:>7}"
            )

        return "\n".join(lines)


# =========================================================================
# Quick self-test
# =========================================================================

if __name__ == "__main__":
    from .gene_structure import build_demo_gene
    from .mutation_classifier import Mutation
    from .strategy_enumerator import StrategyEnumerator

    gene = build_demo_gene()
    m1 = Mutation(exon_number=40, position=5_000_000, ref_allele="G",
                  alt_allele="A", name="c.5000G>A")
    m2 = Mutation(exon_number=60, position=8_000_000, ref_allele="CT",
                  alt_allele="-", name="c.8000delCT")

    enumerator = StrategyEnumerator()
    strategies = enumerator.enumerate_strategies(
        gene, [m1, m2], cell_type="iPSC", nuclease="SpCas9"
    )

    scorer = StrategyScorer()
    ranked = scorer.rank_strategies(strategies, cell_type="iPSC")
    print(scorer.format_ranking(ranked))
