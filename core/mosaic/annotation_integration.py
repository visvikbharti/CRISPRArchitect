"""
Consequence-Aware Score Adjustments
=====================================

This module injects biological consequence information (coding effect,
splice proximity, bystander risk) into the strategy scoring pipeline.

The core idea — and the central claim of the CRISPRArchitect v2 paper — is
that strategies which are technically feasible (PAM available, editing
window compatible) may still be biologically suboptimal if they introduce
deleterious bystander mutations or disrupt splice sites.

Consequence-aware adjustments modify the base score from v1's
StrategyScorer to account for these biological realities.

Penalty / bonus design
-----------------------
All adjustments are additive to the 0-1 overall score, then clamped.

Penalties (reduce score):
  - Bystander creates missense: -0.10 per affected position
  - Bystander creates nonsense: -0.25 (disqualifying if in essential gene)
  - Target or bystander near splice donor (<=2bp): -0.15
  - Target or bystander near splice acceptor (<=2bp): -0.15
  - Target or bystander in splice region (3-8bp): -0.08
  - Multiple simultaneous DSBs in p53-active cells: -0.10

Bonuses (increase score):
  - All bystanders are synonymous: +0.05 (clean edit)
  - PAM-disrupting silent mutation in donor: +0.03 (prevents re-cutting)
  - Short cut-to-edit distance (<10bp): +0.05

Literature basis
-----------------
- Arbab et al., Nature, 2020: bystander editing outcomes are predictable
  and often deleterious
- Kim et al., Nature Biotechnology, 2017: CBE bystander effects depend
  on editing window and local sequence context
- Gaudelli et al., Nature, 2017: ABE has narrower window, fewer bystanders
- Rees & Liu, Nature Reviews Genetics, 2018: comparison of editor fidelity
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from core.models import (
    ConsequenceType,
    EditModality,
    FeasibilityBundle,
    FeasibilityLabel,
    NormalizedVariant,
    RiskLevel,
    ScoredStrategy,
    Strategy,
)


# ═══════════════════════════════════════════════════════════════════════════
# Penalty / Bonus Constants
# ═══════════════════════════════════════════════════════════════════════════

# Consequence-based penalties (additive, negative = worse)
PENALTY_BYSTANDER_MISSENSE = -0.10
PENALTY_BYSTANDER_NONSENSE = -0.25
PENALTY_SPLICE_DONOR = -0.15
PENALTY_SPLICE_ACCEPTOR = -0.15
PENALTY_SPLICE_REGION = -0.08
PENALTY_DUAL_DSB_P53 = -0.10

# Consequence-based bonuses (additive, positive = better)
BONUS_ALL_BYSTANDERS_SYNONYMOUS = 0.05
BONUS_PAM_DISRUPTION = 0.03
BONUS_SHORT_CUT_TO_EDIT = 0.05

# Severity mapping for bystander consequences
CONSEQUENCE_PENALTY_MAP = {
    ConsequenceType.NONSENSE: PENALTY_BYSTANDER_NONSENSE,
    ConsequenceType.MISSENSE: PENALTY_BYSTANDER_MISSENSE,
    ConsequenceType.SPLICE_DONOR: PENALTY_SPLICE_DONOR,
    ConsequenceType.SPLICE_ACCEPTOR: PENALTY_SPLICE_ACCEPTOR,
    ConsequenceType.SPLICE_REGION: PENALTY_SPLICE_REGION,
    ConsequenceType.FRAMESHIFT: PENALTY_BYSTANDER_NONSENSE,
}


class AnnotationIntegrator:
    """Adjust strategy scores based on biological consequences.

    This class bridges the v2 consequence annotation (from CodingAnnotator)
    with the v2 scoring system. It takes a Strategy with its associated
    FeasibilityBundles and computes consequence-aware adjustments.

    The adjusted score is clamped to [0.0, 1.0].

    Parameters
    ----------
    p53_active : bool
        Whether the cell type has active p53 (iPSCs do; HEK293T does not).
        Affects the dual-DSB penalty.
    """

    def __init__(self, p53_active: bool = True):
        self.p53_active = p53_active

    def compute_adjustments(
        self,
        strategy: Strategy,
        bundles: List[FeasibilityBundle],
    ) -> AnnotationAdjustment:
        """Compute consequence-aware score adjustments for a strategy.

        Parameters
        ----------
        strategy : Strategy
            The strategy to evaluate.
        bundles : List[FeasibilityBundle]
            Feasibility bundles for each variant in the strategy.

        Returns
        -------
        AnnotationAdjustment
            Contains total penalty, total bonus, and detailed notes.
        """
        penalties: List[float] = []
        bonuses: List[float] = []
        notes: List[str] = []

        # 1. Bystander penalties for base editing steps
        for step in strategy.steps:
            if step.modality in (EditModality.ABE, EditModality.CBE):
                bundle_idx = step.target_mutation_index
                if bundle_idx < len(bundles):
                    bundle = bundles[bundle_idx]
                    be = bundle.best_base_editing_result()
                    if be and be.bystander_count > 0:
                        # Check bystander consequences
                        for csq in be.bystander_consequences:
                            penalty = CONSEQUENCE_PENALTY_MAP.get(csq, 0.0)
                            if penalty < 0:
                                penalties.append(penalty)
                                notes.append(
                                    f"Bystander {csq.value} penalty: "
                                    f"{penalty:.2f}"
                                )
                        # Bonus if all bystanders are synonymous
                        if (
                            be.bystander_count > 0
                            and all(
                                c == ConsequenceType.SYNONYMOUS
                                for c in be.bystander_consequences
                            )
                        ):
                            bonuses.append(BONUS_ALL_BYSTANDERS_SYNONYMOUS)
                            notes.append(
                                "All bystanders synonymous: "
                                f"+{BONUS_ALL_BYSTANDERS_SYNONYMOUS:.2f}"
                            )

        # 2. Splice proximity penalties for all variants in the strategy
        for bundle in bundles:
            variant = bundle.variant
            splice_prox = variant.coding.splice_proximity
            if splice_prox == "near_exon_start":
                penalties.append(PENALTY_SPLICE_ACCEPTOR)
                notes.append(
                    f"Variant near splice acceptor: "
                    f"{PENALTY_SPLICE_ACCEPTOR:.2f}"
                )
            elif splice_prox == "near_exon_end":
                penalties.append(PENALTY_SPLICE_DONOR)
                notes.append(
                    f"Variant near splice donor: "
                    f"{PENALTY_SPLICE_DONOR:.2f}"
                )

        # 3. Dual DSB penalty in p53-active cells
        if strategy.num_dsbs >= 2 and self.p53_active:
            penalties.append(PENALTY_DUAL_DSB_P53)
            notes.append(
                f"Dual DSB in p53-active cells: "
                f"{PENALTY_DUAL_DSB_P53:.2f}"
            )

        # 4. HDR-specific bonuses
        for step in strategy.steps:
            if step.modality in (
                EditModality.HDR_CSSDNA,
                EditModality.HDR_SSODN,
                EditModality.HDR_LSSDNA,
                EditModality.HDR_DSDNA,
            ):
                bundle_idx = step.target_mutation_index
                if bundle_idx < len(bundles):
                    bundle = bundles[bundle_idx]
                    hdr = bundle.hdr_result
                    if hdr:
                        if hdr.cut_to_edit_distance <= 10:
                            bonuses.append(BONUS_SHORT_CUT_TO_EDIT)
                            notes.append(
                                f"Short cut-to-edit ({hdr.cut_to_edit_distance}bp): "
                                f"+{BONUS_SHORT_CUT_TO_EDIT:.2f}"
                            )
                        if hdr.pam_disruption_possible:
                            bonuses.append(BONUS_PAM_DISRUPTION)
                            notes.append(
                                f"PAM disruption possible: "
                                f"+{BONUS_PAM_DISRUPTION:.2f}"
                            )

        total_penalty = sum(penalties)
        total_bonus = sum(bonuses)

        return AnnotationAdjustment(
            total_penalty=total_penalty,
            total_bonus=total_bonus,
            net_adjustment=total_penalty + total_bonus,
            notes=notes,
        )

    def adjust_score(
        self,
        base_score: float,
        adjustment: AnnotationAdjustment,
    ) -> float:
        """Apply consequence adjustment to a base score.

        Parameters
        ----------
        base_score : float
            The score before consequence adjustment (from v1 scorer or
            modality_prior_score).
        adjustment : AnnotationAdjustment
            The computed adjustment.

        Returns
        -------
        float
            Adjusted score, clamped to [0.0, 1.0].
        """
        adjusted = base_score + adjustment.net_adjustment
        return max(0.0, min(1.0, adjusted))


@dataclass
class AnnotationAdjustment:
    """Result of consequence-aware score adjustment computation."""
    total_penalty: float = 0.0
    total_bonus: float = 0.0
    net_adjustment: float = 0.0
    notes: List[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Testing annotation_integration.py...")

    integrator = AnnotationIntegrator(p53_active=True)

    # Test that a clean base score passes through
    adj = AnnotationAdjustment()
    assert integrator.adjust_score(0.80, adj) == 0.80

    # Test penalty application
    adj = AnnotationAdjustment(total_penalty=-0.15, net_adjustment=-0.15)
    assert abs(integrator.adjust_score(0.80, adj) - 0.65) < 1e-9

    # Test clamping at 0
    adj = AnnotationAdjustment(total_penalty=-1.0, net_adjustment=-1.0)
    assert integrator.adjust_score(0.50, adj) == 0.0

    # Test clamping at 1
    adj = AnnotationAdjustment(total_bonus=0.5, net_adjustment=0.5)
    assert integrator.adjust_score(0.80, adj) == 1.0

    print("PASS")
