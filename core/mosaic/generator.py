"""
Strategy Generator for CRISPRArchitect v2
==========================================

This module takes feasibility bundles (one per variant) and generates
all biologically plausible editing strategies, including:

  - Single-mutation strategies: single BE, single PE, single HDR
  - Two-mutation strategies: dual BE, dual PE, sequential HDR,
    hybrid BE+HDR, hybrid PE+HDR, hybrid BE+PE

Each strategy is represented as a Strategy object with StrategyStep
components, consequence-aware scoring, and explicit rejection reasons
for infeasible combinations.

Design philosophy
------------------
The generator is intentionally exhaustive: it creates every strategy
that could work, then lets the scorer rank them. Strategies that violate
hard constraints (no PAM available, wrong mutation type for BE, etc.)
are tagged with rejection_reasons rather than silently omitted, because
explicit rejection is one of the key contributions of the paper.

Integration with v1
--------------------
This module does NOT call v1's StrategyEnumerator directly. Instead,
it implements its own strategy logic that is informed by v2's PAM-verified
feasibility results. The output Strategy objects are self-contained and
do not depend on v1's EditingStrategy dataclass.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from core.models import (
    EditModality,
    EvidenceTier,
    FeasibilityBundle,
    FeasibilityLabel,
    HDRFeasibility,
    BaseEditingFeasibility,
    NormalizedVariant,
    PrimeEditingFeasibility,
    RiskLevel,
    Strategy,
    StrategyStep,
    ScoredStrategy,
)
from core.mosaic.annotation_integration import AnnotationIntegrator


# ═══════════════════════════════════════════════════════════════════════════
# Modality prior scores — rationale
# ═══════════════════════════════════════════════════════════════════════════
#
# These encode the baseline feasibility of each modality BEFORE locus-specific
# PAM/window checks. They feed into the TOPSIS feasibility dimension.
#
# Evidence basis for the ordering (BE > PE > HDR in iPSCs):
#
#   Base editing (FEASIBLE=0.95, MARGINAL=0.80):
#     BE achieves 30-70% efficiency in iPSCs (Komor et al., Nature, 2016;
#     Gaudelli et al., Nature, 2017). Zero DSBs means no p53 selection.
#     0.95 = normalized to near-maximum because when PAM+window are verified,
#     BE is almost always the preferred choice. MARGINAL (0.80) for cases
#     where the target base is at the edge of the editing window.
#
#   Prime editing (FEASIBLE=0.82, MARGINAL=0.68):
#     PE achieves 5-50% efficiency depending on locus (Anzalone et al.,
#     Nature, 2019; Chen et al., Cell, 2021). Zero DSBs but lower
#     efficiency than BE and more complex delivery (pegRNA + PE protein).
#     0.82 reflects PE's universal mutation-type compatibility offset by
#     its lower and more variable efficiency vs. BE.
#
#   HDR (FEASIBLE=0.72, MARGINAL=0.60):
#     HDR achieves 5-15% in iPSCs (unenhanced), requires a DSB, and
#     has p53 selection concerns. 0.72 reflects the DSB requirement
#     and lower baseline efficiency. MARGINAL (0.60) for guides far
#     from the edit site.
#
# Dual/hybrid strategy priors are derived from single-modality priors
# with multiplicative penalties for additional complexity:
#   Dual BE = 0.92 (two independent BE operations, slightly lower than single)
#   Dual PE = 0.84 (two pegRNAs, moderate complexity)
#   Sequential HDR = 0.68 (two rounds, each with DSB risk)
#   Hybrid BE+HDR = 0.80 (one DSB-free + one DSB operation)
#   Hybrid PE+HDR = 0.76 (similar but PE is less efficient than BE)
#
# NOTE: These priors are SECONDARY to the PAM-verified feasibility scores.
# A strategy with prior 0.95 but no PAM site gets rejected, while a strategy
# with prior 0.72 and a high-quality guide may rank well. The TOPSIS
# sensitivity analysis (10,000 weight permutations) explores how sensitive
# the ranking is to these prior assumptions.
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# Hard capability gates
# ═══════════════════════════════════════════════════════════════════════════
#
# These encode hard biological/technical limits that cannot be overcome by
# favourable PAM, window, or donor geometry. A variant exceeding a modality's
# capability envelope is not a candidate for that modality regardless of
# other scoring factors. This complements — not replaces — the soft priors
# and feasibility results above: the soft priors rank among candidates that
# are capable; the hard gates decide *what is a candidate at all*.
#
# Prime editing (single-pegRNA PE2/PE3):
#   Insertions are bounded by RT-template length (empirically ≤40 bp;
#   Anzalone et al., Nature, 2019). Deletions in the same RT-encoded
#   regime extend to roughly ≤50 bp. Larger structural edits require
#   twin-prime (e.g. PASTE; Anzalone et al., Nat Biotechnol, 2022) — a
#   distinct editor construct, not single-pegRNA PE. The 50-bp threshold
#   below is a conservative single-pegRNA cap.
#
# Base editing (ABE/CBE):
#   Mutation-class (transition-only) and window-position limits are already
#   enforced upstream via FeasibilityLabel.NOT_FEASIBLE, which causes the
#   generator's existing `_is_rankable` guard to skip BE strategies
#   automatically. No additional gate is needed at this layer.
# ═══════════════════════════════════════════════════════════════════════════


PE_MAX_EDIT_SPAN_BP = 50


def _compute_edit_span(variant: NormalizedVariant) -> int:
    """Return the bp span of a variant's edit.

    Prefers the explicit ``GenomicVariantInput.structural_span_bp`` when
    present (for structural variants encoded with placeholder alleles).
    Falls back to ``max(len(ref_allele), len(alt_allele))`` for normal
    alleles. Always returns at least 1.
    """
    gvi = variant.input
    if gvi.structural_span_bp is not None:
        return max(1, int(gvi.structural_span_bp))
    ref_len = 0 if gvi.ref_allele in ("-", "") else len(gvi.ref_allele)
    alt_len = 0 if gvi.alt_allele in ("-", "") else len(gvi.alt_allele)
    return max(ref_len, alt_len, 1)


def _pe_capability_reject(variant: NormalizedVariant) -> Optional[str]:
    """Return a rejection reason string if PE cannot execute this variant.

    A non-None return means the variant's edit span exceeds the single-
    pegRNA PE capability limit. Callers should emit PE as a rejected
    strategy (or skip it in combined strategies).
    """
    span = _compute_edit_span(variant)
    if span > PE_MAX_EDIT_SPAN_BP:
        return (
            f"Edit span {span} bp exceeds the single-pegRNA prime-editing "
            f"capability limit of {PE_MAX_EDIT_SPAN_BP} bp "
            f"(Anzalone et al., Nature, 2019; Chen et al., Cell, 2021). "
            f"Large structural changes require HDR with a donor template "
            f"or twin-prime / PASTE architectures."
        )
    return None


class StrategyGenerator:
    """Generate and score all feasible editing strategies.

    Parameters
    ----------
    p53_active : bool
        Whether the cell type has active p53. Affects safety scoring
        and dual-DSB penalties.
    """

    def __init__(self, p53_active: bool = True):
        self.p53_active = p53_active
        self._integrator = AnnotationIntegrator(p53_active=p53_active)

    def generate(
        self,
        bundles: List[FeasibilityBundle],
    ) -> List[Strategy]:
        """Generate all strategies for the given feasibility bundles.

        For 1 variant: generates single-modality strategies.
        For 2 variants: generates single + dual + hybrid strategies.
        For 3+ variants: generates single strategies for each, plus
        pairwise combinations for adjacent pairs.

        Parameters
        ----------
        bundles : List[FeasibilityBundle]
            One bundle per variant, each containing BE/PE/HDR feasibility.

        Returns
        -------
        List[Strategy]
            All strategies, including rejected ones (check is_rejected).
        """
        strategies: List[Strategy] = []

        if len(bundles) == 1:
            strategies.extend(
                self._generate_single_mutation_strategies(bundles[0])
            )
        elif len(bundles) == 2:
            # Single-mutation strategies for each variant
            for b in bundles:
                strategies.extend(
                    self._generate_single_mutation_strategies(b)
                )
            # Two-mutation combination strategies
            strategies.extend(
                self._generate_two_mutation_strategies(bundles[0], bundles[1])
            )
        else:
            # For 3+ mutations, generate singles and pairwise combos
            for b in bundles:
                strategies.extend(
                    self._generate_single_mutation_strategies(b)
                )
            for i in range(len(bundles) - 1):
                strategies.extend(
                    self._generate_two_mutation_strategies(
                        bundles[i], bundles[i + 1]
                    )
                )

        return self._deduplicate_strategies(strategies)

    # ───────────────────────────────────────────────────────────────────
    # Single-mutation strategies
    # ───────────────────────────────────────────────────────────────────

    def _generate_single_mutation_strategies(
        self,
        bundle: FeasibilityBundle,
    ) -> List[Strategy]:
        """Generate all single-modality strategies for one variant."""
        strategies: List[Strategy] = []

        # --- Base editing ---
        best_be = bundle.best_base_editing_result()
        if best_be and best_be.label != FeasibilityLabel.NOT_FEASIBLE:
            bystander_severity = _extract_bystander_severity(best_be)
            strategies.append(
                self._make_single_edit_strategy(
                    name="Single-step Base Editing",
                    mutation_index=bundle.mutation_index,
                    modality=(
                        EditModality.ABE
                        if best_be.editor_type == "ABE"
                        else EditModality.CBE
                    ),
                    feasibility_results=[best_be],
                    num_dsbs=0,
                    num_rounds=1,
                    num_guides=1,
                    num_proteins=1,
                    num_donors=0,
                    rearrangement_risk=RiskLevel.LOW,
                    evidence_tier=EvidenceTier.A,
                    bystander_severity=bystander_severity,
                    modality_prior_score=(
                        0.95 if best_be.label == FeasibilityLabel.FEASIBLE
                        else 0.80
                    ),
                    included_reasons=[
                        "Single-mutation correction is compatible "
                        "with base editing. BE offers higher efficiency "
                        "(30-70%) and simpler delivery than PE."
                    ],
                    penalties=_warnings_to_penalties(best_be),
                )
            )

        # --- Prime editing ---
        pe = bundle.prime_editing_result
        if pe and pe.label != FeasibilityLabel.NOT_FEASIBLE:
            pe_reject = _pe_capability_reject(bundle.variant)
            if pe_reject is not None:
                strategies.append(
                    Strategy(
                        name="Single-step Prime Editing",
                        steps=[
                            StrategyStep(
                                modality=EditModality.PE,
                                target_mutation_index=bundle.mutation_index,
                                donor_required=False,
                                editor_name=pe.metadata.get(
                                    "editor", "Prime Editor"
                                ),
                            ),
                        ],
                        num_dsbs=0,
                        num_rounds=1,
                        num_distinct_guides=1,
                        num_distinct_proteins=1,
                        num_donors=0,
                        p53_active=self.p53_active,
                        rearrangement_risk=RiskLevel.LOW,
                        evidence_tier=EvidenceTier.A,
                        feasibility_results=[pe],
                        bystander_severity=0.0,
                        donor_feasibility_score=1.0,
                        modality_prior_score=0.0,
                        included_reasons=[],
                        penalties=_warnings_to_penalties(pe),
                        rejection_reasons=[pe_reject],
                    )
                )
            else:
                strategies.append(
                    self._make_single_edit_strategy(
                        name="Single-step Prime Editing",
                        mutation_index=bundle.mutation_index,
                        modality=EditModality.PE,
                        feasibility_results=[pe],
                        num_dsbs=0,
                        num_rounds=1,
                        num_guides=1,
                        num_proteins=1,
                        num_donors=0,
                        rearrangement_risk=RiskLevel.LOW,
                        evidence_tier=EvidenceTier.A,
                        bystander_severity=0.0,
                        modality_prior_score=(
                            0.82 if pe.label == FeasibilityLabel.FEASIBLE
                            else 0.68
                        ),
                        included_reasons=[
                            "Single-mutation correction is compatible "
                            "with prime editing."
                        ],
                        penalties=_warnings_to_penalties(pe),
                    )
                )

        # --- HDR ---
        hdr = bundle.hdr_result
        if hdr and hdr.label != FeasibilityLabel.NOT_FEASIBLE:
            donor_score = _extract_donor_feasibility_score(hdr)
            strategies.append(
                self._make_single_edit_strategy(
                    name="Single-step HDR",
                    mutation_index=bundle.mutation_index,
                    modality=EditModality.HDR_CSSDNA,
                    feasibility_results=[hdr],
                    num_dsbs=1,
                    num_rounds=1,
                    num_guides=1,
                    num_proteins=1,
                    num_donors=1,
                    rearrangement_risk=RiskLevel.LOW,
                    evidence_tier=EvidenceTier.B,
                    bystander_severity=0.0,
                    modality_prior_score=(
                        0.72 if hdr.label == FeasibilityLabel.FEASIBLE
                        else 0.60
                    ),
                    donor_feasibility_score=donor_score,
                    included_reasons=[
                        "Single-mutation correction is compatible "
                        "with HDR design."
                    ],
                    penalties=_warnings_to_penalties(hdr),
                )
            )

        return strategies

    # ───────────────────────────────────────────────────────────────────
    # Two-mutation strategies
    # ───────────────────────────────────────────────────────────────────

    def _generate_two_mutation_strategies(
        self,
        b1: FeasibilityBundle,
        b2: FeasibilityBundle,
    ) -> List[Strategy]:
        """Generate two-mutation combination strategies."""
        strategies: List[Strategy] = []

        be1 = b1.best_base_editing_result()
        be2 = b2.best_base_editing_result()
        pe1 = b1.prime_editing_result
        pe2 = b2.prime_editing_result
        hdr1 = b1.hdr_result
        hdr2 = b2.hdr_result

        # --- Dual base editing ---
        if _is_rankable(be1) and _is_rankable(be2):
            penalties = _warnings_to_penalties(be1) + _warnings_to_penalties(be2)
            bystander = max(
                _extract_bystander_severity(be1),
                _extract_bystander_severity(be2),
            )
            strategies.append(
                Strategy(
                    name="Dual Base Editing",
                    steps=[
                        StrategyStep(
                            modality=(
                                EditModality.ABE
                                if "ABE" in be1.metadata.get("editor", "")
                                else EditModality.CBE
                            ),
                            target_mutation_index=b1.mutation_index,
                            donor_required=False,
                            editor_name=be1.metadata.get("editor", ""),
                        ),
                        StrategyStep(
                            modality=(
                                EditModality.ABE
                                if "ABE" in be2.metadata.get("editor", "")
                                else EditModality.CBE
                            ),
                            target_mutation_index=b2.mutation_index,
                            donor_required=False,
                            editor_name=be2.metadata.get("editor", ""),
                        ),
                    ],
                    num_dsbs=0,
                    simultaneous_dsbs=False,
                    num_rounds=1,
                    num_distinct_guides=2,
                    num_distinct_proteins=_count_distinct_editors([be1, be2]),
                    num_donors=0,
                    requires_selection=False,
                    screening_clones=16,
                    p53_active=self.p53_active,
                    rearrangement_risk=RiskLevel.LOW,
                    evidence_tier=EvidenceTier.A,
                    feasibility_results=[be1, be2],
                    bystander_severity=bystander,
                    donor_feasibility_score=1.0,
                    modality_prior_score=(
                        0.92 if _all_clean([be1, be2]) else 0.78
                    ),
                    estimated_duration_weeks=4,
                    included_reasons=[
                        "Both mutations are addressable by base editing "
                        "without requiring DSBs."
                    ],
                    penalties=penalties,
                    rejection_reasons=[],
                )
            )

        # --- Dual prime editing ---
        pe_gated_b1 = _pe_capability_reject(b1.variant) is not None
        pe_gated_b2 = _pe_capability_reject(b2.variant) is not None
        if (
            _is_rankable(pe1) and _is_rankable(pe2)
            and not pe_gated_b1 and not pe_gated_b2
        ):
            penalties = _warnings_to_penalties(pe1) + _warnings_to_penalties(pe2)
            strategies.append(
                Strategy(
                    name="Dual Prime Editing",
                    steps=[
                        StrategyStep(
                            modality=EditModality.PE,
                            target_mutation_index=b1.mutation_index,
                            donor_required=False,
                            editor_name=pe1.metadata.get("editor", "Prime Editor"),
                        ),
                        StrategyStep(
                            modality=EditModality.PE,
                            target_mutation_index=b2.mutation_index,
                            donor_required=False,
                            editor_name=pe2.metadata.get("editor", "Prime Editor"),
                        ),
                    ],
                    num_dsbs=0,
                    simultaneous_dsbs=False,
                    num_rounds=1,
                    num_distinct_guides=2,
                    num_distinct_proteins=1,
                    num_donors=0,
                    requires_selection=False,
                    screening_clones=20,
                    p53_active=self.p53_active,
                    rearrangement_risk=RiskLevel.LOW,
                    evidence_tier=EvidenceTier.A,
                    feasibility_results=[pe1, pe2],
                    bystander_severity=0.0,
                    donor_feasibility_score=1.0,
                    modality_prior_score=(
                        0.84 if _all_clean([pe1, pe2]) else 0.70
                    ),
                    estimated_duration_weeks=4,
                    included_reasons=[
                        "Both mutations are addressable by prime editing "
                        "without requiring DSBs."
                    ],
                    penalties=penalties,
                    rejection_reasons=[],
                )
            )

        # --- Sequential HDR ---
        if _is_rankable(hdr1) and _is_rankable(hdr2):
            donor_score = min(
                _extract_donor_feasibility_score(hdr1),
                _extract_donor_feasibility_score(hdr2),
            )
            penalties = _warnings_to_penalties(hdr1) + _warnings_to_penalties(hdr2)
            penalties.append("Two editing rounds required.")
            strategies.append(
                Strategy(
                    name="Sequential HDR",
                    steps=[
                        StrategyStep(
                            modality=EditModality.HDR_CSSDNA,
                            target_mutation_index=b1.mutation_index,
                            donor_required=True,
                            editor_name=hdr1.metadata.get("nuclease", ""),
                        ),
                        StrategyStep(
                            modality=EditModality.HDR_CSSDNA,
                            target_mutation_index=b2.mutation_index,
                            donor_required=True,
                            editor_name=hdr2.metadata.get("nuclease", ""),
                        ),
                    ],
                    num_dsbs=1,
                    simultaneous_dsbs=False,
                    num_rounds=2,
                    num_distinct_guides=2,
                    num_distinct_proteins=1,
                    num_donors=2,
                    requires_selection=False,
                    screening_clones=32,
                    p53_active=self.p53_active,
                    rearrangement_risk=RiskLevel.LOW,
                    evidence_tier=EvidenceTier.B,
                    feasibility_results=[hdr1, hdr2],
                    bystander_severity=0.0,
                    donor_feasibility_score=donor_score,
                    modality_prior_score=(
                        0.68 if _all_clean([hdr1, hdr2]) else 0.56
                    ),
                    estimated_duration_weeks=8,
                    included_reasons=[
                        "Both mutations are individually addressable by HDR.",
                        "Sequential design avoids simultaneous DSB burden.",
                    ],
                    penalties=penalties,
                    rejection_reasons=[],
                )
            )

        # --- Hybrid BE + HDR (both directions) ---
        for be_bundle, hdr_bundle in [(b1, b2), (b2, b1)]:
            best_be = be_bundle.best_base_editing_result()
            hdr = hdr_bundle.hdr_result
            if _is_rankable(best_be) and _is_rankable(hdr):
                donor_score = _extract_donor_feasibility_score(hdr)
                penalties = (
                    _warnings_to_penalties(best_be)
                    + _warnings_to_penalties(hdr)
                )
                strategies.append(
                    Strategy(
                        name=f"Hybrid {best_be.metadata.get('editor', 'Base Editing')} + HDR",
                        steps=[
                            StrategyStep(
                                modality=(
                                    EditModality.ABE
                                    if "ABE" in best_be.metadata.get("editor", "")
                                    else EditModality.CBE
                                ),
                                target_mutation_index=be_bundle.mutation_index,
                                donor_required=False,
                                editor_name=best_be.metadata.get("editor", ""),
                            ),
                            StrategyStep(
                                modality=EditModality.HDR_CSSDNA,
                                target_mutation_index=hdr_bundle.mutation_index,
                                donor_required=True,
                                editor_name=hdr.metadata.get("nuclease", ""),
                            ),
                        ],
                        num_dsbs=1,
                        simultaneous_dsbs=False,
                        num_rounds=1,
                        num_distinct_guides=2,
                        num_distinct_proteins=2,
                        num_donors=1,
                        requires_selection=False,
                        screening_clones=20,
                        p53_active=self.p53_active,
                        rearrangement_risk=RiskLevel.LOW,
                        evidence_tier=EvidenceTier.A,
                        feasibility_results=[best_be, hdr],
                        bystander_severity=_extract_bystander_severity(best_be),
                        donor_feasibility_score=donor_score,
                        modality_prior_score=(
                            0.80 if _all_clean([best_be, hdr]) else 0.66
                        ),
                        estimated_duration_weeks=4,
                        included_reasons=[
                            "One mutation is cleanly addressable by base editing.",
                            "The second mutation is addressable by HDR.",
                        ],
                        penalties=penalties,
                        rejection_reasons=[],
                    )
                )

        # --- Hybrid PE + HDR (both directions) ---
        for pe_bundle, hdr_bundle in [(b1, b2), (b2, b1)]:
            pe = pe_bundle.prime_editing_result
            hdr = hdr_bundle.hdr_result
            pe_gated = _pe_capability_reject(pe_bundle.variant) is not None
            if _is_rankable(pe) and _is_rankable(hdr) and not pe_gated:
                donor_score = _extract_donor_feasibility_score(hdr)
                penalties = (
                    _warnings_to_penalties(pe)
                    + _warnings_to_penalties(hdr)
                )
                strategies.append(
                    Strategy(
                        name="Hybrid Prime Editing + HDR",
                        steps=[
                            StrategyStep(
                                modality=EditModality.PE,
                                target_mutation_index=pe_bundle.mutation_index,
                                donor_required=False,
                                editor_name=pe.metadata.get("editor", "Prime Editor"),
                            ),
                            StrategyStep(
                                modality=EditModality.HDR_CSSDNA,
                                target_mutation_index=hdr_bundle.mutation_index,
                                donor_required=True,
                                editor_name=hdr.metadata.get("nuclease", ""),
                            ),
                        ],
                        num_dsbs=1,
                        simultaneous_dsbs=False,
                        num_rounds=1,
                        num_distinct_guides=2,
                        num_distinct_proteins=2,
                        num_donors=1,
                        requires_selection=False,
                        screening_clones=22,
                        p53_active=self.p53_active,
                        rearrangement_risk=RiskLevel.LOW,
                        evidence_tier=EvidenceTier.A,
                        feasibility_results=[pe, hdr],
                        bystander_severity=0.0,
                        donor_feasibility_score=donor_score,
                        modality_prior_score=(
                            0.76 if _all_clean([pe, hdr]) else 0.63
                        ),
                        estimated_duration_weeks=4,
                        included_reasons=[
                            "One mutation is addressable by prime editing.",
                            "The second mutation is addressable by HDR.",
                        ],
                        penalties=penalties,
                        rejection_reasons=[],
                    )
                )

        return self._deduplicate_strategies(strategies)

    # ───────────────────────────────────────────────────────────────────
    # Helper: make a single-edit strategy
    # ───────────────────────────────────────────────────────────────────

    def _make_single_edit_strategy(
        self,
        name: str,
        mutation_index: int,
        modality: EditModality,
        feasibility_results: list,
        num_dsbs: int,
        num_rounds: int,
        num_guides: int,
        num_proteins: int,
        num_donors: int,
        rearrangement_risk: RiskLevel,
        evidence_tier: EvidenceTier,
        bystander_severity: float,
        modality_prior_score: float,
        included_reasons: List[str],
        penalties: List[str],
        donor_feasibility_score: float = 1.0,
    ) -> Strategy:
        return Strategy(
            name=name,
            steps=[
                StrategyStep(
                    modality=modality,
                    target_mutation_index=mutation_index,
                    donor_required=(num_donors > 0),
                    editor_name=(
                        feasibility_results[0].metadata.get("editor", "")
                        if feasibility_results
                        else ""
                    ),
                )
            ],
            num_dsbs=num_dsbs,
            simultaneous_dsbs=False,
            num_rounds=num_rounds,
            num_distinct_guides=num_guides,
            num_distinct_proteins=num_proteins,
            num_donors=num_donors,
            requires_selection=False,
            screening_clones=12 if num_dsbs == 0 else 18,
            p53_active=self.p53_active,
            rearrangement_risk=rearrangement_risk,
            evidence_tier=evidence_tier,
            feasibility_results=feasibility_results,
            bystander_severity=bystander_severity,
            donor_feasibility_score=donor_feasibility_score,
            modality_prior_score=modality_prior_score,
            estimated_duration_weeks=4,
            included_reasons=included_reasons,
            penalties=penalties,
            rejection_reasons=[],
        )

    def _deduplicate_strategies(
        self, strategies: List[Strategy]
    ) -> List[Strategy]:
        """Remove duplicate strategies by (name, step modalities)."""
        seen = set()
        unique: List[Strategy] = []
        for s in strategies:
            key = (
                s.name,
                tuple(
                    sorted(
                        (step.modality.value, step.target_mutation_index)
                        for step in s.steps
                    )
                ),
            )
            if key not in seen:
                seen.add(key)
                unique.append(s)
        return unique


# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════


def _is_rankable(result) -> bool:
    """Check if a feasibility result is usable (not None, not NOT_FEASIBLE)."""
    return result is not None and result.label != FeasibilityLabel.NOT_FEASIBLE


def _all_clean(results: Sequence) -> bool:
    """Check if all results are cleanly FEASIBLE (no warnings)."""
    filtered = [r for r in results if r is not None]
    return bool(filtered) and all(
        r.label == FeasibilityLabel.FEASIBLE for r in filtered
    )


def _extract_bystander_severity(result: BaseEditingFeasibility) -> float:
    """Convert bystander count to a 0-1 severity score."""
    n = result.metadata.get("editable_bystanders", 0)
    return min(1.0, 0.2 * n)


def _extract_donor_feasibility_score(result: HDRFeasibility) -> float:
    """Score donor design quality based on total donor size."""
    design = result.metadata.get("design_summary", {})
    total = design.get("total_donor_size")
    if total is None:
        return result.score
    if total <= 3000:
        return 1.0
    if total <= 6000:
        return 0.8
    if total <= 8000:
        return 0.6
    return 0.4


def _warnings_to_penalties(result) -> List[str]:
    """Extract warnings from a feasibility result as penalty strings."""
    return list(result.warnings)


def _count_distinct_editors(results: Sequence) -> int:
    """Count distinct editor types among feasibility results."""
    editors = set()
    for r in results:
        ed = r.metadata.get("editor")
        if ed:
            editors.add(ed)
    return max(1, len(editors))


# ═══════════════════════════════════════════════════════════════════════════
# Self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Testing generator.py...")

    gen = StrategyGenerator(p53_active=True)

    # Test with empty bundles
    strategies = gen.generate([])
    assert len(strategies) == 0

    print("PASS — generator instantiates and handles empty input")
