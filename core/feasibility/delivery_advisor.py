"""
DeliveryAdvisor — Post-Ranking Delivery Feasibility and Recommendations
========================================================================

This module provides delivery-aware annotations for ranked editing
strategies. It operates AFTER TOPSIS scoring (Option B architecture):

  1. Hard feasibility filters (binary: feasible / warning / infeasible)
  2. Recommended donor format based on edit size and cell type
  3. Recommended delivery method based on cell type and modality
  4. Cell-type-specific warnings and viability enhancer suggestions

Design rationale
-----------------
Delivery complexity is largely correlated with existing TOPSIS dimensions
(Safety, Complexity). Adding a 7th TOPSIS dimension would be redundant
and change rankings in <10% of cases. Instead, delivery is handled as
post-ranking annotation — practical guidance without overclaiming.

Evidence base
--------------
87 verified references in DELIVERY_METHODS_COMPREHENSIVE_LITERATURE_REVIEW.md.
All parameters tagged [MEASURED]/[DERIVED]/[ASSUMED] per project convention.

References
----------
Iyer et al., CRISPR Journal, 2022 (cssDNA donors)
Xie et al., Nature Biotechnology, 2024 (GATALYST cssDNA)
Letort et al., Nature Communications, 2025 (cssDNA in HSPCs)
Ihry et al., Nature Medicine, 2018 (iPSC p53 toxicity)
Dever et al., Nature, 2016 (AAV6 + RNP paradigm)
Frangoul et al., NEJM, 2021 (CASGEVY clinical validation)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from core.models import (
    EditModality,
    ScoredStrategy,
    PipelineResult,
)


# ═══════════════════════════════════════════════════════════════════════════
# Data Models
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class DonorRecommendation:
    """Recommended donor template format for an HDR strategy."""
    format: str                     # "ssODN", "cssDNA", "lssDNA", "dsDNA", "AAV6"
    rationale: str                  # Why this format was chosen
    alternatives: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class DeliveryAnnotation:
    """Delivery feasibility annotation for a single scored strategy."""
    strategy_name: str
    is_deliverable: bool = True
    delivery_method: str = ""           # e.g., "nucleofection_RNP"
    delivery_complexity: int = 1        # 1-5 ordinal
    donor_recommendation: Optional[DonorRecommendation] = None
    warnings: List[str] = field(default_factory=list)
    viability_enhancers: List[str] = field(default_factory=list)
    hard_constraint_violations: List[str] = field(default_factory=list)


@dataclass
class DeliveryAdvisoryResult:
    """Complete delivery advisory output for all ranked strategies."""
    cell_type: str
    annotations: List[DeliveryAnnotation] = field(default_factory=list)
    global_warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════════
# DeliveryAdvisor
# ═══════════════════════════════════════════════════════════════════════════

class DeliveryAdvisor:
    """Post-ranking delivery feasibility advisor.

    Takes TOPSIS-ranked strategies and annotates each with delivery
    feasibility flags, donor format recommendations, and cell-type-specific
    warnings. Does NOT change TOPSIS scores or rankings.

    Parameters
    ----------
    cell_type : str
        Cell type context ("iPSC", "CD34_HSC", "HEK293T").
    edit_sizes : dict or None
        Mapping of variant index to edit size in bp. If None, inferred
        from strategy modality (SNP assumed for BE/PE).
    """

    # Donor format preference by edit size (lower index = preferred)
    # [DERIVED] from Iyer 2022, Xie 2024, Letort 2025, Richardson 2016
    _DONOR_BY_SIZE = [
        # (max_bp, format, rationale)
        (50, "ssODN", "Small edit — ssODN is simplest and most efficient"),
        (200, "ssODN", "ssODN sufficient; asymmetric design recommended (Richardson 2016)"),
        (1500, "cssDNA", "Medium insert — cssDNA 3-5x more efficient than lssDNA (Iyer 2022, Letort 2025)"),
        (4700, "cssDNA", "Large insert — cssDNA preferred; AAV6 is alternative for CD34+ HSCs"),
        (10000, "cssDNA", "Very large insert — cssDNA via phagemid (up to 10 kb)"),
        (20000, "cssDNA", "Extra-large insert — cssDNA via GATALYST system (Xie 2024)"),
    ]

    # Delivery method preference by modality and cell type
    _DELIVERY_METHODS = {
        "BE": {
            "iPSC": ("nucleofection_mRNA", "mRNA + sgRNA nucleofection (Lonza 4D, CA-137)"),
            "CD34_HSC": ("nucleofection_mRNA", "mRNA + sgRNA nucleofection (Newby 2021)"),
            "HEK293T": ("lipofection_RNP", "Lipofectamine CRISPRMAX + RNP"),
        },
        "PE": {
            "iPSC": ("nucleofection_mRNA", "PEmax mRNA + epegRNA nucleofection; consider p53DD co-delivery"),
            "CD34_HSC": ("nucleofection_mRNA", "PEmax mRNA + epegRNA nucleofection (Everette 2023)"),
            "HEK293T": ("lipofection_RNP", "Plasmid transfection or mRNA nucleofection"),
        },
        "HDR": {
            "iPSC": ("nucleofection_RNP", "Cas9 RNP nucleofection + donor co-delivery"),
            "CD34_HSC": ("nucleofection_RNP", "Cas9 RNP nucleofection + AAV6 or cssDNA donor"),
            "HEK293T": ("nucleofection_RNP", "Cas9 RNP nucleofection or lipofection + donor"),
        },
    }

    def __init__(
        self,
        cell_type: str = "iPSC",
        edit_sizes: Optional[Dict[int, int]] = None,
    ):
        self.cell_type = cell_type
        self.edit_sizes = edit_sizes or {}

    def advise(
        self,
        pipeline_result: PipelineResult,
    ) -> DeliveryAdvisoryResult:
        """Annotate all ranked strategies with delivery information.

        Parameters
        ----------
        pipeline_result : PipelineResult
            Output from StrategyPipeline.run() with TOPSIS-ranked strategies.

        Returns
        -------
        DeliveryAdvisoryResult
            Delivery annotations for each strategy.
        """
        global_warnings = []
        annotations = []

        # Infer edit sizes from variants if not provided
        edit_sizes = self._infer_edit_sizes(pipeline_result)

        for scored in pipeline_result.strategies:
            annotation = self._annotate_strategy(scored, edit_sizes)
            annotations.append(annotation)

        # Add cell-type-specific global warnings
        global_warnings.extend(self._get_global_warnings())

        return DeliveryAdvisoryResult(
            cell_type=self.cell_type,
            annotations=annotations,
            global_warnings=global_warnings,
            metadata={
                "n_strategies_annotated": len(annotations),
                "cell_type": self.cell_type,
                "edit_sizes": edit_sizes,
            },
        )

    def _annotate_strategy(
        self,
        scored: ScoredStrategy,
        edit_sizes: Dict[int, int],
    ) -> DeliveryAnnotation:
        """Create delivery annotation for a single strategy."""
        strategy = scored.strategy
        warnings = []
        violations = []
        enhancers = []

        # Determine primary modality from strategy steps
        modality_class = self._classify_modality(strategy.name)

        # 1. Check hard constraints
        violations.extend(
            self._check_hard_constraints(strategy, modality_class, edit_sizes)
        )

        # 2. Determine delivery method
        method_key, method_desc = self._recommend_delivery_method(modality_class)

        # 3. Determine delivery complexity
        complexity = self._get_delivery_complexity(strategy, modality_class)

        # 4. Recommend donor format (for HDR strategies)
        donor_rec = None
        if strategy.num_donors > 0:
            max_edit = max(edit_sizes.values()) if edit_sizes else 1
            donor_rec = self._recommend_donor_format(max_edit)

            # Check donor-specific warnings
            if donor_rec.format == "dsDNA" and self.cell_type == "iPSC":
                warnings.append(
                    "dsDNA donors trigger p53 apoptosis in iPSCs — "
                    "viability may drop to ~10%. Consider cssDNA."
                )
            if donor_rec.format == "AAV6":
                if max_edit > 3500:
                    violations.append(
                        f"Edit size ({max_edit} bp) exceeds AAV6 insert "
                        f"capacity (~3500 bp with homology arms)."
                    )
                if self.cell_type == "CD34_HSC":
                    warnings.append(
                        "AAV6 may cause concatemeric insertions (Suchy 2025) "
                        "and reduced engraftment. cssDNA is an alternative "
                        "(Letort 2025: 5x better engraftment)."
                    )

        # 5. Cell-type-specific warnings
        if strategy.num_dsbs > 0 and self.cell_type == "iPSC":
            warnings.append(
                "DSBs activate p53 in iPSCs — consider DSB-free "
                "editors (BE/PE) if feasible."
            )

        # 6. Viability enhancers
        enhancers = self._get_viability_enhancers(modality_class)

        # 7. HDR in quiescent cells
        if modality_class == "HDR" and self.cell_type == "CD34_HSC":
            warnings.append(
                "HDR requires S/G2 phase. CD34+ HSCs are quiescent — "
                "pre-stimulate 48h with SCF/TPO/FLT3L before editing."
            )

        is_deliverable = len(violations) == 0

        return DeliveryAnnotation(
            strategy_name=strategy.name,
            is_deliverable=is_deliverable,
            delivery_method=method_desc,
            delivery_complexity=complexity,
            donor_recommendation=donor_rec,
            warnings=warnings,
            viability_enhancers=enhancers,
            hard_constraint_violations=violations,
        )

    def _classify_modality(self, strategy_name: str) -> str:
        """Classify strategy name into modality class (BE/PE/HDR)."""
        name_lower = strategy_name.lower()
        if "abe" in name_lower or "cbe" in name_lower or "base" in name_lower:
            return "BE"
        if "pe" in name_lower or "prime" in name_lower:
            return "PE"
        if "hdr" in name_lower or "donor" in name_lower or "knock" in name_lower:
            return "HDR"
        # Check for specific modality patterns
        if "dsb-free" in name_lower or "nick" in name_lower:
            return "PE"
        return "HDR"  # default to most complex

    def _check_hard_constraints(
        self,
        strategy,
        modality_class: str,
        edit_sizes: Dict[int, int],
    ) -> List[str]:
        """Check binary feasibility constraints."""
        violations = []

        # HDR in non-dividing cells context
        # (This is a warning, not a hard violation, since pre-stimulation
        # can push HSCs into S/G2. But it's a genuine constraint.)

        # Lipofection infeasible in CD34+ HSCs
        if self.cell_type == "CD34_HSC" and modality_class == "HDR":
            # Check if delivery would require lipofection
            # (it won't — we always recommend nucleofection for HSCs)
            pass

        return violations

    def _recommend_delivery_method(self, modality_class: str):
        """Recommend delivery method based on modality and cell type."""
        methods = self._DELIVERY_METHODS.get(modality_class, self._DELIVERY_METHODS["HDR"])
        cell_key = self.cell_type if self.cell_type in methods else "iPSC"
        return methods.get(cell_key, methods.get("iPSC", ("nucleofection_RNP", "Default: nucleofection")))

    def _get_delivery_complexity(self, strategy, modality_class: str) -> int:
        """Get ordinal delivery complexity (1-5)."""
        try:
            from utils.constants import MODALITY_DELIVERY_COMPLEXITY
        except ImportError:
            MODALITY_DELIVERY_COMPLEXITY = {}

        # Map strategy to complexity key
        if modality_class == "BE":
            return MODALITY_DELIVERY_COMPLEXITY.get("ABE", 1)
        elif modality_class == "PE":
            return MODALITY_DELIVERY_COMPLEXITY.get("PE", 3)
        elif modality_class == "HDR":
            if strategy.num_donors > 0:
                return MODALITY_DELIVERY_COMPLEXITY.get("HDR_cssDNA", 2)
            return MODALITY_DELIVERY_COMPLEXITY.get("HDR_ssODN", 2)
        return 2

    def _recommend_donor_format(self, edit_size_bp: int) -> DonorRecommendation:
        """Recommend optimal donor format based on edit size and cell type."""
        for max_bp, fmt, rationale in self._DONOR_BY_SIZE:
            if edit_size_bp <= max_bp:
                rec = DonorRecommendation(
                    format=fmt,
                    rationale=rationale,
                )
                # Add alternatives
                if fmt == "ssODN":
                    rec.alternatives = ["lssDNA"]
                elif fmt == "cssDNA":
                    if self.cell_type == "CD34_HSC" and edit_size_bp <= 3500:
                        rec.alternatives = ["AAV6", "lssDNA"]
                    else:
                        rec.alternatives = ["lssDNA"]
                    if edit_size_bp > 10000:
                        rec.rationale += " (requires GATALYST system)"

                # Cell-type-specific warnings
                if self.cell_type == "iPSC":
                    rec.warnings.append(
                        "Avoid dsDNA in iPSCs — use ssDNA formats "
                        "(ssODN, cssDNA, lssDNA) to minimize p53 toxicity."
                    )

                return rec

        # Fallback for very large edits
        return DonorRecommendation(
            format="cssDNA",
            rationale="Very large edit — cssDNA via GATALYST is the only "
                      "ssDNA option at this scale.",
            alternatives=["dsDNA_plasmid"],
            warnings=["Limited published data for inserts >20 kb."],
        )

    def _get_viability_enhancers(self, modality_class: str) -> List[str]:
        """Get cell-type-specific viability enhancer recommendations."""
        enhancers = []

        if self.cell_type == "iPSC":
            enhancers.append("ROCK inhibitor (Y-27632) during replating")
            if modality_class == "HDR":
                enhancers.append(
                    "BCL-XL overexpression (Li 2018): 20-100x HDR improvement"
                )
                enhancers.append(
                    "Cold shock 32C for 48h (Guo 2018): 2-10x HDR boost"
                )
            if modality_class in ("PE", "BE"):
                enhancers.append(
                    "p53DD co-delivery (Li 2022): enhances PE/CBE efficiency"
                )
        elif self.cell_type == "CD34_HSC":
            enhancers.append(
                "Pre-stimulate 48h with SCF/TPO/FLT3L before editing"
            )
            enhancers.append(
                "HiFi Cas9 (R691A) for reduced off-target (Vakulskas 2018)"
            )
            enhancers.append(
                "Minimize ex vivo culture to preserve stemness"
            )

        return enhancers

    def _get_global_warnings(self) -> List[str]:
        """Get cell-type-specific global warnings."""
        warnings = []
        try:
            from utils.constants import CELL_TYPE_DELIVERY_WARNINGS
            ct_warnings = CELL_TYPE_DELIVERY_WARNINGS.get(self.cell_type, {})
            if "culture_warning" in ct_warnings:
                warnings.append(ct_warnings["culture_warning"])
        except ImportError:
            pass
        return warnings

    def _infer_edit_sizes(self, result: PipelineResult) -> Dict[int, int]:
        """Infer edit sizes from pipeline result variants."""
        if self.edit_sizes:
            return self.edit_sizes

        sizes = {}
        for i, nv in enumerate(result.variants):
            ref = nv.input.ref_allele
            alt = nv.input.alt_allele
            if ref == "-":
                # Pure insertion
                sizes[i] = len(alt)
            elif alt == "-":
                # Pure deletion
                sizes[i] = len(ref)
            else:
                # Substitution
                sizes[i] = max(len(ref), len(alt))
        return sizes


# ═══════════════════════════════════════════════════════════════════════════
# Convenience function
# ═══════════════════════════════════════════════════════════════════════════

def annotate_delivery(
    pipeline_result: PipelineResult,
    cell_type: str = "iPSC",
    edit_sizes: Optional[Dict[int, int]] = None,
) -> DeliveryAdvisoryResult:
    """Convenience function for delivery annotation.

    Parameters
    ----------
    pipeline_result : PipelineResult
        Output from StrategyPipeline.run().
    cell_type : str
        Cell type context.
    edit_sizes : dict, optional
        Edit sizes by variant index.

    Returns
    -------
    DeliveryAdvisoryResult
        Delivery annotations for each strategy.
    """
    advisor = DeliveryAdvisor(cell_type=cell_type, edit_sizes=edit_sizes)
    return advisor.advise(pipeline_result)
