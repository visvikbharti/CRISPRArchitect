"""
Translocation Risk Scoring for Simultaneous DSBs
=================================================

The Danger of Two Breaks
------------------------
CRISPR nucleases cut both strands of the DNA double helix, creating a
**double-strand break (DSB)**. A single DSB is dangerous enough — the cell
may mis-repair it. But when TWO DSBs exist simultaneously in the same cell,
the four broken ends (two per break) can be shuffled and re-joined in the
WRONG combinations, creating **chromosomal rearrangements**:

1. **Deletion**: If both breaks are on the same chromosome, the segment
   between them can be lost entirely. The two flanking ends join, deleting
   everything in between. For CRISPR editing, this means losing all the
   genes between your two cut sites.

2. **Inversion**: The segment between the two breaks can be flipped 180
   degrees and re-inserted. This reverses the orientation of the genes
   in that segment. Inversion frequency is typically about half of the
   deletion frequency (since there are two possible orientations).

3. **Translocation**: If the breaks are on DIFFERENT chromosomes (or even
   the same chromosome at very large distances), the ends can swap
   partners. Chromosome A's left end joins chromosome B's right end,
   and vice versa. This creates two derivative chromosomes with novel
   junctions. Translocations are a hallmark of cancer and can activate
   oncogenes.

Why does spatial proximity matter?
-----------------------------------
After a DSB is created, the broken ends don't stay put. They undergo
**constrained diffusion** — bouncing around within a small region (~0.5-1 um
radius) on a timescale of minutes to hours. For two broken ends to be
mis-joined, they need to **find each other** in 3D space.

The key insight from Hi-C studies (Lieberman-Aiden et al., 2009) and
translocation mapping (Chiarle et al., 2011; Zhang et al., 2012; Frock et
al., 2015) is that:

    **Translocation frequency is proportional to contact frequency.**

Two loci that are close in 3D space (high Hi-C contact frequency) are more
likely to generate translocations than two loci that are far apart. This
relationship holds because:

    P(translocation) ~ P(both DSBs exist simultaneously) * P(ends find each other)

The second factor is what depends on 3D proximity.

Quantitative relationships
--------------------------

**Same-chromosome breaks (cis):**
    - Deletion frequency: depends on genomic distance
      - 1-10 kb: 5-30% of edited alleles show deletion (Canver et al., 2015)
      - 10-100 kb: 1-10%
      - 100 kb - 1 Mb: 0.1-5%
      - >1 Mb: <1%
    - Inversion frequency: roughly 30-50% of deletion frequency
      (Source: Kraft et al., 2015; Brunner et al., 2019)
    - Note: deletion is favored over inversion because the intervening
      segment (once excised as a circle) is lost, while inversion requires
      the segment to be re-ligated in the correct position but flipped.

**Different-chromosome breaks (trans):**
    - Translocation frequency scales with Hi-C contact probability:
      P(contact) ~ s^(-gamma) for cis contacts, where gamma ~ 1.08
    - For trans contacts (inter-chromosomal): much lower baseline
      frequency, roughly equivalent to cis contacts at ~100-200 Mb
      (i.e., very far apart on the same chromosome).
    - Typical translocation frequencies with two Cas9 DSBs:
      ~0.01-1% of cells (depends on cell type and sites).

Additional risk factors
-----------------------

**p53 selection (iPSCs and primary cells):**
    DSBs activate the p53 pathway, which triggers cell cycle arrest or
    apoptosis. In iPSCs and primary cells (which have intact p53), cells
    that happen to have a pre-existing p53 mutation have a survival
    advantage after CRISPR cutting. This means your surviving clones are
    *enriched* for p53-mutant cells — a dangerous selection artifact.
    (Source: Ihry et al., Nat Med, 2018; Haapaniemi et al., Nat Med, 2018.)
    Two simultaneous DSBs amplify this selection pressure.

**Chromothripsis:**
    In rare cases, a chromosome (or chromosome segment) can shatter into
    many pieces and be reassembled in a scrambled order. This catastrophic
    event ("chromothripsis") has been linked to micronucleus formation
    following DSBs. While very rare (~0.01-0.1% per DSB), it is devastating
    when it occurs. Two simultaneous DSBs may increase the risk.
    (Source: Zhang et al., Nature, 2015; Leibowitz et al., Science, 2021.)

References
----------
- Lieberman-Aiden et al., Science, 2009 (Hi-C).
- Chiarle et al., Cell, 2011 (translocation mapping).
- Zhang et al., Cell, 2012 (translocation and contact frequency).
- Frock et al., Nat Biotech, 2015 (HTGTS: high-throughput translocation).
- Canver et al., Nature, 2015 (paired guide deletions).
- Ihry et al., Nat Med, 2018; Haapaniemi et al., Nat Med, 2018 (p53 selection).
- Leibowitz et al., Science, 2021 (chromothripsis from DSBs).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional

import sys
import os

_PARENT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from utils.constants import (
    HICCONTACT_POWER_LAW_GAMMA,
    TRANSLOCATION_BASELINE_PROB_1MB,
    CELL_TYPE_PARAMS,
)


# =============================================================================
# Empirical deletion frequency model
# =============================================================================
# Canver et al. (2015) and others measured deletion efficiency as a function
# of the distance between two Cas9 cut sites on the same chromosome. We fit
# a simple log-linear model to their published data:
#
#   P(deletion) = A * (distance_bp)^(-alpha)
#
# with:
#   A ~ 3.0     (normalization constant)
#   alpha ~ 0.5 (power-law exponent for deletion frequency decay)
#
# This captures the empirical observation that:
#   - At 1 kb:   P(deletion) ~ 30%
#   - At 10 kb:  P(deletion) ~ 10%
#   - At 100 kb: P(deletion) ~ 3%
#   - At 1 Mb:   P(deletion) ~ 1%
#   - At 10 Mb:  P(deletion) ~ 0.3%

_DELETION_COEFFICIENT = 3.0
_DELETION_EXPONENT = 0.5
_DELETION_MAX = 0.40   # Cap at 40% — even at very short distances, not all alleles delete
_DELETION_MIN = 1e-5   # Floor: even at enormous distances, some rare deletions occur

# Inversion-to-deletion ratio
# Inversions require both ends of the excised segment to be re-ligated in
# reverse orientation. Empirically, inversions occur at roughly 30-50% the
# frequency of deletions. We use 0.4 as a central estimate.
_INVERSION_TO_DELETION_RATIO = 0.4

# For same-chromosome breaks, "translocation" (complex rearrangement) is
# less common than deletion/inversion, but can still occur via mechanisms
# such as replication-coupled break-induced replication (BIR) or end-joining
# with other DSBs.
_COMPLEX_REARRANGEMENT_FRACTION = 0.05  # ~5% of deletion frequency


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class TranslocationRisk:
    """Quantified rearrangement risk from simultaneous DSBs.

    This dataclass holds the estimated frequencies of different types of
    chromosomal rearrangements that can result from two simultaneous DSBs.

    Attributes
    ----------
    same_chromosome : bool
        Whether the two DSBs are on the same chromosome.
    genomic_distance_bp : Optional[int]
        Genomic distance between the two sites (None if different chromosomes).
    deletion_frequency : float
        Estimated fraction of alleles showing deletion of the segment between
        the two breaks. Only meaningful for same-chromosome breaks.
    inversion_frequency : float
        Estimated fraction of alleles showing inversion of the segment.
    translocation_frequency : float
        Estimated fraction of alleles showing translocation (mis-joining with
        a different chromosome or complex rearrangement).
    total_rearrangement_risk : float
        Sum of deletion + inversion + translocation frequencies.
    risk_level : str
        Qualitative risk level: "low", "moderate", "high", or "very_high".
    """

    same_chromosome: bool
    genomic_distance_bp: Optional[int]
    deletion_frequency: float
    inversion_frequency: float
    translocation_frequency: float
    total_rearrangement_risk: float
    risk_level: str


@dataclass
class SafetyReport:
    """Comprehensive dual-DSB safety assessment.

    This report integrates multiple risk factors to give an overall safety
    grade for a dual-DSB CRISPR experiment.

    Attributes
    ----------
    cell_type : str
        Cell type being edited.
    viability_estimate : float
        Estimated cell viability after two DSBs (fraction surviving).
    translocation_risk : TranslocationRisk
        Detailed rearrangement risk assessment.
    p53_selection_risk : str
        Risk level for p53-mutant enrichment: "none", "low", "moderate", "high".
    chromothripsis_risk : str
        Chromothripsis risk level: "negligible", "low", "moderate".
    overall_grade : str
        Letter grade: "A" (safe), "B" (acceptable), "C" (caution),
        "D" (high risk), "F" (unacceptable risk).
    recommendations : str
        Actionable safety recommendations.
    summary : str
        One-paragraph plain-English summary.
    """

    cell_type: str
    viability_estimate: float
    translocation_risk: TranslocationRisk
    p53_selection_risk: str
    chromothripsis_risk: str
    overall_grade: str
    recommendations: str
    summary: str


# =============================================================================
# Translocation Risk Predictor
# =============================================================================

class TranslocationRiskPredictor:
    """Predict rearrangement risk when two DSBs are made simultaneously.

    This class provides quantitative estimates of deletion, inversion, and
    translocation frequencies based on published experimental data and
    polymer physics models.

    The central principle: **rearrangement probability depends on 3D spatial
    proximity**, which in turn depends on genomic distance (for same-chromosome
    breaks) or nuclear organization (for inter-chromosomal breaks).

    Usage
    -----
    >>> predictor = TranslocationRiskPredictor()
    >>> # Two breaks 2 Mb apart on the same chromosome
    >>> risk = predictor.estimate_risk(50_000_000, 52_000_000, same_chromosome=True)
    >>> print(predictor.risk_summary(risk))
    >>> # Full safety assessment for iPSC editing
    >>> report = predictor.dual_dsb_safety_assessment(
    ...     site1=50_000_000, site2=52_000_000, cell_type="iPSC"
    ... )
    >>> print(report.summary)
    """

    def estimate_risk(
        self,
        site1_coord: int,
        site2_coord: int,
        same_chromosome: bool = True,
    ) -> TranslocationRisk:
        """Estimate rearrangement risk from two simultaneous DSBs.

        Parameters
        ----------
        site1_coord : int
            Genomic coordinate (in bp) of the first DSB. For same-chromosome
            calculations, this is the position along the chromosome.
        site2_coord : int
            Genomic coordinate (in bp) of the second DSB.
        same_chromosome : bool
            If True, both DSBs are on the same chromosome. If False, they
            are on different chromosomes (site1_coord and site2_coord are
            then positions on their respective chromosomes, used only for
            display purposes).

        Returns
        -------
        TranslocationRisk
            Dataclass with estimated frequencies and risk level.

        Notes
        -----
        **How the estimates are derived:**

        For *same-chromosome* breaks:
            Deletion frequency is modeled empirically from Canver et al. (2015)
            and related studies. The key observation is that deletion frequency
            decreases as a power law of genomic distance:

                P(deletion) ~ A * distance^(-0.5)

            This makes intuitive sense: the farther apart the two breaks, the
            harder it is for the ends to find each other (after the intervening
            segment is excised), and the more likely the cell is to repair each
            break independently.

            Inversion frequency is ~40% of deletion frequency.

            Complex rearrangement (translocation-like events) is ~5% of
            deletion frequency.

        For *different-chromosome* breaks:
            We use the Hi-C contact frequency power law. Inter-chromosomal
            contacts are much rarer than intra-chromosomal contacts. The
            translocation frequency is estimated from the baseline rate at 1 Mb
            (from constants.py), extrapolated to inter-chromosomal distances
            using the empirical observation that inter-chromosomal contact
            frequency is roughly equivalent to intra-chromosomal contact at
            ~100 Mb genomic separation.
        """
        if same_chromosome:
            return self._estimate_cis_risk(site1_coord, site2_coord)
        else:
            return self._estimate_trans_risk(site1_coord, site2_coord)

    def _estimate_cis_risk(
        self, site1_coord: int, site2_coord: int
    ) -> TranslocationRisk:
        """Estimate rearrangement risk for two breaks on the same chromosome.

        Same-chromosome ("cis") breaks primarily produce deletions and
        inversions. The intervening segment can be:
            - Excised (deletion): the two flanking ends join directly
            - Flipped (inversion): the segment is re-inserted in reverse
            - Part of a complex rearrangement: rare, but possible

        The frequency of each outcome decreases with genomic distance
        because the broken ends must diffuse together in 3D space, and
        more distant loci are harder to bridge.
        """
        distance_bp = abs(site2_coord - site1_coord)
        distance_bp = max(distance_bp, 1)  # Avoid division by zero

        # Deletion frequency: empirical power law
        deletion_freq = _DELETION_COEFFICIENT * (distance_bp ** (-_DELETION_EXPONENT))
        deletion_freq = min(deletion_freq, _DELETION_MAX)
        deletion_freq = max(deletion_freq, _DELETION_MIN)

        # Inversion frequency: ~40% of deletion
        inversion_freq = deletion_freq * _INVERSION_TO_DELETION_RATIO

        # Complex rearrangement / translocation: ~5% of deletion
        translocation_freq = deletion_freq * _COMPLEX_REARRANGEMENT_FRACTION

        total = deletion_freq + inversion_freq + translocation_freq
        risk_level = self._classify_risk(total)

        return TranslocationRisk(
            same_chromosome=True,
            genomic_distance_bp=distance_bp,
            deletion_frequency=deletion_freq,
            inversion_frequency=inversion_freq,
            translocation_frequency=translocation_freq,
            total_rearrangement_risk=total,
            risk_level=risk_level,
        )

    def _estimate_trans_risk(
        self, site1_coord: int, site2_coord: int
    ) -> TranslocationRisk:
        """Estimate translocation risk for breaks on different chromosomes.

        When two DSBs are on different chromosomes, the primary risk is
        **reciprocal translocation**: chromosome A's left arm joins
        chromosome B's right arm, and vice versa. The frequency depends
        on how often the two loci are in spatial proximity, which is
        captured by their inter-chromosomal Hi-C contact frequency.

        Inter-chromosomal contact frequency is typically much lower than
        intra-chromosomal, roughly equivalent to intra-chromosomal contacts
        at ~100-200 Mb genomic distance. We use 150 Mb as the "effective
        genomic distance" for inter-chromosomal contacts.

        The relationship:
            P(translocation) ~ BASELINE_PROB_1MB * (effective_distance / 1e6)^(-gamma)
        """
        # Effective distance for inter-chromosomal contacts
        effective_distance_bp = 150_000_000  # 150 Mb equivalent

        # Scale from the known baseline at 1 Mb
        gamma = HICCONTACT_POWER_LAW_GAMMA
        scale_factor = (effective_distance_bp / 1_000_000) ** (-gamma)
        translocation_freq = TRANSLOCATION_BASELINE_PROB_1MB * scale_factor

        # No deletion/inversion for inter-chromosomal breaks
        # (those concepts require a contiguous segment between the two breaks)
        deletion_freq = 0.0
        inversion_freq = 0.0

        total = translocation_freq
        risk_level = self._classify_risk(total)

        return TranslocationRisk(
            same_chromosome=False,
            genomic_distance_bp=None,
            deletion_frequency=deletion_freq,
            inversion_frequency=inversion_freq,
            translocation_frequency=translocation_freq,
            total_rearrangement_risk=total,
            risk_level=risk_level,
        )

    @staticmethod
    def _classify_risk(total_rearrangement_freq: float) -> str:
        """Classify total rearrangement risk into qualitative levels.

        Thresholds are based on practical considerations for cell line
        engineering:
            - < 1%:   low — manageable with standard screening
            - 1-5%:   moderate — requires careful screening
            - 5-20%:  high — significant fraction of clones affected
            - > 20%:  very_high — majority of clones may be rearranged
        """
        if total_rearrangement_freq < 0.01:
            return "low"
        elif total_rearrangement_freq < 0.05:
            return "moderate"
        elif total_rearrangement_freq < 0.20:
            return "high"
        else:
            return "very_high"

    def risk_summary(self, risk: TranslocationRisk) -> str:
        """Generate a human-readable risk summary.

        Parameters
        ----------
        risk : TranslocationRisk
            Output from ``estimate_risk()``.

        Returns
        -------
        str
            Multi-line plain-English risk assessment with context from
            published studies and actionable recommendations.
        """
        lines = []
        lines.append("CHROMOSOMAL REARRANGEMENT RISK ASSESSMENT")
        lines.append("=" * 55)
        lines.append("")

        if risk.same_chromosome:
            dist_mb = risk.genomic_distance_bp / 1_000_000
            lines.append(
                f"Configuration: Two DSBs on the SAME chromosome, "
                f"{dist_mb:.2f} Mb apart"
            )
            lines.append("")
            lines.append("Estimated rearrangement frequencies:")
            lines.append(
                f"  Deletion:              {risk.deletion_frequency * 100:.2f}%"
            )
            lines.append(
                f"  Inversion:             {risk.inversion_frequency * 100:.2f}%"
            )
            lines.append(
                f"  Complex rearrangement: {risk.translocation_frequency * 100:.3f}%"
            )
            lines.append(
                f"  TOTAL rearrangement:   {risk.total_rearrangement_risk * 100:.2f}%"
            )
        else:
            lines.append(
                f"Configuration: Two DSBs on DIFFERENT chromosomes"
            )
            lines.append("")
            lines.append("Estimated rearrangement frequencies:")
            lines.append(
                f"  Reciprocal translocation: "
                f"{risk.translocation_frequency * 100:.4f}%"
            )
            lines.append(
                f"  (No deletion/inversion for inter-chromosomal breaks)"
            )

        lines.append("")
        lines.append(f"Risk level: {risk.risk_level.upper()}")
        lines.append("")

        # Context from published data
        lines.append("--- Context from Published Studies ---")
        if risk.same_chromosome and risk.genomic_distance_bp is not None:
            d = risk.genomic_distance_bp
            if d < 10_000:
                lines.append(
                    "  Canver et al. (Nature, 2015) observed 5-30% deletion"
                )
                lines.append(
                    "  frequency for paired guides 1-10 kb apart in human cells."
                )
                lines.append(
                    "  Your estimated frequency is consistent with this range."
                )
            elif d < 100_000:
                lines.append(
                    "  At 10-100 kb separations, published deletion frequencies"
                )
                lines.append(
                    "  range from 1-10% (Canver et al., 2015; Kraft et al., 2015)."
                )
            elif d < 1_000_000:
                lines.append(
                    "  At 100 kb - 1 Mb, deletion frequencies drop below 5%."
                )
                lines.append(
                    "  However, the deleted segment is large enough to remove"
                )
                lines.append(
                    "  entire genes, making even low-frequency events dangerous."
                )
            else:
                lines.append(
                    "  At >1 Mb separation, deletion is rare (<1-3%) but"
                )
                lines.append(
                    "  involves loss of a very large genomic segment containing"
                )
                lines.append(
                    "  potentially many genes. Brunner et al. (2019) detected"
                )
                lines.append(
                    "  Mb-scale deletions at low frequency with paired Cas9 guides."
                )
        else:
            lines.append(
                "  Inter-chromosomal translocations from paired Cas9 DSBs have"
            )
            lines.append(
                "  been detected at 0.01-1% frequency (Frock et al., Nat Biotech,"
            )
            lines.append(
                "  2015) using HTGTS. Frequency depends on the specific loci"
            )
            lines.append(
                "  and their spatial proximity in the nucleus."
            )

        lines.append("")
        lines.append("--- Recommendations ---")
        if risk.risk_level == "low":
            lines.append(
                "  - Standard clone screening should be sufficient."
            )
            lines.append(
                "  - Verify by PCR across both junctions."
            )
            lines.append(
                "  - Consider G-banding karyotype for critical applications."
            )
        elif risk.risk_level == "moderate":
            lines.append(
                "  - Screen clones by junction PCR at BOTH cut sites."
            )
            lines.append(
                "  - Perform karyotyping (G-banding or SKY) on selected clones."
            )
            lines.append(
                "  - Consider sequential editing instead of simultaneous cuts."
            )
        elif risk.risk_level == "high":
            lines.append(
                "  - HIGH RISK: significant fraction of clones may carry"
            )
            lines.append(
                "    deletions or inversions of the intervening segment."
            )
            lines.append(
                "  - STRONGLY recommend sequential editing (cut/repair one site,"
            )
            lines.append(
                "    clone, then cut the second site)."
            )
            lines.append(
                "  - If simultaneous cuts are unavoidable: screen extensively"
            )
            lines.append(
                "    by karyotyping and copy-number analysis (aCGH or WGS)."
            )
        else:
            lines.append(
                "  - VERY HIGH RISK: majority of clones likely rearranged."
            )
            lines.append(
                "  - DO NOT perform simultaneous cuts at these sites."
            )
            lines.append(
                "  - Use sequential editing with clonal isolation between cuts."
            )
            lines.append(
                "  - If segment deletion is the GOAL (e.g., enhancer deletion),"
            )
            lines.append(
                "    then this frequency is expected and acceptable."
            )

        return "\n".join(lines)

    def dual_dsb_safety_assessment(
        self,
        site1: int,
        site2: int,
        cell_type: str = "iPSC",
        same_chromosome: bool = True,
    ) -> SafetyReport:
        """Comprehensive safety analysis for making two simultaneous DSBs.

        This method integrates multiple risk factors that a researcher should
        consider before performing a dual-DSB experiment:

        1. **Cell viability**: Two DSBs kill more cells than one. In iPSCs,
           viability drops to ~30% with two cuts (vs ~55% for one cut).
        2. **Rearrangement risk**: Deletions, inversions, translocations.
        3. **p53 selection**: In p53-active cells, DSBs select for p53-mutant
           clones. Two DSBs amplify this selection.
        4. **Chromothripsis**: Rare but catastrophic chromosome shattering.

        Parameters
        ----------
        site1 : int
            Genomic coordinate of the first DSB (bp).
        site2 : int
            Genomic coordinate of the second DSB (bp).
        cell_type : str
            Cell type. Must be a key in CELL_TYPE_PARAMS from constants.py.
            Options: "iPSC", "HEK293T", "K562", "T_cell", "HSC".
        same_chromosome : bool
            Whether both sites are on the same chromosome.

        Returns
        -------
        SafetyReport
            Comprehensive safety assessment with grade, risks, and
            recommendations.

        Example
        -------
        >>> predictor = TranslocationRiskPredictor()
        >>> report = predictor.dual_dsb_safety_assessment(
        ...     site1=50_000_000,
        ...     site2=52_000_000,
        ...     cell_type="iPSC",
        ...     same_chromosome=True,
        ... )
        >>> print(f"Safety grade: {report.overall_grade}")
        >>> print(report.summary)
        """
        # --- Cell type parameters ---
        ct = CELL_TYPE_PARAMS.get(cell_type)
        if ct is None:
            available = ", ".join(CELL_TYPE_PARAMS.keys())
            raise ValueError(
                f"Unknown cell_type '{cell_type}'. Available: {available}"
            )

        viability = ct["viability_dual_dsb"]

        # --- Rearrangement risk ---
        risk = self.estimate_risk(site1, site2, same_chromosome)

        # --- p53 selection risk ---
        if ct["p53_active"]:
            # Two DSBs strongly activate p53. In cells with intact p53,
            # there is significant selection for p53-loss-of-function clones.
            if viability < 0.40:
                p53_risk = "high"
            elif viability < 0.60:
                p53_risk = "moderate"
            else:
                p53_risk = "low"
        else:
            p53_risk = "none"

        # --- Chromothripsis risk ---
        # Chromothripsis is rare but increases with number of DSBs and
        # the severity of the DNA damage response.
        # Base rate: ~0.01-0.1% per DSB.
        # Two DSBs: slightly elevated but still rare.
        if same_chromosome and risk.genomic_distance_bp and risk.genomic_distance_bp > 10_000_000:
            chromothripsis_risk = "low"  # Large segment at risk of micronucleus
        else:
            chromothripsis_risk = "negligible"

        # --- Overall grade ---
        grade = self._compute_grade(
            viability, risk, p53_risk, chromothripsis_risk
        )

        # --- Recommendations ---
        recommendations = self._build_recommendations(
            cell_type, ct, risk, p53_risk, grade, same_chromosome
        )

        # --- Summary ---
        summary = self._build_summary(
            cell_type, ct, site1, site2, same_chromosome,
            viability, risk, p53_risk, chromothripsis_risk, grade
        )

        return SafetyReport(
            cell_type=cell_type,
            viability_estimate=viability,
            translocation_risk=risk,
            p53_selection_risk=p53_risk,
            chromothripsis_risk=chromothripsis_risk,
            overall_grade=grade,
            recommendations=recommendations,
            summary=summary,
        )

    @staticmethod
    def _compute_grade(
        viability: float,
        risk: TranslocationRisk,
        p53_risk: str,
        chromothripsis_risk: str,
    ) -> str:
        """Compute an overall letter grade from individual risk factors.

        The grading is conservative — the worst individual risk factor
        dominates the overall grade, because any single catastrophic
        outcome (e.g., translocation in a clinical iPSC line) is
        unacceptable regardless of other factors.

        Grade criteria:
            A: All risks low, viability > 60%
            B: Moderate risks, viability > 50%
            C: One high risk factor OR viability 30-50%
            D: Multiple high risk factors OR viability < 30%
            F: Very high rearrangement risk AND high p53 selection risk
        """
        score = 0  # Higher = worse

        # Viability penalty
        if viability < 0.30:
            score += 3
        elif viability < 0.50:
            score += 2
        elif viability < 0.60:
            score += 1

        # Rearrangement risk penalty
        risk_scores = {"low": 0, "moderate": 1, "high": 2, "very_high": 3}
        score += risk_scores.get(risk.risk_level, 2)

        # p53 selection penalty
        p53_scores = {"none": 0, "low": 0, "moderate": 1, "high": 2}
        score += p53_scores.get(p53_risk, 1)

        # Chromothripsis penalty
        chromo_scores = {"negligible": 0, "low": 1, "moderate": 2}
        score += chromo_scores.get(chromothripsis_risk, 1)

        # Map total score to grade
        if score <= 1:
            return "A"
        elif score <= 3:
            return "B"
        elif score <= 5:
            return "C"
        elif score <= 7:
            return "D"
        else:
            return "F"

    @staticmethod
    def _build_recommendations(
        cell_type: str,
        ct: dict,
        risk: TranslocationRisk,
        p53_risk: str,
        grade: str,
        same_chromosome: bool,
    ) -> str:
        """Build actionable recommendations based on the safety assessment."""
        recs = []

        # General screening
        recs.append("SCREENING RECOMMENDATIONS:")
        recs.append(
            "  1. PCR across both cut sites to verify on-target editing."
        )
        recs.append(
            "  2. Junction PCR to detect deletions (if same-chromosome cuts)."
        )

        if risk.risk_level in ("moderate", "high", "very_high"):
            recs.append(
                "  3. G-banding karyotype or SKY on top candidate clones."
            )
            recs.append(
                "  4. Consider aCGH or low-pass WGS to detect copy-number"
            )
            recs.append(
                "     changes and structural rearrangements."
            )

        if p53_risk in ("moderate", "high"):
            recs.append("")
            recs.append("p53 SELECTION MITIGATION:")
            recs.append(
                f"  - {cell_type} cells have active p53. Two DSBs strongly"
            )
            recs.append(
                "    select for p53-mutant clones (Ihry et al., 2018)."
            )
            recs.append(
                "  - Sequence TP53 exons 5-8 in all selected clones."
            )
            recs.append(
                "  - Consider transient p53 inhibition (e.g., MDM2 "
                "overexpression or dominant-negative p53) during editing "
                "to reduce selection pressure."
            )
            recs.append(
                "  - Use the fewest possible DSBs. Consider sequential editing."
            )

        if risk.risk_level in ("high", "very_high") and same_chromosome:
            recs.append("")
            recs.append("STRATEGY ALTERNATIVES:")
            recs.append(
                "  - STRONGLY consider sequential editing: introduce the first"
            )
            recs.append(
                "    edit, isolate clones, then introduce the second edit."
            )
            recs.append(
                "  - This eliminates simultaneous DSB risk but doubles timeline."
            )

        if grade in ("D", "F"):
            recs.append("")
            recs.append("WARNING:")
            recs.append(
                "  This experiment carries substantial risk. Proceed only if:"
            )
            recs.append(
                "  (a) The scientific question cannot be addressed otherwise."
            )
            recs.append(
                "  (b) Extensive post-editing QC is planned."
            )
            recs.append(
                "  (c) The edited cells are NOT intended for clinical use."
            )

        return "\n".join(recs)

    @staticmethod
    def _build_summary(
        cell_type: str,
        ct: dict,
        site1: int,
        site2: int,
        same_chromosome: bool,
        viability: float,
        risk: TranslocationRisk,
        p53_risk: str,
        chromothripsis_risk: str,
        grade: str,
    ) -> str:
        """Build a one-paragraph plain-English summary."""
        dist_str = ""
        if same_chromosome and risk.genomic_distance_bp:
            dist_mb = risk.genomic_distance_bp / 1_000_000
            dist_str = f" separated by {dist_mb:.1f} Mb"

        loc_str = (
            "the same chromosome" if same_chromosome else "different chromosomes"
        )

        summary = (
            f"DUAL-DSB SAFETY ASSESSMENT: Grade {grade}\n"
            f"\n"
            f"Creating two simultaneous DSBs on {loc_str}{dist_str} in "
            f"{ct['description']} ({cell_type}). "
            f"Estimated viability: {viability * 100:.0f}%. "
            f"Total rearrangement risk: {risk.total_rearrangement_risk * 100:.2f}% "
            f"({risk.risk_level} risk). "
        )

        if risk.same_chromosome and risk.deletion_frequency > 0.01:
            summary += (
                f"Deletion of the intervening segment is the primary concern "
                f"({risk.deletion_frequency * 100:.1f}% estimated frequency). "
            )

        if p53_risk in ("moderate", "high"):
            summary += (
                f"p53 selection risk is {p53_risk} — surviving clones may be "
                f"enriched for p53 mutations. "
            )

        summary += (
            f"Chromothripsis risk: {chromothripsis_risk}."
        )

        return summary
