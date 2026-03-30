"""
CRISPRArchitect v3 — Strategy Analysis Page
=============================================

Streamlit page for the v3 multi-nuclease, TOPSIS-scored strategy pipeline.
Features: 5 nucleases, 9 editor profiles, 6D TOPSIS ranking, Monte Carlo
sensitivity analysis, Pareto front, cross-method validation.

Compatible with Python 3.9+.
"""

from __future__ import annotations

import streamlit as st
from typing import List, Optional

# ---------------------------------------------------------------------------
# Try importing v3 modules
# ---------------------------------------------------------------------------
_V3_AVAILABLE = False
_V3_IMPORT_ERROR = ""

try:
    from core.pipeline.strategy_stage import StrategyPipeline
    from core.models import (
        GenomicVariantInput,
        PipelineResult,
        FeasibilityLabel,
        NormalizedVariant,
        FeasibilityBundle,
        ScoredStrategy,
    )
    _V3_AVAILABLE = True
except ImportError as e:
    _V3_IMPORT_ERROR = str(e)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CELL_TYPE_OPTIONS = ["iPSC", "CD34_HSC", "HEK293T", "K562"]
NUCLEASE_OPTIONS = ["SpCas9", "enFnCas9", "SpCas9-NG", "SpRY", "Cas12a"]
NUCLEASE_PAMS = {
    "SpCas9": "NGG",
    "enFnCas9": "NRG",
    "SpCas9-NG": "NG",
    "SpRY": "NNN (near-PAMless)",
    "Cas12a": "TTTV",
}
CHROMOSOME_OPTIONS = [str(i) for i in range(1, 23)] + ["X", "Y"]


# ---------------------------------------------------------------------------
# Styling helpers
# ---------------------------------------------------------------------------
def _card(title: str, content_html: str) -> None:
    st.markdown(
        f'<div class="ca-card"><h3>{title}</h3>{content_html}</div>',
        unsafe_allow_html=True,
    )


def _status_dot(ok: bool) -> str:
    color = "#4CAF50" if ok else "#F44336"
    return f'<span style="color:{color}; font-size:1.2em;">&#9679;</span>'


def _feasibility_badge(label: FeasibilityLabel) -> str:
    color_map = {
        FeasibilityLabel.FEASIBLE: "#4CAF50",
        FeasibilityLabel.MARGINAL: "#FF9800",
        FeasibilityLabel.NOT_FEASIBLE: "#F44336",
    }
    color = color_map.get(label, "#757575")
    text = label.value.replace("_", " ").title()
    return (
        f'<span style="background:{color}; color:white; padding:2px 10px; '
        f'border-radius:12px; font-size:0.85em; font-weight:600;">{text}</span>'
    )


def _tier_badge(tier: str) -> str:
    colors = {"A": "#4CAF50", "B": "#FF9800", "C": "#F44336"}
    color = colors.get(tier, "#757575")
    return (
        f'<span style="background:{color}; color:white; padding:1px 8px; '
        f'border-radius:10px; font-size:0.8em; font-weight:700;">Tier {tier}</span>'
    )


def _stability_bar(stability: float) -> str:
    """Horizontal bar showing rank stability percentage."""
    pct = stability * 100
    if pct >= 80:
        color = "#4CAF50"
    elif pct >= 50:
        color = "#FF9800"
    else:
        color = "#F44336"
    return (
        f'<div style="display:flex; align-items:center; gap:8px;">'
        f'<div style="flex:1; background:#333; border-radius:4px; height:8px; '
        f'max-width:120px;">'
        f'<div style="width:{pct:.0f}%; background:{color}; height:100%; '
        f'border-radius:4px;"></div></div>'
        f'<span style="font-size:0.85em; font-weight:600; color:{color};">'
        f'{pct:.1f}%</span></div>'
    )


def _score_pill(label: str, value: float, is_cost: bool = False) -> str:
    """Small pill showing a dimension score."""
    if is_cost:
        color = "#F44336" if value > 0.5 else ("#FF9800" if value > 0.2 else "#4CAF50")
    else:
        color = "#4CAF50" if value > 0.7 else ("#FF9800" if value > 0.4 else "#F44336")
    return (
        f'<span style="font-size:0.78em; margin-right:6px;">'
        f'{label}: <strong style="color:{color};">{value:.2f}</strong></span>'
    )


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
def render_v3_sidebar() -> Optional[dict]:
    """Render v3 sidebar inputs. Returns params dict on Analyze click."""
    st.sidebar.divider()
    st.sidebar.markdown("### v3: Strategy Analysis")

    input_mode = st.sidebar.radio(
        "Input format",
        ["Genomic coordinates", "HGVS notation"],
        key="v3_input_mode",
        horizontal=True,
    )

    gene = st.sidebar.text_input(
        "Gene symbol", value="NF1",
        help="HGNC gene symbol (e.g., NF1, BRCA1, TP53)",
        key="v3_gene",
    )

    if input_mode == "HGVS notation":
        hgvs = st.sidebar.text_input(
            "HGVS notation", value="",
            help="e.g., NM_000267.3:c.910C>T",
            placeholder="NM_xxx:c.NNNRef>Alt",
            key="v3_hgvs",
        )
        # Still need chromosome and position for GenomicVariantInput
        st.sidebar.caption(
            "Note: HGVS parsing requires genomic coordinates. "
            "Enter them below or leave defaults."
        )

    chromosome = st.sidebar.selectbox(
        "Chromosome", options=CHROMOSOME_OPTIONS,
        index=16, key="v3_chrom",
    )
    position = st.sidebar.number_input(
        "Position (GRCh38)", min_value=1, max_value=300_000_000,
        value=31232193, step=1, key="v3_pos",
    )
    col1, col2 = st.sidebar.columns(2)
    with col1:
        ref_allele = st.text_input(
            "Ref", value="C", max_chars=50, key="v3_ref",
        )
    with col2:
        alt_allele = st.text_input(
            "Alt", value="T", max_chars=50, key="v3_alt",
        )

    st.sidebar.divider()

    cell_type = st.sidebar.selectbox(
        "Cell type", options=CELL_TYPE_OPTIONS, index=0, key="v3_cell",
    )

    nuclease = st.sidebar.selectbox(
        "Primary nuclease", options=NUCLEASE_OPTIONS, index=0,
        key="v3_nuclease",
        help="All 5 nucleases are always evaluated for base editing. "
        "This selects the primary nuclease for HDR/PE guide search.",
    )
    pam = NUCLEASE_PAMS.get(nuclease, "")
    st.sidebar.caption(f"PAM: {pam}")

    variant_name = st.sidebar.text_input(
        "Variant name (optional)", value="",
        placeholder="e.g., c.910C>T", key="v3_name",
    )

    analyze_clicked = st.sidebar.button(
        "Analyze Variant",
        key="v3_analyze_btn",
        type="primary",
        use_container_width=True,
    )

    if analyze_clicked:
        return {
            "gene": gene.strip(),
            "chromosome": str(chromosome).strip(),
            "position": int(position),
            "ref_allele": ref_allele.strip().upper(),
            "alt_allele": alt_allele.strip().upper(),
            "cell_type": cell_type,
            "nuclease": nuclease,
            "name": variant_name.strip() or None,
        }
    return None


# ---------------------------------------------------------------------------
# Display: Variant Annotation
# ---------------------------------------------------------------------------
def _display_variant_annotation(result: PipelineResult) -> None:
    if not result.variants:
        st.warning("No variants were annotated.")
        return

    nv = result.variants[0]
    consequence = nv.coding.consequence.value.replace("_", " ").title()
    hgvs_c = nv.coding.hgvs_c or "N/A"
    hgvs_p = nv.coding.hgvs_p or "N/A"
    ref_valid = nv.ref_validation.is_valid

    rows = [
        ("Consequence", f"<strong>{consequence}</strong>"),
        ("HGVS coding", f'<code>{hgvs_c}</code>'),
        ("HGVS protein", f'<code>{hgvs_p}</code>'),
        ("Transcript", f'<code>{nv.transcript.transcript_id}</code> '
         f'({nv.transcript.gene_symbol})'),
        ("Ref validation", f'{_status_dot(ref_valid)} '
         f'{"Valid" if ref_valid else "Mismatch"}'
         f'{" -- " + nv.ref_validation.message if nv.ref_validation.message else ""}'),
    ]
    if nv.coding.splice_proximity:
        rows.append(("Splice proximity",
                      f'<span style="color:#FF9800; font-weight:600;">'
                      f'{nv.coding.splice_proximity.replace("_", " ").title()}</span>'))
    if nv.coding.message:
        rows.append(("Detail", nv.coding.message))

    html = '<table style="width:100%; border-collapse:collapse;">'
    for label, val in rows:
        html += (
            f'<tr><td style="padding:5px 12px; opacity:0.7; width:140px; '
            f'color:var(--text-color,#E0E0E0);">{label}</td>'
            f'<td style="padding:5px 12px; color:var(--text-color,#E0E0E0);">'
            f'{val}</td></tr>'
        )
    html += '</table>'
    _card("Variant Annotation", html)


# ---------------------------------------------------------------------------
# Display: Strategy Ranking with TOPSIS details
# ---------------------------------------------------------------------------
def _display_strategy_ranking(result: PipelineResult) -> None:
    if not result.strategies:
        st.info("No viable strategies found. Check rejected strategies and warnings below.")
        return

    st.markdown(
        '<div class="section-header">Ranked Strategies (TOPSIS 6D)</div>',
        unsafe_allow_html=True,
    )

    for scored in result.strategies:
        s = scored.strategy
        rank = scored.rank

        if rank == 1:
            border_color = "#4CAF50"
            label_prefix = "RECOMMENDED"
        elif scored.safety_score < 0.35:
            border_color = "#F44336"
            label_prefix = ""
        else:
            border_color = "#455A64"
            label_prefix = ""

        # Rank stability
        stability = getattr(scored, 'rank_stability', 0.0)
        stability_html = _stability_bar(stability) if stability > 0 else ""

        # Pareto status from annotation notes
        pareto_status = ""
        for note in scored.annotation_notes:
            if "pareto" in note.lower():
                if "non-dominated" in note.lower():
                    pareto_status = (
                        '<span style="background:#1B9E77; color:white; '
                        'padding:1px 6px; border-radius:8px; font-size:0.75em; '
                        'margin-left:6px;">Pareto optimal</span>'
                    )
                else:
                    pareto_status = (
                        '<span style="background:#757575; color:white; '
                        'padding:1px 6px; border-radius:8px; font-size:0.75em; '
                        'margin-left:6px;">Dominated</span>'
                    )

        # 6D scores
        scores_html = (
            _score_pill("Safety", scored.safety_score)
            + _score_pill("Feasibility", scored.feasibility_score)
            + _score_pill("Complexity", scored.complexity_score, is_cost=True)
            + _score_pill("Risk", scored.risk_score, is_cost=True)
            + _score_pill("Confidence", scored.confidence_score)
        )

        # Strategy details
        details = []
        if s.num_dsbs > 0:
            details.append(f"DSBs: {s.num_dsbs}")
        if s.num_donors > 0:
            details.append(f"Donors: {s.num_donors}")
        if s.num_rounds > 1:
            details.append(f"Rounds: {s.num_rounds}")
        if s.requires_selection:
            details.append("Selection required")

        # Non-pareto annotation notes
        other_notes = [n for n in scored.annotation_notes
                       if "pareto" not in n.lower()]

        html = (
            f'<div style="border-left:4px solid {border_color}; '
            f'background:var(--card-bg, #233448); padding:12px 16px; '
            f'border-radius:6px; margin-bottom:10px;">'
            # Header row
            f'<div style="display:flex; justify-content:space-between; '
            f'align-items:center; flex-wrap:wrap; gap:6px;">'
            f'<div>'
            f'<strong style="font-size:1.15em; color:var(--text-color,#E0E0E0);">'
            f'#{rank} {scored.strategy_name}</strong>'
            f'{pareto_status}'
        )
        if label_prefix:
            html += (
                f' <span style="background:{border_color}; color:white; '
                f'padding:1px 8px; border-radius:8px; font-size:0.72em; '
                f'font-weight:700; margin-left:6px;">{label_prefix}</span>'
            )
        html += f'</div><div>{_tier_badge(scored.confidence)}</div></div>'

        # TOPSIS score + stability
        html += (
            f'<div style="display:flex; justify-content:space-between; '
            f'align-items:center; margin-top:8px; flex-wrap:wrap; gap:8px;">'
            f'<div style="font-size:0.9em;">'
            f'<span style="opacity:0.7;">TOPSIS score:</span> '
            f'<strong style="font-size:1.3em; color:var(--text-color,#E0E0E0);">'
            f'{scored.overall_score:.3f}</strong></div>'
        )
        if stability_html:
            html += (
                f'<div style="font-size:0.85em;">'
                f'<span style="opacity:0.7;">Rank stability:</span> '
                f'{stability_html}</div>'
            )
        html += '</div>'

        # 6D dimension scores
        html += (
            f'<div style="margin-top:6px; color:var(--text-color,#E0E0E0);">'
            f'{scores_html}</div>'
        )

        # Details
        if details:
            html += (
                f'<div style="font-size:0.8em; opacity:0.6; margin-top:4px;">'
                f'{" &bull; ".join(details)}</div>'
            )

        # Notes
        if other_notes:
            html += (
                f'<div style="font-size:0.78em; opacity:0.55; margin-top:2px; '
                f'font-style:italic;">{" | ".join(other_notes)}</div>'
            )

        html += '</div>'
        st.markdown(html, unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Display: Delivery Advisory (post-ranking annotations)
# ---------------------------------------------------------------------------
def _delivery_complexity_badge(complexity: int) -> str:
    """Badge for delivery complexity (1-5 ordinal)."""
    colors = {1: "#4CAF50", 2: "#8BC34A", 3: "#FF9800", 4: "#F44336", 5: "#B71C1C"}
    labels = {1: "Simple", 2: "Low", 3: "Moderate", 4: "Complex", 5: "Very Complex"}
    color = colors.get(complexity, "#757575")
    label = labels.get(complexity, "?")
    return (
        f'<span style="background:{color}; color:white; padding:1px 8px; '
        f'border-radius:10px; font-size:0.78em; font-weight:600;">'
        f'Delivery: {label} ({complexity}/5)</span>'
    )


def _display_delivery_advisory(result: PipelineResult) -> None:
    """Display delivery advisor annotations from pipeline metadata."""
    delivery = result.metadata.get("delivery_advisory", {})
    if not delivery or "error" in delivery:
        return

    annotations = delivery.get("annotations", [])
    if not annotations:
        return

    st.markdown(
        '<div class="section-header">Delivery Recommendations</div>',
        unsafe_allow_html=True,
    )

    # Global warnings
    global_warns = delivery.get("global_warnings", [])
    for gw in global_warns:
        st.info(gw)

    for ann in annotations:
        strategy_name = ann.get("strategy", "")
        deliverable = ann.get("deliverable", True)
        method = ann.get("method", "")
        complexity = ann.get("complexity", 1)
        donor_format = ann.get("donor_format")
        warnings = ann.get("warnings", [])
        enhancers = ann.get("enhancers", [])
        violations = ann.get("violations", [])

        # Skip if no meaningful delivery info
        if not method and not donor_format and not warnings:
            continue

        border_color = "#4CAF50" if deliverable else "#F44336"
        status_icon = "\u2705" if deliverable else "\u274C"

        html = (
            f'<div style="border-left:3px solid {border_color}; '
            f'background:var(--card-bg, #233448); padding:10px 14px; '
            f'border-radius:5px; margin-bottom:8px;">'
            f'<div style="display:flex; justify-content:space-between; '
            f'align-items:center; flex-wrap:wrap; gap:4px;">'
            f'<strong style="color:var(--text-color,#E0E0E0);">'
            f'{status_icon} {strategy_name}</strong>'
            f'<div>{_delivery_complexity_badge(complexity)}</div></div>'
        )

        # Delivery method
        if method:
            html += (
                f'<div style="font-size:0.85em; margin-top:6px; '
                f'color:var(--text-color,#E0E0E0);">'
                f'<span style="opacity:0.6;">Method:</span> {method}</div>'
            )

        # Donor format
        if donor_format:
            fmt_colors = {
                "ssODN": "#4CAF50", "cssDNA": "#1B9E77",
                "lssDNA": "#FF9800", "dsDNA": "#F44336", "AAV6": "#3498DB",
            }
            fmt_color = fmt_colors.get(donor_format, "#757575")
            html += (
                f'<div style="font-size:0.85em; margin-top:3px; '
                f'color:var(--text-color,#E0E0E0);">'
                f'<span style="opacity:0.6;">Donor format:</span> '
                f'<span style="background:{fmt_color}; color:white; '
                f'padding:1px 8px; border-radius:8px; font-size:0.9em; '
                f'font-weight:600;">{donor_format}</span></div>'
            )

        # Violations
        for v in violations:
            html += (
                f'<div style="font-size:0.82em; margin-top:4px; '
                f'color:#F44336;">\u26A0 {v}</div>'
            )

        # Warnings
        for w in warnings:
            html += (
                f'<div style="font-size:0.82em; margin-top:3px; '
                f'color:#FF9800;">\u26A0 {w}</div>'
            )

        # Viability enhancers (collapsed)
        if enhancers:
            html += (
                f'<div style="font-size:0.78em; margin-top:6px; '
                f'opacity:0.65; color:var(--text-color,#E0E0E0);">'
                f'<strong>Viability tips:</strong> '
                f'{" &bull; ".join(enhancers)}</div>'
            )

        html += '</div>'
        st.markdown(html, unsafe_allow_html=True)

    # Evidence note
    st.caption(
        "Delivery annotations are post-ranking (do not change TOPSIS scores). "
        "Evidence: 87 verified references. "
        "Key: Iyer 2022, Xie 2024, Letort 2025, Ihry 2018, Dever 2016."
    )


# ---------------------------------------------------------------------------
# Display: Feasibility Breakdown
# ---------------------------------------------------------------------------
def _display_feasibility(result: PipelineResult) -> None:
    if not result.bundles:
        return

    st.markdown(
        '<div class="section-header">Feasibility Breakdown</div>',
        unsafe_allow_html=True,
    )

    for bundle in result.bundles:
        nv = bundle.variant
        variant_label = nv.input.name or (
            f"{nv.input.gene_symbol or ''} "
            f"chr{nv.input.chromosome}:{nv.input.position} "
            f"{nv.input.ref_allele}>{nv.input.alt_allele}"
        )

        # Base editing — show all results, not just best
        be_rows = ""
        if bundle.base_editing_results:
            for be in bundle.base_editing_results:
                badge = _feasibility_badge(be.label)
                editor = be.editor_type or "?"
                nuclease_used = getattr(be, 'nuclease', '')
                bystander = f", {be.bystander_count} bystander(s)" if be.bystander_count else ""
                tier = getattr(be, 'evidence_tier', '')
                tier_str = f" [Tier {tier}]" if tier else ""
                nuc_str = f" + {nuclease_used}" if nuclease_used else ""
                be_rows += (
                    f'<tr><td style="padding:4px 12px; opacity:0.7; '
                    f'color:var(--text-color,#E0E0E0); padding-left:24px;">'
                    f'{editor}{nuc_str}</td>'
                    f'<td style="padding:4px 12px; color:var(--text-color,#E0E0E0);">'
                    f'{badge}{bystander}{tier_str}</td></tr>'
                )
        else:
            be_rows = (
                '<tr><td style="padding:4px 12px; opacity:0.7; '
                'color:var(--text-color,#E0E0E0); padding-left:24px;">'
                'All editors</td>'
                f'<td style="padding:4px 12px; color:var(--text-color,#E0E0E0);">'
                f'{_feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)} (not assessed)</td></tr>'
            )

        # Prime editing
        pe = bundle.prime_editing_result
        if pe and pe.label != FeasibilityLabel.NOT_FEASIBLE:
            pe_badge = _feasibility_badge(pe.label)
            pe_detail = f" (PBS {pe.pbs_length}nt, RT {pe.rt_template_length}nt)"
        elif pe:
            pe_badge = _feasibility_badge(pe.label)
            pe_detail = f" ({pe.rejection_reason})" if pe.rejection_reason else ""
        else:
            pe_badge = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
            pe_detail = ""

        # HDR
        hdr = bundle.hdr_result
        if hdr and hdr.label != FeasibilityLabel.NOT_FEASIBLE:
            hdr_badge = _feasibility_badge(hdr.label)
            hdr_detail = (f" ({hdr.recommended_donor_type}, "
                          f"cut-to-edit {hdr.cut_to_edit_distance}bp)")
        elif hdr:
            hdr_badge = _feasibility_badge(hdr.label)
            hdr_detail = f" ({hdr.rejection_reason})" if hdr.rejection_reason else ""
        else:
            hdr_badge = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
            hdr_detail = ""

        html = (
            f'<div style="font-size:0.9em; margin-bottom:4px; font-weight:600; '
            f'color:var(--text-color,#E0E0E0);">{variant_label}</div>'
            f'<table style="width:100%; border-collapse:collapse;">'
            f'<tr><td colspan="2" style="padding:6px 12px; font-weight:600; '
            f'color:#1B9E77;">Base Editing (multi-nuclease)</td></tr>'
            f'{be_rows}'
            f'<tr><td style="padding:6px 12px; opacity:0.7; width:180px; '
            f'color:var(--text-color,#E0E0E0);">Prime Editing</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{pe_badge}{pe_detail}</td></tr>'
            f'<tr><td style="padding:6px 12px; opacity:0.7; '
            f'color:var(--text-color,#E0E0E0);">HDR</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{hdr_badge}{hdr_detail}</td></tr>'
            f'</table>'
        )
        _card("Feasibility", html)


# ---------------------------------------------------------------------------
# Display: Rejected Strategies
# ---------------------------------------------------------------------------
def _display_rejected(result: PipelineResult) -> None:
    if not result.rejected_strategies:
        return

    with st.expander(f"Rejected strategies ({len(result.rejected_strategies)})"):
        for s in result.rejected_strategies:
            reasons = "; ".join(s.rejection_reasons) if s.rejection_reasons else "Unknown"
            st.markdown(
                f"- **{s.name}**: {reasons}",
            )


# ---------------------------------------------------------------------------
# Display: Pipeline Summary
# ---------------------------------------------------------------------------
def _display_summary(result: PipelineResult) -> None:
    meta = result.metadata
    if not meta:
        return

    cols = st.columns(4)
    with cols[0]:
        n_ranked = meta.get("n_strategies_ranked", len(result.strategies))
        st.metric("Strategies ranked", n_ranked)
    with cols[1]:
        n_rejected = meta.get("n_strategies_rejected", len(result.rejected_strategies))
        st.metric("Rejected", n_rejected)
    with cols[2]:
        st.metric("Scoring method", meta.get("scoring_method", "TOPSIS"))
    with cols[3]:
        sens = meta.get("sensitivity_runs", 0)
        st.metric("Sensitivity runs", f"{sens:,}" if sens else "N/A")


# ---------------------------------------------------------------------------
# Display: Warnings & Metadata
# ---------------------------------------------------------------------------
def _display_warnings(result: PipelineResult) -> None:
    if result.warnings:
        with st.expander("Warnings"):
            for w in result.warnings:
                st.warning(w)


def _display_metadata(result: PipelineResult) -> None:
    if result.metadata:
        with st.expander("Pipeline metadata"):
            for k, v in result.metadata.items():
                st.text(f"{k}: {v}")


# ---------------------------------------------------------------------------
# Main page
# ---------------------------------------------------------------------------
def page_v3_analysis() -> None:
    """Render the v3 Strategy Analysis page."""
    st.title("CRISPRArchitect v3: Strategy Analysis")
    st.markdown(
        "Multi-nuclease, TOPSIS-scored strategy recommendation. "
        "Configure variant and parameters in the **sidebar**, then click "
        "**Analyze Variant**."
    )

    if not _V3_AVAILABLE:
        st.error(
            f"v3 pipeline modules are not available.\n\n"
            f"Import error: `{_V3_IMPORT_ERROR}`\n\n"
            f"Ensure the `core/` package is on the Python path."
        )
        return

    if "v3_result" not in st.session_state:
        # Show feature highlights when no result yet
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("**5 Nucleases**")
            st.caption("SpCas9, enFnCas9, SpCas9-NG, SpRY, Cas12a")
        with c2:
            st.markdown("**9 Editor Profiles**")
            st.caption("3 Tier A + 6 Tier B fusions")
        with c3:
            st.markdown("**TOPSIS 6D + Sensitivity**")
            st.caption("10,000 weight permutations")
        st.info("Enter variant details in the sidebar and click **Analyze Variant**.")
        return

    result: PipelineResult = st.session_state["v3_result"]

    # Summary metrics
    _display_summary(result)
    st.markdown("---")

    # Variant annotation
    _display_variant_annotation(result)

    # Strategy ranking
    _display_strategy_ranking(result)

    # Delivery advisor annotations
    _display_delivery_advisory(result)

    # Feasibility
    _display_feasibility(result)

    # Rejected strategies
    _display_rejected(result)

    # Warnings & metadata
    _display_warnings(result)
    _display_metadata(result)


def run_v3_pipeline(params: dict) -> None:
    """Run the v3 pipeline and store results in session state."""
    if not _V3_AVAILABLE:
        st.error("v3 pipeline modules are not available.")
        return

    errors = []
    if not params["gene"]:
        errors.append("Gene symbol is required.")
    if not params["ref_allele"]:
        errors.append("Reference allele is required.")
    if not params["alt_allele"]:
        errors.append("Alternate allele is required.")
    if params["ref_allele"] == params["alt_allele"]:
        errors.append("Ref and Alt alleles must differ.")
    if errors:
        for err in errors:
            st.error(err)
        return

    try:
        with st.spinner(
            f"Running v3 pipeline: {params['gene']} "
            f"chr{params['chromosome']}:{params['position']} "
            f"{params['ref_allele']}>{params['alt_allele']} "
            f"({params['nuclease']}, {params['cell_type']})..."
        ):
            pipeline = StrategyPipeline(
                cell_type=params["cell_type"],
                nuclease=params["nuclease"],
            )
            variant = GenomicVariantInput(
                chromosome=params["chromosome"],
                position=params["position"],
                ref_allele=params["ref_allele"],
                alt_allele=params["alt_allele"],
                gene_symbol=params["gene"],
                name=params.get("name"),
            )
            result = pipeline.run([variant])

        st.session_state["v3_result"] = result

    except Exception as e:
        st.error(
            f"Pipeline failed: {e}\n\n"
            "Check: (1) gene symbol and coordinates are correct, "
            "(2) Ensembl API is reachable, (3) ref allele matches GRCh38."
        )
        if "v3_result" in st.session_state:
            del st.session_state["v3_result"]
