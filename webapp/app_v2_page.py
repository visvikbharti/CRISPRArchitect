"""
CRISPRArchitect v2 — Variant Analysis Page
=============================================

Streamlit page component for the v2 variant analysis pipeline.
Accepts genomic variant inputs and displays consequence-aware
strategy recommendations.

Compatible with Python 3.9+.
"""

from __future__ import annotations

import streamlit as st
from typing import List, Optional

# ---------------------------------------------------------------------------
# Try importing v2 modules — graceful fallback if not available
# ---------------------------------------------------------------------------
_V2_AVAILABLE = False
_V2_IMPORT_ERROR = ""

try:
    from core.pipeline.strategy_stage import StrategyPipeline
    from core.models import (
        GenomicVariantInput,
        PipelineResult,
        FeasibilityLabel,
        NormalizedVariant,
        FeasibilityBundle,
    )
    _V2_AVAILABLE = True
except ImportError as e:
    _V2_IMPORT_ERROR = str(e)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CELL_TYPE_OPTIONS = ["iPSC", "HEK293T", "K562"]
NUCLEASE_OPTIONS = ["SpCas9", "enFnCas9"]

CHROMOSOME_OPTIONS = [str(i) for i in range(1, 23)] + ["X", "Y"]


# ---------------------------------------------------------------------------
# Helper: render a styled card using raw HTML (theme-aware)
# ---------------------------------------------------------------------------
def _card(title: str, content_html: str) -> None:
    """Render a theme-aware card matching the existing ca-card style."""
    st.markdown(
        f'<div class="ca-card">'
        f'<h3>{title}</h3>'
        f'{content_html}'
        f'</div>',
        unsafe_allow_html=True,
    )


def _status_dot(ok: bool) -> str:
    """Green or red status indicator."""
    color = "#4CAF50" if ok else "#F44336"
    return f'<span style="color:{color}; font-size:1.2em;">&#9679;</span>'


def _feasibility_badge(label: FeasibilityLabel) -> str:
    """Colored badge for feasibility status."""
    color_map = {
        FeasibilityLabel.FEASIBLE: "#4CAF50",
        FeasibilityLabel.MARGINAL: "#FF9800",
        FeasibilityLabel.NOT_FEASIBLE: "#F44336",
    }
    color = color_map.get(label, "#757575")
    text = label.value.replace("_", " ").title()
    return (
        f'<span style="background:{color}; color:white; padding:2px 10px; '
        f'border-radius:12px; font-size:0.85em; font-weight:600;">'
        f'{text}</span>'
    )


# ---------------------------------------------------------------------------
# Sidebar inputs
# ---------------------------------------------------------------------------
def render_v2_sidebar() -> Optional[dict]:
    """Render the v2 sidebar inputs and return parameters on Analyze click.

    Returns None if user has not clicked Analyze yet.
    """
    st.sidebar.divider()
    st.sidebar.markdown(
        '<div class="section-header" style="font-size:1.05rem;">'
        'v2: Variant Analysis</div>',
        unsafe_allow_html=True,
    )

    gene = st.sidebar.text_input(
        "Gene symbol",
        value="NF1",
        help="HGNC gene symbol (e.g., NF1, BRCA1, TP53)",
        key="v2_gene",
    )
    chromosome = st.sidebar.selectbox(
        "Chromosome",
        options=CHROMOSOME_OPTIONS,
        index=16,  # default "17" for NF1
        key="v2_chrom",
    )
    position = st.sidebar.number_input(
        "Position (GRCh38)",
        min_value=1,
        max_value=300_000_000,
        value=31232193,
        step=1,
        help="1-based genomic coordinate",
        key="v2_pos",
    )
    ref_allele = st.sidebar.text_input(
        "Ref allele",
        value="C",
        max_chars=50,
        help='Reference allele. Use "-" for insertion.',
        key="v2_ref",
    )
    alt_allele = st.sidebar.text_input(
        "Alt allele",
        value="T",
        max_chars=50,
        help='Alternate (patient) allele. Use "-" for deletion.',
        key="v2_alt",
    )
    cell_type = st.sidebar.selectbox(
        "Cell type",
        options=CELL_TYPE_OPTIONS,
        index=0,
        key="v2_cell",
    )
    nuclease = st.sidebar.selectbox(
        "Nuclease",
        options=NUCLEASE_OPTIONS,
        index=0,
        key="v2_nuclease",
    )

    analyze_clicked = st.sidebar.button(
        "Analyze",
        key="v2_analyze_btn",
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
        }
    return None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _validate_inputs(params: dict) -> List[str]:
    """Return a list of validation error messages (empty if all OK)."""
    errors: List[str] = []
    if not params["gene"]:
        errors.append("Gene symbol is required.")
    if not params["ref_allele"]:
        errors.append("Reference allele is required.")
    if not params["alt_allele"]:
        errors.append("Alternate allele is required.")
    if params["ref_allele"] == params["alt_allele"]:
        errors.append("Ref and Alt alleles must differ.")
    if params["position"] < 1:
        errors.append("Position must be a positive integer.")
    return errors


# ---------------------------------------------------------------------------
# Display functions
# ---------------------------------------------------------------------------
def _display_variant_annotation(result: PipelineResult) -> None:
    """Show variant annotation card for the first variant."""
    if not result.variants:
        st.warning("No variants were annotated.")
        return

    nv = result.variants[0]

    consequence = nv.coding.consequence.value.replace("_", " ").title()
    hgvs_c = nv.coding.hgvs_c or "N/A"
    hgvs_p = nv.coding.hgvs_p or "N/A"
    ref_valid = nv.ref_validation.is_valid
    ref_msg = nv.ref_validation.message or ""

    annotation_html = (
        f'<table style="width:100%; border-collapse:collapse;">'
        f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
        f'opacity:0.7;">Consequence</td>'
        f'<td style="padding:6px 12px; font-weight:600; '
        f'color:var(--text-color,#E0E0E0);">{consequence}</td></tr>'
        f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
        f'opacity:0.7;">HGVS (coding)</td>'
        f'<td style="padding:6px 12px; font-family:monospace; '
        f'color:var(--text-color,#E0E0E0);">{hgvs_c}</td></tr>'
        f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
        f'opacity:0.7;">HGVS (protein)</td>'
        f'<td style="padding:6px 12px; font-family:monospace; '
        f'color:var(--text-color,#E0E0E0);">{hgvs_p}</td></tr>'
        f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
        f'opacity:0.7;">Ref validation</td>'
        f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
        f'{_status_dot(ref_valid)} '
        f'{"Valid" if ref_valid else "Mismatch"}'
        f'{" — " + ref_msg if ref_msg else ""}</td></tr>'
    )

    # Splice proximity
    if nv.coding.splice_proximity:
        splice_info = nv.coding.splice_proximity.replace("_", " ").title()
        annotation_html += (
            f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
            f'opacity:0.7;">Splice proximity</td>'
            f'<td style="padding:6px 12px; color:#FF9800; font-weight:600;">'
            f'{splice_info}</td></tr>'
        )

    # Coding detail
    if nv.coding.message:
        annotation_html += (
            f'<tr><td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
            f'opacity:0.7;">Detail</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{nv.coding.message}</td></tr>'
        )

    annotation_html += '</table>'

    _card("Variant Annotation", annotation_html)


def _display_strategy_table(result: PipelineResult) -> None:
    """Show ranked strategy table."""
    if not result.strategies:
        st.info("No strategies were generated. Check warnings below.")
        return

    st.markdown(
        '<div class="section-header">Ranked Strategies</div>',
        unsafe_allow_html=True,
    )

    for scored in result.strategies:
        s = scored.strategy
        rank = scored.rank
        name = scored.strategy_name

        # Determine CSS class for styling
        if rank == 1:
            css_class = "strategy-best"
        elif scored.safety_score < 0.35:
            css_class = "strategy-risky"
        else:
            css_class = "strategy-alt"

        # Build the row
        tier_colors = {"A": "#4CAF50", "B": "#FF9800", "C": "#F44336"}
        tier = scored.confidence
        tier_color = tier_colors.get(tier, "#757575")

        row_html = (
            f'<div class="{css_class}" style="margin-bottom:0.6rem;">'
            f'<div style="display:flex; justify-content:space-between; '
            f'align-items:center; flex-wrap:wrap; gap:0.5rem;">'
            f'<div>'
            f'<strong style="font-size:1.1em;">#{rank} {name}</strong>'
            f'</div>'
            f'<div style="display:flex; gap:1rem; align-items:center;">'
            f'<span style="font-size:0.85em; opacity:0.8;">Overall</span> '
            f'<strong style="font-size:1.2em;">'
            f'{scored.overall_score:.2f}</strong>'
            f'<span style="font-size:0.85em; opacity:0.8;">Safety</span> '
            f'<strong>{scored.safety_score:.2f}</strong>'
            f'<span style="font-size:0.85em; opacity:0.8;">Feasibility</span> '
            f'<strong>{scored.feasibility_score:.2f}</strong>'
            f'<span style="background:{tier_color}; color:white; '
            f'padding:1px 8px; border-radius:10px; font-size:0.8em; '
            f'font-weight:700;">Tier {tier}</span>'
            f'</div>'
            f'</div>'
        )

        # Annotation notes
        if scored.annotation_notes:
            notes_str = " | ".join(scored.annotation_notes)
            row_html += (
                f'<div style="font-size:0.82em; opacity:0.7; margin-top:4px;">'
                f'{notes_str}</div>'
            )

        # Strategy details
        details = []
        if s.num_dsbs > 0:
            details.append(f"DSBs: {s.num_dsbs}")
        if s.num_donors > 0:
            details.append(f"Donors: {s.num_donors}")
        if s.num_rounds > 1:
            details.append(f"Rounds: {s.num_rounds}")
        if s.simultaneous_dsbs:
            details.append("Simultaneous DSBs")
        if s.requires_selection:
            details.append("Requires selection")
        if details:
            row_html += (
                f'<div style="font-size:0.82em; opacity:0.65; margin-top:2px;">'
                f'{" &bull; ".join(details)}</div>'
            )

        row_html += '</div>'

        st.markdown(row_html, unsafe_allow_html=True)


def _display_feasibility_breakdown(result: PipelineResult) -> None:
    """Show feasibility breakdown per modality for each variant."""
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

        # Base editing
        be_result = bundle.best_base_editing_result()
        if be_result:
            be_label = _feasibility_badge(be_result.label)
            be_detail = ""
            if be_result.editor_type:
                be_detail += f" ({be_result.editor_type}"
                if be_result.bystander_count > 0:
                    be_detail += f", {be_result.bystander_count} bystander(s)"
                be_detail += ")"
        else:
            # Check if any results exist at all
            if bundle.base_editing_results:
                be_label = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
                reason = bundle.base_editing_results[0].rejection_reason
                be_detail = f" ({reason})" if reason else ""
            else:
                be_label = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
                be_detail = " (not assessed)"

        # Prime editing
        pe = bundle.prime_editing_result
        if pe:
            pe_label = _feasibility_badge(pe.label)
            pe_detail = ""
            if pe.label != FeasibilityLabel.NOT_FEASIBLE:
                pe_detail = f" (PBS {pe.pbs_length}nt, RT {pe.rt_template_length}nt)"
            elif pe.rejection_reason:
                pe_detail = f" ({pe.rejection_reason})"
        else:
            pe_label = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
            pe_detail = " (not assessed)"

        # HDR
        hdr = bundle.hdr_result
        if hdr:
            hdr_label = _feasibility_badge(hdr.label)
            hdr_detail = ""
            if hdr.label != FeasibilityLabel.NOT_FEASIBLE:
                hdr_detail = (
                    f" ({hdr.recommended_donor_type}, "
                    f"cut-to-edit {hdr.cut_to_edit_distance}bp)"
                )
            elif hdr.rejection_reason:
                hdr_detail = f" ({hdr.rejection_reason})"
        else:
            hdr_label = _feasibility_badge(FeasibilityLabel.NOT_FEASIBLE)
            hdr_detail = " (not assessed)"

        card_html = (
            f'<div style="font-size:0.9em; margin-bottom:4px; font-weight:600; '
            f'color:var(--text-color,#E0E0E0);">{variant_label}</div>'
            f'<table style="width:100%; border-collapse:collapse;">'
            f'<tr>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
            f'opacity:0.7; width:140px;">Base Editing</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{be_label}{be_detail}</td></tr>'
            f'<tr>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
            f'opacity:0.7;">Prime Editing</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{pe_label}{pe_detail}</td></tr>'
            f'<tr>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0); '
            f'opacity:0.7;">HDR</td>'
            f'<td style="padding:6px 12px; color:var(--text-color,#E0E0E0);">'
            f'{hdr_label}{hdr_detail}</td></tr>'
            f'</table>'
        )

        _card("Feasibility", card_html)


def _display_warnings(result: PipelineResult) -> None:
    """Display pipeline warnings."""
    if result.warnings:
        st.markdown(
            '<div class="section-header">Warnings</div>',
            unsafe_allow_html=True,
        )
        for w in result.warnings:
            st.warning(w)


def _display_metadata(result: PipelineResult) -> None:
    """Show pipeline metadata in an expander."""
    if result.metadata:
        with st.expander("Pipeline metadata"):
            for k, v in result.metadata.items():
                st.text(f"{k}: {v}")


# ---------------------------------------------------------------------------
# Main page renderer
# ---------------------------------------------------------------------------
def page_v2_analysis() -> None:
    """Render the full v2 Variant Analysis page."""
    st.title("v2: Variant Analysis")
    st.markdown(
        "Consequence-aware strategy recommendation using the v2 pipeline. "
        "Enter a genomic variant in the sidebar and click **Analyze**."
    )
    st.markdown("---")

    # Check v2 availability
    if not _V2_AVAILABLE:
        st.error(
            f"v2 pipeline modules are not available.\n\n"
            f"Import error: `{_V2_IMPORT_ERROR}`\n\n"
            f"Ensure the `core/` package is installed and on the Python path."
        )
        return

    # Check if results are in session state (from sidebar Analyze click)
    if "v2_result" not in st.session_state:
        st.info(
            "Configure variant parameters in the **sidebar** under "
            '"v2: Variant Analysis" and click **Analyze** to begin.'
        )
        return

    result: PipelineResult = st.session_state["v2_result"]

    # Display all result sections
    _display_variant_annotation(result)
    _display_strategy_table(result)
    _display_feasibility_breakdown(result)
    _display_warnings(result)
    _display_metadata(result)


def run_v2_pipeline(params: dict) -> None:
    """Run the v2 pipeline and store results in session state.

    Called from the main app when the Analyze button is clicked.
    Shows a spinner while Ensembl API calls execute.
    """
    if not _V2_AVAILABLE:
        st.error("v2 pipeline modules are not available.")
        return

    errors = _validate_inputs(params)
    if errors:
        for err in errors:
            st.error(err)
        return

    try:
        with st.spinner("Running v2 pipeline (fetching from Ensembl API)..."):
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
            )
            result = pipeline.run([variant])

        st.session_state["v2_result"] = result

    except Exception as e:
        st.error(
            f"Pipeline execution failed: {e}\n\n"
            "This may be due to a network issue reaching the Ensembl API, "
            "or the gene/variant may not be found. Please verify your inputs."
        )
        # Clear stale result
        if "v2_result" in st.session_state:
            del st.session_state["v2_result"]
