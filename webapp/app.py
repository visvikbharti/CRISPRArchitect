#!/usr/bin/env python3
"""
CRISPRArchitect — Multi-Site Genome Editing Optimizer
=====================================================

A Streamlit web application providing an interactive interface to the
CRISPRArchitect scientific toolkit.  Supports five analysis modules:

  1. Gene & Mutation Setup   (shared state across all pages)
  2. MOSAIC Strategy Optimizer
  3. ConversionSim Tract Simulator
  4. ChromBridge 3D Distance & Risk
  5. cssDNA-TopoPred Donor Quality Check
  6. Full Analysis Report

Launch
------
    cd crisprarchitect/
    streamlit run webapp/app.py --server.port 8501
"""

from __future__ import annotations

import io
import random
import sys
import os
import textwrap
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Path setup — ensure the crisprarchitect package is importable
# ---------------------------------------------------------------------------
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_REPO_ROOT = os.path.abspath(os.path.join(_PROJECT_ROOT, ".."))
for _p in (_PROJECT_ROOT, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ---------------------------------------------------------------------------
# Import CRISPRArchitect modules
# ---------------------------------------------------------------------------
from crisprarchitect.mosaic.gene_structure import GeneStructure, build_demo_gene
from crisprarchitect.mosaic.mutation_classifier import Mutation, MutationClassifier
from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
from crisprarchitect.mosaic.scorer import StrategyScorer, ScoredStrategy
from crisprarchitect.mosaic.reporter import StrategyReporter

from crisprarchitect.conversion_sim import ConversionSimulator
from crisprarchitect.chrombridge import ChromatinDistancePredictor
from crisprarchitect.chrombridge.translocation import TranslocationRiskPredictor
from crisprarchitect.topopred import DonorAnalyzer

from crisprarchitect.utils.constants import (
    CELL_TYPE_PARAMS,
    NUCLEASE_PARAMS,
)
from crisprarchitect.utils.ensembl import fetch_gene, EXAMPLE_GENES, EnsemblError

# Local styling helpers — try both import paths so the app works whether
# launched from crisprarchitect/ or from the webapp/ directory itself.
try:
    from webapp.style import (
        inject_global_css, section_header, grade_badge, render_grade,
        score_colour, format_bp, format_pct, metric_row_html,
        strategy_css_class, GRADE_COLOURS, SCORE_COLOUR_MAP,
        PRIMARY, SUCCESS, WARNING, DANGER,
    )
except ImportError:
    from style import (
        inject_global_css, section_header, grade_badge, render_grade,
        score_colour, format_bp, format_pct, metric_row_html,
        strategy_css_class, GRADE_COLOURS, SCORE_COLOUR_MAP,
        PRIMARY, SUCCESS, WARNING, DANGER,
    )

# ---------------------------------------------------------------------------
# Plotly / matplotlib imports
# ---------------------------------------------------------------------------
try:
    import plotly.graph_objects as go
    import plotly.express as px

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===================================================================
# Streamlit page configuration (must be first Streamlit call)
# ===================================================================
st.set_page_config(
    page_title="CRISPRArchitect",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

inject_global_css()


# ===================================================================
# Session state defaults
# ===================================================================
_DEFAULTS = {
    "gene_name": "MyGene",
    "num_exons": 60,
    "gene_input_mode": "Generate demo gene",
    "manual_exons_text": "",
    "cell_type": "iPSC",
    "nuclease": "enFnCas9",
    "mutations": [],
    "gene": None,
    "mosaic_results": None,
    "conversion_results": None,
    "chrombridge_results": None,
    "topopred_results": None,
}
for key, val in _DEFAULTS.items():
    if key not in st.session_state:
        st.session_state[key] = val


# ===================================================================
# Sidebar navigation
# ===================================================================
PAGES = {
    "Home": "home",
    "Gene & Mutation Setup": "setup",
    "Strategy Optimizer (MOSAIC)": "mosaic",
    "Conversion Tract Simulator": "conversion",
    "3D Distance & Risk (ChromBridge)": "chrombridge",
    "Donor Quality Check (TopoPred)": "topopred",
    "Full Analysis Report": "report",
}

with st.sidebar:
    st.markdown(
        f'<h1 style="color:{PRIMARY};">CRISPRArchitect</h1>',
        unsafe_allow_html=True,
    )
    st.caption("Multi-Site Genome Editing Optimizer")

    # --- Dark / Light mode toggle ---
    if "dark_mode" not in st.session_state:
        st.session_state["dark_mode"] = False
    dark_mode = st.toggle(
        "Dark Mode",
        value=st.session_state["dark_mode"],
        key="dark_mode_toggle",
    )
    st.session_state["dark_mode"] = dark_mode

    st.divider()
    page_label = st.radio(
        "Navigate",
        list(PAGES.keys()),
        label_visibility="collapsed",
    )
    page = PAGES[page_label]

    # Quick-status indicators
    st.divider()
    gene_ready = st.session_state.get("gene") is not None
    muts_ready = len(st.session_state.get("mutations", [])) >= 2
    st.markdown(
        f"**Gene defined:** {'Yes' if gene_ready else 'No'}  \n"
        f"**Mutations:** {len(st.session_state.get('mutations', []))}  \n"
        f"**Cell type:** {st.session_state.get('cell_type', 'N/A')}  \n"
        f"**Nuclease:** {st.session_state.get('nuclease', 'N/A')}"
    )

# --- Apply dark/light theme CSS based on toggle ---
if st.session_state.get("dark_mode", False):
    st.markdown("""
    <style>
        /* === DARK MODE === */
        [data-testid="stAppViewContainer"],
        .main .block-container,
        [data-testid="stApp"] {
            background-color: #0E1117 !important;
            color: #E0E0E0 !important;
        }
        [data-testid="stHeader"] {
            background-color: #0E1117 !important;
        }
        [data-testid="stSidebar"] {
            background-color: #161B22 !important;
            border-right: 1px solid #30363D !important;
        }
        [data-testid="stSidebar"] * {
            color: #E0E0E0 !important;
        }
        [data-testid="stSidebar"] h1 {
            color: #64B5F6 !important;
        }
        h1, h2, h3, h4, h5, h6, p, span, label, div, li {
            color: #E0E0E0 !important;
        }
        .stMarkdown, .stText, .stCaption {
            color: #E0E0E0 !important;
        }
        .ca-card {
            background: #1E1E2E !important;
            border-color: #30363D !important;
        }
        .ca-card h3 { color: #E0E0E0 !important; }
        .section-header {
            color: #64B5F6 !important;
            border-bottom-color: #64B5F6 !important;
        }
        .metric-item {
            background: #1E1E2E !important;
            border-color: #30363D !important;
        }
        .metric-item .label { color: #9E9E9E !important; }
        .metric-item .value { color: #E0E0E0 !important; }
        .strategy-best { background: rgba(76,175,80,0.2) !important; }
        .strategy-alt  { background: rgba(255,152,0,0.15) !important; }
        .strategy-risky { background: rgba(244,67,54,0.15) !important; }
        /* Streamlit widgets — inputs */
        .stSelectbox > div > div,
        .stTextInput > div > div > input,
        .stNumberInput > div > div > input,
        .stTextArea > div > textarea,
        .stMultiSelect > div > div {
            background-color: #1E1E2E !important;
            color: #E0E0E0 !important;
            border-color: #30363D !important;
        }
        /* Dropdown / selectbox popup list */
        [data-baseweb="popover"],
        [data-baseweb="popover"] > div,
        [data-baseweb="menu"],
        [data-baseweb="list"],
        [data-baseweb="select"] [role="listbox"],
        ul[role="listbox"],
        div[role="listbox"] {
            background-color: #1E1E2E !important;
            border-color: #30363D !important;
        }
        /* Dropdown option items */
        [data-baseweb="menu"] li,
        [data-baseweb="list"] li,
        [role="option"],
        ul[role="listbox"] li,
        div[role="listbox"] > div {
            background-color: #1E1E2E !important;
            color: #E0E0E0 !important;
        }
        [role="option"]:hover,
        [data-baseweb="menu"] li:hover,
        ul[role="listbox"] li:hover {
            background-color: #2D3748 !important;
            color: #FFFFFF !important;
        }
        /* Selected option highlight */
        [aria-selected="true"],
        [data-baseweb="menu"] li[aria-selected="true"] {
            background-color: #1A365D !important;
            color: #90CDF4 !important;
        }
        /* Radio buttons */
        .stRadio > div {
            color: #E0E0E0 !important;
        }
        .stRadio label span {
            color: #E0E0E0 !important;
        }
        /* Slider */
        .stSlider > div > div > div {
            color: #E0E0E0 !important;
        }
        /* Expanders */
        [data-testid="stExpander"] {
            background-color: #161B22 !important;
            border-color: #30363D !important;
        }
        [data-testid="stExpander"] summary span {
            color: #E0E0E0 !important;
        }
        /* Tables */
        .stDataFrame, .stTable {
            background-color: #1E1E2E !important;
        }
        .stDataFrame th, .stDataFrame td {
            color: #E0E0E0 !important;
            border-color: #30363D !important;
        }
        /* Alert boxes */
        .stInfo, .stSuccess, .stWarning, .stError {
            background-color: #1E1E2E !important;
        }
        /* Number input buttons */
        .stNumberInput button {
            background-color: #2D3748 !important;
            color: #E0E0E0 !important;
            border-color: #30363D !important;
        }
        /* Buttons */
        .stButton > button {
            border-color: #30363D !important;
        }
        /* Plotly dark bg */
        .js-plotly-plot .plotly .main-svg {
            background: #0E1117 !important;
        }
    </style>
    """, unsafe_allow_html=True)
else:
    st.markdown("""
    <style>
        /* === LIGHT MODE === */
        [data-testid="stSidebar"] {
            background-color: #FAFBFC !important;
            border-right: 1px solid #E0E0E0 !important;
        }
        [data-testid="stSidebar"] h1 {
            color: #2196F3 !important;
        }
    </style>
    """, unsafe_allow_html=True)


# ===================================================================
# Utility helpers
# ===================================================================

def _build_gene_from_session() -> Optional[GeneStructure]:
    """Reconstruct or retrieve the GeneStructure from session state."""
    return st.session_state.get("gene")


def _get_mutations() -> List[Mutation]:
    return st.session_state.get("mutations", [])


def _gene_viz_plotly(gene: GeneStructure, mutations: List[Mutation]):
    """Build a Plotly figure showing exon/intron structure with mutation markers."""
    fig = go.Figure()

    # Genomic baseline
    total_start = gene.exons[0].start
    total_end = gene.exons[-1].end
    fig.add_shape(
        type="line",
        x0=total_start, x1=total_end, y0=0, y1=0,
        line=dict(color="#B0BEC5", width=2),
    )

    # Exons as blue bars
    for exon in gene.exons:
        fig.add_shape(
            type="rect",
            x0=exon.start, x1=exon.end, y0=-0.3, y1=0.3,
            fillcolor=PRIMARY, line=dict(color=PRIMARY, width=0),
            opacity=0.85,
        )

    # Mutation markers
    mut_positions = []
    for m in mutations:
        pos = m.position
        mut_positions.append(pos)
        fig.add_shape(
            type="line",
            x0=pos, x1=pos, y0=-0.6, y1=0.6,
            line=dict(color=DANGER, width=2.5),
        )
        fig.add_annotation(
            x=pos, y=0.75,
            text=f"Exon {m.exon_number}<br>{m.ref_allele}>{m.alt_allele}",
            showarrow=False, font=dict(size=10, color=DANGER),
        )

    fig.update_layout(
        height=180,
        margin=dict(l=20, r=20, t=10, b=10),
        xaxis=dict(
            title="Genomic position (bp)",
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            visible=False,
            range=[-1, 1.2],
        ),
        showlegend=False,
        plot_bgcolor="white",
    )
    return fig


def _gene_viz_matplotlib(gene: GeneStructure, mutations: List[Mutation]):
    """Fallback matplotlib gene diagram."""
    fig, ax = plt.subplots(figsize=(12, 2))
    total_start = gene.exons[0].start
    total_end = gene.exons[-1].end
    ax.plot([total_start, total_end], [0, 0], color="#B0BEC5", lw=1.5)

    for exon in gene.exons:
        ax.barh(0, exon.end - exon.start, left=exon.start, height=0.5,
                color=PRIMARY, edgecolor="none", alpha=0.85)

    for m in mutations:
        ax.axvline(m.position, color=DANGER, lw=2, ymin=0.1, ymax=0.9)
        ax.text(m.position, 0.45,
                f"Ex{m.exon_number}\n{m.ref_allele}>{m.alt_allele}",
                ha="center", va="bottom", fontsize=7, color=DANGER)

    ax.set_xlim(total_start - (total_end - total_start) * 0.02,
                total_end + (total_end - total_start) * 0.02)
    ax.set_ylim(-0.6, 0.8)
    ax.set_xlabel("Genomic position (bp)", fontsize=9)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    return fig


# ===================================================================
# PAGE: Home
# ===================================================================

def page_home():
    st.title("CRISPRArchitect")
    st.markdown("### Multi-Site Genome Editing Optimizer")
    st.markdown("---")

    st.markdown("""
    **CRISPRArchitect** is a computational toolkit for planning multi-locus
    genome-editing experiments.  It integrates four complementary analysis
    modules to answer the critical question:

    > *What is the safest, most efficient strategy for correcting two (or
    > more) mutations in a patient's cells?*

    Each module addresses a different aspect of the problem:
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.info(
            "**MOSAIC** -- Strategy Optimizer\n\n"
            "Enumerates every feasible editing strategy (base editing, "
            "prime editing, HDR, hybrids) and ranks them by efficiency, "
            "safety, time, and cost."
        )
        st.info(
            "**ConversionSim** -- Tract Simulator\n\n"
            "Monte Carlo simulation of HDR gene conversion tracts.  "
            "Predicts how far donor-encoded edits can propagate from "
            "the cut site."
        )

    with col2:
        st.info(
            "**ChromBridge** -- 3D Distance & Risk\n\n"
            "Polymer-physics prediction of 3D chromatin distances.  "
            "Evaluates donor bridgeability and translocation risk "
            "when two DSBs are introduced simultaneously."
        )
        st.info(
            "**cssDNA-TopoPred** -- Donor Quality\n\n"
            "Secondary-structure analysis of cssDNA donor templates.  "
            "Detects G-quadruplexes, hairpins, and accessibility "
            "problems that could block HDR."
        )

    st.markdown("---")
    st.markdown("### Quick Start")
    st.markdown(
        "1. Go to **Gene & Mutation Setup** to define your gene and mutations.  \n"
        "2. Run any individual module, or jump to **Full Analysis Report** "
        "for the complete picture."
    )


# ===================================================================
# PAGE: Gene & Mutation Setup
# ===================================================================

def page_setup():
    st.title("Gene & Mutation Setup")
    st.markdown(
        "Define the gene structure, mutation sites, cell type, and nuclease. "
        "These inputs are shared across all analysis modules."
    )

    # ---- Gene definition ----
    section_header("Gene Definition")

    gene_name = st.text_input(
        "Gene name",
        value=st.session_state["gene_name"],
        help="HGNC gene symbol (e.g., COL7A1, USH2A, RYR1).  Used for labelling only.",
    )
    st.session_state["gene_name"] = gene_name

    _gene_input_options = [
        "Generate demo gene",
        "Enter exon coordinates manually",
        "Fetch real gene from Ensembl (GRCh38)",
        "Upload CSV",
    ]
    input_mode = st.radio(
        "How would you like to define the gene structure?",
        _gene_input_options,
        index=_gene_input_options.index(
            st.session_state.get("gene_input_mode", "Generate demo gene")
        ),
        help=(
            "**Demo gene**: A realistic 60-exon gene with human-like exon/intron "
            "sizes (inspired by COL7A1, USH2A).  \n"
            "**Manual entry**: Paste exon coordinates (number, start, end) one per line.  \n"
            "**Ensembl fetch**: Fetch real exon coordinates from the Ensembl REST API "
            "(GRCh38/hg38). Requires an internet connection.  \n"
            "**CSV upload**: Upload a CSV with columns: number, start, end."
        ),
    )
    st.session_state["gene_input_mode"] = input_mode

    gene = None

    if input_mode == "Generate demo gene":
        num_exons = st.slider(
            "Number of exons",
            min_value=1, max_value=400, value=st.session_state["num_exons"],
            help=(
                "Human protein-coding genes range from 1 exon (histones) to "
                "363 exons (titin).  The median human gene has ~8 exons, but "
                "many disease-relevant genes have 30-120+ exons."
            ),
        )
        st.session_state["num_exons"] = num_exons

        if st.button("Generate Gene", type="primary"):
            with st.spinner("Building gene model..."):
                rng = np.random.RandomState(42)
                exons_list = []
                pos = 0
                for i in range(1, num_exons + 1):
                    exon_size = int(rng.lognormal(mean=np.log(150), sigma=0.45))
                    exon_size = max(50, min(exon_size, 800))
                    if i == 1:
                        exon_size = int(exon_size * 2.0)
                    if i == num_exons:
                        exon_size = int(exon_size * 2.5)
                    exons_list.append({"number": i, "start": pos, "end": pos + exon_size})
                    pos += exon_size
                    if i < num_exons:
                        if i == 1:
                            intron = int(rng.lognormal(mean=np.log(80000), sigma=0.4))
                        elif rng.random() < 0.15:
                            intron = int(rng.lognormal(mean=np.log(50000), sigma=0.6))
                        else:
                            intron = int(rng.lognormal(mean=np.log(2000), sigma=0.7))
                        intron = max(80, intron)
                        pos += intron

                gene = GeneStructure.from_manual(gene_name, exons_list)
                st.session_state["gene"] = gene
            st.success(f"Gene **{gene_name}** created with {num_exons} exons.")

    elif input_mode == "Enter exon coordinates manually":
        st.markdown(
            '<p class="field-help">'
            "Enter one exon per line: <code>number  start  end</code> "
            "(whitespace- or comma-separated).  Coordinates are 0-based, "
            "half-open (BED convention)."
            "</p>",
            unsafe_allow_html=True,
        )
        manual_text = st.text_area(
            "Exon coordinates",
            value=st.session_state.get("manual_exons_text", ""),
            height=200,
            placeholder="1  1000  1200\n2  3000  3180\n3  5500  5700",
        )
        st.session_state["manual_exons_text"] = manual_text

        if st.button("Build Gene from Coordinates", type="primary"):
            try:
                exons_list = []
                for line in manual_text.strip().splitlines():
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.replace(",", " ").split()
                    if len(parts) < 3:
                        st.error(f"Invalid line (need 3 values): {line}")
                        return
                    num, start, end = int(parts[0]), int(parts[1]), int(parts[2])
                    exons_list.append({"number": num, "start": start, "end": end})

                if not exons_list:
                    st.error("No exon coordinates provided.")
                    return

                gene = GeneStructure.from_manual(gene_name, exons_list)
                st.session_state["gene"] = gene
                st.success(f"Gene **{gene_name}** created with {len(exons_list)} exons.")
            except Exception as e:
                st.error(f"Error parsing exon coordinates: {e}")

    elif input_mode == "Fetch real gene from Ensembl (GRCh38)":
        # -- Quick-pick from example genes --
        example_labels = ["(type your own)"] + [
            f"{sym} — {desc}" for sym, desc in EXAMPLE_GENES.items()
        ]
        example_choice = st.selectbox(
            "Quick-pick an example gene (optional)",
            example_labels,
            index=0,
            help="Select a well-known disease gene, or type your own symbol below.",
        )

        # Pre-fill from example if selected
        prefill = ""
        if example_choice != "(type your own)":
            prefill = example_choice.split(" — ")[0].strip()

        ensembl_symbol = st.text_input(
            "Gene symbol",
            value=prefill or st.session_state.get("ensembl_symbol", ""),
            placeholder="e.g., NF1, DMD, BRCA2",
            help="Enter an HGNC gene symbol. The canonical transcript will be used.",
        )
        st.session_state["ensembl_symbol"] = ensembl_symbol

        if st.button("Fetch Gene", type="primary"):
            symbol = ensembl_symbol.strip()
            if not symbol:
                st.error("Please enter a gene symbol.")
            else:
                try:
                    with st.spinner(f"Fetching **{symbol}** from Ensembl..."):
                        gene_info = fetch_gene(symbol)

                    # Store the raw GeneInfo for potential later use
                    st.session_state["ensembl_gene_info"] = gene_info

                    # Update the gene name field to match the fetched gene
                    st.session_state["gene_name"] = gene_info.symbol

                    # Convert Ensembl exons to the format GeneStructure expects
                    exons_list = [
                        {"number": e["number"], "start": e["start"], "end": e["end"]}
                        for e in gene_info.exons
                    ]
                    chrom = "chr" + gene_info.chromosome
                    strand_str = "+" if gene_info.strand == 1 else "-"
                    gene = GeneStructure.from_manual(
                        gene_info.symbol, exons_list,
                        chromosome=chrom, strand=strand_str,
                    )
                    st.session_state["gene"] = gene

                    st.success(
                        f"**{gene_info.symbol}** loaded from Ensembl "
                        f"(chr{gene_info.chromosome}:{gene_info.start:,}-{gene_info.end:,}, "
                        f"{gene_info.n_exons} exons, "
                        f"{'forward' if gene_info.strand == 1 else 'reverse'} strand, "
                        f"span {gene_info.span_bp:,} bp)"
                    )
                    st.caption(
                        f"Transcript: {gene_info.canonical_transcript_id} | "
                        f"Ensembl ID: {gene_info.ensembl_id} | "
                        f"{gene_info.description}"
                    )
                except EnsemblError as e:
                    st.error(f"Ensembl fetch failed: {e}")
                except Exception as e:
                    st.error(f"Unexpected error: {e}")

    else:  # CSV upload
        uploaded = st.file_uploader(
            "Upload exon CSV",
            type=["csv"],
            help="CSV must have columns: number, start, end. One row per exon.",
        )
        if uploaded is not None:
            try:
                df = pd.read_csv(uploaded)
                needed = {"number", "start", "end"}
                if not needed.issubset(set(df.columns)):
                    st.error(f"CSV must contain columns: {needed}. Found: {set(df.columns)}")
                else:
                    exons_list = df[["number", "start", "end"]].to_dict("records")
                    gene = GeneStructure.from_manual(gene_name, exons_list)
                    st.session_state["gene"] = gene
                    st.success(f"Gene **{gene_name}** loaded with {len(exons_list)} exons.")
            except Exception as e:
                st.error(f"Error reading CSV: {e}")

    # Show gene summary if available
    gene = _build_gene_from_session()
    if gene is not None:
        with st.expander("Gene Structure Summary", expanded=True):
            st.code(gene.summary(), language="text")

    # ---- Mutation definition ----
    section_header("Mutation Definition")
    st.markdown(
        "Define the pathogenic mutations to correct.  At least two mutations "
        "are needed for multi-locus strategy analysis."
    )

    # Existing mutations table
    mutations = st.session_state.get("mutations", [])
    if mutations:
        mut_data = []
        classifier = MutationClassifier()
        for m in mutations:
            mut_data.append({
                "Exon": m.exon_number,
                "Position": m.position,
                "Ref": m.ref_allele,
                "Alt": m.alt_allele,
                "Type": m.mutation_type,
                "Base-editable": "Yes" if classifier.base_editing_amenable(m) else "No",
                "Prime-editable": "Yes" if classifier.prime_editing_amenable(m) else "No",
            })
        st.dataframe(pd.DataFrame(mut_data), use_container_width=True, hide_index=True)

    st.markdown("**Add a mutation:**")
    col_e, col_r, col_a = st.columns(3)
    with col_e:
        max_exon = len(gene.exons) if gene else 400
        mut_exon = st.number_input(
            "Exon number",
            min_value=1,
            max_value=max_exon,
            value=1,
            help="1-based exon number where this mutation resides.",
        )
    with col_r:
        mut_ref = st.text_input(
            "Reference allele",
            value="G",
            help=(
                "The wild-type (reference) allele.  Single nucleotide for "
                "point mutations (e.g., G); multi-nt for MNVs; '-' for "
                "insertions."
            ),
        )
    with col_a:
        mut_alt = st.text_input(
            "Alternate allele",
            value="A",
            help=(
                "The patient's mutant allele.  The correction edit is "
                "alt -> ref.  Use '-' for deletions."
            ),
        )

    if st.button("Add Mutation"):
        if gene is None:
            st.warning("Please define the gene structure first.")
        else:
            try:
                exon_idx = mut_exon - 1
                if 0 <= exon_idx < len(gene.exons):
                    position = gene.exons[exon_idx].start + gene.exons[exon_idx].size // 2
                else:
                    position = 0
                new_mut = Mutation(
                    exon_number=mut_exon,
                    position=position,
                    ref_allele=mut_ref.strip(),
                    alt_allele=mut_alt.strip(),
                )
                st.session_state["mutations"].append(new_mut)
                st.rerun()
            except Exception as e:
                st.error(f"Error adding mutation: {e}")

    if mutations and st.button("Clear All Mutations"):
        st.session_state["mutations"] = []
        st.rerun()

    # ---- Cell type & Nuclease ----
    section_header("Experimental Parameters")

    col_ct, col_nuc = st.columns(2)
    with col_ct:
        cell_type = st.selectbox(
            "Cell type",
            list(CELL_TYPE_PARAMS.keys()),
            index=list(CELL_TYPE_PARAMS.keys()).index(st.session_state.get("cell_type", "iPSC")),
            help=(
                "The cell type being edited.  This affects baseline HDR efficiency, "
                "cell viability after DSBs, p53 activation, and safety scoring.  \n\n"
                "**iPSC**: Patient-derived induced pluripotent stem cells.  "
                "Strong p53 response; low viability after DSBs.  \n"
                "**HEK293T**: Immortalised human kidney cells.  p53 often "
                "inactivated; high HDR efficiency.  \n"
                "**K562**: Chronic myeloid leukaemia cells.  Good HDR.  \n"
                "**T cell**: Activated primary T cells (for CAR-T applications).  \n"
                "**HSC**: Haematopoietic stem cells.  Mostly quiescent; low HDR."
            ),
        )
        st.session_state["cell_type"] = cell_type

        ct_info = CELL_TYPE_PARAMS.get(cell_type, {})
        st.caption(ct_info.get("description", ""))

    with col_nuc:
        nuclease = st.selectbox(
            "Nuclease",
            list(NUCLEASE_PARAMS.keys()),
            index=list(NUCLEASE_PARAMS.keys()).index(st.session_state.get("nuclease", "enFnCas9")),
            help=(
                "The CRISPR nuclease for HDR-based editing steps.  \n\n"
                "**SpCas9**: Standard Cas9.  Blunt cuts.  Moderate specificity.  \n"
                "**enFnCas9**: Engineered FnCas9.  Staggered cuts (~3 bp 5' overhang).  "
                "Higher HDR, single-nucleobase specificity.  \n"
                "**Cas12a**: Also known as Cpf1.  5' overhangs (~5 bp).  TTTV PAM.  \n"
                "**vCas9**: Staggered-cut Cas9 variant.  ~6 bp 5' overhangs; "
                "best HDR enhancement."
            ),
        )
        st.session_state["nuclease"] = nuclease

        nuc_info = NUCLEASE_PARAMS.get(nuclease, {})
        st.caption(nuc_info.get("description", ""))

    # ---- Gene visualisation ----
    if gene is not None:
        section_header("Gene Structure Visualisation")
        try:
            if HAS_PLOTLY:
                fig = _gene_viz_plotly(gene, mutations)
                st.plotly_chart(fig, use_container_width=True)
            else:
                fig = _gene_viz_matplotlib(gene, mutations)
                st.pyplot(fig)
        except Exception as e:
            st.warning(f"Could not render gene diagram: {e}")


# ===================================================================
# PAGE: MOSAIC Strategy Optimizer
# ===================================================================

def page_mosaic():
    st.title("MOSAIC Strategy Optimizer")
    st.markdown(
        "Enumerate and rank all feasible multi-locus editing strategies for "
        "correcting the defined mutations."
    )

    gene = _build_gene_from_session()
    mutations = _get_mutations()

    if gene is None:
        st.warning("Please define the gene structure on the **Gene & Mutation Setup** page.")
        return
    if len(mutations) < 2:
        st.warning("At least 2 mutations are required for strategy analysis.  Add mutations on the Setup page.")
        return

    cell_type = st.session_state.get("cell_type", "iPSC")
    nuclease = st.session_state.get("nuclease", "SpCas9")

    st.info(
        f"**Gene:** {gene.gene_name} ({len(gene.exons)} exons)  |  "
        f"**Cell type:** {cell_type}  |  "
        f"**Nuclease:** {nuclease}  |  "
        f"**Mutations:** {len(mutations)}"
    )

    if st.button("Run MOSAIC Analysis", type="primary"):
        with st.spinner("Classifying mutations, enumerating strategies, scoring..."):
            try:
                classifier = MutationClassifier()
                enumerator = StrategyEnumerator()
                strategies = enumerator.enumerate_strategies(gene, mutations, cell_type, nuclease)
                scorer = StrategyScorer()
                ranked = scorer.rank_strategies(strategies, cell_type)
                reporter = StrategyReporter()
                report_text = reporter.generate_report(
                    gene, mutations, ranked, cell_type, nuclease
                )
                st.session_state["mosaic_results"] = {
                    "ranked": ranked,
                    "strategies": strategies,
                    "report_text": report_text,
                    "classifier": classifier,
                }
            except Exception as e:
                st.error(f"MOSAIC analysis failed: {e}")
                import traceback
                st.code(traceback.format_exc())
                return

    results = st.session_state.get("mosaic_results")
    if results is None:
        return

    ranked: List[ScoredStrategy] = results["ranked"]
    classifier: MutationClassifier = results["classifier"]

    # ---- Mutation classification table ----
    section_header("Mutation Classification")
    mut_rows = []
    for m in mutations:
        mut_rows.append({
            "Mutation": m.name,
            "Exon": m.exon_number,
            "Type": m.mutation_type,
            "Base-editable?": "Yes" if classifier.base_editing_amenable(m) else "No",
            "Prime-editable?": "Yes" if classifier.prime_editing_amenable(m) else "No",
            "HDR required?": "Yes" if classifier.hdr_required(m) else "No",
        })
    st.dataframe(pd.DataFrame(mut_rows), use_container_width=True, hide_index=True)

    # ---- Distance analysis ----
    section_header("Inter-Site Distance Analysis")
    if len(mutations) >= 2:
        ea, eb = mutations[0].exon_number, mutations[1].exon_number
        try:
            gd = gene.genomic_distance(ea, eb)
            ed = gene.exonic_distance(ea, eb)
            intd = gene.intronic_distance(ea, eb)
            c1, c2, c3 = st.columns(3)
            c1.metric("Genomic distance", format_bp(gd))
            c2.metric("Exonic distance", format_bp(ed))
            c3.metric("Intronic distance", format_bp(intd))

            can5 = gene.can_span_with_single_donor(ea, eb, 5000)
            can10 = gene.can_span_with_single_donor(ea, eb, 10000)
            if not can5:
                st.warning(
                    f"The two mutation sites are **{format_bp(gd)}** apart on the chromosome.  "
                    f"A single cssDNA donor (5 kb limit) **cannot** span both sites.  "
                    f"Each must be addressed independently."
                )
            elif not can10:
                st.info(
                    "A 5 kb cssDNA donor can span both sites. Verify donor synthesis feasibility."
                )
            else:
                st.success("Both mutation sites can be spanned by a single donor template.")
        except Exception as e:
            st.error(f"Distance calculation error: {e}")

    # ---- Strategy ranking table ----
    section_header("Strategy Ranking")

    rank_rows = []
    for ss in ranked:
        s = ss.strategy
        rank_rows.append({
            "Rank": ss.rank,
            "Strategy": s.name,
            "Overall": round(ss.overall_score, 3),
            "Safety": round(ss.safety_score, 3),
            "Efficiency": round(ss.efficiency_score, 3),
            "Time": round(ss.time_score, 2),
            "Cost": round(ss.cost_score, 2),
            "DSBs": s.num_dsbs,
            "Rounds": s.num_editing_rounds,
            "Screen": s.screening_burden,
        })
    st.dataframe(pd.DataFrame(rank_rows), use_container_width=True, hide_index=True)

    # ---- Bar chart comparing strategies ----
    section_header("Strategy Comparison")
    if HAS_PLOTLY:
        names = [ss.strategy.name for ss in ranked]
        fig = go.Figure()
        for dim, attr in [
            ("Safety", "safety_score"),
            ("Efficiency", "efficiency_score"),
            ("Time", "time_score"),
            ("Cost", "cost_score"),
        ]:
            fig.add_trace(go.Bar(
                name=dim,
                x=names,
                y=[getattr(ss, attr) for ss in ranked],
                marker_color=SCORE_COLOUR_MAP.get(dim.lower(), "#888"),
            ))
        fig.update_layout(
            barmode="group",
            yaxis_title="Score (0-1, higher = better)",
            height=400,
            legend=dict(orientation="h", y=1.12),
            margin=dict(t=40),
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        fig, ax = plt.subplots(figsize=(10, 5))
        names = [ss.strategy.name for ss in ranked]
        x = np.arange(len(names))
        w = 0.2
        for i, (dim, attr) in enumerate([
            ("Safety", "safety_score"),
            ("Efficiency", "efficiency_score"),
            ("Time", "time_score"),
            ("Cost", "cost_score"),
        ]):
            vals = [getattr(ss, attr) for ss in ranked]
            ax.bar(x + i * w, vals, w, label=dim)
        ax.set_xticks(x + 1.5 * w)
        ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("Score (0-1)")
        ax.legend()
        fig.tight_layout()
        st.pyplot(fig)

    # ---- Expandable details per strategy ----
    section_header("Strategy Details")
    for ss in ranked:
        s = ss.strategy
        css_cls = strategy_css_class(ss.rank, ss.safety_score)
        with st.expander(f"#{ss.rank}  {s.name}  (Overall: {ss.overall_score:.3f})"):
            st.markdown(s.description)
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("DSBs", s.num_dsbs)
            c2.metric("Editing rounds", s.num_editing_rounds)
            c3.metric("Donors needed", s.donor_templates_needed)
            c4.metric("Clones to screen", f"~{s.screening_burden}")

            if s.feasibility_warnings:
                st.warning("**Warnings:**\n" + "\n".join(f"- {w}" for w in s.feasibility_warnings))

    # ---- Download report ----
    section_header("Download Report")
    report_text = results.get("report_text", "")
    st.download_button(
        label="Download Full MOSAIC Report (.txt)",
        data=report_text,
        file_name=f"mosaic_report_{gene.gene_name}_{datetime.now():%Y%m%d}.txt",
        mime="text/plain",
    )


# ===================================================================
# PAGE: ConversionSim
# ===================================================================

def page_conversion():
    st.title("Conversion Tract Simulator")
    st.markdown(
        "Monte Carlo simulation of HDR gene conversion tracts.  Predicts "
        "the distribution of tract lengths and the probability that an "
        "edit at a given distance from the cut site will be incorporated."
    )

    gene = _build_gene_from_session()
    nuclease = st.session_state.get("nuclease", "SpCas9")
    nuc_params = NUCLEASE_PARAMS.get(nuclease, {})

    # ---- Inputs ----
    section_header("Simulation Parameters")

    col1, col2 = st.columns(2)
    with col1:
        cut_type = st.selectbox(
            "Cut type",
            ["blunt", "staggered_5prime"],
            index=0 if nuc_params.get("cut_type", "blunt") == "blunt" else 1,
            help=(
                "**Blunt**: Both DNA strands cut at the same position (SpCas9).  \n"
                "**Staggered 5'**: Non-target strand cut further from PAM, creating "
                "a 5' overhang.  Enhances HDR by facilitating end resection.  "
                "Used by Cas12a, enFnCas9, vCas9."
            ),
        )

        overhang = st.slider(
            "5' overhang length (bp)",
            min_value=0, max_value=10,
            value=nuc_params.get("stagger_bp", 0),
            help=(
                "Length of the 5' overhang in base pairs.  Must be 0 for blunt "
                "cuts.  Typical values: Cas12a = 5 bp, enFnCas9 = 3 bp, vCas9 = 6 bp.  "
                "Longer overhangs enhance HDR ~15% per bp (Chauhan et al., 2023)."
            ),
        )
        if cut_type == "blunt":
            overhang = 0

    with col2:
        donor_topology = st.selectbox(
            "Donor topology",
            ["circular_ssDNA", "linear_ssDNA", "linear_dsDNA", "plasmid"],
            index=0,
            help=(
                "Physical form of the donor template.  \n"
                "**circular_ssDNA (cssDNA)**: ~3x better than linear dsDNA.  "
                "Exonuclease-resistant.  Optimal for iPSCs.  \n"
                "**linear_ssDNA**: ~1.5x better than dsDNA.  \n"
                "**linear_dsDNA**: Standard PCR-amplified donor.  \n"
                "**plasmid**: Circular dsDNA.  Lower efficiency than ssDNA donors."
            ),
        )

        homology_arm = st.slider(
            "Homology arm length (bp)",
            min_value=50, max_value=1500, value=300, step=50,
            help=(
                "Length of each homology arm on the donor.  "
                "Optimal for cssDNA: 300 bp (Iyer et al., 2022).  "
                "For dsDNA donors: 800 bp is typical."
            ),
        )

    n_sims = st.slider(
        "Number of simulations",
        min_value=1000, max_value=100000, value=10000, step=1000,
        help="More simulations yield more precise statistics.  50,000+ recommended for publication.",
    )

    cell_type = st.session_state.get("cell_type", "iPSC")

    if st.button("Run Simulation", type="primary"):
        try:
            with st.spinner(f"Running {n_sims:,} Monte Carlo simulations..."):
                sim = ConversionSimulator(
                    cut_type=cut_type,
                    overhang_length=overhang,
                    donor_topology=donor_topology,
                    homology_arm_length=homology_arm,
                    cell_type=cell_type,
                    n_simulations=n_sims,
                )
                results = sim.run()
                stats = sim.summary()

                st.session_state["conversion_results"] = {
                    "sim": sim,
                    "results": results,
                    "stats": stats,
                }
        except Exception as e:
            st.error(f"Simulation failed: {e}")
            import traceback
            st.code(traceback.format_exc())
            return

    conv = st.session_state.get("conversion_results")
    if conv is None:
        return

    sim = conv["sim"]
    results = conv["results"]
    stats = conv["stats"]

    # ---- Summary statistics ----
    section_header("Summary Statistics")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("HDR success rate", format_pct(stats.get("hdr_success_rate", 0)))
    c2.metric("Mean tract length", f"{stats.get('tract_mean_bp', 0):.0f} bp")
    c3.metric("Median tract length", f"{stats.get('tract_median_bp', 0):.0f} bp")
    c4.metric("95th percentile", f"{stats.get('tract_p95_bp', 0):.0f} bp")

    # ---- Tract distribution histogram ----
    section_header("Tract Length Distribution")
    successful = results.tract_lengths_bp[results.hdr_success]

    if len(successful) > 0:
        if HAS_PLOTLY:
            fig = go.Figure()
            fig.add_trace(go.Histogram(
                x=successful,
                nbinsx=80,
                marker_color=PRIMARY,
                opacity=0.85,
                name="Tract lengths",
            ))
            median_val = float(np.median(successful))
            p95_val = float(np.percentile(successful, 95))
            for val, label, colour, dash in [
                (median_val, f"Median = {median_val:.0f} bp", WARNING, "dash"),
                (p95_val, f"95th pctl = {p95_val:.0f} bp", DANGER, "dot"),
                (homology_arm, f"Homology arm = {homology_arm} bp", "#888", "solid"),
            ]:
                fig.add_vline(x=val, line_dash=dash, line_color=colour,
                              annotation_text=label, annotation_position="top right")

            fig.update_layout(
                xaxis_title="Gene Conversion Tract Length (bp)",
                yaxis_title="Count",
                height=420,
                margin=dict(t=40),
                showlegend=False,
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.hist(successful, bins=80, color=PRIMARY, edgecolor="white", alpha=0.85)
            ax.axvline(np.median(successful), color=WARNING, ls="--", lw=2, label="Median")
            ax.axvline(np.percentile(successful, 95), color=DANGER, ls=":", lw=2, label="95th pctl")
            ax.axvline(homology_arm, color="#888", ls="-", lw=1.5, label="Homology arm")
            ax.set_xlabel("Tract length (bp)")
            ax.set_ylabel("Count")
            ax.legend()
            fig.tight_layout()
            st.pyplot(fig)
    else:
        st.warning("No HDR events occurred in the simulation.  Try increasing the number of simulations.")

    # ---- Survival curve ----
    section_header("Survival Curve (Probability of Reaching Distance X)")
    if len(successful) > 0:
        distances = np.arange(0, min(int(np.percentile(successful, 99.5)), 5000) + 1, 10)
        probs = [float(np.mean(successful >= d)) for d in distances]

        if HAS_PLOTLY:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=distances, y=probs,
                mode="lines", line=dict(color=PRIMARY, width=2.5),
                name="P(tract >= d)",
            ))
            fig.update_layout(
                xaxis_title="Distance from cut site (bp)",
                yaxis_title="P(tract >= distance)",
                height=380,
                margin=dict(t=20),
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(distances, probs, color=PRIMARY, lw=2)
            ax.set_xlabel("Distance from cut (bp)")
            ax.set_ylabel("P(tract >= distance)")
            ax.set_ylim(0, 1.05)
            fig.tight_layout()
            st.pyplot(fig)

    # ---- Probability at key distances ----
    section_header("Probability at Key Distances")
    prob_rows = []
    for d in [100, 200, 300, 500, 800, 1000, 1500, 2000, 3000, 5000]:
        key = f"p_conversion_at_{d}bp"
        p = stats.get(key, None)
        if p is None and len(successful) > 0:
            p = float(np.mean(successful >= d))
        elif p is None:
            p = 0.0
        prob_rows.append({"Distance (bp)": d, "P(tract >= d)": f"{p:.4f}", "P(%)": f"{p*100:.2f}%"})
    st.dataframe(pd.DataFrame(prob_rows), use_container_width=True, hide_index=True)

    # ---- Can HDR reach between exons? ----
    if gene is not None and len(st.session_state.get("mutations", [])) >= 2:
        section_header("Can HDR Reach Between Mutation Sites?")
        mutations = _get_mutations()
        ea, eb = mutations[0].exon_number, mutations[1].exon_number
        try:
            gd = gene.genomic_distance(ea, eb)
            p = sim.probability_at_distance(gd)
            if p < 1e-6:
                st.error(
                    f"HDR from exon {ea} **cannot** reach exon {eb}.  "
                    f"Genomic distance: {format_bp(gd)}.  "
                    f"Probability: essentially zero ({p:.2e})."
                )
            elif p < 0.01:
                st.warning(
                    f"Extremely unlikely.  P(reaching {format_bp(gd)}) = {p:.4f} ({p*100:.3f}%)."
                )
            else:
                st.success(
                    f"P(tract reaching exon {eb} from exon {ea}) = {p:.4f} ({p*100:.2f}%)."
                )
        except Exception:
            pass


# ===================================================================
# PAGE: ChromBridge
# ===================================================================

def page_chrombridge():
    st.title("3D Distance & Translocation Risk")
    st.markdown(
        "Predict the 3D spatial distance between two genomic loci using "
        "polymer physics, and assess whether a donor template can bridge "
        "the gap.  Also estimates translocation risk from simultaneous DSBs."
    )

    gene = _build_gene_from_session()
    mutations = _get_mutations()

    # ---- Inputs ----
    section_header("Input Parameters")

    # Auto-fill genomic distance if gene and mutations are available
    default_gd = 0
    if gene is not None and len(mutations) >= 2:
        try:
            default_gd = gene.genomic_distance(mutations[0].exon_number, mutations[1].exon_number)
        except Exception:
            pass

    genomic_distance = st.number_input(
        "Genomic distance between loci (bp)",
        min_value=0,
        value=default_gd,
        step=1000,
        help=(
            "The 1D distance along the chromosome between the two target loci.  "
            "Auto-filled from your gene and mutation setup if available.  "
            "You can override this with a manual value."
        ),
    )

    col1, col2 = st.columns(2)
    with col1:
        donor_size = st.slider(
            "Donor template size (kb)",
            min_value=1, max_value=20, value=3,
            help=(
                "Total size of the donor template in kilobases.  "
                "GATALYST cssDNA can be up to ~20 kb.  "
                "Typical cssDNA for single-site HDR: 1-5 kb."
            ),
        )

    with col2:
        donor_type_label = st.selectbox(
            "Donor type",
            ["cssDNA (circular ssDNA)", "lssDNA (linear ssDNA)", "dsDNA"],
            index=0,
            help=(
                "Physical form of the donor DNA.  Affects the coil size "
                "(ssDNA is more compact than dsDNA) and thus the bridging capacity."
            ),
        )
        # Map display labels to API-compatible values
        _donor_type_map = {
            "cssDNA (circular ssDNA)": "cssDNA",
            "lssDNA (linear ssDNA)": "lssDNA",
            "dsDNA": "dsDNA",
        }
        donor_type = _donor_type_map.get(donor_type_label, "cssDNA")

    cell_type = st.session_state.get("cell_type", "iPSC")

    if genomic_distance <= 0:
        st.warning("Enter a genomic distance > 0 to run the analysis.")
        return

    if st.button("Run ChromBridge Analysis", type="primary"):
        try:
            with st.spinner("Predicting 3D distance and bridgeability..."):
                predictor = ChromatinDistancePredictor()
                dist_result = predictor.predict_3d_distance(genomic_distance)
                bridge = predictor.can_donor_bridge(
                    genomic_distance_bp=genomic_distance,
                    donor_size_bp=donor_size * 1000,
                    donor_type=donor_type,
                )

                trans_pred = TranslocationRiskPredictor()
                # Use actual exon positions if available
                if gene is not None and len(mutations) >= 2:
                    s1 = mutations[0].position
                    s2 = mutations[1].position
                else:
                    s1 = 0
                    s2 = genomic_distance
                risk = trans_pred.estimate_risk(s1, s2, same_chromosome=True)
                safety = trans_pred.dual_dsb_safety_assessment(
                    s1, s2, cell_type=cell_type, same_chromosome=True
                )

                st.session_state["chrombridge_results"] = {
                    "dist_result": dist_result,
                    "bridge": bridge,
                    "risk": risk,
                    "safety": safety,
                }
        except Exception as e:
            st.error(f"ChromBridge analysis failed: {e}")
            import traceback
            st.code(traceback.format_exc())
            return

    cb = st.session_state.get("chrombridge_results")
    if cb is None:
        return

    dist_result = cb["dist_result"]
    bridge = cb["bridge"]
    risk = cb["risk"]
    safety = cb["safety"]

    # ---- 3D Distance prediction ----
    section_header("3D Distance Prediction")
    c1, c2, c3 = st.columns(3)
    c1.metric("Mean 3D distance", f"{dist_result.mean_3d_distance_nm:.0f} nm")
    c2.metric("5th percentile", f"{dist_result.p5_nm:.0f} nm")
    c3.metric("95th percentile", f"{dist_result.p95_nm:.0f} nm")

    st.markdown(
        f"The two loci are predicted to be **{dist_result.mean_3d_distance_nm:.0f} nm** "
        f"({dist_result.mean_3d_distance_um:.2f} um) apart in 3D nuclear space "
        f"(90% CI: {dist_result.p5_nm:.0f} - {dist_result.p95_nm:.0f} nm).  "
        f"Model: confined polymer, chromatin state: {dist_result.chromatin_state}."
    )

    # ---- Bridgeability analysis ----
    section_header("Donor Bridgeability")
    c1, c2, c3 = st.columns(3)
    c1.metric("Donor coil diameter", f"{bridge.donor_coil_diameter_nm:.0f} nm")
    c2.metric("Inter-locus distance", f"{bridge.inter_locus_distance_nm:.0f} nm")
    c3.metric("Bridgeability ratio", f"{bridge.bridgeability_ratio:.4f}")

    if bridge.feasible:
        st.success(
            f"A {donor_size} kb {donor_type.replace('_', ' ')} donor **may** be able "
            f"to bridge the gap (ratio = {bridge.bridgeability_ratio:.3f}).  "
            f"Verify experimentally."
        )
    else:
        ratio = bridge.bridgeability_ratio
        if ratio < 0.01:
            st.error(
                f"Bridging is **physically implausible**.  The donor coil "
                f"({bridge.donor_coil_diameter_nm:.0f} nm) is {1/ratio:.0f}x smaller "
                f"than the inter-locus distance ({bridge.inter_locus_distance_nm:.0f} nm)."
            )
        else:
            st.warning(
                f"Bridging is **unlikely** (ratio = {ratio:.4f}).  "
                f"Would require rare favourable chromatin conformations."
            )

    with st.expander("Detailed explanation"):
        st.text(bridge.explanation)

    # ---- Visual scale comparison ----
    section_header("Scale Diagram")
    if HAS_PLOTLY:
        fig = go.Figure()
        # Donor coil
        fig.add_shape(
            type="circle",
            x0=50 - bridge.donor_coil_diameter_nm / 2,
            x1=50 + bridge.donor_coil_diameter_nm / 2,
            y0=50 - bridge.donor_coil_diameter_nm / 2,
            y1=50 + bridge.donor_coil_diameter_nm / 2,
            fillcolor=PRIMARY, opacity=0.3,
            line=dict(color=PRIMARY, width=2),
        )
        fig.add_annotation(x=50, y=50, text=f"Donor coil<br>{bridge.donor_coil_diameter_nm:.0f} nm",
                           showarrow=False, font=dict(size=10))
        # Inter-locus distance (scaled)
        max_dim = max(bridge.inter_locus_distance_nm, bridge.donor_coil_diameter_nm) * 1.2
        fig.add_shape(
            type="line",
            x0=50 - bridge.inter_locus_distance_nm / 2,
            x1=50 + bridge.inter_locus_distance_nm / 2,
            y0=-30, y1=-30,
            line=dict(color=DANGER, width=3),
        )
        fig.add_annotation(
            x=50, y=-45,
            text=f"Inter-locus distance: {bridge.inter_locus_distance_nm:.0f} nm",
            showarrow=False, font=dict(size=10, color=DANGER),
        )
        fig.update_layout(
            height=250,
            xaxis=dict(visible=False, range=[50 - max_dim / 2, 50 + max_dim / 2]),
            yaxis=dict(visible=False, range=[-80, 50 + max_dim / 2]),
            showlegend=False,
            margin=dict(l=10, r=10, t=10, b=10),
        )
        st.plotly_chart(fig, use_container_width=True)

    # ---- Translocation risk ----
    section_header("Translocation Risk Assessment")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Deletion frequency", format_pct(risk.deletion_frequency, 2))
    c2.metric("Inversion frequency", format_pct(risk.inversion_frequency, 2))
    c3.metric("Translocation freq.", format_pct(risk.translocation_frequency, 3))
    c4.metric("Total rearrangement", format_pct(risk.total_rearrangement_risk, 2))

    risk_colours = {"low": SUCCESS, "moderate": WARNING, "high": DANGER, "very_high": DANGER}
    st.markdown(
        f"**Risk level:** <span style='color:{risk_colours.get(risk.risk_level, '#888')};font-weight:700;'>"
        f"{risk.risk_level.upper()}</span>",
        unsafe_allow_html=True,
    )

    # ---- Safety grade ----
    section_header("Overall Safety Grade")
    render_grade(safety.overall_grade, f"Safety Grade for dual DSBs in {cell_type}")
    st.markdown(safety.summary)

    with st.expander("Detailed safety breakdown"):
        st.markdown(f"- **Viability estimate:** {format_pct(safety.viability_estimate)}")
        st.markdown(f"- **p53 selection risk:** {safety.p53_selection_risk}")
        st.markdown(f"- **Chromothripsis risk:** {safety.chromothripsis_risk}")
        st.markdown("---")
        st.markdown("**Recommendations:**")
        st.markdown(safety.recommendations)


# ===================================================================
# PAGE: TopoPred (Donor Quality Check)
# ===================================================================

def page_topopred():
    st.title("Donor Quality Check (cssDNA-TopoPred)")
    st.markdown(
        "Analyse a cssDNA donor template for secondary structures that "
        "could impair HDR efficiency: G-quadruplexes, hairpins, and "
        "per-nucleotide accessibility."
    )

    # ---- Auto-generate donor from Ensembl ----
    gene_info_obj = st.session_state.get("ensembl_gene_info")
    gene_obj = st.session_state.get("gene")
    mutations = st.session_state.get("mutations", [])

    if gene_info_obj and mutations:
        section_header("Auto-Generate Donor from Ensembl")
        st.markdown(
            "Your gene was fetched from Ensembl. Select a mutation to auto-design "
            "a cssDNA donor with correct genomic homology arms."
        )
        mut_options = [f"Exon {m.exon_number}: {m.ref_allele}>{m.alt_allele}" for m in mutations]
        selected_mut = st.selectbox("Select mutation for donor design", mut_options)
        arm_len = st.slider("Homology arm length (bp)", 100, 1000, 300, 50,
                           help="Optimal for cssDNA: 300 bp (Iyer et al., 2022)")

        if st.button("Auto-Design cssDNA Donor"):
            mut_idx = mut_options.index(selected_mut)
            mut = mutations[mut_idx]
            try:
                from crisprarchitect.utils.ensembl import design_cssdna_donor
                with st.spinner(f"Fetching genomic sequence for exon {mut.exon_number}..."):
                    donor = design_cssdna_donor(
                        gene_info_obj, mut.exon_number,
                        ref_allele=mut.ref_allele, alt_allele=mut.alt_allele,
                        homology_arm_length=arm_len
                    )
                st.session_state["_auto_donor"] = donor
                st.success(
                    f"Donor designed: {donor['total_length']} bp total "
                    f"({arm_len} bp left arm + {len(donor['exon_corrected'])} bp exon + "
                    f"{arm_len} bp right arm)"
                )
                for note in donor["design_notes"]:
                    st.caption(note)
            except Exception as e:
                st.error(f"Error designing donor: {e}")

        auto_donor = st.session_state.get("_auto_donor")
        if auto_donor:
            st.info(
                f"Auto-designed donor loaded ({auto_donor['total_length']} bp). "
                "It will be used for the analysis below."
            )
        st.divider()

    # ---- Inputs ----
    section_header("Donor Sequence")

    # Pre-fill from auto donor if available
    _auto = st.session_state.get("_auto_donor")
    _default_seq = _auto["donor_sequence"] if _auto else ""

    donor_seq = st.text_area(
        "Paste your cssDNA donor sequence (5' -> 3')",
        value=_default_seq,
        height=200,
        help=(
            "The full circular single-stranded DNA sequence of your donor template.  "
            "Include the left homology arm, insert/correction region, right homology "
            "arm, and any backbone sequence (e.g., phagemid origin).  "
            "Use A, T, C, G only.  Whitespace and line breaks are stripped."
        ),
        placeholder="ATCGATCGATCG...paste your full donor sequence here...",
    )

    if not donor_seq.strip():
        st.info("Paste a donor sequence above, or use the demo button below.")

        if st.button("Load Demo Donor"):
            random.seed(123)
            def _rand_dna(n):
                return "".join(random.choice("ATCG") for _ in range(n))
            demo = _rand_dna(300) + _rand_dna(150) + _rand_dna(300) + _rand_dna(2200)
            st.session_state["_demo_donor"] = demo
            st.rerun()

        demo_seq = st.session_state.get("_demo_donor", "")
        if demo_seq:
            donor_seq = demo_seq
            st.text_area("Demo sequence loaded:", value=demo_seq[:200] + "...", height=80, disabled=True)

    # Clean sequence
    clean_seq = "".join(c for c in donor_seq.upper() if c in "ATCG")
    if clean_seq:
        st.caption(f"Cleaned sequence length: {len(clean_seq)} nt")

    section_header("Region Definitions")
    st.markdown(
        "Define the left and right homology arm boundaries (0-based positions). "
        "Optionally mark coding regions for synonymous codon optimisation."
    )

    # Auto-fill arm positions from auto-donor if available
    _ad = st.session_state.get("_auto_donor")
    _default_la_end = _ad["homology_arm_length"] if _ad else min(300, len(clean_seq))
    _default_ra_start = (_ad["homology_arm_length"] + len(_ad["exon_corrected"])) if _ad else min(450, len(clean_seq))
    _default_ra_end = _ad["total_length"] if _ad else min(750, len(clean_seq))

    col1, col2 = st.columns(2)
    with col1:
        la_start = st.number_input("Left arm start", value=0, min_value=0,
                                   help="Start position of the left homology arm (0-based).")
        la_end = st.number_input("Left arm end", value=_default_la_end,
                                 min_value=0,
                                 help="End position of the left homology arm (exclusive).")
    with col2:
        ra_start = st.number_input("Right arm start", value=_default_ra_start,
                                   min_value=0,
                                   help="Start position of the right homology arm.")
        ra_end = st.number_input("Right arm end", value=min(750, len(clean_seq)),
                                 min_value=0,
                                 help="End position of the right homology arm (exclusive).")

    coding_text = st.text_input(
        "Coding regions (optional)",
        value="",
        help=(
            "Comma-separated start-end pairs for coding regions within the donor.  "
            "Example: 300-450, 500-600.  Coding regions are used for synonymous "
            "codon optimisation if structures are found."
        ),
    )
    coding_regions = []
    if coding_text.strip():
        for part in coding_text.split(","):
            part = part.strip()
            if "-" in part:
                try:
                    s, e = part.split("-")
                    coding_regions.append((int(s.strip()), int(e.strip())))
                except ValueError:
                    st.warning(f"Invalid coding region: {part}")

    if not clean_seq:
        return

    if st.button("Analyse Donor", type="primary"):
        try:
            with st.spinner("Scanning for G-quadruplexes, hairpins, and computing accessibility..."):
                analyzer = DonorAnalyzer()
                report = analyzer.analyze(
                    sequence=clean_seq,
                    left_arm=(la_start, la_end),
                    right_arm=(ra_start, ra_end),
                    coding_regions=coding_regions if coding_regions else None,
                )
                st.session_state["topopred_results"] = report
        except Exception as e:
            st.error(f"Analysis failed: {e}")
            import traceback
            st.code(traceback.format_exc())
            return

    report = st.session_state.get("topopred_results")
    if report is None:
        return

    # ---- G-quadruplex results ----
    section_header("G-Quadruplex Analysis")
    g4s = report.get("g4_motifs", [])
    g4_risk = report.get("g4_risk_score", 0)

    c1, c2 = st.columns(2)
    c1.metric("G4 motifs found", len(g4s))
    c2.metric("G4 risk score", f"{g4_risk:.2f}")

    if g4s:
        g4_rows = []
        for m in g4s:
            g4_rows.append({
                "Position": f"{m.start}-{m.end}",
                "Strand": m.strand,
                "Tetrads": m.num_tetrads,
                "Loops": str(m.loop_lengths),
                "Stability": f"{m.stability_score:.2f}",
            })
        st.dataframe(pd.DataFrame(g4_rows), use_container_width=True, hide_index=True)
    else:
        st.success("No G-quadruplex motifs detected.  Good -- G4s are the most disruptive structures for HDR.")

    # ---- Hairpin results ----
    section_header("Hairpin Analysis")
    hairpins = report.get("hairpins", [])
    st.metric("Stable hairpins found", len(hairpins))

    if hairpins:
        hp_rows = []
        for hp in hairpins:
            hp_rows.append({
                "Position": f"{hp.stem_start}-{hp.stem_end}",
                "Stem (bp)": hp.stem_length,
                "Loop (nt)": hp.loop_length,
                "dG (kcal/mol)": f"{hp.free_energy_kcal_mol:.1f}",
                "GC content": f"{hp.gc_content:.0%}",
            })
        st.dataframe(pd.DataFrame(hp_rows), use_container_width=True, hide_index=True)
    else:
        st.success("No significant hairpins detected in the donor.")

    # ---- Accessibility ----
    section_header("Homology Arm Accessibility")
    acc = report.get("accessibility_report")
    if acc is not None:
        c1, c2, c3 = st.columns(3)
        c1.metric("Left arm accessibility", format_pct(acc.left_arm_accessibility))
        c2.metric("Right arm accessibility", format_pct(acc.right_arm_accessibility))
        c3.metric("Overall score", format_pct(acc.overall_score))

        if acc.warnings:
            for w in acc.warnings:
                st.warning(w)

        st.markdown(f"**Recommendation:** {acc.recommendation}")

    # ---- Accessibility heatmap ----
    acc_map = report.get("accessibility_map")
    if acc_map is not None and len(acc_map) > 0:
        section_header("Per-Nucleotide Accessibility Map")

        if HAS_PLOTLY:
            # Show as a heatmap-like bar/line chart
            positions = np.arange(len(acc_map))
            colours = [SUCCESS if v > 0.6 else (WARNING if v > 0.3 else DANGER) for v in acc_map]

            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=positions, y=acc_map,
                marker_color=colours,
                name="Accessibility",
            ))
            # Mark arm boundaries
            for pos, label in [(la_start, "LA start"), (la_end, "LA end"),
                               (ra_start, "RA start"), (ra_end, "RA end")]:
                fig.add_vline(x=pos, line_dash="dash", line_color="#888",
                              annotation_text=label, annotation_position="top")
            fig.update_layout(
                xaxis_title="Position (nt)",
                yaxis_title="Accessibility (0=structured, 1=free)",
                height=300,
                margin=dict(t=30),
                showlegend=False,
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            fig, ax = plt.subplots(figsize=(12, 2.5))
            ax.bar(range(len(acc_map)), acc_map, width=1.0, color=PRIMARY, alpha=0.7)
            ax.set_xlabel("Position (nt)")
            ax.set_ylabel("Accessibility")
            ax.set_ylim(0, 1.05)
            fig.tight_layout()
            st.pyplot(fig)

    # ---- Optimisation suggestions ----
    opt = report.get("optimization")
    if opt is not None and opt.changes_made:
        section_header("Optimisation Suggestions")
        st.markdown(
            f"**{len(opt.changes_made)} synonymous codon changes** suggested to "
            f"disrupt problematic secondary structures."
        )
        opt_rows = []
        for ch in opt.changes_made:
            opt_rows.append({
                "Position": ch.position,
                "Original codon": ch.original_codon,
                "Suggested codon": ch.suggested_codon,
                "Amino acid": ch.amino_acid,
                "Target structure": ch.target_structure,
            })
        st.dataframe(pd.DataFrame(opt_rows), use_container_width=True, hide_index=True)

        if hasattr(opt, "optimized_sequence") and opt.optimized_sequence:
            st.download_button(
                label="Download Optimised Sequence",
                data=opt.optimized_sequence,
                file_name="optimised_donor.txt",
                mime="text/plain",
            )

    # ---- Full summary download ----
    section_header("Download Full Report")
    summary_text = report.get("summary", "")
    st.download_button(
        label="Download TopoPred Report (.txt)",
        data=summary_text,
        file_name="topopred_report.txt",
        mime="text/plain",
    )


# ===================================================================
# PAGE: Full Analysis Report
# ===================================================================

def page_report():
    st.title("Full Analysis Report")
    st.markdown(
        "Comprehensive summary combining results from all CRISPRArchitect "
        "modules.  Run individual modules first to populate this report."
    )

    gene = _build_gene_from_session()
    mutations = _get_mutations()

    if gene is None:
        st.warning("Define your gene on the Setup page to generate a report.")
        return

    # ---- Gene & Mutation Summary ----
    section_header("1. Gene & Mutation Summary")
    st.code(gene.summary(), language="text")

    if mutations:
        classifier = MutationClassifier()
        for i, m in enumerate(mutations, 1):
            st.markdown(f"**Mutation {i}:** {m.name}")
            st.text(classifier.describe_mutation(m))

    # ---- MOSAIC results ----
    mosaic = st.session_state.get("mosaic_results")
    if mosaic:
        section_header("2. Strategy Analysis (MOSAIC)")
        ranked = mosaic["ranked"]
        if ranked:
            best = ranked[0]
            st.success(
                f"**Recommended strategy:** {best.strategy.name} "
                f"(Overall: {best.overall_score:.3f})"
            )
            st.markdown(best.strategy.description)

            rank_rows = []
            for ss in ranked:
                rank_rows.append({
                    "Rank": ss.rank,
                    "Strategy": ss.strategy.name,
                    "Overall": f"{ss.overall_score:.3f}",
                    "Safety": f"{ss.safety_score:.3f}",
                    "Efficiency": f"{ss.efficiency_score:.3f}",
                })
            st.dataframe(pd.DataFrame(rank_rows), use_container_width=True, hide_index=True)
    else:
        st.info("Run the **Strategy Optimizer (MOSAIC)** module to see results here.")

    # ---- ConversionSim results ----
    conv = st.session_state.get("conversion_results")
    if conv:
        section_header("3. Conversion Tract Analysis")
        stats = conv["stats"]
        c1, c2, c3 = st.columns(3)
        c1.metric("HDR success rate", format_pct(stats.get("hdr_success_rate", 0)))
        c2.metric("Median tract", f"{stats.get('tract_median_bp', 0):.0f} bp")
        c3.metric("95th percentile", f"{stats.get('tract_p95_bp', 0):.0f} bp")
    else:
        st.info("Run the **Conversion Tract Simulator** to see results here.")

    # ---- ChromBridge results ----
    cb = st.session_state.get("chrombridge_results")
    if cb:
        section_header("4. 3D Distance & Risk (ChromBridge)")
        dist = cb["dist_result"]
        bridge = cb["bridge"]
        safety = cb["safety"]

        c1, c2, c3 = st.columns(3)
        c1.metric("Mean 3D distance", f"{dist.mean_3d_distance_nm:.0f} nm")
        c2.metric("Bridgeable?", "Yes" if bridge.feasible else "No")
        c3.metric("Safety grade", safety.overall_grade)

        risk = cb["risk"]
        st.markdown(
            f"Translocation risk: **{risk.risk_level.upper()}** "
            f"(total rearrangement: {format_pct(risk.total_rearrangement_risk, 2)})"
        )
    else:
        st.info("Run the **3D Distance & Risk (ChromBridge)** module to see results here.")

    # ---- TopoPred results ----
    topo = st.session_state.get("topopred_results")
    if topo:
        section_header("5. Donor Quality (TopoPred)")
        acc = topo.get("accessibility_report")
        if acc:
            c1, c2, c3 = st.columns(3)
            c1.metric("Left arm accessibility", format_pct(acc.left_arm_accessibility))
            c2.metric("Right arm accessibility", format_pct(acc.right_arm_accessibility))
            c3.metric("Recommendation", acc.recommendation)
    else:
        st.info("Run the **Donor Quality Check (TopoPred)** module to see results here.")

    # ---- Key findings ----
    section_header("6. Key Findings & Recommendations")

    findings = []
    if mosaic and mosaic["ranked"]:
        best = mosaic["ranked"][0]
        findings.append(
            f"**Best strategy:** {best.strategy.name} "
            f"(overall score {best.overall_score:.3f}, safety {best.safety_score:.3f})"
        )
        if best.strategy.num_dsbs == 0:
            findings.append(
                "This strategy introduces **no double-strand breaks**, making it "
                "the safest option especially for iPSCs."
            )

    if conv:
        stats = conv["stats"]
        findings.append(
            f"**HDR tract reach:** Median {stats.get('tract_median_bp', 0):.0f} bp, "
            f"95th percentile {stats.get('tract_p95_bp', 0):.0f} bp."
        )

    if cb:
        bridge = cb["bridge"]
        if not bridge.feasible:
            findings.append(
                "**Donor bridging:** The donor template **cannot** physically "
                "bridge the two mutation sites in 3D space."
            )

    if findings:
        for f in findings:
            st.markdown(f"- {f}")
    else:
        st.markdown("Run the analysis modules above to generate findings.")

    # ---- Download combined report ----
    section_header("Download Combined Report")
    report_parts = []
    report_parts.append("=" * 78)
    report_parts.append("CRISPRArchitect — Full Analysis Report")
    report_parts.append(f"Generated: {datetime.now():%Y-%m-%d %H:%M}")
    report_parts.append("=" * 78)
    report_parts.append("")

    if gene:
        report_parts.append(gene.summary())
        report_parts.append("")

    if mutations:
        report_parts.append("MUTATIONS:")
        classifier = MutationClassifier()
        for m in mutations:
            report_parts.append(classifier.describe_mutation(m))
            report_parts.append("")

    if mosaic:
        report_parts.append("")
        report_parts.append(mosaic.get("report_text", ""))

    if conv:
        report_parts.append("")
        report_parts.append("CONVERSION TRACT SUMMARY:")
        stats = conv["stats"]
        for k, v in stats.items():
            report_parts.append(f"  {k}: {v}")

    if cb:
        report_parts.append("")
        report_parts.append("CHROMBRIDGE SUMMARY:")
        dist = cb["dist_result"]
        report_parts.append(f"  Mean 3D distance: {dist.mean_3d_distance_nm:.0f} nm")
        report_parts.append(f"  Bridgeable: {cb['bridge'].feasible}")
        report_parts.append(f"  Safety grade: {cb['safety'].overall_grade}")
        report_parts.append(f"  {cb['safety'].summary}")

    if topo:
        report_parts.append("")
        report_parts.append("TOPOPRED SUMMARY:")
        report_parts.append(topo.get("summary", ""))

    full_report = "\n".join(report_parts)

    st.download_button(
        label="Download Full Report (.txt)",
        data=full_report,
        file_name=f"crisprarchitect_report_{datetime.now():%Y%m%d}.txt",
        mime="text/plain",
    )


# ===================================================================
# Main router
# ===================================================================

PAGE_FUNCTIONS = {
    "home": page_home,
    "setup": page_setup,
    "mosaic": page_mosaic,
    "conversion": page_conversion,
    "chrombridge": page_chrombridge,
    "topopred": page_topopred,
    "report": page_report,
}

PAGE_FUNCTIONS[page]()
