#!/usr/bin/env python3
"""
Generate CRISPRArchitect v3 Pipeline Flowchart
===============================================
Creates a professional flowchart of the end-to-end pipeline using Graphviz.
Outputs PNG and PDF for use in README, documentation, and PPT.

Usage:
    python3 generate_pipeline_flowchart.py

Output:
    figures/v3/Fig_Pipeline_Flowchart.png
    figures/v3/Fig_Pipeline_Flowchart.pdf
"""

import subprocess
import os
import sys

# Output paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, "figures", "v3")
os.makedirs(FIG_DIR, exist_ok=True)

DOT_PATH = os.path.join(FIG_DIR, "Fig_Pipeline_Flowchart.dot")

# Color palette matching the PPT dark theme
NAVY = "#1B2838"
DARKER = "#131C28"
CARD = "#233448"
TEAL = "#1B9E77"
CORAL = "#D95F02"
GOLD = "#FFC107"
WHITE = "#FFFFFF"
LIGHT_GRAY = "#CCCCCC"
SOFT_WHITE = "#E8EEF4"
GREEN = "#2ECC71"
BLUE = "#3498DB"
RED = "#E74C3C"

dot_source = f"""
digraph CRISPRArchitect_Pipeline {{
    // Global settings
    graph [
        rankdir=TB
        bgcolor="{NAVY}"
        fontname="Helvetica"
        fontsize=14
        fontcolor="{WHITE}"
        pad="0.5"
        nodesep=0.5
        ranksep=0.6
        dpi=200
        size="10,18"
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="24" COLOR="{WHITE}">CRISPRArchitect v3 Pipeline</FONT></TD></TR>
                <TR><TD><FONT FACE="Helvetica" POINT-SIZE="14" COLOR="{TEAL}">End-to-End Genome Editing Strategy Design</FONT></TD></TR>
            </TABLE>
        >
        labelloc=t
        labeljust=c
    ]

    node [
        fontname="Helvetica"
        fontsize=11
        style="filled,rounded"
        shape=box
        penwidth=1.5
        margin="0.15,0.10"
    ]

    edge [
        color="{TEAL}"
        penwidth=1.5
        arrowsize=0.8
        fontname="Helvetica"
        fontsize=9
        fontcolor="{LIGHT_GRAY}"
    ]

    // ── INPUT ──
    input [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="13" COLOR="{NAVY}">INPUT</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="10" COLOR="{NAVY}">Genomic coordinates / HGVS / ClinVar batch</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{GOLD}"
        color="{GOLD}"
        fontcolor="{NAVY}"
    ]

    // ── STAGE 1: Variant Parsing ──
    s1 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">1. Variant Parsing</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">HGVS parser · ClinVar TSV/VCF · coordinate input</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 2: Transcript Fetch ──
    s2 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">2. Transcript Fetch</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Ensembl REST API · GRCh38 · canonical transcript</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 3: Coordinate Mapping ──
    s3 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">3. Coordinate Mapping</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Genomic → CDS · exon/codon position · strand handling</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 4: Reference Validation + Annotation ──
    s4 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">4. Validation &amp; Annotation</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Reference allele check · consequence classification (ACMG)</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Missense · nonsense · splice · synonymous</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 5: Variant Normalization ──
    s5 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">5. Variant Normalization</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">NormalizedVariant object · local sequence ±200 bp</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 6: Multi-Nuclease PAM Scan ──
    s6 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{CORAL}">6. Multi-Nuclease PAM Scan</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">5 nucleases: SpCas9 · enFnCas9 · SpCas9-NG · SpRY · Cas12a</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">GC filter · poly-T filter · composite guide scoring</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{CORAL}"
    ]

    // ── STAGE 7: Feasibility Assessment (parallel branches) ──
    subgraph cluster_feasibility {{
        graph [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="13" COLOR="{CORAL}">7. Multi-Modality Feasibility Assessment</FONT></TD></TR>
                </TABLE>
            >
            style="rounded,dashed"
            color="{CORAL}"
            bgcolor="{DARKER}"
            fontcolor="{CORAL}"
            penwidth=1.5
        ]

        be [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{TEAL}">Base Editing</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">9 editor profiles (3 Tier A + 6 Tier B)</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">ABE7.10 · ABE8e · BE4max + fusions</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Window check · bystander identification</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{TEAL}"
        ]

        pe [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{CORAL}">Prime Editing</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">pegRNA design (PBS + RT template)</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">PE3 nicking guide search</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{CORAL}"
        ]

        hdr [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{GOLD}">HDR Design</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Donor type: dsDNA · lssDNA · cssDNA · AAV</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">ConversionSim Monte Carlo (10K sims)</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Cut-to-edit distance · HA optimization</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{GOLD}"
        ]
    }}

    // ── STAGE 8: Strategy Generation ──
    s8 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{WHITE}">8. Strategy Generation</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Single-step · dual · hybrid combinations</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">BE+PE · BE+HDR · PE+HDR · dual-BE · dual-PE</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Rejected strategies tagged with reasons</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{TEAL}"
    ]

    // ── STAGE 9: TOPSIS Ranking ──
    s9 [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="12" COLOR="{CORAL}">9. TOPSIS 6D Multi-Criteria Ranking</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Safety (0.28) · Feasibility (0.23) · Complexity (0.19)</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Risk (0.14) · Confidence (0.09) · Consequence (0.07)</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Hwang &amp; Yoon, 1981</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{CARD}"
        color="{CORAL}"
    ]

    // ── STAGE 10: Pareto + Sensitivity ──
    subgraph cluster_validation {{
        graph [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="13" COLOR="{TEAL}">10. Robustness Validation</FONT></TD></TR>
                </TABLE>
            >
            style="rounded,dashed"
            color="{TEAL}"
            bgcolor="{DARKER}"
            fontcolor="{TEAL}"
            penwidth=1.5
        ]

        pareto [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{TEAL}">Pareto Front</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">Weight-independent</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">dominance analysis</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{TEAL}"
        ]

        sensitivity [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{CORAL}">Monte Carlo Sensitivity</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">10,000 Dirichlet weight</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">permutations · rank stability</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{CORAL}"
        ]

        crossmethod [
            label=<
                <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                    <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="11" COLOR="{GOLD}">Cross-Method</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">TOPSIS · VIKOR · WPM</FONT></TD></TR>
                    <TR><TD><FONT POINT-SIZE="9" COLOR="{LIGHT_GRAY}">100% concordance</FONT></TD></TR>
                </TABLE>
            >
            fillcolor="{CARD}"
            color="{GOLD}"
        ]
    }}

    // ── OUTPUT ──
    output [
        label=<
            <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">
                <TR><TD><FONT FACE="Helvetica-Bold" POINT-SIZE="13" COLOR="{NAVY}">OUTPUT: Ranked Strategy Report</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="10" COLOR="{NAVY}">TOPSIS scores · rank stability · Pareto status</FONT></TD></TR>
                <TR><TD><FONT POINT-SIZE="10" COLOR="{NAVY}">rejection reasons · parameter provenance · CIs/SEs</FONT></TD></TR>
            </TABLE>
        >
        fillcolor="{GREEN}"
        color="{GREEN}"
        fontcolor="{NAVY}"
    ]

    // ── EDGES ──
    input -> s1
    s1 -> s2
    s2 -> s3
    s3 -> s4
    s4 -> s5
    s5 -> s6

    s6 -> be [label="  guide\\ncandidates"]
    s6 -> pe [label="  guide\\ncandidates"]
    s6 -> hdr [label="  guide\\ncandidates"]

    be -> s8
    pe -> s8
    hdr -> s8

    s8 -> s9
    s9 -> pareto
    s9 -> sensitivity
    s9 -> crossmethod

    pareto -> output
    sensitivity -> output
    crossmethod -> output
}}
"""

# Write DOT file
with open(DOT_PATH, "w") as f:
    f.write(dot_source)

# Generate PNG and PDF
png_path = os.path.join(FIG_DIR, "Fig_Pipeline_Flowchart.png")
pdf_path = os.path.join(FIG_DIR, "Fig_Pipeline_Flowchart.pdf")
svg_path = os.path.join(FIG_DIR, "Fig_Pipeline_Flowchart.svg")

for fmt, path in [("png", png_path), ("pdf", pdf_path), ("svg", svg_path)]:
    cmd = ["dot", f"-T{fmt}", DOT_PATH, "-o", path]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Generated: {path}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR generating {fmt}: {e.stderr}", file=sys.stderr)
    except FileNotFoundError:
        print("ERROR: Graphviz 'dot' command not found. Install with: brew install graphviz",
              file=sys.stderr)
        sys.exit(1)

print("\nDone! Flowchart generated in PNG, PDF, and SVG formats.")
