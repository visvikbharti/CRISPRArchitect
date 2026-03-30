#!/usr/bin/env python3
"""
Generate CRISPRArchitect v3 Lab Meeting Presentation
=====================================================
Professional 25-slide PPTX for presenting CRISPRArchitect v3
to PI Debojyoti Chakraborty and labmates at CSIR-IGIB.

Design: Dark navy background, teal accents, coral/orange for key findings.
All data is real and verified — no placeholders or dummy content.

Usage:
    python3 generate_presentation_v3.py

Output:
    CRISPRArchitect_v3_LabMeeting.pptx
"""

try:
    from pptx import Presentation
    from pptx.util import Inches, Pt, Emu
    from pptx.dml.color import RGBColor
    from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
    from pptx.enum.shapes import MSO_SHAPE
except ImportError:
    import subprocess
    subprocess.check_call(["pip", "install", "python-pptx"])
    from pptx import Presentation
    from pptx.util import Inches, Pt, Emu
    from pptx.dml.color import RGBColor
    from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
    from pptx.enum.shapes import MSO_SHAPE

import os

# ===========================================================================
# Color palette (dark professional theme)
# ===========================================================================
BG_NAVY    = RGBColor(0x1B, 0x28, 0x38)
BG_DARKER  = RGBColor(0x13, 0x1C, 0x28)
BG_CARD    = RGBColor(0x23, 0x34, 0x48)
WHITE      = RGBColor(0xFF, 0xFF, 0xFF)
TEAL       = RGBColor(0x1B, 0x9E, 0x77)
CORAL      = RGBColor(0xD9, 0x5F, 0x02)  # The orange header you liked
GOLD       = RGBColor(0xFF, 0xC1, 0x07)
LIGHT_GRAY = RGBColor(0xCC, 0xCC, 0xCC)
DARK_GRAY  = RGBColor(0x88, 0x99, 0xAA)
SOFT_WHITE = RGBColor(0xE8, 0xEE, 0xF4)
RED_ACCENT = RGBColor(0xE7, 0x4C, 0x3C)
GREEN_ACC  = RGBColor(0x2E, 0xCC, 0x71)
BLUE_ACC   = RGBColor(0x34, 0x98, 0xDB)

# Slide dimensions (16:9)
SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)
TOTAL_SLIDES = 26


# ===========================================================================
# Helper functions
# ===========================================================================
def set_slide_bg(slide, color=BG_NAVY):
    background = slide.background
    fill = background.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_textbox(slide, left, top, width, height, text, font_size=18,
                color=WHITE, bold=False, alignment=PP_ALIGN.LEFT,
                font_name="Calibri"):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.color.rgb = color
    p.font.bold = bold
    p.font.name = font_name
    p.alignment = alignment
    return txBox


def add_multiline(slide, left, top, width, height, lines, font_size=16,
                  color=WHITE, line_spacing=1.3, font_name="Calibri"):
    """Add multiple lines with formatting tuples: (text, color, bold)."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, line_data in enumerate(lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.line_spacing = line_spacing
        if isinstance(line_data, tuple):
            text, lcolor, lbold = line_data
        else:
            text, lcolor, lbold = line_data, color, False
        run = p.add_run()
        run.text = text
        run.font.size = Pt(font_size)
        run.font.color.rgb = lcolor
        run.font.name = font_name
        run.font.bold = lbold
    return txBox


def add_bullet_list(slide, left, top, width, height, items, font_size=17,
                    color=WHITE, bullet_color=TEAL, line_spacing=1.4,
                    font_name="Calibri"):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        if isinstance(item, tuple):
            text, item_color, item_bold = item
        else:
            text, item_color, item_bold = item, color, False
        p.line_spacing = line_spacing
        p.space_after = Pt(4)
        run_b = p.add_run()
        run_b.text = "  \u25B8  "
        run_b.font.size = Pt(font_size)
        run_b.font.color.rgb = bullet_color
        run_b.font.name = font_name
        run_b.font.bold = True
        run_t = p.add_run()
        run_t.text = text
        run_t.font.size = Pt(font_size)
        run_t.font.color.rgb = item_color
        run_t.font.name = font_name
        run_t.font.bold = item_bold
    return txBox


def add_accent_line(slide, left, top, width, color=TEAL, height=Inches(0.04)):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_card(slide, left, top, width, height, fill_color=BG_CARD):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    shape.line.fill.background()
    return shape


def add_slide_number(slide, num):
    add_textbox(slide, Inches(12.0), Inches(7.05), Inches(1.2), Inches(0.4),
                f"{num}/{TOTAL_SLIDES}", font_size=10, color=DARK_GRAY,
                alignment=PP_ALIGN.RIGHT)


def add_title_bar(slide, title_text, subtitle_text=None):
    add_textbox(slide, Inches(0.7), Inches(0.3), Inches(11.5), Inches(0.7),
                title_text, font_size=30, color=WHITE, bold=True)
    add_accent_line(slide, Inches(0.7), Inches(0.95), Inches(2.5), TEAL)
    if subtitle_text:
        add_textbox(slide, Inches(0.7), Inches(1.05), Inches(11.5), Inches(0.5),
                    subtitle_text, font_size=14, color=DARK_GRAY)


def add_coral_banner(slide, text="KEY FINDING"):
    """The orange/coral full-width banner at the top."""
    add_card(slide, Inches(0.0), Inches(0.0), SLIDE_W, Inches(0.5), CORAL)
    add_textbox(slide, Inches(0.0), Inches(0.05), SLIDE_W, Inches(0.4),
                text, font_size=18, color=WHITE, bold=True,
                alignment=PP_ALIGN.CENTER)


def add_stat_box(slide, left, top, width, height, value, label,
                 value_color=TEAL, bg_color=BG_CARD):
    add_card(slide, left, top, width, height, bg_color)
    add_textbox(slide, left, top + Inches(0.15), width, Inches(0.6),
                value, font_size=32, color=value_color, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, left, top + Inches(0.7), width, Inches(0.4),
                label, font_size=12, color=LIGHT_GRAY,
                alignment=PP_ALIGN.CENTER)


# ===========================================================================
# Figure paths (all generated from real data)
# ===========================================================================
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'figures', 'v3')

def get_fig(name):
    """Return path to a figure PNG if it exists, else None."""
    path = os.path.join(FIG_DIR, name)
    if os.path.exists(path):
        return path
    print(f"  WARNING: Figure not found: {path}")
    return None


def add_figure(slide, fig_name, left, top, width=None, height=None):
    """Add a figure image to a slide if it exists."""
    path = get_fig(fig_name)
    if path:
        kwargs = {}
        if width:
            kwargs['width'] = width
        if height:
            kwargs['height'] = height
        if not kwargs:
            kwargs['width'] = Inches(5.5)
        slide.shapes.add_picture(path, left, top, **kwargs)
        return True
    return False


# ===========================================================================
# Build presentation
# ===========================================================================
prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H
blank_layout = prs.slide_layouts[6]


# ---- SLIDE 1: Title -------------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, BG_DARKER)
add_accent_line(slide, Inches(0), Inches(0), SLIDE_W, TEAL, Inches(0.06))

add_textbox(slide, Inches(1.0), Inches(1.3), Inches(11.3), Inches(1.5),
            "CRISPRArchitect v3", font_size=52, color=WHITE, bold=True)
add_textbox(slide, Inches(1.0), Inches(2.5), Inches(11.3), Inches(1.2),
            "Multi-Nuclease Decision Support for\nGenome Editing Strategy Design",
            font_size=26, color=TEAL)
add_textbox(slide, Inches(1.0), Inches(3.5), Inches(11.3), Inches(0.5),
            "TOPSIS Multi-Criteria Ranking  |  Pareto Analysis  |  Sensitivity Quantification",
            font_size=16, color=CORAL)

add_accent_line(slide, Inches(1.0), Inches(4.3), Inches(3.0), GOLD, Inches(0.03))

add_textbox(slide, Inches(1.0), Inches(4.6), Inches(11.3), Inches(0.5),
            "Vishal Bharti  |  Debojyoti Chakraborty Lab", font_size=20, color=SOFT_WHITE)
add_textbox(slide, Inches(1.0), Inches(5.15), Inches(11.3), Inches(0.5),
            "CSIR-Institute of Genomics and Integrative Biology, New Delhi",
            font_size=16, color=DARK_GRAY)
add_textbox(slide, Inches(1.0), Inches(5.6), Inches(11.3), Inches(0.5),
            "Lab Meeting  |  April 2026", font_size=16, color=DARK_GRAY)

# Version badges
for i, (ver, col) in enumerate([("v1", DARK_GRAY), ("v2", TEAL), ("v3", CORAL)]):
    x = Inches(9.0) + Inches(i * 1.6)
    add_card(slide, x, Inches(5.8), Inches(1.3), Inches(0.5), BG_CARD)
    add_textbox(slide, x, Inches(5.85), Inches(1.3), Inches(0.4),
                ver, font_size=16, color=col, bold=True, alignment=PP_ALIGN.CENTER)

add_slide_number(slide, 1)


# ---- SLIDE 2: The Problem --------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "The Problem", "Why do we need a unified editing framework?")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.5), Inches(5.0), [
    "Multiple modalities (BE, PE, HDR) with distinct constraints",
    "No tool compares all three modalities side by side",
    "PAM + editing window is a hidden bottleneck (our key finding)",
    "iPSCs are p53-active: DSBs kill most cells (Ihry et al., 2018)",
    "Multi-nuclease options expand the design space enormously",
    "Clinical applications demand transparent, auditable decisions",
], font_size=17, line_spacing=1.5)

# Right cards
for i, (title, desc, col) in enumerate([
    ("Base Editing", "A>G or C>T only, narrow window (4-8 nt)", TEAL),
    ("Prime Editing", "Any edit, no DSB, complex pegRNA design", CORAL),
    ("HDR", "DSB required, p53 toxicity in iPSCs", GOLD),
]):
    y = Inches(1.6) + Inches(i * 1.6)
    add_card(slide, Inches(7.0), y, Inches(5.5), Inches(1.3))
    add_textbox(slide, Inches(7.2), y + Inches(0.1), Inches(5.1), Inches(0.4),
                title, font_size=16, color=col, bold=True)
    add_textbox(slide, Inches(7.2), y + Inches(0.55), Inches(5.1), Inches(0.5),
                desc, font_size=13, color=LIGHT_GRAY)

add_slide_number(slide, 2)


# ---- SLIDE 3: What is CRISPRArchitect? ------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "What is CRISPRArchitect?",
              "A decision-support tool, not a predictive optimizer")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(11.5), Inches(5.0), [
    ("Input: pathogenic variant(s) in HGVS notation or genomic coordinates", WHITE, True),
    "Automatically fetches transcript structure from Ensembl (GRCh38)",
    "Maps variant to CDS, annotates coding consequence (ACMG standards)",
    "Scans for PAM sites across 5 nucleases, evaluates 9 base editor profiles",
    "Assesses feasibility of BE, PE, and HDR with locus-specific constraints",
    ("Output: ranked strategies with TOPSIS scores, Pareto analysis,", WHITE, True),
    ("         sensitivity quantification, and explicit rejection reasons", WHITE, True),
    "Transparent: every score, weight, and parameter has documented provenance",
], font_size=16, line_spacing=1.5)

add_slide_number(slide, 3)


# ---- SLIDE 4: Pipeline Flowchart -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Pipeline Overview",
              "10-stage end-to-end workflow from variant input to ranked recommendations")

# Add the flowchart image (generated by generate_pipeline_flowchart.py)
fig_path = get_fig('Fig_Pipeline_Flowchart.png')
if fig_path:
    slide.shapes.add_picture(fig_path, Inches(2.5), Inches(1.3),
                             width=Inches(8.3), height=Inches(6.0))
else:
    # Fallback: text-based pipeline summary
    add_multiline(slide, Inches(0.7), Inches(1.5), Inches(11.5), Inches(5.5), [
        ("1. Variant Parsing  →  2. Transcript Fetch  →  3. Coordinate Mapping", TEAL, True),
        ("4. Validation & Annotation  →  5. Variant Normalization", WHITE, False),
        ("6. Multi-Nuclease PAM Scan (5 nucleases)", CORAL, True),
        ("7. Feasibility: Base Editing (9 profiles) | Prime Editing | HDR", WHITE, False),
        ("8. Strategy Generation (single, dual, hybrid)", WHITE, False),
        ("9. TOPSIS 6D Ranking (Safety, Feasibility, Complexity, Risk, Confidence, Consequence)", CORAL, True),
        ("10. Robustness: Pareto Front | Monte Carlo Sensitivity | Cross-Method", TEAL, True),
        ("", WHITE, False),
        ("Output: Ranked strategies with TOPSIS scores, rank stability, and rejection reasons", GOLD, True),
    ], font_size=15, line_spacing=1.6)

add_slide_number(slide, 4)


# ---- SLIDE 5: v1 Foundation ------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v1 Foundation: 6 Simulation Modules",
              "~24,000 lines of code | Monte Carlo biophysical simulations")

modules = [
    ("ConversionSim", "Gene conversion tract simulation (SDSA pathway)", "Monte Carlo"),
    ("MOSAIC", "Multi-locus strategy optimizer", "Decision engine"),
    ("ChromBridge", "3D chromatin distance predictor", "Polymer physics"),
    ("TopoPred", "cssDNA secondary structure analyzer", "Thermodynamics"),
    ("LoopSim", "Cohesin loop extrusion simulator", "Stochastic 1D"),
    ("WebApp", "Interactive Streamlit interface", "7 analysis pages"),
]
for i, (name, desc, method) in enumerate(modules):
    col_idx = i % 3
    row_idx = i // 3
    x = Inches(0.7) + Inches(col_idx * 4.2)
    y = Inches(1.5) + Inches(row_idx * 2.8)
    add_card(slide, x, y, Inches(3.9), Inches(2.3))
    add_textbox(slide, x + Inches(0.2), y + Inches(0.15), Inches(3.5), Inches(0.4),
                name, font_size=18, color=TEAL, bold=True)
    add_textbox(slide, x + Inches(0.2), y + Inches(0.6), Inches(3.5), Inches(0.8),
                desc, font_size=13, color=LIGHT_GRAY)
    add_textbox(slide, x + Inches(0.2), y + Inches(1.5), Inches(3.5), Inches(0.4),
                f"Method: {method}", font_size=11, color=DARK_GRAY)

add_slide_number(slide, 5)


# ---- SLIDE 6: v2 Recap -----------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v2: Transcript-Aware Pipeline",
              "~9,400 LOC | 30 ClinVar benchmark cases | 123 tests")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.5), Inches(4.5), [
    "Ensembl REST API integration (GRCh38 coordinates)",
    "Reference allele validation against genome sequence",
    "Coding consequence annotation (ACMG standards)",
    "PAM-verified feasibility: BE, PE, HDR engines",
    "Weighted-sum scoring with consequence adjustments",
    "30-case ClinVar benchmark: 86.7% top-1 accuracy",
], font_size=16, line_spacing=1.5)

# v2 result highlight
add_card(slide, Inches(7.0), Inches(1.5), Inches(5.5), Inches(5.0))
add_textbox(slide, Inches(7.2), Inches(1.7), Inches(5.1), Inches(0.4),
            "v2 Benchmark Results", font_size=18, color=TEAL, bold=True)
add_multiline(slide, Inches(7.2), Inches(2.3), Inches(5.1), Inches(3.5), [
    ("Top-1 accuracy: 86.7% (26/30)", WHITE, True),
    ("Top-3 accuracy: 96.7% (29/30)", WHITE, False),
    "",
    ("Strategy distribution:", DARK_GRAY, False),
    ("  PE: 29/30 (96.7%)", CORAL, True),
    ("  HDR: 1/30 (3.3%)", GOLD, False),
    ("  BE: 0/30 (0.0%)", RED_ACCENT, True),
    "",
    ("Problem: PE always wins.", RED_ACCENT, True),
    ("Even at ABE-compatible loci!", RED_ACCENT, False),
], font_size=15, line_spacing=1.3)

add_slide_number(slide, 6)


# ---- SLIDE 7: v2 Key Finding -----------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, RGBColor(0x15, 0x1E, 0x2B))
add_coral_banner(slide, "v2 KEY FINDING")

add_title_bar(slide, "PAM-Window Is the Real Bottleneck",
              "The most important result from v2")

add_card(slide, Inches(0.7), Inches(1.7), Inches(11.8), Inches(2.0),
         RGBColor(0x2A, 0x15, 0x10))
add_textbox(slide, Inches(1.0), Inches(1.9), Inches(11.2), Inches(0.6),
            "At ALL 7 ClinVar loci with ABE-compatible transitions,\n"
            "NO SpCas9 guide placed the target base within ABE positions 4-7.",
            font_size=20, color=CORAL, bold=True, alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1.0), Inches(2.7), Inches(11.2), Inches(0.5),
            "Mutation-type classification alone is INSUFFICIENT for strategy selection.",
            font_size=17, color=WHITE, alignment=PP_ALIGN.CENTER)

add_bullet_list(slide, Inches(0.7), Inches(4.2), Inches(11.5), Inches(2.5), [
    "A>G does NOT mean ABE will work at that locus",
    "The PAM must position the target A at protospacer position 4-7 (ABE7.10)",
    "This motivated v3: broader PAM nucleases + wider editing windows",
], font_size=16, line_spacing=1.5, bullet_color=CORAL)

add_slide_number(slide, 7)


# ---- SLIDE 8: v3 What's New ------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_coral_banner(slide, "VERSION 3.0")

add_title_bar(slide, "What's New in v3?",
              "Three major capabilities + scientific rigor overhaul")

cols = [
    ("Multi-Nuclease\nEngine", "5 nucleases × 9 editor profiles\n(3 Tier A + 6 Tier B fusions)\nrescued BE: 0/30 → 6/30",
     TEAL, Inches(0.7)),
    ("TOPSIS 6D\n+ Pareto", "6-dimension ranking\n+ weight-independent\nPareto dominance analysis",
     CORAL, Inches(4.9)),
    ("Statistical\nRigor", "SEs, 95% CIs, Wilson intervals\nDirichlet sensitivity\nAll params evidence-tagged",
     GOLD, Inches(9.1)),
]
for title, desc, col, x in cols:
    add_card(slide, x, Inches(1.7), Inches(3.8), Inches(4.5))
    add_textbox(slide, x + Inches(0.2), Inches(1.9), Inches(3.4), Inches(1.0),
                title, font_size=22, color=col, bold=True, alignment=PP_ALIGN.CENTER)
    add_accent_line(slide, x + Inches(0.5), Inches(3.0), Inches(2.8), col)
    add_textbox(slide, x + Inches(0.3), Inches(3.2), Inches(3.2), Inches(2.5),
                desc, font_size=14, color=LIGHT_GRAY, alignment=PP_ALIGN.CENTER)

add_textbox(slide, Inches(0.7), Inches(6.5), Inches(11.5), Inches(0.5),
            "Also: HGVS parser, ClinVar batch ingestion, CFD/MIT off-target scoring, "
            "VIKOR/WPM comparison methods",
            font_size=13, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

add_slide_number(slide, 8)


# ---- SLIDE 9: Multi-Nuclease Engine ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Multi-Nuclease Feasibility Engine",
              "Systematic evaluation of 5 nucleases × 9 base editor profiles")

# Nuclease table
nucleases = [
    ("SpCas9", "NGG", "Blunt", "1.0x", "Reference standard"),
    ("enFnCas9", "NRG", "~3 bp*", "1.5x*", "Our lab (Acharya et al., 2024)"),
    ("SpCas9-NG", "NG", "Blunt", "1.0x", "Nishimasu et al., 2018"),
    ("SpRY", "NNN", "Blunt", "0.8x", "Walton et al., 2020"),
    ("Cas12a", "TTTV", "5 bp", "1.4x", "Zetsche et al., 2015"),
]
header_y = Inches(1.5)
add_card(slide, Inches(0.5), header_y, Inches(12.3), Inches(0.5), TEAL)
for j, hdr in enumerate(["Nuclease", "PAM", "Stagger", "HDR mult.", "Reference"]):
    widths = [Inches(1.8), Inches(1.2), Inches(1.2), Inches(1.2), Inches(6.5)]
    x = Inches(0.6) + sum(w for w in [Inches(0)] + widths[:j])
    add_textbox(slide, x, header_y + Inches(0.05), widths[j], Inches(0.4),
                hdr, font_size=13, color=WHITE, bold=True)

for i, (name, pam, stag, mult, ref) in enumerate(nucleases):
    y = header_y + Inches(0.55) + Inches(i * 0.5)
    bg = BG_CARD if i % 2 == 0 else BG_NAVY
    add_card(slide, Inches(0.5), y, Inches(12.3), Inches(0.45), bg)
    vals = [name, pam, stag, mult, ref]
    widths = [Inches(1.8), Inches(1.2), Inches(1.2), Inches(1.2), Inches(6.5)]
    for j, (val, w) in enumerate(zip(vals, widths)):
        x = Inches(0.6) + sum(ww for ww in [Inches(0)] + widths[:j])
        col = TEAL if j == 0 else WHITE
        add_textbox(slide, x, y + Inches(0.05), w, Inches(0.35),
                    val, font_size=12, color=col, bold=(j == 0))

add_textbox(slide, Inches(0.7), Inches(4.5), Inches(11.5), Inches(0.4),
            "* enFnCas9 stagger is [ASSUMED] — not directly measured. "
            "HDR multiplier is [ASSUMED]. See parameter provenance.",
            font_size=11, color=DARK_GRAY)

# Editor summary
add_card(slide, Inches(0.5), Inches(5.0), Inches(12.3), Inches(2.2))
add_textbox(slide, Inches(0.7), Inches(5.1), Inches(11.5), Inches(0.4),
            "Base Editor Profiles (9 profiles: 3 Tier A + 6 Tier B)", font_size=16, color=TEAL, bold=True)
add_multiline(slide, Inches(0.7), Inches(5.55), Inches(4.0), Inches(1.5), [
    ("Tier A (core editors):", TEAL, True),
    ("ABE7.10: window 4-7", WHITE, False),
    ("ABE8e: window 3-9 extended", TEAL, False),
    ("BE4max (CBE): window 4-8", WHITE, False),
], font_size=12, line_spacing=1.2)
add_multiline(slide, Inches(4.8), Inches(5.55), Inches(4.0), Inches(1.5), [
    ("Tier B (nuclease fusions — ABE):", CORAL, True),
    ("ABE8e + enFnCas9: NRG PAM", CORAL, False),
    ("ABE8e + SpCas9-NG: NG PAM", WHITE, False),
    ("ABE8e + SpRY: near-PAMless", WHITE, False),
], font_size=12, line_spacing=1.2)
add_multiline(slide, Inches(8.8), Inches(5.55), Inches(4.0), Inches(1.5), [
    ("Tier B (nuclease fusions — CBE):", CORAL, True),
    ("BE4max + enFnCas9: NRG PAM", CORAL, False),
    ("BE4max + SpCas9-NG: NG PAM", WHITE, False),
    ("BE4max + SpRY: near-PAMless", WHITE, False),
], font_size=12, line_spacing=1.2)

add_slide_number(slide, 9)


# ---- SLIDE 10: enFnCas9 Advantage ------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "enFnCas9: Our Lab's Engineered Nuclease",
              "Acharya et al., Nature Communications 15:5471 (2024)")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.5), Inches(5.0), [
    ("NRG PAM: ~2x more targetable sites than NGG", WHITE, True),
    "Single-nucleobase specificity (lower off-target than SpCas9)",
    "Improved HDR knock-in rates demonstrated",
    "Three engineered variants: en1, en15, en31",
    "RPE65 correction in LCA2 patient iPSCs demonstrated",
], font_size=16, line_spacing=1.5)

# What we know vs assume
add_card(slide, Inches(7.0), Inches(1.5), Inches(5.5), Inches(2.2))
add_textbox(slide, Inches(7.2), Inches(1.6), Inches(5.1), Inches(0.4),
            "[MEASURED] from paper", font_size=14, color=GREEN_ACC, bold=True)
add_multiline(slide, Inches(7.2), Inches(2.1), Inches(5.1), Inches(1.2), [
    ("NRG/NGR PAM recognition", WHITE, False),
    ("Single-nucleobase specificity", WHITE, False),
    ("Improved HDR knock-in", WHITE, False),
], font_size=13, line_spacing=1.3)

add_card(slide, Inches(7.0), Inches(4.0), Inches(5.5), Inches(2.5))
add_textbox(slide, Inches(7.2), Inches(4.1), Inches(5.1), Inches(0.4),
            "[ASSUMED] — need experimental data", font_size=14, color=CORAL, bold=True)
add_multiline(slide, Inches(7.2), Inches(4.6), Inches(5.1), Inches(1.5), [
    ("Stagger: ~3 bp 5' overhang (from FnCas9 structure)", CORAL, False),
    ("HDR multiplier: 1.5x (qualitative from paper)", CORAL, False),
    ("ABE8e-enFnCas9 window: 3-9 (extrapolated)", CORAL, False),
    ("", WHITE, False),
    ("Data request sent to PI", GOLD, True),
], font_size=13, line_spacing=1.3)

add_slide_number(slide, 10)


# ---- SLIDE 11: ABE8e Rescue Mechanism --------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "How Multi-Nuclease Rescues Base Editing",
              "ABE8e + enFnCas9 = broader PAM + wider window")

# Before vs After
add_card(slide, Inches(0.7), Inches(1.5), Inches(5.6), Inches(5.0))
add_textbox(slide, Inches(0.9), Inches(1.6), Inches(5.2), Inches(0.4),
            "v2 (SpCas9 + ABE7.10 only)", font_size=18, color=RED_ACCENT, bold=True)
add_multiline(slide, Inches(0.9), Inches(2.2), Inches(5.2), Inches(3.5), [
    ("PAM: NGG only (~8% of positions)", WHITE, False),
    ("Window: positions 4-7 (4 nt)", WHITE, False),
    ("Result: 0/30 BE top-ranked", RED_ACCENT, True),
    ("", WHITE, False),
    ("Even at 7 ABE-compatible loci,", LIGHT_GRAY, False),
    ("no guide placed target in window", LIGHT_GRAY, False),
], font_size=15, line_spacing=1.4)

add_card(slide, Inches(6.8), Inches(1.5), Inches(5.8), Inches(5.0))
add_textbox(slide, Inches(7.0), Inches(1.6), Inches(5.4), Inches(0.4),
            "v3 (5 nucleases + ABE8e)", font_size=18, color=GREEN_ACC, bold=True)
add_multiline(slide, Inches(7.0), Inches(2.2), Inches(5.4), Inches(3.5), [
    ("PAM: NGG + NRG + NG + NNN + TTTV", WHITE, False),
    ("Window: positions 3-9 (7 nt, ABE8e)", WHITE, False),
    ("Result: 6/30 BE top-ranked", GREEN_ACC, True),
    ("", WHITE, False),
    ("ABE8e window (3-9) is ~3x broader", TEAL, False),
    ("enFnCas9 (NRG) ~2x more PAM sites", TEAL, False),
    ("Combined: rescues BE at 6 loci", TEAL, True),
], font_size=15, line_spacing=1.4)

add_slide_number(slide, 11)


# ---- SLIDE 12: KEY FINDING - BE Rescue -------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, RGBColor(0x15, 0x1E, 0x2B))
add_coral_banner(slide, "v3 KEY FINDING")

add_title_bar(slide, "Base Editing Rescued: 0/30 to 6/30")

add_card(slide, Inches(0.7), Inches(1.7), Inches(11.8), Inches(1.5),
         RGBColor(0x0A, 0x2A, 0x15))
add_textbox(slide, Inches(1.0), Inches(1.9), Inches(11.2), Inches(1.0),
            "v3 strategy distribution: BE 20% | PE 77% | HDR 3%\n"
            "(v2 was: BE 0% | PE 97% | HDR 3%)",
            font_size=20, color=GREEN_ACC, bold=True, alignment=PP_ALIGN.CENTER)

# Strategy distribution figure (real data)
add_figure(slide, 'Fig_StrategyDistribution_v2_v3.png',
           Inches(0.7), Inches(3.3), width=Inches(6.5), height=Inches(3.5))

# Key stats on the right
add_stat_box(slide, Inches(7.8), Inches(3.5), Inches(2.3), Inches(1.0),
             "6/30", "BE v3", GREEN_ACC)
add_stat_box(slide, Inches(10.3), Inches(3.5), Inches(2.3), Inches(1.0),
             "0/30", "BE v2", RED_ACCENT)
add_textbox(slide, Inches(7.8), Inches(5.0), Inches(4.8), Inches(1.5),
            "Primary drivers:\n"
            "  ABE8e broader window (3-9)\n"
            "  enFnCas9 NRG PAM\n"
            "Bug fixed: bystander triple-count",
            font_size=13, color=LIGHT_GRAY)

add_slide_number(slide, 12)


# ---- SLIDE 13: 6D TOPSIS Scoring ------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "6-Dimensional TOPSIS Scoring",
              "Technique for Order Preference by Similarity to Ideal Solution (Hwang & Yoon, 1981)")

# The 6 dimensions
dims = [
    ("Safety", "0.28", "DSBs, p53 risk", TEAL, "Benefit"),
    ("Feasibility", "0.23", "PAM + window verified", TEAL, "Benefit"),
    ("Complexity", "0.19", "Rounds, donors, screening", RED_ACCENT, "Cost"),
    ("Risk", "0.14", "Rearrangement only*", RED_ACCENT, "Cost"),
    ("Confidence", "0.09", "Evidence tier (A/B/C)", TEAL, "Benefit"),
    ("Consequence", "0.07", "Bystander edits, splice", TEAL, "Benefit"),
]
for i, (name, weight, desc, col, dim_type) in enumerate(dims):
    col_idx = i % 3
    row_idx = i // 3
    x = Inches(0.5) + Inches(col_idx * 4.2)
    y = Inches(1.5) + Inches(row_idx * 2.5)
    add_card(slide, x, y, Inches(3.9), Inches(2.1))
    add_textbox(slide, x + Inches(0.2), y + Inches(0.1), Inches(3.5), Inches(0.4),
                f"{name} (w={weight})", font_size=16, color=col, bold=True)
    add_textbox(slide, x + Inches(0.2), y + Inches(0.6), Inches(3.5), Inches(0.5),
                desc, font_size=13, color=LIGHT_GRAY)
    add_textbox(slide, x + Inches(0.2), y + Inches(1.2), Inches(3.5), Inches(0.4),
                f"Type: {dim_type} dimension", font_size=11, color=DARK_GRAY)

add_textbox(slide, Inches(0.5), Inches(6.7), Inches(12.0), Inches(0.5),
            "* Risk now captures ONLY structural rearrangement. Bystander edits moved to "
            "Consequence dimension to prevent double-counting.",
            font_size=11, color=DARK_GRAY)

add_slide_number(slide, 13)


# ---- SLIDE 14: Bystander Bug Fix ------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, RGBColor(0x15, 0x1E, 0x2B))
add_coral_banner(slide, "CRITICAL BUG FIX")

add_title_bar(slide, "Bystander Triple-Counting Eliminated",
              "Root cause of the degenerate 'always PE' pattern in v2")

add_card(slide, Inches(0.7), Inches(1.7), Inches(5.6), Inches(4.5))
add_textbox(slide, Inches(0.9), Inches(1.8), Inches(5.2), Inches(0.4),
            "v2 Bug: 3 Penalty Channels", font_size=18, color=RED_ACCENT, bold=True)
add_multiline(slide, Inches(0.9), Inches(2.4), Inches(5.2), Inches(3.0), [
    ("1. Risk dimension: severity x 0.3 x w=0.15", WHITE, False),
    ("   Penalty for 1 bystander: 0.009", LIGHT_GRAY, False),
    ("2. Consequence penalty: severity x 0.08", WHITE, False),
    ("   Penalty for 1 bystander: 0.016", LIGHT_GRAY, False),
    ("3. Lost clean design bonus: -0.03", WHITE, False),
    ("   (PE gets +0.03, BE loses it)", LIGHT_GRAY, False),
    ("", WHITE, False),
    ("TOTAL penalty for 1 bystander: 0.055", RED_ACCENT, True),
    ("BE feasibility advantage: only 0.033", RED_ACCENT, True),
    ("=> PE always wins with any bystander", RED_ACCENT, True),
], font_size=14, line_spacing=1.3)

add_card(slide, Inches(6.8), Inches(1.7), Inches(5.8), Inches(4.5))
add_textbox(slide, Inches(7.0), Inches(1.8), Inches(5.4), Inches(0.4),
            "v3 Fix: Single Consequence Dimension", font_size=18, color=GREEN_ACC, bold=True)
add_multiline(slide, Inches(7.0), Inches(2.4), Inches(5.4), Inches(3.0), [
    ("Bystander removed from Risk dimension", WHITE, False),
    ("Consequence is now a proper 6th TOPSIS dim", WHITE, False),
    ("No additive post-hoc penalty", WHITE, False),
    ("Clean bonus is a tie-breaker only", WHITE, False),
    ("", WHITE, False),
    ("Result: BE wins with 0-3 bystanders", GREEN_ACC, True),
    ("PE wins only with 4+ bystanders", GREEN_ACC, False),
    ("Ranking reflects true biology", GREEN_ACC, True),
    ("", WHITE, False),
    ("All 3 MCDM methods agree: 100%", TEAL, True),
    ("concordance (TOPSIS, VIKOR, WPM)", TEAL, False),
], font_size=14, line_spacing=1.3)

# Bystander fix figure (inset, lower right area)
add_figure(slide, 'Fig_BystanterFix.png',
           Inches(6.5), Inches(4.5), width=Inches(6.2), height=Inches(2.7))

add_slide_number(slide, 14)


# ---- SLIDE 15: Pareto Front Analysis --------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Pareto Front Analysis",
              "Weight-independent dominance check across all 6 dimensions")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(11.5), Inches(2.5), [
    "Strategy A dominates B if A is at least as good on ALL dimensions and strictly better on at least one",
    "Pareto non-dominated strategies form the 'Pareto front' — defensible under SOME weighting scheme",
    "This analysis requires NO weight assumptions — it is purely structural",
    ("Complement to TOPSIS: shows which strategies COULD be optimal under any weights", TEAL, True),
], font_size=16, line_spacing=1.5)

# Example
add_card(slide, Inches(0.7), Inches(4.5), Inches(11.8), Inches(2.5))
add_textbox(slide, Inches(0.9), Inches(4.6), Inches(11.4), Inches(0.4),
            "Typical Result: BE vs PE vs HDR (1 bystander)", font_size=16, color=TEAL, bold=True)
add_multiline(slide, Inches(0.9), Inches(5.2), Inches(5.5), Inches(1.5), [
    ("BE: Pareto non-dominated", GREEN_ACC, True),
    ("  (better feasibility, slightly worse consequence)", LIGHT_GRAY, False),
    ("PE: Pareto non-dominated", GREEN_ACC, True),
    ("  (better consequence, lower feasibility)", LIGHT_GRAY, False),
], font_size=14, line_spacing=1.3)
add_multiline(slide, Inches(6.8), Inches(5.2), Inches(5.5), Inches(1.5), [
    ("HDR: Pareto DOMINATED", RED_ACCENT, True),
    ("  (worse on safety AND feasibility AND confidence)", LIGHT_GRAY, False),
    ("  => HDR is never optimal regardless of weights", LIGHT_GRAY, False),
], font_size=14, line_spacing=1.3)

add_slide_number(slide, 15)


# ---- SLIDE 16: Sensitivity Analysis ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Monte Carlo Sensitivity Analysis",
              "10,000 Dirichlet-sampled weight permutations (concentration=20)")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(11.5), Inches(2.0), [
    "Instead of reporting a single 'best' strategy, we report: How stable is the ranking?",
    "Rank stability = fraction of 10,000 weight permutations where strategy is #1",
    "Dirichlet distribution with concentration=20, min_alpha=2.0 (prevents pathological weights)",
], font_size=15, line_spacing=1.5)

# Example output
add_card(slide, Inches(0.7), Inches(3.8), Inches(11.8), Inches(3.0))
add_textbox(slide, Inches(0.9), Inches(3.9), Inches(11.4), Inches(0.4),
            "Example Output (typical case):", font_size=16, color=TEAL, bold=True)
add_multiline(slide, Inches(0.9), Inches(4.5), Inches(11.4), Inches(2.0), [
    ("#1 Single-step Base Editing     TOPSIS=0.995  Stability=94.2%  Pareto: non-dominated", GREEN_ACC, False),
    ("#2 Single-step Prime Editing    TOPSIS=0.869  Stability=5.8%   Pareto: non-dominated", WHITE, False),
    ("#3 Single-step HDR              TOPSIS=0.005  Stability=0.0%   Pareto: dominated", DARK_GRAY, False),
    ("", WHITE, False),
    ("Interpretation: BE is top-ranked in 94.2% of weight permutations.", LIGHT_GRAY, False),
    ("The recommendation is ROBUST — it does not depend on the specific weight values.", TEAL, True),
], font_size=14, line_spacing=1.3)

add_slide_number(slide, 16)


# ---- SLIDE 17: Cross-Method Validation -------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Cross-Method Validation: TOPSIS + VIKOR + WPM",
              "Three independent MCDM methods agree on rankings")

# Method descriptions
methods = [
    ("TOPSIS", "Euclidean distance to\nideal/anti-ideal solutions\n(Hwang & Yoon, 1981)", TEAL),
    ("VIKOR", "Compromise solution\nminimizing group utility\n(Opricovic & Tzeng, 2004)", CORAL),
    ("WPM", "Weighted geometric mean\nNon-compensatory\n(Bridgman, 1922)", GOLD),
]
for i, (name, desc, col) in enumerate(methods):
    x = Inches(0.5) + Inches(i * 4.2)
    add_card(slide, x, Inches(1.5), Inches(3.9), Inches(2.5))
    add_textbox(slide, x + Inches(0.2), Inches(1.6), Inches(3.5), Inches(0.5),
                name, font_size=22, color=col, bold=True, alignment=PP_ALIGN.CENTER)
    add_accent_line(slide, x + Inches(0.5), Inches(2.2), Inches(2.9), col)
    add_textbox(slide, x + Inches(0.3), Inches(2.4), Inches(3.3), Inches(1.2),
                desc, font_size=13, color=LIGHT_GRAY, alignment=PP_ALIGN.CENTER)

add_card(slide, Inches(2.0), Inches(4.5), Inches(9.3), Inches(2.5),
         RGBColor(0x0A, 0x2A, 0x15))
add_textbox(slide, Inches(2.5), Inches(4.7), Inches(8.3), Inches(0.5),
            "Rank Concordance: 100%", font_size=28, color=GREEN_ACC, bold=True,
            alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(2.5), Inches(5.4), Inches(8.3), Inches(1.0),
            "All three methods produce identical rankings across test cases.\n"
            "The recommendation is method-robust: it does not depend on the\n"
            "choice of MCDM algorithm.",
            font_size=15, color=WHITE, alignment=PP_ALIGN.CENTER)

# Cross-method figure (replaces the green box)
add_figure(slide, 'Fig_CrossMethod_Agreement.png',
           Inches(2.0), Inches(4.3), width=Inches(9.3), height=Inches(2.8))

add_slide_number(slide, 17)


# ---- SLIDE 18: Statistical Rigor -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Statistical Rigor in ConversionSim",
              "Every output now includes uncertainty quantification")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.5), Inches(5.5), [
    "HDR rate: binomial SE + 95% CI",
    "Mean tract length: SE = std/sqrt(n) + 95% CI",
    "Distance probabilities: Wilson score 95% CI",
    "Sample std with ddof=1 (not population std)",
    "Tract clipping bias removed (lower bound)",
    "SeedSequence spawning (not seed+i)",
    "Numpy vectorized Dirichlet (not stdlib)",
    "Geometric distribution explicitly justified",
], font_size=15, line_spacing=1.5)

# Example output
add_card(slide, Inches(6.5), Inches(1.5), Inches(6.3), Inches(5.5))
add_textbox(slide, Inches(6.7), Inches(1.6), Inches(5.9), Inches(0.4),
            "Example: enFnCas9 + cssDNA + iPSC", font_size=14, color=TEAL, bold=True)
add_multiline(slide, Inches(6.7), Inches(2.2), Inches(5.9), Inches(4.5), [
    ("HDR rate: 3.5% (95% CI: [3.2%, 3.9%])", WHITE, False),
    ("", WHITE, False),
    ("Mean tract: 706 bp (SE=36.9)", WHITE, False),
    ("  95% CI: [634, 778] bp", TEAL, False),
    ("Median: 454 bp", WHITE, False),
    ("", WHITE, False),
    ("P(tract >= 500 bp):", WHITE, False),
    ("  47.2% [42.0%, 52.4%]", TEAL, True),
    ("P(tract >= 1000 bp):", WHITE, False),
    ("  25.9% [21.6%, 30.7%]", TEAL, False),
    ("", WHITE, False),
    ("n = 10,000 simulations, seed = 42", DARK_GRAY, False),
], font_size=13, line_spacing=1.25)

# ConversionSim CI figure (replaces the card with example data)
add_figure(slide, 'Fig_ConversionSim_CIs.png',
           Inches(0.5), Inches(3.8), width=Inches(12.3), height=Inches(3.5))

add_slide_number(slide, 18)


# ---- SLIDE 19: ConversionSim Scope ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "ConversionSim: Honest Scope & Validation",
              "SDSA pathway model for long-donor HDR only")

# Valid / Not valid
add_card(slide, Inches(0.7), Inches(1.5), Inches(5.6), Inches(2.8))
add_textbox(slide, Inches(0.9), Inches(1.6), Inches(5.2), Inches(0.4),
            "Valid for (SDSA pathway):", font_size=16, color=GREEN_ACC, bold=True)
add_bullet_list(slide, Inches(0.9), Inches(2.1), Inches(5.2), Inches(1.8), [
    "cssDNA donors (Iyer et al., 2022: 2.07x predicted vs 1.9x observed)",
    "Staggered cuts (Chauhan et al., 2023: 1.82x pred vs 1.9x obs)",
    "dsDNA/lssDNA donors with HA >= 100 bp",
], font_size=13, line_spacing=1.4, bullet_color=GREEN_ACC)

add_card(slide, Inches(6.8), Inches(1.5), Inches(5.8), Inches(2.8))
add_textbox(slide, Inches(7.0), Inches(1.6), Inches(5.4), Inches(0.4),
            "NOT valid for (SSTR pathway):", font_size=16, color=RED_ACCENT, bold=True)
add_bullet_list(slide, Inches(7.0), Inches(2.1), Inches(5.4), Inches(1.8), [
    "ssODN editing (uses SSTR, not SDSA)",
    "R\u00B2 = -0.56 vs Paquet 2016 (expected mismatch)",
    "Prime editing (RT-based, not strand invasion)",
], font_size=13, line_spacing=1.4, bullet_color=RED_ACCENT)

# Parameter provenance
add_card(slide, Inches(0.7), Inches(4.8), Inches(11.8), Inches(2.2))
add_textbox(slide, Inches(0.9), Inches(4.9), Inches(11.4), Inches(0.4),
            "Central parameter: SDSA displacement probability p = 0.002", font_size=16,
            color=TEAL, bold=True)
add_multiline(slide, Inches(0.9), Inches(5.4), Inches(11.4), Inches(1.5), [
    ("Evidence: [ASSUMED] from functional data (not direct measurement)", CORAL, True),
    ("  Stark lab SDSA assay: >= 350 bp synthesis confirmed in human cells (G3, 2017)", LIGHT_GRAY, False),
    ("  Successful HDR with 300-1000 bp homology arms implies routine incorporation at hundreds of bp", LIGHT_GRAY, False),
    ("  NOT calibrated to: Elliott 1998 (endogenous tracts <58 bp) or Kan 2017 (SSTR tracts ~20 bp)", DARK_GRAY, False),
], font_size=12, line_spacing=1.3)

add_slide_number(slide, 19)


# ---- SLIDE 20: Citation Integrity ------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Citation Integrity: Web-Search Verified",
              "Every reference and calibration value checked against primary sources")

add_card(slide, Inches(0.7), Inches(1.5), Inches(11.8), Inches(5.0))
add_textbox(slide, Inches(0.9), Inches(1.6), Inches(11.4), Inches(0.4),
            "Errors Found and Corrected", font_size=18, color=CORAL, bold=True)

errors = [
    ("Arbab et al., 2020", "Listed as Nature 584:268", "Actually Cell 182:463", "FIXED"),
    ("enFnCas9 (Ref 14)", "Wrong first author & title", "Acharya S et al., correct title", "FIXED"),
    ("Kan et al., 2017", "Listed as Mol Cell", "Actually Genome Research 27:1099", "FIXED"),
    ("Elliott et al., 1998", "Claimed tracts 200-2000 bp", "Actually 80% of tracts <= 58 bp", "REFRAMED"),
    ("HA 300 nt optimal", "Attributed to Iyer 2022", "Paper doesn't specify 300 nt", "REFRAMED"),
    ("Kim et al., 2019", "13,000 BE targets", "Correct: Song et al., Nat Biotech 2020", "DOCUMENTED"),
]
header_y2 = Inches(2.2)
for j, hdr in enumerate(["Paper", "Error", "Correction", "Status"]):
    widths2 = [Inches(2.5), Inches(3.5), Inches(3.8), Inches(1.5)]
    x2 = Inches(0.9) + sum(w for w in [Inches(0)] + widths2[:j])
    add_textbox(slide, x2, header_y2, widths2[j], Inches(0.35),
                hdr, font_size=12, color=TEAL, bold=True)

for i, (paper, error, fix, status) in enumerate(errors):
    y2 = header_y2 + Inches(0.4) + Inches(i * 0.55)
    vals2 = [paper, error, fix, status]
    widths2 = [Inches(2.5), Inches(3.5), Inches(3.8), Inches(1.5)]
    for j, (val, w) in enumerate(zip(vals2, widths2)):
        x2 = Inches(0.9) + sum(ww for ww in [Inches(0)] + widths2[:j])
        col2 = GREEN_ACC if status == "FIXED" and j == 3 else (CORAL if j == 1 else WHITE)
        add_textbox(slide, x2, y2, w, Inches(0.45), val, font_size=11, color=col2)

add_textbox(slide, Inches(0.9), Inches(5.8), Inches(11.4), Inches(0.4),
            "20 references verified via web search (PubMed/DOI). "
            "17 verified, 3 corrected. All PMIDs confirmed.",
            font_size=13, color=DARK_GRAY)

add_slide_number(slide, 20)


# ---- SLIDE 21: Parameter Provenance ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Parameter Provenance System",
              "Every constant tagged [MEASURED], [DERIVED], or [ASSUMED]")

# Evidence level definitions
levels = [
    ("[MEASURED]", "Value directly from a published measurement with citation", GREEN_ACC),
    ("[DERIVED]", "Computed from published data via a stated mathematical procedure", TEAL),
    ("[ASSUMED]", "Modeling assumption with rationale; explored via sensitivity analysis", CORAL),
]
for i, (tag, desc, col) in enumerate(levels):
    y = Inches(1.5) + Inches(i * 0.8)
    add_card(slide, Inches(0.7), y, Inches(11.8), Inches(0.65), BG_CARD)
    add_textbox(slide, Inches(0.9), y + Inches(0.1), Inches(2.0), Inches(0.4),
                tag, font_size=16, color=col, bold=True)
    add_textbox(slide, Inches(3.0), y + Inches(0.1), Inches(9.3), Inches(0.4),
                desc, font_size=13, color=LIGHT_GRAY)

# Key parameters table
add_textbox(slide, Inches(0.7), Inches(4.2), Inches(11.5), Inches(0.4),
            "Key Parameters with Evidence Tags:", font_size=16, color=WHITE, bold=True)

params = [
    ("SpCas9 cut position: -3 bp", "[MEASURED]", "Jinek et al., Science 2012"),
    ("ABE7.10 window: 4-7", "[MEASURED]", "Gaudelli et al., Nature 2017"),
    ("cssDNA HDR: ~1.9x vs lssDNA", "[MEASURED]", "Iyer et al., CRISPR J 2022"),
    ("vCas9 stagger: >= 6 bp, 1.9x HDR", "[MEASURED]", "Chauhan et al., PNAS 2023"),
    ("iPSC viability post-DSB: 55%", "[MEASURED]", "Ihry et al., Nat Med 2018"),
    ("SDSA displacement prob: p=0.002", "[ASSUMED]", "Functional evidence (see scope slide)"),
    ("enFnCas9 stagger: 3 bp", "[ASSUMED]", "FnCas9 structure inference"),
    ("Scoring weights: 0.28/0.23/...", "[ASSUMED]", "iPSC-context rationale + sensitivity"),
]
for i, (param, tag, source) in enumerate(params):
    y = Inches(4.7) + Inches(i * 0.33)
    col3 = GREEN_ACC if "MEASURED" in tag else (TEAL if "DERIVED" in tag else CORAL)
    add_textbox(slide, Inches(0.9), y, Inches(4.5), Inches(0.3),
                param, font_size=11, color=WHITE)
    add_textbox(slide, Inches(5.5), y, Inches(1.8), Inches(0.3),
                tag, font_size=11, color=col3, bold=True)
    add_textbox(slide, Inches(7.4), y, Inches(5.0), Inches(0.3),
                source, font_size=10, color=DARK_GRAY)

# Parameter provenance pie chart (top right)
add_figure(slide, 'Fig_ParameterProvenance.png',
           Inches(8.5), Inches(1.5), width=Inches(4.3), height=Inches(4.3))

add_slide_number(slide, 21)


# ---- SLIDE 22: v3 Benchmark Results ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Benchmark Results: v2 vs v3",
              "30 ClinVar cases, verified GRCh38 coordinates")

# Comparison table
add_card(slide, Inches(0.7), Inches(1.5), Inches(5.6), Inches(5.0))
add_textbox(slide, Inches(0.9), Inches(1.6), Inches(5.2), Inches(0.4),
            "v2 Results (SpCas9 + ABE7.10)", font_size=16, color=TEAL, bold=True)
add_multiline(slide, Inches(0.9), Inches(2.2), Inches(5.2), Inches(3.5), [
    ("Top-1 accuracy: 86.7% (26/30)", WHITE, True),
    ("Top-3 accuracy: 96.7% (29/30)", WHITE, False),
    ("Rejection accuracy: 90.0% (27/30)", WHITE, False),
    ("", WHITE, False),
    ("Strategy: PE=29, HDR=1, BE=0", WHITE, False),
    ("", WHITE, False),
    ("Scoring: 5D weighted sum", DARK_GRAY, False),
    ("Tests: 123 passing", DARK_GRAY, False),
], font_size=14, line_spacing=1.3)

add_card(slide, Inches(6.8), Inches(1.5), Inches(5.8), Inches(5.0))
add_textbox(slide, Inches(7.0), Inches(1.6), Inches(5.4), Inches(0.4),
            "v3 Results (multi-nuclease + TOPSIS 6D)", font_size=16, color=CORAL, bold=True)
add_multiline(slide, Inches(7.0), Inches(2.2), Inches(5.4), Inches(3.5), [
    ("Top-1 accuracy: 86.7% (26/30)", WHITE, True),
    ("Top-3 accuracy: 96.7% (29/30)", WHITE, False),
    ("Rejection accuracy: 86.7% (26/30)", WHITE, False),
    ("", WHITE, False),
    ("Strategy: PE=23, BE=6, HDR=1", GREEN_ACC, True),
    ("", WHITE, False),
    ("Scoring: 6D TOPSIS + Pareto + VIKOR/WPM", CORAL, False),
    ("Tests: 198 passing", CORAL, False),
    ("SEs/CIs on all outputs", CORAL, False),
    ("All 20 references web-verified", CORAL, False),
], font_size=14, line_spacing=1.3)

# Literature benchmark figure (bottom center)
add_figure(slide, 'Fig_LiteratureBenchmark.png',
           Inches(0.5), Inches(5.0), width=Inches(12.3), height=Inches(2.3))

add_slide_number(slide, 22)


# ---- SLIDE 23: Codebase Summary -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v3 Codebase Summary",
              "Version 3.0.0 | MIT License | github.com/visvikbharti/CRISPRArchitect")

stats_data = [
    ("198", "Tests Passing", TEAL),
    ("3.0.0", "Version", CORAL),
    ("~36,500", "Lines of Code", GOLD),
    ("20", "Verified Refs", GREEN_ACC),
    ("6", "TOPSIS Dims", BLUE_ACC),
    ("0", "Test Failures", GREEN_ACC),
]
for i, (val, label, col) in enumerate(stats_data):
    col_idx = i % 3
    row_idx = i // 3
    x = Inches(0.7) + Inches(col_idx * 4.2)
    y = Inches(1.5) + Inches(row_idx * 2.0)
    add_stat_box(slide, x, y, Inches(3.5), Inches(1.5), val, label, col)

add_bullet_list(slide, Inches(0.7), Inches(5.7), Inches(11.5), Inches(1.5), [
    "Python 3.9+ | NumPy, SciPy, Matplotlib | Streamlit web app",
    "Docker + docker-compose ready | GitHub Actions CI (Python 3.9-3.12)",
    "HGVS parser + ClinVar batch ingestion | CLI + web interface",
], font_size=14, line_spacing=1.5)

add_slide_number(slide, 23)


# ---- SLIDE 24: Limitations ------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Limitations (Honest)",
              "What we cannot do and what we need")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(11.5), Inches(5.5), [
    ("No experimental validation yet — benchmark uses self-curated truth labels", RED_ACCENT, True),
    "ConversionSim models SDSA only; ssODN (SSTR) predictions are invalid",
    "SDSA displacement probability p=0.002 is [ASSUMED], not directly measured",
    "enFnCas9 stagger (3 bp) and HDR multiplier (1.5x) are not experimentally verified",
    "Off-target scoring is local only (no genome-wide search)",
    "No chromatin accessibility or replication timing data integration",
    "Large deletion strategies (multi-exon) are poorly handled (0/3 accuracy)",
    "Benchmark has degenerate pattern: PE still dominates in 23/30 cases",
    "Need independent expert labeling of benchmark cases (Cohen's kappa)",
    "Need comparison against empirical editing data (Song 2020, Chen 2021)",
], font_size=14, line_spacing=1.5, bullet_color=CORAL)

add_slide_number(slide, 24)


# ---- SLIDE 25: What's Next ------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "What's Next",
              "Roadmap to PLOS Computational Biology submission")

add_card(slide, Inches(0.7), Inches(1.5), Inches(5.6), Inches(5.0))
add_textbox(slide, Inches(0.9), Inches(1.6), Inches(5.2), Inches(0.4),
            "From the Lab (data request sent)", font_size=16, color=GOLD, bold=True)
add_bullet_list(slide, Inches(0.9), Inches(2.2), Inches(5.2), Inches(3.5), [
    "3-5 iPSC editing cases (any gene, any strategy)",
    "enFnCas9 cut-site stagger measurement",
    "enFnCas9 vs SpCas9 HDR comparison (2-3 loci)",
    "enFnCas9-ABE8e editing window (if tested)",
], font_size=14, line_spacing=1.5, bullet_color=GOLD)

add_card(slide, Inches(6.8), Inches(1.5), Inches(5.8), Inches(5.0))
add_textbox(slide, Inches(7.0), Inches(1.6), Inches(5.4), Inches(0.4),
            "Computational (in progress)", font_size=16, color=TEAL, bold=True)
add_bullet_list(slide, Inches(7.0), Inches(2.2), Inches(5.4), Inches(3.5), [
    "Expand benchmark to 30+ published cases",
    "Validate BE against Song et al. 2020 data",
    "Validate PE against Chen et al. 2021 data",
    "Rebuild manuscript for PLOS Comp Bio",
    "Generate publication-quality figures",
    "Independent expert labeling of cases",
], font_size=14, line_spacing=1.5, bullet_color=TEAL)

add_slide_number(slide, 25)


# ---- SLIDE 26: Thank You --------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, BG_DARKER)
add_accent_line(slide, Inches(0), Inches(0), SLIDE_W, TEAL, Inches(0.06))

add_textbox(slide, Inches(1.0), Inches(1.5), Inches(11.3), Inches(1.0),
            "Thank You", font_size=48, color=WHITE, bold=True)
add_textbox(slide, Inches(1.0), Inches(2.8), Inches(11.3), Inches(0.5),
            "Questions & Discussion", font_size=24, color=TEAL)

add_accent_line(slide, Inches(1.0), Inches(3.6), Inches(3.0), GOLD, Inches(0.03))

add_textbox(slide, Inches(1.0), Inches(4.0), Inches(11.3), Inches(0.5),
            "Vishal Bharti  |  Debojyoti Chakraborty Lab  |  CSIR-IGIB",
            font_size=18, color=SOFT_WHITE)
add_textbox(slide, Inches(1.0), Inches(4.6), Inches(11.3), Inches(0.5),
            "github.com/visvikbharti/CRISPRArchitect  |  MIT License  |  v3.0.0",
            font_size=14, color=DARK_GRAY)

add_textbox(slide, Inches(1.0), Inches(5.5), Inches(11.3), Inches(1.0),
            "Key takeaway: CRISPRArchitect v3 provides transparent,\n"
            "method-robust, uncertainty-quantified strategy recommendations\n"
            "with every parameter traceable to published evidence.",
            font_size=16, color=LIGHT_GRAY, alignment=PP_ALIGN.CENTER)

add_slide_number(slide, 26)


# ===========================================================================
# Save
# ===========================================================================
output_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(output_dir, "CRISPRArchitect_v3_LabMeeting.pptx")
prs.save(output_path)
print(f"Presentation saved to: {output_path}")
print(f"Total slides: {len(prs.slides)}")
