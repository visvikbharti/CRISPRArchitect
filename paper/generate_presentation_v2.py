#!/usr/bin/env python3
"""
Generate CRISPRArchitect Lab Meeting Presentation (v2)
=====================================================
Generates a professional 20-slide PPTX for presenting CRISPRArchitect
to PI Debojyoti Chakraborty and labmates at CSIR-IGIB.

Usage:
    python3 generate_presentation_v2.py

Output:
    CRISPRArchitect_v2_LabMeeting.pptx
"""

try:
    from pptx import Presentation
    from pptx.util import Inches, Pt, Emu
    from pptx.dml.color import RGBColor
    from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
    from pptx.enum.shapes import MSO_SHAPE
    print("python-pptx available")
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
# Color palette
# ===========================================================================
BG_NAVY    = RGBColor(0x1B, 0x28, 0x38)
WHITE      = RGBColor(0xFF, 0xFF, 0xFF)
TEAL       = RGBColor(0x1B, 0x9E, 0x77)
CORAL      = RGBColor(0xD9, 0x5F, 0x02)
GOLD       = RGBColor(0xFF, 0xC1, 0x07)
LIGHT_GRAY = RGBColor(0xCC, 0xCC, 0xCC)
DARK_GRAY  = RGBColor(0x88, 0x99, 0xAA)
SOFT_WHITE = RGBColor(0xE8, 0xEE, 0xF4)
RED_ACCENT = RGBColor(0xE7, 0x4C, 0x3C)
GREEN_ACC  = RGBColor(0x2E, 0xCC, 0x71)
BG_DARKER  = RGBColor(0x13, 0x1C, 0x28)
BG_CARD    = RGBColor(0x23, 0x34, 0x48)

# Slide dimensions (16:9)
SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)


# ===========================================================================
# Helper functions
# ===========================================================================
def set_slide_bg(slide, color=BG_NAVY):
    """Set the background fill of a slide."""
    background = slide.background
    fill = background.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_textbox(slide, left, top, width, height, text, font_size=18,
                color=WHITE, bold=False, alignment=PP_ALIGN.LEFT,
                font_name="Calibri"):
    """Add a text box with specified formatting."""
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


def add_bullet_list(slide, left, top, width, height, items, font_size=18,
                    color=WHITE, bullet_color=TEAL, spacing=Pt(6),
                    font_name="Calibri", line_spacing=1.3):
    """Add a bulleted list to a slide."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True

    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()

        # Handle sub-items (indented with "  - ")
        if isinstance(item, tuple):
            text, item_color, item_bold = item
        else:
            text = item
            item_color = color
            item_bold = False

        p.space_after = spacing
        p.line_spacing = line_spacing

        # Bullet character
        run_bullet = p.add_run()
        run_bullet.text = "  \u25B8  "
        run_bullet.font.size = Pt(font_size)
        run_bullet.font.color.rgb = bullet_color
        run_bullet.font.name = font_name
        run_bullet.font.bold = True

        # Text
        run_text = p.add_run()
        run_text.text = text
        run_text.font.size = Pt(font_size)
        run_text.font.color.rgb = item_color
        run_text.font.name = font_name
        run_text.font.bold = item_bold

    return txBox


def add_accent_line(slide, left, top, width, color=TEAL, height=Inches(0.04)):
    """Add a horizontal accent line."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, left, top, width, height
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_card(slide, left, top, width, height, fill_color=BG_CARD):
    """Add a rounded-rectangle card background."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    shape.line.fill.background()
    # Bring to back -- we'll add text on top
    return shape


def add_slide_number(slide, num, total=20):
    """Add a slide number in bottom-right corner."""
    add_textbox(
        slide, Inches(12.0), Inches(7.05), Inches(1.2), Inches(0.4),
        f"{num}/{total}", font_size=10, color=DARK_GRAY,
        alignment=PP_ALIGN.RIGHT
    )


def add_title_bar(slide, title_text, subtitle_text=None):
    """Add a consistent title bar at the top of content slides."""
    add_textbox(
        slide, Inches(0.7), Inches(0.3), Inches(11.5), Inches(0.7),
        title_text, font_size=30, color=WHITE, bold=True
    )
    add_accent_line(slide, Inches(0.7), Inches(0.95), Inches(2.5), TEAL)
    if subtitle_text:
        add_textbox(
            slide, Inches(0.7), Inches(1.05), Inches(11.5), Inches(0.5),
            subtitle_text, font_size=14, color=DARK_GRAY
        )


def add_stat_box(slide, left, top, width, height, value, label,
                 value_color=TEAL, bg_color=BG_CARD):
    """Add a stat highlight box."""
    add_card(slide, left, top, width, height, bg_color)
    add_textbox(
        slide, left, top + Inches(0.15), width, Inches(0.6),
        value, font_size=32, color=value_color, bold=True,
        alignment=PP_ALIGN.CENTER
    )
    add_textbox(
        slide, left, top + Inches(0.7), width, Inches(0.4),
        label, font_size=12, color=LIGHT_GRAY,
        alignment=PP_ALIGN.CENTER
    )


# ===========================================================================
# Build presentation
# ===========================================================================
prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H

# Use blank layout
blank_layout = prs.slide_layouts[6]  # blank


# ---- SLIDE 1: Title -------------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, BG_DARKER)

# Decorative top bar
add_accent_line(slide, Inches(0), Inches(0), SLIDE_W, TEAL, Inches(0.06))

# Main title
add_textbox(
    slide, Inches(1.0), Inches(1.5), Inches(11.3), Inches(1.5),
    "CRISPRArchitect", font_size=48, color=WHITE, bold=True
)
add_textbox(
    slide, Inches(1.0), Inches(2.5), Inches(11.3), Inches(1.2),
    "Transcript-Aware and Consequence-Guided\nDesign of Genome Editing Strategies",
    font_size=24, color=TEAL, bold=False
)

add_accent_line(slide, Inches(1.0), Inches(3.8), Inches(3.0), GOLD, Inches(0.03))

add_textbox(
    slide, Inches(1.0), Inches(4.1), Inches(11.3), Inches(0.5),
    "Vishal Bharti  |  Debojyoti Chakraborty Lab",
    font_size=20, color=SOFT_WHITE, bold=False
)
add_textbox(
    slide, Inches(1.0), Inches(4.65), Inches(11.3), Inches(0.5),
    "CSIR-Institute of Genomics and Integrative Biology, New Delhi",
    font_size=16, color=DARK_GRAY
)
add_textbox(
    slide, Inches(1.0), Inches(5.1), Inches(11.3), Inches(0.5),
    "Lab Meeting  |  March 2026",
    font_size=16, color=DARK_GRAY
)

# Version badges
add_card(slide, Inches(9.5), Inches(5.5), Inches(1.5), Inches(0.5), BG_CARD)
add_textbox(
    slide, Inches(9.5), Inches(5.55), Inches(1.5), Inches(0.4),
    "v1 + v2", font_size=16, color=TEAL, bold=True,
    alignment=PP_ALIGN.CENTER
)

add_slide_number(slide, 1)


# ---- SLIDE 2: The Problem -------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "The Problem", "Why do we need a unified editing framework?")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.5), Inches(5.0), [
    "Multiple editing modalities (BE, PE, HDR) each with distinct constraints",
    "No unified framework for systematic comparison across strategies",
    "Researchers rely on heuristic, ad hoc decision-making",
    "Risk of selecting suboptimal or infeasible strategies",
    "Feasibility depends on local sequence context (PAM, editing window)",
    "Clinical applications demand explicit rejection reasoning",
], font_size=17, line_spacing=1.4)

# Right side: illustrative cards
add_card(slide, Inches(7.0), Inches(1.6), Inches(5.5), Inches(1.2))
add_textbox(slide, Inches(7.2), Inches(1.7), Inches(5.1), Inches(0.4),
            "Base Editing", font_size=16, color=TEAL, bold=True)
add_textbox(slide, Inches(7.2), Inches(2.1), Inches(5.1), Inches(0.5),
            "Restricted to specific transitions within a narrow window",
            font_size=13, color=LIGHT_GRAY)

add_card(slide, Inches(7.0), Inches(3.1), Inches(5.5), Inches(1.2))
add_textbox(slide, Inches(7.2), Inches(3.2), Inches(5.1), Inches(0.4),
            "Prime Editing", font_size=16, color=CORAL, bold=True)
add_textbox(slide, Inches(7.2), Inches(3.6), Inches(5.1), Inches(0.5),
            "Versatile but complex pegRNA design requirements",
            font_size=13, color=LIGHT_GRAY)

add_card(slide, Inches(7.0), Inches(4.6), Inches(5.5), Inches(1.2))
add_textbox(slide, Inches(7.2), Inches(4.7), Inches(5.1), Inches(0.4),
            "HDR", font_size=16, color=GOLD, bold=True)
add_textbox(slide, Inches(7.2), Inches(5.1), Inches(5.1), Inches(0.5),
            "DSB-dependent, variable efficiency, donor design matters",
            font_size=13, color=LIGHT_GRAY)

add_slide_number(slide, 2)


# ---- SLIDE 3: What is CRISPRArchitect? ------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "What is CRISPRArchitect?",
              "A computational decision framework for genome editing")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.8), Inches(4.5), [
    "Unified evaluation of BE + PE + HDR for any point mutation",
    "Transcript-aware: maps genomic variants to CDS, codons, amino acids",
    "Consequence-guided: penalizes splice damage, bystander nonsense",
    "Explicit rejection of infeasible strategies with reasoning",
    "PAM-verified feasibility (not just mutation-type heuristics)",
    "Supports SpCas9 (NGG) and enFnCas9 (NRG)",
], font_size=17, line_spacing=1.4)

# Key differentiator box
add_card(slide, Inches(7.2), Inches(1.8), Inches(5.3), Inches(3.5))
add_textbox(slide, Inches(7.5), Inches(1.95), Inches(4.8), Inches(0.5),
            "Key Differentiators", font_size=18, color=GOLD, bold=True)
add_accent_line(slide, Inches(7.5), Inches(2.45), Inches(2.0), GOLD, Inches(0.02))

items_diff = [
    ("Not just \"which tool\" -- but \"is it even feasible?\"", WHITE),
    ("Transcript context prevents silent assumption errors", WHITE),
    ("Bystander edits evaluated for coding consequence", WHITE),
    ("Scoring weights tuned for iPSC safety", WHITE),
    ("30-case ClinVar benchmark, zero fabrication", TEAL),
]
y_pos = Inches(2.65)
for text, clr in items_diff:
    add_textbox(slide, Inches(7.5), y_pos, Inches(4.8), Inches(0.45),
                f"\u25B8  {text}", font_size=13, color=clr)
    y_pos += Inches(0.45)

add_slide_number(slide, 3)


# ---- SLIDE 4: v1 Foundation (6 modules) -----------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v1 Foundation: 6 Modules",
              "24,000 lines of code | Validated against 4 datasets | 71% concordance with 14 papers")

modules = [
    ("ConversionSim", "HDR gene conversion tract\nsimulation & cssDNA modeling", TEAL),
    ("MOSAIC", "Strategy enumeration engine\n71% concordance, 14 papers", CORAL),
    ("TopoPred", "cssDNA secondary structure\nprediction", GOLD),
    ("ChromBridge", "3D chromatin distance\nestimation", TEAL),
    ("LoopSim", "Cohesin loop extrusion\nsimulation", CORAL),
    ("WebApp", "Interactive Streamlit\ninterface", GOLD),
]

for i, (name, desc, accent) in enumerate(modules):
    col = i % 3
    row = i // 3
    left = Inches(0.7) + Inches(col * 4.1)
    top = Inches(1.6) + Inches(row * 2.6)
    add_card(slide, left, top, Inches(3.7), Inches(2.2))
    add_accent_line(slide, left + Inches(0.15), top + Inches(0.1),
                    Inches(0.8), accent, Inches(0.04))
    add_textbox(slide, left + Inches(0.2), top + Inches(0.25),
                Inches(3.3), Inches(0.5),
                name, font_size=20, color=accent, bold=True)
    add_textbox(slide, left + Inches(0.2), top + Inches(0.85),
                Inches(3.3), Inches(1.1),
                desc, font_size=14, color=LIGHT_GRAY)

add_slide_number(slide, 4)


# ---- SLIDE 5: v1 Validation Highlights ------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v1 Validation Highlights",
              "Quantitative agreement with published experimental data")

validations = [
    ("cssDNA 2.07x advantage", "Iyer et al. 2022 reported 1.9x", TEAL),
    ("Stagger 1.82x enhancement", "Chauhan et al. 2023 reported 1.9x", TEAL),
    ("MOSAIC 100% accuracy", "Base editing cases fully concordant", CORAL),
    ("4 experimental datasets", "Independent validation sources", GOLD),
]

for i, (metric, reference, accent) in enumerate(validations):
    top = Inches(1.6) + Inches(i * 1.3)
    add_card(slide, Inches(0.7), top, Inches(11.8), Inches(1.1))
    add_textbox(slide, Inches(1.0), top + Inches(0.1), Inches(5.0), Inches(0.5),
                metric, font_size=20, color=accent, bold=True)
    add_textbox(slide, Inches(1.0), top + Inches(0.55), Inches(5.0), Inches(0.4),
                reference, font_size=14, color=LIGHT_GRAY)
    # Checkmark
    add_textbox(slide, Inches(11.0), top + Inches(0.2), Inches(1.0), Inches(0.5),
                "\u2713", font_size=28, color=GREEN_ACC, bold=True,
                alignment=PP_ALIGN.CENTER)

add_textbox(slide, Inches(0.7), Inches(6.5), Inches(11.0), Inches(0.5),
            "Manuscript scope: PLOS Computational Biology (in preparation)",
            font_size=14, color=DARK_GRAY)

add_slide_number(slide, 5)


# ---- SLIDE 6: Why v2? What was missing? -----------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Why v2? What Was Missing?",
              "Critical gaps that motivated the transcript-aware extension")

gaps = [
    ("No transcript-level context",
     "v1 operated on genomic coordinates without mapping to transcripts, CDS, or codons"),
    ("No PAM-verified feasibility",
     "Strategy suggestions assumed feasibility without checking actual guide placement"),
    ("No coding consequence awareness",
     "Could not distinguish synonymous from nonsense bystander edits"),
    ("No explicit rejection reasoning",
     "Infeasible strategies were silently omitted instead of explained"),
    ("No benchmark with real variants",
     "Validation used synthetic examples, not curated ClinVar pathogenic variants"),
]

for i, (gap, detail) in enumerate(gaps):
    top = Inches(1.5) + Inches(i * 1.1)
    # Red X marker
    add_textbox(slide, Inches(0.7), top, Inches(0.5), Inches(0.5),
                "\u2717", font_size=22, color=RED_ACCENT, bold=True)
    add_textbox(slide, Inches(1.3), top, Inches(11.0), Inches(0.45),
                gap, font_size=18, color=WHITE, bold=True)
    add_textbox(slide, Inches(1.3), top + Inches(0.45), Inches(11.0), Inches(0.4),
                detail, font_size=13, color=DARK_GRAY)

add_slide_number(slide, 6)


# ---- SLIDE 7: v2 Architecture (Pipeline) ----------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "v2 Architecture: 7-Stage Pipeline",
              "End-to-end from variant input to ranked strategies")

stages = [
    ("1", "Input\nParsing", TEAL),
    ("2", "Transcript\nMapping", CORAL),
    ("3", "Reference\nValidation", GOLD),
    ("4", "Coding\nAnnotation", TEAL),
    ("5", "Feasibility\n(BE/PE/HDR)", CORAL),
    ("6", "Strategy\nGeneration", GOLD),
    ("7", "Scoring &\nRanking", TEAL),
]

box_w = Inches(1.45)
box_h = Inches(1.8)
start_x = Inches(0.55)
y_top = Inches(2.2)

for i, (num, label, accent) in enumerate(stages):
    x = start_x + Inches(i * 1.8)
    add_card(slide, x, y_top, box_w, box_h, BG_CARD)
    # Number circle-like
    add_textbox(slide, x, y_top + Inches(0.15), box_w, Inches(0.5),
                num, font_size=28, color=accent, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, x, y_top + Inches(0.65), box_w, Inches(0.9),
                label, font_size=13, color=LIGHT_GRAY,
                alignment=PP_ALIGN.CENTER)
    # Arrow between boxes
    if i < len(stages) - 1:
        arrow_x = x + box_w + Inches(0.05)
        add_textbox(slide, arrow_x, y_top + Inches(0.6), Inches(0.3), Inches(0.4),
                    "\u25B6", font_size=16, color=DARK_GRAY,
                    alignment=PP_ALIGN.CENTER)

# Bottom annotation
add_textbox(slide, Inches(0.7), Inches(4.6), Inches(12.0), Inches(0.5),
            "Ensembl REST API (GRCh38)  |  Real-time PAM scanning  |  "
            "Bystander consequence evaluation  |  Multi-nuclease support",
            font_size=13, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Pipeline detail cards
add_card(slide, Inches(0.7), Inches(5.2), Inches(5.8), Inches(1.8))
add_textbox(slide, Inches(0.9), Inches(5.3), Inches(5.4), Inches(0.4),
            "Input \u2192 Annotation", font_size=16, color=TEAL, bold=True)
add_textbox(slide, Inches(0.9), Inches(5.7), Inches(5.4), Inches(1.1),
            "chr:pos:ref>alt \u2192 gene lookup \u2192 transcript selection \u2192\n"
            "CDS mapping \u2192 codon identification \u2192 amino acid change \u2192\n"
            "consequence classification (syn/mis/non/splice/fs)",
            font_size=12, color=LIGHT_GRAY)

add_card(slide, Inches(6.8), Inches(5.2), Inches(5.8), Inches(1.8))
add_textbox(slide, Inches(7.0), Inches(5.3), Inches(5.4), Inches(0.4),
            "Feasibility \u2192 Ranking", font_size=16, color=CORAL, bold=True)
add_textbox(slide, Inches(7.0), Inches(5.7), Inches(5.4), Inches(1.1),
            "PAM scanning (NGG/NRG) \u2192 editing window check \u2192\n"
            "bystander enumeration \u2192 pegRNA/donor design \u2192\n"
            "weighted scoring \u2192 ranked output with rejection reasons",
            font_size=12, color=LIGHT_GRAY)

add_slide_number(slide, 7)


# ---- SLIDE 8: Transcript-Aware Mapping ------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Transcript-Aware Mapping",
              "Ensembl REST API integration for GRCh38 coordinate resolution")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.8), Inches(5.0), [
    "Genomic position \u2192 overlapping gene(s) via Ensembl REST",
    "Gene \u2192 canonical transcript (MANE Select preferred)",
    "Transcript \u2192 exon boundaries, CDS start/end, strand",
    "CDS coordinate \u2192 codon position (1st, 2nd, or 3rd base)",
    "Codon \u2192 reference & mutant amino acid (standard genetic code)",
    "Forward and reverse strand handling verified",
    "Reference allele validation against Ensembl genome sequence",
], font_size=16, line_spacing=1.35)

# Example box
add_card(slide, Inches(7.0), Inches(1.6), Inches(5.5), Inches(4.5))
add_textbox(slide, Inches(7.2), Inches(1.7), Inches(5.1), Inches(0.4),
            "Example: NF1 c.910C>T", font_size=18, color=GOLD, bold=True)
add_accent_line(slide, Inches(7.2), Inches(2.15), Inches(2.0), GOLD, Inches(0.02))

mapping_steps = [
    "chr17:31,232,517  C>T",
    "\u2193",
    "Gene: NF1 (ENSG00000196712)",
    "Transcript: ENST00000358273",
    "Strand: forward (+)",
    "CDS position: 910",
    "Codon: CGA \u2192 TGA",
    "Amino acid: Arg304 \u2192 Ter (Stop)",
    "\u2193",
    "Consequence: NONSENSE",
]
y = Inches(2.35)
for step in mapping_steps:
    if step == "\u2193":
        clr = DARK_GRAY
        sz = 14
    elif "NONSENSE" in step:
        clr = RED_ACCENT
        sz = 15
    elif "\u2192" in step and "Ter" in step:
        clr = CORAL
        sz = 14
    else:
        clr = LIGHT_GRAY
        sz = 13
    add_textbox(slide, Inches(7.4), y, Inches(4.8), Inches(0.35),
                step, font_size=sz, color=clr, alignment=PP_ALIGN.LEFT)
    y += Inches(0.35)

add_slide_number(slide, 8)


# ---- SLIDE 9: Coding & Splice Annotation ----------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Coding & Splice Annotation",
              "Consequence classification following ACMG standards")

# Consequence types
consequences = [
    ("Synonymous", "No amino acid change; silent", TEAL),
    ("Missense", "Different amino acid; functional impact varies", GOLD),
    ("Nonsense", "Premature stop codon; likely loss-of-function", RED_ACCENT),
    ("Splice site", "Within 2bp of exon boundary (ACMG canonical)", CORAL),
    ("Splice region", "Within 3-8bp of exon boundary", CORAL),
    ("Frameshift", "Indel not divisible by 3; reading frame disrupted", RED_ACCENT),
]

for i, (ctype, desc, accent) in enumerate(consequences):
    col = i % 2
    row = i // 2
    left = Inches(0.7) + Inches(col * 6.3)
    top = Inches(1.5) + Inches(row * 1.6)
    add_card(slide, left, top, Inches(5.9), Inches(1.3))
    add_textbox(slide, left + Inches(0.2), top + Inches(0.1),
                Inches(5.5), Inches(0.4),
                ctype, font_size=18, color=accent, bold=True)
    add_textbox(slide, left + Inches(0.2), top + Inches(0.55),
                Inches(5.5), Inches(0.6),
                desc, font_size=13, color=LIGHT_GRAY)

# Bottom note
add_textbox(slide, Inches(0.7), Inches(6.3), Inches(12.0), Inches(0.5),
            "HGVS notation generated for all variants  |  "
            "Splice proximity scored by distance to nearest exon boundary",
            font_size=13, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

add_slide_number(slide, 9)


# ---- SLIDE 10: Feasibility Engines ----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Feasibility Engines",
              "PAM-verified strategy assessment for each editing modality")

# Three columns for BE, PE, HDR
engines = [
    ("Base Editing", TEAL, [
        "ABE window: pos 4-7",
        "CBE window: pos 4-8",
        "PAM-verified guide search",
        "Bystander edit counting",
        "Consequence of each bystander",
        "SpCas9 (NGG) + enFnCas9 (NRG)",
    ]),
    ("Prime Editing", CORAL, [
        "pegRNA design automation",
        "PBS length: 13 nt (default)",
        "RT template: 10-30 nt",
        "PE3 nicking guide search",
        "Nick-to-edit distance scoring",
        "All mutation types supported",
    ]),
    ("HDR", GOLD, [
        "Cut-to-edit distance scoring",
        "Donor type recommendation",
        "ssODN vs. plasmid selection",
        "Conversion probability estimate",
        "Homology arm optimization",
        "DSB-dependent (risk noted)",
    ]),
]

for i, (name, accent, items) in enumerate(engines):
    left = Inches(0.5) + Inches(i * 4.2)
    add_card(slide, left, Inches(1.5), Inches(3.9), Inches(5.3))
    add_textbox(slide, left + Inches(0.15), Inches(1.6), Inches(3.6), Inches(0.5),
                name, font_size=22, color=accent, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_accent_line(slide, left + Inches(0.8), Inches(2.15),
                    Inches(2.3), accent, Inches(0.03))
    y = Inches(2.35)
    for item in items:
        add_textbox(slide, left + Inches(0.2), y, Inches(3.5), Inches(0.42),
                    f"\u25B8  {item}", font_size=13, color=LIGHT_GRAY)
        y += Inches(0.42)

add_slide_number(slide, 10)


# ---- SLIDE 11: Scoring Function -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Scoring Function",
              "Weighted multi-objective scoring tuned for iPSC safety")

# Formula
add_card(slide, Inches(0.7), Inches(1.5), Inches(11.8), Inches(1.2))
add_textbox(
    slide, Inches(1.0), Inches(1.65), Inches(11.2), Inches(0.8),
    "Score = w1\u00B7Safety + w2\u00B7Feasibility - w3\u00B7Complexity - w4\u00B7Risk + w5\u00B7Confidence",
    font_size=22, color=WHITE, bold=True, alignment=PP_ALIGN.CENTER
)

# Weight table
weights = [
    ("Safety", "0.30", "Highest weight -- iPSC context", TEAL),
    ("Feasibility", "0.25", "PAM availability, window placement", CORAL),
    ("Complexity", "0.20", "pegRNA design, donor length", GOLD),
    ("Risk", "0.15", "DSB dependence, off-target potential", RED_ACCENT),
    ("Confidence", "0.10", "Literature support, validated parameters", LIGHT_GRAY),
]

add_card(slide, Inches(0.7), Inches(3.0), Inches(5.5), Inches(4.0))
add_textbox(slide, Inches(0.9), Inches(3.1), Inches(5.1), Inches(0.4),
            "iPSC Weight Configuration", font_size=16, color=TEAL, bold=True)

y = Inches(3.55)
for name, weight, desc, accent in weights:
    add_textbox(slide, Inches(0.9), y, Inches(1.5), Inches(0.35),
                name, font_size=14, color=accent, bold=True)
    add_textbox(slide, Inches(2.5), y, Inches(0.8), Inches(0.35),
                weight, font_size=14, color=WHITE, bold=True)
    add_textbox(slide, Inches(3.4), y, Inches(2.7), Inches(0.35),
                desc, font_size=11, color=DARK_GRAY)
    y += Inches(0.55)

# Consequence penalties
add_card(slide, Inches(6.5), Inches(3.0), Inches(6.0), Inches(4.0))
add_textbox(slide, Inches(6.7), Inches(3.1), Inches(5.6), Inches(0.4),
            "Consequence-Based Penalties", font_size=16, color=CORAL, bold=True)

penalties = [
    ("Splice site disruption", "-0.15", RED_ACCENT),
    ("Bystander nonsense", "-0.25", RED_ACCENT),
    ("Bystander missense", "-0.10", CORAL),
    ("DSB-free bonus (BE/PE)", "+0.10", GREEN_ACC),
    ("Clean edit (no bystander)", "+0.05", GREEN_ACC),
]

y = Inches(3.55)
for name, value, accent in penalties:
    add_textbox(slide, Inches(6.7), y, Inches(3.5), Inches(0.35),
                name, font_size=14, color=LIGHT_GRAY)
    add_textbox(slide, Inches(10.5), y, Inches(1.5), Inches(0.35),
                value, font_size=14, color=accent, bold=True)
    y += Inches(0.55)

add_textbox(slide, Inches(6.7), Inches(6.4), Inches(5.6), Inches(0.4),
            "DSB-free strategies always score higher in iPSC context",
            font_size=12, color=DARK_GRAY)

add_slide_number(slide, 11)


# ---- SLIDE 12: Benchmark Design -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Benchmark Design",
              "30 curated ClinVar cases with tiered truth labels")

# Stats row
stats = [
    ("30", "ClinVar Cases", TEAL),
    ("11", "Mutation Categories", CORAL),
    ("3-Tier", "Truth Labels", GOLD),
    ("GRCh38", "Verified Coords", WHITE),
]
for i, (val, lbl, accent) in enumerate(stats):
    left = Inches(0.7) + Inches(i * 3.1)
    add_stat_box(slide, left, Inches(1.5), Inches(2.7), Inches(1.2),
                 val, lbl, accent)

# Categories
add_card(slide, Inches(0.7), Inches(3.1), Inches(11.8), Inches(3.6))
add_textbox(slide, Inches(0.9), Inches(3.2), Inches(11.4), Inches(0.4),
            "11 Mutation Categories", font_size=16, color=TEAL, bold=True)

categories = [
    "C>T / G>A transitions (CBE-compatible)",
    "A>G / T>C transitions (ABE-compatible)",
    "C>G / C>A transversions (PE required)",
    "Small insertions (1-3 bp)",
    "Small deletions (1-5 bp)",
    "Large deletions (>20 bp)",
    "Splice site variants",
    "Nonsense mutations",
    "Multi-nucleotide variants",
    "Compound heterozygous",
    "Edge cases (non-coding, PAM-less)",
]

for i, cat in enumerate(categories):
    col = i % 3
    row = i // 3
    left = Inches(1.0) + Inches(col * 4.0)
    top = Inches(3.7) + Inches(row * 0.55)
    add_textbox(slide, left, top, Inches(3.8), Inches(0.4),
                f"\u25B8  {cat}", font_size=12, color=LIGHT_GRAY)

# Truth label explanation
add_textbox(slide, Inches(0.7), Inches(6.3), Inches(11.8), Inches(0.5),
            "Truth labels: Preferred (gold standard) | Acceptable (valid alternative) | "
            "Reject (infeasible)   --   All coordinates verified against Ensembl, zero fabrication",
            font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

add_slide_number(slide, 12)


# ---- SLIDE 13: Results -- Overall Performance ------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Results: Overall Performance",
              "30/30 cases executed successfully with zero errors")

# Big stat boxes
metrics = [
    ("86.7%", "Top-1 Accuracy", "26 / 30", TEAL),
    ("96.7%", "Top-3 Accuracy", "29 / 30", CORAL),
    ("90.0%", "Rejection Accuracy", "27 / 30", GOLD),
    ("100%", "Execution Success", "30 / 30", GREEN_ACC),
]

for i, (pct, label, detail, accent) in enumerate(metrics):
    left = Inches(0.5) + Inches(i * 3.2)
    add_card(slide, left, Inches(1.7), Inches(2.9), Inches(3.0))
    add_textbox(slide, left, Inches(1.9), Inches(2.9), Inches(0.8),
                pct, font_size=42, color=accent, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, left, Inches(2.8), Inches(2.9), Inches(0.5),
                label, font_size=16, color=WHITE, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, left, Inches(3.4), Inches(2.9), Inches(0.4),
                detail, font_size=14, color=DARK_GRAY,
                alignment=PP_ALIGN.CENTER)

# Strategy distribution
add_card(slide, Inches(0.7), Inches(5.1), Inches(11.8), Inches(1.8))
add_textbox(slide, Inches(0.9), Inches(5.2), Inches(5.0), Inches(0.4),
            "Top-Ranked Strategy Distribution", font_size=16, color=WHITE, bold=True)

# PE bar (visual)
pe_width = Inches(8.5)  # 29/30
add_card(slide, Inches(1.5), Inches(5.7), pe_width, Inches(0.4), CORAL)
add_textbox(slide, Inches(1.5), Inches(5.7), pe_width, Inches(0.4),
            "PE: 29 / 30", font_size=13, color=WHITE, bold=True,
            alignment=PP_ALIGN.CENTER)

# HDR bar
hdr_width = Inches(0.3)
add_card(slide, Inches(1.5), Inches(6.2), hdr_width, Inches(0.4), GOLD)
add_textbox(slide, Inches(2.0), Inches(6.2), Inches(2.0), Inches(0.4),
            "HDR: 1 / 30", font_size=13, color=GOLD, bold=True)

# BE bar
add_textbox(slide, Inches(5.5), Inches(6.2), Inches(4.0), Inches(0.4),
            "BE: 0 / 30  (PAM window bottleneck -- see Slide 15)",
            font_size=13, color=RED_ACCENT, bold=True)

add_slide_number(slide, 13)


# ---- SLIDE 14: Results -- Per-Category Accuracy ----------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Results: Per-Category Accuracy",
              "Category-level breakdown of Top-1 and Top-3 performance")

# Table-like layout
headers = ["Category", "N", "Top-1", "Top-3", "Notes"]
header_widths = [Inches(3.2), Inches(0.6), Inches(1.2), Inches(1.2), Inches(5.5)]
header_x = [Inches(0.7)]
for w in header_widths[:-1]:
    header_x.append(header_x[-1] + w)

# Header row
y_h = Inches(1.5)
for j, (hdr, x_pos, w) in enumerate(zip(headers, header_x, header_widths)):
    add_textbox(slide, x_pos, y_h, w, Inches(0.4),
                hdr, font_size=14, color=TEAL, bold=True)
add_accent_line(slide, Inches(0.7), Inches(1.85), Inches(11.8), TEAL, Inches(0.02))

rows = [
    ("Clean BE (CBE/ABE)", "4", "100%", "100%", "All correct, but none top-ranked as BE"),
    ("PE transversions", "3", "100%", "100%", "PE correctly top-ranked"),
    ("PE small insertions", "3", "100%", "100%", "PE correctly top-ranked"),
    ("PE small deletions", "3", "100%", "100%", "PE correctly top-ranked"),
    ("Compound heterozygous", "2", "100%", "100%", "Both alleles handled"),
    ("Edge cases", "3", "100%", "100%", "Non-coding, PAM-less correctly rejected"),
    ("Nonsense mutations", "3", "67%", "100%", "1 case: PE ranked over HDR"),
    ("Splice variants", "3", "67%", "100%", "Splice-aware penalties applied"),
    ("Multi-nucleotide", "3", "67%", "100%", "Complex variants handled"),
    ("Large deletions", "3", "0%", "100%", "HDR preferred but PE ranked higher"),
]

y = Inches(1.95)
for row_data in rows:
    cat, n, t1, t3, notes = row_data
    t1_color = GREEN_ACC if t1 == "100%" else (GOLD if t1 == "67%" else RED_ACCENT)
    add_textbox(slide, header_x[0], y, header_widths[0], Inches(0.35),
                cat, font_size=12, color=LIGHT_GRAY)
    add_textbox(slide, header_x[1], y, header_widths[1], Inches(0.35),
                n, font_size=12, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_textbox(slide, header_x[2], y, header_widths[2], Inches(0.35),
                t1, font_size=12, color=t1_color, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, header_x[3], y, header_widths[3], Inches(0.35),
                t3, font_size=12, color=GREEN_ACC, bold=True,
                alignment=PP_ALIGN.CENTER)
    add_textbox(slide, header_x[4], y, header_widths[4], Inches(0.35),
                notes, font_size=11, color=DARK_GRAY)
    y += Inches(0.45)

add_slide_number(slide, 14)


# ---- SLIDE 15: KEY FINDING -- PAM Window Bottleneck -----------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, RGBColor(0x15, 0x1E, 0x2B))  # slightly different bg for emphasis

# IMPORTANT banner
add_card(slide, Inches(0.0), Inches(0.0), SLIDE_W, Inches(0.5), CORAL)
add_textbox(slide, Inches(0.0), Inches(0.05), SLIDE_W, Inches(0.4),
            "KEY FINDING", font_size=18, color=WHITE, bold=True,
            alignment=PP_ALIGN.CENTER)

add_title_bar(slide, "PAM-Editing Window Is the Real Bottleneck",
              "The most important result of this study")

# Main finding box
add_card(slide, Inches(0.7), Inches(1.7), Inches(11.8), Inches(2.2),
         RGBColor(0x2A, 0x15, 0x10))
add_textbox(
    slide, Inches(1.0), Inches(1.9), Inches(11.2), Inches(0.7),
    "At ALL 7 ClinVar loci with ABE-compatible transitions (A>G / T>C),\n"
    "NO SpCas9 guide placed the target base at ABE positions 4-7.",
    font_size=20, color=CORAL, bold=True, alignment=PP_ALIGN.CENTER
)
add_textbox(
    slide, Inches(1.0), Inches(2.8), Inches(11.2), Inches(0.7),
    "The PAM-dependent editing window is MORE restrictive than mutation type.\n"
    "Mutation-type classification alone is INSUFFICIENT for strategy selection.",
    font_size=17, color=WHITE, alignment=PP_ALIGN.CENTER
)

# Implications
add_card(slide, Inches(0.7), Inches(4.2), Inches(5.6), Inches(2.8))
add_textbox(slide, Inches(0.9), Inches(4.3), Inches(5.2), Inches(0.4),
            "Implications", font_size=18, color=TEAL, bold=True)
add_bullet_list(slide, Inches(0.9), Inches(4.8), Inches(5.2), Inches(2.0), [
    "Locus-specific PAM scanning is essential",
    "Mutation classification alone is misleading",
    "PE dominates because it has no editing window",
    "This explains PE = 29/30 top-ranked",
], font_size=14, line_spacing=1.4, bullet_color=CORAL)

# Why this matters
add_card(slide, Inches(6.6), Inches(4.2), Inches(5.9), Inches(2.8))
add_textbox(slide, Inches(6.8), Inches(4.3), Inches(5.5), Inches(0.4),
            "Why This Matters", font_size=18, color=GOLD, bold=True)
add_bullet_list(slide, Inches(6.8), Inches(4.8), Inches(5.5), Inches(2.0), [
    "Most tools classify by mutation type only",
    "A>G does NOT mean ABE will work",
    "PAM context must be checked in silico first",
    "CRISPRArchitect does this automatically",
], font_size=14, line_spacing=1.4, bullet_color=GOLD)

add_slide_number(slide, 15)


# ---- SLIDE 16: enFnCas9 Advantage -----------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "enFnCas9 Advantage",
              "Our lab's engineered nuclease with broader PAM compatibility")

add_bullet_list(slide, Inches(0.7), Inches(1.5), Inches(5.8), Inches(2.5), [
    "Developed in our lab (Chakraborty et al.)",
    "Engineered FnCas9 with enhanced activity",
    "NRG PAM -- broader than SpCas9 NGG",
    "First-class support in CRISPRArchitect v2",
], font_size=17, line_spacing=1.4)

# PAM comparison
add_card(slide, Inches(0.7), Inches(4.0), Inches(5.5), Inches(2.5))
add_textbox(slide, Inches(0.9), Inches(4.1), Inches(5.1), Inches(0.4),
            "PAM Comparison", font_size=16, color=TEAL, bold=True)
add_textbox(slide, Inches(0.9), Inches(4.6), Inches(5.1), Inches(0.4),
            "SpCas9:     NGG  (16 possible dinucleotides)", font_size=14, color=LIGHT_GRAY)
add_textbox(slide, Inches(0.9), Inches(5.05), Inches(5.1), Inches(0.4),
            "enFnCas9:  NRG  (8 possible trinucleotides, 2x coverage)",
            font_size=14, color=TEAL, bold=True)
add_textbox(slide, Inches(0.9), Inches(5.55), Inches(5.1), Inches(0.5),
            "Broader PAM = more guides = better chance of\nplacing target in editing window",
            font_size=13, color=DARK_GRAY)

# Potential impact
add_card(slide, Inches(6.6), Inches(1.6), Inches(6.0), Inches(5.0))
add_textbox(slide, Inches(6.8), Inches(1.7), Inches(5.6), Inches(0.4),
            "Potential Impact on BE Feasibility", font_size=16, color=CORAL, bold=True)
add_accent_line(slide, Inches(6.8), Inches(2.15), Inches(2.5), CORAL, Inches(0.02))
add_bullet_list(slide, Inches(6.8), Inches(2.3), Inches(5.6), Inches(3.5), [
    "7/7 ABE-compatible loci failed with SpCas9 NGG",
    "NRG PAM could rescue some of these loci",
    "2x PAM coverage doubles potential guide sites",
    "CRISPRArchitect already searches both PAMs",
    "Future: quantify rescue rate across full ClinVar",
], font_size=14, line_spacing=1.4, bullet_color=CORAL)

add_textbox(slide, Inches(6.8), Inches(5.8), Inches(5.6), Inches(0.5),
            "enFnCas9 + CRISPRArchitect: a natural synergy\n"
            "between our wet lab and computational tools",
            font_size=12, color=DARK_GRAY)

add_slide_number(slide, 16)


# ---- SLIDE 17: Technical Summary ------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Technical Summary",
              "Codebase, testing, and reproducibility")

# Stats grid
tech_stats = [
    ("24,000", "v1 LOC", TEAL),
    ("11,002", "v2 LOC", CORAL),
    ("~35,000", "Total LOC", GOLD),
    ("123", "Tests Passing", GREEN_ACC),
    ("6", "v1 Modules", TEAL),
    ("26", "v2 New Files", CORAL),
    ("43", "Total Files", GOLD),
    ("0", "Regressions", GREEN_ACC),
]

for i, (val, lbl, accent) in enumerate(tech_stats):
    col = i % 4
    row = i // 4
    left = Inches(0.5) + Inches(col * 3.15)
    top = Inches(1.6) + Inches(row * 1.4)
    add_stat_box(slide, left, top, Inches(2.85), Inches(1.15), val, lbl, accent)

# Reproducibility features
add_card(slide, Inches(0.7), Inches(4.6), Inches(5.5), Inches(2.5))
add_textbox(slide, Inches(0.9), Inches(4.7), Inches(5.1), Inches(0.4),
            "Reproducibility", font_size=16, color=TEAL, bold=True)
add_bullet_list(slide, Inches(0.9), Inches(5.1), Inches(5.1), Inches(1.8), [
    "GRCh38 coordinates verified against Ensembl",
    "API retry logic with exponential backoff",
    "Seed-based simulations for deterministic output",
    "All benchmark results reproducible",
], font_size=13, line_spacing=1.3, bullet_color=TEAL)

# Tech stack
add_card(slide, Inches(6.5), Inches(4.6), Inches(6.0), Inches(2.5))
add_textbox(slide, Inches(6.7), Inches(4.7), Inches(5.6), Inches(0.4),
            "Technology Stack", font_size=16, color=CORAL, bold=True)
add_bullet_list(slide, Inches(6.7), Inches(5.1), Inches(5.6), Inches(1.8), [
    "Python 3.9+  |  NumPy  |  Pandas",
    "Streamlit (WebApp)  |  Biopython",
    "Ensembl REST API  |  pytest",
    "GitHub: github.com/visvikbharti/CRISPRArchitect",
], font_size=13, line_spacing=1.3, bullet_color=CORAL)

add_slide_number(slide, 17)


# ---- SLIDE 18: Limitations (Honest) ---------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Limitations",
              "Honest assessment of current constraints")

limitations = [
    ("No experimental validation",
     "All results are computational; wet-lab confirmation is pending",
     RED_ACCENT),
    ("Simplified CDS model",
     "Uses canonical transcript only; does not handle alternative splicing",
     CORAL),
    ("No chromatin context",
     "Does not integrate ATAC-seq or chromatin accessibility data",
     CORAL),
    ("No off-target prediction",
     "Guide specificity scoring not yet incorporated",
     CORAL),
    ("Large deletions poorly handled",
     "0/3 Top-1 accuracy; PE over-ranked vs. HDR for >20bp deletions",
     RED_ACCENT),
    ("BE applicability lower than expected",
     "PAM constraints reduce practical BE feasibility significantly",
     GOLD),
]

for i, (title, detail, accent) in enumerate(limitations):
    col = i % 2
    row = i // 2
    left = Inches(0.7) + Inches(col * 6.3)
    top = Inches(1.5) + Inches(row * 1.8)
    add_card(slide, left, top, Inches(5.9), Inches(1.5))
    add_textbox(slide, left + Inches(0.15), top + Inches(0.05),
                Inches(0.4), Inches(0.4),
                "\u26A0", font_size=18, color=accent)
    add_textbox(slide, left + Inches(0.5), top + Inches(0.08),
                Inches(5.2), Inches(0.4),
                title, font_size=16, color=accent, bold=True)
    add_textbox(slide, left + Inches(0.5), top + Inches(0.55),
                Inches(5.2), Inches(0.8),
                detail, font_size=13, color=LIGHT_GRAY)

add_slide_number(slide, 18)


# ---- SLIDE 19: Future Directions -------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide)
add_title_bar(slide, "Future Directions",
              "Roadmap for CRISPRArchitect development")

directions = [
    ("Experimental Validation", "Test top-ranked strategies in iPSC lines",
     "Near-term", TEAL),
    ("Chromatin Accessibility", "Integrate ATAC-seq data for guide efficiency",
     "Near-term", TEAL),
    ("ML-Based Scoring", "Train on experimental outcomes to refine weights",
     "Medium-term", CORAL),
    ("HGVS Parser", "Accept clinical variant notation as input",
     "Medium-term", CORAL),
    ("Off-Target Integration", "Incorporate Cas-OFFinder or similar tools",
     "Medium-term", GOLD),
    ("Clinical Reporting", "Generate structured reports for clinical review",
     "Long-term", GOLD),
]

for i, (title, desc, timeline, accent) in enumerate(directions):
    col = i % 2
    row = i // 2
    left = Inches(0.7) + Inches(col * 6.3)
    top = Inches(1.5) + Inches(row * 1.8)
    add_card(slide, left, top, Inches(5.9), Inches(1.5))
    # Timeline badge
    add_textbox(slide, left + Inches(3.8), top + Inches(0.08),
                Inches(1.9), Inches(0.35),
                timeline, font_size=10, color=accent, bold=True,
                alignment=PP_ALIGN.RIGHT)
    add_textbox(slide, left + Inches(0.2), top + Inches(0.1),
                Inches(3.5), Inches(0.4),
                title, font_size=17, color=accent, bold=True)
    add_textbox(slide, left + Inches(0.2), top + Inches(0.6),
                Inches(5.5), Inches(0.7),
                desc, font_size=14, color=LIGHT_GRAY)

add_slide_number(slide, 19)


# ---- SLIDE 20: Thank You --------------------------------------------------
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, BG_DARKER)

# Decorative top bar
add_accent_line(slide, Inches(0), Inches(0), SLIDE_W, TEAL, Inches(0.06))

add_textbox(
    slide, Inches(1.0), Inches(2.0), Inches(11.3), Inches(1.0),
    "Thank You", font_size=48, color=WHITE, bold=True,
    alignment=PP_ALIGN.CENTER
)

add_accent_line(slide, Inches(5.0), Inches(3.2), Inches(3.3), GOLD, Inches(0.03))

add_textbox(
    slide, Inches(1.0), Inches(3.6), Inches(11.3), Inches(0.6),
    "CRISPRArchitect: Transcript-Aware Genome Editing Strategy Design",
    font_size=18, color=TEAL, alignment=PP_ALIGN.CENTER
)

# Links
add_card(slide, Inches(3.5), Inches(4.5), Inches(6.3), Inches(1.6))
add_textbox(
    slide, Inches(3.5), Inches(4.65), Inches(6.3), Inches(0.5),
    "github.com/visvikbharti/CRISPRArchitect",
    font_size=18, color=WHITE, bold=True, alignment=PP_ALIGN.CENTER
)
add_textbox(
    slide, Inches(3.5), Inches(5.2), Inches(6.3), Inches(0.5),
    "Vishal Bharti  |  Debojyoti Chakraborty Lab  |  CSIR-IGIB",
    font_size=14, color=DARK_GRAY, alignment=PP_ALIGN.CENTER
)

add_textbox(
    slide, Inches(1.0), Inches(6.5), Inches(11.3), Inches(0.5),
    "Questions?",
    font_size=28, color=SOFT_WHITE, bold=True, alignment=PP_ALIGN.CENTER
)

add_slide_number(slide, 20)


# ===========================================================================
# Save
# ===========================================================================
output_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(output_dir, "CRISPRArchitect_v2_LabMeeting.pptx")
prs.save(output_path)
print(f"\nPresentation saved to: {output_path}")
print(f"Total slides: {len(prs.slides)}")
