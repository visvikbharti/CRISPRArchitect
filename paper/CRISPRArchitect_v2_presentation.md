# CRISPRArchitect v2 -- Conference Presentation Outline

**Slide-by-slide outline for PowerPoint conversion**
**Style:** Dark theme (charcoal/navy background, white text, accent colors: teal #00B4D8, coral #FF6B6B, gold #FFD166)
**Duration:** 20--25 minutes + 10 minutes Q&A
**Font:** Montserrat (headings), Inter or Open Sans (body), JetBrains Mono (code/numbers)

---

## SLIDE 1: Title Slide

**Background:** Dark navy (#0D1B2A) with subtle double-helix motif, bottom accent bar in teal

**Title (large, white, centered):**
> CRISPRArchitect: Transcript-Aware and Consequence-Guided Design of Genome Editing Strategies

**Subtitle (teal accent):**
> A unified computational framework for base editing, prime editing, and HDR

**Authors:**
> Vishal Bharti, Debojyoti Chakraborty

**Affiliation:**
> CSIR-Institute of Genomics and Integrative Biology, New Delhi, India
> Academy of Scientific and Innovative Research (AcSIR)

**Bottom bar:** GitHub logo + `github.com/visvikbharti/CRISPRArchitect`

---

## SLIDE 2: The Problem

**Background:** Dark charcoal (#1B2838)
**Layout:** Four problem statements as icon-text blocks arranged in a 2x2 grid

| Icon | Statement |
|------|-----------|
| Split arrows | **Multiple editing modalities exist** -- base editing, prime editing, HDR each have distinct mechanisms |
| Warning triangle | **Each has different constraints and consequences** -- PAM requirements, editing windows, DSB burden, bystander effects |
| Broken chain | **No unified framework for comparison** -- existing tools (CRISPOR, PrimeDesign, BE-Designer) operate in isolation |
| Dice | **Researchers rely on heuristic decision-making** -- manual, ad hoc strategy selection without systematic evaluation |

**Bottom callout box (coral border):**
> "Which editing strategy should I use for *this* variant, in *this* cell type?"
> -- A question no current tool answers systematically.

**Speaker note:** Emphasize that the problem is not guide RNA design -- many excellent tools exist for that. The gap is in cross-modality comparison and consequence-aware ranking.

---

## SLIDE 3: What CRISPRArchitect Does

**Background:** Dark charcoal
**Layout:** Five capability pills arranged vertically with teal left-border accents

1. **Unified strategy space** -- Evaluates BE + PE + HDR + hybrid strategies in a single framework
2. **Transcript-aware variant mapping** -- Maps genomic variants to exon structure, CDS position, codon frame, splice proximity via Ensembl
3. **PAM-verified feasibility checking** -- Scans both strands for SpCas9 (NGG) and enFnCas9 (NRG) PAMs; verifies editing window placement
4. **Consequence-aware scoring** -- Multi-objective function penalizing bystander edits, splice disruption, and DSB burden
5. **Explicit rejection of infeasible strategies** -- Strategies that fail feasibility constraints are filtered, not silently omitted

**Right side:** Simplified before/after schematic:
- Before: researcher with multiple browser tabs open (CRISPOR, Benchling, literature)
- After: single CRISPRArchitect output with ranked, annotated strategies

---

## SLIDE 4: Pipeline Architecture

**Background:** Dark navy
**Layout:** Full-width horizontal pipeline diagram (reference manuscript Figure 1)

**Seven stages as connected nodes (teal-to-gold gradient flow):**

```
[Input]  -->  [Transcript  -->  [Reference   -->  [Coding      -->  [Feasibility  -->  [Strategy     -->  [Scoring
 Variant]      Mapping]          Validation]       Annotation]       Engines]           Generation]        & Ranking]
```

**Under each node, brief annotation:**

| Stage | Key Operation |
|-------|--------------|
| 1. Input | Chromosome, position, ref/alt alleles (GRCh38) |
| 2. Transcript Mapping | Ensembl REST API: canonical transcript, exon ID, CDS coordinate |
| 3. Reference Validation | Fetch genomic sequence; confirm ref allele matches genome |
| 4. Coding Annotation | Classify: synonymous, missense, nonsense, splice-proximal |
| 5. Feasibility Engines | Independent evaluation for BE, PE, HDR (see Slide 6) |
| 6. Strategy Generation | All feasible strategies enumerated; infeasible explicitly rejected |
| 7. Scoring & Ranking | Multi-objective score; top strategies with interpretable reasoning |

**Footer:** "15 Python modules, ~9,400 lines of code (v2), building on ~24,000 LOC (v1)"

---

## SLIDE 5: Transcript-Aware Mapping

**Background:** Dark charcoal
**Layout:** Left: schematic of gene structure. Right: annotation table.

**Left panel -- Gene schematic:**
- Multi-exon gene (forward or reverse strand)
- Highlighted variant position in one exon
- Arrows showing: exon identification, CDS position mapping, codon frame determination
- Splice donor/acceptor zones marked (2 bp = red, 3-8 bp = amber)

**Right panel -- What the system extracts:**

| Feature | Source | Example |
|---------|--------|---------|
| Canonical transcript | Ensembl VEP | ENST00000257430 (NF1) |
| Exon number | Exon structure | Exon 37 of 58 |
| CDS position | Coordinate mapping | c.6791 |
| Codon & frame | CDS math | Codon 2264, position 2 |
| Amino acid change | Translation | p.Arg2264Cys |
| Splice distance | Exon boundary | 47 bp from nearest boundary |
| Strand orientation | Ensembl | Reverse strand (auto-complemented) |

**Key point callout (teal box):**
> Reference allele is validated against the Ensembl genome sequence before any downstream analysis.
> In our 30-case benchmark, all variants passed reference validation.

---

## SLIDE 6: Feasibility Engines

**Background:** Dark navy
**Layout:** Three columns, one per modality, with modality-specific color coding

### Column 1: Base Editing (color: green #2EC4B6)

- **ABE window:** positions 4--7 (A>G corrections)
- **CBE window:** positions 4--8 (C>T corrections)
- **PAM-verified:** guide must place target within window
- **Bystander counting:** all C or A residues in window scanned
- **Bystander consequence classification:** synonymous / missense / nonsense
- Nucleases: SpCas9 (NGG), enFnCas9 (NRG)

### Column 2: Prime Editing (color: teal #00B4D8)

- **pegRNA design:** PBS (13 nt default, range 10--17) + RT template (10--30 nt)
- **PE3 nick search:** opposite strand, 40--100 bp from pegRNA nick
- **Scope:** substitutions, insertions up to 40 bp, deletions up to 80 bp
- **No editing window constraint** -- RT template directly encodes the edit
- DSB-free

### Column 3: HDR (color: gold #FFD166)

- **Cut-to-edit distance:** exponential decay model (half-life ~20 bp)
- **Donor type recommendation:**
  - ssODN: distance <=30 bp (90 bp arms)
  - cssDNA: distance <=5,000 bp (300 bp arms)
  - lssDNA/dsDNA: larger spans (800 bp arms)
- **Gene conversion probability** estimated from tract-length distributions
- Requires DSB (safety penalty in iPSC context)

**Bottom banner:**
> enFnCas9 (NRG PAM) -- developed in the Chakraborty lab -- evaluated as a first-class nuclease option alongside SpCas9

---

## SLIDE 7: Scoring Function

**Background:** Dark charcoal
**Layout:** Equation at top, weight table below, consequence penalties on right

**Equation (large, centered, white on dark):**

```
Score = w1*Safety + w2*Feasibility - w3*Complexity - w4*Risk + w5*Confidence
```

**Weight table (iPSC-optimized defaults):**

| Component | Weight | What it captures |
|-----------|--------|-----------------|
| Safety | 0.30 | DSB-free = 1.0; single DSB = 0.5; dual DSB = 0.2 |
| Feasibility | 0.25 | Modality-specific confidence in editing success |
| Complexity | 0.20 | Penalizes multi-step/multi-component designs |
| Risk | 0.15 | Off-target potential, unintended consequence likelihood |
| Confidence | 0.10 | Literature support, design completeness |

**Consequence penalties (coral accent box):**

| Consequence | Penalty |
|-------------|---------|
| Bystander missense | -0.10 per position |
| Bystander nonsense | -0.25 per position |
| Splice donor/acceptor proximity (<=2 bp) | -0.15 |
| Splice region proximity (3--8 bp) | -0.08 |
| Dual DSB in p53-active cells | -0.10 |
| All bystanders synonymous | +0.05 (bonus) |

**Key insight callout:**
> In iPSC context, a typical BE strategy scores 0.650 vs. 0.388 for an equivalent HDR strategy at the same locus -- driven primarily by the safety component (1.00 vs. 0.50).

---

## SLIDE 8: Benchmark Design

**Background:** Dark navy
**Layout:** Left: benchmark overview stats. Right: category breakdown table.

**Left panel -- Key numbers (large, teal):**
- **30** curated ClinVar cases
- **GRCh38** coordinates, all verified against Ensembl
- **11** variant categories
- **Tiered truth labels:** preferred / acceptable / reject
- **All reference alleles validated** against genome sequence

**Right panel -- Category breakdown:**

| Category | n | Examples |
|----------|---|---------|
| Clean BE substitutions | 7 | Transitions at PAM-accessible loci |
| BE negative control | 1 | Transversion falsely resembling BE target |
| PE transversions | 2 | HBB sickle cell (c.20A>T) |
| PE small indels | 5 | Insertions/deletions within PE scope |
| HDR large deletions | 3 | Multi-exon structural deletions |
| HDR/PE small deletions | 1 | Deletion feasible by both modalities |
| Compound het (hybrid) | 5 | Two variants requiring different modalities |
| Sequential HDR | 1 | Two distant variants requiring staged HDR |
| Dual base editing | 3 | Two variants both BE-correctable |
| Edge cases | 2 | Distant variants, non-coding positions |

**Footer:** "Truth labels defined by expert judgment based on published editing parameters and clinical precedent"

---

## SLIDE 9: Results -- Overall Performance

**Background:** Dark charcoal
**Layout:** Three large metric cards across the top; bar chart below (reference Figure 7)

**Metric cards (large numbers, teal background):**

| Metric | Value | Detail |
|--------|-------|--------|
| **Top-1 Accuracy** | **86.7%** | 26/30 cases: top-ranked strategy matches preferred or acceptable label |
| **Top-3 Accuracy** | **96.7%** | 29/30 cases: at least one appropriate strategy in top 3 |
| **Rejection Accuracy** | **90.0%** | 27/30 cases: infeasible strategies correctly excluded |

**Below cards:** Grouped bar chart from Figure 7

**Additional detail (small text, bottom):**
- 30/30 cases ran successfully -- zero pipeline failures
- 4 of 120 API calls initially failed (transient Ensembl errors); all recovered via automatic retry with exponential backoff
- PE top-ranked in 29/30 cases; HDR top-ranked in 1/30; BE top-ranked in 0/30

---

## SLIDE 10: Results -- Feasibility Heatmap

**Background:** Dark navy
**Layout:** Full-width heatmap (reference Figure 3)

**Heatmap description:**
- Rows: 30 benchmark cases (labeled by gene/variant)
- Columns: Base Editing, Prime Editing, HDR
- Color coding:
  - Green (+): modality produced the top-ranked strategy
  - Amber (?): modality was feasible but not top-ranked
  - White/dark: not feasible for this variant

**Three key takeaways (right side, bullet points):**

1. **Not all variants are editable by all methods** -- the heatmap shows substantial gaps, particularly for BE
2. **PE is most broadly applicable** -- feasible for nearly all 30 cases
3. **BE restricted by PAM-window constraints** -- even nominally compatible transitions frequently lack a guide placing the target at positions 4--7

**Callout (coral border):**
> Feasibility is locus-specific and PAM-dependent. Mutation-type alone is insufficient to determine editability.

---

## SLIDE 11: Results -- Per-Category Accuracy

**Background:** Dark charcoal
**Layout:** Grouped bar chart (reference Figure 5) with category labels

**Chart:** Two bars per category (blue = Top-1, green = Top-3)

| Category | Top-1 | Top-3 |
|----------|-------|-------|
| Clean BE substitutions (n=7) | 7/7 (100%) | 7/7 (100%) |
| BE negative control (n=1) | 1/1 (100%) | 1/1 (100%) |
| PE transversions (n=2) | 2/2 (100%) | 2/2 (100%) |
| PE small indels (n=5) | 5/5 (100%) | 5/5 (100%) |
| HDR large deletions (n=3) | 0/3 (0%) | 3/3 (100%) |
| HDR/PE small deletions (n=1) | 1/1 (100%) | 1/1 (100%) |
| Compound het hybrid (n=5) | 5/5 (100%) | 5/5 (100%) |
| Sequential HDR (n=1) | 0/1 (0%) | 0/1 (0%) |
| Dual base editing (n=3) | 3/3 (100%) | 3/3 (100%) |
| Edge cases (n=2) | 2/2 (100%) | 2/2 (100%) |

**Key observation (teal box):**
> 100% top-1 accuracy for all categories except HDR large deletions (0/3 top-1, but 3/3 top-3) and sequential HDR (0/1). These represent structural editing scenarios not yet fully modeled.

---

## SLIDE 12: Key Finding -- PAM Window Bottleneck

**Background:** Dark navy with subtle DNA helix watermark
**Layout:** Central finding statement with supporting evidence below

**Central statement (large, gold text on dark):**
> PAM-dependent editing window constraints are MORE restrictive than mutation-type classification alone.

**Evidence panel (three numbered points, white text):**

1. **7 ClinVar loci** with ABE-compatible transitions were tested
2. At **all 7 loci**, no SpCas9 (NGG) guide placed the target nucleotide at ABE positions 4--7
3. This is a **genuine biological constraint**, not a software limitation

**Implication box (teal border):**
> Estimates of base editing applicability based solely on the fraction of ClinVar pathogenic variants that are transitions may substantially **overestimate** the fraction that is practically editable with current base editors and canonical SpCas9.

**Bottom note:**
> Expanded-PAM nucleases (e.g., enFnCas9 with NRG PAM) partially alleviate this constraint but do not eliminate it.
> Locus-specific PAM scanning is essential for realistic BE feasibility assessment.

**Speaker note:** This is the most novel finding of the paper. Emphasize that this emerged from systematic evaluation, not from an a priori hypothesis. The community often quotes "~60% of ClinVar pathogenic SNVs are transitions" as evidence that BE has broad applicability -- our data suggest the practical fraction is substantially lower.

---

## SLIDE 13: Technical Innovation

**Background:** Dark charcoal
**Layout:** Two-column comparison (v1 vs. v2) with shared innovation highlights

### Left column: v1 Foundation

- **~24,000 lines of code**
- 5 validated modules:
  - ConversionSim (gene conversion tract modeling)
  - MOSAIC (multi-objective strategy analysis)
  - TopoPred (topological domain prediction)
  - ChromBridge (chromatin interaction modeling)
  - LoopSim (loop extrusion simulation)

### Right column: v2 Extension

- **~9,400 lines of code** (24 files)
- 4 new capability layers:
  - Transcript-aware mapping (Ensembl integration)
  - Modality-specific feasibility engines (BE, PE, HDR)
  - Consequence-aware scoring framework
  - Systematic benchmark with evaluation pipeline

### Bottom banner (shared highlights):

| Metric | Value |
|--------|-------|
| Total codebase | ~33,400 LOC |
| Test suite | **121 tests passing** |
| Nuclease support | SpCas9 (NGG) + enFnCas9 (NRG) as first-class options |
| Dependencies | Python 3.9+, NumPy, SciPy, Matplotlib only |
| Validation | All parameters validated against published experimental data |

---

## SLIDE 14: Limitations (Stated Honestly)

**Background:** Dark navy
**Layout:** Six limitation cards with severity indicators (amber/coral dots)

| # | Limitation | Severity | Detail |
|---|-----------|----------|--------|
| 1 | No experimental validation | Amber | Computational predictions only; no prospective design-outcome loop |
| 2 | Simplified CDS model | Amber | Treats spliced exonic transcript as surrogate for full coding sequence; adequate for exonic variants but not non-canonical architectures |
| 3 | No chromatin accessibility context | Amber | Does not incorporate ATAC-seq, replication timing, or epigenomic data |
| 4 | No off-target prediction | Amber | Defers to CRISPOR or Cas-OFFinder; evaluates on-target only |
| 5 | Large deletion handling | Coral | 3/3 HDR large deletion cases missed at top-1; pipeline does not model multi-exon structural constraints |
| 6 | iPSC-specific weights | Amber | Current scoring optimized for iPSC; other cell types would need retuning |

**Callout (bottom, teal text):**
> "Honest reporting of limitations is not a weakness -- it defines the actionable scope of the tool and builds trust with users and reviewers."

**Speaker note:** Present this slide confidently. Reviewers and PIs respect honesty. Each limitation maps directly to a future direction.

---

## SLIDE 15: Future Directions

**Background:** Dark charcoal
**Layout:** Five future directions as a roadmap timeline (left to right), with icons

| Priority | Direction | Description |
|----------|-----------|-------------|
| Near-term | **Chromatin accessibility integration** | Incorporate ENCODE/ATAC-seq data for locus-specific efficiency adjustment |
| Near-term | **HGVS parser** | Accept clinical variant nomenclature (e.g., NM_000546.6:c.215C>G) directly |
| Medium-term | **ML-based scoring refinement** | Train on experimental editing outcomes to supplement or replace heuristic weights |
| Medium-term | **Off-target integration** | Interface with Cas-OFFinder for comprehensive guide safety assessment |
| Long-term | **Experimental feedback loop** | Prospective validation: design experiments with CRISPRArchitect, use outcomes to refine the model |

**Bottom note:**
> The framework is designed to be modular and extensible. New editing modalities, scoring data, and nuclease variants can be integrated without restructuring the pipeline.

---

## SLIDE 16: Acknowledgements

**Background:** Dark navy
**Layout:** Clean, centered text with institutional logo placeholder

**Acknowledgements:**

- **CSIR-Institute of Genomics and Integrative Biology**, New Delhi
- **Debojyoti Chakraborty Lab** -- supervision, enFnCas9 development, iPSC editing expertise
- **enFnCas9 development team** -- for the expanded-PAM nuclease that CRISPRArchitect supports as a first-class option
- **Academy of Scientific and Innovative Research (AcSIR)** -- doctoral program
- **CSIR** -- fellowship support

**Logo row:** CSIR-IGIB logo | AcSIR logo | CSIR logo

---

## SLIDE 17: Thank You / Questions

**Background:** Dark navy with subtle double-helix motif (matching Slide 1)
**Layout:** Centered, minimal

**Large text (white):**
> Thank you

**Contact information (teal accent):**

| | |
|---|---|
| GitHub | `github.com/visvikbharti/CRISPRArchitect` |
| Code | Open source, Python 3.9+ |
| Manuscript | In preparation for PLOS Computational Biology |

**QR code placeholder:** Link to GitHub repository

**Bottom text (small, gray):**
> All benchmark data, pipeline code, scoring parameters, and evaluation scripts are publicly available for independent reproduction.

---

## APPENDIX: Design Notes for PowerPoint Conversion

### Color Palette

| Use | Color | Hex |
|-----|-------|-----|
| Background (primary) | Dark navy | #0D1B2A |
| Background (alternate) | Dark charcoal | #1B2838 |
| Text (primary) | White | #FFFFFF |
| Text (secondary) | Light gray | #A8B2C1 |
| Accent 1 (highlights, links) | Teal | #00B4D8 |
| Accent 2 (warnings, negatives) | Coral | #FF6B6B |
| Accent 3 (key numbers, findings) | Gold | #FFD166 |
| Base editing | Green-teal | #2EC4B6 |
| Prime editing | Teal | #00B4D8 |
| HDR | Gold | #FFD166 |

### Typography

| Element | Font | Size | Weight |
|---------|------|------|--------|
| Slide title | Montserrat | 36--40 pt | Bold |
| Body text | Inter or Open Sans | 18--22 pt | Regular |
| Code/numbers | JetBrains Mono | 16--20 pt | Regular |
| Footnotes | Inter | 12--14 pt | Light |
| Large metric numbers | Montserrat | 72--96 pt | Bold |

### General Guidelines

- Maximum 6 bullet points per slide
- No full sentences on slides -- use fragments and keywords
- All numbers must match verified values (see header of this document)
- Figures referenced (Fig. 1, 3, 5, 7) correspond to manuscript figures
- Animations: use simple fade-in for sequential reveals; no spinning or bouncing effects
- Aspect ratio: 16:9 widescreen
- Export: .pptx for editing, .pdf for distribution

### Timing Guide

| Slide | Duration | Cumulative |
|-------|----------|------------|
| 1. Title | 30 sec | 0:30 |
| 2. The Problem | 2 min | 2:30 |
| 3. What CRISPRArchitect Does | 2 min | 4:30 |
| 4. Pipeline Architecture | 2 min | 6:30 |
| 5. Transcript-Aware Mapping | 2 min | 8:30 |
| 6. Feasibility Engines | 3 min | 11:30 |
| 7. Scoring Function | 2 min | 13:30 |
| 8. Benchmark Design | 1.5 min | 15:00 |
| 9. Results -- Overall | 2 min | 17:00 |
| 10. Results -- Heatmap | 1.5 min | 18:30 |
| 11. Results -- Per-Category | 1.5 min | 20:00 |
| 12. Key Finding -- PAM | 2 min | 22:00 |
| 13. Technical Innovation | 1 min | 23:00 |
| 14. Limitations | 1.5 min | 24:30 |
| 15. Future Directions | 1 min | 25:30 |
| 16. Acknowledgements | 30 sec | 26:00 |
| 17. Thank You | -- | -- |
| **Total** | **~26 min** | + Q&A |
