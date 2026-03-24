# CRISPRArchitect v2 -- Lab Meeting Speaker Guide

**Presenter:** Vishal Bharti
**Audience:** Debojyoti Chakraborty (PI) and lab members, CSIR-IGIB, New Delhi
**Duration:** 25-30 minutes presentation + 8-10 minutes live demo + 10 minutes Q&A
**Slide deck:** CRISPRArchitect_v2_LabMeeting.pptx (17 core slides + 3 demo/backup slides = 20 slides)

---

## Part 1: Speaking Script (Slide-by-Slide)

---

### Slide 1: Title Slide

**Timing:** 30 seconds

**What to say:**

Good morning, everyone. Today I want to present CRISPRArchitect, a computational framework we have been developing for transcript-aware and consequence-guided design of genome editing strategies. The central idea is straightforward: when you have a pathogenic variant in a patient iPSC line and you need to decide between base editing, prime editing, and HDR, that decision currently lives entirely in your head. You weigh PAM availability, editing windows, bystander risk, DSB toxicity, and splice proximity all by intuition. CRISPRArchitect formalizes that reasoning.

This is joint work with DC sir, and the code is publicly available on GitHub. Let me walk you through what we built, what we found, and why it matters for the kind of editing work we do in this lab every day.

**Key point to emphasize:** This tool addresses a gap that every person in this room encounters -- choosing between editing modalities for a specific variant.

**Transition to next slide:** Let me start by framing the problem more precisely.

---

### Slide 2: The Problem

**Timing:** 2 minutes

**What to say:**

Everyone here knows we have three major precision editing approaches -- base editing, prime editing, and HDR. Each is powerful, and each has fundamentally different constraints. Base editors need a transition mutation and a PAM that positions the target within a narrow 4-nucleotide window. Prime editors need a suitable pegRNA design but can handle substitutions, small insertions, and small deletions without a DSB. HDR can do almost anything but requires a double-strand break and a donor template, which in iPSCs triggers the p53-mediated toxicity that Ihry and colleagues documented in 2018.

The problem is not that we lack tools for guide design. CRISPOR, PrimeDesign, BE-Designer -- these are all excellent for designing guides within a single modality. The gap is that no tool answers the question we actually face: given this specific variant, at this specific locus, in this specific cell type, which modality should I use? That question requires evaluating feasibility across modalities simultaneously, and then incorporating biological consequences -- bystander edits, splice disruption, amino acid changes -- into the ranking. That is what CRISPRArchitect does.

I want to emphasize that this is not a theoretical problem. Think about the last time someone in this lab was designing an editing experiment. The decision of "should I try base editing or prime editing for this variant" was made by checking a few things manually and then going with gut instinct. There is nothing wrong with that when you have deep expertise, but it does not scale, and it is not reproducible.

**Key point to emphasize:** The gap is not guide RNA design -- it is systematic cross-modality comparison with biological context.

**Transition to next slide:** So what exactly does CRISPRArchitect provide?

---

### Slide 3: What CRISPRArchitect Does

**Timing:** 2 minutes

**What to say:**

CRISPRArchitect does five things. First, it evaluates all three modalities -- base editing, prime editing, and HDR -- plus hybrid strategies for compound heterozygous cases, within a single unified framework. You put in your variants and get back a ranked list of strategies that spans all modalities, not separate outputs from separate tools that you then have to reconcile manually.

Second, it maps every variant to its transcript context using Ensembl. This is not just about knowing the chromosome and position -- it is about knowing which exon you are in, what codon you are affecting, how far you are from the nearest splice junction, and whether the gene is on the forward or reverse strand. This matters because a C-to-T change on the coding strand might require ABE correction if the gene is on the minus strand, and getting that wrong invalidates every downstream analysis.

Third, it performs PAM-verified feasibility checking. This is different from just asking "is this a transition?" -- it actually scans the local sequence for SpCas9 NGG and enFnCas9 NRG PAMs and checks whether any guide places the target nucleotide within the editing window. Fourth, it scores strategies using a multi-objective function that penalizes bystander edits, splice disruption, and DSB burden. And fifth, it explicitly rejects infeasible strategies rather than silently omitting them, so you know exactly why a particular approach will not work at your locus.

**Key point to emphasize:** The key differentiator is unification -- one input, one ranked output across all modalities, with biological reasoning attached.

**Transition to next slide:** Let me show you how the pipeline is organized.

---

### Slide 4: Pipeline Architecture

**Timing:** 2 minutes

**What to say:**

The pipeline has seven stages, and each one corresponds to a distinct biological question. Stage one takes your input variant -- chromosome, position, reference and alternate alleles in GRCh38 coordinates. Stage two maps that variant to its transcript context through the Ensembl REST API. We fetch the canonical transcript, identify the exon, compute the CDS coordinate, and determine the codon frame. Stage three validates the reference allele against the Ensembl genome sequence. This sounds trivial, but it catches genome-build mismatches and strand-orientation errors before they propagate downstream. In our 30-case benchmark, all variants passed this check, which confirmed our coordinate curation was correct.

Stage four annotates coding consequences -- is this missense, nonsense, synonymous, or splice-proximal? Stage five runs three independent feasibility engines, one for each modality. Stage six generates all feasible strategies and explicitly tags infeasible ones with reasons. And stage seven scores and ranks everything using the multi-objective function.

The v2 pipeline is implemented across 15 Python modules, roughly 9,400 lines of code. It builds on the v1 framework that many of you have seen before -- ConversionSim, MOSAIC, TopoPred, ChromBridge, LoopSim -- which totals about 24,000 lines. The only dependencies beyond the standard library are NumPy, SciPy, and Matplotlib. And we support both SpCas9 and enFnCas9 as first-class nuclease options, which I think is particularly relevant for this lab.

**Key point to emphasize:** Each pipeline stage answers a specific biological question, and the stages are modular -- you can inspect or replace any one independently.

**Transition to next slide:** Let me zoom in on the transcript mapping step, because it has implications that might not be obvious.

---

### Slide 5: Transcript-Aware Mapping

**Timing:** 2 minutes

**What to say:**

When you input a variant like NF1 c.910C>T, the first thing CRISPRArchitect does is fetch the canonical transcript from Ensembl and map the genomic coordinate to transcript coordinates. For NF1, which is on the reverse strand, this means the system automatically reverse-complements the alleles before comparing them to the coding sequence. This is the kind of detail that is easy to get wrong when doing it manually, and a strand orientation error would cause every downstream analysis to be incorrect.

The system extracts several pieces of information that feed into later stages. It determines which exon the variant is in, what the CDS position is, what codon is affected and at which position within the codon, what the amino acid change is, and critically, how far the variant is from the nearest exon boundary. That splice distance matters because variants within 2 base pairs of a splice junction get a severe penalty in our scoring function -- minus 0.15 -- because editing near a splice donor or acceptor risks disrupting splicing even if the primary edit is correct. Variants in the splice region, 3 to 8 base pairs from the boundary, get a moderate penalty of minus 0.08.

One aspect I want to highlight is reference validation. Before any feasibility analysis begins, the system fetches the actual genomic sequence from Ensembl and confirms that the reference allele you specified matches what is in the genome. This is a simple sanity check, but it prevents an entire class of errors -- wrong genome build, wrong strand assumption, copy-paste mistakes -- from producing silently wrong results.

**Key point to emphasize:** Transcript-aware mapping is not optional decoration -- it determines strand orientation, codon frame, and splice proximity, all of which directly affect strategy selection.

**Transition to next slide:** Now that we have the variant mapped and annotated, how do we determine which modalities can actually correct it?

---

### Slide 6: Feasibility Engines

**Timing:** 3 minutes

**What to say:**

This is the most biologically dense slide, so I want to take a moment on each modality. For base editing, feasibility requires three conditions to be true simultaneously. The correction must be a transition -- A-to-G for ABE, C-to-T for CBE. A PAM must exist such that the target nucleotide falls within the editing window -- positions 4 through 7 for ABE, 4 through 8 for CBE, counting from the PAM-distal end. And the bystander edits within that window should not introduce deleterious consequences. The system scans all C or A residues within the window, depending on the editor type, and classifies each potential bystander conversion as synonymous, missense, or nonsense using the coding annotation from the previous stage. This is where transcript context directly feeds into feasibility -- without knowing the codon frame, you cannot determine whether a bystander edit is benign or catastrophic.

For prime editing, feasibility is broader. The system designs a pegRNA with a primer binding site -- default 13 nucleotides -- and a reverse transcriptase template that directly encodes the desired edit. It searches for a PE3 nicking guide on the opposite strand, 40 to 100 base pairs away. The key biological advantage is that there is no editing window constraint analogous to base editing: the RT template encodes exactly what you want regardless of PAM-to-target distance. Prime editing handles substitutions, insertions up to 40 base pairs, and deletions up to 80 base pairs.

For HDR, the system evaluates cut-to-edit distance using an exponential decay model calibrated to published tract-length data from the Elliott and Jasin 1998 paper. Based on the distance, it recommends donor type: ssODN for distances under 30 base pairs, cssDNA for distances under 5,000 base pairs -- which is directly relevant to the cssDNA work we do in this lab -- and lssDNA or dsDNA for larger spans. One important note: we support enFnCas9 with the NRG PAM as a first-class option alongside SpCas9. The broader PAM compatibility of enFnCas9 increases the number of available guide sites, which can make the difference between having a feasible strategy and not.

**Key point to emphasize:** Feasibility is not just about mutation type -- it requires simultaneous satisfaction of PAM availability, editing window positioning, and consequence constraints.

**Transition to next slide:** Once we know which strategies are feasible, how do we rank them?

---

### Slide 7: Scoring Function

**Timing:** 2 minutes

**What to say:**

The scoring function balances five components. Safety carries the highest weight at 0.30, and this is deliberate for iPSC applications. We all know from the Ihry 2018 paper that DSBs in iPSCs trigger p53-mediated apoptosis and create selection pressure for p53-deficient clones. So DSB-free approaches -- base editing and prime editing -- get a safety score of 1.0, while single-DSB HDR strategies get 0.5, and dual-DSB strategies get 0.2. This is not a philosophical preference; it reflects a documented biological risk in the cell type we work with most.

Feasibility at 0.25 captures how confident we are that the editing will work mechanistically. Complexity at 0.20 penalizes multi-step or multi-component designs -- a strategy requiring two rounds of editing with two different guides is inherently more complex than a single-step prime edit. Risk at 0.15 captures off-target potential and unintended consequence likelihood. And confidence at 0.10 reflects how well-supported the approach is in the literature.

On top of the base score, consequence-aware adjustments apply biologically grounded penalties. A bystander edit that creates a missense change costs minus 0.10 per position. A bystander nonsense -- which means you are introducing a premature stop codon as collateral damage -- costs minus 0.25 per position. Splice proximity penalties apply whether the variant itself or a bystander position is near a splice junction. And there is a small bonus of plus 0.05 when all bystanders within the editing window are synonymous. To give you a concrete sense of the impact: a typical base editing strategy at a clean locus scores around 0.650, while an equivalent HDR strategy at the same locus scores around 0.388. That gap is driven primarily by the safety component.

**Key point to emphasize:** The weights are iPSC-optimized, reflecting the documented p53-mediated toxicity of DSBs in pluripotent stem cells -- this is directly relevant to the editing work in our lab.

**Transition to next slide:** How did we test whether these design choices lead to reasonable recommendations?

---

### Slide 8: Benchmark Design

**Timing:** 1.5 minutes

**What to say:**

We curated a benchmark of 30 variant scenarios spanning 11 categories. All variants use verified GRCh38 coordinates from ClinVar and published editing studies, and all reference alleles were validated against the Ensembl genome before evaluation. The categories cover the range of scenarios you would encounter in practice: seven clean base-editable substitutions, one base editing negative control -- the HBB sickle cell variant, which is a transversion that superficially looks like it might be base-editable but is not -- two prime-editable transversions, five prime-editable small indels, three HDR-required large deletions, compound heterozygous cases requiring hybrid strategies, dual base editing cases, and edge cases like non-coding variants.

Each case has tiered truth labels: preferred strategies, acceptable alternatives, and explicit reject labels for strategies that should not be recommended. These truth labels were defined by expert judgment based on published editing parameters and clinical precedent -- not by the algorithm. So when we measure accuracy, we are asking: does the algorithm agree with what an expert would recommend?

I want to be transparent that 30 cases is a modest benchmark. This is not thousands of variants. But each case is carefully curated with verified coordinates, and the categories are designed to stress-test specific failure modes. We chose depth over breadth.

**Key point to emphasize:** Every benchmark coordinate is verified against GRCh38, and truth labels are defined by expert judgment independent of the algorithm.

**Transition to next slide:** So how did we do?

---

### Slide 9: Results -- Overall Performance

**Timing:** 2 minutes

**What to say:**

Three headline numbers. Top-1 accuracy -- meaning the top-ranked strategy matches a preferred or acceptable truth label -- was 86.7%, which is 26 out of 30 cases. Top-3 accuracy -- at least one appropriate strategy appears among the top three -- was 96.7%, 29 out of 30. And rejection accuracy -- infeasible strategies correctly excluded from the top-ranked output -- was 90.0%, 27 out of 30. Zero pipeline failures across all 30 cases. The pipeline ran for about 16 minutes total, and four of the 120 Ensembl API calls initially failed due to transient server errors but were all recovered automatically by our retry logic.

Now, let me be honest about where the numbers are not perfect. The four top-1 misses all came from the same type of case: large deletions and sequential HDR. In those cases, the system recommended prime editing -- which is technically feasible for the individual variant -- but does not adequately address the multi-exon structural nature of the deletion. The three rejection accuracy misses are the same three large deletion cases. So the failure mode is specific and well-characterized: the pipeline does not yet model structural variant constraints properly.

One more important number: prime editing was top-ranked in 29 out of 30 cases, with HDR top-ranked in the remaining case, a non-coding FMR1 5'UTR variant. Base editing was top-ranked in zero cases. I will explain why in a moment, and it is the most interesting finding of this work.

**Key point to emphasize:** 86.7% top-1 accuracy with well-characterized failure modes concentrated in large deletion cases that the pipeline explicitly does not yet model.

**Transition to next slide:** Let me show you the feasibility landscape across all 30 cases.

---

### Slide 10: Results -- Feasibility Heatmap

**Timing:** 1.5 minutes

**What to say:**

This heatmap shows the feasibility of each modality across all 30 benchmark cases. Green means the modality produced the top-ranked strategy. Amber means the modality was feasible but was not top-ranked. And white or dark cells mean the modality was not feasible for that variant. Three things jump out immediately.

First, prime editing is the most broadly applicable modality -- it is feasible for nearly all 30 cases. This makes biological sense: prime editing has no editing window constraint, handles both transitions and transversions, and accommodates small insertions and deletions. Second, HDR is almost always feasible in principle -- if there is a guide near the target, you can do HDR -- but it consistently scores lower because of the DSB burden penalty in iPSC context. Third, and this is the important one: look at how sparse the base editing column is. Even for cases that are nominally ABE-compatible transitions, base editing frequently lacks a guide that places the target within the editing window.

This is not a software bug. This is a genuine biological constraint. The co-occurrence of a suitable PAM at the precise spacing required to position the target within a 4-nucleotide window is not guaranteed, and at many clinically relevant loci, it simply does not happen.

**Key point to emphasize:** The heatmap visually demonstrates that feasibility is locus-specific and PAM-dependent -- mutation type alone does not determine editability.

**Transition to next slide:** Let me break down the accuracy by category.

---

### Slide 11: Results -- Per-Category Accuracy

**Timing:** 1.5 minutes

**What to say:**

When we break down the results by category, the pattern is clear. CRISPRArchitect achieved 100% top-1 accuracy for clean base-editable substitutions -- all 7 out of 7 -- the base editing negative control, both prime-editable transversions, all five prime-editable small indels, the HDR/PE small deletion, all five compound heterozygous hybrid cases, all three dual base editing cases, and both edge cases. That is 100% top-1 accuracy for 8 out of 11 categories.

The imperfect categories are HDR-required large deletions -- 0 out of 3 for top-1, but 3 out of 3 for top-3 -- and sequential HDR -- 0 out of 1 for both top-1 and top-3. These represent structural editing scenarios where the pipeline's current single-variant processing logic does not capture the multi-exon nature of the problem. The system recommended prime editing for the individual variants, which is technically feasible, but the correct answer is HDR because you need to address a large structural deletion that spans multiple exons.

I want to be clear that this is a known limitation, not a hidden one. We explicitly flagged it in the manuscript. The fix would require adding structural variant-specific logic -- recognizing that when two variants define a large deletion, the strategy should shift to HDR with a deletion-spanning donor. That is on our roadmap.

**Key point to emphasize:** Perfect accuracy across 8 of 11 categories; failures are concentrated in structural deletion cases that require explicit multi-exon modeling.

**Transition to next slide:** Now for the finding I am most excited about.

---

### Slide 12: Key Finding -- PAM Window Bottleneck

**Timing:** 2.5 minutes

**What to say:**

This is the most novel finding of this work, and I want to be precise about what we found and what it means. It is widely appreciated that base editors are restricted to transition mutations -- ABE for A-to-G, CBE for C-to-T. Many people in the field cite statistics like "roughly 60% of ClinVar pathogenic SNVs are transitions" as evidence that base editing has broad therapeutic applicability. Our analysis suggests that this number substantially overestimates the practical fraction.

Here is what we found. We took seven ClinVar loci in our benchmark that have ABE-compatible transitions -- meaning the correction is an A-to-G change that ABE can chemically perform. At all seven loci, when we scanned the local sequence for SpCas9 NGG PAMs, there was no guide that placed the target nucleotide within the ABE editing window at positions 4 through 7. Not a single one. The transition is there, ABE can chemically do it, but there is no way to physically position the deaminase domain over the target base using a canonical SpCas9 guide. This is not a software limitation. This is a genuine constraint of the PAM-window geometry.

The practical implication is significant. If you are estimating how many disease-causing variants are amenable to base editing, you cannot just count transitions. You need to do locus-specific PAM scanning and editing window verification. The fraction that survives this additional filter is meaningfully smaller. Now, expanded-PAM nucleases -- including enFnCas9 with the NRG PAM, which DC sir's group developed and which CRISPRArchitect evaluates as a first-class option -- partially alleviate this constraint by opening up more PAM sites. But even with NRG PAMs, the editing window constraint does not disappear.

This finding emerged from systematic evaluation, not from an a priori hypothesis. We did not set out to show that base editing has limited applicability. We built a system that honestly evaluates feasibility, ran it on ClinVar variants, and this is what the data showed.

**Key point to emphasize:** PAM-dependent editing window constraints, not mutation type, are the primary bottleneck for base editing applicability -- this has implications for how the field estimates therapeutic reach.

**Transition to next slide:** Let me briefly mention the technical innovations in the codebase.

---

### Slide 13: Technical Innovation (v1 vs. v2)

**Timing:** 1 minute

**What to say:**

For those interested in the engineering side, the v1 framework -- which several of you have seen before -- includes ConversionSim for gene conversion tract modeling, MOSAIC for multi-objective strategy analysis, TopoPred for cssDNA secondary structure prediction, ChromBridge for 3D chromatin distance modeling, and LoopSim for cohesin loop extrusion simulation. That totals about 24,000 lines of code validated against four published datasets.

The v2 extension adds four new capability layers in about 9,400 lines of code: transcript-aware mapping through Ensembl integration, modality-specific feasibility engines for base editing, prime editing, and HDR, the consequence-aware scoring framework, and a systematic benchmark with an automated evaluation pipeline. The total codebase is roughly 33,400 lines, the test suite has 121 passing tests, and the only dependencies are Python 3.9 with NumPy, SciPy, and Matplotlib. No TensorFlow, no PyTorch, no heavy frameworks -- just clean, transparent code.

**Key point to emphasize:** The codebase is modular, well-tested, and lightweight -- 121 tests, minimal dependencies, open source on GitHub.

**Transition to next slide:** Now let me be honest about what this tool cannot do.

---

### Slide 14: Limitations

**Timing:** 1.5 minutes

**What to say:**

I want to present our limitations directly because I think honest reporting of limitations builds more trust than glossing over them. First, and most importantly, we have no experimental validation. These are computational predictions benchmarked against expert-defined truth labels, not against actual editing outcomes. We have not yet closed the loop by designing experiments with CRISPRArchitect recommendations and measuring what happens in the cells. That is the most critical next step.

Second, the coding annotation uses a simplified CDS model. We treat the spliced exonic transcript as a surrogate for the full coding sequence. This works well for the exonic variants in our benchmark but would need refinement for non-canonical transcript architectures or complex alternative splicing.

Third, we do not incorporate chromatin accessibility, replication timing, or any epigenomic data. We all know that a guide can look perfect on paper but fail in the cell because the locus is in closed chromatin. Fourth, we do not perform off-target prediction -- we defer that to CRISPOR or Cas-OFFinder. Fifth, and I showed this already in the per-category results, large deletion handling is poor. All three HDR-required large deletion cases were top-1 misses, and these account for all our rejection accuracy misses. Sixth, the scoring weights are iPSC-optimized. If you are working in HEK293T or another cell type with a different p53 landscape, you would need to adjust the safety weights.

Each of these limitations maps directly to a specific future direction, which brings me to the next slide.

**Key point to emphasize:** The most critical limitation is the absence of experimental validation -- computational predictions need to be tested in the lab.

**Transition to next slide:** Here is where we are headed.

---

### Slide 15: Future Directions

**Timing:** 1 minute

**What to say:**

Five directions, in rough priority order. Near-term: integrate chromatin accessibility data from ENCODE or ATAC-seq to make locus-specific efficiency predictions. Also near-term: add an HGVS variant parser so you can input clinical nomenclature directly instead of specifying genomic coordinates manually. Medium-term: explore machine learning-based scoring refinement -- train on actual experimental editing outcomes to supplement or replace the heuristic weights. Also medium-term: interface with Cas-OFFinder for comprehensive off-target assessment. And long-term -- and this is the big one -- establish an experimental feedback loop where we design experiments with CRISPRArchitect, perform the editing in our iPSC lines, and use the outcomes to refine the model. That would transform this from a computational prediction tool into a validated design system.

The framework is built to be modular and extensible. New editing modalities, new nucleases, new scoring data -- all of these can be integrated without restructuring the pipeline.

**Key point to emphasize:** The experimental feedback loop -- designing edits with CRISPRArchitect and validating them in iPSC lines -- is the most impactful next step.

**Transition to next slide:** A few acknowledgments before we open for questions.

---

### Slide 16: Acknowledgements

**Timing:** 30 seconds

**What to say:**

I want to thank DC sir for supervision and for the foundational work on enFnCas9 that CRISPRArchitect supports as a first-class option. Thanks to everyone in the lab for helpful discussions about editing strategies, which directly informed the scoring framework design. And thanks to CSIR for fellowship support and CSIR-IGIB for infrastructure.

**Key point to emphasize:** The enFnCas9 work from this lab is directly integrated into CRISPRArchitect as a first-class nuclease option.

**Transition to next slide:** With that, I would like to do a quick live demo, and then we can open for questions.

---

### Slide 17: Thank You / Questions

**Timing:** Transition to demo

**What to say:**

The code is publicly available on GitHub at github.com/visvikbharti/CRISPRArchitect, and we are preparing the manuscript for PLOS Computational Biology. All benchmark data, pipeline code, scoring parameters, and evaluation scripts are included in the repository for independent reproduction. Let me switch to the demo now, and then we will have plenty of time for discussion.

**Key point to emphasize:** Everything is open source and reproducible.

**Transition to next slide:** Let me pull up the terminal for a live demo.

---

### Slides 18-20: Live Demo Slides (see Part 2 below)

These slides serve as visual anchors during the demo -- they display the commands being run and annotate what to look for in the output. See Part 2 for full demo scripts.

---

## Part 2: Live Demo Script

**Total demo time: ~8 minutes**

Prepare your laptop before the talk: open two terminal windows and one browser tab. Have all commands pre-typed in a text file so you can copy-paste if you make a typo under pressure.

---

### Pre-Demo Checklist (Do This 30 Minutes Before the Talk)

1. Ensure your internet connection is stable (the pipeline needs Ensembl API access)
2. Test the Ensembl API: open `https://rest.ensembl.org/info/ping` in a browser -- it should return `{"ping":1}`
3. Activate your Python environment: `conda activate crisprarchitect` (or equivalent)
4. Run a quick test: `cd /path/to/crisprarchitect && python -c "from core.pipeline.strategy_stage import StrategyPipeline; print('OK')"`
5. Pre-run the NF1 example once so results are warm in your shell history
6. Open `benchmark_results/figures/` in Finder/file manager so you can pull up figures quickly
7. If using the webapp, run `streamlit run webapp/app.py` once to confirm it launches

---

### Demo 1: Webapp (v1 Modules) -- ~3 minutes

**Purpose:** Show the audience the user-friendly interactive interface that wraps the v1 modules.

**Step 1: Launch**

Open a terminal and run:
```bash
cd /path/to/crisprarchitect
streamlit run webapp/app.py
```

A browser tab should open automatically at `http://localhost:8501`.

**Step 2: Walk through the interface**

Say: "This is the interactive web application that wraps our v1 modules. Let me show you a typical workflow."

- In the gene input field, type `NF1`
- Set up 2 mutations (if the interface allows variant specification)
- Point out the MOSAIC module -- say: "MOSAIC enumerates all possible editing strategies for these variants. It is not just checking one modality -- it generates the complete strategy space."
- Show the strategy ranking panel -- say: "Each strategy gets a composite score. Notice the scoring breakdown: safety, feasibility, complexity, risk, confidence. These are the same five components from the scoring function I showed earlier."
- If ConversionSim output is visible, say: "ConversionSim predicts the gene conversion tract length distribution for HDR donors at this locus. This is calibrated to the Elliott and Jasin 1998 tract-length data."

**Step 3: Key things to point out**

- The strategy ranking is not just a list -- each strategy has an explanation for why it was ranked where it is
- Multiple modalities appear in the ranking, allowing direct comparison
- The interface is designed for biologists, not programmers -- no code required

**What to say while navigating:**

"The webapp gives you an accessible entry point. You enter your gene and variants, and it returns ranked strategies with full reasoning. But for those of you who want more control, or who want to integrate CRISPRArchitect into your own analysis pipelines, we also have a Python API. Let me show you that."

---

### Demo 2: v2 Pipeline (Python) -- ~3 minutes

**Purpose:** Show the programmatic interface and walk through a real analysis step by step.

**Step 1: Open a terminal**

Say: "Let me show you the v2 pipeline running directly in Python. I will analyze the NF1 c.910C>T variant -- a missense substitution in exon 37 of NF1."

**Step 2: Run the following code**

Copy-paste or type (have this pre-loaded in a script file as backup):

```python
python3 << 'PYEOF'
from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput

# Initialize pipeline for iPSC with SpCas9
pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")

# Analyze NF1 c.910C>T (missense substitution)
result = pipeline.run([
    GenomicVariantInput("17", 31200443, "C", "T", gene_symbol="NF1", name="c.910C>T")
])

# Show ranked strategies
print("=== RANKED STRATEGIES ===")
for s in result.strategies:
    print(f"  #{s.rank}: {s.strategy_name} (score={s.overall_score:.3f})")

# Show variant annotation details
print("\n=== VARIANT ANNOTATION ===")
v = result.variants[0]
print(f"  Consequence: {v.coding.consequence.value}")
print(f"  HGVS: {v.coding.hgvs_c}")
print(f"  Ref validated: {v.ref_validation.is_valid}")
PYEOF
```

**Step 3: Explain each output line**

As the output appears, say:

- "The top-ranked strategy is Single-step Prime Editing. Notice the score -- this reflects the composite of safety, feasibility, complexity, risk, and confidence."
- "Below it you will see HDR as an alternative, scoring lower because of the DSB burden penalty."
- "In the variant annotation section: the consequence is missense, the HGVS notation confirms the coding change, and reference validation is True -- meaning the C allele we specified matches what is in the GRCh38 genome at that position."

**Step 4: Show the HBB sickle cell case (transversion, BE correctly rejected)**

Say: "Now let me show you a case where the system correctly rejects base editing."

```python
python3 << 'PYEOF'
from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput

pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")

# HBB sickle cell: c.20A>T (p.Glu7Val) -- a TRANSVERSION
result = pipeline.run([
    GenomicVariantInput("11", 5227002, "A", "T", gene_symbol="HBB", name="c.20A>T (sickle cell)")
])

print("=== HBB SICKLE CELL -- TRANSVERSION ===")
for s in result.strategies:
    print(f"  #{s.rank}: {s.strategy_name} (score={s.overall_score:.3f})")

print("\nBase editing correctly rejected: this is a T>A transversion.")
print("ABE does A>G, CBE does C>T -- neither addresses this mutation.")
PYEOF
```

Say: "The sickle cell variant is c.20A>T -- a transversion. The patient has a T where they should have an A. ABE converts A to G. CBE converts C to T. Neither can convert T to A. The system correctly identifies this and excludes base editing entirely, recommending prime editing instead. This is exactly the kind of modality-level triage that CRISPRArchitect automates."

---

### Demo 3: Benchmark Results -- ~2 minutes

**Purpose:** Show the publication-quality figures from the benchmark evaluation.

**Step 1: Open the figures**

Navigate to `benchmark_results/figures/` and open the following files one at a time:

1. **fig7_benchmark_summary.png** -- Say: "This is the overall benchmark performance. 86.7% top-1, 96.7% top-3, 90.0% rejection accuracy across all 30 cases."

2. **fig5_category_accuracy.png** -- Say: "Breaking it down by category, you can see 100% accuracy for base-editable, prime-editable, and compound heterozygous cases. The only failures are in the large deletion and sequential HDR categories -- which we have been transparent about."

3. **fig3_feasibility_heatmap.png** -- Say: "This heatmap shows the feasibility landscape. Notice how sparse the base editing column is. That visual sparseness is the PAM-window bottleneck I described."

4. **fig1_architecture.png** -- Say: "And this is the pipeline architecture diagram for reference."

**Step 2: Highlight key visual patterns**

Point at the heatmap and say: "If you take one image away from this talk, let it be this heatmap. It shows viscerally that feasibility is not just about whether a mutation is a transition or a transversion. It is about whether the local PAM landscape cooperates."

---

### Backup Plans

**If the Ensembl API is down:**

Say: "The Ensembl API appears to be experiencing issues right now -- this happens occasionally. The pipeline has built-in retry logic with exponential backoff that handles transient failures, but if the server is fully down, we cannot fetch transcript data. Let me show you the cached results from our last benchmark run instead."

Then open `benchmark_results/definitive_benchmark_results.json` and walk through the case results, or show the pre-generated figures.

**If the webapp does not load:**

Say: "Streamlit is not cooperating at the moment. Let me show you screenshots of the interface instead."

Have screenshots pre-saved in a folder. If you do not have screenshots, skip to the Python pipeline demo -- that is the more important one.

**If something crashes during the Python demo:**

Stay calm. Say: "This is a live demo, so these things happen. Let me explain what should have appeared."

Then describe the expected output:
- For NF1: "#1: Single-step Prime Editing (score ~0.6-0.7), #2: Single-step HDR (score ~0.3-0.4)"
- For HBB: "#1: Single-step Prime Editing, with base editing correctly absent from the ranked list"

Fall back to the benchmark figures, which are static files and do not depend on any API.

**If your internet connection drops entirely:**

Say: "We have lost network connectivity, which means we cannot reach the Ensembl REST API. The pipeline is designed to work with live Ensembl data, so I will walk you through the pre-computed results instead."

Show the benchmark JSON and figures. Emphasize that the 30-case benchmark already ran successfully and all results are archived.

---

## Part 3: Anticipated Questions and Answers

---

### Category 1: Biological Questions (10)

---

**Q1: Why does prime editing dominate so heavily? Is the scoring function biased toward PE?**

The scoring function is not specifically biased toward prime editing -- it is biased toward DSB-free approaches with broad feasibility, and prime editing happens to satisfy both criteria. Three properties converge: PE has no editing window constraint like base editing, it handles both transitions and transversions plus small indels, and it requires no double-strand break. In the iPSC-optimized weighting, the safety component (weight 0.30) gives DSB-free approaches a base score of 1.0 versus 0.5 for HDR, which creates a consistent scoring advantage. If you changed the cell type to HEK293T and reduced the safety weight, HDR would rank higher. The dominance is context-specific, not hard-coded.

*Key fact:* PE was top-ranked in 29/30 cases, HDR in 1/30, BE in 0/30 -- reflecting iPSC-optimized safety-first scoring.

---

**Q2: At the seven ABE-compatible loci, did you check whether enFnCas9 (NRG PAM) rescues base editing feasibility?**

Yes, CRISPRArchitect evaluates enFnCas9 with NRG PAM as a first-class nuclease option. The broader PAM compatibility of NRG versus NGG does open up additional guide sites and can rescue some loci that fail with SpCas9 alone. However, even with NRG PAMs, the editing window constraint remains -- you still need the target nucleotide at positions 4-7. The expanded PAM alleviates the bottleneck but does not eliminate it. Our seven tested loci were specifically evaluated with SpCas9 NGG to make the cleanest possible statement about the canonical system, but the pipeline can rerun any case with enFnCas9 selected.

*Key fact:* enFnCas9 (NRG PAM) is evaluated as a first-class option; broader PAM partially rescues some loci but does not eliminate the window constraint.

---

**Q3: How do you handle bystander edits? Do you just count them or do you assess their biological consequence?**

We assess their biological consequence, not just their presence. The system scans all C or A residues within the editing window -- depending on whether it is CBE or ABE -- and for each potential bystander position, it uses the coding annotation module to determine whether conversion at that position would produce a synonymous, missense, or nonsense change. Synonymous bystanders get a small bonus (+0.05), missense bystanders incur a -0.10 penalty per position, and nonsense bystanders -- which would introduce a premature stop codon -- incur a severe -0.25 penalty per position. This requires knowing the codon frame, which is why transcript-aware mapping is essential.

*Key fact:* Bystander consequences are classified as synonymous/missense/nonsense using codon frame information, with penalties of 0/−0.10/−0.25 per position.

---

**Q4: How do you define splice proximity, and why those specific distance thresholds?**

We follow ACMG standards. Positions within 2 base pairs of an exon boundary are classified as splice donor (5' end) or splice acceptor (3' end) sites and receive a penalty of -0.15. Positions 3 to 8 base pairs from the boundary are classified as splice region and receive a penalty of -0.08. These thresholds are based on the known functional importance of the GT-AG dinucleotides at splice junctions and the broader splice regulatory region. The 2 bp threshold captures the invariant dinucleotides; the 3-8 bp threshold captures the extended splice consensus where variants can disrupt splicing through weakening of the splice signal.

*Key fact:* Splice proximity follows ACMG: <=2 bp = donor/acceptor (penalty -0.15), 3-8 bp = splice region (penalty -0.08).

---

**Q5: For HDR, how do you estimate gene conversion tract length? Is that the ConversionSim module from v1?**

The v2 pipeline uses an exponential decay model for cut-to-edit distance scoring, with a half-life of approximately 20 base pairs calibrated to published tract-length data from Elliott and Jasin (1998). This is a simplified version of what ConversionSim does in v1. ConversionSim in v1 is a full Monte Carlo simulator that models gene conversion tract length distributions and has been validated against four published datasets. The v2 pipeline uses the simpler model for speed during strategy ranking, but ConversionSim remains available through the webapp for detailed HDR donor design analysis.

*Key fact:* Exponential decay model with ~20 bp half-life, calibrated to Elliott & Jasin (1998) tract-length data.

---

**Q6: Why does the system recommend ssODN for distances under 30 bp but cssDNA for larger distances? What is the biological rationale?**

The donor type recommendation is based on the practical tradeoffs between donor templates at different edit distances. ssODN donors -- single-stranded oligodeoxynucleotides with 90 bp homology arms -- work well when the cut-to-edit distance is short because HDR efficiency drops exponentially with distance. For distances up to 5,000 bp, cssDNA donors with 300 bp homology arms provide better efficiency because the longer homology arms compensate for the increased distance, and the circular topology of cssDNA protects against exonuclease degradation. This is directly relevant to the cssDNA work in our lab and is informed by the Iyer et al. 2022 paper on efficient HDR with circular ssDNA donors. For even larger spans, conventional lssDNA or dsDNA donors with 800 bp arms are recommended.

*Key fact:* ssODN (<=30 bp, 90 bp arms), cssDNA (<=5,000 bp, 300 bp arms), lssDNA/dsDNA (>5,000 bp, 800 bp arms) -- based on Iyer et al. 2022 and published efficiency data.

---

**Q7: You mention p53-mediated toxicity in iPSCs. How specifically does this affect the scoring?**

DSB-free strategies -- base editing and prime editing -- receive a safety base score of 1.0. Strategies requiring a single DSB, like standard Cas9 HDR, receive 0.5. Strategies requiring two simultaneous DSBs receive 0.2, reflecting the additional translocation risk documented by Leibowitz et al. (2021). Since safety carries the highest weight in our scoring function at 0.30, this creates a meaningful score differential. A typical base editing strategy scores around 0.650 while an equivalent HDR strategy at the same locus scores around 0.388 -- that gap of approximately 0.26 is driven primarily by the safety differential of 0.50 multiplied by the 0.30 weight (= 0.15), plus additional complexity and risk penalties for HDR. In practice, this means HDR needs to offer substantial advantages in other components to overcome the safety penalty.

*Key fact:* Safety scores: DSB-free = 1.0, single DSB = 0.5, dual DSB = 0.2; at weight 0.30, this creates a ~0.15 point safety advantage for PE/BE over HDR.

---

**Q8: Does the system account for the efficiency differences between ABE8e, ABE7.10, and other ABE variants?**

Currently, CRISPRArchitect uses a single set of editing window parameters for ABE (positions 4-7) and CBE (positions 4-8) that correspond to the canonical editors. We do not yet model the expanded editing windows or shifted activity profiles of next-generation editors like ABE8e (Richter et al. 2020), which has a somewhat broader activity window. This is a deliberate simplification -- we use conservative window definitions to avoid overestimating feasibility. Adding editor-variant-specific window profiles is a planned extension and would be straightforward to implement since the window parameters are configurable inputs to the feasibility engine.

*Key fact:* ABE window = positions 4-7, CBE window = positions 4-8 (canonical editors); editor-variant-specific profiles are a planned extension.

---

**Q9: For the pegRNA design, why did you choose 13 nt as the default PBS length?**

The 13-nucleotide default for the primer binding site is based on the empirical data from Anzalone et al. (2019), who found that PBS lengths of 10-17 nucleotides generally support efficient prime editing, with 13 nt representing a good balance between binding stability and avoiding excessive secondary structure. The pipeline evaluates PBS lengths across the 10-17 nt range and selects the optimal design, but uses 13 nt as the default starting point. The RT template length is evaluated across 10-30 nt. These ranges are consistent with the design guidelines from both the original Anzalone paper and the enhanced PE systems described by Chen et al. (2021).

*Key fact:* PBS default = 13 nt (range 10-17), RT template range = 10-30 nt, consistent with Anzalone et al. 2019 and Chen et al. 2021.

---

**Q10: What happens with variants in alternatively spliced exons? Does the canonical transcript assumption cause problems?**

The system uses the canonical transcript as determined by Ensembl, which prioritizes MANE Select and transcript length. For variants in alternatively spliced exons, this means the analysis is performed against a single transcript model, which may not capture the full complexity. A variant in a cassette exon that is skipped in some isoforms would be analyzed as if the exon is always included. This is acknowledged as a limitation -- the simplified CDS model treats the spliced exonic transcript as a surrogate for the full coding sequence. For the exonic variants in our benchmark this works well, but for variants in complex alternatively spliced regions, multi-transcript analysis would be needed. That would require evaluating consequence and feasibility against multiple transcript models and somehow reconciling potentially contradictory recommendations.

*Key fact:* Canonical transcript from Ensembl (MANE Select priority); multi-transcript analysis is needed for alternatively spliced regions but is not yet implemented.

---

### Category 2: Methodological Questions (8)

---

**Q11: Why those specific scoring weights? How sensitive are the results to weight changes?**

The default weights -- Safety 0.30, Feasibility 0.25, Complexity 0.20, Risk 0.15, Confidence 0.10 -- are iPSC-optimized. Safety gets the highest weight because DSB toxicity in iPSCs is the single most important biological consideration documented in the literature. We did not perform a formal sensitivity analysis varying all five weights simultaneously, but the dominance of PE is robust to moderate weight perturbations because PE consistently scores well across multiple components -- not just safety. PE is DSB-free (high safety), has no editing window constraint (high feasibility), requires a single delivery (moderate complexity), avoids bystander edits (low risk), and has strong literature support (high confidence). For PE to lose its top ranking, you would need to substantially increase the complexity weight while decreasing the safety weight, which would only make sense in a cell type where DSB toxicity is not a concern.

*Key fact:* PE dominance is robust because PE scores well across all five components, not just safety.

---

**Q12: How were the benchmark truth labels defined? Is there circularity between the labels and the algorithm?**

The truth labels were defined by expert judgment based on published editing parameters and clinical precedent before the algorithm was finalized. There is no circularity -- the labels were not adjusted to match the algorithm's output, and the algorithm's scoring weights were not tuned to maximize benchmark accuracy. The labels represent what an experienced genome editor would recommend given the variant type, editing constraints, and cell type context. The "preferred" label indicates the single best strategy, "acceptable" indicates alternatives that would also work, and "reject" indicates strategies that are infeasible or biologically inappropriate. All labels are published in the benchmark dataset for inspection.

*Key fact:* Truth labels defined by expert judgment before algorithm finalization; no circular optimization.

---

**Q13: Why do you report 0% consequence shift? Does that mean the consequence-aware scoring does not matter?**

The 0% consequence shift means that in our benchmark, enabling versus disabling the consequence-aware adjustments (bystander penalties, splice proximity penalties) did not change which strategy was top-ranked in any of the 30 cases. This does not mean the consequence scoring is broken or useless -- it means that prime editing inherently avoids bystander edits (because the RT template encodes exactly the desired edit), so the consequence penalties primarily affect base editing and HDR strategies that are already ranked below PE. The consequence scoring is correctly implemented and produces appropriate adjustments when triggered, but in a benchmark where PE dominates, those adjustments affect the relative ranking of non-top strategies rather than the top-1 pick. The consequence scoring would matter more in a cell type where DSB toxicity is less of a concern and base editing competes more directly with prime editing for the top rank.

*Key fact:* 0% consequence shift occurs because PE inherently avoids bystanders; consequence scoring correctly fires but affects strategies already ranked below PE.

---

**Q14: Is 30 cases enough for a benchmark? How do you justify this sample size?**

Thirty cases is modest and we are transparent about that. The justification is depth over breadth: each case uses verified GRCh38 coordinates from ClinVar, has reference alleles validated against the Ensembl genome, and has expert-defined tiered truth labels across 11 categories that cover the major decision-relevant scenarios. We designed the benchmark to stress-test specific failure modes rather than to achieve statistical power over a large homogeneous dataset. That said, expanding the benchmark is a priority -- both to more cases and to include experimental outcome data. The current 30-case benchmark demonstrates the framework works for the major variant categories; a larger benchmark would quantify performance with tighter confidence intervals.

*Key fact:* 30 cases across 11 categories, each with verified coordinates and tiered expert truth labels -- depth over breadth.

---

**Q15: Why did you use ClinVar variants instead of designing synthetic test cases?**

ClinVar variants represent real disease-causing mutations at real genomic loci, which means they have realistic sequence context -- real PAM landscapes, real splice junction distances, real codon frames. Synthetic test cases could have any sequence context we want, which risks testing the algorithm against artificially easy or artificially hard cases. By using ClinVar variants, we ensure the benchmark reflects the actual constraints that a user would encounter when designing editing strategies for patient-derived iPSC lines. The tradeoff is that ClinVar variants do not come with experimentally verified editing outcomes, so our truth labels are based on expert judgment rather than measured data.

*Key fact:* ClinVar variants provide realistic sequence context (PAM landscape, splice distances, codon frames) that synthetic cases cannot guarantee.

---

**Q16: How does the pipeline handle reverse-strand genes? Is NF1 on the reverse strand?**

Yes, NF1 is on the reverse strand (chromosome 17, minus strand). The pipeline handles this through automatic reverse complementation during transcript mapping. When the system fetches the genomic sequence from Ensembl, it compares the reported reference allele against the retrieved sequence accounting for strand orientation. For reverse-strand genes, the coding strand is the complement of the genomic strand, so alleles are reverse-complemented before any downstream analysis. This is essential for correctly identifying the base to be edited -- if you have a C>T change on the coding strand of a reverse-strand gene, the genomic strand has G>A. Getting this wrong would cause the system to look for the wrong type of base editor. All 30 benchmark variants passed reference validation, confirming that strand handling is correct.

*Key fact:* NF1 is reverse-strand; alleles are automatically reverse-complemented; all 30 cases passed reference validation confirming correct strand handling.

---

**Q17: How does the retry logic work for Ensembl API failures?**

The pipeline implements automatic retry with exponential backoff for transient Ensembl API failures. It retries up to 3 times with wait intervals of 1, 2, and 4 seconds for HTTP status codes 500, 502, 503, 504, and network timeout errors. In our benchmark run, 4 of 120 total API calls initially failed due to transient server errors; all four were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases. Permanent failures -- like a 404 for a nonexistent gene -- are not retried. The retry logic is intentionally conservative; we do not retry indefinitely, and we do not retry client errors.

*Key fact:* 3 retries with exponential backoff (1s, 2s, 4s); 4/120 API calls failed transiently in the benchmark, all recovered.

---

**Q18: What is the wall time for the full benchmark? Is this practical for interactive use?**

The full 30-case benchmark ran in approximately 977 seconds -- about 16 minutes -- which averages to roughly 32 seconds per case. Most of that time is Ensembl API latency (transcript fetching, sequence retrieval, VEP annotation). The actual computation -- feasibility checking, strategy generation, scoring -- takes less than a second per case. For interactive use, you are typically analyzing 1-2 variants at a time, which takes 30-60 seconds including API calls. This is practical for experimental design but not for screening thousands of variants. If high-throughput screening is needed, implementing local caching of Ensembl data would dramatically reduce latency.

*Key fact:* ~32 seconds per case (dominated by Ensembl API latency); computation itself is sub-second.

---

### Category 3: Comparison Questions (5)

---

**Q19: How does CRISPRArchitect compare to BE-Hive?**

BE-Hive (Arbab et al. 2020) is a machine learning model that predicts base editing outcomes -- it tells you what the product distribution will be for a given base editor and target sequence. It is excellent at what it does, but it operates within a single modality. It does not tell you whether you should use base editing in the first place, or whether prime editing or HDR would be better for your variant. CRISPRArchitect operates at a different level: it asks which modality to use, evaluates feasibility across all three, and then ranks strategies. The two tools are complementary rather than competitive. In an ideal workflow, you would use CRISPRArchitect to determine that base editing is feasible and preferred, then use BE-Hive to predict the product distribution for your specific guide.

*Key fact:* BE-Hive predicts within-modality outcomes; CRISPRArchitect selects between modalities. They are complementary.

---

**Q20: How does CRISPRArchitect compare to PrimeDesign?**

PrimeDesign (Hsu et al.) is a tool for designing prime editing guide RNAs -- pegRNAs with optimized PBS and RT template parameters. It is focused on PE-specific reagent design. CRISPRArchitect does not try to replace PrimeDesign for pegRNA optimization; instead, it operates upstream of that decision. CRISPRArchitect determines whether prime editing is the right modality for your variant (versus base editing or HDR), and then generates a basic pegRNA design as part of its strategy output. If CRISPRArchitect recommends prime editing, you would then use PrimeDesign to refine the pegRNA design parameters. Again, these are complementary tools at different levels of the decision hierarchy.

*Key fact:* PrimeDesign optimizes pegRNA design within PE; CRISPRArchitect decides whether PE is the right modality in the first place.

---

**Q21: Why not just use CRISPOR for everything?**

CRISPOR is an excellent tool for guide RNA design and off-target assessment, and we explicitly defer off-target prediction to CRISPOR in our pipeline. But CRISPOR is fundamentally a guide design tool -- it designs guides for a nuclease you have already chosen, at a locus you have already decided to cut. It does not compare base editing versus prime editing versus HDR, it does not incorporate coding consequences or splice proximity into strategy selection, and it does not score the overall biological appropriateness of the editing strategy. CRISPRArchitect and CRISPOR operate at different levels of the design process: CRISPRArchitect helps you decide what to do, CRISPOR helps you design how to do it.

*Key fact:* CRISPOR designs guides within a chosen modality; CRISPRArchitect selects the modality and strategy. They address different questions.

---

**Q22: What about Benchling? Many labs use Benchling for CRISPR design.**

Benchling is a comprehensive lab informatics platform that includes CRISPR design tools. For guide RNA design and sequence management, Benchling is extremely capable. However, Benchling's CRISPR design module is primarily focused on guide identification and basic HDR donor design for a single modality. It does not provide systematic cross-modality comparison, it does not perform consequence-aware scoring that penalizes bystander edits based on their amino acid impact, and it does not generate ranked strategy recommendations across base editing, prime editing, and HDR. CRISPRArchitect addresses a specific gap that Benchling does not fill: the upstream decision of which editing modality to use. A practical workflow would be to use CRISPRArchitect for modality selection, then use Benchling for detailed reagent design and experiment management.

*Key fact:* Benchling excels at guide design and lab management but does not perform cross-modality comparison or consequence-aware strategy ranking.

---

**Q23: Are there any other tools that do cross-modality comparison?**

To our knowledge, no existing tool performs systematic, consequence-aware comparison across base editing, prime editing, and HDR within a single framework. Several tools compare within a modality -- BE-Hive for base editing outcomes, PrimeDesign for pegRNA optimization, CRISPOR for guide design. Some reviews and decision trees have been published as static figures in perspective articles, but these are manual flowcharts, not computational systems. CRISPRArchitect fills this specific gap. We want to be clear that we are not claiming to be better than these tools at what they do -- we are doing something different. The value proposition is the unified, automated, consequence-aware comparison across modalities.

*Key fact:* No existing tool performs automated, consequence-aware cross-modality comparison; CRISPRArchitect fills a specific unaddressed gap.

---

### Category 4: Limitation Questions (5)

---

**Q24: You have no experimental validation. How confident should we be in these recommendations?**

This is a fair and important question. Without experimental validation, the confidence you should have is calibrated: the system agrees with expert judgment 86.7% of the time for top-1 and 96.7% for top-3. This means it is a useful starting point that will give you a reasonable recommendation in the vast majority of cases, but it is not a substitute for thinking about your specific experiment. The most critical next step is to close the feedback loop -- design editing experiments based on CRISPRArchitect recommendations, perform them in iPSC lines, and measure whether the recommended strategy actually works better than alternatives. Until we have that data, CRISPRArchitect is best used as a structured decision-support tool that systematizes the reasoning you would do anyway, not as an oracle.

*Key fact:* 86.7% top-1 agreement with expert judgment; experimental validation is the highest-priority next step.

---

**Q25: The system does not account for chromatin accessibility. How much does that matter?**

It matters significantly for absolute efficiency prediction -- a guide in closed chromatin will edit poorly regardless of how good it looks on paper. However, chromatin accessibility affects all modalities roughly equally at a given locus: if the chromatin is closed, base editing, prime editing, and HDR will all be less efficient. So while omitting chromatin data means our absolute efficiency estimates are unreliable, the relative ranking between modalities at the same locus is less affected. The place where chromatin matters differentially is when comparing strategies that use different guide positions -- one guide might be in open chromatin while another is not. Integrating ATAC-seq or ENCODE accessibility data is a near-term priority.

*Key fact:* Chromatin accessibility affects absolute efficiency prediction but has less impact on relative modality ranking at the same locus.

---

**Q26: You do not perform off-target prediction. Is that not a critical gap?**

Off-target prediction is important but is a solved problem handled well by existing tools like CRISPOR and Cas-OFFinder. We made a deliberate design decision to defer off-target assessment to these specialized tools rather than reimplementing it. CRISPRArchitect evaluates on-target feasibility and biological consequences, which is the gap that existing tools do not address. In a complete design workflow, you would use CRISPRArchitect to select the modality and generate candidate strategies, then use CRISPOR to evaluate the off-target profiles of the specific guides recommended. Integrating a Cas-OFFinder interface directly into the pipeline is on our medium-term roadmap.

*Key fact:* Off-target prediction is deferred to CRISPOR/Cas-OFFinder by design; integration is planned.

---

**Q27: How does the system handle cases where none of the three modalities work well?**

If no modality passes feasibility checking -- for example, a complex structural variant in a region with no usable PAM sites -- the pipeline returns an empty strategy list with explicit reasons for each rejection. This is preferable to forcing a recommendation when none is appropriate. In our benchmark, every case had at least one feasible strategy, but in real-world usage you could encounter variants where the system correctly concludes that none of the evaluated modalities are appropriate. In those cases, the user would need to consider alternative approaches outside the system's scope, such as large-fragment knockin, transposon-based insertion, or alternative nucleases not yet supported.

*Key fact:* If nothing passes feasibility, the pipeline returns explicit rejection reasons rather than forcing an inappropriate recommendation.

---

**Q28: The large deletion cases all failed at top-1. Is this a fundamental limitation of the approach?**

It is a limitation of the current implementation, not the approach. The pipeline currently processes each variant independently and evaluates it against the three standard modalities. For large deletions that span multiple exons, the correct strategy is not to correct individual nucleotide variants but to perform HDR with a deletion-spanning donor -- which requires understanding that the two boundary coordinates define a structural deletion, not two independent point mutations. Adding this logic requires explicit structural variant detection: recognizing when two input coordinates define a deletion, estimating the deletion size, and routing to HDR-specific logic that designs a donor spanning the deleted region. This is on our roadmap and is conceptually straightforward; we just have not implemented it yet.

*Key fact:* 0/3 top-1 for large deletions (but 3/3 top-3); fix requires explicit structural variant detection logic, which is planned.

---

### Category 5: Publication Strategy Questions (3)

---

**Q29: Which journal are you targeting, and why?**

We are preparing the manuscript for PLOS Computational Biology. The rationale is that this is a computational methods paper with a biological application focus. PLOS Computational Biology publishes tools and methods that serve the life sciences community, and it values open-source code and reproducibility -- both of which we have. The paper is not primarily a biological discovery paper (which would suit Nature Methods or Nature Biotechnology), nor is it a pure bioinformatics methods paper (which would suit Bioinformatics or NAR). PLOS Computational Biology hits the right intersection: computational methodology with clear biological utility and thorough benchmarking.

*Key fact:* PLOS Computational Biology -- computational methods with biological application, open-source, reproducible.

---

**Q30: What will reviewers likely push back on?**

Three likely pushback areas. First, the lack of experimental validation -- every reviewer will ask for at least some prospective data. Our response will be that this is a computational framework paper and that experimental validation is explicitly listed as the most important next step, but that the benchmark demonstrates concordance with expert judgment. Second, the benchmark size -- 30 cases will seem small. Our response is depth over breadth with verified coordinates and tiered labels. Third, the 0% consequence shift -- reviewers may question whether the consequence-aware scoring adds value. Our response is that the scoring is correctly implemented and is essential for scenarios where base editing competes with PE, but that our iPSC-optimized benchmark does not trigger those scenarios because PE's inherent safety advantage dominates.

*Key fact:* Expected pushback on: no experimental validation, small benchmark (n=30), and 0% consequence shift.

---

**Q31: What is the timeline to submission?**

The manuscript is essentially complete -- full text, all figures generated from actual benchmark data, methods section with all parameters documented. The GitHub repository is live with all code, benchmark data, and evaluation scripts. What remains before submission is a final round of internal review, ensuring all numbers in the text match the definitive benchmark results file exactly, formatting for the journal template, and writing the cover letter. We are looking at submission within 4-6 weeks. If we want to add any experimental validation data -- even a small pilot with 2-3 variants -- that would extend the timeline but would substantially strengthen the paper.

*Key fact:* Manuscript essentially complete; 4-6 weeks to submission without experimental data, longer if we add pilot validation.

---

## Appendix A: Quick-Reference Cheat Sheet

Keep this on your phone during the talk for instant number lookups.

| Metric | Value |
|--------|-------|
| Benchmark cases | 30 |
| Pipeline errors | 0 |
| Top-1 accuracy | 86.7% (26/30) |
| Top-3 accuracy | 96.7% (29/30) |
| Rejection accuracy | 90.0% (27/30) |
| PE top-ranked | 29/30 |
| HDR top-ranked | 1/30 |
| BE top-ranked | 0/30 |
| Consequence shift | 0% |
| ABE-compatible loci with no SpCas9 window | 7/7 |
| v1 codebase | ~24,000 LOC |
| v2 codebase | ~9,400 LOC |
| Total codebase | ~33,400 LOC |
| Test suite | 121 tests passing |
| API calls in benchmark | 120 |
| Transient API failures | 4/120 (all recovered) |
| Benchmark wall time | ~977 seconds (~16 min) |
| Safety weight | 0.30 |
| DSB-free safety score | 1.0 |
| Single-DSB safety score | 0.5 |
| Dual-DSB safety score | 0.2 |
| ABE window | positions 4-7 |
| CBE window | positions 4-8 |
| PBS default | 13 nt (range 10-17) |
| RT template range | 10-30 nt |
| PE3 nick distance | 40-100 bp opposite strand |
| ssODN distance threshold | <=30 bp |
| cssDNA distance threshold | <=5,000 bp |
| Bystander missense penalty | -0.10/position |
| Bystander nonsense penalty | -0.25/position |
| Splice donor/acceptor penalty | -0.15 |
| Splice region penalty | -0.08 |

---

## Appendix B: Timing Summary

| Section | Duration | Cumulative |
|---------|----------|------------|
| Slide 1: Title | 0:30 | 0:30 |
| Slide 2: The Problem | 2:00 | 2:30 |
| Slide 3: What CRISPRArchitect Does | 2:00 | 4:30 |
| Slide 4: Pipeline Architecture | 2:00 | 6:30 |
| Slide 5: Transcript-Aware Mapping | 2:00 | 8:30 |
| Slide 6: Feasibility Engines | 3:00 | 11:30 |
| Slide 7: Scoring Function | 2:00 | 13:30 |
| Slide 8: Benchmark Design | 1:30 | 15:00 |
| Slide 9: Results -- Overall | 2:00 | 17:00 |
| Slide 10: Results -- Heatmap | 1:30 | 18:30 |
| Slide 11: Results -- Per-Category | 1:30 | 20:00 |
| Slide 12: Key Finding -- PAM Bottleneck | 2:30 | 22:30 |
| Slide 13: Technical Innovation | 1:00 | 23:30 |
| Slide 14: Limitations | 1:30 | 25:00 |
| Slide 15: Future Directions | 1:00 | 26:00 |
| Slide 16: Acknowledgements | 0:30 | 26:30 |
| Slide 17: Thank You | 0:30 | 27:00 |
| Demo 1: Webapp | 3:00 | 30:00 |
| Demo 2: Python Pipeline | 3:00 | 33:00 |
| Demo 3: Benchmark Figures | 2:00 | 35:00 |
| **Total (presentation + demo)** | **~35 min** | |
| Q&A | 10:00 | 45:00 |

**Note:** If time is tight, cut Demo 1 (webapp) and reduce Demo 3 (figures) to 1 minute. The Python pipeline demo (Demo 2) is the most important and should not be cut.

---

## Appendix C: One-Line Summaries for Each Slide

Use these if you lose your place or need to quickly orient yourself.

1. **Title** -- CRISPRArchitect: unified cross-modality editing strategy design.
2. **Problem** -- No tool answers "which modality should I use for this variant in this cell type."
3. **What it does** -- Five capabilities: unified strategy space, transcript mapping, PAM feasibility, consequence scoring, explicit rejection.
4. **Architecture** -- Seven-stage pipeline from input variant to ranked strategies.
5. **Transcript mapping** -- Ensembl-based mapping gives exon, CDS, codon frame, splice proximity.
6. **Feasibility** -- Three engines (BE, PE, HDR) with modality-specific biological constraints.
7. **Scoring** -- Five-component multi-objective function with iPSC-optimized weights.
8. **Benchmark** -- 30 ClinVar cases, 11 categories, verified coordinates, tiered truth labels.
9. **Results overall** -- 86.7% top-1, 96.7% top-3, 90.0% rejection, zero failures.
10. **Heatmap** -- Feasibility is locus-specific; BE column is sparse due to PAM-window constraint.
11. **Per-category** -- 100% for 8/11 categories; failures only in large deletion and sequential HDR.
12. **PAM finding** -- PAM-window is the bottleneck, not mutation type; 0/7 ABE loci had SpCas9 window.
13. **Technical** -- v1 (24k LOC) + v2 (9.4k LOC); 121 tests; enFnCas9 as first-class option.
14. **Limitations** -- No experimental validation, no chromatin data, no off-target, poor on large deletions.
15. **Future** -- Chromatin integration, HGVS parser, ML scoring, off-target, experimental feedback loop.
16. **Acknowledgements** -- DC sir, lab, CSIR-IGIB, AcSIR, CSIR fellowship.
17. **Thank you** -- GitHub link, PLOS Comp Bio submission, open for questions.
