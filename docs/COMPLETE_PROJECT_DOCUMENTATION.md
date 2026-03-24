# CRISPRArchitect: Complete Project Documentation

**Transcript-aware, consequence-guided design of genome editing strategies across modalities**

Authors: Vishal Bharti and Debojyoti Chakraborty
Institution: CSIR-Institute of Genomics and Integrative Biology (CSIR-IGIB), New Delhi, India
Repository: https://github.com/visvikbharti/CRISPRArchitect
License: MIT
Version: v2 (March 2026)

---

## Table of Contents

1. [Project Genesis and Motivation](#1-project-genesis-and-motivation)
2. [Scientific Background](#2-scientific-background)
3. [Project Architecture (v1 + v2)](#3-project-architecture-v1--v2)
4. [Complete File Structure](#4-complete-file-structure)
5. [Data Flow (End-to-End)](#5-data-flow-end-to-end)
6. [Scoring System (Detailed)](#6-scoring-system-detailed)
7. [Benchmark Design and Execution](#7-benchmark-design-and-execution)
8. [Results and Interpretation](#8-results-and-interpretation)
9. [Validation Summary](#9-validation-summary)
10. [Limitations (Honest and Detailed)](#10-limitations-honest-and-detailed)
11. [Future Directions](#11-future-directions)
12. [Technical Details](#12-technical-details)

---

## 1. Project Genesis and Motivation

### 1.1 Where the Idea Came From

CRISPRArchitect was born in the Debojyoti Chakraborty laboratory at CSIR-IGIB, New Delhi, a group with deep expertise in CRISPR-based genome editing in human induced pluripotent stem cells (iPSCs). The Chakraborty lab is best known for developing **enFnCas9**, a variant of FnCas9 engineered for broadened PAM recognition (NRG instead of NGG), which substantially expands the targetable space of the genome. The lab routinely corrects pathogenic mutations in patient-derived iPSCs as part of disease modeling and therapeutic development pipelines.

The practical challenge that motivated CRISPRArchitect arose from a recurring experimental scenario: a patient-derived iPSC line carries **compound heterozygous mutations** -- two different pathogenic variants on the two alleles of the same gene. To create an isogenic control or a therapeutic cell product, both mutations must be corrected. But how? Should each mutation be corrected by base editing (if it is a compatible transition)? Should prime editing be used (if base editing windows are unavailable)? Should HDR be used with a cssDNA donor (if the mutations are close enough for single-template correction)? What if one mutation is a transition and the other is a transversion -- should a hybrid strategy (base editing for one, prime editing for the other) be employed?

### 1.2 The Specific Problem

When correcting compound heterozygous mutations in iPSCs, the researcher faces a combinatorial decision space: for n mutations, each potentially correctable by base editing (ABE or CBE), prime editing, or HDR (with ssODN, cssDNA, lssDNA, or dsDNA donors), the number of possible strategies grows rapidly. Moreover, the choice is not purely technical -- it depends on:

- **Sequence context**: Is a PAM available that positions the target base within the editing window? For base editing, the target nucleotide must fall within a narrow 4-nucleotide window (positions 4-7 for ABE, 4-8 for CBE) relative to the protospacer. This is a hard constraint that cannot be overcome by better reagents.
- **Biological consequences**: Will bystander edits introduce missense or nonsense changes? Is the variant near a splice site? Will the correction create unintended coding changes?
- **Safety profile**: iPSCs are exquisitely sensitive to double-strand breaks (DSBs) because they have active p53 pathways. DSBs trigger p53-mediated apoptosis, and surviving clones may have acquired p53 mutations or chromosomal rearrangements. DSB-free approaches (base editing, prime editing) are therefore strongly preferred in this cell type.
- **Practical complexity**: How many editing rounds, donors, and screening colonies are required? Sequential approaches are safer but slower; simultaneous approaches are faster but riskier.

No single existing tool addresses this decision space.

### 1.3 Why No Existing Tool Does This

The genome editing computational tool landscape in 2024-2026 is fragmented by modality:

- **BE-Hive** (Arbab et al., Cell, 2020): Predicts base editing outcomes (bystander editing profiles) for a given guide, but evaluates only base editing. It does not compare BE to PE or HDR, does not perform transcript-aware annotation, and does not score biological consequences.
- **PrimeDesign** (Hsu et al., Nature Biotechnology, 2021): Designs pegRNAs for prime editing targets, but evaluates only prime editing. It does not assess whether base editing might be simpler for a compatible transition, or whether HDR with a cssDNA donor might be preferable for a multi-site correction.
- **CRISPOR** (Concordet and Haeussler, Nucleic Acids Research, 2018): Excellent guide RNA design tool with off-target prediction, but it does not recommend editing strategies. It helps you design a guide once you have decided on a strategy, but does not help you decide between BE, PE, and HDR.
- **CRISPick** (Doench et al.): Guide scoring for CRISPR knockout screens, not editing strategy recommendation.

None of these tools provides a **unified framework** that (a) evaluates all three major editing modalities side by side, (b) incorporates transcript-level context (exon structure, coding frame, splice proximity), (c) scores downstream biological consequences (bystander mutations, splice disruption), and (d) produces ranked recommendations with transparent reasoning.

### 1.4 The Hypothesis

CRISPRArchitect tests a specific hypothesis:

> **A unified, transcript-aware, consequence-guided computational framework will make better editing strategy recommendations than modality-specific tools used in isolation.**

"Better" is operationalized as: when the framework's top-ranked strategy is compared against expert-defined truth labels (determined by biological reasoning from the published literature), the agreement rate (top-1 accuracy) exceeds that achievable by simple heuristic rules (e.g., "always use PE" or "use BE if it's a transition, otherwise PE").

The framework is designed to be transparent (every score component is interpretable), extensible (new nucleases, cell types, or modalities can be added), and grounded in published biology (every parameter has a literature citation).

---

## 2. Scientific Background

### 2.1 Base Editing Biology

Base editors are fusion proteins that combine a catalytically impaired Cas protein (nickase or dead Cas9) with a nucleotide deaminase enzyme. They convert one base to another without introducing a double-strand break.

**Adenine Base Editors (ABE)**: Developed by Gaudelli et al. (Nature, 2017). ABE converts adenine (A) to inosine (I), which is read as guanine (G) by the cellular machinery. Net effect: A-to-G conversion (or T-to-C on the complementary strand). The most widely used variant is ABE8e (Richter et al., Nature Biotechnology, 2020).

- Editing window: positions 4-7 within the 20-nt protospacer (1-indexed from the PAM-distal end)
- The target adenine must be within this window for efficient editing
- Other adenines within the window may also be edited (bystander editing)

**Cytosine Base Editors (CBE)**: Developed by Komor et al. (Nature, 2016). CBE converts cytosine (C) to uracil (U), which is read as thymine (T). Net effect: C-to-T conversion (or G-to-A on the complementary strand).

- Editing window: positions 4-8 within the 20-nt protospacer
- Slightly wider window than ABE
- Other cytosines within the window are bystander risks

**Key constraint**: Base editors can only perform transition mutations (purine-to-purine or pyrimidine-to-pyrimidine). ABE does A>G; CBE does C>T. Transversions (e.g., A>C, G>T) cannot be corrected by base editing. Furthermore, even for compatible transitions, a suitable PAM must exist at the precise spacing to place the target base within the editing window -- a constraint that CRISPRArchitect reveals is more restrictive than commonly appreciated.

**Bystander risk**: Any same-type base (A for ABE, C for CBE) within the editing window may be edited along with the target. If a bystander edit falls in a coding region, it may introduce a missense or nonsense change. CRISPRArchitect classifies all bystander consequences and penalizes strategies with deleterious bystanders.

### 2.2 Prime Editing Biology

Prime editing was developed by Anzalone et al. (Nature, 2019) and represents a fundamentally different approach: instead of chemically converting a base, prime editing uses a reverse transcriptase fused to a Cas9 nickase to directly write new genetic information into the genome.

**pegRNA design**: The prime editing guide RNA (pegRNA) contains three functional elements:
1. **Spacer** (20 nt): Directs Cas9 nickase to the target site
2. **Primer Binding Site (PBS)** (10-17 nt, default 13 nt): Hybridizes to the nicked strand to prime reverse transcription
3. **Reverse Transcriptase (RT) template** (10-30 nt): Encodes the desired edit plus flanking sequence

**PE3 nicking**: To improve efficiency, a second nicking guide is placed 40-100 bp away on the opposite strand. This creates a nick that biases mismatch repair toward incorporating the edit. CRISPRArchitect searches for suitable PE3 nicking guides automatically.

**Advantages over base editing**:
- Can install any substitution (transitions AND transversions)
- Can install small insertions (up to ~40 bp) and deletions (up to ~80 bp)
- No editing window constraint -- the RT template directly encodes the edit
- No bystander risk (only the intended edit is written)
- No DSB (only nicks), making it safer than HDR

**Limitations**: Lower efficiency than base editing for compatible transitions; pegRNA design complexity; sensitivity to PBS/RT template length optimization.

### 2.3 HDR Biology

Homology-Directed Repair (HDR) is the classical approach for precise genome editing. It requires:

1. **A DSB at or near the target site**: Created by Cas9 (or another nuclease). The DSB activates the cellular DNA repair machinery.
2. **A donor template**: Provides the desired sequence flanked by homology arms. The cell copies information from the donor into the genome.

**The repair process (SDSA pathway)**:
1. **End resection**: After the DSB, 5'-to-3' exonucleases (MRE11, EXO1) chew back the 5' ends, creating 3' single-stranded overhangs (~200-2,000 bp)
2. **RAD51 filament formation**: RAD51 protein coats the single-stranded DNA, forming a nucleoprotein filament
3. **Strand invasion**: The RAD51 filament searches for and invades the donor template at the homology arm, forming a D-loop
4. **DNA synthesis**: DNA polymerase extends the invading strand, copying the donor sequence (this copied region is the "gene conversion tract")
5. **Displacement (SDSA)**: The newly synthesized strand is displaced from the donor and re-anneals to the other side of the break

**Gene conversion tract**: The length of copied sequence follows a right-skewed, approximately geometric distribution. Median tract lengths for exogenous donors are ~500 bp (range 200-2,000 bp). The probability of incorporating an edit decreases with distance from the cut site, following an exponential decay with a half-life of ~20 bp.

**Donor types and their properties**:

| Donor Type | Homology Arms | Best For | Efficiency |
|------------|--------------|----------|------------|
| ssODN | 30-90 bp each | Edits within 30 bp of cut | Highest for proximal edits |
| cssDNA | 300 bp each | Edits within 5,000 bp of cut | ~2x better than lssDNA |
| lssDNA | 300 bp each | Moderate distance edits | Baseline |
| dsDNA | 800 bp each | Large edits, gene insertion | Lowest per-template |

**Cut-to-edit distance**: CRISPRArchitect uses an exponential decay model calibrated to published tract-length distributions (Elliott et al., MCB, 1998) to estimate the probability that the gene conversion tract will reach from the cut site to the desired edit position.

### 2.4 Why iPSCs Are Special

Human iPSCs present unique challenges for genome editing that directly influence strategy selection:

**p53 sensitivity**: iPSCs have active, wild-type p53 pathways. DSBs trigger p53-dependent apoptosis, killing the majority of edited cells (Ihry et al., Nature Medicine, 2018). This creates two problems:
1. **Low survival**: Most cells that receive a DSB die, reducing the number of correctly edited clones available for screening
2. **Selection for p53 mutations**: The surviving clones are enriched for cells that have acquired p53 loss-of-function mutations -- precisely the cells you do not want for disease modeling or therapy

**Karyotype concerns**: DSBs can cause chromosomal rearrangements (Leibowitz et al., Nature Genetics, 2021). Simultaneous DSBs at two loci on different chromosomes create a risk of reciprocal translocations. Even single DSBs can cause large deletions (>10 kb) around the cut site (Kosicki et al., Nature Biotechnology, 2018).

**Implication for strategy scoring**: CRISPRArchitect assigns the Safety component the highest weight (0.30) in its scoring function, directly reflecting the premium on DSB-free approaches in iPSC work. Strategies requiring zero DSBs (base editing, prime editing) receive a safety score of 1.0; strategies requiring one DSB (HDR) receive 0.5; strategies requiring two simultaneous DSBs receive 0.2.

### 2.5 enFnCas9: Broadened PAM from the Chakraborty Lab

enFnCas9 (engineered FnCas9) was developed in the Chakraborty laboratory at CSIR-IGIB. Unlike SpCas9, which requires an NGG PAM (where N is any nucleotide), enFnCas9 recognizes the broader NRG PAM (where R is A or G). This means enFnCas9 can target approximately 50% more sites in the genome compared to SpCas9.

In CRISPRArchitect, enFnCas9 is supported as a first-class nuclease alongside SpCas9. When the user specifies SpCas9 as the primary nuclease, the pipeline automatically also checks enFnCas9 compatibility for base editing, providing a broader view of the editable sequence space.

The broadened PAM partially alleviates the PAM-window bottleneck identified in this work -- at some loci where no SpCas9 guide places the target within the base editing window, an enFnCas9 guide may succeed. However, as the benchmark results demonstrate, broadened PAM does not fully eliminate this constraint.

---

## 3. Project Architecture (v1 + v2)

CRISPRArchitect is organized in two layers: the **v1 foundation** (~24,000 lines of code, 6 modules) and the **v2 extension** (~11,000 lines of code, 26 new files). The v1 modules provide validated biological simulations and strategy enumeration logic; the v2 modules add transcript-aware mapping, PAM-verified feasibility, consequence-aware scoring, and a unified pipeline.

### 3.1 v1 Modules (~24,000 LOC)

#### ConversionSim -- Monte Carlo HDR Gene Conversion Tract Simulator

**What it does**: Simulates the HDR process step by step using Monte Carlo methods. Given a DSB position, donor type, nuclease cut geometry, and cell type, ConversionSim runs thousands of virtual repair events to produce distributions of gene conversion tract lengths and incorporation probabilities at specified distances from the cut.

**Key classes**:
- `ResectionSimulator`: Models 5'-to-3' end resection (short-range by MRE11, long-range by EXO1)
- `FilamentModel`: Models RAD51 filament formation on the resected single-stranded DNA
- `SynthesisSimulator`: Models DNA polymerase extension along the donor template with stochastic displacement
- `ConversionSimulator`: Orchestrates the full simulation pipeline

**Biological basis**: Parameters are calibrated to published experimental data:
- Resection lengths from Symington (Annual Review of Genetics, 2011) and Cejka (Annual Review of Genetics, 2015)
- Tract length distributions from Elliott et al. (MCB, 1998)
- SDSA displacement probability fitted to tract data (~0.002 per bp)
- Donor topology multipliers from Iyer et al. (CRISPR Journal, 2022): cssDNA = 3.0, lssDNA = 1.5
- Staggered cut enhancement from Chauhan et al. (PNAS, 2023): +15% per bp of overhang

**Validation status**: Validated against 4 published datasets:
- Elliott et al. (1998): Tract length distribution shape -- PASS (qualitative shape match)
- Paquet et al. (2016): Distance-dependent incorporation -- POOR FIT (model over-predicts because it models SDSA, not SSTR which dominates for ssODN donors)
- Iyer et al. (2022): cssDNA vs lssDNA enhancement -- GOOD MATCH (predicted 2.07x, observed 1.9x, within experimental range 1.5-2.1x)
- Chauhan et al. (2023): Staggered cut enhancement -- GOOD MATCH (predicted 1.82x, observed 1.9x, within experimental range 1.4-2.8x)

#### MOSAIC -- Multi-locus Optimized Strategy for Allele-specific Integrated Correction

**What it does**: Given a gene structure (exon/intron architecture), mutation positions, cell type, and nuclease choice, MOSAIC enumerates every feasible editing strategy and scores them on efficiency, safety, time, and cost.

**Key classes**:
- `GeneStructure`: Represents the gene with exon/intron boundaries and lengths
- `Mutation`: Represents a single mutation with type classification (transition, transversion, insertion, deletion)
- `StrategyEnumerator`: Generates all possible strategies including SINGLE_TEMPLATE_HDR, DUAL_TEMPLATE_SIMULTANEOUS_HDR, SEQUENTIAL_HDR, BASE_EDITING, PRIME_EDITING, HYBRID approaches, and EXON_DELETION
- `StrategyScorer`: Scores each strategy on 4 axes (efficiency, safety, time, cost)
- `Reporter`: Generates human-readable reports with strategy comparisons

**Biological basis**:
- HDR efficiency estimates from cell-type-specific published rates
- Translocation risk from the polymer chain model (ChromBridge)
- p53-mediated selection risk for iPSCs (Ihry et al., 2018)
- Gene conversion probability from ConversionSim

**Validation status**: Benchmarked against 14 published papers. Author strategy appeared in MOSAIC top-3 in 10/14 cases (71.4% concordance). Perfect 100% accuracy on base editing cases (5/5). Disagreements occurred exclusively where authors used HDR or exon deletion while MOSAIC recommended DSB-free approaches -- defensible recommendations that reflect MOSAIC's safety-first scoring for iPSC context.

#### TopoPred -- cssDNA Secondary Structure Analyzer

**What it does**: Analyzes circular single-stranded DNA donor templates for secondary structures that could interfere with HDR. ssDNA is flexible and folds onto itself, forming G-quadruplexes and hairpins that block RAD51 binding and strand invasion.

**Key classes**:
- `GQuadruplexScanner`: Identifies G-quadruplex-forming sequences using the G3+N1-7G3+N1-7G3+N1-7G3+ pattern
- `HairpinPredictor`: Finds self-complementary regions and calculates thermodynamic stability using nearest-neighbor energy parameters
- `AccessibilityScorer`: Computes per-nucleotide accessibility scores for RAD51 binding
- `DonorOptimizer`: Suggests synonymous codon changes to disrupt deleterious structures

**Biological basis**: G-quadruplex stability from established thermodynamic models; hairpin energy calculations from the nearest-neighbor model; accessibility scoring based on the premise that structured regions are inaccessible to RAD51.

**Validation status**: Qualitative validation -- structures predicted by TopoPred match those predicted by ViennaRNA for test sequences. No systematic quantitative validation against experimental HDR data.

#### ChromBridge -- 3D Chromatin Distance Predictor

**What it does**: Calculates the physical 3D distance between two genomic loci using a polymer chain model of chromatin, and assesses whether a donor template can physically bridge two target sites.

**Key classes**:
- `PolymerModel`: Implements the Gaussian chain model for chromatin (mean squared displacement = N * b^2, where b = Kuhn length ~300 nm)
- `DistanceCalculator`: Converts genomic distance (bp) to physical distance (nm)
- `TADAnalyzer`: Assesses whether two loci are within the same topologically associating domain
- `TranslocationRiskModel`: Estimates translocation probability from 3D proximity using the empirical power law P(translocation) ~ s^(-1.08)

**Biological basis**: Polymer physics calibrated to FISH and Hi-C data; translocation risk model from Leibowitz et al. (Nature Genetics, 2021).

**Validation status**: Qualitative validation against published Hi-C contact frequency maps. The polymer model gives order-of-magnitude estimates; cell-type-specific Hi-C data would improve precision.

#### LoopSim -- Cohesin Loop Extrusion Simulator

**What it does**: Simulates the cohesin loop extrusion process that organizes chromatin into loops, and models how loop extrusion affects the spatial proximity between genomic loci.

**Key classes**:
- `ChromatinFiber`: Represents a chromatin region with nucleosome positions
- `CohesinExtruder`: Simulates the cohesin ring sliding along chromatin and forming loops
- `HomologySearchModel`: Models how RAD51-mediated homology search is affected by 3D chromatin organization
- `LoopSimulator`: Orchestrates the extrusion simulation

**Biological basis**: Cohesin extrusion dynamics from Fudenberg et al. (Cell Reports, 2016); CTCF boundary effects; loop formation kinetics.

**Validation status**: Qualitative validation. The loop extrusion simulator produces contact probability curves consistent with published Hi-C data patterns.

#### WebApp -- Streamlit Interactive Interface

**What it does**: Provides a web-based graphical interface for interacting with CRISPRArchitect. Users can input gene names, mutation positions, cell type, and nuclease choice, and receive strategy recommendations with interactive visualizations.

**Key files**:
- `app.py`: Main Streamlit application with v1 module interfaces
- `app_v2_page.py`: v2 pipeline interface page
- `style.py`: Custom styling for the web interface
- `run.sh`: Launch script

**Validation status**: Functional testing only. The web app is a presentation layer over the validated pipeline modules.

### 3.2 v2 Modules (~11,000 LOC)

#### core/models.py -- Central Data Models

**What it does**: Defines all shared dataclasses and enums used across the v2 pipeline. This is the single source of truth for data structures.

**Major dataclasses** (in dependency order):

| Dataclass | Purpose |
|-----------|---------|
| `ConsequenceType` (Enum) | Functional consequence classification: synonymous, missense, nonsense, splice_donor, splice_acceptor, splice_region, frameshift, inframe_insertion, inframe_deletion, intronic, 5'UTR, 3'UTR, non_coding, unknown |
| `EditModality` (Enum) | Editing modality: ABE, CBE, PE, HDR_ssODN, HDR_cssDNA, HDR_lssDNA, HDR_dsDNA, exon_deletion |
| `FeasibilityLabel` (Enum) | Hard verdict: feasible, marginal, not_feasible |
| `EvidenceTier` (Enum) | Confidence: A (PAM-verified), B (feasible with caveats), C (theoretical) |
| `RiskLevel` (Enum) | Rearrangement risk: low, moderate, high, very_high |
| `ExonRecord` | One exon: ID, number, start, end, strand, chromosome |
| `TranscriptInfo` | Full transcript: ID, gene symbol, chromosome, strand, exon list |
| `TranscriptCoordinate` | Mapped position: genomic position, exon number, CDS position, codon index/position, splice distance |
| `GenomicVariantInput` | Raw user input: chromosome, position, ref/alt alleles, gene symbol |
| `CodingAnnotation` | Consequence: type, HGVS notation, ref/alt codon and amino acid |
| `ReferenceValidation` | Ref allele check: valid/invalid, expected vs provided |
| `NormalizedVariant` | Fully annotated variant: combines input, transcript, coordinate, coding, ref validation |
| `GuideCandidate` | A candidate sgRNA: 20-mer, PAM, strand, cut position, GC content |
| `BaseEditingFeasibility` | BE result: label, editor type, best guide, bystander info |
| `PrimeEditingFeasibility` | PE result: label, best guide, PBS/RT lengths, PE3 nick guide |
| `HDRFeasibility` | HDR result: label, best guide, cut-to-edit distance, donor type |
| `FeasibilityBundle` | All modality results for one variant |
| `StrategyStep` | One step in a multi-step strategy |
| `Strategy` | Complete strategy: name, steps, DSB count, rounds, risk, evidence tier |
| `ScoredStrategy` | Strategy with scores: safety, feasibility, complexity, risk, confidence, overall |
| `PipelineResult` | Full pipeline output: transcript, variants, bundles, strategies |
| `BenchmarkCase` | One benchmark case: variants, truth labels, category |
| `BenchmarkResult` | One case's evaluation result |
| `BenchmarkSummary` | Aggregate metrics across all cases |

#### core/sequence/ -- Transcript Mapping Layer (5 modules)

**fetcher.py** (`TranscriptFetcher`): Connects to the Ensembl REST API to retrieve canonical transcript information for a gene. Fetches transcript ID, exon structure, coding sequence coordinates, strand, and biotype. Includes automatic retry logic with exponential backoff (3 retries with waits of 1, 2, 4 seconds) for transient server errors (HTTP 500, 502, 503, 504).

**transcript_mapper.py** (`TranscriptMapper`): Maps a genomic position to transcript context. Given a genomic coordinate and a TranscriptInfo object, determines which exon the position falls in, computes the CDS position, codon index, codon position (1st, 2nd, or 3rd position), and distance to exon boundaries. Handles both forward-strand and reverse-strand genes correctly.

**reference_validator.py** (`ReferenceValidator`): Validates that the user-provided reference allele matches the Ensembl genome sequence at the specified position. Fetches a window of genomic sequence from Ensembl and compares. For reverse-strand genes, reverse-complements appropriately. Catches genome-build mismatches and strand-orientation errors before they propagate.

**coding_annotation.py** (`CodingAnnotator`): Determines the coding consequence of a variant. Translates reference and alternate codons using the standard genetic code and classifies the consequence as synonymous, missense, nonsense, or frameshift. Annotates splice proximity following ACMG standards: within 2 bp of exon boundary = splice donor/acceptor, 3-8 bp = splice region.

**variant_normalizer.py** (`VariantNormalizer`): Orchestrates the full normalization pipeline: fetches transcript, maps position, validates reference allele, annotates coding consequence, fetches local genomic sequence context. Produces a NormalizedVariant dataclass that contains everything needed for downstream feasibility assessment.

#### core/feasibility/ -- Feasibility Engines (4 modules)

**pam_scan.py** (`EnhancedPAMScanner`): Scans both strands of a genomic sequence window for PAM sequences. Supports SpCas9 (NGG) and enFnCas9 (NRG). For each PAM found, extracts the 20-nt protospacer, calculates the cut position (3 bp upstream of PAM on protospacer strand), computes GC content, and checks for poly-T stretches (>=4 consecutive T's, which terminate Pol III transcription). Returns a ranked list of GuideCandidate objects.

**base_editing.py** (`BaseEditingEngine`): Evaluates base editing feasibility for a given variant. Three-step process: (1) determine if the correction is a compatible transition (A>G for ABE, C>T for CBE), (2) scan for PAM sites that position the target nucleotide within the editing window (positions 4-7 for ABE, 4-8 for CBE), (3) identify and classify bystander edits within the window. Returns BaseEditingFeasibility with label, best guide, bystander count, and bystander consequences.

**prime_editing.py** (`PrimeEditingEngine`): Evaluates prime editing feasibility. For each candidate guide near the target, designs a pegRNA with PBS (default 13 nt) and RT template (10-30 nt) that encodes the desired edit. Searches for a PE3 nicking guide 40-100 bp away on the opposite strand. Feasible for substitutions, insertions up to 40 bp, and deletions up to 80 bp. Returns PrimeEditingFeasibility with guide, PBS/RT lengths, and PE3 nick details.

**hdr_design.py** (`HDRDesignEngine`): Evaluates HDR feasibility. For each candidate guide, calculates cut-to-edit distance and scores it using an exponential decay function. Recommends donor type based on distance: ssODN for <=30 bp, cssDNA for <=5,000 bp, lssDNA or dsDNA for larger spans. Estimates gene conversion probability using the ConversionSim model. Recommends homology arm lengths: 90 bp for ssODN, 300 bp for cssDNA/lssDNA, 800 bp for dsDNA. Supports asymmetric donor design. Returns HDRFeasibility with all design parameters.

#### core/mosaic/ -- Strategy Generation and Scoring (2 modules)

**generator.py** (`StrategyGenerator`): Takes a list of FeasibilityBundle objects (one per variant) and generates all possible editing strategies. For single-variant cases, generates separate BE, PE, and HDR strategies. For multi-variant cases, additionally generates dual/hybrid strategies (dual BE, dual PE, hybrid BE+PE, hybrid BE+HDR, hybrid PE+HDR, sequential HDR) and evaluates which modality combinations are feasible. Marks infeasible strategies with rejection reasons.

**annotation_integration.py**: Integrates consequence annotations from the coding annotation module into strategy objects. Propagates bystander severity scores, splice proximity flags, and consequence penalties to the strategy layer so the scorer can apply biologically grounded adjustments.

#### core/pipeline/ -- Orchestrator (1 module)

**strategy_stage.py** (`StrategyPipeline`, `StrategyScorer`): This is the single entry point for the v2 pipeline. Contains two classes:

- `StrategyScorer`: Implements the multi-objective scoring function with 5 components (Safety, Feasibility, Complexity, Risk, Confidence) and consequence-aware adjustments (bystander penalties, splice proximity penalties, clean design bonuses). Ranks strategies by overall score with safety as tiebreaker.

- `StrategyPipeline`: Orchestrates the 8-stage pipeline from raw GenomicVariantInput to ranked ScoredStrategy list. Handles errors gracefully at each stage, falls back to minimal normalization if full annotation is unavailable, and automatically checks both SpCas9 and enFnCas9 for base editing.

#### benchmarks/ -- Evaluation Framework (4 files)

**dataset_v1.json**: The curated benchmark dataset of 30 variant scenarios with verified GRCh38 coordinates from ClinVar and published editing studies. Each case includes: case_id, gene_symbol, variants (with chromosome, position, ref/alt alleles), category, truth labels (preferred, acceptable, reject), disease context, source PMID, and rationale.

**evaluator.py** (`BenchmarkEvaluator`): Runs each benchmark case through the pipeline and evaluates the output against truth labels. Computes top-1 accuracy, top-3 accuracy, and rejection accuracy. Handles pipeline errors gracefully and reports them as failures.

**run_benchmark.py**: CLI runner for the benchmark. Accepts the dataset JSON as input, executes all cases, writes results to JSON, and prints a summary table. Includes wall-time measurement.

**plotting.py**: Generates publication-quality figures: architecture diagram (Fig. 1), feasibility heatmap (Fig. 3), category-level accuracy bar chart (Fig. 5), and overall benchmark summary (Fig. 7). Outputs both PDF and PNG formats.

#### tests/ -- Test Suite (4 files, 93 v2 tests)

**test_v2_models.py** (259 lines): Tests instantiation, default values, and serialization of all core dataclasses. Verifies enum values and property computations.

**test_feasibility.py** (752 lines): Tests PAM scanning, base editing feasibility, prime editing feasibility, and HDR design feasibility. Includes edge cases: transversions correctly rejected by BE engine, poly-T guides filtered, bystander consequences classified correctly.

**test_strategy_generation.py** (668 lines): Tests strategy generation for single-variant and multi-variant cases. Verifies that the correct strategy types are generated, infeasible strategies are properly rejected, and scoring produces expected rankings (BE/PE > HDR in iPSC context).

**test_conversion_sim.py** (309 lines): Tests the v1 ConversionSim module (30 tests for v1 validation).

---

## 4. Complete File Structure

Below is every file in the project with a one-line description.

### Root Files

| File | Description |
|------|-------------|
| `__init__.py` | Package initializer for the crisprarchitect package |
| `README.md` | Project overview, installation instructions, and quick start guide |
| `LICENSE` | MIT License text |
| `requirements.txt` | Python dependencies: numpy, scipy, matplotlib, seaborn, pandas, requests |
| `pyproject.toml` | Python project metadata and build configuration |
| `MANIFEST.in` | Specifies files to include in source distribution |
| `CITATION.cff` | Citation metadata in Citation File Format |
| `CONTRIBUTING.md` | Contributor guidelines |
| `RELEASE_NOTES.md` | Version history and release notes |
| `Dockerfile` | Docker container definition for reproducible deployment |
| `docker-compose.yml` | Docker Compose configuration for multi-service deployment |
| `cli.py` | Command-line interface entry point |
| `.gitignore` | Git ignore patterns |

### core/ -- v2 Pipeline Modules

| File | Description |
|------|-------------|
| `core/__init__.py` | Package initializer for v2 core modules |
| `core/models.py` | Central data models: 22 dataclasses and 5 enums shared across the pipeline |
| `core/sequence/__init__.py` | Package initializer for the sequence layer |
| `core/sequence/fetcher.py` | Ensembl REST API client for transcript retrieval with retry logic |
| `core/sequence/transcript_mapper.py` | Genomic-to-transcript coordinate mapper with CDS position resolution |
| `core/sequence/reference_validator.py` | Reference allele validator against Ensembl genome sequence |
| `core/sequence/coding_annotation.py` | Coding consequence annotator (synonymous, missense, nonsense, splice) |
| `core/sequence/variant_normalizer.py` | Full variant normalization orchestrator combining all sequence modules |
| `core/feasibility/__init__.py` | Package initializer for the feasibility layer |
| `core/feasibility/pam_scan.py` | Enhanced PAM scanner for SpCas9 (NGG) and enFnCas9 (NRG) |
| `core/feasibility/base_editing.py` | Base editing feasibility engine with ABE/CBE window and bystander analysis |
| `core/feasibility/prime_editing.py` | Prime editing feasibility engine with pegRNA and PE3 nicking guide design |
| `core/feasibility/hdr_design.py` | HDR feasibility engine with donor type recommendation and conversion probability |
| `core/mosaic/__init__.py` | Package initializer for the strategy layer |
| `core/mosaic/generator.py` | Strategy generator: enumerates BE, PE, HDR, and hybrid strategies |
| `core/mosaic/annotation_integration.py` | Integrates consequence annotations into strategy objects for scoring |
| `core/pipeline/__init__.py` | Package initializer for the pipeline |
| `core/pipeline/strategy_stage.py` | Pipeline orchestrator and multi-objective scoring engine (single entry point) |

### conversion_sim/ -- v1 Gene Conversion Tract Simulator

| File | Description |
|------|-------------|
| `conversion_sim/__init__.py` | Package initializer exposing ConversionSimulator |
| `conversion_sim/models.py` | Data models for simulation parameters and results |
| `conversion_sim/resection.py` | End resection simulator (MRE11 short-range + EXO1 long-range) |
| `conversion_sim/filament.py` | RAD51 filament formation model |
| `conversion_sim/synthesis.py` | DNA synthesis and SDSA displacement simulator |
| `conversion_sim/simulator.py` | Main simulation orchestrator with Monte Carlo loop |

### mosaic/ -- v1 Strategy Optimizer

| File | Description |
|------|-------------|
| `mosaic/__init__.py` | Package initializer exposing MOSAIC classes |
| `mosaic/gene_structure.py` | Gene structure representation with exon/intron architecture |
| `mosaic/mutation_classifier.py` | Mutation type classifier (transition, transversion, indel) |
| `mosaic/strategy_enumerator.py` | Enumerates all feasible multi-locus editing strategies |
| `mosaic/scorer.py` | Multi-axis strategy scorer (efficiency, safety, time, cost) |
| `mosaic/reporter.py` | Human-readable report generator with strategy comparisons |

### topopred/ -- v1 cssDNA Structure Analyzer

| File | Description |
|------|-------------|
| `topopred/__init__.py` | Package initializer exposing TopoPred classes |
| `topopred/g_quadruplex.py` | G-quadruplex sequence scanner and stability predictor |
| `topopred/hairpin.py` | Hairpin/stem-loop predictor with nearest-neighbor thermodynamics |
| `topopred/accessibility.py` | Per-nucleotide accessibility scorer for RAD51 binding |
| `topopred/optimizer.py` | Donor sequence optimizer suggesting synonymous changes to disrupt structures |

### chrombridge/ -- v1 3D Chromatin Distance Predictor

| File | Description |
|------|-------------|
| `chrombridge/__init__.py` | Package initializer exposing ChromBridge classes |
| `chrombridge/polymer_model.py` | Gaussian chain polymer model for chromatin |
| `chrombridge/distance.py` | Genomic-to-physical distance converter |
| `chrombridge/tad_analysis.py` | Topologically associating domain (TAD) boundary analyzer |
| `chrombridge/translocation.py` | Translocation risk estimator from 3D proximity |

### loopsim/ -- v1 Cohesin Loop Extrusion Simulator

| File | Description |
|------|-------------|
| `loopsim/__init__.py` | Package initializer exposing LoopSim classes |
| `loopsim/chromatin_fiber.py` | Chromatin fiber model with nucleosome positions |
| `loopsim/cohesin_extruder.py` | Cohesin ring loop extrusion dynamics |
| `loopsim/homology_search.py` | RAD51-mediated homology search in 3D chromatin context |
| `loopsim/simulator.py` | Main loop extrusion simulation orchestrator |
| `loopsim/visualize.py` | Visualization functions for loop extrusion results |

### utils/ -- Shared Utilities

| File | Description |
|------|-------------|
| `utils/__init__.py` | Package initializer for shared utilities |
| `utils/constants.py` | All biological parameters with literature citations (nuclease PAMs, cell type params, donor multipliers) |
| `utils/sequence.py` | DNA sequence tools (complement, reverse_complement, PAM finding, GC content) |
| `utils/plotting.py` | General-purpose visualization functions |
| `utils/ensembl.py` | Ensembl REST API helper with retry logic and error handling |

### benchmarks/ -- Evaluation Framework

| File | Description |
|------|-------------|
| `benchmarks/__init__.py` | Package initializer for benchmark modules |
| `benchmarks/dataset_v1.json` | 30 curated ClinVar variant scenarios with truth labels |
| `benchmarks/evaluator.py` | Benchmark evaluator computing top-1, top-3, and rejection accuracy |
| `benchmarks/run_benchmark.py` | CLI runner with JSON output and summary table |
| `benchmarks/plotting.py` | Publication-quality figure generator (Figs. 1, 3, 5, 7) |

### benchmark_results/ -- Stored Results

| File | Description |
|------|-------------|
| `benchmark_results/definitive_benchmark_results.json` | Full results: 30 cases, per-case outcomes, accuracy metrics |
| `benchmark_results/consequence_shift_analysis.json` | Consequence-aware vs naive scoring comparison (0% shift rate) |
| `benchmark_results/figures/fig1_architecture.png` | Pipeline architecture diagram |
| `benchmark_results/figures/fig1_architecture.pdf` | Pipeline architecture diagram (vector) |
| `benchmark_results/figures/fig3_feasibility_heatmap.png` | Feasibility heatmap across modalities and cases |
| `benchmark_results/figures/fig3_feasibility_heatmap.pdf` | Feasibility heatmap (vector) |
| `benchmark_results/figures/fig5_category_accuracy.png` | Per-category accuracy bar chart |
| `benchmark_results/figures/fig5_category_accuracy.pdf` | Per-category accuracy (vector) |
| `benchmark_results/figures/fig7_benchmark_summary.png` | Overall benchmark summary panel |
| `benchmark_results/figures/fig7_benchmark_summary.pdf` | Overall benchmark summary (vector) |

### tests/ -- Test Suite

| File | Description |
|------|-------------|
| `tests/__init__.py` | Package initializer for test modules |
| `tests/test_conversion_sim.py` | 30 tests for v1 ConversionSim module validation |
| `tests/test_v2_models.py` | Tests for all core dataclass instantiation and serialization |
| `tests/test_feasibility.py` | Tests for PAM scanning, BE/PE/HDR feasibility engines |
| `tests/test_strategy_generation.py` | Tests for strategy generation, scoring, and ranking |

### webapp/ -- Streamlit Interactive Interface

| File | Description |
|------|-------------|
| `webapp/app.py` | Main Streamlit application for v1 module interfaces |
| `webapp/app_v2_page.py` | v2 pipeline interface page with variant input and strategy display |
| `webapp/style.py` | Custom CSS and styling for the web interface |
| `webapp/run.sh` | Bash script to launch the Streamlit application |
| `webapp/requirements.txt` | Web app-specific dependencies (streamlit) |

### validation/ -- v1 Validation Reports and Scripts

| File | Description |
|------|-------------|
| `validation/VALIDATION_REPORT.md` | ConversionSim validation against 4 published datasets |
| `validation/MOSAIC_BENCHMARK_REPORT.md` | MOSAIC benchmark against 14 published papers |
| `validation/validate_conversionsim.py` | Script to reproduce ConversionSim validation |
| `validation/benchmark_mosaic.py` | Script to reproduce MOSAIC benchmark |
| `validation/figures/validation1_tract_distribution.png` | Tract length distribution plot |
| `validation/figures/validation2_distance_incorporation.png` | Distance-dependent incorporation plot |
| `validation/figures/validation3_cssdna_vs_lssdna.png` | cssDNA vs lssDNA comparison plot |
| `validation/figures/validation4_staggered_enhancement.png` | Staggered cut enhancement plot |

### paper/ -- Manuscript and Figures

| File | Description |
|------|-------------|
| `paper/CRISPRArchitect_v2_manuscript.md` | Full manuscript text (Markdown source) |
| `paper/CRISPRArchitect_PLOS_manuscript.docx` | PLOS Computational Biology formatted manuscript |
| `paper/CRISPRArchitect_PLOS_manuscript_with_figures.docx` | Manuscript with embedded figures |
| `paper/CRISPRArchitect_Supplementary.docx` | Supplementary materials |
| `paper/CRISPRArchitect_Cover_Letter.docx` | Cover letter for journal submission |
| `paper/cover_letter.md` | Cover letter (Markdown source) |
| `paper/cover_letter_v2_nature_methods.md` | Cover letter variant for Nature Methods |
| `paper/PLOS_CompBio_manuscript.md` | PLOS-formatted manuscript (Markdown) |
| `paper/manuscript_draft.md` | Earlier manuscript draft |
| `paper/supplementary_materials.md` | Supplementary materials (Markdown) |
| `paper/REFERENCE_VERIFICATION.md` | Verification log for all cited references |
| `paper/CRISPRArchitect_v2_presentation.md` | Presentation script (Markdown) |
| `paper/CRISPRArchitect_v2_LabMeeting.pptx` | Lab meeting PowerPoint presentation |
| `paper/generate_figures.py` | Script to generate all manuscript figures |
| `paper/generate_presentation_v2.py` | Script to generate presentation slides |
| `paper/figures/Fig1_ConversionSim.png` | Figure 1: ConversionSim overview |
| `paper/figures/Fig1_ConversionSim.pdf` | Figure 1 (vector) |
| `paper/figures/Fig1.tif` | Figure 1 (TIFF for journal) |
| `paper/figures/Fig2_Validation.png` | Figure 2: Validation results |
| `paper/figures/Fig2_Validation.pdf` | Figure 2 (vector) |
| `paper/figures/Fig2.tif` | Figure 2 (TIFF for journal) |
| `paper/figures/Fig3_MOSAIC_Benchmark.png` | Figure 3: MOSAIC benchmark |
| `paper/figures/Fig3_MOSAIC_Benchmark.pdf` | Figure 3 (vector) |
| `paper/figures/Fig3.tif` | Figure 3 (TIFF for journal) |
| `paper/figures/FigS1_Sensitivity.png` | Supplementary Figure 1: Parameter sensitivity |
| `paper/figures/FigS1_Sensitivity.pdf` | Supplementary Figure 1 (vector) |
| `paper/figures/FigS1.tif` | Supplementary Figure 1 (TIFF) |
| `paper/figures/FigS2_Weight_Sensitivity.png` | Supplementary Figure 2: Scoring weight sensitivity |
| `paper/figures/FigS2_Weight_Sensitivity.pdf` | Supplementary Figure 2 (vector) |
| `paper/figures/FigS2.tif` | Supplementary Figure 2 (TIFF) |

### docs/ -- Documentation

| File | Description |
|------|-------------|
| `docs/ARCHITECTURE.md` | System architecture and design document with biological explanations |
| `docs/HONEST_ASSESSMENT.md` | Candid assessment of project strengths and weaknesses |
| `docs/LIVE_DEMO_SCRIPT.md` | Step-by-step script for live demonstrations |
| `docs/PRESENTATION_GUIDE.md` | Guide for presenting CRISPRArchitect at conferences |
| `docs/ROADMAP.md` | Development roadmap and planned features |
| `docs/SESSION_CONTEXT.md` | Development session context and notes |
| `docs/USER_GUIDE.md` | User guide for running the pipeline |

### examples/ -- Example Scripts

| File | Description |
|------|-------------|
| `examples/fetch_real_genes.py` | Example: fetching real gene structures from Ensembl |
| `examples/your_scenario_analysis.py` | Template for analyzing your own editing scenario |

### Configuration Files

| File | Description |
|------|-------------|
| `.github/workflows/ci.yml` | GitHub Actions CI/CD workflow for automated testing |
| `.streamlit/config.toml` | Streamlit configuration for the web app |

---

## 5. Data Flow (End-to-End)

This section walks through exactly what happens when a user inputs a variant such as "NF1 c.910C>T" into CRISPRArchitect, from the raw input to the final ranked strategy list.

### Step 1: GenomicVariantInput Created

The user creates a `GenomicVariantInput` object specifying the variant in GRCh38 coordinates:

```python
variant = GenomicVariantInput(
    chromosome="17",
    position=31200443,
    ref_allele="C",
    alt_allele="T",
    gene_symbol="NF1",
    name="c.910C>T",
)
```

All coordinates are 1-based (Ensembl convention). The gene_symbol is required for transcript lookup. The name field is optional but improves readability in output reports.

### Step 2: TranscriptFetcher Calls Ensembl

The `TranscriptFetcher` makes an HTTP request to the Ensembl REST API (`https://rest.ensembl.org/lookup/symbol/homo_sapiens/NF1?expand=1`). It retrieves:

- Canonical transcript ID (e.g., ENST00000358273)
- Gene ID (e.g., ENSG00000196712)
- Chromosome, start, end positions
- Strand (-1 for NF1, which is on the reverse strand)
- Complete exon list with coordinates

Output: a `TranscriptInfo` object containing the full exon structure.

If the API call fails (HTTP 500, 502, 503, 504, or network timeout), automatic retry logic kicks in: up to 3 retries with exponential backoff (waits of 1, 2, and 4 seconds).

### Step 3: TranscriptMapper Maps Position

The `TranscriptMapper` receives the genomic position (31200443) and the TranscriptInfo object. It:

1. Iterates through exons (in transcript order, accounting for strand) to find which exon contains the position
2. Calculates the transcript position (position within the spliced exonic sequence)
3. Calculates the CDS position (position within the coding sequence, adjusting for the start codon)
4. Determines codon_index (which codon number) and codon_position (1st, 2nd, or 3rd position within the codon)
5. Extracts the reference codon and looks up the reference amino acid
6. Calculates distance to both exon boundaries (distance_to_exon_start and distance_to_exon_end) for splice proximity annotation

Output: a `TranscriptCoordinate` object.

### Step 4: ReferenceValidator Checks Ref Allele

The `ReferenceValidator` fetches a window of genomic sequence centered on position 31200443 from the Ensembl sequence API (`/sequence/region/homo_sapiens/17:31200343..31200543:1`). It then:

1. Extracts the nucleotide at the exact position
2. For reverse-strand genes like NF1, takes the complement (because the user specifies the ref allele on the coding strand, but Ensembl returns the forward strand)
3. Compares the expected reference allele against the user-provided "C"

If they match, the variant passes validation. If not, the variant is flagged and excluded from downstream analysis.

Output: a `ReferenceValidation` object with `is_valid=True` and a message confirming the match.

### Step 5: CodingAnnotator Determines Consequence

The `CodingAnnotator` uses the TranscriptCoordinate to:

1. Construct the reference codon from the CDS position and codon_position
2. Construct the alternate codon by substituting the alternate allele at the appropriate position
3. Translate both codons using the standard genetic code
4. Compare the amino acids to classify the consequence:
   - Same amino acid: synonymous
   - Different amino acid: missense
   - Stop codon introduced: nonsense
   - Insertion/deletion not divisible by 3: frameshift
5. Check splice proximity: distance_to_exon_start or distance_to_exon_end <=2 bp = splice_donor/splice_acceptor; 3-8 bp = splice_region

For NF1 c.910C>T, the consequence is p.Arg304Ter (nonsense -- a premature stop codon is introduced).

Output: a `CodingAnnotation` object.

### Step 6: VariantNormalizer Produces NormalizedVariant

The `VariantNormalizer` also fetches the local genomic sequence (+-200 bp flanking the variant position) and identifies the edit index within this local sequence. This local sequence is needed for PAM scanning and feasibility assessment.

Output: a `NormalizedVariant` object that bundles together: the original GenomicVariantInput, TranscriptInfo, TranscriptCoordinate, CodingAnnotation, ReferenceValidation, local_sequence, and local_seq_edit_index.

### Step 7: PAM Scanner Finds Guides

The `EnhancedPAMScanner` takes the 401-bp local sequence and scans both strands for PAM sequences:

- **SpCas9 (NGG)**: Scans for any two-base sequence followed by GG on the forward strand, or CC followed by any two bases on the reverse strand
- **enFnCas9 (NRG)**: Scans for any nucleotide, then A or G, then G (forward strand) -- broader recognition

For each PAM found, the scanner:
1. Extracts the 20-nt protospacer upstream of the PAM
2. Calculates the cut position (3 bp upstream of the PAM on the protospacer strand)
3. Computes distance from cut to the edit position
4. Calculates GC content (accepted range: 30-70%)
5. Checks for poly-T (>=4 consecutive T's -- rejected)
6. Computes a composite ranking score

Output: a list of `GuideCandidate` objects, sorted by score.

### Step 8: BE/PE/HDR Engines Evaluate Feasibility

Each feasibility engine receives the NormalizedVariant, local sequence, and edit index:

**Base Editing Engine**:
1. Determines if C>T correction requires CBE or T>C requires ABE (or the reverse). For NF1 c.910C>T, correction is T>C, which is ABE-compatible.
2. For each guide candidate, checks if the target base falls within the ABE window (positions 4-7) or CBE window (positions 4-8)
3. If a window-compatible guide exists, scans for bystander A's (ABE) or C's (CBE) within the window
4. Classifies each bystander's coding consequence
5. Returns BaseEditingFeasibility with label (FEASIBLE, MARGINAL, or NOT_FEASIBLE)

**Prime Editing Engine**:
1. For each guide near the edit, designs pegRNA parameters: PBS length (default 13 nt), RT template length (15 nt)
2. Searches for PE3 nicking guide on opposite strand, 40-100 bp away
3. Checks that the edit is within PE range (substitution, insertion <=40 bp, deletion <=80 bp)
4. Returns PrimeEditingFeasibility

**HDR Design Engine**:
1. For each guide, calculates cut-to-edit distance
2. Scores incorporation probability using exponential decay (half-life ~20 bp)
3. Recommends donor type: ssODN if <=30 bp, cssDNA if <=5,000 bp
4. Estimates required homology arm lengths
5. Returns HDRFeasibility

Output: a `FeasibilityBundle` containing results from all three engines.

### Step 9: StrategyGenerator Produces Strategies

The `StrategyGenerator` receives the list of FeasibilityBundle objects (one per variant) and generates all feasible strategies:

For a single variant:
- **Single-step Base Editing** (if BE is feasible): One ABE or CBE guide, no DSB, no donor
- **Single-step Prime Editing** (if PE is feasible): One pegRNA + PE3 nick, no DSB, no donor
- **Single-step HDR** (if HDR is feasible): One guide + one donor template, one DSB

For multiple variants, additional strategies are generated:
- **Dual Base Editing**: Both variants corrected by BE (if both are BE-feasible)
- **Dual Prime Editing**: Both variants corrected by PE
- **Hybrid BE+PE**, **Hybrid BE+HDR**, **Hybrid PE+HDR**: Mixed modality strategies
- **Sequential HDR**: Edit one site, clone, validate, edit the other

Strategies that fail hard constraints (e.g., two simultaneous DSBs with p53 active) are marked with rejection reasons.

Output: `List[Strategy]` including both viable and rejected strategies.

### Step 10: StrategyScorer Ranks Strategies

The `StrategyScorer` applies the multi-objective scoring function to each non-rejected strategy:

```
Score = w1*Safety + w2*Feasibility - w3*Complexity - w4*Risk + w5*Confidence
       - consequence_penalty + consequence_bonus
```

Where w1=0.30, w2=0.25, w3=0.20, w4=0.15, w5=0.10 (iPSC defaults).

Each component is computed as described in Section 6.

Strategies are sorted by overall score (descending), with safety score as tiebreaker. Each strategy is assigned a rank (1 = best).

Output: `List[ScoredStrategy]` sorted by rank.

### Step 11: PipelineResult Returned

The pipeline assembles everything into a `PipelineResult` object containing:

- `transcript`: The TranscriptInfo object
- `variants`: List of NormalizedVariant objects
- `bundles`: List of FeasibilityBundle objects
- `strategies`: List of ScoredStrategy objects (ranked, best first)
- `rejected_strategies`: List of Strategy objects that were rejected
- `metadata`: Pipeline execution metadata (cell type, nuclease, counts)
- `warnings`: Any warnings generated during execution

The caller can access the top strategy via `result.top_strategy`, iterate over all ranked strategies, or inspect rejected strategies to understand why they were excluded.

---

## 6. Scoring System (Detailed)

### 6.1 The Multi-Objective Scoring Equation

The overall score for each strategy is computed as:

```
Score = w1 * S_safety + w2 * S_feasibility - w3 * S_complexity - w4 * S_risk + w5 * S_confidence
        - consequence_penalty + consequence_bonus
```

The score is clamped to the range [0.0, 1.0].

**Default weights** (optimized for iPSC applications):

| Component | Symbol | Weight | Rationale |
|-----------|--------|--------|-----------|
| Safety | w1 | 0.30 | Highest weight because iPSCs have active p53; DSBs are toxic and select for p53-null clones |
| Feasibility | w2 | 0.25 | PAM availability and editing window compatibility are hard constraints |
| Complexity | w3 | 0.20 | Fewer rounds, donors, and screening colonies improve practical success |
| Risk | w4 | 0.15 | Rearrangement risk and bystander consequences affect biological outcome |
| Confidence | w5 | 0.10 | Evidence tier reflects design quality and literature support |

Note: Complexity and Risk are subtracted (they are penalties), while Safety, Feasibility, and Confidence are added (they are rewards). The weights are normalized to sum to 1.0 internally.

### 6.2 How Each Component Is Computed

#### Safety Score (0-1, higher is better)

DSB-free strategies are safest for iPSCs:

| Condition | Score | Biological Justification |
|-----------|-------|-------------------------|
| 0 DSBs (BE, PE) | 1.0 | No DSB means no p53 activation, no karyotype risk |
| 1 DSB, sequential, p53 active | 0.5 | Single DSB triggers p53 apoptosis (Ihry 2018) but no translocation risk |
| 1 DSB, sequential, p53 inactive | 0.6 | Less toxicity in p53-null cells |
| 2+ DSBs, simultaneous, p53 active | 0.1 | Dual DSBs risk translocation + severe p53 selection (Leibowitz 2021) |
| 2+ DSBs, simultaneous, p53 inactive | 0.2 | Translocation risk without p53 selection |

#### Feasibility Score (0-1, higher is better)

Product of two factors:
- `modality_prior_score`: Base confidence in the modality's ability to achieve the edit. PE = ~0.85, BE = ~0.90 (when feasible), HDR = ~0.72.
- `donor_feasibility_score`: Quality of the donor design (1.0 for no-donor strategies, discounted for large or complex donors).

Feasibility = modality_prior_score * donor_feasibility_score

#### Complexity Score (0-1, higher is worse -- penalizes complex designs)

Weighted combination of four penalty terms:

```
Complexity = 0.35 * rounds_penalty + 0.25 * donor_penalty + 0.20 * guide_penalty + 0.20 * screening_penalty
```

Where:
- `rounds_penalty` = min(1.0, (num_rounds - 1) * 0.3). Single-round strategies get 0; two-round strategies get 0.3.
- `donor_penalty` = min(1.0, num_donors * 0.15). No-donor strategies (BE, PE) get 0; one-donor HDR gets 0.15.
- `guide_penalty` = min(1.0, (num_distinct_guides - 1) * 0.1). Single-guide strategies get 0.
- `screening_penalty` = min(1.0, screening_clones / 100.0). 12 clones = 0.12; 48 clones = 0.48.

#### Risk Score (0-1, higher is worse -- penalizes risky designs)

Two components:

```
Risk = rearrangement_risk_score + bystander_severity * 0.3
```

Where rearrangement_risk_score maps from the RiskLevel enum:
- LOW: 0.0
- MODERATE: 0.3
- HIGH: 0.6
- VERY_HIGH: 0.9

And bystander_severity (0-1) reflects the worst-case bystander consequence within the editing window.

#### Confidence Score (0-1, higher is better)

Maps from the EvidenceTier enum:

| Tier | Score | Meaning |
|------|-------|---------|
| A | 1.0 | All components PAM-verified, editing window confirmed |
| B | 0.7 | Feasible but with caveats (bystanders, large donor) |
| C | 0.4 | Theoretical only (no PAM found, extrapolated efficiency) |

### 6.3 Consequence Penalties and Bonuses

Consequence-aware adjustments are applied after the main scoring equation:

#### Penalties (subtracted from score)

| Condition | Penalty | Biological Justification |
|-----------|---------|-------------------------|
| Bystander edit creating a missense change | -0.10 per position | Missense bystanders may alter protein function unpredictably |
| Bystander edit creating a nonsense change | -0.25 per position | Nonsense bystanders create premature stop codons; devastating for protein function |
| Variant within splice donor/acceptor (<=2 bp from exon boundary) | -0.15 | Splice site disruption typically causes exon skipping or intron retention |
| Variant within splice region (3-8 bp from exon boundary) | -0.08 | Splice region variants have moderate risk of affecting splicing |
| DSB burden: >=2 simultaneous DSBs with p53 active | -0.10 | Elevated risk of p53-mediated selection and chromosomal rearrangement |

The total consequence penalty is capped at 0.30 to prevent a single bad bystander from completely zeroing out an otherwise good strategy.

#### Bonuses (added to score)

| Condition | Bonus | Biological Justification |
|-----------|-------|-------------------------|
| All bystander edits are synonymous | +0.05 | Synonymous bystanders do not affect protein sequence; clean design |
| No DSBs AND no bystander risk | +0.03 | Cleanest possible design; maximally safe |

The total consequence bonus is capped at 0.10.

### 6.4 Why Safety Gets the Highest Weight for iPSC Work

The weight distribution is not arbitrary. It reflects a specific prioritization for iPSC-based disease modeling and therapeutic development:

1. **Patient safety**: iPSC-derived cells may eventually be transplanted into patients. Karyotypic abnormalities or p53 mutations acquired during editing could lead to tumorigenesis.
2. **Data integrity**: iPSC disease models are used to study gene function. If surviving clones are enriched for p53 mutations, the observed phenotype may reflect p53 loss rather than the intended correction.
3. **Regulatory requirements**: Clinical-grade iPSC products require extensive karyotype and whole-genome sequencing quality control. DSB-based approaches create more QC burden.

For other cell types (e.g., HEK293T, which has inactive p53), the weights can be adjusted to reduce safety emphasis and increase efficiency emphasis.

### 6.5 Example Calculation: BE vs PE vs HDR for a Specific Case

Consider a ClinVar variant that is a G>A transition in a coding exon, 50 bp from the nearest exon boundary. Suppose PAM scanning finds:
- An ABE-compatible guide that places the target A at position 5 (within window), with one bystander A at position 6 (predicted synonymous)
- A PE-compatible guide with a PE3 nicking guide 60 bp away
- An HDR-compatible guide with cut-to-edit distance of 8 bp

**Base Editing (ABE)**:
```
Safety = 1.0 (no DSB)
Feasibility = 0.90 * 1.0 = 0.90
Complexity = 0.35*0 + 0.25*0 + 0.20*0 + 0.20*0.12 = 0.024
Risk = 0.0 + 0.1*0.3 = 0.03 (one synonymous bystander, low severity)
Confidence = 1.0 (Tier A, PAM-verified)

Base score = 0.30*1.0 + 0.25*0.90 - 0.20*0.024 - 0.15*0.03 + 0.10*1.0
           = 0.300 + 0.225 - 0.005 - 0.005 + 0.100 = 0.615
Consequence bonus = +0.05 (all bystanders synonymous) + 0.03 (no DSB, clean) = +0.08
Final = 0.615 + 0.08 = 0.695
```

**Prime Editing**:
```
Safety = 1.0 (no DSB)
Feasibility = 0.85 * 1.0 = 0.85
Complexity = 0.35*0 + 0.25*0 + 0.20*0.1 + 0.20*0.12 = 0.044
Risk = 0.0 (no bystander risk)
Confidence = 1.0 (Tier A)

Base score = 0.30*1.0 + 0.25*0.85 - 0.20*0.044 - 0.15*0.0 + 0.10*1.0
           = 0.300 + 0.213 - 0.009 - 0.000 + 0.100 = 0.604
Consequence bonus = +0.03 (no DSB, no bystander)
Final = 0.604 + 0.03 = 0.634
```

**HDR (cssDNA)**:
```
Safety = 0.5 (one DSB, p53 active)
Feasibility = 0.72 * 0.95 = 0.684
Complexity = 0.35*0 + 0.25*0.15 + 0.20*0 + 0.20*0.12 = 0.062
Risk = 0.0 (low rearrangement risk) + 0.0 (no bystander)
Confidence = 0.7 (Tier B)

Base score = 0.30*0.5 + 0.25*0.684 - 0.20*0.062 - 0.15*0.0 + 0.10*0.7
           = 0.150 + 0.171 - 0.012 - 0.000 + 0.070 = 0.379
Consequence bonus = 0.0
Final = 0.379
```

**Ranking**: ABE (0.695) > PE (0.634) > HDR (0.379)

This illustrates how the safety premium drives DSB-free approaches to the top. In this hypothetical case, ABE beats PE because it has a slightly higher feasibility prior. In practice, CRISPRArchitect's benchmark revealed that ABE rarely achieves window placement at real ClinVar loci, so PE typically wins.

---

## 7. Benchmark Design and Execution

### 7.1 Why 30 Cases

The choice of 30 cases balances several considerations:

- **Statistical power**: 30 cases provide enough data to compute meaningful accuracy percentages (each case contributes ~3.3 percentage points). Smaller benchmarks would have high variance; a single case's outcome would swing the metric by >5%.
- **Category coverage**: 30 cases allow 11 categories with 1-7 cases each, covering the full spectrum of editing scenarios a researcher might encounter.
- **API feasibility**: Each case requires multiple Ensembl API calls (transcript fetch, reference validation, VEP annotation). With retry logic, the full benchmark takes ~16 minutes of wall time. Scaling to hundreds of cases would require rate-limiting management.
- **Manual curation quality**: Every case has manually verified GRCh38 coordinates, manually assigned truth labels with biological rationale, and manually confirmed reference alleles. This level of curation is labor-intensive and does not scale easily.

### 7.2 The 11 Categories and Why Each Matters

| Category | n | Why It Matters |
|----------|---|----------------|
| Clean base editable | 7 | Core use case: ClinVar transitions where BE should be feasible. Tests PAM-window verification. |
| Base editing negative | 1 | Sanity check: HBB sickle cell is a transversion (A>T). BE must be correctly rejected. |
| PE transversion | 2 | Transversions where only PE (not BE) can correct without DSB. Tests PE-specific design. |
| PE small indel | 5 | Small insertions/deletions within PE range (<=40 bp ins, <=80 bp del). Tests PE indel handling. |
| HDR large deletion | 3 | Multi-exon deletions too large for PE. HDR (or exon deletion) is the only option. Tests structural variant handling. |
| HDR/PE small deletion | 1 | Small deletion addressable by either HDR or PE. Tests modality comparison for ambiguous cases. |
| Compound het hybrid | 5 | Two mutations requiring different modalities. Tests hybrid strategy generation. |
| Sequential HDR | 1 | Compound het requiring two rounds of HDR. Tests multi-round strategy handling. |
| Dual base editing | 3 | Two transitions in the same gene. Tests dual BE strategy generation. |
| Edge: distant variants | 1 | Two variants very far apart in the same gene. Tests distance-aware strategy selection. |
| Edge: non-coding | 1 | Variant in 5'UTR (FMR1). Tests non-coding variant handling. |

### 7.3 How Truth Labels Were Assigned

Each case has three tiers of truth labels, assigned by biological reasoning (NOT by running the pipeline and accepting its output):

- **Preferred**: The strategy a domain expert would recommend as first-line. Based on: mutation type compatibility, PAM availability (verified by manual inspection), DSB-free preference for iPSCs, and published experimental precedent.
- **Acceptable**: Strategies that are biologically sound but suboptimal. For example, HDR is acceptable for a base-editable transition, but not preferred due to DSB risk in iPSCs.
- **Reject**: Strategies that are infeasible or dangerous. For example, base editing for a transversion, or simultaneous dual DSBs in iPSCs.

Example truth label for HBB sickle cell (c.20A>T, p.Glu7Val):
- Preferred: Prime Editing (transversion, DSB-free)
- Acceptable: HDR with ssODN (standard approach, but requires DSB)
- Reject: Base Editing (transversion cannot be corrected by ABE or CBE)

### 7.4 How Positions Were Verified

All 30 cases use GRCh38 genomic coordinates from ClinVar or published editing studies. Verification steps:

1. **ClinVar lookup**: For each variant, the ClinVar accession was checked to confirm the genomic coordinates, gene symbol, and clinical significance.
2. **Ensembl VEP**: Each position was submitted to the Ensembl Variant Effect Predictor (VEP) to verify that (a) the position maps to the expected gene, (b) the reference allele matches, and (c) the predicted consequence matches expectations.
3. **Reference allele validation**: The pipeline's own ReferenceValidator independently checks each position against the Ensembl genome. All 30 cases passed.

### 7.5 Retry Logic for API Robustness

Ensembl REST API calls occasionally fail with transient server errors. CRISPRArchitect handles this with:

- **Retry count**: Up to 3 retries per API call
- **Backoff schedule**: 1 second, 2 seconds, 4 seconds (exponential)
- **Retried errors**: HTTP 500 (Internal Server Error), 502 (Bad Gateway), 503 (Service Unavailable), 504 (Gateway Timeout), and network timeouts
- **Non-retried errors**: HTTP 400 (Bad Request), 404 (Not Found), 429 (Rate Limited -- these should be handled differently)

In the definitive benchmark run, 4 of 120 total API calls initially failed due to transient server errors. All 4 were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases.

---

## 8. Results and Interpretation

### 8.1 Raw Results

The definitive benchmark run (version 1.2) produced the following aggregate results across 30 cases:

| Metric | Value | Count |
|--------|-------|-------|
| **Top-1 Accuracy** | **86.7%** | 26/30 correct |
| **Top-3 Accuracy** | **96.7%** | 29/30 correct |
| **Rejection Accuracy** | **90.0%** | 27/30 correct |
| Pipeline Errors | 0% | 0/30 |
| Wall Time | 977 seconds | ~16 minutes |

**Strategy distribution** in top-1 rankings:

| Modality | Count | Percentage |
|----------|-------|------------|
| Prime Editing | 29 | 96.7% |
| HDR | 1 | 3.3% |
| Base Editing | 0 | 0.0% |

### 8.2 Per-Category Breakdown

| Category | n | Top-1 | Top-3 | Rejection |
|----------|---|-------|-------|-----------|
| Clean base editable | 7 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) |
| Base editing negative | 1 | 1/1 (100%) | 1/1 (100%) | 1/1 (100%) |
| PE transversion | 2 | 2/2 (100%) | 2/2 (100%) | 2/2 (100%) |
| PE small indel | 5 | 5/5 (100%) | 5/5 (100%) | 5/5 (100%) |
| HDR large deletion | 3 | 0/3 (0%) | 3/3 (100%) | 0/3 (0%) |
| HDR/PE small deletion | 1 | 1/1 (100%) | 1/1 (100%) | 1/1 (100%) |
| Compound het hybrid | 5 | 5/5 (100%) | 5/5 (100%) | 5/5 (100%) |
| Sequential HDR | 1 | 0/1 (0%) | 0/1 (0%) | 1/1 (100%) |
| Dual base editing | 3 | 3/3 (100%) | 3/3 (100%) | 3/3 (100%) |
| Edge: distant variants | 1 | 1/1 (100%) | 1/1 (100%) | 1/1 (100%) |
| Edge: non-coding | 1 | 1/1 (100%) | 1/1 (100%) | 1/1 (100%) |

### 8.3 Why PE Dominates

Prime editing emerged as the top-ranked modality in 29 of 30 cases. This is not a software bias but reflects the convergence of three genuine advantages within the iPSC-optimized scoring framework:

1. **No editing window constraint**: Unlike base editing, which requires the target base at positions 4-7 (ABE) or 4-8 (CBE) within the protospacer, prime editing encodes the desired edit directly in the RT template. Any edit within the PE range is feasible regardless of PAM-to-target distance.

2. **Broad mutation compatibility**: PE can correct transitions, transversions, and small indels. Base editing is restricted to transitions only (ABE: A>G, CBE: C>T).

3. **No DSB**: PE uses a nickase, not a nuclease. The safety score is 1.0 (same as BE), giving PE a massive advantage over HDR (safety score 0.5 in iPSC context). With Safety weighted at 0.30, this single factor creates a score gap of 0.30 * (1.0 - 0.5) = 0.15 -- often larger than all other score differences combined.

### 8.4 Why BE Never Tops

Despite 7 cases being designed as "clean base editable" transitions (ClinVar variants where ABE or CBE should be applicable), base editing never emerged as the top-1 strategy. The reason is the PAM-window bottleneck:

At all seven tested ClinVar loci with ABE-compatible transitions, **no SpCas9 NGG PAM site positioned the target nucleotide within the ABE editing window (positions 4-7)**. This means base editing was either infeasible or marginal at every one of these loci. Since PE has no such window constraint and shares the same safety advantage (no DSB), PE consistently outscored BE.

This finding -- that PAM-dependent editing window constraints are more restrictive than mutation-type classification -- is the most important scientific result of the benchmark. It suggests that published estimates of base editing applicability based on the fraction of ClinVar pathogenic variants that are transitions substantially overestimate the fraction that are practically editable with current base editors and canonical SpCas9.

### 8.5 The 4 Top-1 Misses

All four top-1 misses occurred in categories involving large deletions or sequential HDR:

1. **HDR_DMD_016** (DMD large deletion): The pipeline recommended PE, but the deletion spans multiple exons and is too large for PE. HDR or exon deletion is required. The pipeline treats the variant as a point deletion at the boundary, not as a structural rearrangement.

2. **HDR_NF1_017** (NF1 large deletion): Same issue. Multi-exon deletion mishandled as a single-variant problem.

3. **HDR_FBN1_022** (FBN1 large deletion): Same pattern. The pipeline lacks explicit multi-exon deletion handling logic.

4. **HDR_COL7A1_020** (COL7A1 compound het, sequential HDR): The expected strategy was sequential HDR (edit one allele, clone, then edit the other). The pipeline recommended PE for each variant independently, which is technically feasible for the individual variants but does not address the compound heterozygous coordination requirement.

The common thread: the pipeline currently processes variants independently and does not model structural variant-specific constraints or multi-variant coordination beyond strategy enumeration.

### 8.6 The 1 Top-3 Miss

The single top-3 miss was **HDR_COL7A1_020**: the COL7A1 compound heterozygous case requiring sequential HDR. The expected strategy (sequential HDR with two rounds of editing) was not generated because the pipeline processes each variant independently rather than as a coordinated pair requiring temporal ordering.

### 8.7 Consequence-Shift Analysis

A critical question: does consequence-aware scoring actually change strategy rankings compared to naive (consequence-unaware) scoring?

The answer from the benchmark is: **no, not in this dataset**. The consequence shift analysis shows:

- 28 of 30 cases had valid strategies for comparison (2 cases produced no strategies for one or both scoring modes)
- 27 of 28 cases (96.4%) had consequence adjustments applied (splice proximity penalties, bystander penalties, or clean design bonuses)
- **0 of 28 cases (0.0%)** had their top-1 ranking shifted by consequence adjustments

This means that while the consequence-aware scoring machinery is correctly implemented and consistently produces adjustments, those adjustments never change which strategy is ranked first. The reason is PE dominance: PE inherently avoids bystander edits (no editing window) and avoids DSBs (nick only), so the consequence penalties that differentiate BE and HDR strategies rarely apply to PE. Since PE already wins on safety, adding consequence penalties to BE/HDR only widens the gap.

This is a genuine limitation of the current benchmark composition and the iPSC weight configuration. In a dataset with more cases where BE and PE are both feasible and closely scored, consequence adjustments would be more likely to shift rankings.

### 8.8 The Key Finding: PAM-Window Bottleneck

The most important scientific finding from this work deserves emphasis:

**PAM-dependent editing window constraints are a more significant bottleneck for base editing applicability than mutation-type classification alone.**

Why this matters for the field:

1. **Overestimation of BE applicability**: Many papers cite statistics like "~60% of ClinVar pathogenic SNVs are transitions, and therefore potentially correctable by base editing." This number is based purely on mutation-type classification. The actual fraction that are practically editable requires checking that (a) a suitable PAM exists AND (b) the PAM positions the target within the 4-nt editing window. Our data suggest the practical fraction is substantially lower.

2. **Guide design is not trivial**: Finding a guide is necessary but not sufficient. The guide must position the target at exactly the right place. This is a constraint that cannot be appreciated from mutation databases alone -- it requires locus-specific sequence analysis.

3. **Expanded-PAM nucleases help but do not eliminate the problem**: enFnCas9 (NRG PAM) approximately doubles the number of targetable sites compared to SpCas9 (NGG), but even with NRG PAM, window placement is not guaranteed at every locus.

4. **Implications for therapeutic development**: For clinical applications of base editing, each target locus should be individually verified for PAM-window compatibility. CRISPRArchitect's automated PAM-window checking provides exactly this capability.

---

## 9. Validation Summary

### 9.1 v1 Validation

#### ConversionSim vs 4 Published Datasets

| # | Reference | Metric | Model Prediction | Published Value | Verdict |
|---|-----------|--------|-----------------|-----------------|---------|
| 1 | Elliott et al., MCB, 1998 | Tract length distribution shape | Right-skewed, geometric-like | Right-skewed, 80% <= 58 bp | PASS (qualitative shape match) |
| 2 | Paquet et al., Nature, 2016 | Distance-dependent incorporation | RMSE=0.41, R2=-0.56 | Monotonic decline | POOR FIT (model over-predicts; see note) |
| 3 | Iyer et al., CRISPR J, 2022 | cssDNA vs lssDNA ratio | 2.07x | 1.9x (range 1.5-2.1x) | GOOD MATCH |
| 4 | Chauhan et al., PNAS, 2023 | Staggered cut enhancement | 1.82x | 1.9x (range 1.4-2.8x) | GOOD MATCH |

**Note on Validation 2 (Paquet)**: The poor fit is expected and understood. Paquet used ssODN donors (~100-200 nt), which are incorporated largely through SSTR (single-strand template repair), a RAD51-independent pathway with much shorter tracts (~50 bp). ConversionSim models SDSA (mean ~500 bp tracts), which is the dominant pathway for longer donors (cssDNA, dsDNA with 300+ bp arms). The model is designed for cssDNA/dsDNA donors, not ssODN donors. A separate SSTR sub-model would be needed for ssODN validation.

**Simulations per validation**: 50,000 Monte Carlo runs per test case.
**Random seed**: 42 (for reproducibility).

#### MOSAIC vs 14 Published Papers

| Metric | Value |
|--------|-------|
| Papers benchmarked | 14 |
| Author strategy in MOSAIC top-3 | 10/14 (71.4%) |
| MOSAIC rank-1 matches | 7 papers |
| MOSAIC rank-2 matches | 1 paper |
| MOSAIC rank-3 matches | 2 papers |
| Misses (rank 4+) | 4 papers |

**Accuracy by strategy type**:

| Author Strategy | Papers | Hits | Accuracy |
|----------------|--------|------|----------|
| Base editing (ABE) | 5 | 5 | 100% |
| HDR (ssODN) | 5 | 3 | 60% |
| Prime editing | 1 | 1 | 100% |
| Exon skip via base editing | 1 | 1 | 100% |
| HDR (dsDNA) | 1 | 0 | 0% |
| Exon deletion (NHEJ) | 1 | 0 | 0% |

**Key pattern**: All 4 disagreements involved MOSAIC recommending a DSB-free approach (PE or BE) while authors used HDR or exon deletion. In each case, MOSAIC's recommendation is defensible as the safer approach for iPSCs. The disagreements reflect MOSAIC's safety-first prioritization versus authors who achieved high HDR efficiency through aggressive pharmacological optimization (p53 inhibition, HDR enhancers) or who published before PE was widely available.

### 9.2 v2 Validation

#### 30-Case ClinVar Benchmark on GRCh38

- **86.7% Top-1 accuracy** (26/30): The system's top-ranked strategy matches an expert-defined preferred or acceptable label
- **96.7% Top-3 accuracy** (29/30): At least one appropriate strategy appears in the top 3
- **90.0% Rejection accuracy** (27/30): Infeasible strategies correctly excluded from top-ranked output
- **0% pipeline failures**: All 30 cases processed without errors
- **100% reference allele validation**: All variants confirmed against the Ensembl genome

#### Biological Sanity Checks

The following sanity checks passed across all 30 cases:

- ABE is never recommended for a transversion mutation (verified: HBB sickle cell correctly gets PE, not ABE)
- Splice proximity is correctly detected and penalized when variants fall near exon boundaries
- HDR correctly penalized for DSB risk in iPSC context
- PE correctly assessed as feasible for all substitutions, small insertions, and small deletions within range
- Non-coding variants (FMR1 5'UTR) correctly identified as non-coding, with HDR recommended (consequence-aware scoring has limited applicability for non-coding regions)

---

## 10. Limitations (Honest and Detailed)

### Limitation 1: No Experimental Validation

**What it is**: CRISPRArchitect's recommendations are computational predictions that have not been validated in a prospective experimental design-outcome loop. The benchmark evaluates against expert-defined truth labels, not against actual editing outcomes in cells.

**Why it matters**: Computational feasibility does not guarantee experimental success. Many factors that influence editing efficiency (chromatin accessibility, cell cycle stage, delivery efficiency, RNA secondary structure) are not modeled. A strategy ranked first by the pipeline might yield lower efficiency than a strategy ranked third.

**How it could be fixed**: Design editing experiments based on CRISPRArchitect recommendations and compare predicted rankings against observed editing efficiencies. Use the experimental outcomes to refine the scoring weights. This is the most critical next step for the project.

### Limitation 2: Simplified CDS Model

**What it is**: The coding annotation module treats the spliced exonic transcript as a surrogate for the full coding sequence. It does not model alternative splicing, non-canonical reading frames, or overlapping genes.

**Why it matters**: For genes with complex transcript architectures (multiple isoforms, retained introns, non-AUG start codons), the simplified model may miss consequences that affect specific isoforms. A variant classified as synonymous in the canonical transcript might be missense or splice-disrupting in an alternative isoform.

**How it could be fixed**: Integrate MANE Select transcripts (which represent the consensus clinical transcript) and optionally evaluate multiple isoforms. Flag variants where consequence differs across isoforms.

### Limitation 3: No Chromatin Accessibility

**What it is**: CRISPRArchitect does not incorporate chromatin state (open vs. closed chromatin), replication timing, or epigenomic context. All loci are treated as equally accessible.

**Why it matters**: Editing efficiency varies dramatically by locus. Cas9 binding is strongly influenced by chromatin accessibility (open chromatin enables binding; heterochromatin blocks it). Two loci with identical PAM and guide scores may differ by 10-fold in actual editing efficiency. HDR rates are also influenced by replication timing (S/G2 phase is required for HDR).

**How it could be fixed**: Integrate ATAC-seq or DNase-seq data from iPSCs to compute per-locus accessibility scores. Use ENCODE chromatin state annotations (e.g., ChromHMM). Adjust feasibility scores by accessibility.

### Limitation 4: No Off-Target Prediction

**What it is**: CRISPRArchitect evaluates on-target feasibility and consequence but does not predict off-target editing sites for the recommended guides.

**Why it matters**: Off-target editing is a major safety concern, especially for therapeutic applications. A guide with perfect on-target design but numerous off-target sites would be a poor choice. Existing tools (Cas-OFFinder, CRISPOR, CRISPRscan) provide this analysis, but it is not integrated into CRISPRArchitect's scoring.

**How it could be fixed**: Integrate Cas-OFFinder or a similar alignment-based off-target predictor. Compute off-target scores for each candidate guide and incorporate them into the Risk component of the scoring function.

### Limitation 5: BE Applicability Lower Than Expected

**What it is**: At all seven tested ClinVar loci with ABE-compatible transitions, no SpCas9 guide placed the target within the editing window. Base editing therefore never emerged as the top-1 strategy.

**Why it matters**: This may give the impression that base editing is rarely useful, which is not the biological reality. The benchmark's 7 loci may not be representative of all ClinVar transitions. Furthermore, with expanded-PAM nucleases (enFnCas9), some of these loci might become BE-feasible.

**How it could be fixed**: Expand the benchmark to include ClinVar loci where BE feasibility has been experimentally confirmed (e.g., loci used in published BE studies). Systematically test enFnCas9 PAM availability at all loci. This limitation is a genuine biological constraint, not a software deficiency, but the benchmark composition amplifies its apparent severity.

### Limitation 6: Large Deletions Poorly Handled

**What it is**: All three HDR-required large deletion cases were top-1 misses. The pipeline recommended PE for individual variants at the deletion boundaries, which is technically feasible for the point mutations but does not address the multi-exon structural nature of the deletion.

**Why it matters**: Large deletions (>100 bp spanning multiple exons) require fundamentally different approaches: dual-guide excision followed by HDR, or exon skipping/deletion. The pipeline's current architecture processes variants as point mutations, not as structural rearrangements.

**How it could be fixed**: Add explicit structural variant detection: if the deletion spans multiple exons or exceeds PE range (80 bp), automatically flag it as requiring HDR or exon deletion. Generate dual-guide excision strategies for large deletions.

### Limitation 7: Consequence Scoring Doesn't Shift Rankings

**What it is**: In the 30-case benchmark, consequence-aware scoring adjustments were applied in 96.4% of cases but shifted the top-1 ranking in 0% of cases. The consequence-aware machinery is implemented and tested but had no practical impact on strategy selection.

**Why it matters**: This undermines the claim that consequence-aware scoring "influences strategy prioritization." While the machinery is correctly implemented, PE dominance in the iPSC scoring context means consequence penalties (which primarily affect BE and HDR) cannot overcome PE's inherent advantages.

**How it could be fixed**: (1) Design benchmark cases specifically to test consequence-shift scenarios (e.g., a locus where BE is feasible but has a severe bystander nonsense consequence, making PE preferable). (2) Test with non-iPSC weight configurations where safety has lower weight, allowing consequence differences to matter more. (3) Acknowledge that PE dominance makes consequence adjustments less impactful in iPSC context, which is itself an informative result.

### Limitation 8: Single-Gene Focus

**What it is**: CRISPRArchitect processes variants within a single gene. It cannot handle multi-gene editing scenarios (e.g., correcting mutations in two different genes on different chromosomes simultaneously).

**Why it matters**: Some disease models require multi-gene editing (e.g., creating a double-knockout, or correcting mutations in a gene and its regulatory element on a different chromosome). The pipeline's architecture assumes all variants share a single transcript.

**How it could be fixed**: Extend the pipeline to accept variants from multiple genes, fetch separate transcripts for each, and generate strategies that span genes. Multi-gene strategies would need special translocation risk assessment (ChromBridge integration).

### Limitation 9: No Indel Outcome Prediction

**What it is**: For HDR strategies, the pipeline does not predict the distribution of indel outcomes at the cut site. When HDR fails, the DSB is typically repaired by NHEJ, producing indels. The ratio of HDR to NHEJ and the indel spectrum are not predicted.

**Why it matters**: NHEJ outcomes determine the "background" of unedited cells. If NHEJ produces mainly in-frame deletions, the unedited cells may retain partial gene function. If NHEJ produces frameshifts, the unedited cells will have knockout alleles. This affects screening strategy and clone selection.

**How it could be fixed**: Integrate inDelphi (Shen et al., Nature, 2018) or a similar NHEJ outcome predictor. Compute the expected HDR:NHEJ ratio and incorporate the indel spectrum into the complexity/risk scoring.

### Limitation 10: Ensembl API Dependency

**What it is**: The v2 pipeline requires internet connectivity to access the Ensembl REST API for transcript information, reference validation, and VEP annotation. Offline use is not fully supported.

**Why it matters**: The Ensembl API is a shared resource with rate limits and occasional downtime. In a clinical setting, internet dependency introduces a point of failure. Ensembl database versions change over time, potentially affecting reproducibility.

**How it could be fixed**: Cache transcript and genome data locally. Provide a pre-built SQLite database with canonical transcripts for common genes. Allow offline mode that uses cached data, falling back to API only for novel genes. Pin the Ensembl release version in results metadata for reproducibility.

---

## 11. Future Directions

### 11.1 Experimental Validation Plan

The most critical next step. Proposed approach:
1. Select 10-15 ClinVar variants from the benchmark where the pipeline makes strong recommendations
2. Design editing reagents (BE, PE, HDR) according to CRISPRArchitect's recommendations AND according to alternative strategies
3. Edit iPSC lines with both the recommended and alternative strategies
4. Compare editing efficiency, specificity (bystander profile), and clone quality
5. Use the experimental outcomes to refine scoring weights

### 11.2 Chromatin Integration

Incorporate ATAC-seq and Hi-C data to make locus-specific editing efficiency predictions:
- ATAC-seq from iPSCs (ENCODE) for chromatin accessibility at guide target sites
- Hi-C for translocation risk refinement (replacing the generic polymer model)
- ChromHMM state annotations for context-dependent efficiency adjustment

### 11.3 ML-Based Scoring Refinement

Replace the heuristic scoring weights with a learned model:
- Train a gradient-boosted tree or neural network on experimental editing outcomes
- Features: safety score, feasibility score, complexity score, guide quality metrics, chromatin accessibility
- Label: editing efficiency and/or clone quality
- Requires experimental validation data (see 11.1)

### 11.4 HGVS Parser

Allow direct input of clinical variant nomenclature (e.g., "NM_000267.3:c.910C>T" or "NF1 p.Arg304Ter") without requiring the user to look up genomic coordinates. This requires:
- HGVS string parsing (using biocommons/hgvs library or custom parser)
- Transcript-to-genomic coordinate conversion
- Support for protein-level notation back-translation

### 11.5 Off-Target Integration

Integrate off-target prediction into the scoring framework:
- Run Cas-OFFinder for each candidate guide
- Compute off-target score based on number and severity of off-target sites
- Incorporate into the Risk component of the scoring function
- Flag guides with off-targets in known oncogenes or tumor suppressors

### 11.6 Clinical-Grade Reporting

Generate reports suitable for clinical review boards:
- Standardized format following ACMG/AMP guidelines
- Include variant classification, strategy rationale, risk assessment
- Comprehensive guide QC (GC content, poly-T, off-targets, specificity scores)
- Donor template sequences with annotated features
- Screening protocol recommendations (colony number, genotyping primers)

### 11.7 Multi-Gene Support

Extend the pipeline to handle variants across multiple genes:
- Separate transcript fetching per gene
- Cross-gene translocation risk assessment
- Multi-gene strategy enumeration (sequential editing of different genes)

---

## 12. Technical Details

### 12.1 Requirements

**Python version**: 3.9 or higher

**Core dependencies** (requirements.txt):
```
numpy>=1.24.0
scipy>=1.10.0
matplotlib>=3.7.0
seaborn>=0.12.0
pandas>=2.0.0
requests>=2.28.0
```

**Web app dependencies** (webapp/requirements.txt):
```
streamlit>=1.30.0
```

**No GPU required**. All computations are CPU-based. The Monte Carlo simulations in ConversionSim use NumPy vectorized operations for efficiency.

### 12.2 Installation

```bash
# Clone the repository
git clone https://github.com/visvikbharti/CRISPRArchitect.git
cd CRISPRArchitect/crisprarchitect

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "from core.models import GenomicVariantInput; print('OK')"
```

**Docker installation** (alternative):
```bash
docker-compose up --build
```

### 12.3 Running the Pipeline

```python
from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput

# Initialize pipeline for iPSC editing with SpCas9
pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")

# Define variant(s)
result = pipeline.run([
    GenomicVariantInput(
        chromosome="17",
        position=31200443,
        ref_allele="C",
        alt_allele="T",
        gene_symbol="NF1",
        name="c.910C>T",
    ),
])

# Access results
print(f"Top strategy: {result.top_strategy.strategy_name}")
print(f"Score: {result.top_strategy.overall_score:.3f}")

for strategy in result.strategies:
    print(f"  #{strategy.rank}: {strategy.strategy_name} "
          f"(score={strategy.overall_score:.3f}, "
          f"confidence={strategy.confidence})")

# Access rejected strategies
for rej in result.rejected_strategies:
    print(f"  REJECTED: {rej.name} — {rej.rejection_reasons}")
```

### 12.4 Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test files
python -m pytest tests/test_v2_models.py -v
python -m pytest tests/test_feasibility.py -v
python -m pytest tests/test_strategy_generation.py -v
python -m pytest tests/test_conversion_sim.py -v

# Run with coverage
python -m pytest tests/ --cov=core --cov=conversion_sim --cov-report=html
```

Test suite summary:
- **test_conversion_sim.py**: 30 tests for v1 ConversionSim validation
- **test_v2_models.py**: Tests for all core dataclass instantiation
- **test_feasibility.py**: Tests for PAM scanning and BE/PE/HDR engines
- **test_strategy_generation.py**: Tests for strategy generation and scoring
- **Total**: 93 v2 tests + 30 v1 tests = 123 total tests

### 12.5 Running the Benchmark

```bash
# Run the full 30-case benchmark (requires internet for Ensembl API)
python -m benchmarks.run_benchmark --input benchmarks/dataset_v1.json

# Results are saved to benchmark_results/
# Figures are generated in benchmark_results/figures/
```

The benchmark takes approximately 16 minutes to complete (977 seconds in the definitive run) due to Ensembl API calls. Each case requires 3-4 API calls (transcript fetch, reference validation, sequence fetch, and optionally VEP).

### 12.6 Running the Web App

```bash
# Launch the Streamlit interface
streamlit run webapp/app.py

# Or use the convenience script
bash webapp/run.sh
```

The web app provides:
- Gene name input with Ensembl lookup
- Interactive mutation definition (position, ref, alt alleles)
- Cell type and nuclease selection
- Real-time strategy ranking display
- Visualization of feasibility across modalities

### 12.7 Adding a New Nuclease

Edit `utils/constants.py` and add an entry to the `NUCLEASE_PARAMS` dictionary:

```python
"NewCas9": {
    "pam": "NNGG",              # PAM sequence (IUPAC notation)
    "cut_type": "staggered_5prime",  # "blunt" or "staggered_5prime" or "staggered_3prime"
    "stagger_bp": 3,            # Overhang length in bp (0 for blunt)
    "hdr_multiplier": 1.3,      # HDR efficiency multiplier relative to SpCas9
    "specificity": "high",      # "high", "medium", or "low"
    "description": "My new Cas9 variant with NNGG PAM",
}
```

Then update the PAM scanner in `core/feasibility/pam_scan.py` to handle the new PAM pattern.

### 12.8 Adding a New Cell Type

Edit `utils/constants.py` and add an entry to the `CELL_TYPE_PARAMS` dictionary:

```python
"my_cell_type": {
    "hdr_base_efficiency": 0.15,       # Baseline HDR efficiency (0-1)
    "cell_cycle_s_g2_fraction": 0.40,  # Fraction of cells in S/G2 (HDR-competent)
    "p53_active": True,                # Whether p53 pathway is active
    "viability_single_dsb": 0.65,      # Survival after single DSB
    "viability_dual_dsb": 0.40,        # Survival after dual simultaneous DSBs
    "description": "My custom cell type",
}
```

The scoring engine reads these parameters automatically to adjust safety scores and HDR feasibility estimates.

### 12.9 Adding New Benchmark Cases

Edit `benchmarks/dataset_v1.json` and add a new case entry:

```json
{
    "case_id": "NEW_GENE_031",
    "gene_symbol": "GENE",
    "category": "your_category",
    "variants": [
        {
            "chromosome": "1",
            "position": 12345678,
            "ref_allele": "A",
            "alt_allele": "G",
            "gene_symbol": "GENE",
            "name": "c.100A>G"
        }
    ],
    "cell_type": "iPSC",
    "nuclease": "SpCas9",
    "truth_label": {
        "preferred": ["Single-step Base Editing"],
        "acceptable": ["Single-step Prime Editing"],
        "reject": ["Simultaneous Dual HDR"]
    },
    "disease_context": "Disease name",
    "source_pmid": "12345678",
    "rationale": [
        "A>G transition is ABE-compatible",
        "PAM available at position X"
    ],
    "notes": "Verified against ClinVar accession VCV000012345"
}
```

Ensure:
1. GRCh38 coordinates are correct (verify against ClinVar and Ensembl VEP)
2. Reference allele matches the Ensembl genome
3. Truth labels are assigned by biological reasoning, not by running the pipeline
4. Rationale explains why each truth label tier was chosen

---

## Appendix A: Verified Numbers

These numbers have been verified against the actual code and benchmark results. They should be cited exactly as listed here.

| Metric | Value | Source |
|--------|-------|--------|
| v1 codebase | ~24,000 LOC | 6 modules across conversion_sim, mosaic, topopred, chrombridge, loopsim, webapp |
| v2 codebase | ~11,000 LOC | 26 new files across core/, benchmarks/, tests/ |
| Total v1 tests | 30 | test_conversion_sim.py |
| Total v2 tests | 93 | test_v2_models.py + test_feasibility.py + test_strategy_generation.py |
| Total tests | 123 | 30 + 93 |
| Benchmark cases | 30 | benchmarks/dataset_v1.json |
| Top-1 accuracy | 86.7% (26/30) | benchmark_results/definitive_benchmark_results.json |
| Top-3 accuracy | 96.7% (29/30) | benchmark_results/definitive_benchmark_results.json |
| Rejection accuracy | 90.0% (27/30) | benchmark_results/definitive_benchmark_results.json |
| Strategy distribution | PE=29, HDR=1, BE=0 | benchmark_results/definitive_benchmark_results.json |
| Consequence shift rate | 0.0% (0/28) | benchmark_results/consequence_shift_analysis.json |
| Cases with adjustments | 96.4% (27/28) | benchmark_results/consequence_shift_analysis.json |
| ConversionSim cssDNA ratio | 2.07x (Iyer observed: 1.9x) | validation/VALIDATION_REPORT.md |
| ConversionSim stagger ratio | 1.82x (Chauhan observed: 1.9x) | validation/VALIDATION_REPORT.md |
| MOSAIC concordance | 71.4% (10/14 papers) | validation/MOSAIC_BENCHMARK_REPORT.md |
| MOSAIC BE accuracy | 100% (5/5 papers) | validation/MOSAIC_BENCHMARK_REPORT.md |
| Scoring weight: Safety | 0.30 | core/pipeline/strategy_stage.py |
| Scoring weight: Feasibility | 0.25 | core/pipeline/strategy_stage.py |
| Scoring weight: Complexity | 0.20 | core/pipeline/strategy_stage.py |
| Scoring weight: Risk | 0.15 | core/pipeline/strategy_stage.py |
| Scoring weight: Confidence | 0.10 | core/pipeline/strategy_stage.py |
| Consequence penalty: splice donor/acceptor | -0.15 | paper/CRISPRArchitect_v2_manuscript.md |
| Consequence penalty: splice region | -0.08 | paper/CRISPRArchitect_v2_manuscript.md |
| Consequence penalty: bystander missense | -0.10 per position | paper/CRISPRArchitect_v2_manuscript.md |
| Consequence penalty: bystander nonsense | -0.25 per position | paper/CRISPRArchitect_v2_manuscript.md |
| Consequence bonus: all bystanders synonymous | +0.05 | paper/CRISPRArchitect_v2_manuscript.md |
| Benchmark wall time | 977 seconds (~16 min) | benchmark_results/definitive_benchmark_results.json |
| API retries in benchmark | 4/120 calls (all recovered) | paper/CRISPRArchitect_v2_manuscript.md |

---

## Appendix B: Key References

1. Komor, A. C. et al. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature* **533**, 420-424 (2016). [CBE]
2. Gaudelli, N. M. et al. Programmable base editing of A*T to G*C in genomic DNA without DNA cleavage. *Nature* **551**, 464-471 (2017). [ABE]
3. Anzalone, A. V. et al. Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature* **576**, 149-157 (2019). [Prime editing]
4. Paquet, D. et al. Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. *Nature* **533**, 125-129 (2016). [Cut-to-edit distance]
5. Ihry, R. J. et al. p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nat. Med.* **24**, 939-946 (2018). [p53 in iPSCs]
6. Kosicki, M. et al. Repair of double-strand breaks induced by CRISPR-Cas9 leads to large deletions and complex rearrangements. *Nat. Biotechnol.* **36**, 765-771 (2018). [DSB risks]
7. Iyer, S. et al. Efficient homology-directed repair with circular single-stranded DNA donors. *CRISPR J.* **5**, 685-701 (2022). [cssDNA]
8. Richards, S. et al. Standards and guidelines for the interpretation of sequence variants. *Genet. Med.* **17**, 405-424 (2015). [ACMG standards]
9. Arbab, M. et al. Determinants of base editing outcomes from target library analysis and machine learning. *Cell* **182**, 463-480 (2020). [BE outcomes, BE-Hive]
10. Elliott, B. et al. Gene conversion tracts from double-strand break repair in mammalian cells. *Mol. Cell. Biol.* **18**, 93-101 (1998). [Tract lengths]
11. Richter, M. F. et al. Phage-assisted evolution of an adenine base editor with improved Cas domain compatibility and activity. *Nat. Biotechnol.* **38**, 883-891 (2020). [ABE8e]
12. Chen, P. J. et al. Enhanced prime editing systems by manipulating cellular determinants of editing outcomes. *Cell* **184**, 5635-5652 (2021). [PE optimization]
13. Nelson, J. W. et al. Engineered pegRNAs improve prime editing efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022). [epegRNA]
14. Chauhan, V. P. et al. Altered DNA repair pathway engagement by engineered CRISPR-Cas9 nucleases. *PNAS* **120**, e2300605120 (2023). [Staggered cuts]
15. Leibowitz, M. L. et al. Chromothripsis as an on-target consequence of CRISPR-Cas9 genome editing. *Nat. Genet.* **53**, 895-905 (2021). [Translocation risk]
16. Richardson, C. D. et al. Enhancing homology-directed genome editing by catalytically active and inactive CRISPR-Cas9 using asymmetric donor DNA. *Nat. Biotechnol.* **34**, 339-344 (2016). [Asymmetric donors]
17. McLaren, W. et al. The Ensembl Variant Effect Predictor. *Genome Biol.* **17**, 122 (2016). [VEP]
18. Rees, H. A. & Liu, D. R. Base editing: precision chemistry on the genome and transcriptome of living cells. *Nat. Rev. Genet.* **19**, 770-788 (2018). [BE review]

---

*Document generated: March 2026*
*CRISPRArchitect version: v2*
*Total project size: ~35,000 LOC (v1 + v2), 43+ files*
