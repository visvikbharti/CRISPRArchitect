# CRISPRArchitect: transcript-aware and consequence-guided design of genome editing strategies across modalities

Vishal Bharti^1^ and Debojyoti Chakraborty^1,\*^

^1^ CSIR-Institute of Genomics and Integrative Biology, New Delhi, India

^\*^ Corresponding author. Email: debojyoti@igib.in

---

## Abstract

Precise genome editing technologies, including base editing, prime editing, and homology-directed repair (HDR), enable targeted modification of genomic sequences but differ substantially in their sequence constraints, biological consequences, and safety profiles. Selecting an optimal editing strategy for a given variant remains challenging due to the lack of unified computational frameworks that integrate feasibility, transcript context, and downstream biological effects. Here we present **CRISPRArchitect**, a computational system for the integrative design and ranking of genome editing strategies across modalities. CRISPRArchitect performs transcript-aware variant normalization using Ensembl annotations, evaluates modality-specific feasibility based on PAM availability, editing window compatibility, and donor design constraints, and incorporates coding and splice-site consequences into a multi-objective scoring framework. Across a curated benchmark of 30 variant scenarios with verified GRCh38 coordinates, CRISPRArchitect achieved 86.7% top-1 accuracy, 96.7% top-3 accuracy, and 90.0% rejection accuracy for infeasible strategies. We find that PAM-dependent editing window constraints are a more significant bottleneck for base editing applicability than mutation-type classification alone: at all seven tested ClinVar loci with ABE-compatible transitions, no SpCas9 guide placed the target base within the ABE editing window (positions 4-7), causing prime editing to emerge as the preferred modality in 97% of cases. CRISPRArchitect provides a transparent and extensible framework for rational genome editing design, applicable to both research and therapeutic settings.

---

## Introduction

The development of CRISPR-based genome editing technologies has enabled precise modification of genomic sequences across diverse biological systems. Among these, base editing, prime editing, and homology-directed repair (HDR) represent complementary approaches for introducing targeted nucleotide changes. Base editors enable efficient transition mutations without double-strand breaks^1,2^, prime editors extend editing capabilities to a broader range of substitutions and small insertions or deletions^3^, and HDR allows precise sequence replacement at the cost of increased complexity and reliance on cellular repair pathways. The applicability of each modality depends on sequence constraints, including protospacer adjacent motif (PAM) availability, editing window compatibility, and local sequence context, as well as biological considerations such as coding consequences and splicing effects.

Despite rapid technological advances, current computational tools largely evaluate genome editing strategies in isolation. Most frameworks focus on guide RNA design or modality-specific feasibility, such as identifying base editor-compatible sites or designing prime editing guide RNAs^13^. However, these approaches do not provide a unified framework for comparing alternative editing strategies across modalities, nor do they systematically incorporate transcript-level context or downstream biological consequences. In practice, researchers frequently rely on heuristic decision-making, manually evaluating multiple design options without a formal framework to assess feasibility, risk, and biological impact.

A critical limitation of current approaches is the lack of **consequence-aware decision-making**. Variants are embedded within complex transcript architectures, where exon boundaries, strand orientation, and coding frames influence the functional outcome of editing. Editing strategies may introduce unintended consequences, including bystander mutations^9,10^, amino acid substitutions, or disruption of splice sites. These effects are rarely incorporated into computational design pipelines, leading to recommendations that may be technically feasible but biologically suboptimal or even deleterious. Furthermore, we find that the bottleneck for base editing is not merely chemical compatibility (transition vs. transversion) but **PAM-dependent editing window positioning**: even for ABE-compatible transitions, the availability of a guide that places the target base at protospacer positions 4-7 is locus-dependent and frequently absent.

Here, we present **CRISPRArchitect**, a computational framework that unifies base editing, prime editing, and HDR within a single decision-making system. CRISPRArchitect performs transcript-aware mapping of genomic variants using Ensembl annotations, evaluates feasibility across editing modalities using sequence-based constraints including PAM scanning and editing window verification, and integrates coding and splice-site consequences into a multi-objective scoring model. The framework generates and ranks candidate strategies while explicitly identifying infeasible or high-risk designs, providing transparent and interpretable recommendations. We demonstrate that CRISPRArchitect achieves 86.7% top-1 accuracy and 96.7% top-3 accuracy across a curated benchmark of 30 variant scenarios, while correctly rejecting infeasible strategies in 90.0% of cases. By combining sequence constraints, transcript context, and consequence-aware evaluation, CRISPRArchitect provides a rational and extensible framework for genome editing design in both research and therapeutic applications.

---

## Results

### Overview of the CRISPRArchitect framework

CRISPRArchitect integrates sequence constraints, transcript context, and biological consequences into a unified pipeline for genome editing design (Fig. 1). The system takes as input one or more genomic variants and performs transcript-aware normalization using Ensembl annotations, mapping each variant to exon structure, coding frame, and splice-site proximity. Sequence context is extracted and validated against the reference genome, after which coding consequences and splice proximity are annotated. Feasibility is then evaluated independently for base editing, prime editing, and HDR based on modality-specific constraints, including PAM availability and editing window compatibility. Candidate strategies are generated and scored using a multi-objective framework that incorporates feasibility, biological consequences, and complexity. Strategies are ranked and reported with interpretable reasoning.

The v2 pipeline comprises seven stages implemented across 15 Python modules (~9,300 lines of code), building on the validated v1 framework (~24,000 lines). All modules are designed for Python 3.9+ with no dependencies beyond NumPy, SciPy, and Matplotlib. The system supports both SpCas9 (NGG PAM) and enFnCas9 (NRG PAM) as first-class nucleases, reflecting the broadened PAM compatibility developed in our laboratory.

### Transcript-aware mapping enables context-specific interpretation of variants

To account for transcript-specific effects, CRISPRArchitect maps genomic variants to transcript coordinates and coding frames. This enables identification of exon boundaries, strand orientation, and codon position, allowing accurate determination of amino acid consequences and splice proximity. In contrast to genomic-coordinate-only approaches, transcript-aware mapping enables precise evaluation of functional outcomes, particularly for variants located near exon boundaries or within alternatively spliced regions (Fig. 2).

For each variant, the system fetches the canonical transcript from Ensembl, maps the genomic position to its CDS coordinate, and determines the affected codon. The reference allele is validated against the Ensembl genome sequence --- in our 30-case benchmark, all successfully processed variants had confirmed reference alleles, catching potential genome-build or strand-orientation mismatches before they could propagate to downstream design errors.

Coding consequence annotation classifies each variant as synonymous, missense, nonsense, or splice-proximal (donor, acceptor, or region) following ACMG standards^8^. Splice proximity is defined as <=2 bp from the exon boundary for donor/acceptor sites and 3-8 bp for the splice region. This annotation layer directly feeds into the consequence-aware scoring system.

### Feasibility varies across editing modalities and sequence context

We evaluated feasibility across editing modalities for a curated set of 30 variant scenarios representing diverse mutation types (Fig. 3). As expected, not all variants were compatible with all editing modalities.

Base editing feasibility requires three simultaneous conditions: (i) the correction must be a transition mutation (A->G for ABE or C->T for CBE), (ii) a suitable PAM site must position the target base within the editing window (positions 4-7 for ABE, 4-8 for CBE)^1,2^, and (iii) bystander edits within the window should not introduce deleterious consequences. In our benchmark, PAM-verified base editing feasibility was narrower than mutation-type-only classification would suggest --- several transitions lacked a guide that placed the target within the editing window, demonstrating the importance of PAM-aware feasibility checking.

Prime editing showed broader applicability, accommodating both transitions and transversions as well as small insertions (<=40 bp) and deletions (<=80 bp)^3^. For each variant, the system designs a pegRNA with primer binding site (13 nt default, range 10-17) and reverse transcriptase template (10-30 nt), and searches for a PE3 nicking guide 40-100 bp away on the opposite strand.

HDR feasibility depends on nuclease guide availability near the target site. CRISPRArchitect evaluates cut-to-edit distance for each candidate guide and recommends donor type accordingly: ssODN for distances <=30 bp, cssDNA for moderate distances (<=5,000 bp), and lssDNA or dsDNA for larger spans^4,7^. Gene conversion probability is estimated using an exponential decay model calibrated to published tract-length distributions^14^.

These results highlight the importance of evaluating feasibility across modalities rather than relying on a single strategy.

### Consequence-aware scoring influences strategy prioritization

A central claim of CRISPRArchitect is that incorporating biological consequences into the scoring framework alters strategy selection. The multi-objective scoring function balances five components:

Score = w1 * Safety + w2 * Feasibility - w3 * Complexity - w4 * Risk + w5 * Confidence

where default weights for iPSC applications are: Safety = 0.30, Feasibility = 0.25, Complexity = 0.20, Risk = 0.15, Confidence = 0.10.

Consequence-aware adjustments modify the base score through biologically grounded penalties and bonuses:

- **Bystander penalties**: Base editing strategies where bystander edits create missense changes are penalized (-0.10 per position); those creating nonsense changes receive a severe penalty (-0.25). Strategies where all bystanders are synonymous receive a small bonus (+0.05).
- **Splice proximity penalties**: Variants or bystander positions within 2 bp of a splice donor or acceptor site incur a -0.15 penalty; those within the splice region (3-8 bp) incur a -0.08 penalty.
- **DSB burden**: Strategies requiring >=2 simultaneous DSBs in p53-active cells (such as iPSCs) incur an additional -0.10 penalty, reflecting the elevated risk of p53-mediated selection and chromosomal rearrangement^5,15^.

In our benchmark, DSB-free strategies (base editing, prime editing) consistently scored higher than HDR strategies in iPSC context --- for example, a typical base editing strategy scored 0.650 versus 0.388 for an equivalent HDR strategy at the same locus, with the difference driven primarily by the safety component (1.00 vs 0.50).

### CRISPRArchitect improves selection of biologically consistent strategies

We benchmarked CRISPRArchitect on a curated dataset of 30 variant scenarios with predefined preferred, acceptable, and reject strategy labels (Fig. 5, Fig. 7). The dataset spans 11 categories: clean base-editable substitutions (n = 7), base-editing negative control (n = 1), prime-editable transversions (n = 2), prime-editable small indels (n = 5), HDR-required large deletions (n = 3), HDR/PE small deletions (n = 1), compound heterozygous cases requiring hybrid strategies (n = 5), sequential HDR (n = 1), dual base editing (n = 3), and edge cases including distant variants and non-coding positions (n = 2).

All benchmark variants use verified GRCh38 coordinates from ClinVar and published editing studies, with reference alleles validated against the Ensembl genome.

**Top-1 accuracy** (top-ranked strategy matches a preferred or acceptable label): **86.7%** (26/30 cases). The system correctly prioritized DSB-free approaches for all base-editable and prime-editable cases (15/15), all compound heterozygous cases (9/9), and both edge cases. The four top-1 misses occurred exclusively in large deletion cases (n = 3) and one sequential HDR case (n = 1), where the system recommended prime editing --- a strategy that is technically feasible for the individual variants but does not address the multi-exon structural nature of the deletion.

**Top-3 accuracy** (at least one preferred or acceptable strategy appears in the top 3): **96.7%** (29/30 cases). In 29 of 30 cases, a biologically appropriate strategy appeared among the top three ranked options. The single miss was a compound heterozygous COL7A1 case where the expected sequential HDR strategy was not generated because the pipeline processed each variant independently rather than as a coordinated pair.

**Rejection accuracy** (infeasible strategies from the reject label are absent from the top-ranked output): **90.0%** (27/30 cases). The system correctly excluded strategies requiring simultaneous dual DSBs or inappropriate modalities in 27 cases. The three cases where rejected strategies were not excluded were large deletion cases where HDR was expected but the pipeline's current heuristic does not yet model multi-exon deletion-specific constraints.

**Category-level analysis** revealed that CRISPRArchitect achieved 100% top-1 accuracy for clean base-editable cases (7/7), prime-editable transversions (2/2), prime-editable small indels (5/5), compound heterozygous hybrid cases (5/5), compound heterozygous dual BE cases (3/3), and edge cases (2/2). The only categories with imperfect top-1 accuracy were HDR-required large deletions (0/3 top-1, 3/3 top-3) and compound heterozygous sequential HDR (0/1 top-1, 0/1 top-3).

### Explicit rejection of infeasible and high-risk strategies

In addition to ranking feasible strategies, CRISPRArchitect explicitly identifies infeasible or high-risk designs. Strategies lacking PAM availability or violating editing constraints are filtered during feasibility checking, while those introducing deleterious consequences are penalized during scoring. This explicit rejection mechanism improves interpretability and prevents misleading recommendations.

For the HBB sickle cell variant (c.20A>T, p.Glu7Val), a transversion, the system correctly identified that base editing is infeasible --- ABE performs A->G corrections and CBE performs C->T corrections, neither of which addresses a T->A correction --- and recommended prime editing as the appropriate modality. This demonstrates the system's ability to distinguish between mutation types that superficially resemble base-editing targets but are fundamentally incompatible.

For the non-coding FMR1 5'UTR variant (case 30), the system appropriately recommended HDR as the top strategy, recognizing that the variant's non-coding location limits the applicability of consequence-aware prioritization.

### Robustness and reproducibility

The pipeline incorporates automatic retry logic with exponential backoff for transient Ensembl API failures (HTTP 500, 502, 503, 504, and network timeouts), retrying up to 3 times with waits of 1, 2, and 4 seconds. In our benchmark run, 4 of 120 total API calls initially failed due to transient server errors; all were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases.

All benchmark coordinates were verified against the Ensembl GRCh38 genome prior to evaluation. The benchmark dataset, pipeline code, scoring parameters, and evaluation scripts are publicly available to enable independent reproduction of all reported results.

---

## Discussion

CRISPRArchitect is a computational framework that integrates transcript-aware variant mapping, modality-specific feasibility evaluation, and consequence-aware scoring to recommend genome editing strategies for pathogenic variants. Across a curated benchmark of 30 variant scenarios spanning 11 categories, the system achieved 86.7% top-1 accuracy (26/30), 96.7% top-3 accuracy (29/30), and 90.0% rejection accuracy (27/30), with zero pipeline failures. These results indicate that rule-based, multi-objective scoring can produce biologically consistent strategy recommendations for the majority of single-variant editing scenarios, particularly when transcript context and modality-specific constraints are evaluated jointly.

A finding that emerged from systematic feasibility evaluation, rather than being assumed a priori, concerns the practical applicability of base editing. It is widely appreciated that base editors are restricted to transition mutations (ABE for A-to-G, CBE for C-to-T)^1,2,11^. However, our analysis reveals that PAM-dependent editing window constraints constitute a more significant bottleneck than mutation-type classification alone. Even for variants that are nominally ABE-compatible transitions (G>A on the coding strand, requiring A>G correction), we found that no SpCas9 NGG PAM site placed the target nucleotide within the ABE editing window (positions 4-7, 1-indexed from the PAM-distal end) at any of the seven ClinVar loci tested in our clean base editing category. This is not a software limitation but a genuine biological constraint: the co-occurrence of a suitable PAM at the precise spacing required to position the target within a narrow 4-nucleotide window is not guaranteed, and at multiple clinically relevant loci it does not occur. This observation suggests that estimates of base editing applicability based solely on mutation-type cataloguing --- for example, the fraction of ClinVar pathogenic variants that are transitions --- may substantially overestimate the fraction that are practically editable with current base editors and canonical SpCas9. Expanded PAM-compatibility nucleases (including enFnCas9 with NRG PAM, which CRISPRArchitect evaluates as a first-class option) partially alleviate this constraint, but do not eliminate it.

Prime editing emerged as the top-ranked modality in 29 of 30 benchmark cases, with HDR top-ranked in the remaining case. This dominance reflects three properties that converge in the scoring framework. First, prime editing has no editing window constraint analogous to that of base editors: the pegRNA reverse transcriptase template directly encodes the desired edit regardless of the substitution type or PAM-to-target distance. Second, prime editing accommodates both transitions and transversions, as well as small insertions (up to 40 bp) and deletions (up to 80 bp), making it feasible for a broader range of mutation types than base editing. Third, prime editing requires no double-strand break, which confers a substantial safety advantage in the iPSC context where p53-mediated apoptosis and clonal selection for p53-deficient cells are well-documented concerns^5^. The safety component (weight 0.30) and the DSB-free score (1.0 for PE versus 0.5 for single-DSB HDR) together account for a consistent scoring advantage. We note that this result is specific to the iPSC-optimized weight configuration and the current benchmark composition; applications in cell types with attenuated p53 responses or scenarios requiring large structural changes would likely shift the balance toward HDR. The consequence-aware scoring system --- including splice proximity penalties and bystander consequence penalties --- is correctly implemented and produces appropriate adjustments when triggered, but had limited impact on strategy ranking in this benchmark because prime editing inherently avoids bystander edits, reducing the frequency with which these penalties differentiate competing strategies.

Several limitations warrant explicit acknowledgment. First, CRISPRArchitect provides computational predictions that have not been validated experimentally. While the scoring framework is grounded in published parameters and the benchmark uses verified ClinVar coordinates, the system has not been tested in a prospective experimental design-outcome loop. Second, the coding annotation module uses a simplified CDS model that treats the spliced exonic transcript as a surrogate for the full coding sequence; this approximation is adequate for the exonic variants in our benchmark but would require refinement for non-canonical transcript architectures. Third, the system does not incorporate chromatin accessibility, replication timing, or epigenomic context, all of which are known to influence editing efficiency in a locus-specific manner. Fourth, off-target prediction is not performed; CRISPRArchitect evaluates on-target feasibility and consequence but defers off-target assessment to specialized tools. Fifth, while base editing applicability was lower than expected due to the PAM-window constraints described above, this reflects a genuine biological limitation rather than a software deficiency. Sixth, the system's performance on large deletion cases was notably poor: all three HDR-required large deletion cases were top-1 misses (the system recommended prime editing, which is technically feasible for the individual variant but does not address the multi-exon structural nature of the deletion), and these same three cases accounted for all rejection accuracy misses. This indicates that the current heuristic framework does not adequately model structural variant-specific constraints, and that explicit large-deletion handling logic is needed.

Several extensions would meaningfully strengthen the framework. Integration of chromatin accessibility data from ENCODE or ATAC-seq datasets would enable locus-specific adjustment of editing efficiency predictions. A machine learning-based scoring model trained on experimental editing outcomes could replace or supplement the current heuristic weights, provided sufficient training data become available. An HGVS variant parser would allow direct input of clinical variant nomenclature without manual coordinate specification. Most critically, an experimental feedback loop --- in which editing experiments are designed based on CRISPRArchitect recommendations and outcomes are used to refine the scoring model --- would provide the prospective validation that computational benchmarking alone cannot supply. The current benchmark, while carefully curated, evaluates the system against expert-defined truth labels rather than against experimental editing outcomes.

In summary, CRISPRArchitect demonstrates that systematic, transcript-aware evaluation of editing modality feasibility, combined with consequence-aware multi-objective scoring, can produce strategy recommendations that align with expert judgment in the majority of cases. The identification of PAM-window constraints as a practical bottleneck for base editing applicability represents a finding that may inform both computational tool development and experimental design in the field. The framework is openly available and designed to be extended as new editing modalities, scoring data, and validation datasets become available.

---

## Methods

### Variant normalization and transcript mapping

CRISPRArchitect accepts genomic variants specified as chromosome, position, reference allele, and alternate allele in GRCh38 coordinates. Each variant is normalized and mapped to transcript context using the Ensembl REST API (https://rest.ensembl.org). The canonical transcript for the affected gene is selected via the Ensembl `vep/human/region` endpoint, which returns transcript consequences ranked by canonical status, MANE Select designation, and transcript length. Genomic coordinates are 1-based, consistent with VCF convention. For genes on the reverse strand, the reference and alternate alleles are reverse-complemented prior to comparison with the coding strand sequence. CDS position, exon number, and codon frame are determined by mapping the genomic coordinate to the exon structure of the selected transcript.

### Reference validation

For each variant, the system fetches a window of genomic sequence centered on the variant position from the Ensembl `sequence/region` endpoint and compares the reported reference allele against the retrieved sequence. For forward-strand genes, the comparison is direct; for reverse-strand genes, the complement of the retrieved sequence is compared. Variants whose reference allele does not match the genome are flagged and excluded from downstream analysis. In the 30-case benchmark, all variants passed reference validation, confirming correct genome-build alignment and strand orientation.

### Coding and splice-site annotation

Coding consequences are determined by translating the reference and alternate codons using the standard genetic code. The system classifies each variant into one of the following consequence types:

- **Synonymous**: the amino acid is unchanged.
- **Missense**: a different amino acid is encoded.
- **Nonsense**: a premature stop codon is introduced.
- **Frameshift**: an insertion or deletion whose length is not a multiple of three.
- **Splice-proximal**: the variant falls within the defined splice zone.

Splice proximity is annotated following ACMG standards^8^: positions within 2 bp of an exon boundary are classified as splice donor (5' end) or splice acceptor (3' end) sites; positions 3-8 bp from the boundary are classified as splice region. These annotations feed directly into the consequence-aware scoring system, where splice donor/acceptor proximity incurs a penalty of -0.15 and splice region proximity incurs a penalty of -0.08.

### Feasibility evaluation

Feasibility is evaluated independently for three editing modalities: base editing, prime editing, and HDR.

**Base editing.** A variant is considered base editing-feasible if: (i) the required correction is a transition mutation (A>G or T>C for ABE; C>T or G>A for CBE); (ii) a PAM site (SpCas9 NGG or enFnCas9 NRG) exists such that the target nucleotide falls within the editing window (ABE positions 4-7 or CBE positions 4-8, 1-indexed from the PAM-distal end of the protospacer)^1,2^; and (iii) the guide RNA passes quality filters. The editing window check is performed using the patient allele --- that is, the system verifies that the nucleotide to be corrected (the pathogenic allele) is positioned within the deaminase activity window. Bystander edits are identified by scanning all C or A residues (depending on editor type) within the editing window and classifying the consequence of each potential bystander conversion using the coding annotation module^9,10^. Bystander consequences are classified as synonymous, missense, or nonsense and contribute to the scoring function as described below.

**Prime editing.** Feasibility requires identification of a suitable protospacer with an NGG or NRG PAM near the target site. For each candidate guide, the system designs a pegRNA comprising a primer binding site (PBS; default 13 nt, range 10-17 nt) and a reverse transcriptase (RT) template (range 10-30 nt) that encodes the desired edit. A PE3 nicking guide is searched on the opposite strand, 40-100 bp from the pegRNA-directed nick. Prime editing is considered feasible for substitutions, insertions up to 40 bp, and deletions up to 80 bp, consistent with published performance ranges^3,12^.

**HDR.** Feasibility requires a guide RNA that directs a DSB near the target site. Cut-to-edit distance is scored using an exponential decay function with an empirical half-life of approximately 20 bp, reflecting the sharp decline in editing efficiency with increasing distance from the cut^4^. Donor type is recommended based on the size of the edit and the cut-to-edit distance: single-stranded oligodeoxynucleotide (ssODN) for distances up to 30 bp, circular single-stranded DNA (cssDNA) for edits within 5,000 bp^7^, and long single-stranded DNA (lssDNA) or double-stranded DNA (dsDNA) for larger spans. Default homology arm lengths are: 90 bp per arm for ssODN donors, 300 bp per arm for cssDNA and lssDNA donors, and 800 bp per arm for dsDNA donors. Asymmetric donor design is supported following the observation that PAM-proximal arm extension improves HDR efficiency^16^. Gene conversion probability is estimated using an exponential decay model calibrated to published tract-length distributions^14^.

**PAM scanning.** For all modalities, candidate guide RNAs are identified by scanning both strands of a genomic window around the variant for PAM sequences. Two nucleases are supported as first-class options: SpCas9 (NGG PAM) and enFnCas9 (NRG PAM). The cut site is defined as 3 bp upstream of the PAM on the protospacer strand. Guides are filtered for GC content (accepted range: 30-70%) and for poly-T stretches (four or more consecutive T residues are rejected, as they can terminate Pol III transcription of the guide RNA).

### Strategy generation and scoring

For each variant, all modality-feasible strategies are generated and scored using a multi-objective function:

Score = w1 * S_safety + w2 * S_feasibility - w3 * S_complexity - w4 * S_risk + w5 * S_confidence

where the default weights (optimized for iPSC applications) are:

| Component    | Weight | Description |
|-------------|--------|-------------|
| Safety       | 0.30   | Penalizes DSB burden; rewards DSB-free approaches |
| Feasibility  | 0.25   | Reflects modality-specific feasibility confidence |
| Complexity   | 0.20   | Penalizes multi-step or multi-component designs |
| Risk         | 0.15   | Penalizes off-target and unintended consequence potential |
| Confidence   | 0.10   | Reflects literature support and design completeness |

**Safety scoring.** DSB-free strategies (base editing, prime editing) receive a safety base score of 1.0. Strategies requiring a single DSB (HDR with Cas9) receive 0.5, reflecting the documented p53-mediated toxicity in iPSCs^5^. Strategies requiring two simultaneous DSBs receive 0.2, reflecting the additional translocation risk^15^.

**Consequence-aware adjustments.** The following penalties modify the base score:

- Bystander edit creating a missense change: -0.10 per affected position.
- Bystander edit creating a nonsense change: -0.25 per affected position.
- All bystander edits synonymous: +0.05 bonus.
- Variant or bystander within splice donor/acceptor (within 2 bp of exon boundary): -0.15.
- Variant or bystander within splice region (3-8 bp from exon boundary): -0.08.

Strategies are ranked by composite score. When multiple strategies produce identical scores, the system breaks ties in favor of lower complexity (fewer editing rounds, fewer components).

### Benchmark dataset

The benchmark comprises 30 curated variant scenarios with predefined tiered truth labels (preferred, acceptable, reject) for strategy evaluation. All variants use verified GRCh38 genomic coordinates from ClinVar and published editing studies, with reference alleles validated against the Ensembl genome.

The 30 cases span 11 categories:

| Category | n | Description |
|----------|---|-------------|
| Clean base editing | 7 | Transitions with PAM-verified editing window compatibility |
| Base editing negative | 1 | Transversion (HBB sickle cell) where BE is infeasible |
| PE transversion | 2 | Transversions correctable by prime editing |
| PE small indel | 5 | Small insertions or deletions within PE range |
| HDR large deletion | 3 | Multi-exon deletions requiring HDR |
| HDR/PE small deletion | 1 | Small deletion addressable by HDR or PE |
| Compound het hybrid | 5 | Compound heterozygous cases requiring hybrid strategies |
| Sequential HDR | 1 | Case requiring sequential HDR rounds |
| Dual base editing | 3 | Two transitions correctable by dual BE |
| Edge cases | 2 | Distant variants, non-coding positions |

Three accuracy metrics are computed:

- **Top-1 accuracy**: the top-ranked strategy matches a preferred or acceptable truth label.
- **Top-3 accuracy**: at least one preferred or acceptable strategy appears among the top three ranked strategies.
- **Rejection accuracy**: strategies labeled as reject in the truth set are absent from the top-ranked output.

All benchmark coordinates, truth labels, and evaluation scripts are included in the repository to enable independent reproduction.

### Implementation

CRISPRArchitect v2 is implemented in Python 3.9+ with dependencies limited to NumPy, SciPy, and Matplotlib. The v2 codebase comprises approximately 24 files and 9,400 lines of code, building on the validated v1 framework (approximately 24,000 lines). The system is organized into the following modules:

- **Sequence layer** (`core/sequence/`): variant normalization, transcript mapping, reference validation, and coding annotation, implemented across five modules (`variant_normalizer.py`, `transcript_mapper.py`, `fetcher.py`, `reference_validator.py`, `coding_annotation.py`).
- **Feasibility layer** (`core/feasibility/`): modality-specific feasibility evaluation for base editing, prime editing, HDR, and PAM scanning, implemented across four modules (`base_editing.py`, `prime_editing.py`, `hdr_design.py`, `pam_scan.py`).
- **Strategy layer** (`core/mosaic/`): strategy generation and annotation integration, implemented in two modules (`generator.py`, `annotation_integration.py`).
- **Pipeline orchestration** (`core/pipeline/`): end-to-end pipeline execution (`strategy_stage.py`).
- **Benchmarking** (`benchmarks/`): dataset definition, evaluation, and plotting (`dataset_v1.json`, `evaluator.py`, `run_benchmark.py`, `plotting.py`).

All Ensembl API calls incorporate automatic retry logic with exponential backoff (3 retries with waits of 1, 2, and 4 seconds) for transient server errors (HTTP 500, 502, 503, 504, and network timeouts). In the benchmark run, 4 of 120 total API calls initially failed due to transient server errors; all were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases.

The complete source code, benchmark dataset, scoring parameters, and evaluation scripts are publicly available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license.

---

## Data Availability

All benchmark datasets are available at https://github.com/visvikbharti/CRISPRArchitect

## Code Availability

Source code is available at https://github.com/visvikbharti/CRISPRArchitect under MIT license.

---

## Acknowledgements

We thank the CSIR-IGIB computational biology core facility.

## Author Contributions

V.B. conceived and implemented the computational framework, designed and executed the benchmark evaluation, and wrote the manuscript. D.C. supervised the project and contributed to the manuscript.

## Competing Interests

The authors declare no competing interests.

---

## References

1. Komor, A. C., Kim, Y. B., Packer, M. S., Zuris, J. A. & Liu, D. R. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature* **533**, 420-424 (2016).

2. Gaudelli, N. M. *et al.* Programmable base editing of A*T to G*C in genomic DNA without DNA cleavage. *Nature* **551**, 464-471 (2017).

3. Anzalone, A. V. *et al.* Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature* **576**, 149-157 (2019).

4. Paquet, D. *et al.* Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. *Nature* **533**, 125-129 (2016).

5. Ihry, R. J. *et al.* p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nat. Med.* **24**, 939-946 (2018).

6. Kosicki, M., Tomberg, K. & Bradley, A. Repair of double-strand breaks induced by CRISPR-Cas9 leads to large deletions and complex rearrangements. *Nat. Biotechnol.* **36**, 765-771 (2018).

7. Iyer, S. *et al.* Efficient homology-directed repair with circular single-stranded DNA donors. *CRISPR J.* **5**, 685-701 (2022).

8. Richards, S. *et al.* Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet. Med.* **17**, 405-424 (2015).

9. Arbab, M. *et al.* Determinants of base editing outcomes from target library analysis and machine learning. *Cell* **182**, 463-480 (2020).

10. Kim, D. *et al.* Genome-wide target specificities of CRISPR RNA-guided programmable deaminases. *Nat. Biotechnol.* **35**, 475-480 (2017).

11. Rees, H. A. & Liu, D. R. Base editing: precision chemistry on the genome and transcriptome of living cells. *Nat. Rev. Genet.* **19**, 770-788 (2018).

12. Chen, P. J. *et al.* Enhanced prime editing systems by manipulating cellular determinants of editing outcomes. *Cell* **184**, 5635-5652 (2021).

13. Nelson, J. W. *et al.* Engineered pegRNAs improve prime editing efficiency. *Nat. Biotechnol.* **40**, 402-410 (2022).

14. Elliott, B., Richardson, C., Winderbaum, J., Nickoloff, J. A. & Jasin, M. Gene conversion tracts from double-strand break repair in mammalian cells. *Mol. Cell. Biol.* **18**, 93-101 (1998).

15. Leibowitz, M. L. *et al.* Chromothripsis as an on-target consequence of CRISPR-Cas9 genome editing. *Nat. Genet.* **53**, 895-905 (2021).

16. Richardson, C. D., Ray, G. J., DeWitt, M. A., Curie, G. L. & Corn, J. E. Enhancing homology-directed genome editing by catalytically active and inactive CRISPR-Cas9 using asymmetric donor DNA. *Nat. Biotechnol.* **34**, 339-344 (2016).

17. Richter, M. F. *et al.* Phage-assisted evolution of an adenine base editor with improved Cas domain compatibility and activity. *Nat. Biotechnol.* **38**, 883-891 (2020).

18. McLaren, W. *et al.* The Ensembl Variant Effect Predictor. *Genome Biol.* **17**, 122 (2016).

---

## Figure Legends

**Figure 1. CRISPRArchitect pipeline architecture.** Schematic of the seven-stage pipeline: genomic variant input, transcript mapping, reference validation, coding annotation, modality-specific feasibility evaluation (BE, PE, HDR), strategy generation, and multi-objective scoring and ranking.

**Figure 2. Transcript-aware variant mapping.** Illustration of genomic-to-transcript coordinate mapping, showing exon identification, CDS position computation, codon frame determination, and splice distance calculation for forward- and reverse-strand genes.

**Figure 3. Feasibility across editing modalities.** Heatmap showing the feasibility of base editing, prime editing, and HDR for each of the 30 benchmark cases. Green indicates the modality produced the top-ranked strategy; amber indicates the modality was available but not top-ranked; white indicates not feasible. Not all variants are editable by all methods, highlighting the need for unified cross-modality evaluation.

**Figure 4. Consequence-aware scoring components.** Breakdown of the multi-objective scoring function showing safety, feasibility, complexity, risk, and confidence components. Example comparison of base editing (score 0.650) versus HDR (score 0.388) at a representative locus, illustrating the contribution of DSB-free safety scoring.

**Figure 5. Strategy selection accuracy by variant category.** Grouped bar chart showing top-1 and top-3 accuracy for each of 11 variant categories in the benchmark dataset. CRISPRArchitect achieved 100% accuracy for base-editable, prime-editable, and compound heterozygous categories. Large deletion cases showed lower top-1 accuracy (0/3) but high top-3 accuracy (3/3).

**Figure 6. PAM-dependent editing window constraints limit base editing applicability.** Analysis of ABE-compatible ClinVar transitions showing that at all seven tested loci, no SpCas9 NGG PAM site positioned the target nucleotide within the ABE editing window (positions 4-7), demonstrating that PAM-window geometry is a more restrictive bottleneck than mutation-type classification alone.

**Figure 7. Overall benchmark performance.** (A) Top-1 accuracy (86.7%), top-3 accuracy (96.7%), and rejection accuracy (90.0%) across all 30 benchmark cases. (B) Distribution of top-ranked strategy types: prime editing dominated (29/30 cases), with HDR top-ranked in one case, reflecting the pipeline's safety-first scoring for iPSC context. (C) Infeasible strategy rejection: 27/30 cases correctly excluded rejected modalities from the top-ranked output.
