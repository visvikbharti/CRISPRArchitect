# CRISPRArchitect v2: Discussion and Methods

**For:** Nature Methods manuscript
**Authors:** Vishal Bharti and Debojyoti Chakraborty
**Date:** March 2026

---

## Discussion

CRISPRArchitect is a computational framework that integrates transcript-aware variant mapping, modality-specific feasibility evaluation, and consequence-aware scoring to recommend genome editing strategies for pathogenic variants. Across a curated benchmark of 30 variant scenarios spanning 11 categories, the system achieved 86.7% top-1 accuracy (26/30), 96.7% top-3 accuracy (29/30), and 90.0% rejection accuracy (27/30), with zero pipeline failures. These results indicate that rule-based, multi-objective scoring can produce biologically consistent strategy recommendations for the majority of single-variant editing scenarios, particularly when transcript context and modality-specific constraints are evaluated jointly.

A finding that emerged from systematic feasibility evaluation, rather than being assumed a priori, concerns the practical applicability of base editing. It is widely appreciated that base editors are restricted to transition mutations (ABE for A-to-G, CBE for C-to-T). However, our analysis reveals that PAM-dependent editing window constraints constitute a more significant bottleneck than mutation-type classification alone. Even for variants that are nominally ABE-compatible transitions (G>A on the coding strand, requiring A>G correction), we found that no SpCas9 NGG PAM site placed the target nucleotide within the ABE editing window (positions 4--7, 1-indexed from the PAM-distal end) at any of the seven ClinVar loci tested in our clean base editing category. This is not a software limitation but a genuine biological constraint: the co-occurrence of a suitable PAM at the precise spacing required to position the target within a narrow 4-nucleotide window is not guaranteed, and at multiple clinically relevant loci it does not occur. This observation suggests that estimates of base editing applicability based solely on mutation-type cataloguing --- for example, the fraction of ClinVar pathogenic variants that are transitions --- may substantially overestimate the fraction that are practically editable with current base editors and canonical SpCas9. Expanded PAM-compatibility nucleases (including enFnCas9 with NRG PAM, which CRISPRArchitect evaluates as a first-class option) partially alleviate this constraint, but do not eliminate it.

Prime editing emerged as the top-ranked modality in 29 of 30 benchmark cases, with HDR top-ranked in the remaining case. This dominance reflects three properties that converge in the scoring framework. First, prime editing has no editing window constraint analogous to that of base editors: the pegRNA reverse transcriptase template directly encodes the desired edit regardless of the substitution type or PAM-to-target distance. Second, prime editing accommodates both transitions and transversions, as well as small insertions (up to 40 bp) and deletions (up to 80 bp), making it feasible for a broader range of mutation types than base editing. Third, prime editing requires no double-strand break, which confers a substantial safety advantage in the iPSC context where p53-mediated apoptosis and clonal selection for p53-deficient cells are well-documented concerns (Ihry et al., 2018; Haapaniemi et al., 2018). The safety component (weight 0.30) and the DSB-free score (1.0 for PE versus 0.5 for single-DSB HDR) together account for a consistent scoring advantage. We note that this result is specific to the iPSC-optimized weight configuration and the current benchmark composition; applications in cell types with attenuated p53 responses or scenarios requiring large structural changes would likely shift the balance toward HDR. The consequence-aware scoring system --- including splice proximity penalties and bystander consequence penalties --- is correctly implemented and produces appropriate adjustments when triggered, but had limited impact on strategy ranking in this benchmark because prime editing inherently avoids bystander edits, reducing the frequency with which these penalties differentiate competing strategies.

Several limitations warrant explicit acknowledgment. First, CRISPRArchitect provides computational predictions that have not been validated experimentally. While the scoring framework is grounded in published parameters and the benchmark uses verified ClinVar coordinates, the system has not been tested in a prospective experimental design-outcome loop. Second, the coding annotation module uses a simplified CDS model that treats the spliced exonic transcript as a surrogate for the full coding sequence; this approximation is adequate for the exonic variants in our benchmark but would require refinement for non-canonical transcript architectures. Third, the system does not incorporate chromatin accessibility, replication timing, or epigenomic context, all of which are known to influence editing efficiency in a locus-specific manner (Aymard et al., 2014). Fourth, off-target prediction is not performed; CRISPRArchitect evaluates on-target feasibility and consequence but defers off-target assessment to specialized tools such as CRISPOR (Concordet and Haeussler, 2018) or Cas-OFFinder (Bae et al., 2014). Fifth, while base editing applicability was lower than expected due to the PAM-window constraints described above, this reflects a genuine biological limitation rather than a software deficiency. Sixth, the system's performance on large deletion cases was notably poor: all three HDR-required large deletion cases were top-1 misses (the system recommended prime editing, which is technically feasible for the individual variant but does not address the multi-exon structural nature of the deletion), and these same three cases accounted for all rejection accuracy misses. This indicates that the current heuristic framework does not adequately model structural variant-specific constraints, and that explicit large-deletion handling logic is needed.

Several extensions would meaningfully strengthen the framework. Integration of chromatin accessibility data from ENCODE or ATAC-seq datasets would enable locus-specific adjustment of editing efficiency predictions. A machine learning-based scoring model trained on experimental editing outcomes could replace or supplement the current heuristic weights, provided sufficient training data become available. An HGVS variant parser would allow direct input of clinical variant nomenclature without manual coordinate specification. Most critically, an experimental feedback loop --- in which editing experiments are designed based on CRISPRArchitect recommendations and outcomes are used to refine the scoring model --- would provide the prospective validation that computational benchmarking alone cannot supply. The current benchmark, while carefully curated, evaluates the system against expert-defined truth labels rather than against experimental editing outcomes.

In summary, CRISPRArchitect demonstrates that systematic, transcript-aware evaluation of editing modality feasibility, combined with consequence-aware multi-objective scoring, can produce strategy recommendations that align with expert judgment in the majority of cases. The identification of PAM-window constraints as a practical bottleneck for base editing applicability represents a finding that may inform both computational tool development and experimental design in the field. The framework is openly available and designed to be extended as new editing modalities, scoring data, and validation datasets become available.

---

## Methods

### 1. Variant normalization and transcript mapping

CRISPRArchitect accepts genomic variants specified as chromosome, position, reference allele, and alternate allele in GRCh38 coordinates. Each variant is normalized and mapped to transcript context using the Ensembl REST API (https://rest.ensembl.org). The canonical transcript for the affected gene is selected via the Ensembl `vep/human/region` endpoint, which returns transcript consequences ranked by canonical status, MANE Select designation, and transcript length. Genomic coordinates are 1-based, consistent with VCF convention. For genes on the reverse strand, the reference and alternate alleles are reverse-complemented prior to comparison with the coding strand sequence. CDS position, exon number, and codon frame are determined by mapping the genomic coordinate to the exon structure of the selected transcript.

### 2. Reference validation

For each variant, the system fetches a window of genomic sequence centered on the variant position from the Ensembl `sequence/region` endpoint and compares the reported reference allele against the retrieved sequence. For forward-strand genes, the comparison is direct; for reverse-strand genes, the complement of the retrieved sequence is compared. Variants whose reference allele does not match the genome are flagged and excluded from downstream analysis. In the 30-case benchmark, all variants passed reference validation, confirming correct genome-build alignment and strand orientation.

### 3. Coding and splice-site annotation

Coding consequences are determined by translating the reference and alternate codons using the standard genetic code. The system classifies each variant into one of the following consequence types:

- **Synonymous**: the amino acid is unchanged.
- **Missense**: a different amino acid is encoded.
- **Nonsense**: a premature stop codon is introduced.
- **Frameshift**: an insertion or deletion whose length is not a multiple of three.
- **Splice-proximal**: the variant falls within the defined splice zone.

Splice proximity is annotated following ACMG standards (Richards et al., 2015): positions within 2 bp of an exon boundary are classified as splice donor (5' end) or splice acceptor (3' end) sites; positions 3--8 bp from the boundary are classified as splice region. These annotations feed directly into the consequence-aware scoring system, where splice donor/acceptor proximity incurs a penalty of -0.15 and splice region proximity incurs a penalty of -0.08.

### 4. Feasibility evaluation

Feasibility is evaluated independently for three editing modalities: base editing, prime editing, and HDR.

**Base editing.** A variant is considered base editing-feasible if: (i) the required correction is a transition mutation (A>G or T>C for ABE; C>T or G>A for CBE); (ii) a PAM site (SpCas9 NGG or enFnCas9 NRG) exists such that the target nucleotide falls within the editing window (ABE positions 4--7 or CBE positions 4--8, 1-indexed from the PAM-distal end of the protospacer); and (iii) the guide RNA passes quality filters. The editing window check is performed using the patient allele --- that is, the system verifies that the nucleotide to be corrected (the pathogenic allele) is positioned within the deaminase activity window. Bystander edits are identified by scanning all C or A residues (depending on editor type) within the editing window and classifying the consequence of each potential bystander conversion using the coding annotation module. Bystander consequences are classified as synonymous, missense, or nonsense and contribute to the scoring function as described below.

**Prime editing.** Feasibility requires identification of a suitable protospacer with an NGG or NRG PAM near the target site. For each candidate guide, the system designs a pegRNA comprising a primer binding site (PBS; default 13 nt, range 10--17 nt) and a reverse transcriptase (RT) template (range 10--30 nt) that encodes the desired edit. A PE3 nicking guide is searched on the opposite strand, 40--100 bp from the pegRNA-directed nick. Prime editing is considered feasible for substitutions, insertions up to 40 bp, and deletions up to 80 bp, consistent with published performance ranges (Anzalone et al., 2019; Chen et al., 2021).

**HDR.** Feasibility requires a guide RNA that directs a DSB near the target site. Cut-to-edit distance is scored using an exponential decay function with an empirical half-life of approximately 20 bp, reflecting the sharp decline in editing efficiency with increasing distance from the cut (Paquet et al., 2016). Donor type is recommended based on the size of the edit and the cut-to-edit distance: single-stranded oligodeoxynucleotide (ssODN) for distances up to 30 bp, circular single-stranded DNA (cssDNA) for edits within 5,000 bp (Iyer et al., 2022), and long single-stranded DNA (lssDNA) or double-stranded DNA (dsDNA) for larger spans. Default homology arm lengths are: 90 bp per arm for ssODN donors, 300 bp per arm for cssDNA and lssDNA donors, and 800 bp per arm for dsDNA donors.

**PAM scanning.** For all modalities, candidate guide RNAs are identified by scanning both strands of a genomic window around the variant for PAM sequences. Two nucleases are supported as first-class options: SpCas9 (NGG PAM) and enFnCas9 (NRG PAM; Kulcsár et al., 2017). The cut site is defined as 3 bp upstream of the PAM on the protospacer strand. Guides are filtered for GC content (accepted range: 30--70%) and for poly-T stretches (four or more consecutive T residues are rejected, as they can terminate Pol III transcription of the guide RNA).

### 5. Strategy generation and scoring

For each variant, all modality-feasible strategies are generated and scored using a multi-objective function:

$$\text{Score} = w_1 \cdot S_{\text{safety}} + w_2 \cdot S_{\text{feasibility}} - w_3 \cdot S_{\text{complexity}} - w_4 \cdot S_{\text{risk}} + w_5 \cdot S_{\text{confidence}}$$

where the default weights (optimized for iPSC applications) are:

| Component    | Weight ($w$) | Description |
|-------------|-------------|-------------|
| Safety       | 0.30        | Penalizes DSB burden; rewards DSB-free approaches |
| Feasibility  | 0.25        | Reflects modality-specific feasibility confidence |
| Complexity   | 0.20        | Penalizes multi-step or multi-component designs |
| Risk         | 0.15        | Penalizes off-target and unintended consequence potential |
| Confidence   | 0.10        | Reflects literature support and design completeness |

**Safety scoring.** DSB-free strategies (base editing, prime editing) receive a safety base score of 1.0. Strategies requiring a single DSB (HDR with Cas9) receive 0.5, reflecting the documented p53-mediated toxicity in iPSCs (Ihry et al., 2018). Strategies requiring two simultaneous DSBs receive 0.2, reflecting the additional translocation risk (Leibowitz et al., 2021).

**Consequence-aware adjustments.** The following penalties modify the base score:

- Bystander edit creating a missense change: -0.10 per affected position.
- Bystander edit creating a nonsense change: -0.25 per affected position.
- All bystander edits synonymous: +0.05 bonus.
- Variant or bystander within splice donor/acceptor (within 2 bp of exon boundary): -0.15.
- Variant or bystander within splice region (3--8 bp from exon boundary): -0.08.

Strategies are ranked by composite score. When multiple strategies produce identical scores, the system breaks ties in favor of lower complexity (fewer editing rounds, fewer components).

### 6. Benchmark dataset

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

### 7. Implementation

CRISPRArchitect v2 is implemented in Python 3.9+ with dependencies limited to NumPy, SciPy, and Matplotlib. The v2 codebase comprises approximately 24 files and 9,400 lines of code, building on the validated v1 framework (approximately 24,000 lines). The system is organized into the following modules:

- **Sequence layer** (`core/sequence/`): variant normalization, transcript mapping, reference validation, and coding annotation, implemented across five modules (`variant_normalizer.py`, `transcript_mapper.py`, `fetcher.py`, `reference_validator.py`, `coding_annotation.py`).
- **Feasibility layer** (`core/feasibility/`): modality-specific feasibility evaluation for base editing, prime editing, HDR, and PAM scanning, implemented across four modules (`base_editing.py`, `prime_editing.py`, `hdr_design.py`, `pam_scan.py`).
- **Strategy layer** (`core/mosaic/`): strategy generation and annotation integration, implemented in two modules (`generator.py`, `annotation_integration.py`).
- **Pipeline orchestration** (`core/pipeline/`): end-to-end pipeline execution (`strategy_stage.py`).
- **Benchmarking** (`benchmarks/`): dataset definition, evaluation, and plotting (`dataset_v1.json`, `evaluator.py`, `run_benchmark.py`, `plotting.py`).

All Ensembl API calls incorporate automatic retry logic with exponential backoff (3 retries with waits of 1, 2, and 4 seconds) for transient server errors (HTTP 500, 502, 503, 504, and network timeouts). In the benchmark run, 4 of 120 total API calls initially failed due to transient server errors; all were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases.

The complete source code, benchmark dataset, scoring parameters, and evaluation scripts are publicly available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license.

---

## References

Anzalone AV, Randolph PB, Davis JR, Sousa AA, Koblan LW, Levy JM, et al. Search-and-replace genome editing without double-strand breaks or donor DNA. Nature. 2019;576:149-157.

Aymard F, Bugler B, Schmidt CK, Guillou E, Caron P, et al. Transcriptionally active chromatin recruits homologous recombination at DNA double-strand breaks. Nat Struct Mol Biol. 2014;21:366-374.

Bae S, Park J, Kim JS. Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics. 2014;30:1473-1475.

Chen PJ, Hussmann JA, Yan J, Knott GJ, Myber P, Qi LS, et al. Enhanced prime editing systems by manipulating cellular determinants of editing outcomes. Cell. 2021;184:5635-5652.

Concordet JP, Haeussler M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. Nucleic Acids Res. 2018;46:W242-W245.

Elliott B, Richardson C, Winderbaum J, Nickoloff JA, Jasin M. Gene conversion tracts from double-strand break repair in mammalian cells. Mol Cell Biol. 1998;18:93-101.

Gaudelli NM, Komor AC, Rees HA, Packer MS, Badran AH, Bryson DI, et al. Programmable base editing of A*T to G*C in genomic DNA without DNA cleavage. Nature. 2017;551:464-471.

Haapaniemi E, Botla S, Persson J, Schmierer B, Taipale J. CRISPR-Cas9 genome editing induces a p53-mediated DNA damage response. Nat Med. 2018;24:927-930.

Ihry RJ, Worringer KA, Salick MR, Frias E, Ho D, Theriault K, et al. p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. Nat Med. 2018;24:939-946.

Iyer S, Mir A, Vega-Badillo J, Roscoe BP, Ibraheim R, Zhu LJ, et al. Efficient homology-directed repair with circular single-stranded DNA donors. CRISPR J. 2022;5:685-701.

Komor AC, Kim YB, Packer MS, Zuris JA, Liu DR. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. Nature. 2016;533:420-424.

Kulcsar PI, Talas A, Huszar K, Ligeti Z, Toth E, Weinhardt N, et al. Crossing enhanced and high fidelity SpCas9 nucleases to optimize specificity and cleavage. Genome Biol. 2017;18:190.

Leibowitz ML, Papathanasiou S, Dober SA, Blaine LJ, Sun L, Yao Y, et al. Chromothripsis as an on-target consequence of CRISPR-Cas9 genome editing. Nat Genet. 2021;53:895-905.

Paquet D, Kwart D, Chen A, Sproul A, Jacob S, Teo S, et al. Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. Nature. 2016;533:125-129.

Richards S, Aziz N, Bale S, Bick D, Das S, Gastier-Foster J, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17:405-424.

---

*Prepared for Nature Methods. All benchmark numbers are from verified computational results. No experimental data were generated. All parameter values and their sources are documented in the repository.*
