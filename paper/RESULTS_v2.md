# Results

## 1. Overview of the CRISPRArchitect framework

CRISPRArchitect integrates sequence constraints, transcript context, and biological consequences into a unified pipeline for genome editing design (Fig. 1). The system takes as input one or more genomic variants and performs transcript-aware normalization using Ensembl annotations, mapping each variant to exon structure, coding frame, and splice-site proximity. Sequence context is extracted and validated against the reference genome, after which coding consequences and splice proximity are annotated. Feasibility is then evaluated independently for base editing, prime editing, and HDR based on modality-specific constraints, including PAM availability and editing window compatibility. Candidate strategies are generated and scored using a multi-objective framework that incorporates feasibility, biological consequences, and complexity. Strategies are ranked and reported with interpretable reasoning.

The v2 pipeline comprises seven stages implemented across 15 Python modules (~9,300 lines of code), building on the validated v1 framework (~24,000 lines). All modules are designed for Python 3.9+ with no dependencies beyond NumPy, SciPy, and Matplotlib. The system supports both SpCas9 (NGG PAM) and enFnCas9 (NRG PAM) as first-class nucleases, reflecting the broadened PAM compatibility developed in our laboratory.


## 2. Transcript-aware mapping enables context-specific interpretation of variants

To account for transcript-specific effects, CRISPRArchitect maps genomic variants to transcript coordinates and coding frames. This enables identification of exon boundaries, strand orientation, and codon position, allowing accurate determination of amino acid consequences and splice proximity. In contrast to genomic-coordinate-only approaches, transcript-aware mapping enables precise evaluation of functional outcomes, particularly for variants located near exon boundaries or within alternatively spliced regions (Fig. 2).

For each variant, the system fetches the canonical transcript from Ensembl, maps the genomic position to its CDS coordinate, and determines the affected codon. The reference allele is validated against the Ensembl genome sequence — in our 30-case benchmark, all successfully processed variants had confirmed reference alleles, catching potential genome-build or strand-orientation mismatches before they could propagate to downstream design errors.

Coding consequence annotation classifies each variant as synonymous, missense, nonsense, or splice-proximal (donor, acceptor, or region) following ACMG standards (Richards et al., 2015). Splice proximity is defined as ≤2 bp from the exon boundary for donor/acceptor sites and 3–8 bp for the splice region. This annotation layer directly feeds into the consequence-aware scoring system.


## 3. Feasibility varies across editing modalities and sequence context

We evaluated feasibility across editing modalities for a curated set of 30 variant scenarios representing diverse mutation types (Fig. 3). As expected, not all variants were compatible with all editing modalities.

Base editing feasibility requires three simultaneous conditions: (i) the correction must be a transition mutation (A→G for ABE or C→T for CBE), (ii) a suitable PAM site must position the target base within the editing window (positions 4–7 for ABE, 4–8 for CBE; Gaudelli et al., 2017; Komor et al., 2016), and (iii) bystander edits within the window should not introduce deleterious consequences. In our benchmark, PAM-verified base editing feasibility was narrower than mutation-type-only classification would suggest — several transitions lacked a guide that placed the target within the editing window, demonstrating the importance of PAM-aware feasibility checking.

Prime editing showed broader applicability, accommodating both transitions and transversions as well as small insertions (≤40 bp) and deletions (≤80 bp) (Anzalone et al., 2019). For each variant, the system designs a pegRNA with primer binding site (13 nt default, range 10–17) and reverse transcriptase template (10–30 nt), and searches for a PE3 nicking guide 40–100 bp away on the opposite strand.

HDR feasibility depends on nuclease guide availability near the target site. CRISPRArchitect evaluates cut-to-edit distance for each candidate guide and recommends donor type accordingly: ssODN for distances ≤30 bp, cssDNA for moderate distances (≤5,000 bp), and lssDNA or dsDNA for larger spans (Paquet et al., 2016; Iyer et al., 2022). Gene conversion probability is estimated using an exponential decay model calibrated to published tract-length distributions (Elliott et al., 1998).

These results highlight the importance of evaluating feasibility across modalities rather than relying on a single strategy.


## 4. Consequence-aware scoring influences strategy prioritization

A central claim of CRISPRArchitect is that incorporating biological consequences into the scoring framework alters strategy selection. The multi-objective scoring function balances five components:

$$\text{Score} = w_1 \cdot \text{Safety} + w_2 \cdot \text{Feasibility} - w_3 \cdot \text{Complexity} - w_4 \cdot \text{Risk} + w_5 \cdot \text{Confidence}$$

where default weights for iPSC applications are: Safety = 0.30, Feasibility = 0.25, Complexity = 0.20, Risk = 0.15, Confidence = 0.10.

Consequence-aware adjustments modify the base score through biologically grounded penalties and bonuses:

- **Bystander penalties**: Base editing strategies where bystander edits create missense changes are penalized (−0.10 per position); those creating nonsense changes receive a severe penalty (−0.25). Strategies where all bystanders are synonymous receive a small bonus (+0.05).
- **Splice proximity penalties**: Variants or bystander positions within 2 bp of a splice donor or acceptor site incur a −0.15 penalty; those within the splice region (3–8 bp) incur a −0.08 penalty.
- **DSB burden**: Strategies requiring ≥2 simultaneous DSBs in p53-active cells (such as iPSCs) incur an additional −0.10 penalty, reflecting the elevated risk of p53-mediated selection and chromosomal rearrangement (Ihry et al., 2018; Leibowitz et al., 2021).

In our benchmark, DSB-free strategies (base editing, prime editing) consistently scored higher than HDR strategies in iPSC context — for example, a typical base editing strategy scored 0.650 versus 0.388 for an equivalent HDR strategy at the same locus, with the difference driven primarily by the safety component (1.00 vs 0.50).


## 5. CRISPRArchitect improves selection of biologically consistent strategies

We benchmarked CRISPRArchitect on a curated dataset of 30 variant scenarios with predefined preferred, acceptable, and reject strategy labels (Fig. 5, Fig. 7). The dataset spans 11 categories: clean base-editable substitutions (n = 7), base-editing negative control (n = 1), prime-editable transversions (n = 2), prime-editable small indels (n = 5), HDR-required large deletions (n = 3), HDR/PE small deletions (n = 1), compound heterozygous cases requiring hybrid strategies (n = 5), sequential HDR (n = 1), dual base editing (n = 3), and edge cases including distant variants and non-coding positions (n = 2).

All benchmark variants use verified GRCh38 coordinates from ClinVar and published editing studies, with reference alleles validated against the Ensembl genome.

**Top-1 accuracy** (top-ranked strategy matches a preferred or acceptable label): **86.7%** (26/30 cases). The system correctly prioritized DSB-free approaches for all base-editable and prime-editable cases (15/15), all compound heterozygous cases (9/9), and both edge cases. The four top-1 misses occurred exclusively in large deletion cases (n = 3) and one sequential HDR case (n = 1), where the system recommended prime editing — a strategy that is technically feasible for the individual variants but does not address the multi-exon structural nature of the deletion.

**Top-3 accuracy** (at least one preferred or acceptable strategy appears in the top 3): **96.7%** (29/30 cases). In 29 of 30 cases, a biologically appropriate strategy appeared among the top three ranked options. The single miss was a compound heterozygous COL7A1 case where the expected sequential HDR strategy was not generated because the pipeline processed each variant independently rather than as a coordinated pair.

**Rejection accuracy** (infeasible strategies from the reject label are absent from the top-ranked output): **90.0%** (27/30 cases). The system correctly excluded strategies requiring simultaneous dual DSBs or inappropriate modalities in 27 cases. The three cases where rejected strategies were not excluded were large deletion cases where HDR was expected but the pipeline's current heuristic does not yet model multi-exon deletion-specific constraints.

**Category-level analysis** revealed that CRISPRArchitect achieved 100% top-1 accuracy for clean base-editable cases (7/7), prime-editable transversions (2/2), prime-editable small indels (5/5), compound heterozygous hybrid cases (5/5), compound heterozygous dual BE cases (3/3), and edge cases (2/2). The only categories with imperfect top-1 accuracy were HDR-required large deletions (0/3 top-1, 3/3 top-3) and compound heterozygous sequential HDR (0/1 top-1, 0/1 top-3).


## 6. Explicit rejection of infeasible and high-risk strategies

In addition to ranking feasible strategies, CRISPRArchitect explicitly identifies infeasible or high-risk designs. Strategies lacking PAM availability or violating editing constraints are filtered during feasibility checking, while those introducing deleterious consequences are penalized during scoring. This explicit rejection mechanism improves interpretability and prevents misleading recommendations.

For the HBB sickle cell variant (c.20A>T, p.Glu7Val), a transversion, the system correctly identified that base editing is infeasible — ABE performs A→G corrections and CBE performs C→T corrections, neither of which addresses a T→A correction — and recommended prime editing as the appropriate modality. This demonstrates the system's ability to distinguish between mutation types that superficially resemble base-editing targets but are fundamentally incompatible.

For the non-coding FMR1 5′UTR variant (case 30), the system appropriately recommended HDR as the top strategy, recognizing that the variant's non-coding location limits the applicability of consequence-aware prioritization.


## 7. Robustness and reproducibility

The pipeline incorporates automatic retry logic with exponential backoff for transient Ensembl API failures (HTTP 500, 502, 503, 504, and network timeouts), retrying up to 3 times with waits of 1, 2, and 4 seconds. In our benchmark run, 4 of 120 total API calls initially failed due to transient server errors; all were recovered automatically on retry, resulting in zero pipeline failures across all 30 cases.

All benchmark coordinates were verified against the Ensembl GRCh38 genome prior to evaluation. The benchmark dataset, pipeline code, scoring parameters, and evaluation scripts are publicly available to enable independent reproduction of all reported results.


---

## Figure Legends

**Figure 1. CRISPRArchitect pipeline architecture.** Schematic of the seven-stage pipeline: genomic variant input → transcript mapping → reference validation → coding annotation → modality-specific feasibility (BE, PE, HDR) → strategy generation → multi-objective scoring and ranking.

**Figure 2. Transcript-aware variant mapping.** Illustration of genomic-to-transcript coordinate mapping, showing exon identification, CDS position computation, codon frame determination, and splice distance calculation for forward- and reverse-strand genes.

**Figure 3. Feasibility across editing modalities.** Heatmap showing the feasibility of base editing, prime editing, and HDR for each of the 30 benchmark cases. Green (+) indicates the modality produced the top-ranked strategy; amber (?) indicates the modality was available but not top-ranked; white indicates not feasible. Not all variants are editable by all methods, highlighting the need for unified cross-modality evaluation.

**Figure 5. Strategy selection accuracy by variant category.** Grouped bar chart showing top-1 (blue) and top-3 (green) accuracy for each of 11 variant categories in the benchmark dataset. CRISPRArchitect achieved 100% accuracy for base-editable, prime-editable, and compound heterozygous categories. Large deletion cases showed lower top-1 accuracy (0/3) but high top-3 accuracy (3/3).

**Figure 7. Overall benchmark performance.** (A) Top-1 accuracy (86.7%), top-3 accuracy (96.7%), and rejection accuracy (90.0%) across all 30 benchmark cases. (B) Distribution of top-ranked strategy types: prime editing dominated (28/30 cases), reflecting the pipeline's safety-first scoring for iPSC context. (C) Infeasible strategy rejection: 27/30 cases correctly excluded rejected modalities from the top-ranked output.
