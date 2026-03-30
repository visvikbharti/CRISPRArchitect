# CRISPRArchitect v3: multi-nuclease decision support for genome editing strategy design with TOPSIS ranking and sensitivity analysis

Vishal Bharti^1^ and Debojyoti Chakraborty^1,\*^

^1^ CSIR-Institute of Genomics and Integrative Biology, New Delhi, India

^\*^ Corresponding author. Email: debojyoti@igib.in

---

## Abstract

Selecting an optimal genome editing strategy for a clinically relevant variant requires simultaneous evaluation of nuclease availability, editor compatibility, biological consequences, and off-target risk --- constraints that span multiple editing modalities and are rarely assessed within a unified framework. Here we present CRISPRArchitect v3, a computational decision-support platform that integrates multi-nuclease feasibility evaluation, consequence-aware scoring, and formal multi-criteria ranking to recommend editing strategies across base editing, prime editing, and homology-directed repair (HDR). Building on our transcript-aware pipeline (v2), v3 introduces three capabilities: (i) a multi-nuclease feasibility engine that evaluates five Cas nucleases (SpCas9, enFnCas9, SpCas9-NG, SpRY, Cas12a) paired with nine base editor profiles (ABE7.10, ABE8e, BE4max, and six nuclease-specific fusions: ABE8e-enFnCas9, ABE8e-SpCas9-NG, ABE8e-SpRY, BE4max-enFnCas9, BE4max-SpCas9-NG, BE4max-SpRY), rescuing base editing from 0/30 to 6/30 top-ranked cases through broader PAM access and expanded editing windows; (ii) TOPSIS-based multi-criteria ranking with Monte Carlo sensitivity analysis (10,000 Dirichlet-sampled weight permutations) that provides rank stability reporting alongside strategy scores; and (iii) CFD and MIT off-target specificity scoring integrated into guide selection. We also implement an HGVS notation parser and ClinVar batch ingestion for clinical variant input. Across 30 ClinVar benchmark cases with verified GRCh38 coordinates, CRISPRArchitect v3 achieved 86.7% top-1 accuracy (26/30), 96.7% top-3 accuracy (29/30), and 86.7% rejection accuracy (26/30), with a strategy distribution of 20% base editing, 77% prime editing, and 3% HDR. The codebase comprises approximately 36,500 lines of Python across the v3 core (~12,500 lines) and v1 modules (~24,000 lines), with 224 passing tests. CRISPRArchitect is positioned as a decision-support tool rather than a predictive optimizer, providing transparent, interpretable recommendations to inform --- not replace --- experimental judgment.

---

## Introduction

The genome editing toolkit has expanded rapidly in the past decade, with base editing^1,2^, prime editing^3^, and homology-directed repair (HDR) now offering complementary approaches to correct pathogenic variants. Cytosine base editors (CBEs) convert C-to-T within a defined activity window^1^, adenine base editors (ABEs) convert A-to-G^2^, and the more recently developed ABE8e extends this capability with a broader window (positions 3--9) and higher on-target activity^11^. Prime editors use a reverse transcriptase template to install a wider range of edits without requiring double-strand breaks (DSBs)^3^, while HDR enables precise sequence replacement but depends on cellular repair pathways and introduces DSB-associated risks including p53-mediated toxicity in iPSCs^5^ and chromosomal rearrangements. Each modality is constrained by distinct sequence requirements: base editing demands a PAM site that positions the target nucleotide within a narrow deaminase activity window; prime editing requires pegRNA design with primer binding site and RT template optimization; and HDR efficiency decays sharply with increasing cut-to-edit distance^4^.

A further dimension of complexity arises from nuclease diversity. While SpCas9 (NGG PAM) remains the most widely used nuclease, its PAM restriction limits the fraction of genomic positions accessible for editing. Engineered variants with relaxed PAM requirements --- including enFnCas9 (NRG PAM)^14^, SpCas9-NG (NG PAM)^12^, and the near-PAMless SpRY (NNN PAM)^13^ --- substantially expand targetable sequence space but differ in on-target activity, off-target propensity, and compatibility with specific base editors. The combinatorial space of nuclease-editor-modality-variant configurations is too large for manual evaluation, yet existing computational tools largely address individual modalities or single nucleases in isolation. CRISPOR provides guide RNA design and off-target scoring for nuclease selection^20^ but does not evaluate base editing or prime editing feasibility. Modality-specific design tools exist for base editing and prime editing but do not compare across modalities or account for the expanded nuclease landscape.

The absence of a unified decision framework has practical consequences. Our previous analysis (v2) revealed that PAM-dependent editing window constraints are a more significant bottleneck for base editing than mutation-type classification alone: at ABE-compatible transitions, SpCas9 NGG PAM sites frequently failed to position the target nucleotide within the narrow ABE7.10 window (positions 4--7). This observation motivated the multi-nuclease approach implemented in v3, where broader PAM nucleases and wider-window editors can rescue otherwise infeasible base editing designs. Additionally, strategy ranking requires balancing multiple objectives --- safety, feasibility, complexity, off-target risk, and confidence --- that are inherently incommensurable. Simple weighted sums, while intuitive, do not penalize strategies that are catastrophically poor on a single axis. Formal multi-criteria decision methods, such as TOPSIS (Technique for Order Preference by Similarity to Ideal Solution)^6^, address this limitation by ranking alternatives based on their geometric distance to ideal and anti-ideal solutions in normalized decision space.

Several computational tools address individual aspects of editing design. CRISPOR^20^ provides guide RNA scoring and off-target prediction for a range of nucleases but does not evaluate base editing or prime editing feasibility, nor does it compare across modalities. PrimeDesign (Hsu et al., 2021) generates pegRNA designs for prime editing but does not assess whether base editing might be simpler at the same locus. BE-Designer (Hwang et al., 2018) identifies base editing-compatible sites for a single nuclease-editor pair but does not consider alternative nucleases, expanded editing windows, or bystander consequence severity. CRISPick (Doench & Root, 2021) provides on-target activity scoring but is limited to guide design, not strategy-level decision-making. To our knowledge, no existing tool systematically compares base editing, prime editing, and HDR across multiple nuclease platforms, evaluates biological consequences of bystander edits, and provides formal uncertainty quantification for strategy rankings.

Here we present CRISPRArchitect v3, which extends our transcript-aware, consequence-guided framework with three capabilities: (i) a multi-nuclease feasibility engine spanning five nucleases and nine base editor profiles; (ii) TOPSIS-based strategy ranking with Monte Carlo sensitivity analysis; and (iii) off-target specificity scoring using CFD^15^ and MIT^16^ frameworks. We also introduce an HGVS notation parser^18^ and ClinVar batch ingestion for clinical variant input. We evaluate the platform on 30 ClinVar benchmark cases and demonstrate that multi-nuclease evaluation rescues base editing from 0/30 to 6/30 top-ranked cases, while TOPSIS sensitivity analysis provides rank stability reporting that is, to our knowledge, not available in any existing CRISPR design tool. We frame CRISPRArchitect explicitly as a decision-support platform: its recommendations are intended to inform experimental design, not to substitute for empirical validation.

---

## Results

### Overview of the v3 framework

CRISPRArchitect v3 is a modular pipeline that takes genomic variants as input and produces ranked, annotated editing strategy recommendations across base editing, prime editing, and HDR (Fig. 1). The pipeline operates in seven stages: (1) variant normalization and HGVS parsing; (2) transcript mapping via the Ensembl REST API; (3) reference allele validation against GRCh38; (4) coding consequence and splice-site annotation following ACMG standards^8^; (5) multi-nuclease, multi-editor feasibility evaluation; (6) TOPSIS-based scoring with consequence adjustments; and (7) Monte Carlo sensitivity analysis and reporting.

The v3 codebase comprises approximately 12,500 lines of new Python code building on the v1 framework (~24,000 lines), with 224 passing tests. All modules are Python 3.9 compatible with dependencies limited to NumPy, SciPy, and Matplotlib. The system supports five nucleases (SpCas9, enFnCas9, SpCas9-NG, SpRY, Cas12a) as first-class options with nine base editor profiles (ABE7.10, ABE8e, ABE8e-enFnCas9, ABE8e-SpCas9-NG, ABE8e-SpRY, BE4max, BE4max-enFnCas9, BE4max-SpCas9-NG, BE4max-SpRY), each with defined editing windows and evidence tiers. Clinical variant input is supported via an HGVS parser that accepts NM_xxx:c.NNNRef>Alt notation, and batch analysis is enabled through TSV and VCF ingestion from ClinVar exports.

### Multi-nuclease feasibility engine

A central limitation identified in our v2 analysis was that base editing feasibility, when restricted to SpCas9 (NGG PAM) and ABE7.10 (window positions 4--7), was zero across all 30 benchmark cases despite seven cases carrying ABE-compatible transitions. The v3 multi-nuclease engine addresses this by systematically evaluating each variant against all compatible nuclease-editor combinations (Fig. 2, Fig. 3).

The key changes are twofold. First, enFnCas9 (NRG PAM, where R = A or G)^14^ approximately doubles the number of candidate PAM sites relative to SpCas9 (NGG), because NRG encompasses both NAG and NGG. Second, ABE8e^11^ expands the editing window from positions 4--7 (ABE7.10) to positions 3--9, nearly tripling the number of protospacer positions at which a target adenine is within the deaminase activity window. The combination of broader PAM and wider window substantially increases the probability that at least one guide places the target within the editing window.

In the 30-case benchmark, this multi-nuclease evaluation rescued base editing from 0/30 top-ranked cases (v2, SpCas9 + ABE7.10 only) to 6/30 top-ranked cases (v3, multi-nuclease) (Fig. 2). The overall strategy distribution shifted from effectively 97% prime editing in v2 to a more balanced profile: base editing 20% (6/30), prime editing 77% (23/30), and HDR 3% (1/30). This redistribution reflects the genuine expansion of the base editing-accessible target space, not a change in scoring weights.

However, prime editing remains dominant (23/30 top-ranked), and we note that the PAM-window constraint is not fully eliminated even with broader nucleases. At many clinically relevant loci, no nuclease-editor combination places the target within the editing window. This is a genuine biological limitation, not a software deficiency: the co-occurrence of a suitable PAM at the precise spacing required for a narrow editing window is locus-dependent and not guaranteed even with near-PAMless nucleases such as SpRY, which have reduced on-target activity^13^ that lowers their TOPSIS scores.

To illustrate the mechanism of rescue, consider case BE_BRCA2_004, a G>A transition in BRCA2 at a ClinVar-verified pathogenic locus. In v2, the ABE7.10 window (positions 4--7) with SpCas9 (NGG PAM) found no guide placing the target adenine within the editing window --- the closest guide positioned the target at protospacer position 3, one position outside the window boundary. In v3, ABE8e's expanded window (positions 3--9) captures this position, and additionally, an enFnCas9 (NRG PAM) guide that was invisible to the NGG-only scanner provides an alternative placement at position 5. The result is that two independent editor-nuclease combinations now support base editing at this locus (ABE8e + SpCas9 and ABE8e + enFnCas9), both with Tier A and Tier B evidence respectively. The TOPSIS scorer ranks base editing above prime editing for this case because the safety component (1.0 for DSB-free BE versus 1.0 for DSB-free PE, but with ABE8e's higher modality prior of 0.95 versus PE's 0.85) and the reduced complexity (single guide, no donor) tip the balance. Cases BE_TSC2_006 and BE_FBN1_007 follow a similar pattern, where the ABE8e window expansion alone was sufficient to capture the target.

For the three remaining rescued cases (HDR_BRCA2_019, COMP_TSC_025, COMP_PKD1_027), the mechanism differed: these involved variants where the correction direction maps to CBE (C>T on the protospacer), and BE4max paired with enFnCas9 provided the rescuing PAM site. In all six cases, the primary driver was the combination of ABE8e's broader window and enFnCas9's NRG PAM, rather than the near-PAMless nucleases (SpCas9-NG, SpRY), which contributed additional options but at lower efficiency modifiers.

Each editor-nuclease combination is assigned an evidence tier. Combinations with direct experimental validation (e.g., ABE8e + SpCas9) are Tier A; combinations inferred from component-level data (e.g., ABE8e + SpRY, BE4max + enFnCas9) are Tier B. Evidence tiers are reported alongside strategy recommendations and are factored into the confidence component of the TOPSIS scoring. This transparency ensures that users can distinguish between well-validated and extrapolated recommendations when making experimental decisions.

### TOPSIS multi-criteria scoring with sensitivity analysis

Strategy ranking in v3 replaces the simple weighted-sum scorer (v2) with TOPSIS^6^, a formal multi-criteria decision analysis method (Fig. 4). TOPSIS ranks alternatives based on their Euclidean distance to the ideal solution (best score on every criterion) and anti-ideal solution (worst score on every criterion) in normalized, weighted decision space. The relative closeness score C = D^-^ / (D^+^ + D^-^) ranges from 0 (identical to anti-ideal) to 1 (identical to ideal).

TOPSIS is more principled than a weighted sum for two reasons. First, it penalizes strategies that are catastrophically poor on any single axis, even if they score well on others --- a strategy with high feasibility but zero safety is pushed toward the anti-ideal, whereas a weighted sum might still produce a moderate composite score. Second, it operates in normalized space, making it invariant to the absolute scale of individual criteria.

The five scoring criteria and their default weights (optimized for iPSC applications) are: Safety (0.30), Feasibility (0.25), Complexity (0.20), Risk (0.15), and Confidence (0.10). Safety scores reflect DSB burden: DSB-free approaches (base editing, prime editing) receive 1.0; single-DSB approaches (HDR) receive 0.5; dual-DSB approaches receive 0.2, reflecting documented p53-mediated toxicity in iPSCs^5^. Consequence-aware adjustments are applied after TOPSIS scoring: bystander missense penalties (-0.10 per position), bystander nonsense penalties (-0.25), and splice donor/acceptor penalties (-0.15) modify the TOPSIS score within the [0, 1] range.

To assess rank stability under weight uncertainty, v3 implements Monte Carlo sensitivity analysis using 10,000 Dirichlet-sampled weight permutations. For each permutation, a weight vector is drawn from a Dirichlet distribution centered on the default weights (concentration parameter = 10), TOPSIS is re-run, and the rank of each strategy is recorded. The output includes, for each strategy: the fraction of permutations in which it is top-ranked (rank stability), the mean rank across permutations, and the full rank distribution. This provides an explicit measure of how sensitive the recommendation is to weight assumptions --- to our knowledge, no existing CRISPR design tool provides this form of rank stability reporting. Strategies with high rank stability (e.g., >0.90 fraction top-ranked across 10,000 permutations) represent robust recommendations; strategies with low rank stability signal that the choice between alternatives is sensitive to weight assumptions and warrants additional experimental consideration.

### Off-target specificity scoring

CRISPRArchitect v3 integrates off-target specificity scoring into guide evaluation using two complementary frameworks (Tier B: literature-informed heuristic). The CFD (Cutting Frequency Determination) score^15^ quantifies the fraction of cutting activity retained at a mismatched off-target site, based on position-specific mismatch tolerance derived from large-scale guide activity profiling. For each candidate guide, the system computes an aggregate specificity score by scanning a local sequence window for potential off-target sites, scoring each mismatch pattern using the CFD position-weight matrix, and aggregating into a specificity score where higher values indicate fewer predicted off-targets.

The MIT specificity score^16^ provides an alternative framework based on the Hsu et al. (2013) mismatch tolerance model, which uses a different position-weight matrix and aggregation scheme. Both scores are computed for each candidate guide and reported in the strategy output; the CFD score is used as the default input to the Risk component of the TOPSIS scoring.

We acknowledge that this approach uses a simplified off-target enumeration limited to a local sequence window, not genome-wide search. For clinical applications, we recommend supplementing CRISPRArchitect off-target scores with genome-wide predictions from dedicated tools such as Cas-OFFinder or CRISPOR^20^. Furthermore, the CFD and MIT matrices are trained on SpCas9 data; their accuracy for alternative nucleases (enFnCas9, SpCas9-NG, SpRY) is extrapolated and should be interpreted with caution. The position tolerance profile, showing decreasing mismatch tolerance from PAM-distal to PAM-proximal positions, is visualized in Fig. S2.

### Benchmark evaluation

We evaluated CRISPRArchitect v3 on 30 ClinVar variant scenarios with verified GRCh38 coordinates and predefined truth labels (preferred, acceptable, reject) spanning diverse mutation types and clinical contexts (Fig. 5, Fig. 6).

**Top-1 accuracy** (top-ranked strategy matches a preferred or acceptable label): **86.7%** (26/30). The system correctly prioritized DSB-free approaches for base-editable and prime-editable cases, compound heterozygous cases, and edge cases. The four top-1 misses occurred in large deletion cases (n = 3) and one sequential HDR case (n = 1), where the system recommended prime editing --- technically feasible for the individual variant but inappropriate for the multi-exon structural nature of the deletion.

**Top-3 accuracy** (at least one preferred or acceptable strategy in the top 3): **96.7%** (29/30). The single miss was a compound heterozygous COL7A1 case where the expected sequential HDR strategy was not generated because the pipeline processes each variant independently.

**Rejection accuracy** (infeasible strategies from the reject label absent from top-ranked output): **86.7%** (26/30). The four cases where rejected strategies were not excluded include the three large deletion cases and one case where a multi-step strategy was ranked above the expected single-step alternative.

**Strategy distribution** across the 30 benchmark cases: base editing 20% (6/30), prime editing 77% (23/30), HDR 3% (1/30). Compared to v2, where prime editing was top-ranked in 29/30 cases and base editing in 0/30, the v3 multi-nuclease engine achieved a meaningful redistribution by rescuing six cases where ABE8e paired with broader-PAM nucleases (primarily enFnCas9) placed the target within the expanded editing window.

**Comparison to v2.** The overall top-1 and top-3 accuracies are identical between v2 and v3 (86.7% and 96.7%, respectively), because the multi-nuclease rescue affects strategy *type* (BE vs. PE) but not correctness --- in all six rescued cases, both base editing and prime editing are acceptable strategies in the truth labels. The substantive improvement is that v3 provides a more diverse and informative set of recommendations: where v2 effectively defaulted to prime editing for nearly all cases, v3 identifies base editing as a viable or preferred option at loci where the biology supports it. Rejection accuracy changed from 90.0% in v2 (using different rejection criteria) to 86.7% in v3 (with stricter evaluation of multi-step strategy handling).

### Explicit rejection of infeasible strategies

CRISPRArchitect explicitly identifies and reports infeasible strategies rather than silently omitting them. For each variant, modalities that fail feasibility checks (no compatible PAM, target outside editing window, edit size exceeding modality limits) are listed as rejected with a human-readable reason. This is particularly important for base editing, where the difference between "nominally ABE-compatible transition" and "PAM-verified, window-verified base editing" can determine whether a strategy is feasible.

For the HBB sickle cell variant (c.20A>T, p.Glu7Val), the system correctly identified that base editing is infeasible (the required T-to-A correction is a transversion, incompatible with both ABE and CBE) and recommended prime editing. For non-coding variants, the system appropriately defaults to HDR when consequence-aware scoring has limited applicability. Rejected strategies carry Tier C evidence labels, signaling that they should not be pursued without independent experimental justification.

---

## Discussion

CRISPRArchitect v3 demonstrates that systematic multi-nuclease evaluation, combined with formal multi-criteria ranking and sensitivity analysis, can provide informative decision support for genome editing strategy selection. The multi-nuclease feasibility engine rescued base editing from 0/30 to 6/30 top-ranked cases by evaluating ABE8e with broader-window editors and enFnCas9/SpCas9-NG/SpRY with relaxed PAM requirements. This finding has implications beyond the software itself: it suggests that estimates of base editing applicability based solely on mutation-type cataloguing --- for example, the fraction of ClinVar pathogenic variants that are transitions --- may overestimate the fraction correctable with canonical ABE7.10 + SpCas9, while simultaneously underestimating the fraction accessible through next-generation editors and alternative nucleases. The practical message is that nuclease-editor combinatorics matter as much as mutation chemistry.

The adoption of TOPSIS scoring with Monte Carlo sensitivity analysis represents a methodological contribution to the CRISPR computational tool landscape. While TOPSIS is well-established in operations research^6^, its application to genome editing strategy ranking --- and specifically the use of Dirichlet-sampled weight perturbations to report rank stability --- provides a form of uncertainty quantification that we have not encountered in existing CRISPR design tools. When a strategy is top-ranked in >90% of 10,000 weight permutations, the recommendation is robust to reasonable differences in how users prioritize safety versus feasibility versus complexity. When rank stability is low, this explicitly signals that the choice is weight-sensitive and warrants careful experimental comparison. We consider this transparency valuable, particularly in therapeutic contexts where the consequences of suboptimal strategy selection can be severe.

Several findings from the benchmark merit discussion. First, prime editing remains the dominant modality (77% of top-ranked strategies), reflecting its inherent advantages: no editing window constraint, broad mutation-type compatibility, and DSB-free operation. The safety-first weight configuration (Safety = 0.30) amplifies this advantage in the iPSC context; applications in cell types with attenuated p53 responses (e.g., HEK293T, where p53 pathways are often disrupted) would shift the balance toward HDR, as the safety penalty for DSBs would be substantially reduced. Users can adjust the weight configuration through the TOPSIS scorer to match their specific cell-type context and experimental priorities — and the sensitivity analysis will report how such changes affect the ranking.

Second, the six base editing rescues all involved ABE8e rather than ABE7.10, consistent with the broader window (positions 3--9 versus 4--7) being the primary driver rather than PAM diversity alone. This has a practical implication: for researchers planning base editing experiments, the choice of editor version (ABE8e versus ABE7.10) may matter more than the choice of nuclease for many loci. The enFnCas9 NRG PAM contributes additional rescues, particularly at loci where SpCas9 NGG sites are sparse. Given that enFnCas9 was developed with enhanced specificity and improved HDR characteristics^14^, it represents a particularly attractive option for iPSC editing workflows.

Third, HDR remained top-ranked in only one case (3%), reflecting the substantial safety penalties in the iPSC-optimized configuration. This distribution is context-specific and would change with different weight configurations or cell types. We note that CRISPRArchitect's value is not in producing a single "correct" answer but in systematically enumerating the option space and making the trade-offs explicit. A researcher who prioritizes speed over safety, or who works in a p53-null cell line, would obtain different rankings from the same pipeline by adjusting the weight profile.

The comparison with existing tools highlights CRISPRArchitect's distinct positioning. CRISPOR^20^ excels at guide design and off-target prediction but evaluates one nuclease at a time without cross-modality comparison. PrimeDesign generates pegRNA designs but cannot assess whether base editing might be simpler and safer at the same locus. BE-Designer identifies base editing-compatible sites but does not consider consequence severity of bystander edits or compare with HDR alternatives. CRISPRArchitect does not aim to replace any of these tools for their specific design functions; rather, it operates upstream — at the strategy selection level — to help researchers decide which modality to pursue before committing to detailed guide or donor design. The recommended workflow is to use CRISPRArchitect for strategy triage, then employ modality-specific tools (CRISPOR for guide optimization, PrimeDesign for pegRNA design, donor synthesis tools for HDR templates) for detailed design of the selected strategy.

We acknowledge several limitations explicitly. First and most importantly, CRISPRArchitect has not been validated experimentally. The benchmark evaluates the system against expert-defined truth labels, not against experimental editing outcomes. The tool is designed as decision support: it narrows the search space and provides structured reasoning, but it does not predict editing efficiency or clinical success. Second, the off-target scoring uses simplified local enumeration rather than genome-wide search, and the CFD/MIT matrices are SpCas9-derived; their applicability to enFnCas9, SpCas9-NG, SpRY, and Cas12a is extrapolated (Tier B). Third, the system processes variants independently and does not model multi-variant coordination (e.g., compound heterozygous cases requiring sequential editing rounds). This accounts for the single top-3 miss in the benchmark. Fourth, the coding annotation module uses a spliced-exon surrogate model that is adequate for canonical exonic variants but would require refinement for non-canonical transcript architectures. Fifth, large structural deletions are handled poorly (0/3 top-1 accuracy for HDR-required large deletion cases), indicating that explicit deletion-aware logic is needed. Sixth, chromatin accessibility, replication timing, and epigenomic context are not incorporated; all of these are known to influence editing efficiency in a locus-specific manner.

Future extensions that would meaningfully strengthen the framework include: (i) prospective experimental validation in iPSC disease models to establish the concordance between computational rankings and observed editing outcomes; (ii) genome-wide off-target prediction via integration with Cas-OFFinder or CRISPOR; (iii) chromatin accessibility integration from ENCODE or ATAC-seq datasets for locus-specific efficiency adjustment; (iv) multi-variant coordination for compound heterozygous and complex structural variant cases; and (v) machine learning-based scoring refinement as experimental outcome data accumulates. We release CRISPRArchitect as an open-source tool and invite community contributions toward these goals.

In summary, CRISPRArchitect v3 provides a transparent, extensible decision-support platform for genome editing strategy design. The multi-nuclease engine expands the accessible target space for base editing, TOPSIS scoring with sensitivity analysis provides principled ranking with uncertainty quantification, and off-target specificity scoring adds a safety dimension to guide selection. We emphasize that CRISPRArchitect is a decision-support tool: it is designed to structure and inform the strategy selection process, not to replace the experimental judgment that remains essential for safe and effective genome editing.

---

## Methods

### Variant normalization and transcript mapping

CRISPRArchitect accepts genomic variants in three formats: (i) direct genomic coordinates (chromosome, position, reference allele, alternate allele in GRCh38), (ii) HGVS coding DNA notation (NM_xxx:c.NNNRef>Alt)^18^, and (iii) ClinVar batch files (TSV or VCF format). The HGVS parser uses regular expression matching to extract transcript accession, CDS position, and allele information, then resolves transcript-relative coordinates to genomic coordinates via the Ensembl REST API. Supported HGVS formats include substitutions (c.910C>T), single-base and multi-base deletions (c.1234del, c.1234_1236del), insertions (c.1234_1235insATG), and deletion-insertions (c.1234delinsATG). Gene-symbol prefixed notation (e.g., NF1:c.910C>T) is also supported, with transcript resolution via Ensembl.

Each variant is normalized and mapped to transcript context using the Ensembl REST API (https://rest.ensembl.org). The canonical transcript for the affected gene is selected via the Ensembl `vep/human/region` endpoint, which returns transcript consequences ranked by canonical status, MANE Select designation, and transcript length. Genomic coordinates are 1-based, consistent with VCF convention. For genes on the reverse strand, reference and alternate alleles are reverse-complemented prior to comparison with the coding strand sequence. CDS position, exon number, and codon frame are determined by mapping the genomic coordinate to the exon structure of the selected transcript.

Reference validation is performed by fetching a window of genomic sequence from the Ensembl `sequence/region` endpoint and comparing the reported reference allele against the retrieved sequence, with appropriate complementation for reverse-strand genes. Variants whose reference allele does not match the genome are flagged and excluded from downstream analysis. In the 30-case benchmark, all variants passed reference validation.

### Coding and splice-site annotation

Coding consequences are determined by translating reference and alternate codons using the standard genetic code. Variants are classified as synonymous (amino acid unchanged), missense (different amino acid), nonsense (premature stop codon), frameshift (insertion or deletion not a multiple of three), or splice-proximal. Splice proximity follows ACMG standards^8^: positions within 2 bp of an exon boundary are classified as splice donor (5' end) or splice acceptor (3' end); positions 3--8 bp from the boundary are classified as splice region.

### Multi-nuclease base editing feasibility

Base editing feasibility is evaluated through systematic enumeration of all nuclease-editor combinations. For each variant, the system determines whether the required correction is a transition mutation compatible with ABE (A>G or T>C on the protospacer strand) or CBE (C>T or G>A). If the mutation type is compatible, the system scans for PAM sites across all supported nucleases:

- **SpCas9**: NGG PAM (Jinek et al., 2012)
- **enFnCas9**: NRG PAM (R = A or G)^14^
- **SpCas9-NG**: NG PAM^12^
- **SpRY**: NNN PAM (near-PAMless)^13^
- **Cas12a**: TTTV PAM (V = A, C, or G)

For each PAM site, the system evaluates whether the target nucleotide falls within the editing window of each compatible editor:

- **ABE7.10**: positions 4--7 (1-indexed, PAM-distal = 1)^2^
- **ABE8e**: positions 3--9^11^
- **BE4max**: positions 4--8^17^

The editing window check uses the patient allele: the system verifies that the nucleotide to be corrected (the pathogenic allele) is positioned within the deaminase activity window of the protospacer. Bystander edits are identified by scanning all C residues (for CBE) or A residues (for ABE) within the editing window and classifying the consequence of each potential bystander conversion using the coding annotation module^9,10^. Bystander consequences (synonymous, missense, or nonsense) directly feed into the consequence-aware scoring system.

Editor-nuclease compatibility is governed by a defined matrix. ABE8e is compatible with SpCas9 (efficiency modifier 1.0), enFnCas9 (0.8), SpCas9-NG (0.7), and SpRY (0.5). BE4max is compatible with SpCas9 (1.0), enFnCas9 (0.8), SpCas9-NG (0.7), and SpRY (0.5). Efficiency modifiers scale the feasibility score and are incorporated into the TOPSIS input matrix. Each combination carries an evidence tier: Tier A for combinations with direct experimental support (ABE8e + SpCas9, BE4max + SpCas9, ABE7.10 + SpCas9), and Tier B for combinations inferred from component-level data.

### Prime editing and HDR feasibility

Prime editing feasibility requires identification of a suitable protospacer with a PAM site near the target. For each candidate guide, the system designs a pegRNA comprising a primer binding site (PBS; default 13 nt, range 10--17 nt, GC content constrained to 40--60%) and a reverse transcriptase (RT) template (range 10--30 nt) encoding the desired edit. A PE3 nicking guide is searched on the opposite strand, 40--100 bp from the pegRNA-directed nick. Prime editing is considered feasible for substitutions, insertions up to 40 bp, and deletions up to 80 bp^3^.

HDR feasibility requires a guide RNA directing a DSB near the target site. Cut-to-edit distance is scored using an exponential decay function with an empirical half-life of approximately 20 bp^4^. Donor type is recommended based on edit size and distance: ssODN for distances up to 30 bp (90 bp homology arms), cssDNA for edits within 5,000 bp (300 bp arms)^7^, lssDNA (300 bp arms) or dsDNA (800 bp arms) for larger spans. Gene conversion probability is estimated using an exponential decay model calibrated to published tract-length distributions^19^.

### TOPSIS scoring algorithm

Strategies are scored using the TOPSIS (Technique for Order Preference by Similarity to Ideal Solution) method^6^. For each variant, all feasible strategies form the alternatives; the criteria are Safety, Feasibility, Complexity, Risk, and Confidence.

The algorithm proceeds in five steps:

1. **Construct the decision matrix** X where x_ij is the score of alternative i on criterion j. Safety: 1.0 (DSB-free), 0.5 (single DSB), 0.2 (dual DSB). Feasibility: modality-specific confidence incorporating PAM quality, window position, and editor-nuclease efficiency. Complexity: inverse of the number of editing steps and components. Risk: based on off-target specificity (CFD score) and consequence penalties. Confidence: based on evidence tier and literature support.

2. **Normalize** each column by its Euclidean norm: r_ij = x_ij / sqrt(sum_k x_kj^2).

3. **Apply weights**: v_ij = w_j * r_ij, where default weights are Safety = 0.30, Feasibility = 0.25, Complexity = 0.20, Risk = 0.15, Confidence = 0.10.

4. **Determine ideal (A+) and anti-ideal (A-) solutions**: A+ = max_i(v_ij) for benefit criteria (Safety, Feasibility, Confidence); A- = min_i(v_ij) for benefit criteria. Complexity and Risk are cost criteria (lower is better), so their ideal/anti-ideal assignments are reversed.

5. **Compute relative closeness**: D_i+ = sqrt(sum_j (v_ij - A_j+)^2), D_i- = sqrt(sum_j (v_ij - A_j-)^2), C_i = D_i- / (D_i+ + D_i-).

Consequence-aware adjustments are applied after TOPSIS scoring. Bystander missense changes incur a penalty of -0.10 per affected position. Bystander nonsense changes incur -0.25. Splice donor/acceptor proximity (within 2 bp of exon boundary) incurs -0.15. Splice region proximity (3--8 bp) incurs -0.08. All-synonymous bystanders receive a +0.05 bonus. Adjusted scores are clamped to [0, 1]. Strategies are ranked by adjusted TOPSIS score, with ties broken by safety score.

### Sensitivity analysis

To quantify the robustness of strategy rankings under weight uncertainty, the system performs Monte Carlo sensitivity analysis with 10,000 Dirichlet-sampled weight permutations. For each iteration, a weight vector is drawn from a Dirichlet distribution with concentration parameter alpha = 10 * w_default (where w_default is the default weight vector), ensuring that sampled weights are centered on the defaults but vary moderately. TOPSIS is re-run with the perturbed weights, and the rank of each strategy is recorded.

The output for each strategy includes: (i) rank stability (fraction of 10,000 runs in which the strategy is top-ranked), (ii) mean rank across all runs, and (iii) the full rank distribution (fraction of runs at each rank position). This analysis is performed using Python's standard library `random.gammavariate` for Gamma sampling (the building block of Dirichlet sampling), requiring no additional dependencies beyond the standard library.

### Off-target scoring (CFD and MIT)

Off-target specificity is scored using two complementary frameworks.

The **CFD score**^15^ uses a position-specific mismatch tolerance matrix derived from Table S19 of Doench et al. (2016). For each candidate 20-mer guide, potential off-target sites are enumerated within a local genomic window by identifying sequences with up to 4 mismatches. Each mismatch is scored by position (PAM-distal position 1 is most tolerant; PAM-proximal position 20 is least tolerant) and mismatch type (rG:dT wobble pairs are better tolerated; rC:dC mismatches are poorly tolerated). The per-site CFD score is the product of individual mismatch penalties. An aggregate specificity score is computed as 1 / (1 + sum of off-target CFD scores), where higher values indicate fewer and less severe predicted off-targets.

The **MIT specificity score**^16^ uses the Hsu et al. (2013) position-weight matrix with a different aggregation scheme. Both scores are computed for each guide and reported in the output. The CFD aggregate score is used as the default input to the Risk component of TOPSIS scoring.

We note that this off-target scoring uses local sequence enumeration, not genome-wide off-target search. The CFD and MIT matrices are trained on SpCas9 data; their extension to alternative nucleases is an extrapolation (Tier B evidence). For clinical applications, genome-wide off-target prediction using dedicated tools is recommended.

### HGVS parser and ClinVar integration

The HGVS parser^18^ accepts clinical variant nomenclature in coding DNA format and converts it to genomic coordinates for pipeline input. The parser uses four regular expression patterns to match substitutions, deletions, insertions, and deletion-insertions. Transcript accessions are resolved to genomic coordinates via the Ensembl REST API. The parser handles both transcript-prefixed (NM_000267.3:c.910C>T) and gene-prefixed (NF1:c.910C>T) notation.

ClinVar batch ingestion accepts TSV files (ClinVar variant summary format) and VCF files, extracting chromosome, position, reference allele, alternate allele, gene symbol, and clinical significance. This enables analysis of user-defined variant panels or entire ClinVar gene-level exports.

### Delivery-aware post-ranking annotations

After TOPSIS ranking, CRISPRArchitect appends delivery feasibility annotations to each recommended strategy. Rather than incorporating delivery as an additional TOPSIS dimension --- which would be redundant with Safety and Complexity in approximately 90% of cases --- the system applies a two-tier annotation scheme. First, hard feasibility filters flag biologically incompatible strategy-cell type combinations (e.g., dsDNA donors in iPSCs, where p53-mediated toxicity is a documented concern^5^). Second, practical delivery recommendations are generated, including donor format selection by edit size and cell-type-specific guidance. These annotations are informed by a comprehensive literature review of 87 references covering cssDNA, lssDNA, dsDNA, and AAV donor formats. This module is implemented in `core/feasibility/delivery_advisor.py` and does not modify the TOPSIS-derived strategy rankings.

### Benchmark dataset

The benchmark comprises 30 curated variant scenarios with predefined tiered truth labels (preferred, acceptable, reject). All variants use verified GRCh38 genomic coordinates from ClinVar, with reference alleles validated against the Ensembl genome. The dataset spans 11 categories: clean base editing (n = 7), base editing negative control (n = 1), prime-editable transversions (n = 2), prime-editable small indels (n = 5), HDR-required large deletions (n = 3), HDR/PE small deletion (n = 1), compound heterozygous hybrid (n = 5), sequential HDR (n = 1), dual base editing (n = 3), and edge cases (n = 2).

Three accuracy metrics are computed: (i) top-1 accuracy: the top-ranked strategy matches a preferred or acceptable truth label; (ii) top-3 accuracy: at least one preferred or acceptable strategy appears among the top three; (iii) rejection accuracy: strategies labeled as reject are absent from the top-ranked output.

### Implementation

CRISPRArchitect v3 is implemented in Python 3.9 with dependencies limited to NumPy, SciPy, and Matplotlib. The v3 core codebase comprises approximately 12,500 lines of code across the following modules:

- **Sequence layer** (`core/sequence/`): variant normalization (`variant_normalizer.py`), transcript mapping (`transcript_mapper.py`), Ensembl client (`fetcher.py`), reference validation (`reference_validator.py`), coding annotation (`coding_annotation.py`), and HGVS parser (`hgvs_parser.py`).
- **Feasibility layer** (`core/feasibility/`): multi-nuclease base editing (`base_editing.py`), prime editing (`prime_editing.py`), HDR design (`hdr_design.py`), PAM scanning (`pam_scan.py`), off-target scoring (`off_target.py`), and delivery-aware annotations (`delivery_advisor.py`).
- **Strategy layer** (`core/mosaic/`): strategy generation (`generator.py`) and annotation integration (`annotation_integration.py`).
- **Pipeline orchestration** (`core/pipeline/`): end-to-end execution with TOPSIS scoring (`strategy_stage.py`).
- **Data models** (`core/models.py`): shared dataclasses and enums for variants, strategies, feasibility, and evidence tiers.
- **Benchmarking** (`benchmarks/`): dataset definition, evaluation, and plotting.
- **v1 modules** (~24,000 lines): ConversionSim, MOSAIC, TopoPred, ChromBridge, LoopSim, and WebApp.

All Ensembl API calls incorporate automatic retry logic with exponential backoff (3 retries with waits of 1, 2, and 4 seconds) for transient server errors (HTTP 500, 502, 503, 504, and network timeouts). The full test suite comprises 224 tests covering unit tests, integration tests, the 30-case benchmark, and delivery advisor tests.

---

## Data Availability

All benchmark datasets, including the 30-case ClinVar variant panel with verified GRCh38 coordinates and truth labels, are available at https://github.com/visvikbharti/CRISPRArchitect.

## Code Availability

Source code is available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license.

---

## Acknowledgements

We thank the CSIR-IGIB computational biology core facility for computing resources and the Ensembl team for maintaining the REST API infrastructure used for transcript mapping and reference validation.

## Author Contributions

V.B. conceived the computational framework, implemented all modules including the multi-nuclease engine, TOPSIS scorer, off-target scoring, and HGVS parser, designed and executed the benchmark evaluation, and wrote the manuscript. D.C. supervised the project, provided guidance on nuclease biology and clinical relevance, and contributed to the manuscript.

## Competing Interests

The authors declare no competing interests.

---

## References

1. Komor, A. C., Kim, Y. B., Packer, M. S., Zuris, J. A. & Liu, D. R. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature* **533**, 420--424 (2016).

2. Gaudelli, N. M. *et al.* Programmable base editing of A-T to G-C in genomic DNA without DNA cleavage. *Nature* **551**, 464--471 (2017).

3. Anzalone, A. V. *et al.* Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature* **576**, 149--157 (2019).

4. Paquet, D. *et al.* Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. *Nature* **533**, 125--129 (2016).

5. Ihry, R. J. *et al.* p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nat. Med.* **24**, 939--946 (2018).

6. Hwang, C.-L. & Yoon, K. *Multiple Attribute Decision Making: Methods and Applications*. Springer-Verlag, Berlin (1981).

7. Iyer, S. *et al.* Efficient homology-directed repair with circular single-stranded DNA donors. *CRISPR J.* **5**, 685--701 (2022).

8. Richards, S. *et al.* Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet. Med.* **17**, 405--424 (2015).

9. Arbab, M. *et al.* Determinants of base editing outcomes from target library analysis and machine learning. *Cell* **182**, 463--480.e30 (2020).

10. Rees, H. A. & Liu, D. R. Base editing: precision chemistry on the genome and transcriptome of living cells. *Nat. Rev. Genet.* **19**, 770--788 (2018).

11. Richter, M. F. *et al.* Phage-assisted evolution of an adenine base editor with improved Cas domain compatibility and activity. *Nat. Biotechnol.* **38**, 883--891 (2020).

12. Nishimasu, H. *et al.* Engineered CRISPR-Cas9 nuclease with expanded targeting space. *Science* **361**, 1259--1262 (2018).

13. Walton, R. T., Christie, K. A., Whittaker, M. N. & Kleinstiver, B. P. Unconstrained genome targeting with near-PAMless engineered CRISPR-Cas9 variants. *Science* **368**, 290--296 (2020).

14. Acharya, S. *et al.* PAM-flexible Engineered FnCas9 variants for robust and ultra-precise genome editing and diagnostics. *Nat. Commun.* **15**, 5471 (2024).

15. Doench, J. G. *et al.* Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat. Biotechnol.* **34**, 184--191 (2016).

16. Hsu, P. D. *et al.* DNA targeting specificity of RNA-guided Cas9 nucleases. *Nat. Biotechnol.* **31**, 827--832 (2013).

17. Koblan, L. W. *et al.* Improving cytidine and adenine base editors by expression optimization and ancestral reconstruction. *Nat. Biotechnol.* **36**, 843--846 (2018).

18. den Dunnen, J. T. *et al.* HGVS recommendations for the description of sequence variants: 2016 update. *Hum. Mutat.* **37**, 564--569 (2016).

19. Elliott, B., Richardson, C., Winderbaum, J., Nickoloff, J. A. & Jasin, M. Gene conversion tracts from double-strand break repair in mammalian cells. *Mol. Cell. Biol.* **18**, 93--101 (1998).

20. Concordet, J.-P. & Haeussler, M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Res.* **46**, W242--W245 (2018).

---

## Figure Legends

**Figure 1. CRISPRArchitect v3 system architecture.** Schematic of the seven-stage pipeline: variant input (genomic coordinates, HGVS notation, or ClinVar batch), transcript mapping via Ensembl, reference validation, coding and splice-site annotation, multi-nuclease feasibility evaluation across base editing (five nucleases, seven editors), prime editing, and HDR, TOPSIS scoring with consequence adjustments, and Monte Carlo sensitivity analysis. New v3 components are highlighted: multi-nuclease engine, TOPSIS scorer, off-target scoring (CFD/MIT), and HGVS parser.

**Figure 2. Multi-nuclease impact on strategy distribution.** Comparison of v2 (SpCas9 + ABE7.10 only) and v3 (multi-nuclease) strategy distributions across the 30-case benchmark. v2: BE = 0/30 (0%), PE = 29/30 (97%), HDR = 1/30 (3%). v3: BE = 6/30 (20%), PE = 23/30 (77%), HDR = 1/30 (3%). The six base editing rescues are attributable to ABE8e (broader window, positions 3--9) paired with enFnCas9 (NRG PAM) or SpCas9-NG (NG PAM).

**Figure 3. Editor-nuclease feasibility heatmap.** Heatmap showing base editing feasibility across the 30 benchmark cases for each editor-nuclease combination. Rows: editor-nuclease pairs (ABE7.10+SpCas9, ABE8e+SpCas9, ABE8e+enFnCas9, ABE8e+SpCas9-NG, ABE8e+SpRY, BE4max+SpCas9, BE4max+enFnCas9, BE4max+SpCas9-NG, BE4max+SpRY). Columns: benchmark cases. Color indicates feasibility (green = feasible with target in window; amber = PAM available but target outside window; red = no PAM; white = wrong mutation type).

**Figure 4. TOPSIS scoring and sensitivity analysis.** (A) Schematic of the TOPSIS algorithm showing normalization, weighting, ideal/anti-ideal solution determination, and relative closeness computation. (B) Example TOPSIS scores for a representative variant with base editing, prime editing, and HDR alternatives. (C) Monte Carlo sensitivity analysis output: rank distribution across 10,000 Dirichlet-sampled weight permutations for a case where prime editing has high rank stability (>0.95) versus a case where base editing and prime editing are close competitors (rank stability ~0.55).

**Figure 5. Benchmark accuracy results.** (A) Overall accuracy: top-1 = 86.7% (26/30), top-3 = 96.7% (29/30), rejection = 86.7% (26/30). (B) Accuracy by variant category: 100% top-1 for clean base editing, PE transversion, PE small indel, compound heterozygous, dual base editing, and edge cases. Large deletion cases show 0/3 top-1 but 3/3 top-3 accuracy.

**Figure 6. Per-case benchmark results.** Heatmap showing top-1 correctness, top-3 correctness, and rejection correctness for each of the 30 benchmark cases, grouped by category. Green = correct, red = incorrect, white = not applicable.

**Figure S3. SDSA displacement probability sensitivity analysis.** Rank stability of strategy recommendations across SDSA displacement probability values from p=0.001 to p=0.005 (corresponding to mean gene conversion tract lengths of 1000 bp to 200 bp). The default p=0.002 produces robust rankings: the relative ordering of strategies and the predicted cssDNA/lssDNA donor ratio remain within published experimental ranges across the full tested parameter range, confirming that the assumed displacement probability does not materially affect strategy recommendations.

**Supplementary Figure 1. Editor activity windows.** Comparison of editing windows for ABE7.10 (positions 4--7), ABE8e (positions 3--9), and BE4max (positions 4--8) on a 20-nt protospacer, illustrating how ABE8e's broader window increases the probability of placing a target nucleotide within the active region.

**Supplementary Figure 2. CFD position tolerance profile.** Position-dependent mismatch tolerance from the CFD scoring matrix (Doench et al., 2016), showing decreasing tolerance from PAM-distal (position 1, tolerance ~0.85) to PAM-proximal (position 20, tolerance ~0.05). This gradient underlies the off-target specificity scores computed for each candidate guide.
