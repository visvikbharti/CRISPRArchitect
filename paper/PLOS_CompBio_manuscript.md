# TITLE PAGE

**Full Title:** A Monte Carlo framework for predicting HDR gene conversion outcomes and optimizing multi-site genome editing strategies

**Short Title:** Computational framework for multi-site genome editing

**Authors:**

Vishal Bharti ^1,*^ and Debojyoti Chakraborty ^1,2,*^

^1^ CSIR-Institute of Genomics and Integrative Biology (IGIB), New Delhi, 110025, India

^2^ Academy of Scientific and Innovative Research (AcSIR), Ghaziabad, 201002, India

^*^ Corresponding authors

E-mail: vishalvikashbharti@gmail.com (VB); debojyoti.chakraborty@igib.in (DC)

---

# ABSTRACT

Correcting multiple pathogenic mutations in large multi-exon genes using CRISPR-based technologies requires choosing among several editing modalities, including homology-directed repair with various donor formats, base editing, and prime editing. This decision is consequential in human induced pluripotent stem cells, where double-strand breaks trigger p53-dependent apoptosis, yet it is currently guided by intuition rather than quantitative analysis. Here we present a computational framework with two components. The first, ConversionSim, models homology-directed repair as a four-step stochastic process — end resection, RAD51 filament formation, strand invasion, and synthesis-dependent strand annealing — to predict gene conversion tract length distributions and the probability of incorporating donor-encoded edits at various distances from the double-strand break. The second, MOSAIC, enumerates all feasible editing strategies for a given multi-site correction problem and ranks them by predicted efficiency, safety, time, and cost. We calibrated ConversionSim using parameters from independent published studies and validated it against four experimental datasets. The model correctly predicted the approximately two-fold enhancement of circular over linear single-stranded DNA donors (predicted 2.07-fold versus observed 1.5-2.1-fold) and the enhancement from staggered nuclease cuts (predicted 1.82-fold versus observed 1.4-2.8-fold). It also revealed a systematic discrepancy for short oligonucleotide donors, attributable to the single-strand template repair pathway not captured by our synthesis-dependent strand annealing model. We benchmarked MOSAIC against strategies employed in 14 published genome editing studies, finding 71 percent concordance with expert decisions and complete agreement for base editing applications. The framework, available as open-source software, provides the first quantitative tool for comparing editing modalities in multi-site correction scenarios.

---

# AUTHOR SUMMARY

When patients carry two or more mutations in a single gene, correcting them using genome editing is not straightforward. Researchers must decide: should each mutation be fixed separately or simultaneously? Should they use traditional cut-and-paste repair, or newer approaches like base editing that change single DNA letters without cutting both strands? These choices matter enormously in patient-derived stem cells, where DNA cuts can trigger cell death.

Currently, researchers make these decisions based on experience and intuition. We developed a computational tool that formalizes this decision. One component simulates the molecular repair process after a DNA cut to predict how far a correction template can reach — typically only a few hundred base pairs, far too short to bridge two mutations separated by the vast stretches of non-coding DNA between exons. The other component evaluates all possible editing strategies and ranks them. When we compared its recommendations against what expert researchers actually chose in 14 published studies, it agreed 71 percent of the time. The tool is freely available and requires no specialized computational expertise.

---

# INTRODUCTION

Precise correction of disease-causing mutations using CRISPR-Cas nucleases holds considerable therapeutic promise, yet the path from identifying a pathogenic variant to designing an optimal correction strategy is neither straightforward nor well-supported by existing computational tools. This challenge is particularly acute when a patient carries two or more mutations in a large multi-exon gene — a common scenario in autosomal recessive disorders where compound heterozygosity is the rule rather than the exception.

The researcher confronting a multi-site correction problem faces a combinatorial decision space. At minimum, the following choices must be made: which nuclease to employ, which DNA repair pathway to exploit, whether to use a donor template (and in what format), whether to edit both sites simultaneously or in sequential rounds, and how to balance efficiency against the considerable safety constraints imposed by certain cell types. In human induced pluripotent stem cells (iPSCs), these safety constraints are particularly stringent. Double-strand breaks (DSBs) activate the p53 tumor suppressor pathway, leading to apoptosis and, more concerningly, selective expansion of p53-deficient cells that may harbor oncogenic potential [1,2]. Each additional DSB introduced into an iPSC genome exacerbates these risks.

The expanding toolkit of editing modalities provides alternatives that partially or completely avoid DSBs. Cytosine base editors (CBEs) and adenine base editors (ABEs) can install transition mutations without cleaving both DNA strands [3,4]. Prime editors use a Cas9 nickase fused to a reverse transcriptase to write short edits at a target site, again without a DSB [5]. Circular single-stranded DNA (cssDNA) donors offer improved homology-directed repair (HDR) efficiency relative to linear donors, with reduced random integration [6,7]. Engineered nuclease variants that generate staggered rather than blunt DNA ends can shift repair outcomes toward HDR [8]. Yet despite this expanding repertoire, no existing computational tool helps researchers compare these modalities or select the most appropriate strategy for a given multi-site problem.

Current CRISPR design tools address individual steps in the editing workflow. Guide RNA design tools such as CRISPOR [9] predict on-target efficiency and off-target sites. The IDT Alt-R HDR Design Tool assists with single-site donor construction. PrimeDesign [10] and BE-Designer [11] support prime editing and base editing guide design, respectively. However, none of these tools considers the multi-site problem, none compares across editing modalities, and none models the mechanistic constraints that determine whether a given editing strategy is physically feasible.

Among these constraints, the gene conversion tract length — the distance from a DSB over which donor-encoded sequence is incorporated during HDR — is perhaps the least appreciated yet most consequential. Published measurements indicate that conversion tracts in mammalian cells are typically short, with most events incorporating donor sequence over fewer than 1,000 base pairs [12,13]. This constraint means that a single donor template cannot correct two mutations separated by large intronic distances, even if the exonic sequences are close in the messenger RNA. To our knowledge, no computational tool models this process quantitatively.

Here we present a framework that addresses these gaps. ConversionSim is a Monte Carlo simulator that models HDR as a sequential stochastic process — end resection, RAD51 nucleoprotein filament formation, strand invasion into a donor template, and synthesis-dependent strand annealing (SDSA) — to predict gene conversion tract length distributions. MOSAIC (Multi-locus Optimized Strategy for Allele-specific Integrated Correction) takes a gene structure, mutation positions, cell type, and nuclease as input, enumerates all feasible editing strategies, and ranks them by composite scores reflecting efficiency, safety, time, and cost.

We validated ConversionSim against four independent published datasets and benchmarked MOSAIC against strategies employed in 14 published genome editing studies. We report where the models succeed, where they fail, and what the failures reveal about the underlying biology.

---

# RESULTS

## ConversionSim: a stochastic model of HDR gene conversion

We modeled HDR as a sequential stochastic process comprising four mechanistic steps (Fig 1A). Each step is parameterized using values from independent published biochemical and cell biological studies (S1 Table). No parameters were fit to the validation datasets described below.

**End resection.** Following a DSB, the MRN complex and CtIP initiate short-range 5'-to-3' resection, and EXO1 or BLM-DNA2 perform long-range resection, generating 3' ssDNA overhangs [14,15]. We model resection as the sum of a normally distributed short-range component (mean 200 bp, standard deviation 80 bp) and a log-normally distributed long-range component (median approximately 1,500 bp). For nucleases that generate staggered cuts with 5' overhangs, the pre-existing overhang is added to the resected length and long-range resection is increased by 20%, reflecting evidence that 5' overhangs promote resection initiation [8].

**RAD51 filament formation.** RPA coats the resected ssDNA and is subsequently displaced by RAD51 with the assistance of BRCA2 [16]. We model filament coverage as a fraction of the resected length, drawn from a Beta distribution (mean approximately 85%). Strand invasion requires a minimum filament of approximately 15 nucleotides, based on the 8-nucleotide microhomology sampling mechanism characterized in single-molecule studies [17].

**Strand invasion.** We model the probability of productive strand invasion as a function of cell-type-specific HDR competence, donor topology (circular ssDNA, linear ssDNA, linear dsDNA), nuclease-induced stagger, cell cycle phase, and homology arm length.

**SDSA synthesis.** Following strand invasion, DNA polymerase extends the invading strand using the donor as template. In mitotic cells, SDSA is the dominant sub-pathway, in which the extended strand is displaced and reanneals to the other resected end [18]. We model synthesis as a geometric process: at each base pair of extension, there is a probability *p* of D-loop collapse and strand displacement. This produces a right-skewed tract length distribution consistent with the observation that most conversion tracts are short [12]. For circular ssDNA donors, *p* is reduced by 20% (reflecting greater donor stability), and for staggered cuts, *p* is reduced by 15% (reflecting longer resection-mediated invasion stability).

For a configuration typical of iPSC editing with a staggered-cut nuclease and cssDNA donor bearing 300 bp homology arms, ConversionSim predicts a median gene conversion tract of 542 bp (interquartile range: 230–950 bp) and an HDR success rate of approximately 3.8% of cells (50,000 simulations; Fig 1B). The probability of incorporating a donor-encoded edit decreases with distance from the DSB: approximately 89% at 100 bp, 53% at 500 bp, 28% at 1,000 bp, and 7% at 2,000 bp (Fig 1C).

## Validation of ConversionSim

We compared ConversionSim predictions against four published datasets that were not used for parameter calibration (Fig 2, Table 1).

**Circular versus linear ssDNA donor enhancement.** Iyer et al. reported that cssDNA donors produced approximately 1.9-fold higher HDR frequency than equimolar linear ssDNA donors across multiple loci and cell types (range 1.5–2.1-fold) [6]. ConversionSim predicts a ratio of 2.07-fold. This value falls within the experimentally observed range (Fig 2A).

**Staggered-cut enhancement.** Chauhan et al. showed that engineered Cas9 variants producing 5' overhangs of six or more base pairs achieved a mean 1.9-fold increase in precise editing frequency (range 1.4–2.8-fold) compared to wild-type Cas9 [8]. ConversionSim predicts a ratio of 1.82-fold for six-base-pair staggered versus blunt cuts, within the reported range (Fig 2B).

**Tract length distribution shape.** Elliott et al. measured gene conversion tracts in mouse embryonic stem cells, finding a strongly right-skewed distribution with a maximum tract of 511 bp [12]. ConversionSim produces a qualitatively similar right-skewed distribution but predicts longer tracts (median approximately 350 bp versus the majority of tracts below 58 bp in Elliott et al.). We attribute this quantitative difference to the distinct experimental systems: Elliott et al. used an endogenous chromosomal donor, whereas ConversionSim models exogenous donor-mediated repair, for which longer SDSA tracts have been reported in subsequent studies (Fig 2C) [19].

**Distance-dependent incorporation with ssODN donors.** Paquet et al. demonstrated that mutation incorporation frequency decreases sharply with distance from the Cas9 cut site when using short ssODN donors, with edits beyond approximately 30-50 bp from the cut incorporated at markedly reduced frequency [13]. ConversionSim shows poor agreement with this dataset, systematically over-predicting incorporation at all distances (Fig 2D).

This discrepancy is, we believe, informative rather than simply a model failure. Short ssODN donors (less than approximately 200 nucleotides) are incorporated primarily through single-strand template repair (SSTR), a RAD51-independent pathway with very short conversion tracts [20]. ConversionSim models the SDSA pathway, which dominates when longer donor templates — cssDNA, long ssDNA, dsDNA, or AAV vectors — are used. The poor fit to the Paquet et al. data thus delineates a meaningful biological boundary: ConversionSim is appropriate for long-donor HDR but not for ssODN-mediated editing.

## MOSAIC: strategy optimization for multi-site editing

MOSAIC addresses the question: given two or more mutations in a gene, which editing strategy should be used? MOSAIC takes as input a gene structure, mutation positions and types, cell type, and nuclease, and enumerates up to eight strategy classes (Fig 3A):

(1) Single-template HDR, feasible only when mutations are separated by less than approximately 5,000 bp of genomic distance; (2) dual-template simultaneous HDR; (3) sequential HDR in separate rounds; (4) dual base editing, when both mutations are transitions amenable to CBE or ABE; (5) dual prime editing; (6) hybrid approaches combining DSB-free editing at one site with HDR at the other; (7) NHEJ-mediated exon deletion; and (8) sequential mixed-modality editing.

Each strategy is scored along four dimensions: predicted efficiency, safety, time, and cost. The safety dimension penalizes DSBs, with particular severity in cell types with active p53 pathways. Two simultaneous DSBs incur an additional penalty reflecting translocation risk, estimated from the power-law relationship between genomic distance and rearrangement frequency [21]. The four dimension scores are combined into a weighted composite (default weights: safety 40%, efficiency 30%, time 15%, cost 15%).

We note that these weights are not derived from data but reflect a judgment that safety is the primary concern in iPSC applications. They are adjustable by the user, and we do not claim that the default weights are universally optimal.

## Benchmarking MOSAIC against published studies

We benchmarked MOSAIC against editing strategies employed in 14 published studies spanning 8 genes and 6 distinct editing modalities (Table 2).

**Overall concordance.** MOSAIC's top-three ranked strategies included the study authors' chosen strategy in 10 of 14 cases (71.4%). When considering only the top-ranked recommendation, concordance was 7 of 14 (50%).

**Concordance by modality.** Agreement was highest for base editing studies (5 of 5, 100%) and prime editing (1 of 1, 100%). For HDR-based studies, concordance was 3 of 5 (60%).

**Analysis of discordant cases.** The four discordant cases share a common pattern: MOSAIC recommended a DSB-free approach (base editing or prime editing) while the study authors used HDR or NHEJ-based exon deletion. In two of these cases, the studies were published before the relevant base or prime editors were available for the target mutation type, suggesting that MOSAIC's recommendation may be prospectively appropriate even when retrospectively discordant. In one case, the mutation was a large Alu insertion that we could not accurately model in our benchmark. In one case, MOSAIC lacks a "single-cut exon reframing" strategy type that would match the approach used.

We emphasize that concordance with published studies is an imperfect metric. A published strategy reflects the tools available at the time, laboratory expertise and preferences, and practical considerations that MOSAIC cannot model. Conversely, MOSAIC may recommend strategies that the authors would have preferred had they been available.

## Application to clinically relevant human genes

We applied the framework to three genes using exon coordinates from Ensembl (GRCh38 assembly): NF1 (58 exons, neurofibromatosis type 1), DMD (79 exons, Duchenne muscular dystrophy), and BRCA2 (27 exons, hereditary breast and ovarian cancer).

For hypothetical dual-mutation scenarios in each gene, CRISPRArchitect consistently found that single-template dual-site HDR is infeasible when the two mutation sites are separated by tens of kilobases or more of genomic distance. Even for BRCA2, where the tested mutations were separated by only 47 kb, the predicted three-dimensional nuclear distance between the two loci (846 nm) far exceeds the random-coil diameter of a 3 kb cssDNA donor (261 nm), precluding simultaneous engagement of both target sites by a single donor molecule.

---

# DISCUSSION

We have presented a computational framework for predicting HDR gene conversion outcomes and optimizing multi-site genome editing strategies. Here we discuss what the framework contributes, what it does not, and what would be needed to strengthen it.

## Contributions

The primary contributions are twofold. First, ConversionSim provides, to our knowledge, the first publicly available quantitative model of the SDSA-mediated gene conversion process, parameterized from published biochemical data and validated against independent experimental datasets. The model correctly predicts the relative enhancements conferred by circular donor topology and staggered nuclease cuts, and produces tract length distributions consistent with the expected SDSA behavior.

Second, MOSAIC formalizes a decision that is currently made informally: when multiple mutations require correction, which combination of editing modalities is most appropriate? By enumerating all feasible strategies and scoring them on multiple dimensions, MOSAIC ensures that researchers consider options they might otherwise overlook — particularly DSB-free approaches whose advantages in iPSCs are substantial but may not be immediately apparent.

## What we do not claim

We do not claim that ConversionSim predicts absolute HDR efficiency at specific genomic loci. The model uses cell-type averages and does not account for the well-documented locus-to-locus variation driven by chromatin state, transcription, and replication timing [22]. The model provides relative comparisons (circular versus linear donors, staggered versus blunt cuts) with greater confidence than absolute predictions.

We do not claim that MOSAIC's recommendations are always superior to expert judgment. A 71% concordance rate means that in nearly a third of cases, the model and published experts disagree. We have analyzed these disagreements and found them to be largely attributable to the temporal availability of editing tools and to missing strategy types in our framework, rather than to fundamental scoring errors. However, we cannot rule out that MOSAIC's safety-first weighting may be overly conservative for some applications.

We do not claim novelty for the supporting modules (ChromBridge, cssDNA-TopoPred, LoopSim), which apply established methods — polymer physics, secondary structure prediction, and loop extrusion simulation — to the genome editing context. These modules provide useful supporting analyses but do not represent methodological advances.

## Limitations

Several limitations should be acknowledged. (1) All model parameters carry uncertainty that is not formally propagated. Future versions should incorporate Bayesian parameter estimation. (2) The SSTR pathway, used by short ssODN donors, is not modeled; extending ConversionSim to include SSTR would substantially broaden its applicability. (3) MOSAIC does not check PAM site availability, which can render a theoretically optimal strategy impractical. (4) The strategy scoring weights are heuristic. (5) We have not performed prospective experimental validation, in which editing experiments would be designed based on CRISPRArchitect's recommendations and outcomes compared to control strategies.

## Relation to existing work

No existing tool addresses the multi-site editing strategy optimization problem. The closest conceptual precedent is the manual decision trees described in review articles on iPSC genome editing, which informally weigh similar considerations to MOSAIC's scoring system. ConversionSim has no direct precedent; while analytical models of gene conversion exist in the population genetics literature, Monte Carlo simulators of the cell biological process that are parameterized for genome editing applications have not been described.

## Future directions

Three extensions would most meaningfully strengthen the framework: incorporating locus-specific chromatin features from public datasets (ENCODE, ATAC-seq) to adjust HDR efficiency predictions; adding an SSTR sub-model for ssODN donors; and prospective experimental validation of MOSAIC's strategy recommendations in iPSCs.

---

# MATERIALS AND METHODS

## ConversionSim: model specification

ConversionSim is implemented in Python using NumPy for vectorized computation. Each of *n* simulations (default: 10,000) comprises four stochastic steps executed as parallel array operations.

**End resection.** For each simulation *i*, left-side and right-side resection lengths are computed as:

*R_short* ~ N(μ = 200, σ = 80), clipped to [50, 500] bp

*R_long* ~ LogNormal(μ = ln(1500), σ = 0.8), clipped to [300, 10000] bp

*R_total* = *R_short* + *R_long*

For staggered cuts with overhang *h* bp: the overhang length is added to the left-side total, and *R_long* is scaled by 1.2 on that side.

**RAD51 filament formation.** Coverage fraction *f* ~ Beta(α = 8, β = 2), rescaled to [0.70, 0.95]. Filament length (nt) = *R_total* × *f*. Simulations with filament length less than 15 nt are marked as invasion-incompetent.

**Strand invasion.** Invasion probability *P_inv* is computed as:

*P_inv* = *η_cell* × *m_topology* × (1 + *k* × *h*) × *f_S/G2* × *φ(L_arm, L_fil)*

where *η_cell* is the cell-type base HDR rate, *m_topology* is a donor topology multiplier, *k* is the HDR enhancement per base pair of overhang, *h* is the overhang length, *f_S/G2* is the S/G2 cell cycle fraction, and *φ* is a homology-dependent function that decays exponentially when the filament length exceeds the homology arm length.

**SDSA synthesis.** Tract length *T* ~ Geometric(*p_eff*), where:

*p_eff* = *p_base* × *m_topo* × *m_stagger*

*p_base* = 0.002 per bp (producing baseline median tracts of approximately 350 bp before donor topology and stagger modifiers are applied; with cssDNA and staggered cuts, the effective median increases to approximately 540 bp). For circular ssDNA donors, *m_topo* = 0.80; for staggered cuts, *m_stagger* = 0.85. Tracts are clipped to [50, 5000] bp.

All parameter values and their literature sources are provided in S1 Table.

## MOSAIC: strategy enumeration and scoring

**Mutation classification.** Each mutation is classified by its amenability to different editing modalities. Single-nucleotide transitions (A↔G, C↔T) are flagged as amenable to adenine or cytosine base editing. All substitutions and small insertions/deletions less than 50 bp are flagged as amenable to prime editing. Larger mutations or those requiring donor-mediated correction are flagged as requiring HDR.

**Strategy enumeration.** MOSAIC generates all strategy classes whose prerequisites are met (mutation type compatibility, genomic distance constraints). Single-template HDR is offered only when the inter-mutation genomic distance is less than 5,000 bp.

**Scoring.** Each strategy receives scores on four dimensions (range 0-1):

*Efficiency*: square root transform of the predicted combined editing rate, where combined rate is the product of per-site rates for independent events.

*Safety*: starts at 1.0, with penalties for each DSB (iPSC penalty: 0.25 per DSB), an additional penalty for simultaneous dual DSBs (translocation risk), and a penalty for p53-active cell types.

*Time*: inverse of the number of editing rounds, normalized.

*Cost*: inverse of a composite of donors needed, screening burden, and selection requirements.

*Overall* = *w_s* × Safety + *w_e* × Efficiency + *w_t* × Time + *w_c* × Cost

Default weights: *w_s* = 0.40, *w_e* = 0.30, *w_t* = 0.15, *w_c* = 0.15.

## Validation and benchmarking

**ConversionSim validation.** We compared model predictions with four published datasets: Elliott et al. [12] (tract lengths), Paquet et al. [13] (distance-dependent incorporation), Iyer et al. [6] (cssDNA enhancement), and Chauhan et al. [8] (staggered-cut enhancement). Model parameters were set a priori from independent literature sources; no fitting to these validation datasets was performed.

**MOSAIC benchmarking.** We identified 14 published studies in which mutations were corrected in iPSCs or related cell types. For each study, we reconstructed the gene structure using approximate exon counts and inter-exon distances, defined mutations matching those corrected, and ran MOSAIC with the corresponding cell type and nuclease. We recorded whether the authors' chosen strategy appeared in MOSAIC's top-three ranked strategies.

## Software availability

CRISPRArchitect (version 0.1.0) is implemented in Python (compatible with versions 3.9 through 3.12) and depends only on NumPy, SciPy, and Matplotlib. An interactive web interface is provided using Streamlit. The complete source code, documentation, validation scripts, and a Docker container for reproducible deployment are available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license.

---

# ACKNOWLEDGMENTS

[To be added]

---

# REFERENCES

1. Ihry RJ, Worringer KA, Salick MR, Frias E, Ho D, Theriault K, et al. p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. Nat Med. 2018;24:939-946.

2. Haapaniemi E, Botla S, Persson J, Schmierer B, Taipale J. CRISPR-Cas9 genome editing induces a p53-mediated DNA damage response. Nat Med. 2018;24:927-930.

3. Komor AC, Kim YB, Packer MS, Zuris JA, Liu DR. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. Nature. 2016;533:420-424.

4. Gaudelli NM, Komor AC, Rees HA, Packer MS, Badran AH, Bryson DI, et al. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage. Nature. 2017;551:464-471.

5. Anzalone AV, Randolph PB, Davis JR, Sousa AA, Koblan LW, Levy JM, et al. Search-and-replace genome editing without double-strand breaks or donor DNA. Nature. 2019;576:149-157.

6. Iyer S, Mir A, Vega-Badillo J, Roscoe BP, Ibraheim R, Zhu LJ, et al. Efficient homology-directed repair with circular single-stranded DNA donors. CRISPR J. 2022;5:685-701.

7. Xie K, Starzyk J, Majumdar I, Xiao K, Pham T, Wang YH, et al. Efficient non-viral immune cell engineering using circular single-stranded DNA-mediated genomic integration. Nat Biotechnol. 2025;43:1821-1832.

8. Chauhan VP, Sharp PA, Langer R. Altered DNA repair pathway engagement by engineered CRISPR-Cas9 nucleases. Proc Natl Acad Sci USA. 2023;120:e2300605120.

9. Concordet JP, Haeussler M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. Nucleic Acids Res. 2018;46:W242-W245.

10. Hsu JY, Grünewald J, Szalay R, Shih J, Anzalone AV, Lam KC, et al. PrimeDesign software for rapid and simplified design of prime editing guide RNAs. Nat Commun. 2021;12:1034.

11. Hwang GH, Park J, Lim K, Kim S, Yu J, Kim ST, et al. Web-based design and analysis tools for CRISPR base editing. BMC Bioinformatics. 2018;19:542.

12. Elliott B, Richardson C, Winderbaum J, Nickoloff JA, Jasin M. Gene conversion tracts from double-strand break repair in mammalian cells. Mol Cell Biol. 1998;18:93-101.

13. Paquet D, Kwart D, Chen A, Sproul A, Jacob S, Teo S, et al. Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. Nature. 2016;533:125-129.

14. Symington LS. Mechanism and regulation of DNA end resection in eukaryotes. Crit Rev Biochem Mol Biol. 2016;51:195-212.

15. Cejka P. DNA end resection: nucleases team up with the right partners to initiate homologous recombination. J Biol Chem. 2015;290:22931-22938.

16. Jensen RB, Carreira A, Kowalczykowski SC. Purified human BRCA2 stimulates RAD51-mediated recombination. Nature. 2010;467:678-683.

17. Qi Z, Redding S, Lee JY, Gibb B, Kwon Y, Niu H, et al. DNA sequence alignment by microhomology sampling during homologous recombination. Cell. 2015;160:856-869.

18. Renkawitz J, Lademann CA, Jentsch S. Mechanisms and principles of homology search during recombination. Nat Rev Mol Cell Biol. 2014;15:369-383.

19. Kan Y, Ruis B, Takasugi T, Hendrickson EA. Mechanisms of precise genome editing using oligonucleotide donors. Genome Res. 2017;27:1099-1111.

20. Gallagher DN, Haber JE. Repair of a site-specific DNA cleavage: old-school lessons for Cas9-mediated gene editing. ACS Chem Biol. 2018;13:397-405.

21. Frock RL, Hu J, Meyers RM, Ho YJ, Kii E, Alt FW. Genome-wide detection of DNA double-stranded breaks induced by engineered nucleases. Nat Biotechnol. 2015;33:179-186.

22. Aymard F, Bugler B, Schmidt CK, Guillou E, Caron P, et al. Transcriptionally active chromatin recruits homologous recombination at DNA double-strand breaks. Nat Struct Mol Biol. 2014;21:366-374.

---

# SUPPORTING INFORMATION

**S1 Table. ConversionSim parameter values and literature sources.** Complete listing of all model parameters, default values, distributions, and the published studies from which they were derived.

**S1 Fig. ConversionSim sensitivity analysis.** Effect of varying individual parameters on predicted tract length distributions, illustrating which parameters most strongly influence model output.

**S2 Fig. MOSAIC scoring weight sensitivity.** Effect of varying the four dimension weights on strategy rankings, showing robustness of the top-ranked strategy across a range of weighting schemes.

**S1 Text. Detailed description of supporting modules.** ChromBridge polymer model, cssDNA-TopoPred secondary structure analysis, and LoopSim cohesin simulation — including their methodological basis, implementation, and limitations.

**S2 Text. Complete MOSAIC benchmarking results.** Detailed case-by-case analysis of all 14 benchmarked studies, including discordance analysis.

---

# FIGURE LEGENDS

**Fig 1. ConversionSim model and predictions.** (A) Schematic of the four-step stochastic HDR model: end resection generates 3' ssDNA overhangs; RAD51 filament assembles on resected ssDNA; the filament invades the donor template; SDSA-mediated synthesis copies donor sequence before D-loop collapse. (B) Predicted gene conversion tract length distribution (50,000 simulations; staggered-cut nuclease, cssDNA donor, iPSC cell type). Vertical lines indicate median (542 bp) and 95th percentile (2,258 bp). (C) Probability that a donor-encoded edit is incorporated as a function of distance from the DSB. Key distances annotated.

**Fig 2. Validation against published datasets.** (A) cssDNA versus lssDNA HDR enhancement ratio: ConversionSim predicts 2.07-fold; Iyer et al. [6] observed 1.9-fold (range 1.5-2.1). (B) Staggered-cut enhancement: ConversionSim predicts 1.82-fold; Chauhan et al. [8] observed 1.9-fold (range 1.4-2.8). (C) Tract length distribution shape compared with Elliott et al. [12]. (D) Distance-dependent incorporation: ConversionSim (SDSA model) versus Paquet et al. [13] ssODN data (SSTR pathway), illustrating the pathway-specific boundary of model applicability.

**Fig 3. MOSAIC benchmarking.** (A) Schematic of strategy enumeration for a hypothetical dual-mutation scenario. (B) Concordance between MOSAIC top-three recommendations and published study strategies, by editing modality.

---

# AUTHOR CONTRIBUTIONS

Conceptualization: VB. Software: VB. Validation: VB. Formal analysis: VB. Writing – original draft: VB. Writing – review & editing: VB, DC. Supervision: DC. Funding acquisition: DC.

---

# FUNDING

This work was supported by [CSIR grant number(s) to be added by authors]. VB acknowledges [fellowship details to be added]. The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

---

# COMPETING INTERESTS

The authors declare no competing interests.

---

# DATA AVAILABILITY STATEMENT

All source code, validation scripts, parameter files, and analysis pipelines are available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license. The repository includes a Docker container for reproducible deployment. All gene coordinates were obtained from the Ensembl REST API (GRCh38 assembly). No novel experimental data were generated in this study; all validation comparisons use data from previously published studies cited in the manuscript.

---

*Manuscript prepared for PLOS Computational Biology, Research Article format.*
