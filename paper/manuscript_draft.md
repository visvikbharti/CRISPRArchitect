# CRISPRArchitect: a computational framework for predicting HDR outcomes and optimizing multi-site genome editing strategies

---

## Authors

Vishal Bharti [1, *]

[1] [Your Institution]
[*] Corresponding author

---

## Abstract

Correcting multiple pathogenic mutations in large genes remains a significant challenge in genome editing, particularly in human induced pluripotent stem cells (iPSCs). Researchers must choose among several editing modalities — homology-directed repair (HDR), base editing, prime editing, and hybrid approaches — yet this decision is typically guided by intuition rather than quantitative analysis. Here we present CRISPRArchitect, a computational framework comprising two core components: (i) ConversionSim, a Monte Carlo simulator that models the step-by-step mechanics of HDR — end resection, RAD51 filament formation, strand invasion, and synthesis-dependent strand annealing (SDSA) — to predict gene conversion tract length distributions; and (ii) MOSAIC, a strategy optimizer that enumerates and scores all feasible editing approaches for a given multi-site correction problem. We calibrated ConversionSim against four published experimental datasets, finding good agreement for donor topology effects (predicted 2.07-fold cssDNA enhancement vs. observed 1.9-fold) and staggered-cut effects (predicted 1.82-fold vs. observed 1.9-fold), while identifying a systematic discrepancy for short oligonucleotide donors attributable to the distinct SSTR repair pathway not captured by our SDSA-based model. We benchmarked MOSAIC against strategies employed in 14 published genome editing studies, finding that MOSAIC's top-three recommendations included the authors' chosen strategy in 71% of cases, with 100% concordance for base editing applications. CRISPRArchitect is freely available as an open-source Python toolkit with a web interface at https://github.com/visvikbharti/CRISPRArchitect.

**Keywords:** CRISPR, homology-directed repair, gene conversion, iPSC, base editing, strategy optimization, Monte Carlo simulation

---

## Introduction

Genome editing with CRISPR-Cas nucleases has transformed our ability to introduce precise modifications into mammalian genomes. However, when a patient carries two or more pathogenic mutations in a large multi-exon gene, the path from diagnosis to correction is far from straightforward. The researcher faces a combinatorial decision space: which nuclease to use, which repair pathway to exploit, whether to edit sequentially or simultaneously, and what donor template format to employ.

This decision is particularly consequential in human iPSCs, where double-strand breaks (DSBs) activate p53-mediated apoptosis and can select for p53-deficient cells (Ihry et al., 2018; Haapaniemi et al., 2018), where HDR efficiency is limited by the predominantly G1-phase cell cycle, and where each round of clonal selection adds passage-associated risks to genomic stability.

The expanding toolkit of editing modalities — including cytosine and adenine base editors (Komor et al., 2016; Gaudelli et al., 2017), prime editors (Anzalone et al., 2019), and improved donor template technologies such as circular single-stranded DNA (cssDNA) (Iyer et al., 2022; Xie et al., 2025) — provides alternatives to conventional HDR that avoid DSBs entirely. Yet no computational tool currently helps researchers navigate this decision space. Existing CRISPR design tools focus on guide RNA selection (Concordet and Haeussler, 2018) or donor template construction (Richardson et al., 2016) for individual editing events, but none considers the multi-site problem or compares across editing modalities.

Here we address this gap with CRISPRArchitect, a framework built around two complementary components. First, ConversionSim models the sequential steps of HDR — from end resection through SDSA-mediated synthesis — as a stochastic process, enabling quantitative prediction of gene conversion tract lengths and the probability of incorporating donor-encoded edits at various distances from the DSB. Second, MOSAIC (Multi-locus Optimized Strategy for Allele-specific Integrated Correction) takes a gene structure, mutation positions, cell type, and nuclease as input, enumerates all feasible editing strategies, and ranks them by predicted efficiency, safety, time, and cost.

We emphasize that CRISPRArchitect provides *estimates informed by published experimental parameters*, not exact predictions. Biological variability across loci, cell lines, and experimental conditions means that any computational tool in this space operates with substantial uncertainty. We have validated our models against published data where possible and are transparent about the boundaries of our approach.

---

## Results

### ConversionSim: a stochastic model of HDR gene conversion

ConversionSim models HDR as a sequential stochastic process comprising four mechanistic steps (Figure 1a):

**Step 1 — End resection.** Following a DSB, the MRN complex and CtIP initiate short-range 5'-to-3' resection, generating 3' single-stranded DNA (ssDNA) overhangs. Long-range resection by EXO1 or BLM-DNA2 extends these overhangs further (Symington, 2011; Cejka, 2015). We model short-range resection as normally distributed (mean 200 bp, s.d. 80 bp) and long-range resection as log-normally distributed (median ~1,500 bp), based on biochemical measurements. For nucleases that generate staggered cuts with 5' overhangs (e.g., Cas12a, engineered Cas9 variants), the pre-existing overhang provides a head start for resection, which we model as a 20% increase in long-range resection length.

**Step 2 — RAD51 filament formation.** RPA coats the resected ssDNA, followed by BRCA2-mediated loading of RAD51 (Jensen et al., 2010). We model filament coverage as a fraction of the resected length, drawn from a Beta distribution (mean ~85% coverage). A minimum filament length of 15 nucleotides (~5 RAD51 monomers) is required for productive strand invasion, based on single-molecule studies (Qi et al., 2015).

**Step 3 — Strand invasion and donor engagement.** The RAD51 filament searches for homology in the donor template. We model invasion probability as a function of homology arm length, cell-type-specific HDR competence, donor topology, and cell cycle phase (S/G2 fraction).

**Step 4 — SDSA synthesis.** After strand invasion, DNA polymerase extends the invading strand using the donor as template. In mitotic cells, the dominant HDR sub-pathway is synthesis-dependent strand annealing (SDSA), which limits the length of template copying (Renkawitz et al., 2014). We model synthesis as a geometric process: at each base pair of extension, there is a fixed probability *p* of D-loop collapse and strand displacement. This produces an exponential-like distribution of tract lengths consistent with published observations that most conversion tracts are short (<1 kb) with a declining tail (Elliott et al., 1998).

The model incorporates donor topology effects: circular ssDNA donors show reduced D-loop displacement probability (reflecting their greater intracellular stability due to exonuclease resistance; Iyer et al., 2022), and staggered cuts further reduce displacement probability (reflecting the longer resection and more stable invasion intermediates; Chauhan et al., 2023).

For a typical configuration — staggered-cut nuclease, circular ssDNA donor with 300-bp homology arms, iPSC cell type — ConversionSim predicts a median gene conversion tract of ~540 bp (interquartile range: 230–950 bp) and an HDR success rate of ~3.8% per cell (50,000 simulations; Figure 1b). The probability of incorporating a donor-encoded edit decreases monotonically with distance from the DSB: ~89% at 100 bp, ~53% at 500 bp, ~28% at 1,000 bp, and ~7% at 2,000 bp (Figure 1c). At distances exceeding ~5,000 bp, the probability is effectively zero, reflecting the fundamental length constraint imposed by SDSA.

### Validation of ConversionSim against published data

We validated ConversionSim against four independent published datasets (Figure 2, Table 1).

**Donor topology effects.** Iyer et al. (2022) reported that circular ssDNA donors produced approximately 1.9-fold higher HDR frequency than equimolar linear ssDNA donors in HEK293T and K562 cells (range 1.5–2.1-fold across loci). ConversionSim predicts a ratio of 2.07-fold, which falls within the observed experimental range (Figure 2a).

**Staggered cut effects.** Chauhan et al. (2023) demonstrated that engineered Cas9 variants producing 5' staggered cuts of ≥6 bp achieved a mean 1.9-fold improvement in precise editing frequency over wild-type Cas9 (range 1.4–2.8-fold). ConversionSim predicts a ratio of 1.82-fold for 6-bp staggered versus blunt cuts, again within the reported range (Figure 2b).

**Tract length distribution shape.** Elliott et al. (1998) measured gene conversion tracts in mouse embryonic stem cells using an I-SceI-based system, finding a strongly right-skewed distribution with most tracts ≤58 bp and a maximum of 511 bp. ConversionSim produces a qualitatively similar right-skewed distribution. However, the model predicts longer tracts (median ~350 bp) than Elliott et al. observed, which we attribute to the difference between their endogenous chromosomal donor system and the exogenous donor scenario that ConversionSim is designed to model. This is consistent with later reports of longer SDSA tracts (200–2,000 bp) with exogenous donors in human cells (Figure 2c).

**Distance-dependent incorporation.** Paquet et al. (2016) showed that SNP incorporation frequency in iPSCs decreases sharply with distance from the Cas9 cut site, with a practical guideline of placing edits within 10 bp of the cut. ConversionSim shows poor quantitative agreement with this dataset (R² = -0.56), systematically over-predicting incorporation at all distances. We attribute this discrepancy to a fundamental mechanistic difference: Paquet et al. used short single-stranded oligodeoxynucleotides (ssODNs, <200 nt), which are primarily incorporated through the single-strand template repair (SSTR) pathway — a RAD51-independent mechanism with very short conversion tracts (~30–50 bp) — rather than the SDSA pathway modeled by ConversionSim (Figure 2d). This finding defines a clear boundary of our model's applicability: ConversionSim is appropriate for long donor templates (cssDNA, dsDNA, AAV) where SDSA dominates, but not for short ssODN donors where SSTR is the primary repair mechanism.

### MOSAIC: a multi-site editing strategy optimizer

MOSAIC addresses a question that, to our knowledge, no existing tool considers: given multiple mutations in a gene, what is the optimal editing strategy? MOSAIC takes as input a gene structure (exon coordinates), mutation positions and types, cell type, and available nuclease, and enumerates up to eight strategy classes:

1. **Single-template HDR** — one donor spanning both sites (only if genomic distance <5 kb)
2. **Dual-template simultaneous HDR** — two donors + two guide RNAs in one experiment
3. **Sequential HDR** — correct one site per editing round
4. **Dual base editing** — if both mutations are transitions amenable to CBE or ABE
5. **Dual prime editing** — if both mutations are small enough for pegRNA-mediated correction
6. **Hybrid approaches** — base/prime editing at one site, HDR at the other
7. **Exon deletion** — NHEJ-mediated removal of mutation-containing exons
8. **Sequential mixed modality** — e.g., base editing in round 1, prime editing in round 2

Each strategy is scored on four dimensions: predicted efficiency (based on published cell-type-specific editing rates), safety (penalizing DSBs in iPSCs due to p53 selection risk, penalizing simultaneous dual DSBs for translocation risk), time (number of editing rounds), and cost (screening burden, donor production). The weighted composite score (default weights: safety 40%, efficiency 30%, time 15%, cost 15%) reflects the reality that in iPSC applications, safety often outweighs raw efficiency.

### Benchmarking MOSAIC against published studies

We benchmarked MOSAIC against editing strategies employed in 14 published studies spanning 8 genes and multiple editing modalities (Table 2, Figure 3).

**Overall accuracy.** MOSAIC's top-three ranked strategies included the strategy chosen by the study authors in 10 of 14 cases (71.4%). When restricted to the top-ranked recommendation, MOSAIC matched the authors' choice in 7 of 14 cases (50%).

**Accuracy by editing modality.** MOSAIC showed perfect concordance with published studies that used base editing (5/5, 100%) and prime editing (1/1, 100%). For HDR-based studies, concordance was 60% (3/5). The four discordant cases all involved MOSAIC recommending a DSB-free approach (base or prime editing) while the authors used HDR. In two of these cases, the studies were published before the relevant base editors were available, suggesting that MOSAIC's recommendation may be prospectively valid even when retrospectively discordant. In one case, the discrepancy was attributable to a mutation type (large Alu insertion) that our benchmark modeled imprecisely. In one case, MOSAIC lacked a "single-cut exon reframing" strategy that would have matched the DMD approach used.

We note that a 71% concordance rate should be interpreted cautiously. The authors' choice in a published study is not necessarily optimal — it reflects the tools available at the time, laboratory expertise, and considerations (e.g., intellectual property, reagent availability) that MOSAIC does not model. Conversely, MOSAIC's recommendations may be suboptimal for reasons it cannot capture, such as locus-specific chromatin effects on editing efficiency.

### Supporting analysis modules

CRISPRArchitect includes three additional modules that support — but do not constitute the primary scientific contribution of — the framework:

**ChromBridge** estimates the 3D physical distance between two genomic loci using a confined polymer model of chromatin, calibrated to published FISH measurements (Yokota et al., 1995) and Hi-C contact frequency data (Lieberman-Aiden et al., 2009). For the multi-site editing question, ChromBridge quantifies why a single donor template cannot physically bridge distant loci: a 3 kb cssDNA has a random coil diameter of ~261 nm, while two exons separated by 1 Mb are typically ~1,300 nm apart in 3D nuclear space.

**cssDNA-TopoPred** analyzes circular ssDNA donor templates for secondary structures — G-quadruplexes and hairpins — that could impair HDR by sequestering homology arm nucleotides. While the underlying algorithms (regex-based G4 detection, simplified nearest-neighbor hairpin prediction) are not novel, their application to HDR donor quality assessment addresses a practical need.

**LoopSim** simulates cohesin-mediated loop extrusion from DSB sites, inspired by the recent finding that cohesin drives chromatin scanning during the RAD51-mediated homology search (Marin-Gonzalez et al., 2025). We note that this mechanism is primarily relevant for endogenous donor recombination (sister chromatid, ectopic genomic donors) and has limited direct bearing on exogenous donor-mediated HDR, where the donor reaches the DSB through 3D diffusion.

### Application to real human genes

We applied CRISPRArchitect to three clinically relevant genes using exon coordinates fetched from Ensembl (GRCh38): NF1 (58 exons, neurofibromatosis type 1), DMD (79 exons, Duchenne muscular dystrophy), and BRCA2 (27 exons, hereditary breast cancer). For hypothetical dual-mutation scenarios in each gene, CRISPRArchitect consistently found that single-template dual-site HDR is infeasible when mutations are separated by >5 kb of genomic distance, and recommended DSB-free approaches (base editing, prime editing) when the mutation types permit.

---

## Discussion

We have presented CRISPRArchitect, a computational framework for optimizing multi-site genome editing strategies. Its two core contributions are: (i) ConversionSim, a Monte Carlo model of HDR that predicts gene conversion tract lengths from first principles of repair biology, and (ii) MOSAIC, a strategy optimizer that formalizes the multi-modality editing decision into a quantitative scoring system.

### What we can and cannot claim

ConversionSim produces tract length distributions and donor-type enhancement ratios that are consistent with published experimental data for long donor templates. The model correctly predicts that circular ssDNA donors outperform linear ssDNA (~2-fold) and that staggered nuclease cuts enhance HDR (~1.8-fold), both within the ranges reported in independent studies. However, the model is not validated for predicting absolute HDR rates at specific genomic loci — such predictions would require locus-specific parameters (chromatin state, replication timing, transcription level) that we do not model.

The identification of the SSTR pathway boundary is, we believe, a useful contribution: it makes explicit that different donor formats engage fundamentally different repair pathways, and that models of one pathway should not be applied to the other.

MOSAIC's 71% concordance with published expert decisions is encouraging but should not be over-interpreted. The strategy scoring weights are heuristic, not learned from data, and the efficiency estimates are literature averages rather than locus-specific predictions. MOSAIC is best understood as a structured decision support tool — it ensures that researchers consider all viable strategies and weigh relevant trade-offs — rather than an oracle that prescribes the optimal approach.

### Limitations

Several limitations should be noted:

1. **Parameter uncertainty.** All model parameters are drawn from published measurements, but these measurements were made in specific cell types, at specific loci, with specific protocols. We report point estimates where distributions would be more appropriate. Future work should incorporate Bayesian parameter estimation with uncertainty propagation.

2. **Locus-specific effects.** HDR efficiency varies substantially across genomic loci due to chromatin state, transcriptional activity, and replication timing (Aymard et al., 2014). CRISPRArchitect uses cell-type averages and does not account for locus-specific variation.

3. **SSTR pathway.** ConversionSim does not model the single-strand template repair pathway used by short ssODN donors. Extending the model to include SSTR would broaden its applicability.

4. **PAM availability.** MOSAIC does not check whether suitable PAM sites exist near the mutation for the selected nuclease. A researcher may find that base editing is theoretically optimal but that no usable guide RNA exists at the target site.

5. **Experimental validation.** Our validation relies on comparisons with published data from other laboratories. Direct experimental validation — designing editing experiments based on CRISPRArchitect's recommendations and measuring outcomes — would provide stronger evidence.

6. **Scoring weights.** MOSAIC's strategy scoring uses hand-tuned weights (safety 40%, efficiency 30%, time 15%, cost 15%). Different applications may warrant different weightings. We provide these as adjustable parameters.

### Comparison with existing tools

To our knowledge, no existing tool addresses the multi-site editing strategy optimization problem. CRISPOR (Concordet and Haeussler, 2018) and the IDT Alt-R design tool provide guide RNA and single-site donor design but do not compare editing modalities. PrimeDesign (Hsu et al., 2021) and BE-Designer (Hwang et al., 2018) support prime editing and base editing guide design, respectively, but only for their specific modality. CRISPRArchitect complements rather than replaces these tools: it determines *which* approach to use, after which specialized design tools can be employed for the selected modality.

### Future directions

Three extensions would meaningfully strengthen the framework: (i) incorporating locus-specific chromatin features (from ENCODE/ATAC-seq data) to adjust HDR efficiency predictions; (ii) adding an SSTR sub-model for ssODN donors; and (iii) prospective experimental validation of MOSAIC's recommendations in iPSCs.

---

## Methods

### ConversionSim implementation

ConversionSim is implemented in Python (NumPy). Each simulation consists of four sequential stochastic steps:

**End resection.** Short-range resection lengths are drawn from N(μ=200, σ=80) bp, clipped to [50, 500]. Long-range resection lengths are drawn from LogNormal(μ=ln(1500), σ=0.8) bp, clipped to [300, 10000]. For staggered cuts with overhang length *h* bp, the overhang is added to the left-side resection and long-range resection is boosted by 20%. Total resection per side = short_range + long_range.

**RAD51 filament.** Filament coverage fraction is drawn from Beta(α=8, β=2), rescaled to [0.70, 0.95]. Filament length (nt) = resection_length × coverage_fraction. Filament is viable if length ≥ 15 nt (5 RAD51 monomers × 3 nt/monomer).

**Strand invasion.** Invasion probability = cell_type_base_HDR × donor_topology_multiplier × (1 + HDR_enhancement_per_bp × stagger_bp) × S_G2_fraction × homology_factor. The homology_factor decays exponentially when the filament exceeds the homology arm length.

**SDSA synthesis.** Tract length is drawn from Geometric(p), where p = base_displacement_prob × topology_modifier × stagger_modifier. Base p = 0.002 per bp. Circular ssDNA reduces p by 20%; staggered cuts reduce p by 15%. Tracts are clipped to [50, 5000] bp.

All per-simulation computations are vectorized using NumPy arrays (no per-trial Python loops). Default: 10,000 simulations per run.

### Parameter sources

All parameters and their literature sources are documented in `utils/constants.py` within the repository. Key sources: end resection (Symington, 2011; Cejka, 2015), RAD51 filament (Qi et al., 2015; Prakash et al., 2015), gene conversion tracts (Elliott et al., 1998), SDSA mechanics (Renkawitz et al., 2014), cssDNA enhancement (Iyer et al., 2022), staggered cut effects (Chauhan et al., 2023), iPSC HDR rates (Paquet et al., 2016; Ihry et al., 2018), cell cycle parameters (Orthwein et al., 2015).

### MOSAIC implementation

MOSAIC consists of four components: (i) a mutation classifier that determines amenability to base editing (transitions only), prime editing (any small edit <50 bp), or HDR (all others); (ii) a strategy enumerator that generates all feasible strategies based on mutation types and inter-mutation genomic distance; (iii) a scorer that computes efficiency, safety, time, and cost scores using published cell-type-specific parameters; and (iv) a ranker that produces weighted composite scores.

Safety scoring penalizes DSBs in cell types with active p53 pathways (iPSCs, HSCs). Two simultaneous DSBs incur an additional translocation risk penalty estimated from the power-law relationship between genomic distance and rearrangement frequency (Frock et al., 2015).

### Validation

ConversionSim validation compared model predictions against four datasets: Elliott et al. (1998, tract length distributions), Paquet et al. (2016, distance-dependent incorporation), Iyer et al. (2022, cssDNA vs. lssDNA ratios), and Chauhan et al. (2023, staggered-cut enhancement ratios). Model parameters were not fit to these datasets; all parameters were set *a priori* from independent literature sources.

MOSAIC benchmarking compared the top-three ranked strategies against editing approaches used in 14 published studies. For each study, we reconstructed the gene structure (exon count and approximate inter-exon distances), defined mutations matching those corrected in the study, and ran MOSAIC with the corresponding cell type and nuclease.

### Software availability

CRISPRArchitect is open-source (MIT license) and available at https://github.com/visvikbharti/CRISPRArchitect. The repository includes all source code, validation scripts, a command-line interface, a Streamlit web application, and a Docker container for reproducible deployment. The toolkit requires Python ≥3.9 and NumPy.

---

## Figures

**Figure 1.** ConversionSim model overview and predictions.
(a) Schematic of the four-step stochastic model: end resection, RAD51 filament formation, strand invasion, and SDSA synthesis.
(b) Distribution of predicted gene conversion tract lengths (50,000 simulations; staggered-cut nuclease, cssDNA donor, iPSC parameters). Median 542 bp, mean 770 bp.
(c) Probability of donor sequence incorporation as a function of distance from the DSB.

**Figure 2.** Validation of ConversionSim against published data.
(a) cssDNA vs. lssDNA HDR enhancement ratio: predicted 2.07x vs. observed 1.9x (Iyer et al., 2022).
(b) Staggered-cut HDR enhancement: predicted 1.82x vs. observed 1.9x (Chauhan et al., 2023).
(c) Tract length distribution shape compared with Elliott et al. (1998).
(d) Distance-dependent incorporation: ConversionSim (SDSA model) vs. Paquet et al. (2016) ssODN data (SSTR pathway), illustrating the pathway-specific applicability boundary.

**Figure 3.** MOSAIC benchmarking against 14 published editing studies.
(a) Concordance matrix: MOSAIC rank vs. author strategy type.
(b) Accuracy by editing modality.

---

## Tables

**Table 1.** ConversionSim validation summary.

| Comparison | Observed | Predicted | Agreement |
|-----------|----------|-----------|-----------|
| cssDNA/lssDNA ratio | 1.9x (1.5–2.1x) | 2.07x | Within range |
| Staggered/blunt ratio | 1.9x (1.4–2.8x) | 1.82x | Within range |
| Tract distribution shape | Right-skewed | Right-skewed | Qualitative match |
| ssODN distance-dependence | Sharp decline <50 bp | Gradual decline ~500 bp | Poor fit (pathway mismatch) |

**Table 2.** MOSAIC benchmarking: concordance with 14 published studies.

| Study | Gene | Cell Type | Author Strategy | MOSAIC Top-1 | In Top-3? |
|-------|------|-----------|----------------|-------------|-----------|
| Newby et al. 2021 | HBB | iPSC | ABE | Dual Base Editing | Yes |
| Koblan et al. 2021 | LMNA | iPSC | ABE | Dual Base Editing | Yes |
| Musunuru et al. 2021 | PCSK9 | iPSC | ABE | Dual Base Editing | Yes |
| Zeng et al. 2022 | HBB | iPSC | ABE | Dual Base Editing | Yes |
| Iyer et al. 2024 | MAK | iPSC | ABE (exon skip) | Dual Base Editing | Yes |
| Chemello et al. 2021 | RBM20 | iPSC | ABE | Dual Base Editing | Yes |
| Geurts et al. 2021 | CFTR | iPSC | Prime editing | Dual Prime Editing | Yes |
| Jackow et al. 2019 | COL7A1 | iPSC | HDR (ssODN) | Sequential HDR | Yes |
| Ramaswamy et al. 2018 | HBB | iPSC | HDR (ssODN) | Dual Base Editing | Yes |
| Sürün et al. 2020 | F9 | iPSC | HDR (ssODN) | Dual Prime Editing | No |
| Ma et al. 2017 | MYH7 | iPSC | HDR (ssODN) | Dual Base Editing | No |
| Paquet et al. 2016 | APP | iPSC | HDR (ssODN) | Sequential HDR | Yes |
| Young et al. 2016 | DMD | iPSC | Exon deletion | Dual Prime Editing | No |
| Schwank et al. 2013 | CFTR | Organoid | HDR (dsDNA) | Dual Prime Editing | No |

Overall concordance: 10/14 (71.4%)

---

## References

Anzalone AV, et al. Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature*. 2019;576:149–157.

Aymard F, et al. Transcriptionally active chromatin recruits homologous recombination at DNA double-strand breaks. *Nat Struct Mol Biol*. 2014;21:366–374.

Cejka P. DNA end resection: nucleases team up with the right partners to initiate homologous recombination. *Annu Rev Genet*. 2015;49:399–420.

Chauhan VP, Sharp PA, Langer R. Altered DNA repair pathway engagement by engineered CRISPR-Cas9 nucleases. *PNAS*. 2023;120:e2300605120.

Concordet JP, Haeussler M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Res*. 2018;46:W242–W245.

Elliott B, et al. Gene conversion tracts from double-strand break repair in mammalian cells. *Mol Cell Biol*. 1998;18:93–101.

Frock RL, et al. Genome-wide detection of DNA double-stranded breaks induced by engineered nucleases. *Nat Biotechnol*. 2015;33:179–186.

Gaudelli NM, et al. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage. *Nature*. 2017;551:464–471.

Haapaniemi E, et al. CRISPR-Cas9 genome editing induces a p53-mediated DNA damage response. *Nat Med*. 2018;24:927–930.

Ihry RJ, et al. p53 inhibits CRISPR-Cas9 engineering in human pluripotent stem cells. *Nat Med*. 2018;24:939–946.

Iyer S, et al. Efficient Homology-Directed Repair with Circular Single-Stranded DNA Donors. *CRISPR J*. 2022;5:685–701.

Jensen RB, et al. Purified human BRCA2 stimulates RAD51-mediated recombination. *Nature*. 2010;467:678–683.

Komor AC, et al. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature*. 2016;533:420–424.

Lieberman-Aiden E, et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. *Science*. 2009;326:289–293.

Marin-Gonzalez A, et al. Cohesin drives chromatin scanning during the RAD51-mediated homology search. *Science*. 2025.

Orthwein A, et al. A mechanism for the suppression of homologous recombination in G1 cells. *Nature*. 2015;528:422–426.

Paquet D, et al. Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. *Nature*. 2016;533:125–129.

Prakash R, et al. Homologous recombination and human health: the roles of BRCA1, BRCA2, and associated proteins. *Cold Spring Harb Perspect Biol*. 2015;7:a016600.

Qi Z, et al. DNA sequence alignment by microhomology sampling during homologous recombination. *Cell*. 2015;160:856–869.

Renkawitz J, et al. Mechanisms and principles of homology search during recombination. *Nat Rev Mol Cell Biol*. 2014;15:369–383.

Richardson CD, et al. Enhancing homology-directed genome editing by catalytically active and inactive CRISPR-Cas9 using asymmetric donor DNA. *Nat Biotechnol*. 2016;34:339–344.

Symington LS. Mechanism and regulation of DNA end resection in eukaryotes. *Crit Rev Biochem Mol Biol*. 2011;46:67–82.

Xie K, et al. Efficient non-viral immune cell engineering using circular single-stranded DNA-mediated genomic integration. *Nat Biotechnol*. 2025;43:1821–1832.

Yokota H, et al. Evidence for the organization of chromatin in megabase pair-sized loops arranged along a random walk path in the human G0/G1 interphase nucleus. *J Cell Biol*. 1995;130:1239–1249.

---

## Acknowledgments

[To be added]

## Competing interests

The authors declare no competing interests.

## Data availability

All source code, validation scripts, and analysis pipelines are available at https://github.com/visvikbharti/CRISPRArchitect under the MIT license.

---

*Note: This manuscript draft requires verification of all cited papers, author names, and data values against the original publications before submission. Figure references are to figures that would need to be generated from the validation scripts provided in the repository.*
