# Supporting Information

## S1 Table. ConversionSim parameter values and literature sources.

All parameters used in the ConversionSim Monte Carlo model, their default values, distributions, biological basis, and the published studies from which they were derived. No parameters were fit to the validation datasets reported in the main text.

### End Resection Parameters

| Parameter | Value | Distribution | Biological Basis | Source |
|-----------|-------|-------------|-----------------|--------|
| Short-range resection mean | 200 bp | Normal(200, 80) | MRN/CtIP-initiated resection generates initial 3' ssDNA overhangs of ~100-300 bp | Symington, 2011 [14] |
| Short-range resection s.d. | 80 bp | — | Cell-to-cell variability in MRN recruitment | Estimated |
| Short-range resection range | [50, 500] bp | Clipped | Physical limits: minimum ~50 bp for MRN engagement; CtIP rarely exceeds 500 bp | Cejka, 2015 [15] |
| Long-range resection median | ~1,500 bp | LogNormal(ln(1500), 0.8) | EXO1/BLM-DNA2 extends resection to 1-5+ kb | Symington, 2011 [14] |
| Long-range resection range | [300, 10000] bp | Clipped | Minimum: must exceed short-range; Maximum: rarely >10 kb observed | Cejka, 2015 [15] |
| Staggered cut resection boost | 20% | Multiplicative on long-range | 5' overhangs provide entry point for EXO1; promote resection initiation | Chauhan et al., 2023 [8] |

### RAD51 Filament Parameters

| Parameter | Value | Distribution | Biological Basis | Source |
|-----------|-------|-------------|-----------------|--------|
| RAD51 monomer footprint | 3 nt | Fixed | Each RAD51 monomer covers 3 nucleotides of ssDNA | Ogawa et al., 1993 |
| Minimum viable filament | 15 nt (5 monomers) | Fixed | Minimum nucleus for stable strand invasion; based on 8-nt sampling mechanism | Qi et al., 2015 [17] |
| Filament coverage fraction | Mean ~85% | Beta(8, 2) rescaled to [0.70, 0.95] | Not all resected ssDNA is coated by RAD51; competition with RPA, secondary structures | Estimated from in vitro data |

### Strand Invasion Parameters

| Parameter | Value | Biological Basis | Source |
|-----------|-------|-----------------|--------|
| iPSC base HDR rate | 0.08 | 5-15% HDR in iPSCs before optimization | Multiple studies; Paquet et al., 2016 [13] |
| HEK293T base HDR rate | 0.25 | 20-40% HDR in HEK293T | Multiple studies |
| K562 base HDR rate | 0.20 | 15-30% HDR in K562 | Iyer et al., 2022 [6] |
| iPSC S/G2 fraction | 0.35 | iPSCs spend ~35% of time in S/G2, when HDR is active | Orthwein et al., 2015 |
| HEK293T S/G2 fraction | 0.55 | Fast-cycling transformed cells | Estimated |
| HDR enhancement per bp overhang | 0.15 | ~15% HDR increase per bp of 5' overhang, based on vCas9 data | Chauhan et al., 2023 [8] |

### Donor Topology Multipliers

| Donor Type | Multiplier | Basis | Source |
|-----------|-----------|-------|--------|
| Linear dsDNA | 1.0 | Baseline | — |
| Plasmid dsDNA | 0.8 | Slightly reduced due to supercoiling | Estimated |
| Linear ssDNA | 1.5 | ~1.5x better than linear dsDNA | Iyer et al., 2022 [6] |
| Circular ssDNA (cssDNA) | 3.0 | ~2x better than linear ssDNA; ~3x vs linear dsDNA | Iyer et al., 2022 [6] |
| AAV ssDNA | 4.0 | Highest efficiency for large knock-ins in iPSCs | Martin et al., 2019 |

### SDSA Synthesis Parameters

| Parameter | Value | Distribution | Biological Basis | Source |
|-----------|-------|-------------|-----------------|--------|
| Base displacement probability | 0.002 per bp | Geometric(p) | Probability of D-loop collapse at each bp of synthesis; produces mean tract ~500 bp | Fitted to Elliott et al., 1998 [12] tract distributions |
| Circular ssDNA modifier | 0.80 (reduces p by 20%) | Multiplicative | Circular topology stabilizes D-loop; exonuclease resistance increases effective concentration | Iyer et al., 2022 [6] |
| Staggered cut modifier | 0.85 (reduces p by 15%) | Multiplicative | Longer resection/filament stabilizes strand invasion intermediate | Chauhan et al., 2023 [8] |
| Tract length range | [50, 5000] bp | Clipped | Physical limits on SDSA synthesis | Elliott et al., 1998 [12] |

### MOSAIC Scoring Parameters

| Parameter | Default Value | Adjustable? | Basis |
|-----------|--------------|-------------|-------|
| Safety weight | 0.40 | Yes | Safety prioritized in iPSC applications due to p53 risk | Ihry et al., 2018 [1] |
| Efficiency weight | 0.30 | Yes | Important but secondary to safety | — |
| Time weight | 0.15 | Yes | Practical consideration | — |
| Cost weight | 0.15 | Yes | Practical consideration | — |
| DSB safety penalty (iPSC) | 0.25 per DSB | No | Strong p53 response to DSBs in iPSCs | Ihry et al., 2018 [1] |
| Simultaneous dual-DSB penalty | Additional 0.10 | No | Translocation risk between concurrent DSBs | Frock et al., 2015 [21] |
| Base editing efficiency | ~29% per site | No | Typical ABE/CBE efficiency in iPSCs | Multiple studies |
| Prime editing efficiency | ~13% per site | No | Typical PE efficiency in iPSCs (lower than BE) | Anzalone et al., 2019 [5] |

---

## S1 Fig. ConversionSim sensitivity analysis.

We assessed the sensitivity of ConversionSim predictions to variation in three key parameters by running 10,000 simulations for each parameter value while holding all other parameters at their defaults (staggered-cut nuclease, cssDNA donor, iPSC cell type, 300 bp homology arms).

**(A) Sensitivity to SDSA displacement probability.** The base displacement probability *p* (default: 0.002 per bp) is the primary determinant of tract length. Varying *p* from 0.001 to 0.003 shifts the median tract length from approximately 700 bp to approximately 230 bp. This parameter has the strongest influence on model output and is the most uncertain, as direct measurements of SDSA displacement rates in mammalian cells are limited.

**(B) Sensitivity to long-range resection median.** Varying the median of the log-normal resection distribution from 1,000 to 3,000 bp primarily affects HDR success rate (through the filament length → invasion probability pathway) rather than tract length. The tract length distribution is relatively insensitive to resection length, because SDSA synthesis is decoupled from the initial resection step.

**(C) Sensitivity to donor topology multiplier.** The donor topology multiplier affects both the invasion probability and the SDSA displacement rate. Varying the cssDNA multiplier from 1.0 (equivalent to linear dsDNA) to 4.0 (equivalent to AAV) shifts the predicted cssDNA/lssDNA enhancement ratio from 1.0-fold to approximately 2.7-fold. The default value of 3.0 produces a ratio of 2.07-fold, consistent with the observed range of 1.5-2.1-fold reported by Iyer et al. [6].

**Conclusion:** The model is most sensitive to the SDSA displacement probability, which directly controls tract length. This parameter should be the priority target for future experimental calibration, ideally through direct measurement of conversion tract lengths with long donor templates in iPSCs.

---

## S2 Fig. MOSAIC scoring weight sensitivity.

We assessed how the top-ranked editing strategy changes as the safety weight is varied from 0.10 to 0.80, with the remaining weight distributed proportionally among efficiency (3/6), time (1.5/6), and cost (1.5/6).

**Test case:** NF1 gene, mutations in exon 20 (G>A, transition) and exon 50 (C>T, transition), iPSC cell type, SpCas9 nuclease.

**Results:**

| Safety Weight | Efficiency Weight | Top-Ranked Strategy |
|--------------|-------------------|-------------------|
| 0.10 | 0.53 | Sequential HDR |
| 0.20 | 0.47 | Sequential HDR |
| 0.30 | 0.41 | Dual Base Editing |
| 0.40 (default) | 0.35 | Dual Base Editing |
| 0.50 | 0.29 | Dual Base Editing |
| 0.60 | 0.24 | Dual Base Editing |
| 0.70 | 0.18 | Dual Base Editing |
| 0.80 | 0.12 | Dual Base Editing |

**Interpretation:** Dual base editing is the top-ranked strategy across a wide range of safety weights (0.30-0.80), reflecting the fundamental advantage of zero DSBs in iPSCs. Only when the safety weight drops below 0.30 — placing disproportionate emphasis on raw efficiency — does sequential HDR overtake base editing. This is because sequential HDR has higher per-round efficiency (~20% per site) but introduces DSBs. The robustness of the base editing recommendation across weight configurations provides confidence that it is not an artifact of specific weight choices.

For mutation types not amenable to base editing (transversions, indels), the crossover between sequential HDR and dual prime editing occurs at a safety weight of approximately 0.35-0.40.

---

## S1 Text. Supporting analysis modules.

### ChromBridge: 3D chromatin distance estimation

ChromBridge estimates the mean three-dimensional distance between two genomic loci using a confined polymer model of chromatin. The chromatin fiber is modeled as a Gaussian chain with persistence length ~300 nm (approximately 30 kb in euchromatin), incorporating confinement within chromosome territories of radius ~2 μm.

The mean-squared end-to-end distance is:

⟨R²⟩_free = N × b²

where N is the number of Kuhn segments (genomic distance / 30,000 bp) and b is the Kuhn length (300 nm). Confinement is applied as:

⟨R²⟩_confined = R²_max × (1 - exp(-⟨R²⟩_free / R²_max))

where R²_max is the squared chromosome territory radius.

The model was calibrated against published FISH measurements of inter-locus distances in human cells (Yokota et al., 1995; Mateos-Langerak et al., 2009) and is consistent with the power-law decay of Hi-C contact frequencies (Lieberman-Aiden et al., 2009).

**Bridgeability analysis.** For a donor template of length L nucleotides (ssDNA), the radius of gyration is estimated as:

R_g = √(L × l_k / 6)

where l_k is the Kuhn length for ssDNA (approximately 1.5 nm). The effective coil diameter (2 × R_g) is compared to the predicted inter-locus distance. If the coil diameter is less than the inter-locus distance, bridging by a single donor molecule is physically implausible.

**Translocation risk.** When two DSBs are made simultaneously, the probability of inter-break joining (translocation, deletion, or inversion) correlates with 3D proximity. We estimate this from the power-law relationship: P(rearrangement) ~ d^(-γ) where d is genomic distance and γ ≈ 0.5, calibrated to published paired-sgRNA deletion frequencies (Canver et al., 2015).

**Limitations.** ChromBridge uses a generic polymer model and does not incorporate cell-type-specific or locus-specific Hi-C data. The model provides order-of-magnitude estimates suitable for feasibility assessment but not precise distance predictions.

### cssDNA-TopoPred: donor template secondary structure analysis

cssDNA-TopoPred analyzes circular ssDNA donor templates for secondary structures that could impair HDR by sequestering homology arm nucleotides.

**G-quadruplex detection.** G4 motifs are identified using regex-based scanning for the pattern G₃₊N₁₋₇G₃₊N₁₋₇G₃₊N₁₋₇G₃₊ on both sense and antisense orientations. Each motif is scored by stability (tetrad count and loop lengths).

**Hairpin prediction.** Potential hairpins are identified by scanning for complementary regions within the sequence. Free energy is estimated using a simplified nearest-neighbor model (GC stack: -2.0 kcal/mol, AT stack: -1.0 kcal/mol, loop penalty: 4.0 + 1.4 × ln(loop_length) kcal/mol). Hairpins with ΔG < -3.0 kcal/mol at 37°C are classified as stable.

**Accessibility scoring.** For each homology arm, the fraction of nucleotides not involved in predicted secondary structures (G4 stems or hairpin stems) is computed. Arms with less than 60% accessibility are flagged as potentially impaired.

**Limitations.** The secondary structure predictions are simplified relative to full thermodynamic algorithms (e.g., ViennaRNA). The relationship between predicted accessibility and actual HDR efficiency has not been experimentally validated.

### LoopSim: cohesin loop extrusion simulator

LoopSim simulates cohesin-mediated chromatin loop extrusion from DSB sites, inspired by the finding that cohesin drives chromatin scanning during the RAD51-mediated homology search (Marin-Gonzalez et al., 2025).

**Model.** The chromatin fiber is represented as a one-dimensional lattice (default: 1 kb per bead). After a DSB, a cohesin complex is loaded at the break site and extrudes chromatin bidirectionally at approximately 1 kb/s. Extrusion stalls probabilistically at CTCF sites in convergent orientation and slows through heterochromatin.

**Relevance and limitations.** The cohesin-mediated search domain is relevant for endogenous donor recombination (sister chromatid repair, loss of heterozygosity risk) but has limited direct bearing on exogenous donor-mediated HDR, where the donor reaches the DSB through three-dimensional diffusion in the nucleoplasm. We include LoopSim for completeness but emphasize that its predictions do not directly inform exogenous donor experimental design.

---

## S2 Text. Detailed MOSAIC benchmarking results.

### Concordant cases (10/14)

**Case 1: Newby et al. (2021) — HBB sickle cell correction with ABE.**
Gene: HBB (3 exons). Mutation: E6V (GAG>GTG), transition. Authors used ABE in iPSCs. MOSAIC ranked Dual Base Editing #1 (score 0.756). Match.

**Case 2: Koblan et al. (2021) — LMNA progeria correction with ABE.**
Gene: LMNA (12 exons). Mutation: C608T splice site, transition. Authors used ABE. MOSAIC ranked Dual Base Editing #1. Match.

**Case 3: Musunuru et al. (2021) — PCSK9 base editing.**
Gene: PCSK9. Mutation: splice site disruption, transition. Authors used ABE. MOSAIC ranked Dual Base Editing #1. Match.

**Case 4: Zeng et al. (2022) — HBB correction with ABE.**
Gene: HBB. Mutation: sickle cell variant, transition. Authors used ABE in iPSCs. MOSAIC ranked Dual Base Editing #1. Match.

**Case 5: Iyer et al. (2024) — MAK exon skipping with ABE.**
Gene: MAK. Mutation: Alu insertion causing mis-splicing. Authors used ABE to disrupt splice site. MOSAIC ranked Dual Base Editing #1 (recognizing the mutation as transition-amenable). Match.

**Case 6: Chemello et al. (2021) — RBM20 correction with ABE.**
Gene: RBM20. Mutation: R636S missense, transition. Authors used ABE in iPSC-derived cardiomyocytes. MOSAIC ranked Dual Base Editing #1. Match.

**Case 7: Geurts et al. (2021) — CFTR correction with prime editing.**
Gene: CFTR (27 exons). Mutation: F508del (3-bp deletion). Authors used prime editing in organoids. MOSAIC ranked Dual Prime Editing #1. Match.

**Case 8: Jackow et al. (2019) — COL7A1 HDR correction.**
Gene: COL7A1 (118 exons). Mutation: point mutation. Authors used HDR with ssODN in iPSCs. MOSAIC ranked Sequential HDR in top 3. Match.

**Case 9: Ramaswamy et al. (2018) — HBB HDR correction.**
Gene: HBB. Mutation: sickle cell, transition. Authors used HDR with ssODN. MOSAIC ranked Dual Base Editing #1 (reflecting that ABE was available by this time). Match within top 3 (HDR ranked #3).

**Case 10: Paquet et al. (2016) — APP Alzheimer's mutations with HDR.**
Gene: APP (18 exons). Mutations: A673T and Swedish mutation. Authors used ssODN HDR. MOSAIC ranked Sequential HDR in top 3. Match.

### Discordant cases (4/14)

**Case 11: Sürün et al. (2020) — F9 hemophilia B correction with HDR.**
Gene: F9 (8 exons). Mutation: missense. Authors used HDR with ssODN. MOSAIC ranked Dual Prime Editing #1, Sequential HDR #3 (outside top 3 for HDR). Disagreement analysis: The mutation is a transversion, so base editing is not available. MOSAIC prefers prime editing (no DSB) over HDR (requires DSB). The authors may have chosen HDR due to higher efficiency at this specific locus or because prime editing efficiency at the time was lower than current estimates.

**Case 12: Ma et al. (2017) — MYH7 hypertrophic cardiomyopathy correction.**
Gene: MYH7 (40 exons). Mutation: R403Q missense, transition (G>A). Authors used HDR with ssODN in human embryos. MOSAIC ranked Dual Base Editing #1. Disagreement analysis: This study was conducted before ABE was widely available for this application. MOSAIC's preference for ABE is prospectively reasonable given current tools.

**Case 13: Young et al. (2016) — DMD exon deletion.**
Gene: DMD (79 exons). Authors deleted exons 45-55 using paired sgRNAs (NHEJ). MOSAIC ranked Dual Prime Editing #1. Disagreement analysis: MOSAIC lacks a "single-cut exon reframing" or "multi-exon deletion" strategy optimized for the DMD reading frame restoration approach. This is a genuine gap in MOSAIC's strategy enumeration.

**Case 14: Schwank et al. (2013) — CFTR correction with HDR in organoids.**
Gene: CFTR. Mutation: F508del. Authors used dsDNA donor HDR. MOSAIC ranked Dual Prime Editing #1. Disagreement analysis: This study predated prime editing. MOSAIC's recommendation reflects current tool availability and is prospectively appropriate.

### Summary of disagreement patterns

Three of four disagreements (Cases 11, 12, 14) involve MOSAIC recommending a DSB-free approach that was not available or not yet established when the study was conducted. This suggests that MOSAIC's forward-looking recommendations may be more appropriate than the historical concordance rate implies. The remaining disagreement (Case 13) reflects a genuine gap in MOSAIC's strategy enumeration that should be addressed in future versions.
