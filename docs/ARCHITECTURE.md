# CRISPRArchitect — Architecture & Design Document

## For Scientists: Understanding the Code From Scratch

This document explains the entire CRISPRArchitect system from first principles. You don't need a computer science background to understand it — we explain every concept as it relates to the biology.

---

## What This Software Does

CRISPRArchitect answers one fundamental question:

> **"I have mutations at two (or more) positions in a large gene. What is the best way to correct them using CRISPR in my iPSCs?"**

It does this by combining four types of analysis:

1. **ConversionSim** — Simulates what happens at the molecular level when HDR occurs at a single DSB
2. **MOSAIC** — Takes your gene structure and mutations as input, considers all possible editing strategies, and recommends the best one
3. **cssDNA-TopoPred** — Checks if your cssDNA donor template has secondary structures that could block HDR
4. **ChromBridge** — Calculates the physical 3D distance between your target sites and determines if single-template editing is physically possible

---

## How the Code is Organized

```
crisprarchitect/
├── utils/               # Shared building blocks
│   ├── constants.py     # All biological parameters (with citations)
│   ├── sequence.py      # DNA sequence tools (complement, PAM finding, etc.)
│   └── plotting.py      # Visualization functions
├── conversion_sim/      # Module 1: Gene Conversion Tract Simulator
├── mosaic/              # Module 2: Strategy Optimizer
├── topopred/            # Module 3: cssDNA Structure Analyzer
├── chrombridge/         # Module 4: 3D Distance Predictor
├── tests/               # Automated tests
├── examples/            # Ready-to-run example scripts
└── docs/                # Documentation (you are here)
```

### Why this structure?

Each module is **independent** — you can use any one without the others. But they also **integrate**: MOSAIC can call ConversionSim to check if a single-template strategy is viable, and ChromBridge to assess translocation risk.

---

## Module 1: ConversionSim — Gene Conversion Tract Simulator

### What it simulates

When CRISPR/Cas9 makes a cut in your DNA and you provide a donor template, the cell uses **Homology-Directed Repair (HDR)** to copy information from the donor into the genome. But HDR doesn't copy the entire donor — it only copies a limited stretch called the **gene conversion tract**.

ConversionSim simulates this process step by step:

```
Step 1: END RESECTION
        After the DSB, enzymes chew back the 5' ends of the DNA,
        creating 3' single-stranded overhangs.

        5'------>     <------5'    (DSB)
             ↓ Resection ↓
        5'----->      <-----5'
                ======        ======    (3' ssDNA overhangs)

Step 2: RAD51 FILAMENT
        RAD51 protein coats the ssDNA, forming a helical filament.
        This filament is the "search engine" that finds the donor template.

        ========RAD51=RAD51=RAD51========

Step 3: STRAND INVASION
        The RAD51 filament invades the donor template at the homology arm,
        creating a D-loop structure.

        Donor:  =========================================
                     /--- RAD51 filament invades ---\
        Genome: ====                                  ====

Step 4: DNA SYNTHESIS (SDSA)
        A polymerase extends the invading strand, copying donor sequence.
        This is where the gene conversion tract is created.

        Donor:  =========================================
                     |-----> copying ------->|
        Genome: ====                                  ====

        The LENGTH of this copying step = the GENE CONVERSION TRACT LENGTH
        Typical range: 200 bp to 2,000 bp (most are <1 kb)

Step 5: DISPLACEMENT (SDSA)
        The newly synthesized strand is kicked off the donor and re-anneals
        to the other side of the break. This is the SDSA pathway.

        This displacement is what LIMITS tract length.
```

### How the simulation works (Monte Carlo)

Instead of solving equations, we simulate the process thousands of times with random variation (Monte Carlo method). Each simulation is like running the experiment once in a virtual cell:

```python
for i in range(10000):  # Run 10,000 virtual experiments
    resection_length = random(mean=2000, std=1000)  # Variable resection
    filament_ok = resection_length > minimum_threshold
    if filament_ok:
        synthesis_length = geometric_random(p=0.002)  # Stochastic synthesis
        record(synthesis_length)
```

The distributions of outcomes tell us:
- **Median tract length** — the "typical" result
- **95th percentile** — how long the tract could be in lucky cases
- **P(reaching distance X)** — the probability that an edit at distance X from the cut gets incorporated

### Key parameters and where they come from

| Parameter | Value | Source |
|-----------|-------|--------|
| Short-range resection | ~200 bp | Symington, Ann Rev Genet, 2011 |
| Long-range resection | ~2,000 bp (median) | Cejka, Ann Rev Genet, 2015 |
| Conversion tract (median) | ~500 bp | Elliott et al., MCB, 1998 |
| SDSA displacement rate | ~0.2% per bp | Fitted to tract data |
| Stagger boost (enFnCas9) | +20% resection | Chauhan et al., PNAS, 2023 |

---

## Module 2: MOSAIC — Strategy Optimizer

### What it does

Given:
- Your gene structure (exons, introns, sizes)
- Mutation positions (which exons, what type of change)
- Your cell type (iPSC, HEK293T, etc.)
- Your nuclease (enFnCas9, SpCas9, etc.)

MOSAIC enumerates **every feasible editing strategy** and scores them.

### Strategy types it considers

```
1. SINGLE_TEMPLATE_HDR
   - One donor template spanning both mutations
   - Only feasible if mutations are <5 kb apart in the genome
   - Scored: efficiency × feasibility

2. DUAL_TEMPLATE_SIMULTANEOUS_HDR
   - Two donors + two sgRNAs delivered together
   - Efficiency = (single-site efficiency)²
   - Risk: translocations, p53 selection

3. SEQUENTIAL_HDR
   - Edit one site, clone, validate, then edit the other
   - Safest approach for iPSCs
   - Takes longer (8-16 weeks total)

4. BASE_EDITING (if mutation is a transition)
   - No DSB needed! CBE or ABE corrects the base directly
   - Highest safety score
   - Only works for C>T, G>A, A>G, T>C changes

5. PRIME_EDITING
   - No DSB needed (uses a nickase)
   - Can correct any small mutation
   - Lower efficiency than base editing

6. HYBRID approaches
   - Base edit one site + HDR the other
   - Prime edit one site + HDR the other
   - Only one DSB total → better safety in iPSCs

7. EXON_DELETION
   - Delete the exon(s) containing mutations
   - Only if the protein tolerates the deletion
```

### Scoring system

Each strategy gets scored on 4 axes:

```
EFFICIENCY (0-1): How likely are you to get a correctly edited clone?
  - Based on published HDR rates for your cell type
  - Multiplied by donor type effectiveness
  - For dual-site: efficiency₁ × efficiency₂

SAFETY (0-1): How safe is this for your cells?
  - Number of DSBs (more DSBs = lower safety)
  - Translocation risk (from ChromBridge)
  - p53 selection risk (critical for iPSCs!)
  - Base/prime editing gets highest safety scores

TIME (0-1): How fast can you get results?
  - Sequential: 2× time per round
  - Simultaneous: 1 round but more screening
  - Base editing: fastest (no cloning/screening needed if high efficiency)

COST (0-1): How expensive?
  - Number of donors to produce
  - Screening burden (clones to pick)
  - Number of editing rounds
```

---

## Module 3: cssDNA-TopoPred — Structure Analyzer

### Why DNA structure matters for HDR

Your cssDNA donor is a single-stranded circle. Unlike double-stranded DNA (which is a rigid helix), ssDNA is flexible and **folds onto itself** forming secondary structures:

```
Linear representation:
5'---AGCTTTT-GGGG-ATCG-GGGG-TTTT-GGGG-AAAC-GGGG---3'

What actually happens (G-quadruplex):
              G-G
             / | \
            G  G  G
            |  |  |
            G  G  G
             \ | /
              G-G
    (four G-runs stack into a square)

This structure is EXTREMELY stable and blocks RAD51 binding!
```

### What TopoPred analyzes

1. **G-quadruplex scanning**: Finds G-rich sequences that form G4 structures
   - Pattern: G₃₊N₁₋₇G₃₊N₁₋₇G₃₊N₁₋₇G₃₊
   - Flags those overlapping with homology arms

2. **Hairpin prediction**: Finds self-complementary regions
   - Uses energy-based scoring (nearest-neighbor model)
   - Stable hairpins (ΔG < -3 kcal/mol) are flagged

3. **Accessibility scoring**: For each nucleotide in the homology arms:
   - Is it free (accessible for RAD51)? → Score = 1
   - Is it in a hairpin stem? → Score = 0
   - Is it in a G4? → Score = 0

4. **Optimization**: Suggests synonymous codon changes to disrupt structures

---

## Module 4: ChromBridge — 3D Distance Predictor

### The physics of chromatin

Chromatin (DNA wrapped around histones) behaves like a **polymer chain**. We can predict how far apart two points on the chain are in 3D space using polymer physics:

```
Imagine a ball of yarn (= chromosome in the nucleus)

If you pick two points on the yarn:
  - If they're close along the yarn (small genomic distance):
    → They're likely close in 3D space too

  - If they're far along the yarn (large genomic distance):
    → They could be anywhere in the ball
    → Average 3D distance increases, but slower than linear

Mathematical model (Gaussian chain):
  <R²> = N × b²

  Where:
  - R = 3D distance between the two points
  - N = number of chain segments between them
  - b = size of each segment (Kuhn length ≈ 300 nm for chromatin)
```

### What ChromBridge calculates

1. **3D distance between two genomic loci**
   - Converts genomic distance (bp) to physical distance (nm)
   - Uses polymer chain model calibrated to FISH and Hi-C data

2. **Donor bridgeability**
   - Calculates the random coil diameter of your cssDNA donor
   - Compares it to the inter-locus distance
   - Result: can the donor physically reach both sites simultaneously?
   - (Spoiler: almost never, unless sites are <5 kb apart)

3. **Translocation risk**
   - When two DSBs exist, broken ends can be mis-joined
   - Risk correlates with 3D proximity (Hi-C contact frequency)
   - Empirical model: P(translocation) ~ s^(-1.08) where s = genomic distance

---

## Data Flow: How the Modules Connect

```
USER INPUT:
  Gene name, exon coordinates, mutation positions,
  cell type, nuclease choice, donor sequence

                    ↓

MOSAIC (Strategy Optimizer)
  ├── Fetches gene structure → calculates distances
  ├── Classifies mutations → determines which tools are applicable
  ├── Calls ConversionSim → checks if single-template HDR is viable
  ├── Calls ChromBridge → assesses translocation risk
  ├── Enumerates all strategies
  ├── Scores and ranks
  └── Generates recommendation report

                    ↓

If HDR strategy is recommended:
  ├── User designs donor template
  └── cssDNA-TopoPred → analyzes donor for structural problems
                    ↓

FINAL OUTPUT:
  - Ranked strategies with scores
  - Recommended approach with rationale
  - Donor template quality assessment
  - Risk warnings (translocation, p53, etc.)
  - Experimental design suggestions
```

---

## Key Design Decisions

### 1. Why Monte Carlo instead of analytical formulas?

The HDR process involves multiple stochastic steps (resection, filament formation, synthesis). While we could derive approximate formulas, Monte Carlo simulation:
- Captures the **full distribution** of outcomes, not just the mean
- Easily incorporates correlated parameters
- Is easy to extend with new biology
- Produces intuitive outputs (histograms, survival curves)

### 2. Why no external database dependencies for basic use?

We designed CRISPRArchitect to work with just `numpy` and `scipy`. This means:
- No API calls needed (works offline)
- No database setup
- Reproducible results
- Easy to install

For advanced use, you CAN connect to Ensembl, UCSC, and ENCODE, but the core simulations work without internet.

### 3. Why separate modules?

Each module answers a distinct biological question. A researcher working on donor design needs TopoPred but not ChromBridge. A researcher comparing strategies needs MOSAIC but not ConversionSim's internals. Modularity lets you use what you need.

---

## Extending CRISPRArchitect

### Adding a new nuclease

Edit `utils/constants.py` and add to `NUCLEASE_PARAMS`:

```python
"NewCas9": {
    "pam": "NNGG",
    "cut_type": "staggered_5prime",
    "stagger_bp": 3,
    "hdr_multiplier": 1.3,
    "specificity": "high",
    "description": "My new Cas9 variant",
}
```

### Adding a new cell type

Edit `utils/constants.py` and add to `CELL_TYPE_PARAMS`:

```python
"my_cell_type": {
    "hdr_base_efficiency": 0.15,
    "cell_cycle_s_g2_fraction": 0.40,
    "p53_active": True,
    "viability_single_dsb": 0.65,
    "viability_dual_dsb": 0.40,
    "description": "My custom cell type",
}
```

### Adding a new donor type

Edit `DONOR_TOPOLOGY_MULTIPLIER` in `utils/constants.py`.

---

## Limitations and Caveats

1. **All models are approximations.** Real HDR involves thousands of proteins, chromatin remodeling, cell cycle regulation, and stochastic processes we cannot fully capture. These simulations give useful estimates, not exact predictions.

2. **Parameters vary by locus.** The same cell type can have very different HDR rates at different genomic loci due to chromatin state, transcription, replication timing, etc.

3. **No off-target analysis.** CRISPRArchitect does not predict Cas9 off-target sites. Use dedicated tools (Cas-OFFinder, CRISPOR, etc.) for that.

4. **Translocation risk is approximate.** Without cell-type-specific Hi-C data, ChromBridge uses a general polymer model. For precise translocation risk, experimental Hi-C data would be needed.

5. **Secondary structure prediction is simplified.** TopoPred uses a basic energy model for hairpins. For rigorous structure prediction, consider ViennaRNA or Mfold.
