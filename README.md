# CRISPRArchitect

**A computational toolkit for multi-site genome editing strategy optimization**

CRISPRArchitect is a suite of bioinformatics tools that helps researchers plan complex genome editing experiments, particularly those involving multiple mutation sites in large multi-exon genes. It integrates knowledge from DNA repair biology, 3D genome organization, polymer physics, and donor template chemistry.

## Motivation

When correcting multiple mutations in a large gene (e.g., a gene with 60 exons where mutations exist in exon 40 and exon 60), researchers face critical questions:

- Can a single donor template correct both sites? (Usually no — but why, quantitatively?)
- What is the optimal editing strategy? (Sequential HDR? Simultaneous? Base editing?)
- What is the translocation risk of making two simultaneous DSBs?
- Is my cssDNA donor folding into structures that block HDR?
- How far from the cut site will my donor sequence be incorporated?

**No existing tool answers these questions.** CRISPRArchitect fills this gap.

## Tool Suite

### 1. ConversionSim — Gene Conversion Tract Simulator
Monte Carlo simulation of HDR mechanics: end resection, RAD51 filament formation, strand invasion, and SDSA-mediated synthesis. Predicts the distribution of gene conversion tract lengths and the probability of incorporating edits at various distances from the DSB.

### 2. MOSAIC — Multi-locus Optimized Strategy for Allele-specific Integrated Correction
Given a gene structure and mutation positions, MOSAIC enumerates all feasible editing strategies and ranks them by predicted efficiency, safety, time, and cost.

### 3. cssDNA-TopoPred — Secondary Structure Analyzer for ssDNA Donors
Analyzes cssDNA donor templates for problematic secondary structures (hairpins, G-quadruplexes) that could impair HDR. Scores homology arm accessibility and suggests sequence optimizations.

### 4. ChromBridge — 3D Chromatin-Aware Feasibility Predictor
Uses polymer physics models calibrated to Hi-C data to estimate 3D distances between loci, translocation risk, and the physical feasibility of single-template multi-site editing.

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from crisprarchitect.conversion_sim import ConversionSimulator
from crisprarchitect.mosaic import StrategyOptimizer
from crisprarchitect.topopred import DonorAnalyzer
from crisprarchitect.chrombridge import ChromatinDistancePredictor

# Example: Simulate gene conversion tract lengths
sim = ConversionSimulator(
    cut_type="staggered_5prime",  # enFnCas9-like
    overhang_length=4,
    donor_topology="circular_ss",
    n_simulations=10000
)
results = sim.run()
results.plot_tract_distribution()
results.probability_at_distance(500)  # P(edit incorporated) at 500bp from cut
```

## Project Structure

```
crisprarchitect/
├── conversion_sim/      # Gene conversion tract simulator
│   ├── __init__.py
│   ├── resection.py     # End resection models
│   ├── filament.py      # RAD51 filament formation
│   ├── synthesis.py     # Template-directed synthesis (SDSA)
│   ├── simulator.py     # Main Monte Carlo engine
│   └── models.py        # Parameter sets and calibration data
├── mosaic/              # Multi-site strategy optimizer
│   ├── __init__.py
│   ├── gene_structure.py    # Gene annotation fetcher
│   ├── mutation_classifier.py  # Classify mutation types
│   ├── strategy_enumerator.py  # Enumerate editing strategies
│   ├── scorer.py        # Score and rank strategies
│   └── reporter.py      # Generate recommendation reports
├── topopred/            # cssDNA secondary structure analyzer
│   ├── __init__.py
│   ├── g_quadruplex.py  # G-quadruplex motif scanner
│   ├── hairpin.py       # Hairpin/stem-loop predictor
│   ├── accessibility.py # Homology arm accessibility scorer
│   └── optimizer.py     # Sequence optimization suggestions
├── chrombridge/         # 3D chromatin distance predictor
│   ├── __init__.py
│   ├── polymer_model.py # Polymer physics (WLC/Gaussian chain)
│   ├── distance.py      # 3D distance estimation
│   ├── translocation.py # Translocation risk scoring
│   └── tad_analysis.py  # TAD boundary analysis
├── utils/               # Shared utilities
│   ├── __init__.py
│   ├── constants.py     # Biological constants
│   ├── plotting.py      # Visualization helpers
│   └── sequence.py      # Sequence manipulation utilities
├── tests/               # Unit tests
├── examples/            # Example scripts and notebooks
└── docs/                # Documentation
```

## Scientific Basis

Every parameter and model in CRISPRArchitect is calibrated to published experimental data. Key references:

- **Gene conversion tracts:** Elliott et al., Mol Cell Biol, 1998; Paquet et al., Nature, 2016
- **cssDNA donors:** Iyer et al., CRISPR Journal, 2022; Xie et al., Nature Biotechnology, 2025
- **Staggered cuts & HDR:** Chauhan et al., PNAS, 2023
- **enFnCas9:** Chakraborty et al., Nature Communications, 2024
- **3D genome & HDR:** Marin-Gonzalez et al., Science, 2025
- **Polymer models:** Lieberman-Aiden et al., Science, 2009; Rao et al., Cell, 2014

## License

MIT License

## Citation

If you use CRISPRArchitect in your research, please cite:
[Publication pending]
