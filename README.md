# CRISPRArchitect

**Transcript-aware, consequence-guided design of genome editing strategies across modalities**

CRISPRArchitect is a computational framework for designing and ranking genome editing strategies that unifies base editing (BE), prime editing (PE), and homology-directed repair (HDR) within a single decision-making system.

## Key Features

- **Unified strategy space**: Evaluates BE, PE, HDR, and hybrid combinations together
- **Transcript-aware mapping**: Maps variants to exon structure, coding frame, and splice sites via Ensembl
- **Consequence-aware scoring**: Integrates amino acid changes, splice proximity, and bystander effects
- **PAM-verified feasibility**: Checks actual PAM availability and editing window compatibility
- **Explicit rejection**: Identifies and explains why infeasible strategies are excluded
- **Multi-nuclease support**: SpCas9 (NGG) and enFnCas9 (NRG) as first-class citizens

## Benchmark Results

Evaluated on 30 curated variant scenarios with verified GRCh38 coordinates:

| Metric | Value |
|--------|-------|
| Top-1 Accuracy | 86.7% (26/30) |
| Top-3 Accuracy | 96.7% (29/30) |
| Rejection Accuracy | 90.0% (27/30) |

## Architecture

CRISPRArchitect comprises two layers:

### v1 Modules (~24,000 LOC)
- **ConversionSim**: Monte Carlo HDR gene conversion tract simulator (validated against 4 published datasets)
- **MOSAIC**: Multi-locus editing strategy optimizer (benchmarked against 14 published papers)
- **TopoPred**: cssDNA secondary structure analyzer (G-quadruplex, hairpin, accessibility)
- **ChromBridge**: 3D chromatin distance and translocation risk predictor
- **LoopSim**: Cohesin loop extrusion simulator
- **WebApp**: Interactive Streamlit application

### v2 Modules (~9,400 LOC)
- **Sequence Layer**: Ensembl transcript fetcher, genomic-to-transcript mapper, reference validator, coding annotator, variant normalizer
- **Feasibility Layer**: Enhanced PAM scanner, base editing engine (ABE/CBE), prime editing engine, HDR design engine
- **Strategy Layer**: Consequence-aware strategy generator, annotation integrator
- **Pipeline**: End-to-end orchestrator (single entry point)
- **Benchmarks**: 30-case evaluation framework with publication-quality plotting

## Quick Start

### Installation

```bash
git clone https://github.com/visvikbharti/CRISPRArchitect.git
cd CRISPRArchitect/crisprarchitect
pip install -r requirements.txt
```

### Run the v2 pipeline

```python
from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput

pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")

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

for strategy in result.strategies:
    print(f"#{strategy.rank}: {strategy.strategy_name} "
          f"(score={strategy.overall_score:.3f})")
```

### Run the benchmark

```bash
python -m benchmarks.run_benchmark --input benchmarks/dataset_v1.json
```

### Run the interactive web app

```bash
streamlit run webapp/app.py
```

### Run tests

```bash
python -m pytest tests/ -v
```

## Project Structure

```
crisprarchitect/
    core/                          # v2 modules
        models.py                  # Central data models
        sequence/                  # Transcript mapping and annotation
        feasibility/               # BE, PE, HDR feasibility engines
        mosaic/                    # Strategy generation and scoring
        pipeline/                  # End-to-end orchestrator
    mosaic/                        # v1 strategy optimizer
    conversion_sim/                # v1 gene conversion simulator
    topopred/                      # v1 cssDNA structure analyzer
    chrombridge/                   # v1 3D chromatin predictor
    loopsim/                       # v1 cohesin loop simulator
    utils/                         # Shared utilities and constants
    webapp/                        # Streamlit interactive app
    benchmarks/                    # Evaluation framework
        dataset_v1.json            # 30 curated ClinVar cases
        evaluator.py               # Benchmark evaluator
        run_benchmark.py           # CLI runner
        plotting.py                # Publication figures
    tests/                         # Test suite
    paper/                         # Manuscript and figures
    validation/                    # v1 validation reports
```

## Scoring Function

```
Score = w1*Safety + w2*Feasibility - w3*Complexity - w4*Risk + w5*Confidence
```

Default weights (iPSC-optimized):
- Safety: 0.30 (DSB-free approaches preferred)
- Feasibility: 0.25 (PAM availability, window compatibility)
- Complexity: 0.20 (rounds, donors, screening burden)
- Risk: 0.15 (rearrangement, bystanders, splice proximity)
- Confidence: 0.10 (evidence tier)

## Key Finding

PAM-dependent editing window constraints are a more significant bottleneck for base editing applicability than mutation-type classification alone. At all seven tested ClinVar loci with ABE-compatible transitions, no SpCas9 guide placed the target base within the ABE editing window (positions 4-7), causing prime editing to emerge as the preferred modality in 97% of cases.

## Requirements

- Python 3.9+
- NumPy >= 1.24.0
- SciPy >= 1.10.0
- Matplotlib >= 3.7.0
- Streamlit >= 1.30.0 (for web app)

## Citation

If you use CRISPRArchitect in your research, please cite:

> Bharti V, Chakraborty D. CRISPRArchitect: transcript-aware and consequence-guided design of genome editing strategies across modalities. (2026). *In preparation.*

## License

MIT License

## Authors

- **Vishal Bharti** — CSIR-Institute of Genomics and Integrative Biology, New Delhi
- **Debojyoti Chakraborty** — CSIR-Institute of Genomics and Integrative Biology, New Delhi
