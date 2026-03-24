# CRISPRArchitect — Next Session Context

## How to use this document

Copy the prompt at the bottom of this document and paste it as your first message in the next Claude Code session. It contains everything Claude needs to pick up where we left off.

---

## Project State Summary (as of March 25, 2026)

### What's Built and Working
- **v1** (24,000 LOC): ConversionSim, MOSAIC, TopoPred, ChromBridge, LoopSim, WebApp — validated, PLOS-ready
- **v2** (11,002 LOC): Transcript-aware pipeline with Ensembl integration, PAM-verified feasibility (BE/PE/HDR), consequence-aware scoring, 30-case benchmark
- **Results**: Top-1=86.7%, Top-3=96.7%, Rejection=90.0% on 30 real ClinVar cases
- **Key Finding**: PAM-window constraints are a bigger bottleneck for base editing than mutation-type classification
- **All code committed to GitHub**: https://github.com/visvikbharti/CRISPRArchitect (branch: main)
- **123 tests passing**, zero v1 regressions

### Key Files to Read First
1. `docs/COMPLETE_PROJECT_DOCUMENTATION.md` — comprehensive project documentation (1,586 lines)
2. `README.md` — project overview with quick start
3. `benchmark_results/definitive_benchmark_results.json` — real benchmark results
4. `paper/CRISPRArchitect_v2_manuscript.md` — the Nature Methods manuscript

### Architecture
```
crisprarchitect/
  core/                    # v2 pipeline (18 modules)
    models.py              # 22 dataclasses, 6 enums
    sequence/              # Ensembl fetcher, transcript mapper, validator, annotator, normalizer
    feasibility/           # PAM scanner, BE engine, PE engine, HDR engine
    mosaic/                # Strategy generator, annotation integrator
    pipeline/              # StrategyPipeline orchestrator + StrategyScorer
  benchmarks/              # 30-case evaluation framework
  mosaic/                  # v1 strategy optimizer
  conversion_sim/          # v1 HDR simulator
  topopred/                # v1 cssDNA structure analyzer
  chrombridge/             # v1 3D chromatin predictor
  loopsim/                 # v1 cohesin loop simulator
  webapp/                  # Streamlit app (v1 + v2 tab)
  tests/                   # 5 test files, 123 tests
  paper/                   # Manuscript, PPTX, cover letter
  docs/                    # 3 docs (project doc, speaker guide, user guide)
```

### Entry Point for v2 Pipeline
```python
from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput

pipeline = StrategyPipeline(cell_type="iPSC", nuclease="SpCas9")
result = pipeline.run([
    GenomicVariantInput("17", 31200443, "C", "T", gene_symbol="NF1", name="c.910C>T")
])
```

### Verified Numbers (never fabricate)
- Top-1: 86.7% (26/30), Top-3: 96.7% (29/30), Rejection: 90.0% (27/30)
- PE top-ranked: 29/30, HDR: 1/30, BE: 0/30
- Consequence shift: 0% (96.4% had adjustments applied, 0% changed ranking)
- v1 validation: cssDNA 2.07x (Iyer: 1.9x), stagger 1.82x (Chauhan: 1.9x), MOSAIC 71% concordance
- Scoring weights: Safety=0.30, Feasibility=0.25, Complexity=0.20, Risk=0.15, Confidence=0.10

---

## What Needs To Be Done Next

### Priority 1: Enhance BE applicability
The biggest issue: BE never tops the ranking because PAM-window constraints prevent it. Options:
- Expand the editing window to include positions 3-9 (newer ABE8e variants have broader windows per Richter et al. 2020)
- Add enFnCas9 NRG PAM as default second check (broader PAM = more BE opportunities)
- Consider both sense and antisense strand guides more thoroughly

### Priority 2: HGVS parser
Allow clinical-format input like `NM_000267.3:c.910C>T` instead of requiring chromosome/position/ref/alt separately.

### Priority 3: Better large deletion handling
The 3 top-1 misses are all large deletions. The pipeline currently treats them as point variants. Need:
- Exon-deletion detection (if deletion spans entire exon)
- Exon-skipping strategy suggestion (antisense oligonucleotides as alternative)
- Multi-exon donor design for HDR

### Priority 4: Off-target integration
Add off-target scoring using seed sequence matching or Cas-OFFinder integration.

### Priority 5: Chromatin context
Integrate ATAC-seq or DNase-seq data to score guide accessibility in specific cell types.

### Priority 6: Full CDS annotation
Replace the current spliced-exon surrogate with proper Ensembl CDS start/end tracking for accurate UTR vs coding distinction.

### Priority 7: ML-based scoring
Train a model on published editing outcomes to refine the scoring weights beyond the current literature-based defaults.

---

## Prompt for Next Session

Copy everything below this line and paste as your first message:

---

```
I'm continuing work on CRISPRArchitect v2. Please read these files first to get full context:

1. docs/NEXT_SESSION_CONTEXT.md — this context document
2. docs/COMPLETE_PROJECT_DOCUMENTATION.md — comprehensive project docs
3. README.md — project overview

The project is at /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect/

Quick summary: CRISPRArchitect is a computational framework for genome editing strategy design that unifies base editing, prime editing, and HDR. v1 (24K LOC) is validated. v2 (11K LOC) adds transcript-aware mapping, PAM-verified feasibility, and consequence-aware scoring. Benchmark: 86.7% Top-1, 96.7% Top-3, 90.0% Rejection on 30 real ClinVar cases.

Key finding: PAM-window constraints are a bigger bottleneck for base editing than mutation-type classification — at all 7 tested ClinVar loci with ABE-compatible transitions, no SpCas9 guide placed the target at ABE positions 4-7.

Today I want to work on the remaining enhancements listed in docs/NEXT_SESSION_CONTEXT.md. Let's start with [SPECIFY WHICH PRIORITY].

Important constraints:
- Python 3.9 compatibility
- No data fabrication — all numbers must come from real pipeline output
- v2 code lives in core/ and must not break v1
- Scientific integrity: honest limitations, humble claims
- My PI is Debojyoti Chakraborty (developer of enFnCas9) at CSIR-IGIB, New Delhi
```
