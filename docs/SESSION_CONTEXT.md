# CRISPRArchitect — Session Context for Future Reference

## What This Document Is

This document captures the complete state of the CRISPRArchitect project as of March 20, 2026. It is written so that anyone (including a future AI assistant or the authors themselves) can pick up exactly where we left off.

---

## THE STORY IN 60 SECONDS

Vishal Bharti (PhD student/researcher at CSIR-IGIB, New Delhi) asked a question: "Can a single cssDNA donor correct mutations in both exon 40 and exon 60 of a large gene using enFnCas9 in iPSCs?"

The answer is NO — for quantifiable biological reasons (gene conversion tracts are too short, 3D distances too large, HDR is local). But answering this required building tools that don't exist. So we built CRISPRArchitect — a computational framework with a Monte Carlo HDR simulator (ConversionSim) and a multi-site editing strategy optimizer (MOSAIC). We validated both against published data, wrote a manuscript for PLOS Computational Biology, and published the code on GitHub.

---

## WHAT EXISTS RIGHT NOW

### Code (GitHub: visvikbharti/CRISPRArchitect)

```
crisprarchitect/
├── conversion_sim/    # Monte Carlo HDR tract simulator (6 files)
├── mosaic/            # Strategy optimizer (6 files)
├── topopred/          # cssDNA structure analyzer (5 files)
├── chrombridge/       # 3D distance predictor (5 files)
├── loopsim/           # Cohesin loop extrusion sim (6 files)
├── webapp/            # Streamlit web app (2 files + style)
├── utils/             # Constants, Ensembl API, plotting (4 files)
├── tests/             # 30 unit tests (1 file)
├── validation/        # ConversionSim + MOSAIC validation (4 files + figures)
├── examples/          # Demo scripts (2 files)
├── paper/             # PLOS manuscript + figures + supplementary
├── docs/              # Architecture, user guide, roadmap, Q&A, honest assessment
├── cli.py             # Command-line interface
├── Dockerfile         # Container deployment
├── pyproject.toml     # Package config
└── .github/workflows/ # CI pipeline
```

### Paper (paper/ directory)

| File | What |
|------|------|
| PLOS_CompBio_manuscript.md | Full manuscript in PLOS format |
| CRISPRArchitect_PLOS_manuscript.docx | Word version for submission |
| supplementary_materials.md | S1 Table, S1-S2 Figs, S1-S2 Text |
| CRISPRArchitect_Supplementary.docx | Word version |
| cover_letter.md | Cover letter for PLOS |
| CRISPRArchitect_Cover_Letter.docx | Word version |
| REFERENCE_VERIFICATION.md | All 22 refs verified against PubMed |
| generate_figures.py | Script to regenerate all figures |
| figures/Fig1_ConversionSim.png/.pdf | Model + tracts + survival curve |
| figures/Fig2_Validation.png/.pdf | 4-panel validation |
| figures/Fig3_MOSAIC_Benchmark.png/.pdf | 14-study benchmarking |
| figures/FigS1_Sensitivity.png/.pdf | Parameter sensitivity |
| figures/FigS2_Weight_Sensitivity.png/.pdf | Scoring weight robustness |

### Presentation

CRISPRArchitect_Presentation.pptx (18 slides, in parent directory)
- Conference/DAC format: Background → Gap → Objectives → Methods → Results → Validation → Discussion → Conclusions
- Includes all publication figures embedded
- Live demo guide (6 steps with NF1)

### Documentation

| File | Purpose |
|------|---------|
| docs/PRESENTATION_GUIDE.md | 32 Q&As, slide-by-slide speaker notes, elevator pitch |
| docs/ARCHITECTURE.md | How the code works, for non-CS scientists |
| docs/USER_GUIDE.md | Step-by-step usage (web + CLI + Python) |
| docs/ROADMAP.md | 6-phase development plan |
| docs/HONEST_ASSESSMENT.md | What is/isn't novel, publication strategy |

---

## KEY DECISIONS THAT WERE MADE

### Why these tools?
- **ConversionSim**: No public Monte Carlo HDR tract simulator exists. This is the core novel contribution.
- **MOSAIC**: No tool compares HDR vs base editing vs prime editing for multi-site scenarios. Second novel contribution.
- **ChromBridge/TopoPred/LoopSim**: Supporting utilities, NOT claimed as novel methods.

### Why PLOS Computational Biology?
- Research Article format (no word limit, full methods)
- Values reproducibility and open source
- Publishes both methods and software papers
- Our work sits between the two

### Why Streamlit (not Flask/React)?
- Pragmatic: science is identical regardless of UI framework
- 1,700 lines vs weeks of frontend work
- CLI + Docker + Python API also provided
- UI framework adds zero scientific value

### Why no machine learning?
- Insufficient training data for HDR tract lengths
- Mechanistic model is interpretable — failure at ssODN told us about SSTR
- Better extrapolation to novel conditions

### Why safety weight = 40%?
- iPSCs have active p53 → DSBs trigger apoptosis (Ihry 2018, Haapaniemi 2018)
- Weight is heuristic but robust (base editing stays #1 across 0.30-0.80 range)
- User-adjustable

---

## VALIDATION RESULTS TO REMEMBER

### ConversionSim (4 datasets)

| Test | Predicted | Observed | Source | Match? |
|------|-----------|----------|--------|--------|
| cssDNA/lssDNA ratio | 2.07x | 1.9x (1.5-2.1) | Iyer 2022 CRISPR J | YES |
| Stagger/blunt ratio | 1.82x | 1.9x (1.4-2.8) | Chauhan 2023 PNAS | YES |
| Tract shape | Right-skewed | Right-skewed | Elliott 1998 MCB | Qualitative YES |
| ssODN distance | Gradual decline | Sharp decline <50bp | Paquet 2016 Nature | NO (SSTR gap) |

### MOSAIC (14 papers)

- Overall: 10/14 = 71.4%
- Base editing: 5/5 = 100%
- Prime editing: 1/1 = 100%
- HDR (ssODN): 3/5 = 60%
- 4 disagreements: 3 due to tools not existing when paper was published; 1 genuine gap (DMD exon reframing)

---

## WHAT NEEDS TO HAPPEN NEXT

### Before lab meeting (Wednesday):
- [ ] Review the 18-slide presentation
- [ ] Practice with the presentation guide (32 Q&As)
- [ ] Test the live demo: `cd crisprarchitect && streamlit run webapp/app.py`
- [ ] Fetch NF1, add mutations, run MOSAIC + ConversionSim + ChromBridge

### Before paper submission:
- [ ] Get DC's review and approval of manuscript
- [ ] Fill in CSIR grant numbers in funding statement
- [ ] Add 3 reviewer suggestions to cover letter
- [ ] Final read-through of .docx files
- [ ] Submit via PLOS editorial manager

### Future development (post-submission):
- Add SSTR sub-model for ssODN donors (highest priority)
- Add exon-reframing strategy to MOSAIC
- Integrate locus-specific chromatin features (ENCODE)
- PAM availability check in MOSAIC
- Bayesian parameter estimation
- Prospective experimental validation

### When colleague shares the real gene:
- Just need: gene symbol, mutation 1 (exon, ref>alt), mutation 2 (exon, ref>alt)
- Run: `python cli.py analyze --gene SYMBOL --exon1 X --mut1 R>A --exon2 Y --mut2 R>A`
- Or use web app: Gene Setup → Fetch from Ensembl → Add mutations → Run all modules

---

## IMPORTANT TECHNICAL DETAILS

### Python version: 3.9 (Anaconda)
- No `int | float` union syntax (use Optional from typing)
- No `list[str]` (use List[str] from typing)
- This was fixed once already — don't reintroduce

### Git config
- Author: Vishal Bharti <vishalbharti@Vishals-MacBook-Air-2.local>
- GitHub shows this correctly now
- Co-Authored-By lines were removed via filter-branch + force push

### Web app
- Launch: `cd crisprarchitect && streamlit run webapp/app.py`
- Dark mode toggle in sidebar
- Ensembl fetch on Gene Setup page (requires internet)
- Auto donor design on TopoPred page (requires Ensembl gene to be fetched first)

### Running simulations
```python
# Quick test
from conversion_sim import ConversionSimulator
sim = ConversionSimulator(cut_type='staggered_5prime', overhang_length=3,
    donor_topology='circular_ssDNA', homology_arm_length=300,
    cell_type='iPSC', n_simulations=10000)
sim.run()
sim.summary()
```

### Running tests
```bash
cd crisprarchitect
python -m pytest tests/ -v
```

---

## WHAT NOT TO DO

- Don't claim ChromBridge/TopoPred/LoopSim are novel methods (they're utilities)
- Don't claim absolute locus-specific HDR prediction
- Don't claim MOSAIC replaces expert judgment (it supplements it)
- Don't use Python 3.10+ syntax
- Don't add AI/LLM features to the tool
- Don't over-engineer the UI
- Don't fabricate citations or data — ever
