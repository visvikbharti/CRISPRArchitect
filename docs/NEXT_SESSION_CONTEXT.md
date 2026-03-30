# CRISPRArchitect — Session Context & Continuation Guide

**Last updated:** 2026-03-31
**Version:** 3.0.0
**Status:** Lab meeting ready (Wednesday). Manuscript rebuild pending.
**Repository:** github.com/visvikbharti/CRISPRArchitect (branch: main, clean)
**Tests:** 224 passing, 0 failures

---

## What This Document Is

This is the single source of truth for continuing CRISPRArchitect development.
It records everything that was done, what remains, what decisions were made and
why, and the exact state of every component. Read this before making any changes.

---

## Project Summary

CRISPRArchitect is a computational decision-support tool that recommends the
optimal genome editing strategy (base editing vs. prime editing vs. HDR) for
any given pathogenic variant. It evaluates 5 nucleases (SpCas9, enFnCas9,
SpCas9-NG, SpRY, Cas12a) paired with 9 base editor profiles, uses 6-dimensional
TOPSIS multi-criteria decision analysis with Pareto front analysis and Monte
Carlo sensitivity quantification, and provides transparent, uncertainty-
quantified recommendations. Developed by Vishal Bharti under PI Debojyoti
Chakraborty at CSIR-IGIB, New Delhi. Target journal: PLOS Computational
Biology.

---

## Complete Session Log (2026-03-30)

### All Changes Made (in order)

1. **Codebase cleanup**: Unified version to 3.0.0 (5 files), deleted 7 dead files, updated CLI to v3, deprecated v1 MOSAIC
2. **Parameter audit**: All 10 critical params tagged [MEASURED]/[DERIVED]/[ASSUMED] with evidence chains
3. **Monte Carlo rigor**: SEs/CIs on all outputs, Wilson intervals, ddof=1, SeedSequence, vectorized Dirichlet
4. **References verified**: 20 web-search verified, 3 corrected (Arbab: Cell not Nature; enFnCas9: Acharya first author; Walton: pages)
5. **Bystander triple-counting bug FIXED**: Root cause of degenerate PE pattern. Bystander now in consequence dimension only.
6. **6D TOPSIS**: Consequence as proper 6th dimension (not post-hoc additive)
7. **Pareto front analysis**: Weight-independent dominance check
8. **VIKOR + WPM**: Cross-method comparison (100% rank concordance)
9. **Safety-Risk independence**: Verified after bystander removal
10. **ConversionSim scope restricted**: SDSA only, ssODN/SSTR explicitly excluded
11. **Citation integrity audit**: Elliott 1998 tracts <58bp (not 200-2000bp), Kan 2017 journal corrected, HA 300nt unsupported by Iyer, Kim 2019→Song 2020
12. **Literature benchmark**: 33 cases built (20 genes), 10 run through pipeline (30% top-1, 80% top-3)
13. **6 figures generated**: All from real data, embedded in PPTX
14. **28-slide presentation**: Dark theme, coral banners, embedded figures (includes SDSA sensitivity + delivery advisor slides)
15. **Speaker guide**: 27 Q&A pairs, figure map, benchmark interpretation
16. **Meeting prep**: Demo scripts, PI anticipation, contingency plans
17. **Documentation**: 1,449-line comprehensive doc (16 sections, 61 evidence tags)
18. **File cleanup**: Deleted 24 stale figures, consolidated directories, fixed all remaining citation errors
19. **PI data request**: Formal letter + 2 TSV collection sheets
20. **SDSA p sensitivity analysis**: p=0.002 robustness confirmed across 0.001-0.005 range (Fig_SDSA_Sensitivity.png)
21. **Delivery-aware post-ranking annotations**: Option B architecture (hard filters + post-ranking annotations, NOT a TOPSIS dimension). 87 verified references, 34 tests. File: `core/feasibility/delivery_advisor.py`

---

## Key Technical Decisions

| Decision | Rationale |
|----------|-----------|
| TOPSIS 6D over weighted sum | Penalizes extremes; consequence as proper dimension |
| Bystander in consequence only | Prevents triple-counting that made PE unbeatable |
| SDSA scope restriction | ssODN uses SSTR pathway; wrong model gives R²=-0.56 |
| SDSA p=0.002 as [ASSUMED] | No direct measurement exists; functional evidence only |
| enFnCas9 stagger = 3 bp [ASSUMED] | Not published; inferred from crystal structure |
| PLOS Comp Bio target | Nature Methods requires experimental validation |
| Dirichlet conc=20, min_alpha=2.0 | Prevents pathological weights from small dimensions |
| ABE8e window 3-9 (extended) | Canonical 4-8; we use extended to maximize rescue |

---

## Known Limitations (14 items)

1. No experimental validation — self-curated + literature labels only
2. SDSA p=0.002 [ASSUMED] — no direct measurement
3. enFnCas9 stagger 3 bp [ASSUMED] — not published
4. enFnCas9 HDR multiplier 1.5x [ASSUMED]
5. ConversionSim invalid for ssODN (SSTR pathway)
6. Off-target scoring local only (no genome-wide)
7. No chromatin accessibility (ATAC-seq)
8. PE dominates 23/30 (genuine biology)
9. Literature benchmark 30% top-1 (explainable)
10. Large deletion handling 0/3 accuracy
11. ABE8e-enFnCas9 window [ASSUMED]
12. SpCas9-NRCH not in nuclease set
13. Self-curated truth labels need independent expert validation
14. HA 300 nt [ASSUMED]

---

## TODO List (Prioritized for Next Sessions)

### Immediate (Before Submission)
- [ ] Rebuild manuscript for PLOS Comp Bio (Methods, Results, Supplementary)
- [ ] Expand runnable benchmark cases (resolve HGVS→GRCh38 for 22 skipped cases)
- [ ] Independent expert labeling (3+ PIs, Cohen's kappa)
- [ ] Integrate lab data when available
- [ ] Write PLOS Comp Bio cover letter
- [x] Update README.md for v3 (done 2026-03-30)

### Strengthening
- [ ] Validate BE against Song et al. 2020 (13,504 ABE targets)
- [ ] Validate PE against Chen et al. 2021 (iPSC data)
- [ ] Add SpCas9-NRCH nuclease
- [ ] Update Streamlit webapp for v3
- [ ] Run VIKOR/WPM on all 30 ClinVar cases
- [x] SDSA p sensitivity analysis (0.001-0.005) (done 2026-03-31)
- [x] Delivery-aware post-ranking annotations (87 verified refs, 34 tests) (done 2026-03-31)
- [ ] Cell-type-specific modality priors

### Future
- [ ] SSTR sub-model for ssODN
- [ ] Genome-wide off-target (Cas-OFFinder)
- [ ] ATAC-seq integration
- [ ] SINGLE_CUT_EXON_REFRAMING for DMD

---

## Key Numbers (All Verified)

| Metric | Value |
|--------|-------|
| Version | 3.0.0 |
| Tests | 224 pass, 0 fail |
| TOPSIS dimensions | 6 |
| Nucleases | 5 |
| Base editors | 9 (3 Tier A + 6 Tier B) |
| ClinVar benchmark | 30 cases, 86.7% top-1 |
| Literature benchmark | 33 cases, 30% top-1, 80% top-3 |
| v3 strategy distribution | PE=23, BE=6, HDR=1 |
| Cross-method concordance | 100% |
| Manuscript references | 20 (all verified with PMIDs) |
| Sensitivity runs | 10,000 Dirichlet permutations |
| SDSA displacement prob | 0.002 [ASSUMED] |
| Delivery advisor refs | 87 verified |
| Delivery advisor tests | 34 |

### New Files (2026-03-31)

| File | Description |
|------|-------------|
| `core/feasibility/delivery_advisor.py` | Delivery-aware post-ranking annotations (Option B architecture) |
| `tests/test_delivery_advisor.py` | 34 tests for delivery advisor module |
| `paper/figures/v3_results/Fig_SDSA_Sensitivity.png` | SDSA p sensitivity analysis figure (p=0.001-0.005) |

### Parent Directory Files (Delivery Literature Review)

| File | Description |
|------|-------------|
| `../DELIVERY_METHODS_COMPREHENSIVE_LITERATURE_REVIEW.md` | Comprehensive delivery methods literature review (87 verified references) |

---

## Prompt for Next Session

Copy and paste this as your first message:

```
I'm continuing work on CRISPRArchitect v3. Please read these files for context:
1. crisprarchitect/docs/NEXT_SESSION_CONTEXT.md (project state and TODOs)
2. crisprarchitect/docs/COMPLETE_PROJECT_DOCUMENTATION.md (full technical details)

Key facts:
- Version 3.0.0, 224 tests passing, all pushed to GitHub
- Lab meeting was on Wednesday — [DESCRIBE: how it went, PI feedback, data decisions]
- 6D TOPSIS + Pareto + VIKOR/WPM scoring (bystander bug fixed)
- 33-case literature benchmark, 10 run through pipeline
- All citations web-verified, all parameters evidence-tagged
- Target: PLOS Computational Biology

Today I want to: [CHOOSE ONE OR MORE]
1. Rebuild the manuscript for PLOS Comp Bio submission
2. Expand and run more benchmark cases
3. Update the Streamlit webapp for v3
4. Integrate experimental data from the lab [ATTACH DATA IF AVAILABLE]
5. Address reviewer feedback [IF APPLICABLE]
6. Other: [DESCRIBE]
```

---

*Last verified: 2026-03-31*
*GitHub: github.com/visvikbharti/CRISPRArchitect (main, clean)*
