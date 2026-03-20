# CRISPRArchitect — Roadmap to a Legacy Tool

## Vision Statement

CRISPRArchitect will be the go-to computational platform that any researcher can use to
answer: **"I have multiple mutations in a gene — what is the safest, most efficient
strategy to correct them in my cells?"**

No existing tool answers this question. We will be the first.

---

## Phase 1: FOUNDATION (Week 1-2) — Make It Real

### 1.1 Ensembl API Integration in Web App
**Priority: CRITICAL — this is the #1 barrier to usability**

- [ ] Add "Fetch from Ensembl" option on Gene Setup page
- [ ] User types gene symbol (e.g., NF1, DMD, BRCA2) → auto-populate exon coordinates
- [ ] Show gene metadata: chromosome, strand, span, transcript ID
- [ ] Cache fetched genes in session state
- [ ] Handle errors gracefully (gene not found, API timeout)
- [ ] Support custom transcript selection (not just canonical)

### 1.2 Auto Donor Sequence Generation
**Priority: HIGH — researchers need actual sequences, not just strategy advice**

- [ ] Fetch genomic sequence for target region from Ensembl REST API
- [ ] Auto-generate cssDNA donor with:
  - Correct genomic homology arms (including intronic sequence flanking the exon)
  - Corrected exon sequence (with the mutation fixed)
  - Optional silent PAM-blocking mutations (Paquet et al., 2016)
  - Configurable arm lengths (default 300 bp)
- [ ] Output: downloadable FASTA of the donor template
- [ ] Feed generated sequence directly into TopoPred for quality check

### 1.3 sgRNA Suggestion Engine
**Priority: MEDIUM — complement, don't replace, existing tools**

- [ ] Scan target region for PAM sites (NGG for SpCas9, NRG for enFnCas9, TTTV for Cas12a)
- [ ] Rank by proximity to mutation site (closer = better for HDR)
- [ ] Flag guides that cut within the donor homology arms
- [ ] Link out to CRISPOR for full off-target analysis (we don't replicate that)
- [ ] Display cut position relative to the mutation

### 1.4 Code Quality
- [ ] Add `__init__.py` for the top-level crisprarchitect package
- [ ] Add `setup.py` or `pyproject.toml` for pip-installable package
- [ ] Write unit tests for each module (pytest)
- [ ] Add input validation with clear error messages
- [ ] Type checking with mypy (at least for public APIs)

---

## Phase 2: VALIDATION (Week 3-6) — Prove It Works

### 2.1 ConversionSim Validation
**Priority: CRITICAL — this is the paper's backbone**

**Published datasets to validate against:**

- [ ] Elliott et al., MCB, 1998 — original mammalian gene conversion tract measurements
  - Compare our simulated tract distribution with their observed distribution
  - They used I-SceI endonuclease in mouse cells

- [ ] Paquet et al., Nature, 2016 — distance-dependent mutation incorporation
  - They measured how incorporation frequency decreases with distance from cut
  - Our "probability_at_distance()" function should match their curve

- [ ] Iyer et al., CRISPR Journal, 2022 — cssDNA vs lssDNA HDR rates
  - Compare predicted HDR enhancement from circular topology with their observed 2x improvement

- [ ] Chauhan et al., PNAS, 2023 — staggered cut HDR enhancement
  - Compare predicted effect of 5' overhangs with their vCas9 data (1.9x improvement)

**Validation approach:**
- [ ] Create `validation/` directory with scripts for each comparison
- [ ] Generate comparison plots: predicted vs observed
- [ ] Calculate goodness-of-fit metrics (R², RMSE)
- [ ] Document any parameter adjustments needed to match data
- [ ] Write up as supplementary validation in the paper

### 2.2 MOSAIC Benchmarking
**Priority: HIGH — proves the strategy optimizer makes correct recommendations**

**Approach:** Find 15-20 published papers where researchers corrected multiple mutations
in iPSCs and check if MOSAIC would have recommended the same strategy.

Papers to benchmark against:

- [ ] Young et al., Cell Stem Cell, 2016 — DMD exon 45-55 deletion (NHEJ)
- [ ] Paquet et al., Nature, 2016 — Alzheimer's mutations in APP (ssODN HDR)
- [ ] Mandal et al., Cell Stem Cell, 2014 — B2M + CCR5 dual knockout (NHEJ)
- [ ] Jackow et al., Mol Therapy, 2019 — COL7A1 correction (HDR)
- [ ] Roth et al., Nature, 2018 — TCR replacement in T cells (AAV HDR)
- [ ] 10-15 additional papers from PubMed search

**For each paper:**
- [ ] Extract: gene, mutations, cell type, strategy used, efficiency achieved
- [ ] Run MOSAIC with same inputs
- [ ] Record: did MOSAIC's #1 recommendation match what they did?
- [ ] If not: was their choice actually better? Or was MOSAIC right to suggest differently?

**Success criterion:** MOSAIC's top-3 recommendations should include the published strategy
in ≥80% of cases. If not, adjust scoring weights.

### 2.3 ChromBridge Validation
- [ ] Compare predicted 3D distances with published FISH data (Yokota et al., 1995; Mateos-Langerak et al., 2009)
- [ ] Compare translocation risk predictions with HTGTS data (Frock et al., 2015)
- [ ] Show that bridgeability predictions are consistent with known HDR donor size limits

### 2.4 TopoPred Validation
- [ ] Find papers reporting HDR efficiency differences between donors at same locus
- [ ] Run TopoPred on the donor sequences
- [ ] Check: do lower-accessibility donors have lower HDR efficiency?
- [ ] This is harder to validate — may be best presented as "hypothesis-generating"

---

## Phase 3: REFINEMENT (Week 5-8) — Make It Robust

### 3.1 Improved Parameter Estimation
- [ ] Replace fixed constants with uncertainty ranges
- [ ] Run ConversionSim with parameter sampling (sensitivity analysis)
- [ ] Report confidence intervals on all predictions, not just point estimates
- [ ] Add cell-type-specific parameter tuning based on published HDR data

### 3.2 Enhanced MOSAIC Strategies
- [ ] Add "twinPE" (twin prime editing for larger edits)
- [ ] Add "PASTE" (programmable addition via site-specific targeting elements)
- [ ] Add "HITI" (homology-independent targeted insertion) for non-dividing cells
- [ ] Weight strategies differently based on edit size (not just mutation type)
- [ ] Account for PAM availability at each mutation site

### 3.3 Web App Polish
- [ ] Add progress bars for long computations
- [ ] Add "Save/Load Session" (export/import analysis as JSON)
- [ ] Add comparative mode: run same mutations with different nucleases side-by-side
- [ ] Add batch mode: analyze multiple genes from a CSV
- [ ] Responsive design for tablet/mobile
- [ ] Add "Example Analyses" page with pre-loaded results for NF1, DMD, BRCA2

### 3.4 Visualization Improvements
- [ ] Interactive gene structure diagram (click exon → show details)
- [ ] Animated ConversionSim showing tract growth in real-time
- [ ] 3D chromosome territory visualization for ChromBridge
- [ ] Side-by-side strategy comparison radar charts
- [ ] Exportable publication-quality figures (SVG/PDF)

---

## Phase 4: DEPLOYMENT (Week 7-10) — Make It Public

### 4.1 GitHub Repository
- [ ] Create public GitHub repo: `github.com/[username]/CRISPRArchitect`
- [ ] Clean commit history
- [ ] Add LICENSE (MIT)
- [ ] Add CONTRIBUTING.md
- [ ] Add CITATION.cff
- [ ] Add GitHub Actions CI (pytest on push)
- [ ] Add badges (build status, Python version, license)

### 4.2 Streamlit Cloud Deployment
- [ ] Deploy to Streamlit Community Cloud (free, public URL)
- [ ] URL: `crisprarchitect.streamlit.app` (or similar)
- [ ] Add to Streamlit gallery
- [ ] Monitor usage analytics

### 4.3 PyPI Package
- [ ] Make pip-installable: `pip install crisprarchitect`
- [ ] CLI entry point: `crisprarchitect --gene NF1 --mutations "exon20:G>A,exon50:C>T"`
- [ ] API documentation (Sphinx or MkDocs)

### 4.4 Documentation Website
- [ ] Host docs on Read the Docs or GitHub Pages
- [ ] API reference (auto-generated from docstrings)
- [ ] Tutorials with Jupyter notebooks
- [ ] FAQ

---

## Phase 5: PUBLICATION (Week 8-14) — Make It Citable

### 5.1 Paper Strategy

**Primary paper: Application Note**
- **Journal:** NAR Genomics and Bioinformatics, or Bioinformatics (Application Note)
- **Title:** "CRISPRArchitect: a computational platform for optimizing multi-site genome editing strategies"
- **Structure:**
  1. Introduction: the multi-site editing strategy problem
  2. Methods: ConversionSim, MOSAIC, ChromBridge, TopoPred (brief)
  3. Results: validation against published data, benchmarking against published strategies
  4. Case studies: NF1, DMD, BRCA2 with real coordinates
  5. Web server description
  6. Discussion

**Follow-up paper (if validation is strong): Methods Paper**
- **Journal:** PLOS Computational Biology or Genome Research
- **Focus:** ConversionSim — a Monte Carlo framework for predicting HDR outcomes
- **Deeper treatment of the simulation model, parameter sensitivity, validation**

### 5.2 Paper Writing Timeline
- [ ] Week 8-9: Draft introduction and methods
- [ ] Week 9-10: Generate all validation figures
- [ ] Week 10-11: Write results and case studies
- [ ] Week 11-12: Write discussion
- [ ] Week 12-13: Internal review, revise
- [ ] Week 13-14: Submit

### 5.3 Preprint
- [ ] Post to bioRxiv simultaneously with journal submission
- [ ] Share on Twitter/X, LinkedIn, Reddit r/bioinformatics
- [ ] Email to 5-10 labs doing iPSC editing for feedback

---

## Phase 6: COMMUNITY (Ongoing) — Make It a Legacy

### 6.1 User Feedback
- [ ] GitHub Issues for bug reports and feature requests
- [ ] Google Form for user survey
- [ ] Track which features are most used

### 6.2 Collaborations
- [ ] Partner with 2-3 wet-lab groups to validate predictions experimentally
- [ ] If ConversionSim predictions match real data → co-authored validation paper
- [ ] Offer as a resource to iPSC core facilities

### 6.3 Maintenance
- [ ] Update parameters as new data is published
- [ ] Add new editing modalities as they emerge (e.g., MAVE-based editors, retrons)
- [ ] Keep Ensembl API compatibility current
- [ ] Respond to community issues

### 6.4 Long-term Vision
- [ ] Integration with laboratory information management systems (LIMS)
- [ ] API for programmatic access by other tools
- [ ] Machine learning model trained on accumulated user feedback and outcomes
- [ ] Mobile-friendly version for quick lookups at the bench

---

## Success Metrics

| Metric | Target (6 months) | Target (1 year) |
|--------|-------------------|-----------------|
| GitHub stars | 50+ | 200+ |
| Monthly active users (web) | 100+ | 500+ |
| Citations | Preprint posted | 5-10 |
| Validated predictions | ConversionSim against 3+ datasets | All modules validated |
| Community contributions | 2-3 external issues/PRs | Active contributor base |

---

## What NOT To Build

To stay focused, we explicitly will NOT:

- Build an AI/LLM chatbot interface
- Replicate CRISPOR's off-target prediction (link out instead)
- Build a LIMS or sample tracker
- Support non-CRISPR editing (meganucleases, ZFNs, TALENs)
- Build a sequence editor/viewer (Benchling, SnapGene already do this)
- Over-optimize the UI before validating the science

---

## Priority Matrix

```
                    HIGH IMPACT
                        |
    Ensembl API     ConversionSim
    Integration     Validation
    (Phase 1)       (Phase 2)
                        |
  LOW EFFORT -----+------+----- HIGH EFFORT
                        |
    TopoPred        LoopSim
    Validation      Refinement
    (Phase 2)       (Phase 3)
                        |
                    LOW IMPACT
```

**Do first:** Top-left quadrant (high impact, low effort)
**Do second:** Top-right quadrant (high impact, high effort)
**Do later:** Bottom-left (low effort but lower impact)
**Deprioritize:** Bottom-right (high effort, low impact)
