# CRISPRArchitect — Honest Scientific Assessment

## A researcher-to-researcher evaluation of novelty, gaps, and next steps

---

## 1. Is This App Complete?

**No. It is a working prototype, not a finished product.**

### What IS done:
- Core algorithms for all 5 modules
- Web interface with dark/light mode
- Demo gene functionality
- Literature-backed parameter calibration

### What IS MISSING for a real tool:

| Missing Feature | Why It Matters | Effort |
|----------------|---------------|--------|
| **Ensembl API integration in web app** | Users shouldn't have to manually paste exon coordinates | 1-2 days |
| **Actual sgRNA design** | A real workflow needs guide RNA sequences, not just strategy advice | Would need to integrate CRISPOR/Cas-OFFinder |
| **Real donor sequence generation** | Should auto-generate cssDNA donor sequence with correct genomic homology arms | 2-3 days |
| **Experimental validation** | Without comparing predictions to real experimental outcomes, the numbers are educated guesses | Weeks-months of wet lab |
| **Proper statistical framework** | Current parameter uncertainties are not propagated through the models | 1-2 weeks |
| **Unit tests** | No automated test suite | 1 week |

### Should we add genome data behind it?

**Yes, but minimally.** We should:
- Add Ensembl REST API fetch directly into the web app (type gene name, auto-populate exons)
- Add sequence fetch for the target region (for sgRNA design and donor generation)
- Do NOT download/store the whole genome — the REST API is sufficient

---

## 2. Should We Add an AI Assistant?

**No. Not necessary, and it would hurt more than help.**

### Why NOT:
- The tool's value is in its **quantitative models**, not conversation
- An AI chatbot adds complexity, API costs, and a dependency on external services
- Users (scientists) want numbers and plots, not chat responses
- It would blur the line between "computed result" and "AI opinion" — dangerous for scientific tools
- Reviewers would question whether the outputs are deterministic or stochastic AI hallucinations

### The one exception where AI could help:
- A **report generation** step that takes all the computed results and writes a human-readable summary paragraph — but we already do this with template-based text, which is more trustworthy than LLM output for scientific claims

### Bottom line:
Keep the tool deterministic and transparent. Scientists trust tools where they can trace every number to a formula and a citation.

---

## 3. The Hard Question: Is This Tool Genuinely Novel?

### Honest comparison with existing tools:

| Existing Tool | What It Does | What CRISPRArchitect Adds |
|--------------|-------------|--------------------------|
| **CRISPOR** | sgRNA design, off-target prediction | CRISPOR does NOT compare editing strategies or simulate HDR mechanics |
| **IDT Alt-R HDR Design** | Single-site donor template design | Does NOT handle multi-site scenarios or compare HDR vs base editing |
| **GenScript HDR Designer** | Single-site donor design up to 20 kb | Does NOT score strategies or predict tract lengths |
| **Benchling** | Full CRISPR workflow (sgRNA + donor) | Does NOT simulate HDR, predict 3D distances, or optimize multi-site strategies |
| **PrimeDesign** | Prime editing pegRNA design | Only does prime editing — does NOT compare with HDR or base editing |
| **BE-Designer** | Base editing guide design | Only does base editing — no comparison with other modalities |
| **CRISPRscan** | sgRNA efficiency scoring | No HDR modeling at all |

### What IS genuinely novel:

#### 1. ConversionSim — Gene Conversion Tract Simulator
**NOVEL. Nothing like this exists publicly.**

- No existing tool simulates the HDR process step-by-step (resection → RAD51 → synthesis → SDSA)
- No tool predicts "probability of incorporating an edit at distance X from the cut"
- Labs currently GUESS based on rough literature values ("tracts are usually <2 kb")
- This simulation gives quantitative, locus-aware predictions

**Novelty level: HIGH**
**Publication potential: YES — as a methods paper**

#### 2. MOSAIC — Multi-Site Strategy Optimizer
**NOVEL CONCEPT. The framework is new, but the scoring is heuristic.**

- No existing tool takes two mutations in a gene and compares ALL possible strategies (HDR, base editing, prime editing, exon deletion, hybrids, sequential)
- The concept of a "strategy score" combining efficiency, safety, time, and cost for iPSC editing is new
- BUT: the scoring weights are hand-tuned, not learned from data
- AND: the efficiency estimates are literature averages, not locus-specific predictions

**Novelty level: MODERATE-HIGH (concept is new; implementation needs validation)**
**Publication potential: YES — but needs experimental validation of at least a few strategy recommendations**

#### 3. ChromBridge — 3D Distance Predictor
**INCREMENTAL. The physics is textbook; the application to CRISPR is new-ish.**

- The polymer models (Gaussian chain, WLC) are standard biophysics
- The "bridgeability ratio" concept (can a donor coil span the inter-locus distance?) is a useful pedagogical tool but not a deep contribution
- The translocation risk scoring from Hi-C contact frequency is known (Frock et al., 2015)
- The packaging of this into a CRISPR-focused tool is useful but not deeply novel

**Novelty level: LOW-MODERATE**
**Publication potential: Only as part of a larger platform paper, not standalone**

#### 4. cssDNA-TopoPred — Structure Analyzer
**INCREMENTAL. Applies known algorithms to a new context.**

- G-quadruplex regex scanning is well-established
- Hairpin prediction with nearest-neighbor energy is simplified Mfold/ViennaRNA
- The "accessibility score" for homology arms is a new application concept
- BUT: without experimental validation (do low-accessibility donors actually fail?), it's speculative

**Novelty level: LOW-MODERATE**
**Publication potential: Only as part of platform; the G4/hairpin analysis alone is not new**

#### 5. LoopSim — Cohesin Loop Extrusion Simulator
**POTENTIALLY NOVEL but with a critical caveat.**

- Simulating cohesin loop extrusion from a DSB is new — extends the 2025 Science paper into a predictive tool
- BUT: as we noted in our own analysis, loop extrusion is relevant for ENDOGENOUS donors (sister chromatid), NOT for exogenous donors (cssDNA)
- For the user's actual question (cssDNA-mediated HDR), LoopSim is intellectually interesting but practically irrelevant
- It would be more useful for predicting loss-of-heterozygosity or sister chromatid exchange distances

**Novelty level: MODERATE (algorithm) but LOW relevance to the cssDNA use case**
**Publication potential: YES — but as a standalone tool for studying HR search domains, not as part of the cssDNA story**

---

## 4. What Would Make This ACTUALLY Publishable?

### Option A: Platform/Application Note Paper (Moderate effort)
**Target journals:** NAR Genomics & Bioinformatics, Bioinformatics (Application Note), BMC Bioinformatics

**Requirements:**
1. Integrate real gene fetching in the web app (Ensembl API)
2. Run analysis on 10-20 clinically relevant genes
3. Compare MOSAIC recommendations with what published studies actually did
4. Show that ConversionSim predictions are consistent with published tract length data
5. Deploy web app publicly (Streamlit Cloud or similar)
6. Write clear documentation

**This is achievable in 2-3 months.**

### Option B: Methods Paper Focused on ConversionSim (Higher novelty)
**Target journals:** PLOS Computational Biology, Genome Research, NAR

**Requirements:**
1. Calibrate ConversionSim against multiple published datasets of gene conversion tracts
2. Show it predicts tract lengths more accurately than the current "rule of thumb"
3. Validate with experimental data (collaborate with a lab doing HDR experiments)
4. Include the MOSAIC strategy optimizer as a downstream application
5. Compare with any other quantitative HDR models in the literature

**This would be a stronger paper but needs experimental validation data.**

### Option C: Full Platform Paper (Most ambitious)
**Target journals:** Nature Methods, Genome Biology

**Requirements:**
- Everything in Option B, plus
- Experimental validation of at least 3-5 strategy recommendations
- Demonstration that MOSAIC's top-ranked strategy outperforms naive choices
- ChromBridge validated against actual Hi-C data
- TopoPred validated by comparing donor structures with HDR efficiency data
- LoopSim validated against the Marin-Gonzalez et al. data
- Community beta testing

**This would take 6-12 months and wet-lab collaboration.**

---

## 5. What Should We Do Next?

### Immediate (this week):
1. **Integrate Ensembl API into the web app** — so users type a gene name and get real coordinates
2. **Add auto-donor design** — generate cssDNA donor sequence with correct genomic homology arms
3. **Clean up and test** — make sure all edge cases are handled

### Short-term (1-2 months):
4. **Validate ConversionSim** — find published datasets of gene conversion tract lengths and compare
5. **Benchmark MOSAIC** — check 10-20 published multi-site editing papers and see if MOSAIC recommends the same strategy the authors used
6. **Deploy publicly** — Streamlit Cloud (free for public apps)

### Medium-term (3-6 months):
7. **Write the paper** — likely Option A (application note) first, Option B if validation goes well
8. **Get feedback** — share with 3-5 labs doing iPSC editing and collect user feedback
9. **Add sgRNA design** — integrate with an existing tool (CRISPOR API) rather than building from scratch

### What NOT to do:
- Don't add AI/LLM features — keeps the tool trustworthy
- Don't over-engineer — a simple, well-validated tool beats a complex unvalidated one
- Don't try to replace Benchling/CRISPOR — complement them instead

---

## 6. Final Honest Take

**The CONCEPT is solid and fills a real gap.** No existing tool helps researchers choose between HDR, base editing, prime editing, and hybrid strategies for multi-site editing. That decision is currently made by reading papers and asking colleagues. MOSAIC formalizes this.

**ConversionSim is the most novel component.** A Monte Carlo simulator of HDR gene conversion tracts does not exist publicly. If validated, this alone could be a useful contribution.

**The current IMPLEMENTATION is a good prototype** but not publication-ready. The main weakness is the lack of experimental validation — all the numbers are "educated estimates from literature," not predictions that have been tested against ground truth.

**The honest question to ask yourself:** Would I, as a reviewer, accept a paper describing this tool?

My answer: **Yes for an application note in NAR/Bioinformatics**, with the caveat that the predictions are framed as "literature-calibrated estimates" rather than "validated predictions." **Not yet for a methods paper in a top journal** — that would require experimental validation of at least ConversionSim.

**Is it worth making public?** Yes — even without a paper, releasing it on GitHub with good documentation would be useful to the community. Many bioinformatics tools start as GitHub repos and gain users before being formally published.
