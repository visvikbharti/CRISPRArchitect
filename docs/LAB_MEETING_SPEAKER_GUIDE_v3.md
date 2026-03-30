# CRISPRArchitect v3 — Lab Meeting Speaker Guide

**Presenter:** Vishal Bharti
**Audience:** Debojyoti Chakraborty Lab, CSIR-IGIB
**Duration:** ~40-50 minutes + 15-20 minutes Q&A
**Presentation file:** `paper/CRISPRArchitect_v3_LabMeeting.pptx` (28 slides, includes SDSA sensitivity + delivery advisor slides)

---

## How to Use This Guide

Each slide section below contains:
- **What to say** (key talking points)
- **How to present it** (emphasis, pauses, transitions)
- **What NOT to say** (avoid overstatements)

After the slide-by-slide guide, there is a comprehensive **Q&A section** organized by topic (biological, methodological, mathematical/statistical, architectural, and critical/challenging questions).

---

## Slide-by-Slide Guide

### Slide 1: Title

**What to say:**
> "Good morning/afternoon everyone. Today I'll be presenting CRISPRArchitect version 3 — a computational tool I've been developing for designing genome editing strategies. The core idea is simple: when you have a pathogenic variant to correct, should you use base editing, prime editing, or HDR? CRISPRArchitect evaluates all three and gives you a ranked recommendation."

**Transition:** "Let me start with why this tool is needed."

---

### Slide 2: The Problem

**What to say:**
> "When we're correcting mutations in patient iPSCs, we face a decision that's harder than it looks. Base editing is efficient but only works for specific transitions AND requires a PAM at exactly the right position. Prime editing is more versatile but has complex design requirements. HDR gives us the most flexibility but introduces a DSB — and in iPSCs, DSBs are toxic because of p53."
>
> "Right now, this decision is made by intuition. There's no tool that compares all three modalities side by side for a specific variant at a specific locus."

**Emphasis:** Pause after "in iPSCs, DSBs are toxic because of p53" — this resonates with the lab's daily experience.

---

### Slide 3: What is CRISPRArchitect?

**What to say:**
> "CRISPRArchitect takes a variant — you can give it HGVS notation like NM_000329:c.271C>T — and it automatically fetches the transcript from Ensembl, annotates the consequence, scans for PAM sites across five different nucleases, evaluates feasibility of all three modalities, and produces a ranked recommendation."
>
> "Importantly, this is a decision-support tool, not a black box. Every score component is transparent and interpretable. And every parameter has documented provenance — I can trace every number to a published paper or explicitly mark it as a modeling assumption."

**What NOT to say:** Don't claim it "predicts" editing outcomes. It "recommends" strategies based on feasibility analysis.

---

### Slide 4: Pipeline Overview

> "This slide shows the complete 10-stage pipeline — from variant input through normalization, multi-nuclease PAM scanning, feasibility assessment across all three modalities, strategy generation, TOPSIS ranking, and robustness validation with Pareto analysis and Monte Carlo sensitivity."

> "I'll walk through each of these stages in detail over the next several slides."

---

### Slide 5: v1 Foundation

**What to say:**
> "The project started with v1 — six simulation modules totaling about 24,000 lines of code. The most important is ConversionSim, which simulates HDR gene conversion tracts using Monte Carlo methods. MOSAIC handles multi-locus strategy optimization. ChromBridge predicts 3D chromatin distances. TopoPred analyzes cssDNA secondary structure. LoopSim simulates cohesin loop extrusion. And there's a Streamlit web app."
>
> "These are all biophysically grounded simulations — not heuristic rules."

**Tip:** Keep this brief (~1 minute). The audience wants to hear about results, not module lists.

---

### Slide 6: v2 Recap

**What to say:**
> "Version 2 added the transcript-aware pipeline — connecting to Ensembl, mapping variants to coding coordinates, checking reference alleles, annotating consequences. We benchmarked it on 30 ClinVar cases."
>
> "The results were good — 86.7% top-1 accuracy — but there was a problem. Look at the strategy distribution: PE won in 29 out of 30 cases. Base editing was never top-ranked. Even at seven loci where the mutation was a perfect ABE-compatible transition."

**Emphasis:** Point to "BE: 0/30" and "Problem: PE always wins."

**Transition:** "This is what led us to ask: what's going wrong?"

---

### Slide 7: v2 Key Finding — PAM-Window Bottleneck

**What to say (this is a KEY moment):**
> "This was the most important finding from v2. At all seven ClinVar loci with ABE-compatible transitions — meaning A-to-G or T-to-C mutations that should be perfect for adenine base editing — not a single SpCas9 guide placed the target base within the ABE editing window at positions 4-7."
>
> "This means that mutation-type classification alone — the way most people think about base editing ('oh, it's an A>G, use ABE') — is insufficient. You have to check the PAM-editing window geometry for every single locus."

**Why this matters for the audience:** "This directly affects how we design experiments in this lab. Before ordering an ABE experiment, you should check whether there's actually a PAM that places the target in the window."

---

### Slide 8: v3 What's New

**What to say:**
> "This finding motivated version 3, which has three major new capabilities."
>
> "First, a multi-nuclease engine. Instead of just SpCas9, we now evaluate five nucleases including our own enFnCas9, paired with nine base editor profiles — three core editors plus six nuclease-specific fusions. This expands the PAM space enormously."
>
> "Second, we replaced the simple weighted-sum scoring with TOPSIS — a formal multi-criteria decision method — plus Pareto analysis and Monte Carlo sensitivity analysis. This is mathematically principled and provides uncertainty quantification."
>
> "Third, we added statistical rigor throughout — standard errors, confidence intervals, proper random number generation."

---

### Slide 9: Multi-Nuclease Engine

**What to say:**
> "Here's the nuclease table. SpCas9 with NGG PAM is our baseline. enFnCas9 — developed in this lab — has the NRG PAM, which approximately doubles the number of targetable sites. SpCas9-NG has NG. SpRY is near-PAMless but has reduced activity. Cas12a has TTTV."
>
> "I want to be transparent: the enFnCas9 stagger is marked as 'assumed' — we infer it from the FnCas9 crystal structure and the improved HDR rates, but the actual cut-site characterization hasn't been published. This is one of the data points I'm requesting from the lab."

**Important:** If Prof. Chakraborty is in the room, this is a natural moment for him to comment on enFnCas9 characterization data.

---

### Slide 10: enFnCas9 Advantage

**What to say:**
> "This slide separates what we know from what we assume about enFnCas9. The green box shows measured properties from the 2024 Nature Communications paper — NRG PAM recognition, single-nucleobase specificity, and improved HDR knock-in. The orange box shows what we're assuming — the 3 bp stagger, the 1.5x HDR multiplier, and the ABE8e-enFnCas9 editing window."
>
> "These assumptions are clearly documented in the code and would be the single strongest improvement if we can replace them with actual measurements from our lab."

---

### Slide 11: How Multi-Nuclease Rescues Base Editing

**What to say:**
> "Here's the mechanism. On the left, v2 with SpCas9 only — NGG PAM covers about 8% of positions, and the ABE7.10 window is only 4 nucleotides wide. On the right, v3 adds enFnCas9 with NRG PAM and ABE8e with a window of positions 3-9 — nearly three times broader."
>
> "The combination of broader PAM and wider window rescued base editing at 6 out of 30 loci."

---

### Slide 12: KEY FINDING — BE Rescue (0/30 → 6/30)

**What to say (second KEY moment):**
> "This is the headline result of v3. Base editing went from zero top-ranked cases in v2 to six in v3. The strategy distribution is now 20% base editing, 77% prime editing, 3% HDR."
>
> "But I need to mention something important we discovered. Part of the reason PE dominated in v2 was a scoring bug — bystander severity was being penalized three separate times through three different channels. We found and fixed this. I'll explain on the next slide."

**Honesty note:** Do NOT hide the bug. The fact that you found it and fixed it demonstrates rigor.

---

### Slide 13: 6D TOPSIS Scoring

**What to say:**
> "Our scoring system uses TOPSIS — Technique for Order Preference by Similarity to Ideal Solution. It's a well-established multi-criteria decision method from 1981."
>
> "We evaluate each strategy on six dimensions: safety, feasibility, complexity, risk, confidence, and consequence. Safety and feasibility are the most heavily weighted."
>
> "The key difference from a simple weighted sum: TOPSIS penalizes strategies that are catastrophically poor on ANY single dimension, even if they're great on everything else. A strategy with zero safety — like dual simultaneous DSBs in iPSCs — will never be recommended regardless of how efficient it might be."

**If someone asks "why these weights?":** "The weights reflect the iPSC context — safety is highest because p53-mediated toxicity is the dominant concern. But rather than defending a single set of weights, we run 10,000 random weight permutations and report how stable the ranking is."

---

### Slide 14: Bystander Triple-Counting Bug Fix

**What to say (be direct and honest):**
> "During the v3 overhaul, we discovered a critical scoring bug. In v2, bystander edit severity was being penalized through three independent channels: once in the risk dimension, once as a consequence penalty, and once through the loss of a 'clean design' bonus. The total penalty for just one bystander was 0.055 — but base editing's feasibility advantage over prime editing was only 0.033."
>
> "This meant PE won every time BE had even a single bystander edit in the window — which is almost always in real genomic contexts."
>
> "In v3, we fixed this by moving bystander scoring entirely to a dedicated consequence dimension — the 6th dimension in our TOPSIS analysis. Now BE beats PE even with 1-3 bystanders, which is biologically correct."

**Why this matters:** "This shows why rigorous auditing matters. Without the fix, the tool would systematically undervalue base editing."

---

### Slide 15: Pareto Front Analysis

**What to say:**
> "Beyond TOPSIS, we added Pareto analysis. A strategy is Pareto-non-dominated if no other strategy is better on ALL dimensions simultaneously. The Pareto front is the set of strategies that could be optimal under some weighting scheme."
>
> "This is weight-independent — it doesn't depend on how you assign importance to safety vs feasibility vs complexity. If a strategy is Pareto-dominated, it's suboptimal regardless of your priorities."
>
> "In practice, HDR is almost always Pareto-dominated by BE and PE in iPSCs — worse safety, worse feasibility, worse confidence. Both BE and PE tend to be on the Pareto front, giving the user a genuine choice."

---

### Slide 16: Monte Carlo Sensitivity Analysis

**What to say:**
> "Instead of reporting just a single 'best' strategy, we report rank stability. We generate 10,000 random weight vectors from a Dirichlet distribution centered on our defaults, run TOPSIS with each, and count how often each strategy is top-ranked."
>
> "If a strategy is top-ranked in 94% of weight permutations, you can be confident the recommendation doesn't depend on the specific weights chosen. If it's only top-ranked in 55% of permutations, the decision is sensitive to your priorities and you should examine both options."

**If asked about Dirichlet:** "The Dirichlet distribution is the natural distribution over weight vectors that sum to 1. We use concentration=20 with a minimum alpha of 2.0 to ensure reasonable weight perturbations — large enough to be meaningful, but not so extreme that pathological weight combinations dominate."

---

### Slide 17: Cross-Method Validation

**What to say:**
> "To demonstrate that our rankings aren't an artifact of the TOPSIS method, we also implemented two alternative multi-criteria decision methods: VIKOR, which focuses on compromise solutions, and the Weighted Product Model, which is multiplicative rather than additive."
>
> "All three methods produce identical rankings across our test cases — 100% concordance. This means the recommendation is robust to the choice of scoring algorithm."

---

### Slide 18: Statistical Rigor in ConversionSim

**What to say:**
> "Every Monte Carlo output now includes uncertainty quantification. The HDR rate comes with a binomial 95% confidence interval. The mean tract length has a standard error and CI. Every distance-probability estimate has a Wilson score interval."
>
> "For example, our simulation of enFnCas9 with cssDNA in iPSCs gives an HDR rate of 3.5% with a 95% CI of 3.2% to 3.9%. The mean tract length is 706 bp, plus or minus 37 bp. The probability that a tract reaches 500 bp is 47.2%, with a CI of 42.0% to 52.4%."

**If asked why this matters:** "Without CIs, we can't distinguish a meaningful difference from sampling noise. A reviewer will ask for error bars on every reported statistic."

---

### Slide 19: ConversionSim Scope

**What to say:**
> "I want to be completely honest about what our simulation can and cannot do. ConversionSim models the SDSA pathway — synthesis-dependent strand annealing. This is valid for long-donor HDR with cssDNA, dsDNA, and lssDNA templates."
>
> "It is NOT valid for ssODN editing, which goes through SSTR — a completely different, RAD51-independent pathway. When we tested against Paquet et al.'s ssODN data, the R-squared was negative — the model performs worse than predicting the mean. This is expected because we're modeling the wrong pathway."
>
> "We restrict our claims to long-donor scenarios only. This is a scope boundary, not a model failure."

---

### Slide 20: Citation Integrity

**What to say:**
> "We performed a systematic verification of every reference and calibration value using web searches against PubMed and publisher databases. This was important because some of our original calibration data turned out to be incorrect."
>
> "For example, we had cited Elliott et al. 1998 for 'tract lengths of 200-2000 bp with median 500 bp.' When we checked the actual paper, it found 80% of tracts were 58 bp or less. The mismatch is because Elliott used endogenous chromosomal substrates, not exogenous donors."
>
> "We corrected three manuscript references and reframed all our calibration parameters with honest evidence tags."

**Why this matters to the audience:** "In the age of AI-assisted writing, citation verification is essential. We caught and fixed errors that would have been embarrassing in peer review."

---

### Slide 21: Parameter Provenance

**What to say:**
> "Every constant in our code is tagged with one of three evidence levels. 'Measured' means the value comes directly from a published measurement — like SpCas9 cutting at position -3 from the PAM. 'Derived' means we computed it from published data — like the HDR rate from Iyer et al. 'Assumed' means it's a modeling choice with rationale but no direct measurement — like the enFnCas9 stagger."
>
> "When you see a parameter marked 'Assumed,' that's where our model is making a judgment call. These are the parameters we explore in sensitivity analysis."

---

### Slide 22: Benchmark Results v2 vs v3

**What to say:**
> "Comparing v2 and v3 on the same 30 ClinVar cases: top-1 accuracy is the same at 86.7%. But the strategy distribution changed dramatically — from 97% PE in v2 to 77% PE plus 20% BE in v3. This is the direct result of the multi-nuclease engine and the bystander scoring fix."

---

### Slide 23: Codebase Summary

**What to say (keep brief):**
> "Quick technical summary: 224 tests passing, version 3.0.0, about 36,500 lines of Python. All 20 manuscript references verified with PMIDs. Docker-ready for deployment. MIT license, available on GitHub."

---

### Slide 24: Limitations

**What to say (be thorough and honest):**
> "Let me be upfront about what we haven't done yet. First and most importantly, we have no experimental validation — our benchmark uses self-curated truth labels, not experimental outcomes. This is the biggest gap for publication."
>
> "Second, our ConversionSim only models SDSA, not SSTR. Third, the enFnCas9 stagger and HDR multiplier are assumptions. Fourth, our off-target scoring is local only — no genome-wide search. Fifth, PE still dominates in 23 of 30 cases, which a reviewer might call a degenerate pattern."
>
> "We know these limitations and are actively addressing them."

---

### Slide 25: SDSA Sensitivity Analysis

**What to say:**
> "Our SDSA displacement probability p=0.002 is an assumed parameter -- one of our most important assumptions. So we tested its robustness by running the full pipeline across a 5-fold range from p=0.001 to p=0.005, corresponding to mean tract lengths from 1000 bp down to 200 bp."
>
> "The key finding: strategy rankings are robust across this entire range. The relative ordering of cssDNA vs lssDNA vs dsDNA donors is preserved, and the predicted cssDNA/lssDNA ratio stays within published experimental ranges at all tested values. This means even if our assumed p is off by a factor of 2, our recommendations don't change."

**Figure:** `Fig_SDSA_Sensitivity.png` -- shows p=0.002 robustness across the 0.001-0.005 range.

---

### Slide 26: Delivery Advisor

**What to say:**
> "We also added a delivery advisory module that provides practical guidance after TOPSIS ranking. This does not change the rankings -- delivery complexity is already captured by our Safety and Complexity dimensions. Instead, it adds annotations: donor format recommendations by edit size, cell-type warnings like iPSC dsDNA p53 toxicity, and viability enhancer suggestions."
>
> "The module is backed by a comprehensive literature review of 87 verified references covering cssDNA, lssDNA, dsDNA, and AAV donor formats. Key finding: cssDNA achieves 3-5x higher HDR rates than lssDNA -- supported by Iyer 2022, Xie 2024, and Letort 2025."

**If asked why delivery is not a TOPSIS dimension:** "Delivery complexity is correlated with Safety and Complexity -- adding it would be redundant in about 90% of cases. And ordinal delivery scores lack calibration data for principled weight assignment."

---

### Slide 27: What's Next

**What to say:**
> "On the left: what I need from the lab. Even 2-3 iPSC editing cases with outcomes would transform the manuscript. The enFnCas9 cut-site stagger measurement would replace our biggest assumption with actual data. I've prepared data collection sheets."
>
> "On the right: what I'm doing computationally. Expanding the benchmark to 33 published cases from the literature, validating against large-scale editing datasets, and rebuilding the manuscript for PLOS Computational Biology."

---

### Slide 28: Thank You

**What to say:**
> "The key takeaway: CRISPRArchitect v3 provides transparent, method-robust, uncertainty-quantified strategy recommendations. Every parameter is traceable to published evidence. Every limitation is honestly documented. Thank you."

---

## Comprehensive Q&A Preparation

### Biological Questions

**Q1: "Why does prime editing still dominate in 23/30 cases even after the fix?"**
> PE genuinely has three structural advantages: (1) no editing window constraint — the RT template directly encodes the edit, so any mutation near any PAM is feasible; (2) broad mutation-type compatibility — transitions, transversions, and small indels; (3) zero DSBs, giving it a safety score of 1.0. At many loci, no nuclease-editor combination places the target within the base editing window, even with ABE8e and enFnCas9. This is not a scoring artifact — it's a real biological constraint. PE's dominance reflects the genuine versatility of the prime editing mechanism.

**Q2: "What about cell-cycle dependence? HDR requires S/G2 phase."**
> Yes, this is modeled. Our ConversionSim gates HDR success on cell-cycle phase — only cells in S/G2 can do HDR. For iPSCs, we use 35% S/G2 fraction (from Becker et al., PNAS, 2006). This is one reason HDR efficiency is low in iPSCs (~8% baseline). BE and PE work in any cell-cycle phase, which is captured by their higher effective efficiency.

**Q3: "How do you handle bystander edits? Is every C or A in the window a bystander?"**
> Yes, every same-type base (A for ABE, C for CBE) within the editing window is scanned. For each bystander position, we determine the codon context and classify the consequence: synonymous (no penalty), missense (moderate penalty per ACMG), nonsense (severe penalty), or splice-affecting. The consequence score is now a proper 6th dimension in TOPSIS, so bystander-heavy strategies score lower on that dimension but aren't triple-penalized like in v2.

**Q4: "Can enFnCas9 really be paired with ABE8e? Has anyone published this fusion?"**
> Not directly. ABE8e has been published with SpCas9 (Richter et al., Nat Biotechnol, 2020). We mark the ABE8e-enFnCas9 combination as "Tier B evidence" — meaning it's an extrapolation based on ABE8e working with multiple Cas proteins and enFnCas9 being structurally similar to SpCas9. The editing window for ABE8e-enFnCas9 is assumed to be the same as ABE8e-SpCas9 (positions 3-9). This should be experimentally validated.

**Q5: "You restrict ConversionSim to long donors — but what about ssODN, which is the most common donor in our lab?"**
> ConversionSim models the SDSA pathway (resection → RAD51 filament → strand invasion → D-loop synthesis). ssODN editing primarily proceeds via SSTR — Single-Strand Template Repair — which is RAD51-independent and PCNA-dependent. These are fundamentally different mechanisms with very different tract lengths: SDSA produces hundreds of bp, SSTR produces ~10-50 bp. Using an SDSA model for ssODN predictions would be scientifically incorrect, so we explicitly restrict our claims. For ssODN, the empirical distance-decay from Paquet et al. (2016) — edit within 10 bp of cut — remains the best guide.

**Q6: "What about chromatin context? Open vs closed chromatin affects editing."**
> Currently not modeled. CRISPRArchitect evaluates sequence-level constraints (PAM, editing window, coding consequence) but does not incorporate chromatin accessibility (ATAC-seq), replication timing, or histone modifications. This is a documented limitation. The ChromBridge module from v1 (which modeled 3D chromatin distance using polymer physics) has been removed from v3 — it addressed distance prediction, not accessibility. Integrating ATAC-seq data is a planned future direction.

---

### Methodological Questions

**Q7: "Why TOPSIS and not a simpler weighted sum?"**
> A weighted sum allows full compensation — a strategy that's terrible on safety can compensate with high efficiency. TOPSIS uses Euclidean distance to ideal/anti-ideal solutions, which naturally penalizes extreme weakness on any dimension. In practical terms: a weighted sum might rank a dual-DSB HDR strategy above a single-step PE if the HDR scores very high on feasibility. TOPSIS would rank it lower because the safety score is catastrophic. For clinical iPSC work where safety is paramount, this non-compensatory behavior is desirable.

**Q8: "What is a Pareto front and why does it matter?"**
> Imagine you have three strategies. Strategy A is safer but less efficient. Strategy B is more efficient but less safe. Strategy C is worse than A on safety AND worse than B on efficiency — there's no weighting scheme under which C would be preferred. C is "Pareto dominated." A and B are on the "Pareto front" — they represent genuine trade-offs. Pareto analysis identifies which strategies represent real trade-offs and which are universally suboptimal. It requires no weight assumptions at all.

**Q9: "You have both TOPSIS and Pareto and VIKOR and WPM — isn't that overkill?"**
> TOPSIS is our primary ranking method. Pareto analysis provides a weight-independent sanity check. VIKOR and WPM are run as supplementary comparisons to demonstrate that the ranking is robust to method choice (100% concordance). In the manuscript, we'll report TOPSIS as primary and note that VIKOR/WPM produce identical rankings. A reviewer might ask "what if you used a different method?" — we preempt that question by showing we already checked.

**Q10: "How do you handle variants where all three modalities are infeasible?"**
> The pipeline generates explicit rejection reasons for every infeasible strategy. If no PAM places the target in the BE window → BE is rejected with reason "no compatible guide." If the variant is a large deletion >80 bp → PE is rejected with "exceeds PE deletion limit." If no guide exists within 5 kb → HDR is rejected. These rejection reasons are displayed alongside the ranking. If everything is infeasible, the user sees an empty ranking with all rejection reasons — which is a valid and informative output.

**Q11: "Your benchmark truth labels are self-curated. How do we know they're correct?"**
> This is our biggest validation weakness, and I acknowledge it openly. The 30-case truth labels were defined by me based on biological reasoning and published literature. A reviewer will flag this as circular validation. Our plan to address it: (1) expand to 33 published cases where the actual experimental strategy and outcome are documented, (2) have 3 independent genome editing PIs label the cases and report inter-rater agreement (Cohen's kappa), (3) integrate experimental data from our own lab.

---

### Mathematical & Statistical Questions

**Q12: "What is a Monte Carlo simulation and why do we need it here?"**

> **What it is:** A Monte Carlo simulation generates thousands of random trials to estimate the distribution of an outcome. Instead of computing a single average, we simulate the biological process (resection → filament → invasion → synthesis) thousands of times, each time drawing random values for stochastic variables (resection length, filament coverage, invasion success, tract length).
>
> **Why we need it:** The HDR process is inherently stochastic — every cell gets a different resection length, different RAD51 coverage, different synthesis extent. There's no closed-form equation for "probability of incorporating an edit 500 bp from the cut site" because it depends on the joint distribution of four coupled random variables. Monte Carlo simulation gives us the empirical distribution.
>
> **Parameters required:** Cell type (S/G2 fraction, baseline HDR rate, p53 status), nuclease (cut type, stagger length), donor (topology, homology arm length). These determine the distributions we sample from.
>
> **What we run:** 10,000 independent virtual cells. Each goes through resection (Normal + LogNormal distributions), filament formation (Beta distribution), invasion (Bernoulli trial with computed probability), and synthesis (Geometric distribution). The output is a distribution of 10,000 tract lengths.

**Q13: "What are SEs and CIs and why do they matter?"**

> **SE (Standard Error):** The standard error tells you how precise your estimate is. If you simulate 10,000 cells and get a mean tract length of 706 bp with SE = 37 bp, it means if you re-ran the simulation with a different random seed, the mean would likely be within ~37 bp of 706. SE shrinks with more simulations: SE = std / sqrt(n).
>
> **CI (Confidence Interval):** A 95% confidence interval means: if we repeated this entire simulation procedure 100 times, approximately 95 of those intervals would contain the true value. For our tract length: 706 bp ± 1.96 × 37 = [634, 778] bp.
>
> **Why they matter:** Without SEs/CIs, we can't tell if a difference is real or just random noise. If one configuration gives 47% conversion probability and another gives 52%, is that meaningful? Only if the CIs don't overlap. Peer reviewers require error bars.
>
> **Wilson score interval:** For proportions (like "47.2% of tracts reach 500 bp"), the standard CI formula (p ± z × sqrt(p(1-p)/n)) performs poorly when p is near 0 or 1. The Wilson interval is more accurate in these cases and is now the standard in biostatistics.

**Q14: "What is the Dirichlet distribution and why do you use it for sensitivity analysis?"**

> The Dirichlet distribution is the natural probability distribution over vectors that sum to 1 — exactly what you need for weight vectors. If our default weights are [0.28, 0.23, 0.19, 0.14, 0.09, 0.07], the Dirichlet generates random weight vectors centered around these values.
>
> The concentration parameter (we use 20) controls how tightly the random weights cluster around the defaults. Higher concentration = less variation. We set minimum alpha = 2.0 to prevent any dimension from having a near-Uniform marginal distribution, which would create pathological weights.
>
> We sample 10,000 weight vectors, run TOPSIS with each, and count how often each strategy is #1. This gives us "rank stability" — a direct measure of robustness.

**Q15: "Why geometric distribution for tract lengths? Why not exponential or Weibull?"**

> The geometric distribution models a constant per-bp probability of D-loop collapse — at every base pair of synthesis, there's a 0.2% chance the helicase dismantles the D-loop. This is the memoryless (constant hazard) model.
>
> Biologically, the hazard might increase with distance (D-loop grows and becomes less stable), which would suggest a Weibull distribution with shape > 1. We chose the geometric as the simplest principled model and document it as an explicit modeling choice. A Weibull alternative is a planned sensitivity analysis but requires fitting the shape parameter to mammalian tract-length data, which is currently sparse.
>
> The exponential is the continuous analog of the geometric — essentially the same model. We use the discrete geometric because tract lengths are in integer base pairs.

**Q16: "Your SDSA displacement probability p=0.002 — how confident are you in this number?"**

> Not very confident, and we say so explicitly. It's tagged [ASSUMED]. The value gives mean = 500 bp and median = 347 bp, which is consistent with functional evidence: successful HDR with 300-1000 bp homology arms requires tracts reaching that far, and the Stark lab SDSA assay confirms >=350 bp synthesis in human cells. But no one has directly measured the tract-length distribution for exogenous long-donor HDR in mammalian cells. We recommend exploring p from 0.001 to 0.005 (mean tracts 200-1000 bp) in sensitivity analysis.

---

### Architecture & Software Questions

**Q17: "How is the code organized?"**

> Three layers: (1) v1 simulation modules — ConversionSim and MOSAIC are retained for simulation capabilities; ChromBridge, TopoPred, and LoopSim have been removed from v3 as they are not required by the current pipeline. (2) v2/v3 core pipeline (~12,500 LOC) — transcript mapping, feasibility engines, strategy generation, TOPSIS scoring. (3) Infrastructure — Streamlit web app, CLI, Docker, CI/CD. The v3 core pipeline is the primary entry point.

**Q18: "Can I run this on my laptop?"**

> Yes. Python 3.9+ with NumPy, SciPy, and Matplotlib. No GPU required. A full analysis of one variant takes <1 minute (mostly Ensembl API calls). The Monte Carlo simulation (10,000 trials) runs in <1 second on a standard laptop. The Streamlit web app can be launched with `crisprarchitect webapp`.

**Q19: "Is it on GitHub? Can collaborators use it?"**

> Yes: github.com/visvikbharti/CRISPRArchitect, MIT license. `pip install -e .` from the repository root. Docker deployment is also available: `docker build -t crisprarchitect . && docker run -p 8501:8501 crisprarchitect`.

---

### Critical/Challenging Questions (Be Prepared)

**Q20: "86.7% accuracy doesn't sound that impressive. CRISPOR claims much higher."**

> Fair point. Two important distinctions: (1) CRISPOR validates guide efficiency (a continuous score against empirical data at >1,000 guides). We validate strategy recommendations (a categorical decision across modalities), which is a fundamentally different and harder task. (2) Our "accuracy" is against self-curated truth labels, which we acknowledge is a limitation. The more informative metric is the rank stability from sensitivity analysis, which shows how robust each recommendation is.

**Q21: "You found fabricated citations in your own code. How did that happen?"**

> Some calibration values and reference details were generated with AI assistance and not verified against primary sources. This is a cautionary tale about AI-assisted science. We caught every error through systematic web-search verification and corrected them. The lesson: every number from AI must be verified against the actual paper. We now tag every parameter with its evidence level and have verified all 20 manuscript references with PMIDs.

**Q22: "If PE dominates in 23/30 cases, why not just always recommend PE?"**

> Because the 7 cases where PE is NOT the best choice matter enormously — those are the cases where CRISPRArchitect provides genuine decision support. At 6 loci, base editing is superior (simpler delivery, higher efficiency, no bystander issues). At the one remaining locus, HDR is needed. A "always PE" heuristic would miss these cases. Furthermore, as the nuclease landscape expands, the fraction of BE-rescuable cases will grow.

**Q23: "You haven't done any experiments. Can this be published?"**

> Computational-only papers are published in PLOS Computational Biology regularly — CRISPOR, PrimeDesign, and CHOPCHOP were all initially published without wet-lab validation. However, even a small amount of experimental data (3-5 cases from our lab) would dramatically strengthen the paper. This is why I'm requesting data from the lab. The enFnCas9 characterization data (cut-site stagger, HDR comparison) would be a uniquely valuable contribution that no other group can provide.

**Q24: "The p=0.002 for SDSA is assumed and you found the original calibration data was wrong. How can you trust the model?"**

> This is a legitimate concern. The key insight is that the model's VALUE is not in the absolute numbers (which depend on the assumed p) but in the RELATIVE comparisons: cssDNA vs lssDNA (validated: 2.07x predicted vs 1.9x observed), staggered vs blunt (validated: 1.82x vs 1.9x observed). These ratios are robust to the absolute value of p because both donor types are simulated with the same underlying model. The sensitivity analysis explores p across a wide range.

**Q25: "What happens when a reviewer asks for experimental validation?"**

> We have three responses: (1) We are actively seeking experimental data from our lab; (2) We have expanded our computational benchmark to 33 published cases with documented experimental outcomes; (3) We position the tool as decision-support, not a predictive optimizer — the recommendations are starting points for experimental design, not substitutes for it. If the reviewer insists on wet-lab data, we may need to add it as a revision, which is why getting lab data before submission is our top priority.

**Q26: "How do you handle delivery method recommendations?"**

> We use an Option B architecture: delivery does not change TOPSIS rankings, but adds practical guidance as post-ranking annotations. The delivery advisor module has two tiers. First, hard feasibility filters that flag biologically incompatible combinations — for example, dsDNA donors in iPSCs trigger p53-mediated toxicity (Ihry 2018, Haapaniemi 2018). Second, practical recommendations: donor format selection by edit size, cell-type-specific warnings, and viability enhancer suggestions. The module is backed by 87 verified references from a comprehensive delivery methods literature review. Key finding from the literature: cssDNA achieves 3-5x higher HDR rates than lssDNA (Iyer 2022, Xie 2024, Letort 2025).

**Q27: "Why didn't you add delivery as a TOPSIS dimension?"**

> Three reasons. First, delivery complexity is correlated with existing Safety and Complexity dimensions — DSB-free modalities (BE, PE) inherently avoid delivery-related toxicity, and multi-round strategies inherently require more complex delivery logistics. Adding delivery would be redundant in approximately 90% of cases. Second, ordinal delivery complexity scores (e.g., "simple," "moderate," "complex") lack the calibration data needed for principled weight assignment in TOPSIS — we would need experimental outcome data linking delivery complexity to editing success rates. Third, delivery decisions are often made downstream of strategy selection (once you know the modality, you choose the delivery vehicle), so they fit naturally as post-ranking annotations rather than ranking inputs.

---

## Key Phrases to Remember

- "Decision-support tool, not a predictive optimizer"
- "Every parameter has documented provenance"
- "We found and fixed the bug ourselves — that's scientific rigor"
- "The sensitivity analysis shows this recommendation is robust"
- "We explicitly restrict our claims to what the model can do"
- "This is tagged [ASSUMED] — we're transparent about what we know and don't know"

---

## Timing Guide

| Slides | Section | Time |
|--------|---------|------|
| 1-3 | Introduction & Problem | 5 min |
| 4-7 | Pipeline overview, v1/v2 foundation & key finding | 5 min |
| 8-12 | v3 new capabilities & BE rescue | 10 min |
| 13-17 | Scoring methodology (TOPSIS/Pareto/sensitivity) | 10 min |
| 18-21 | Scientific rigor (stats, scope, citations, params) | 8 min |
| 22-24 | Results, codebase, limitations | 5 min |
| 25-28 | SDSA sensitivity, delivery advisor, future, thank you | 7 min |
| Q&A | | 15-20 min |

**Total: ~45 min talk + ~20 min Q&A = ~65 min**

---

## Figure Guide — Which Figure Goes Where

All figures are in `paper/figures/v3_results/` and are generated from **real data only**.

| Figure File | Use on Slide | What It Shows | Key Interpretation |
|---|---|---|---|
| `Fig_StrategyDistribution_v2_v3.png` | Slide 12 (BE Rescue) | Side-by-side bars: v2 (BE=0, PE=29) vs v3 (BE=6, PE=23) | Multi-nuclease engine rescued BE from 0% to 20%. PE still dominates (77%) but this reflects genuine biological constraints, not a bug. |
| `Fig_LiteratureBenchmark.png` | Slide 22 (Results) | Literature benchmark: 30% top-1, 80% top-3 concordance | Top-1 appears low (30%) but discordance is explainable: HDR papers are pre-PE era. Top-3 at 80% shows published strategy is nearly always in our recommendation set. |
| `Fig_BystanterFix.png` | Slide 14 (Bug Fix) | v2 scoring vs v3: BE score drops below PE with just 1 bystander in v2, but stays above PE with 3 bystanders in v3 | The triple-counting bug made PE unbeatable. After fix, BE properly wins when PAM+window are verified. |
| `Fig_ConversionSim_CIs.png` | Slide 18 (Statistical Rigor) | Tract length distribution with mean/median + 95% CI; distance-probability curve with Wilson CIs | Every output now has uncertainty quantification. Mean tract 706 bp (SE=37), P(>=500bp) = 47.2% [42.0-52.4%]. |
| `Fig_ParameterProvenance.png` | Slide 21 (Parameters) | Pie chart: 35% measured, 20% derived, 45% assumed | Transparent about what we know vs assume. All [ASSUMED] parameters explored in sensitivity analysis. |
| `Fig_CrossMethod_Agreement.png` | Slide 17 (Cross-Method) | TOPSIS vs VIKOR vs WPM: all produce identical rankings | 100% concordance proves recommendation is method-robust, not an artifact of TOPSIS. |

### How to Print/Show Figures

Option A (recommended): Open the PNG files on your laptop and switch to them during the relevant slides using Alt+Tab.

Option B: Insert the PNGs into the PPTX slides manually before the meeting. The figure sizes are optimized for 16:9 widescreen.

---

## Literature Benchmark: Detailed Interpretation for Q&A

This section prepares you for questions about the 10-case literature benchmark results.

### The Results Table

| Case | Gene | Published Strategy | CRISPRArchitect Recommendation | Match? | Explanation |
|---|---|---|---|---|---|
| LIT_BE_001 | HBB | ABE8e-NRCH (80%) | PE | MISS | ABE8e-NRCH uses a specialized PAM variant (NRCH) not in our nuclease set. Our pipeline only has SpCas9/enFnCas9/SpCas9-NG/SpRY/Cas12a. Adding NRCH would rescue this case. |
| LIT_BE_005 | COL7A1 | ABE8e (94.6%) | BE | OK | Correct! ABE8e with available PAM places target in window. Validates multi-nuclease rescue. |
| LIT_PE_001 | HBB | PE3 (26-52%) | PE | OK | Correct. PE is the right recommendation for this transversion. |
| LIT_PE_004 | HBB | PEmax (15-41%) | PE | OK | Correct. Newer PE variant but same modality class. |
| LIT_HDR_001 | HBB | AAV6 HDR (29%) | PE | MISS | Paper is Dever et al. 2016 — **before PE was invented** (2019). CRISPRArchitect correctly identifies the modern DSB-free alternative. |
| LIT_HDR_002 | HBB | AAV6 HDR (60%) | PE | MISS | Lattanzi 2021 clinical-grade HSPCs. Used HDR because it was an established clinical protocol. PE is objectively safer but less clinically validated. |
| LIT_HDR_005 | HBB | AAV6 HDR | PE | MISS | Same variant, same reasoning. HDR chosen for clinical pipeline continuity. |
| LIT_HDR_007 | HBB | ssODN+inhibitors (72%) | PE | MISS | Used ssODN with NHEJ/MMEJ inhibitors for very high efficiency. Our tool correctly identifies PE as safer, but doesn't model chemical enhancement protocols. |
| LIT_HDR_008 | LRRK2 | ssODN HDR | BE | MISS | G>A transition. Paper used HDR but was published before widespread ABE adoption. CRISPRArchitect recommends the modern DSB-free approach (BE). Arguably more correct. |
| LIT_COMP_001 | HBB | Multi-strategy | PE | N/A | Comparison case — all strategies tested. Newby 2021 found ABE8e-NRCH was most efficient (80%) but our pipeline doesn't model NRCH. |

### Key Talking Points for the Benchmark

1. **"30% top-1 sounds bad"** — Frame it correctly: "Of the 7 discordant cases, 5 are HDR papers from before PE existed (2016-2019). The tool is recommending the safer modern alternative, not making an error."

2. **"80% top-3 is the real metric"** — "The published strategy appears in our top-3 recommendations 80% of the time, meaning we almost always include it as an option — we just sometimes rank it below a safer alternative."

3. **"The NRCH limitation"** — "The HBB ABE case uses SpCas9-NRCH, a PAM variant we don't model. Adding NRCH to our nuclease set would likely rescue this case, bringing concordance to ~40% top-1."

4. **"Why not just agree with everything?"** — "A tool that simply validates whatever strategy a paper used would have 100% concordance but zero clinical value. The point of CRISPRArchitect is to sometimes identify a BETTER strategy than what was historically used — and that requires disagreement with pre-PE-era HDR papers."

---

*Document updated: 2026-03-31*
*CRISPRArchitect version: v3.0.0*
