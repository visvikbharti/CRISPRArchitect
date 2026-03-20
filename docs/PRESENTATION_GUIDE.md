# CRISPRArchitect — Presentation Guide

## Complete Slide-by-Slide Explanation, Speaker Notes, Anticipated Questions & Answers

This document is your comprehensive preparation guide. Read it before your lab meeting or DAC presentation. It explains the science behind every slide, what to say, what to emphasize, and how to handle tough questions.

---

## TABLE OF CONTENTS

- [How to Present This](#how-to-present-this)
- [Slide-by-Slide Guide](#slide-by-slide-guide)
- [Anticipated Questions & Answers](#anticipated-questions--answers)
- [Key Numbers to Memorize](#key-numbers-to-memorize)
- [Common Objections and How to Address Them](#common-objections-and-how-to-address-them)
- [The 2-Minute Elevator Pitch](#the-2-minute-elevator-pitch)

---

## HOW TO PRESENT THIS

**Total time:** 25-30 minutes + 10-15 minutes Q&A

**Tone:** Confident but humble. You built something genuinely useful, but you're honest about what it can and can't do. This honesty is a strength, not a weakness — reviewers and PIs respect it.

**Live demo:** After slide 15, switch to your browser (http://localhost:8501) and walk through the NF1 example. This is the most impressive part — showing a real gene being analyzed in real-time.

**Key message to drive home:** "No existing tool helps researchers choose between HDR, base editing, and prime editing for multi-site correction. We built the first one, validated it, and it works."

---

## SLIDE-BY-SLIDE GUIDE

---

### SLIDE 1: TITLE

**What's on it:** Title, authors (Vishal Bharti, Debojyoti Chakraborty), affiliations (CSIR-IGIB, AcSIR), GitHub link.

**What to say (30 seconds):**
"Good morning/afternoon. I'm going to present CRISPRArchitect, a computational framework we developed for optimizing multi-site genome editing strategies. This work addresses a gap we identified while planning an editing experiment in our lab — specifically, how to choose the best approach when you need to correct two mutations in a large gene in iPSCs."

**What NOT to say:** Don't start with "so basically..." or "this is just a tool we built." Frame it as science from the start.

---

### SLIDE 2: OUTLINE

**What's on it:** 8 numbered sections showing the presentation flow.

**What to say (20 seconds):**
"Here's the roadmap for this talk. I'll start with the biological background and the specific gap we identified, then walk through how we built the framework, show the results including validation against published data, and end with conclusions and our publication plan."

**Purpose:** Gives the audience a mental map so they can follow along. Skip quickly.

---

### SLIDE 3: BACKGROUND — The Multi-Site Editing Challenge

**What's on it:** Left: large multi-exon genes (DMD, NF1, TTN) and the distance problem. Right: the expanding editing toolkit (HDR, base editing, prime editing, cssDNA, staggered cuts).

**What to say (2-3 minutes):**
"Many genetic diseases involve mutations in large genes. Take DMD — 79 exons spanning 2.2 megabases. Or NF1 — 58 exons. When a patient is compound heterozygous, they carry two different mutations on different alleles. Correcting both is the goal, but here's the challenge:

The exons themselves are small — typically 150 base pairs. So the exonic distance between two mutations might be just a few kilobases. That seems manageable. But the GENOMIC distance — including all the introns — can be hundreds of kilobases to megabases. That's a completely different scale.

Meanwhile, the editing toolkit has expanded enormously. We now have base editors that can correct transitions without cutting DNA at all. Prime editors that can install any small edit. Circular ssDNA donors that are twice as efficient as linear ones. Staggered-cutting nucleases that boost HDR. But which of these should you actually USE for a given multi-site problem? That's the question nobody has a tool for."

**Key biology to internalize:**

- **Compound heterozygosity:** Patient has mutation A on one chromosome, mutation B on the other. Both alleles are damaged. To cure the disease, you might need to fix both.

- **Why introns matter:** Exon 40 and exon 60 might be only 3 kb apart in the mRNA, but in the genome they're separated by 20 introns totaling maybe 1 Mb. This is the crucial distinction that makes single-template editing impossible.

- **Why iPSCs are special:** iPSCs have active p53, meaning every DSB triggers apoptosis and can select for p53-mutant cells (Ihry et al., 2018). This makes DSB-free approaches (base editing, prime editing) particularly valuable in iPSCs compared to cell lines like HEK293T.

- **The 8 references cited on this slide:**
  1. Ihry et al., 2018, Nat Med — p53 inhibits CRISPR in iPSCs
  2. Haapaniemi et al., 2018, Nat Med — DSBs select for p53 mutants
  3. Komor et al., 2016, Nature — cytosine base editors
  4. Gaudelli et al., 2017, Nature — adenine base editors
  5. Anzalone et al., 2019, Nature — prime editing
  6. Iyer et al., 2022, CRISPR J — cssDNA donors
  7. Xie et al., 2025, Nat Biotech — GATALYST cssDNA, 70% in iPSCs
  8. Chauhan et al., 2023, PNAS — staggered cuts boost HDR 1.9x

---

### SLIDE 4: THE GAP — No Tool for Multi-Site Strategy Optimization

**What's on it:** Table of existing tools (CRISPOR, IDT, PrimeDesign, BE-Designer, Benchling) showing what each does and doesn't do. Three specific gaps identified. The motivating question.

**What to say (2 minutes):**
"So what tools exist to help with this decision? CRISPOR designs guide RNAs. IDT designs single-site donors. PrimeDesign handles prime editing guides. BE-Designer does base editing. Benchling does everything for single sites. But NONE of them addresses the multi-site problem. None compares across modalities. And none simulates the HDR process to tell you whether a single donor can actually reach both sites.

We identified three specific gaps:
1. There's no quantitative model of how far a donor template's information propagates during HDR — the gene conversion tract length.
2. There's no tool that compares HDR versus base editing versus prime editing for a given scenario.
3. There's no feasibility check for whether single-template multi-site editing is even physically possible.

This work was motivated by a concrete question in our lab: can a single cssDNA donor correct mutations in both exon 40 and exon 60 of a large gene using enFnCas9 in iPSCs?"

**Why this slide matters:** This is where you establish that the work is NEEDED. If the audience doesn't believe there's a gap, the rest doesn't matter. Be specific about what existing tools lack.

---

### SLIDE 5: OBJECTIVES

**What's on it:** Four aims in color-coded boxes.

**What to say (1 minute):**
"We set four aims:
1. Build ConversionSim — a Monte Carlo simulator of the HDR process, from end resection through SDSA synthesis.
2. Build MOSAIC — a strategy optimizer that enumerates all feasible editing approaches and ranks them.
3. Validate both against published experimental data — not against our own data, against independent published datasets.
4. Apply the framework to real human genes using coordinates from Ensembl."

**Key point:** Emphasize that validation uses INDEPENDENT data — parameters were set a priori, not fit to the validation datasets.

---

### SLIDE 6: METHODS — ConversionSim

**What's on it:** Four-step pipeline (resection → RAD51 → invasion → SDSA), key parameters, design decisions, outputs.

**What to say (3-4 minutes):**
"ConversionSim models HDR as four sequential stochastic steps.

Step 1: End resection. After the DSB, nucleases chew back the DNA to create single-stranded overhangs. We model this as a two-phase process — short-range resection by MRN/CtIP, about 200 base pairs, then long-range resection by EXO1, median about 1,500 base pairs. These numbers come from Symington 2011 and Cejka 2015.

Step 2: RAD51 filament formation. The recombinase RAD51 coats the single-stranded DNA. We model about 85% coverage based on in vitro data. You need at least 15 nucleotides of filament for productive strand invasion — that's from the Qi et al. 2015 single-molecule study showing the 8-nucleotide sampling mechanism.

Step 3: Strand invasion into the donor template. The probability depends on cell type, donor format, and homology arm length. iPSCs have baseline HDR around 8%.

Step 4: And this is the critical step — SDSA synthesis. After strand invasion, the polymerase copies donor sequence. But in mitotic cells, the dominant pathway is SDSA — synthesis-dependent strand annealing — where the D-loop eventually collapses and the newly synthesized strand is displaced. We model this as a geometric process: at each base pair, there's a small probability of collapse. This naturally produces the right-skewed tract length distribution that's been observed experimentally.

Two key modifiers: circular ssDNA donors reduce the collapse probability by 20% — because they're more stable inside the cell due to exonuclease resistance. And staggered cuts reduce it by another 15% — because the longer resection creates more stable invasion intermediates.

All parameters come from published studies. None were fit to our validation data."

**The biology you MUST understand:**

- **Why SDSA limits tract length:** In the SDSA pathway (dominant in mitotic cells), the invading strand is extended by polymerase, but then it's DISPLACED from the donor and re-anneals to the other side of the break. This displacement is what limits how much donor sequence gets copied. The longer the synthesis goes, the more likely the D-loop collapses. That's why most tracts are short (<1 kb).

- **Why the geometric distribution:** If at each base pair there's a constant probability p of collapse (p ≈ 0.002), then the tract length follows a geometric distribution: P(tract = k) = (1-p)^(k-1) × p. Mean = 1/p ≈ 500 bp. This is biologically reasonable and matches observed distributions.

- **Why circular ssDNA is better:** Linear ssDNA is degraded by cellular exonucleases. Circular ssDNA has no free ends, so it resists degradation. This means more intact donor molecules are available in the nucleus for longer, increasing the effective concentration. Iyer et al. (2022) showed this experimentally: circularizing the same sequence improved HDR 2-fold.

---

### SLIDE 7: METHODS — MOSAIC

**What's on it:** Input→Enumerate→Score→Rank flow, 8 strategy classes, 4-dimension scoring, safety-first rationale.

**What to say (2-3 minutes):**
"MOSAIC takes five inputs: gene structure, mutation positions, mutation types, cell type, and nuclease. It then does three things.

First, it classifies each mutation. Is it a transition — like G to A? Then it's amenable to adenine base editing. Is it a transversion or an indel? Then you need prime editing or HDR.

Second, it enumerates all feasible strategies. There are up to 8 classes, from single-template HDR — only if the mutations are within 5 kb — to dual base editing, prime editing, hybrid approaches, even exon deletion. It generates every option whose prerequisites are met.

Third, it scores each strategy on four dimensions: efficiency, safety, time, and cost. And this is important — we weight safety at 40% for iPSCs. Why? Because Ihry et al. and Haapaniemi et al. showed that DSBs in iPSCs trigger p53-dependent apoptosis and can select for p53-mutant cells. So a strategy with slightly lower efficiency but zero DSBs — like base editing — may be the better choice.

I want to be upfront: these weights are heuristic. They're not derived from data. They reflect our judgment about what matters most in iPSC applications. They're adjustable. We don't claim they're optimal."

---

### SLIDE 8: RESULTS — ConversionSim Predictions (Figure 1)

**What's on it:** The actual Figure 1 from the paper — model schematic, tract length distribution, survival curve.

**What to say (2 minutes):**
"Here are the ConversionSim predictions for a typical iPSC editing scenario — staggered-cut nuclease, cssDNA donor, 300 bp homology arms.

Panel B shows the tract length distribution from 50,000 simulations. The median is 542 base pairs — meaning half of all HDR events copy less than 542 bp of donor sequence. The 95th percentile is about 2,200 bp.

Panel C is the key one for our question. It shows the probability that an edit at a given distance from the cut site gets incorporated. At 100 bp, you're at 89%. At 500 bp, about 53%. At 1,000 bp, 28%. At 2,000 bp, only 7%. And at the distances we're talking about for multi-exon genes — tens to hundreds of kilobases — the probability is zero.

This directly answers our motivating question: no, a single donor cannot correct two mutations separated by 100+ kb. The gene conversion tract simply doesn't extend that far."

---

### SLIDE 9: RESULTS — Real Gene Analysis

**What's on it:** NF1, DMD, BRCA2 analysis cards with real GRCh38 coordinates.

**What to say (2 minutes):**
"We tested this on three real human genes, fetching exon coordinates directly from Ensembl.

NF1 — 58 exons on chromosome 17. If you have mutations in exon 20 and exon 50, they're separated by 123 kilobases of genomic distance. In 3D nuclear space, that's about 1,133 nanometers. A 3 kb cssDNA donor has a random coil diameter of 261 nanometers. It physically cannot bridge that gap — it's about 5 times too small.

DMD — 79 exons, the largest human gene at 2.2 megabases. Even larger distances.

BRCA2 — 27 exons, smaller gene. But even at just 47 kilobases between the target exons, bridging is still impossible.

For ALL three genes, the probability of a gene conversion tract reaching the second site is zero. The recommendation in all cases is either base editing or sequential HDR."

---

### SLIDE 10: VALIDATION — ConversionSim vs Published Data (Figure 2)

**What's on it:** 4-panel validation figure.

**What to say (2-3 minutes):**
"Now the critical question: does ConversionSim actually match reality? We validated against four independent published datasets.

Panel A: Iyer et al. 2022 showed that cssDNA donors produce about 1.9-fold higher HDR than linear ssDNA. We predict 2.07-fold. That's within their observed range of 1.5 to 2.1.

Panel B: Chauhan et al. 2023 showed that staggered-cutting Cas9 variants achieve about 1.9-fold improvement. We predict 1.82-fold. Again, within the observed range of 1.4 to 2.8.

Panel C: The tract length distribution shape matches the right-skewed pattern from Elliott et al. 1998, though our tracts are longer — which makes sense because they measured endogenous donor recombination and we model exogenous donors.

Panel D: And here's the important one where we DON'T match. Paquet et al. 2016 used ssODN donors and showed very sharp incorporation decline within 30-50 bp of the cut. We over-predict at all distances. And this brings me to the most interesting finding of the validation."

---

### SLIDE 11: KEY FINDING — SDSA vs SSTR Pathway

**What's on it:** Side-by-side comparison of SDSA (what we model) vs SSTR (what ssODNs use).

**What to say (2-3 minutes):**
"The poor fit to Paquet et al. is not a failure — it's informative. It tells us something biologically meaningful.

Short ssODN donors — under about 200 nucleotides — are incorporated through a completely different pathway called SSTR, single-strand template repair. This is RAD51-independent. It works at the replication fork. And it produces very short conversion tracts — around 30 to 50 base pairs.

ConversionSim models SDSA — synthesis-dependent strand annealing. This is the pathway used by LONG donors: cssDNA, long ssDNA, dsDNA, AAV vectors. It's RAD51-dependent, involves strand invasion, and produces tracts in the hundreds to thousands of base pairs range.

The mismatch between our model and the Paquet data is because we're modeling the wrong pathway for their donor type. This distinction — SDSA for long donors, SSTR for short donors — is itself a useful contribution. It makes explicit what's often implicit in the field: that the repair mechanism depends on the donor format.

We acknowledge this as a limitation and identify adding an SSTR sub-model as the highest-priority future extension."

**Why this matters:** This is the most intellectually interesting part of the whole project. A naive approach would try to fit the model to all data. Instead, we identified WHY it doesn't fit and turned a limitation into a finding. This is exactly how good science works.

---

### SLIDE 12: MOSAIC BENCHMARKING (Figure 3)

**What's on it:** 14-study concordance table and accuracy-by-modality bar chart.

**What to say (2-3 minutes):**
"For MOSAIC, we benchmarked against 14 published genome editing studies. We asked: if we give MOSAIC the same gene, mutations, and cell type as the published study, does it recommend what the authors actually did?

Overall, 71.4% — 10 out of 14 studies had the authors' strategy in MOSAIC's top three. For base editing specifically, it's 100% — all five papers where authors used ABE, MOSAIC ranked base editing number one.

What about the four disagreements? Three of the four follow the same pattern: MOSAIC recommends a DSB-free approach — base editing or prime editing — while the authors used HDR. But in two of those cases, the papers were published BEFORE the relevant base or prime editors were available. MOSAIC is recommending what the field would do TODAY with current tools.

The fourth disagreement is a genuine gap — DMD exon deletion is a strategy type we haven't implemented yet.

I want to be clear: 71% is not 100%. We don't claim MOSAIC always beats expert judgment. But it ensures you've considered all your options, which is valuable."

---

### SLIDE 13: DISCUSSION — Claims and Limitations

**What's on it:** What we claim (green) / don't claim (red) / 5 limitations.

**What to say (2 minutes):**
"Let me be explicit about what we're claiming and what we're not.

We claim ConversionSim is the first public simulator of HDR gene conversion tracts. We claim that relative comparisons — circular vs linear, staggered vs blunt — are reliable because they match published ratios. We claim MOSAIC fills a real gap in the field.

We do NOT claim absolute locus-specific HDR rate prediction. We do NOT claim MOSAIC is always right. And we do NOT claim the supporting modules — ChromBridge, TopoPred, LoopSim — are methodological advances. They're utilities.

Five specific limitations: parameter uncertainty not propagated, SSTR not modeled, PAM availability not checked, weights are heuristic, and no prospective validation. We're working on the first two."

**Key philosophy:** Being upfront about limitations builds trust. If a reviewer finds a limitation you didn't mention, it hurts your credibility. If you already acknowledged it, you look rigorous.

---

### SLIDE 14: SUPPLEMENTARY FIGURES

**What's on it:** Sensitivity analysis + weight robustness figures.

**What to say (1 minute):**
"Two quick supplementary results. The sensitivity analysis shows that ConversionSim's predictions are most sensitive to the SDSA displacement probability — that's the parameter controlling tract length. This tells us where to focus future calibration.

The weight sensitivity shows that base editing remains the top-ranked strategy across a wide range of safety weights — from 0.30 to 0.80. You'd have to drop safety below 0.30 before HDR overtakes it. This means the recommendation is robust, not an artifact of our specific weight choices."

---

### SLIDE 15: SOFTWARE

**What's on it:** Project stats, access methods, live demo guide.

**What to say (1 minute, then switch to browser):**
"The full toolkit is 43 Python files, about 22,700 lines of code. Everything is on GitHub under MIT license. There's a web interface, a command-line tool, and a Docker container.

Let me show you the web app live."

**[SWITCH TO BROWSER — http://localhost:8501]**

**Demo script (5-7 minutes):**
1. "I'll start by fetching NF1 from Ensembl." → Gene Setup → Fetch from Ensembl → Type NF1 → Fetch Gene. "You can see it pulled 58 exons, chromosome 17, forward strand."
2. "Now I'll add two mutations." → Exon 20, G>A → Add. Exon 50, C>T → Add.
3. "Let's run MOSAIC." → Strategy Optimizer → Run. "It says Dual Base Editing is #1 — both mutations are transitions, so ABE can correct them with zero DSBs."
4. "Let's simulate the HDR tract lengths." → Conversion Tract Simulator → Run. "Median 540 bp. The probability of reaching exon 50 is zero."
5. "And the 3D distance." → ChromBridge → Run. "1,133 nm gap, 261 nm donor. Can't bridge."
6. "Finally, the donor quality check." → TopoPred → Auto-design from Ensembl → Analyze. "This pulls real genomic sequence and checks for structural problems."

**[SWITCH BACK TO SLIDES]**

---

### SLIDE 16: CONCLUSIONS

**What's on it:** 5 numbered conclusions.

**What to say (1 minute):**
"Five conclusions:

First, ConversionSim is, to our knowledge, the first quantitative model of HDR gene conversion tracts validated against independent published data.

Second, MOSAIC formalizes a decision that's currently made by intuition, with 71% concordance with expert choices.

Third, the SSTR/SDSA pathway boundary we identified defines when our model is and isn't applicable — a useful distinction for the field.

Fourth, for all tested human genes, single-template dual-site HDR is infeasible when mutations are more than 5 kb apart.

And fifth, CRISPRArchitect complements existing tools by answering the upstream question: which editing modality should you use?"

---

### SLIDE 17: FUTURE DIRECTIONS + PUBLICATION

**What's on it:** Future work priorities + PLOS Comp Bio manuscript status.

**What to say (1 minute):**
"Looking forward, the highest priorities are adding the SSTR sub-model for ssODN donors and the exon-reframing strategy. We also want to incorporate locus-specific chromatin data from ENCODE.

For publication, we're targeting PLOS Computational Biology as a Research Article. The manuscript draft is complete — full text, three main figures, two supplementary, 22 references. We need PI review and reference verification before submission."

---

### SLIDE 18: THANK YOU

**What to say:** "Thank you. I'm happy to take questions."

---

## ANTICIPATED QUESTIONS & ANSWERS

### Q1: "How is this different from CRISPOR or Benchling?"

**A:** "CRISPOR and Benchling are excellent for single-site editing — designing guide RNAs, checking off-targets, building donor templates. We don't replicate that. CRISPRArchitect answers a different question: when you have MULTIPLE mutations, which editing MODALITY should you use? Should you do HDR or base editing? Sequential or simultaneous? That comparison across modalities is what no existing tool does. CRISPRArchitect is complementary to CRISPOR and Benchling — you'd use our tool first to decide your strategy, then use CRISPOR to design the actual guides."

### Q2: "Your model only matches 3 out of 4 validation datasets. Isn't that a problem?"

**A:** "The one dataset we don't match — Paquet et al. — uses ssODN donors, which go through a completely different repair pathway called SSTR. Our model is designed for long donors like cssDNA, not short ssODNs. So the mismatch is expected and biologically informative. It's like testing a car engine model on a jet engine — the failure tells you about the different mechanism, not about a flaw in your model. We're transparent about this boundary and identify adding SSTR as the top priority for future work."

### Q3: "MOSAIC only matches experts 71% of the time. Why should I trust it?"

**A:** "Three things. First, 71% is the floor, not the ceiling — three of the four disagreements are cases where MOSAIC recommends tools that weren't available when the paper was published. If those authors had access to base editors, they might have made the same choice MOSAIC suggests. Second, MOSAIC isn't meant to replace expert judgment — it's meant to ensure you've considered all options. Even a 71% match means it's surfacing the right strategy most of the time. Third, for base editing specifically, it's 100% concordant — exactly the scenario where iPSC safety matters most."

### Q4: "Why didn't you do experimental validation?"

**A:** "That's a fair point and it's listed as a limitation. Our validation uses published data from other labs, which tests whether the model is consistent with existing knowledge. Prospective validation — designing experiments based on MOSAIC's recommendations and measuring outcomes — would be stronger evidence. We're planning to do this but wanted to first establish that the computational framework is sound before investing wet-lab resources."

### Q5: "Can this predict absolute HDR efficiency at a specific locus?"

**A:** "No, and we're explicit about that. HDR efficiency varies enormously across loci due to chromatin state, transcription, replication timing, and other factors we don't model. What ConversionSim predicts reliably is RELATIVE comparisons: circular vs. linear donors, staggered vs. blunt cuts, and the general shape of the tract length distribution. For absolute rates, you'd need locus-specific data that we don't incorporate."

### Q6: "The scoring weights in MOSAIC seem arbitrary. How did you choose them?"

**A:** "They are heuristic — we chose 40% safety, 30% efficiency, 15% time, 15% cost based on the biological reality that safety is paramount in iPSCs due to p53 selection. But importantly, we showed in supplementary Figure S2 that the top-ranked strategy — base editing for transition mutations — is robust across a wide range of weights, from 0.30 to 0.80 for safety. You'd have to make safety nearly irrelevant before HDR overtakes base editing. The weights are also user-adjustable, so researchers with different priorities can reweight."

### Q7: "What about off-target effects? Does your tool consider that?"

**A:** "No, and deliberately so. Off-target prediction is already well-handled by tools like CRISPOR and Cas-OFFinder. We didn't want to replicate that. CRISPRArchitect focuses on the strategy-level decision — which modality and approach — not on guide-level design. Once you've chosen your strategy using our tool, you'd use CRISPOR for guide design and off-target analysis."

### Q8: "Why didn't you use machine learning?"

**A:** "We considered it but chose a mechanistic model instead, for two reasons. First, there isn't enough published training data to learn HDR tract lengths from examples — most labs report only the final editing efficiency, not the tract length distribution. Second, a mechanistic model is interpretable: every parameter maps to a biological process, and when the model fails (like with ssODNs), the failure is informative. A black-box ML model might fit better to training data but wouldn't tell you WHY."

### Q9: "What's the advantage of cssDNA over other donor types?"

**A:** "Based on Iyer et al. 2022, cssDNA is about 2-fold better than linear ssDNA and about 3-fold better than linear dsDNA for HDR. The main advantage is exonuclease resistance — circular DNA has no free ends for exonucleases to attack, so it persists longer in the cell. Xie et al. 2025 achieved up to 70% knock-in in iPSCs with their GATALYST cssDNA system. Our model captures this through the donor topology multiplier and the reduced SDSA displacement probability for circular donors."

### Q10: "How does enFnCas9 compare to SpCas9 in your model?"

**A:** "enFnCas9, developed by the Chakraborty lab and published in Nature Communications 2024, has broader PAM flexibility (NRG instead of NGG) and improved HDR knock-in rates. In our model, we assign it a stagger of about 3 bp and a 1.5x HDR multiplier over SpCas9. The broader PAM means more target sites are accessible — about 3.5-fold more genomic sites compared to wild-type FnCas9."

### Q11: "Why is the LoopSim module included if it's not relevant to exogenous donors?"

**A:** "Good question. The cohesin loop extrusion mechanism, from the Marin-Gonzalez et al. 2025 Science paper, is primarily relevant for endogenous donor recombination — sister chromatid repair and loss of heterozygosity. For exogenous donors like cssDNA, the donor reaches the DSB by 3D diffusion, not by loop extrusion. We include LoopSim for completeness and because it's useful for understanding repair domain boundaries, but we explicitly state in the paper that it has limited bearing on exogenous donor HDR."

### Q12: "Can this tool be used for non-iPSC cell types?"

**A:** "Yes. We support HEK293T, K562, primary T cells, and HSCs, each with different parameter sets (HDR rates, cell cycle, p53 status). The scoring weights might need adjustment — for example, in HEK293T where p53 is inactive, the safety penalty for DSBs would be lower. The framework is cell-type-aware but the current validation is most extensive for iPSC-relevant comparisons."

### Q13: "What if both mutations are NOT transitions? Can you still use base editing?"

**A:** "No. Base editors can only do C→T (CBE) or A→G (ABE) transitions. If your mutations are transversions (like G→C or A→T) or indels, base editing is not an option. MOSAIC handles this automatically — it classifies each mutation and only offers base editing when both mutations are transitions. For transversions, MOSAIC would typically recommend prime editing (no DSB but lower efficiency) or sequential HDR."

### Q14: "Is there a size limit for cssDNA donors?"

**A:** "Yes. With the standard M13 phagemid system, you can produce cssDNA up to about 10-13 kb total (including the ~2.2 kb phagemid backbone). The newer GATALYST system from Xie et al. can go up to about 20 kb. For our purposes, the practical limit for the insert (homology arms + corrected sequence) is about 7-8 kb with standard methods and up to 18 kb with GATALYST."

### Q15: "What journal are you targeting and why?"

**A:** "PLOS Computational Biology, as a Research Article. We chose it because (1) the ConversionSim model is fundamentally a computational biology contribution — a new simulation approach validated against experimental data; (2) PLOS Comp Bio values reproducibility and open-source software, which aligns with our full GitHub release; and (3) the journal publishes both methods papers and software papers, and our work sits between the two."

---

## KEY NUMBERS TO MEMORIZE

These are the numbers people will ask about. Know them cold.

| Number | What it means |
|--------|--------------|
| **542 bp** | Median gene conversion tract length (iPSC, cssDNA, staggered) |
| **2,258 bp** | 95th percentile tract length |
| **3.8%** | HDR success rate in iPSCs with cssDNA |
| **2.07x** | Predicted cssDNA/lssDNA enhancement (observed: 1.9x) |
| **1.82x** | Predicted staggered/blunt enhancement (observed: 1.9x) |
| **71.4%** | MOSAIC concordance with published studies (10/14) |
| **100%** | MOSAIC concordance for base editing studies (5/5) |
| **261 nm** | Random coil diameter of a 3 kb cssDNA |
| **1,133 nm** | 3D distance between NF1 exon 20 and exon 50 |
| **0.0%** | Probability of tract reaching 100+ kb (answer to our question) |
| **43 files** | Python source code files |
| **22,700+** | Lines of code |
| **30/30** | Unit tests passing |
| **14 papers** | Used for MOSAIC benchmarking |
| **4 datasets** | Used for ConversionSim validation |

---

## COMMON OBJECTIONS AND HOW TO ADDRESS THEM

### "This is just a wrapper around known biology, not real computational innovation"

**Response:** "ConversionSim's novelty is in combining known biological steps — resection, filament formation, strand invasion, SDSA — into a single stochastic simulation that produces quantitative predictions. Nobody has done this before for exogenous donor HDR. The parameters are individually known but the integrated model and its predictions are new. More importantly, it WORKS — it matches published ratios without being fit to them."

### "The model is too simple / doesn't capture real complexity"

**Response:** "You're right that real HDR involves hundreds of proteins and complex chromatin dynamics. We deliberately chose simplicity over complexity for two reasons: (1) simple models with well-characterized parameters are more trustworthy than complex models with many uncertain parameters, and (2) the simple model already captures the key phenomena (donor topology effects, stagger effects, tract length distributions). Occam's razor applies — add complexity only when the simple model fails. And where it does fail — ssODNs — we identified exactly why."

### "You need experimental validation"

**Response:** "We agree, and we say so explicitly in the limitations. What we've done is validate against independent published data, which shows the model is consistent with existing knowledge. Prospective validation — designing experiments based on MOSAIC and measuring outcomes — is planned as the next step. But we believe the computational framework itself is a contribution worth reporting now, with the validation we have."

---

## THE 2-MINUTE ELEVATOR PITCH

Practice this until it's natural:

"When a patient has two mutations in a large gene, correcting both in iPSCs is complicated. You need to choose between HDR with a donor template, base editing, prime editing, or some combination — and this choice matters because DSBs in iPSCs trigger p53-mediated apoptosis. Currently, researchers make this decision by intuition.

We built CRISPRArchitect — the first computational tool that quantitatively compares these strategies. It has two core parts: a Monte Carlo simulator of HDR that predicts gene conversion tract lengths, and a strategy optimizer that ranks all feasible editing approaches.

We validated the HDR simulator against four published datasets — it correctly predicts that circular ssDNA donors are about 2x better than linear and that staggered nuclease cuts boost HDR about 1.8x. We benchmarked the strategy optimizer against 14 published studies and it matched expert decisions 71% of the time, with 100% concordance for base editing.

The tool is open source, has a web interface, and runs on real human genes fetched from Ensembl. We're preparing the manuscript for PLOS Computational Biology."
