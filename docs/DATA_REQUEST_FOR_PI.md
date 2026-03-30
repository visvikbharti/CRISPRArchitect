# Data Request for CRISPRArchitect v3 Validation

**To:** Prof. Debojyoti Chakraborty
**From:** Vishal Bharti
**Date:** March 2026
**Re:** Experimental data needed for CRISPRArchitect manuscript (targeting PLOS Computational Biology)

---

## Context

Sir, we have built CRISPRArchitect v3 — a computational tool that recommends
the optimal editing strategy (base editing vs. prime editing vs. HDR) for any
given pathogenic variant. The tool evaluates multiple nucleases including
enFnCas9 and uses TOPSIS multi-criteria decision analysis with sensitivity
analysis.

For the manuscript, we need **experimental validation data** to demonstrate that
the tool's recommendations match real-world editing outcomes. Data from our lab
would be uniquely valuable because we developed enFnCas9, and no other group
can provide direct comparisons.

---

## What We Need (in order of priority)

### Priority 1: Any iPSC Editing Cases (3-5 cases would transform the paper)

For any locus where the lab has performed genome editing in iPSCs (or other
cell lines), we need the information listed in the attached spreadsheet
(`CRISPRArchitect_DataCollection.tsv`). Even 2-3 cases would be sufficient.

**Ideal scenario:** Cases where multiple strategies were compared at the same
locus (e.g., enFnCas9 vs SpCas9, or base editing vs HDR).

### Priority 2: enFnCas9 Characterization Data

These would directly replace assumptions in our model with measurements:

1. **Cut-site stagger:** Does enFnCas9 produce blunt or staggered cuts?
   If staggered, what is the overhang length (in bp)?
   - Currently we assume 3 bp based on structural inference from FnCas9
     (Hirano et al., Cell, 2016). A direct measurement would be much stronger.
   - Method: Run-off sequencing or adapter ligation at cut sites.

2. **HDR efficiency comparison:** enFnCas9 vs SpCas9 at the same locus with
   the same donor in iPSCs. Even at 2-3 loci, this calibrates our
   `hdr_multiplier` parameter (currently assumed at 1.5x).

3. **enFnCas9-ABE8e editing window:** If the lab has tested ABE8e fused with
   enFnCas9, what editing window was observed? We currently extrapolate from
   SpCas9-ABE8e data (Tier B evidence).

### Priority 3: Published Data Pointers

If any of the above has been published in the Nat Commun 2024 paper
(Acharya et al., 15:5471) or other publications, please point me to the
specific supplementary tables and I will extract the data myself.

---

## What We Do NOT Need

- Raw sequencing files (FASTQ/BAM) — just the analyzed numbers
- Detailed protocols — just the key parameters and outcomes
- Unpublished data you prefer to keep confidential — published data is fine

---

## Why This Matters

Without experimental validation, the manuscript relies entirely on
computational benchmarks with self-defined "correct" answers. Reviewers at
PLOS Comp Bio will flag this as circular validation. Even a small amount of
experimental data fundamentally changes the narrative from "we predicted what
should work" to "we predicted what should work, and it matched what actually
worked in the enFnCas9-developing laboratory."

---

## Timeline

The manuscript is in advanced preparation. Data collection can proceed in
parallel with manuscript writing. Any data available within 2-4 weeks would
be incorporated before submission.

Thank you, Sir.

Vishal
