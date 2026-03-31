# Wednesday Lab Meeting — Complete Preparation Guide

## The Night Before

### 1. Test the presentation
- Open `paper/CRISPRArchitect_v3_LabMeeting.pptx` on the laptop you'll present from
- Check that all 26 slides render correctly (fonts, colors, layout)
- If using lab projector: test the PPTX on the projector — dark backgrounds can look washed out on old projectors. If so, increase room darkness.

### 2. Have figures ready in a separate window
Open these 6 PNG files in Preview/image viewer, ready to Alt+Tab:
```
paper/figures/v3_results/Fig_StrategyDistribution_v2_v3.png    → for slide 12
paper/figures/v3_results/Fig_BystanterFix.png                  → for slide 14
paper/figures/v3_results/Fig_CrossMethod_Agreement.png          → for slide 17
paper/figures/v3_results/Fig_ConversionSim_CIs.png             → for slide 18
paper/figures/v3_results/Fig_ParameterProvenance.png            → for slide 21
paper/figures/v3_results/Fig_LiteratureBenchmark.png            → for slide 22
```

### 3. Have a live demo ready (in case PI asks "show me")
In terminal, navigate to the crisprarchitect directory and prepare:
```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
# Test the CLI:
python cli.py fetch --gene COL7A1
# If asked for a full analysis, this takes ~10 seconds:
# python cli.py analyze --gene COL7A1 --hgvs "c.5047C>T"
```
Keep this terminal window open but minimized. Only show it if asked.

### 4. Print the speaker guide
Print `docs/LAB_MEETING_SPEAKER_GUIDE_v3.md` or keep it open on your phone/second screen. You won't read from it, but having the Q&A section accessible is reassuring.

---

## The 5 Numbers You Must Know By Heart

| Number | What It Means | When to Use It |
|---|---|---|
| **6 out of 30** | BE cases rescued by multi-nuclease engine | Slide 12, the headline result |
| **86.7%** | Top-1 accuracy on ClinVar benchmark | Slide 22 |
| **80%** | Top-3 concordance on literature benchmark | When discussing real-world validation |
| **100%** | Cross-method concordance (TOPSIS = VIKOR = WPM) | When defending scoring methodology |
| **224** | Tests passing, 0 failures | When discussing code quality |

---

## Anticipating Your PI's Thinking

Prof. Chakraborty said "discuss objectives first." This means he wants to understand:

### What he'll want to know:
1. **"What is this tool FOR?"** → It recommends the optimal editing strategy for any pathogenic variant, comparing BE/PE/HDR across multiple nucleases including enFnCas9.

2. **"How does this help OUR lab?"** → Two ways: (a) Before starting any iPSC editing experiment, run CRISPRArchitect to check if a simpler/safer strategy exists. (b) enFnCas9 is highlighted as expanding the base editing-accessible space — this is a unique selling point for the lab's nuclease.

3. **"Is this publishable?"** → Yes, targeting PLOS Computational Biology. Needs experimental validation (which the lab can provide) and independent expert labeling. The computational framework is complete.

4. **"What do you need from us?"** → Slide 24: (a) 3-5 iPSC editing cases with outcomes, (b) enFnCas9 cut-site stagger measurement, (c) his input on benchmark case labeling.

5. **"What's the timeline?"** → Computational work is 90% done. Manuscript can be drafted in 2-3 weeks. Experimental validation (if done) adds 1-2 months. Submission by May-June 2026 is realistic.

### What might concern him:
- **"30% top-1 concordance on literature benchmark"** → Be ready: "80% top-3, and most disagreements are with pre-PE HDR papers. The tool recommends the safer modern alternative."
- **"You found fabricated citations?"** → Be honest: "Yes, from AI-assisted writing. We caught every error through systematic verification. This demonstrates the rigorous auditing process."
- **"The enFnCas9 stagger is assumed"** → "Yes, and that's exactly why we need the lab's data. One measurement replaces our biggest assumption."

---

## Making the Project Appealing and Useful

### The Pitch (30-second version for the PI):
> "CRISPRArchitect is a tool that compares base editing, prime editing, and HDR side-by-side for any variant, using our own enFnCas9 as a first-class nuclease. For every iPSC editing experiment this lab starts, running CRISPRArchitect first takes 30 seconds and might reveal a simpler or safer strategy. The tool found that PAM-window constraints — not mutation type — are the real bottleneck for base editing, and that enFnCas9 helps overcome this bottleneck."

### Why the Lab Should Care:
1. **Saves experimental time**: Before committing to a 3-month HDR experiment, check if BE or PE would work at that locus.
2. **Showcases enFnCas9**: The tool demonstrates quantitatively how NRG PAM expands the editable space — this is a story that strengthens enFnCas9 publications.
3. **Publication value**: Computational tool papers in PLOS Comp Bio are well-cited (CRISPOR: >4,000 citations).
4. **Practical for the whole lab**: Anyone can use the Streamlit webapp or CLI.

### The "Wow" Moment in the Presentation
Slide 12 (BE rescue) is your strongest moment. The visual of going from 0/30 to 6/30 base editing cases, with enFnCas9 as a primary driver, directly connects to the lab's identity. Pause here and say: **"This rescue is possible because of the broader PAM — because of enFnCas9, developed in this lab."**

---

## If Things Go Wrong

### "The projector doesn't show the dark slides well"
→ Increase room darkness. The slides are designed for dark rooms.

### "Someone asks about a gene you haven't tested"
→ "Let me run it right now." Open terminal, type:
```
python cli.py fetch --gene [GENE_NAME]
```
This shows exon structure in real-time. If they want a full analysis:
```
python cli.py analyze --gene [GENE_NAME] --hgvs "c.[POSITION][REF]>[ALT]"
```

### "PI asks a question you can't answer"
→ "That's a great question and I want to give you an accurate answer. Let me check the documentation and get back to you." This is better than guessing.

### "Someone challenges the 30% concordance"
→ "You're right that 30% top-1 looks low. But look at WHY the cases disagree — 5 of 7 misses are HDR papers published before 2019, before prime editing existed. Our tool is recommending the safer modern alternative. The more informative metric is 80% top-3: the published strategy is almost always in our recommendation set."

### "PI asks 'why not just use CRISPOR?'"
→ "CRISPOR is excellent for guide RNA design — we don't replace it. But CRISPOR doesn't compare modalities. It doesn't tell you whether to use base editing vs prime editing vs HDR at a given locus. That's what CRISPRArchitect does."

---

## After the Meeting

1. **Take notes on every question and suggestion** — these become paper revisions
2. **If PI asks for specific experiments** — note exactly what, with gene names and variants
3. **If PI suggests collaborators for independent labeling** — follow up within 24 hours
4. **Send a follow-up email summarizing action items** within the same day

---

## Quick Reference: File Locations

| What | Where |
|---|---|
| Presentation | `paper/CRISPRArchitect_v3_LabMeeting.pptx` |
| Speaker guide | `docs/LAB_MEETING_SPEAKER_GUIDE_v3.md` |
| Figures | `paper/figures/v3_results/*.png` |
| Live demo | `python cli.py analyze --gene COL7A1 --hgvs "c.5047C>T"` |
| Webapp | `python cli.py webapp` (opens localhost:8501) |
| Data request | `docs/DATA_REQUEST_FOR_PI.md` |
| Data sheets | `docs/CRISPRArchitect_DataCollection.tsv` |
| Full documentation | `docs/COMPLETE_PROJECT_DOCUMENTATION.md` |
| GitHub | github.com/visvikbharti/CRISPRArchitect |
