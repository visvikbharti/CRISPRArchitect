# CRISPRArchitect — Live Demo Script for Lab Meeting

## Before the Meeting (5 minutes before)

Open Terminal and run:
```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
streamlit run webapp/app.py
```

Open your browser to **http://localhost:8501** and verify it loads.

Keep Terminal running in the background. Do NOT close it.

---

## STEP 1: Set the Stage (Home Page)

**What you see:** The Home page with CRISPRArchitect title and 4 module cards.

**What to say:**

"This is CRISPRArchitect — a web-based toolkit we built. It has four analysis modules. Let me show you how it works with a real human gene."

**What to click:**
- Toggle **Dark Mode** ON in the sidebar (looks better on projector)
- Click **"Gene & Mutation Setup"** in the sidebar

---

## STEP 2: Fetch a Real Gene (Gene & Mutation Setup Page)

**What you see:** Gene & Mutation Setup page with input options.

**What to say:**

"I'm going to fetch NF1 — a 58-exon gene on chromosome 17 that causes neurofibromatosis. This is very close to the 60-exon scenario we were investigating. The gene coordinates come directly from Ensembl GRCh38 — no genome download needed."

**What to click:**
1. Under gene input mode, select **"Fetch real gene from Ensembl (GRCh38)"**
2. From the dropdown, pick **NF1**
3. Click **"Fetch Gene"**

**What appears:** A green success box saying:

> NF1 fetched: chr17, 58 exons, forward strand, span 287 kb

**Point at the screen and say:**

"58 exons. Median exon size 135 base pairs. This took 2 seconds — it's a live API call to Ensembl."

---

## STEP 3: Add Two Mutations (Same Page, Scroll Down)

**Stay on the same page.** Scroll down to "Mutation Definition."

**What to say:**

"Now I'll add two mutations. Exon 20 — a G to A change, that's a transition. And exon 50 — C to T, also a transition. Both are the kind of mutation a patient might carry."

**What to click:**
1. Set **Exon number: 20**, **Reference: G**, **Alternate: A** → Click **"Add Mutation"**
2. Set **Exon number: 50**, **Reference: C**, **Alternate: T** → Click **"Add Mutation"**

**Check the sidebar** — it should now say:

> Gene defined: Yes
> Mutations: 2

**What to say:**

"The sidebar confirms we have a gene and two mutations. The genomic distance between exon 20 and exon 50 is about 123 kilobases — that's all introns. The exonic distance is only 5 kb. This is the fundamental problem: the donor template sees the genomic distance, not the exonic distance. Now let's ask the tool: what's the best way to correct both?"

---

## STEP 4: Run MOSAIC Strategy Optimizer

**What to click:** Click **"Strategy Optimizer (MOSAIC)"** in the sidebar.

**What to say:**

"MOSAIC enumerates every possible editing strategy and ranks them by efficiency, safety, time, and cost."

**What to click:** Click **"Run MOSAIC Analysis"** button.

**What appears:** A ranked table of strategies. The #1 strategy will be **DUAL_BASE_EDITING**.

**Point at the #1 ranked strategy and say:**

"Dual Base Editing is ranked number one with a score of 0.756. Why? Because both mutations are transitions — G to A can be fixed by adenine base editing, C to T by cytosine base editing. And the crucial point — zero double-strand breaks. In iPSCs, that means no p53 selection risk. Ihry et al. and Haapaniemi et al. in 2018 showed that DSBs in iPSCs select for p53-mutant cells. Base editing avoids this entirely.

Look at the bottom of the list — Simultaneous Dual HDR is ranked last. It requires two DSBs, two donors, and the combined efficiency drops to about 4%. MOSAIC correctly identifies this as the worst option for iPSCs."

**If there's a bar chart**, point to it:

"This chart compares all strategies across four dimensions. Base editing wins on safety — that bar is at 100%."

---

## STEP 5: Run ConversionSim Tract Simulator

**What to click:** Click **"Conversion Tract Simulator"** in the sidebar.

**What to say:**

"But what if the mutations were transversions — where base editing doesn't work? Then you'd need HDR with a donor template. The question becomes: how far does the donor information actually propagate during repair? That's what ConversionSim simulates. This is a Monte Carlo simulation of the HDR process — end resection, RAD51 filament formation, strand invasion, and SDSA synthesis — step by step."

**What to click:**
- Keep all defaults (staggered 3bp, cssDNA, 300bp arms, iPSC)
- Set simulations to **10,000**
- Click **"Run Simulation"**

**What appears:** A histogram of tract lengths, a survival curve, and summary statistics.

**Point at the key numbers and say:**

"The median gene conversion tract is about 500 base pairs. That means half of all HDR events copy less than 500 bp of donor sequence.

At 1,000 bp, only about 25% of events reach that far.
At 2,000 bp, less than 7%.

Our two exons are 123 THOUSAND base pairs apart. The probability of the tract reaching from exon 20 to exon 50 is zero — not low, literally zero. This is a fundamental biological constraint. The SDSA repair pathway collapses the D-loop after a few hundred base pairs of synthesis. No donor design, no concentration optimization, can overcome this."

---

## STEP 6: Run ChromBridge 3D Distance

**What to click:** Click **"3D Distance & Risk (ChromBridge)"** in the sidebar.

**What to say:**

"There's a second, completely independent reason this can't work — the physical distance in 3D nuclear space. This uses polymer physics — the same models used to interpret FISH and Hi-C data."

**What to click:**
- Genomic distance should auto-fill (~122,795 bp). If not, type it.
- Set **Donor size: 3000** bp
- Set **Donor type: Circular ssDNA**
- Click **"Run ChromBridge Analysis"**

**What appears:** 3D distance prediction and bridgeability analysis.

**Point at the numbers and say:**

"In the nucleus, these two exons are about 1,133 nanometers apart — that's 1.1 micrometers. A 3 kilobase cssDNA donor folds into a random coil of just 261 nanometers diameter. It's about 4 to 5 times too small to physically touch both sites at once.

Think of it this way — the donor is a marble, and it's trying to span a hallway. No amount of concentration will make the marble bigger."

**If translocation risk numbers appear:**

"ChromBridge also estimates that if you did make two simultaneous cuts at these positions, there's about a 0.3% chance of a large deletion between the sites. Small in absolute terms, but in clinical iPSC work, any chromosomal rearrangement risk matters."

---

## STEP 7 (Optional): Donor Quality Check

**Only do this step if you have time and the audience is engaged. If they look ready for Q&A, skip directly to the Closing below.**

**What to click:** Click **"Donor Quality Check (TopoPred)"** in the sidebar.

**What to say:**

"One last thing — if we design a cssDNA donor for just one of these mutations — say exon 20 — we can check whether the donor sequence has structural problems that would impair HDR."

**What to click:**
- If you see "Auto-Design cssDNA Donor" at the top → select the exon 20 mutation → click the button
- Then click **"Analyze Donor"**

**What appears:** G-quadruplex scan results, hairpin analysis, accessibility scores per homology arm.

**What to say:**

"This pulled the actual genomic DNA sequence from Ensembl — the real intronic sequences flanking exon 20 — and checked for G-quadruplexes and hairpins. These are secondary structures in single-stranded DNA that can block RAD51 from loading onto the homology arms. If the accessibility score drops below 60%, the donor should be redesigned with synonymous codon substitutions to break those structures."

---

## CLOSING: Switch Back to PowerPoint Slides

**Switch from browser back to your presentation. Go to Slide 16 (Conclusions).**

**What to say (this is your punchline — deliver it clearly):**

"So in 5 minutes we just:

1. Fetched a real 58-exon gene from the human genome in 2 seconds
2. Asked MOSAIC for the optimal strategy — it said base editing, zero DSBs, safest for iPSCs
3. Simulated HDR and showed the gene conversion tract reaches only 500 base pairs — not 123 kilobases
4. Showed the cssDNA donor is physically 4 times too small to bridge the 3D nuclear gap

The answer to our original question is definitively no — a single cssDNA cannot correct two distant mutations simultaneously. But now we have a quantitative, validated tool that tells researchers what they SHOULD do instead.

When we benchmarked MOSAIC against 14 published editing studies, it agreed with expert decisions 71% of the time — and 100% of the time for base editing, which is exactly the scenario where iPSC safety matters most.

This is open source on GitHub. The manuscript is drafted for PLOS Computational Biology."

**Then proceed through:**
- **Slide 17** (Future Directions + Publication Plan) — brief, 30 seconds
- **Slide 18** (Thank You) — "Questions?"

---

## EMERGENCY BACKUP PLANS

### If Ensembl fetch fails (no wifi / API down):

Don't panic. Say: "The API seems to be slow today — let me use a pre-built demo gene instead."

**Do this:**
- Select **"Generate demo gene"** instead of Ensembl
- Set number of exons to 60
- Click Generate
- Everything else works exactly the same

### If the web app crashes or freezes:

Say: "Let me show you the pre-computed results instead."

Switch to these slides which have the same data:
- **Slide 8** = ConversionSim figure (tract distributions)
- **Slide 10** = Validation figure (4-panel comparison)
- **Slide 12** = MOSAIC benchmarking figure

### If a computation takes more than 10 seconds:

Say: "This normally runs in under a second. Let me show the pre-computed version while it loads."

Switch to the relevant slide.

### If someone asks to try a different gene:

Say: "Absolutely — type any gene symbol from Ensembl. Try DMD, BRCA2, or CFTR."

This works live. Let them see it fetch in real time — it's impressive.

---

## TIMING GUIDE

| Section | Time |
|---------|------|
| Slides 1-7 (Background → Methods) | 10 min |
| Slides 8-9 (Results + Real Genes) | 3 min |
| Slides 10-12 (Validation + Benchmarking) | 5 min |
| Slides 13-14 (Discussion + Supplementary) | 3 min |
| Slide 15 → LIVE DEMO (Steps 1-6 above) | 7 min |
| Slides 16-18 (Conclusions → Thank You) | 2 min |
| **TOTAL** | **30 min** |
| Q&A | 10-15 min |

---

## KEY NUMBERS TO HAVE READY

If anyone asks for a specific number, here they are:

| Number | What it is |
|--------|-----------|
| 542 bp | Median gene conversion tract (iPSC, cssDNA, staggered) |
| 2,258 bp | 95th percentile tract length |
| 3.8% | HDR success rate in iPSCs |
| 2.07x | Our predicted cssDNA/lssDNA ratio (observed: 1.9x) |
| 1.82x | Our predicted stagger/blunt ratio (observed: 1.9x) |
| 71.4% | MOSAIC concordance with 14 published papers |
| 100% | MOSAIC concordance for base editing papers (5/5) |
| 261 nm | cssDNA random coil diameter (3 kb donor) |
| 1,133 nm | 3D nuclear distance between NF1 exon 20 and 50 |
| 0.0% | Probability of tract reaching 123 kb |
| 43 files | Python source files |
| 22,700+ | Lines of code |
| 30/30 | Unit tests passing |
