# CRISPRArchitect Webapp — Step-by-Step Demo Walkthrough

**For:** Wednesday Lab Meeting live demonstration
**Time needed:** 5-8 minutes (or shorter segments per page)
**Launch command:** `python cli.py webapp` (opens http://localhost:8501)

---

## Before You Start

1. Open terminal, navigate to the project:
   ```bash
   cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
   python cli.py webapp
   ```
2. Wait for `You can now view your Streamlit app in your browser` message
3. Open http://localhost:8501 in Chrome/Safari
4. The app opens on the **Home** page

---

## Page 1: Home (show for 15 seconds)

**What you see:** Title "CRISPRArchitect" with four module cards (MOSAIC, ConversionSim, ChromBridge, cssDNA-TopoPred).

**What to say:**
> "This is the web interface. It has seven analysis pages accessible from the sidebar on the left. Let me walk through a real example — I'll use a real gene fetched from Ensembl."

**What to click:** Nothing. Move to the next page via the sidebar.

---

## Page 2: Gene & Mutation Setup (main demo page — spend 2-3 minutes here)

**Step-by-step:**

### Step 2a: Select "Fetch real gene from Ensembl (GRCh38)"
- In the sidebar, click **Gene & Mutation Setup**
- You'll see a radio button: "How would you like to define the gene structure?"
- Select **"Fetch real gene from Ensembl (GRCh38)"**

### Step 2b: Type a gene name
- In the "Gene symbol" text box, type: **NF1**
- Click **"Fetch from Ensembl"**
- Wait ~5 seconds for the API call

**What to say while waiting:**
> "The tool is fetching the real exon structure of NF1 from the Ensembl database — GRCh38 coordinates. NF1 has 58 exons spanning about 290 kilobases on chromosome 17."

**What you see:** A gene structure visualization with exon boxes on a line, plus a summary showing exon count, gene span, and coordinates.

### Step 2c: Set cell type and nuclease
- Under "Cell type," select **iPSC**
- Under "Nuclease," select **enFnCas9** (this is the default)

**What to say:**
> "I'm selecting iPSC as the cell type — this activates the p53-aware safety scoring. And enFnCas9 as the primary nuclease, because that's our lab's engineered nuclease with the broader NRG PAM."

### Step 2d: Add mutations
- Under "Add Mutation," set:
  - Exon number: **20**
  - Position within exon: leave as default or type **50**
  - Reference allele: **G**
  - Alternate allele: **A**
- Click **"Add Mutation"**
- Add a second mutation:
  - Exon number: **50**
  - Reference allele: **C**
  - Alternate allele: **T**
- Click **"Add Mutation"**

**What to say:**
> "I'm defining two mutations — a G>A transition in exon 20 and a C>T transition in exon 50. In a real scenario, these would be compound heterozygous pathogenic variants from a patient."

**What you see:** The gene visualization now shows two red markers at the mutation positions. Below, a mutation table shows each mutation with its classification (transition, base-editable: Yes).

---

## Page 3: MOSAIC Strategy Optimizer (the core demo — 2 minutes)

**What to click:** In the sidebar, select **"Strategy Optimizer (MOSAIC)"**

### Step 3a: Run the analysis
- Click the **"Run MOSAIC Analysis"** button
- Wait ~2 seconds for enumeration and scoring

**What you see:**
1. **Mutation Classification Table** — shows each mutation's type, base-editable status, prime-editable status
2. **Inter-Site Distance Analysis** — genomic distance between exon 20 and 50, 3D nuclear distance
3. **Strategy Ranking Table** — the key output! Shows 6-8 strategies ranked by score

**What to say:**
> "Here's the strategy ranking. You can see the tool has enumerated every feasible combination — dual base editing, dual prime editing, sequential HDR, hybrid approaches. Each strategy is scored on four dimensions: efficiency, safety, time, and cost."
>
> "Dual base editing ranks first because both mutations are transitions, and base editing doesn't require any DSBs — which is critical in iPSCs where p53 kills most DSB-bearing cells."

### Step 3b: Expand a strategy (if time allows)
- Click on any strategy row to expand its details
- Shows: description, efficiency estimate, safety assessment, required reagents

**What to say:**
> "Each strategy comes with a full description including what reagents you need, how many rounds of editing, and the estimated screening effort. The tool also tells you explicitly why rejected strategies were rejected — for example, if single-template HDR was rejected because the two mutations are too far apart."

---

## Page 4: ConversionSim (show if someone asks about Monte Carlo — 1 minute)

**What to click:** In sidebar, select **"Conversion Tract Simulator"**

### Step 4a: Configure simulation
- Cut type: **Staggered 5'** (already selected for enFnCas9)
- Overhang: **3 bp**
- Donor: **Circular ssDNA**
- Homology arms: **300 bp**
- Simulations: **10,000** (default)

### Step 4b: Run
- Click **"Run Simulation"**
- Takes ~1 second

**What you see:**
1. **Tract length histogram** — right-skewed distribution with median ~450 bp
2. **Survival curve** — P(tract >= distance) vs distance
3. **Summary statistics** with mean, median, percentiles
4. **Conversion probability bar chart** at key distances (100, 200, 500, 1000 bp)

**What to say:**
> "This is a Monte Carlo simulation of 10,000 virtual cells undergoing HDR. Each cell gets a random resection length, random RAD51 filament coverage, and a random synthesis tract from a geometric distribution. The output tells us: if I place an edit 500 bp from the cut site, about 47% of HDR events will incorporate it."
>
> "In version 3, every number comes with a standard error and 95% confidence interval — you can see them in the summary."

---

## Page 5: ChromBridge 3D Distance (show if someone asks about spatial biology — 1 minute)

**What to click:** Sidebar → **"3D Distance & Risk (ChromBridge)"**

**What you see:**
1. **3D distance prediction** — physical distance between the two mutation sites in nanometers
2. **Donor bridgeability** — can a 3 kb cssDNA bridge the gap? (answer: No, for distant exons)
3. **Scale diagram** — visual showing donor size vs. inter-locus distance
4. **Translocation risk assessment** — probability of deletion, inversion, translocation

**What to say:**
> "ChromBridge uses polymer physics to predict the 3D nuclear distance between the two sites. For exon 20 and exon 50 in NF1 — about 123 kilobases apart — the 3D distance is about 1,100 nanometers. A 3 kb cssDNA donor has a random coil diameter of only 260 nanometers. So a single template physically cannot bridge these two sites."
>
> "This is the quantitative evidence for why single-template dual-site HDR doesn't work — the donor is literally too small. This was actually the founding question that started the entire CRISPRArchitect project."

---

## Page 6: TopoPred Donor Quality (show briefly — 30 seconds)

**What to click:** Sidebar → **"Donor Quality Check (TopoPred)"**

**What to say:**
> "This module analyzes the secondary structure of a cssDNA donor template. It scans for G-quadruplexes and hairpins that could block RAD51 filament formation. If a homology arm is more than 40% folded, the tool warns you and suggests synonymous codon substitutions to disrupt the structures."

---

## Page 7: v2 Variant Analysis (the v3 pipeline page — show if specifically asked)

**What to click:** Sidebar → **"v2: Variant Analysis"**

This page uses the v2/v3 core pipeline with transcript-aware annotation.

### Step 7a: Enter a variant
- In the sidebar (which changes for this page):
  - Gene: **COL7A1**
  - Chromosome: **3**
  - Position: **48580586**
  - Ref: **C**
  - Alt: **T**
  - Cell type: **iPSC**
  - Nuclease: **SpCas9**

### Step 7b: Run
- Click **"Run v2 Pipeline"**
- Takes ~10 seconds (Ensembl API call)

**What you see:** Variant annotation card showing:
- HGVS c./p. notation
- Consequence: stop_gained (nonsense)
- Splice proximity
- Reference validation: PASS

Then strategy ranking cards showing the TOPSIS-scored strategies.

**What to say:**
> "This page uses the full v3 pipeline — transcript mapping, consequence annotation, PAM scanning across all nucleases, and TOPSIS ranking. For this COL7A1 nonsense variant, the tool recommends base editing as the top strategy, with a score of 0.98. This is because ABE can correct the reverse-complement G>A, and there's a guide with the right PAM that places the target in the editing window."

---

## When to Show Which Page

| If someone asks... | Show this page |
|---|---|
| "What does the tool do?" | Home (15 sec) |
| "How do I use it?" | Gene & Mutation Setup (2 min) |
| "What strategies does it recommend?" | MOSAIC Strategy Optimizer (2 min) |
| "How does the Monte Carlo work?" | ConversionSim (1 min) |
| "What about 3D distance?" | ChromBridge (1 min) |
| "What about donor design?" | TopoPred (30 sec) |
| "Show me the v3 pipeline" | v2: Variant Analysis (1-2 min) |
| "Can you run it on a specific gene?" | Gene & Mutation Setup → fetch from Ensembl |

---

## Troubleshooting During the Demo

| Problem | Solution |
|---|---|
| App won't start | `pip install streamlit plotly pandas` then retry |
| Ensembl fetch times out | "API is slow today. Let me use the demo gene." → Select "Generate demo gene" |
| No strategies generated | Probably need 2 mutations. Add a second mutation. |
| Page looks weird | Toggle "Dark Mode" in sidebar. Or resize browser window. |
| v2 page not available | v2 modules may not be on path. Use MOSAIC page instead. |

---

## Pro Tips

1. **Don't navigate randomly.** Follow the flow: Home → Setup → MOSAIC → (others if asked)
2. **Talk while things load.** Ensembl calls take 5-10 seconds. Explain what's happening.
3. **Use NF1 or COL7A1 for the demo.** These are well-characterized, multi-exon genes.
4. **If the PI suggests a gene, try it!** Type it into the Ensembl fetch. Most HGNC symbols work.
5. **Don't apologize for the interface.** It's a research tool, not a commercial product. The science matters more than the design.
6. **End the demo with the MOSAIC ranking.** That's the money shot — the ranked strategy table.
