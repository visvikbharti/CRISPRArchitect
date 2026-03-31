# CRISPRArchitect Webapp — Step-by-Step Demo Walkthrough

**For:** Wednesday Lab Meeting live demonstration
**Version:** v3.1 (updated 2026-03-31)
**Time needed:** 8-12 minutes (or shorter segments per page)
**Features shown:** MOSAIC strategy optimizer, ConversionSim Monte Carlo, v3 multi-nuclease TOPSIS pipeline with delivery-aware recommendations

---

## Before You Start

1. Open terminal, navigate to the project:
   ```bash
   cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
   streamlit run webapp/app.py
   ```
   Alternative: `python cli.py webapp`

2. Wait for `You can now view your Streamlit app in your browser` message
3. Open http://localhost:8501 in Chrome/Safari
4. The app opens on the **Home** page
5. **Pre-flight check:** Ensure internet access (Ensembl API calls need it)

---

## Page 1: Home (15 seconds)

**What you see:** Title "CRISPRArchitect" with module cards.

**What to say:**
> "This is the web interface for CRISPRArchitect. It has four main analysis pages in the sidebar: Gene & Mutation Setup, MOSAIC strategy optimization, ConversionSim Monte Carlo simulation, and the Strategy Analysis page with multi-nuclease TOPSIS scoring and delivery-aware recommendations. Let me walk through a real example."

**What to click:** Move to the next page via the sidebar.

---

## Page 2: Gene & Mutation Setup (2-3 minutes)

### Step 2a: Select "Fetch real gene from Ensembl (GRCh38)"
- In the sidebar, click **Gene & Mutation Setup**
- Select **"Fetch real gene from Ensembl (GRCh38)"**

### Step 2b: Fetch a gene
- In the "Gene symbol" text box, type: **NF1**
- Click **"Fetch from Ensembl"**
- Wait ~5 seconds for the API call

**What to say while waiting:**
> "The tool is fetching the real exon structure of NF1 from the Ensembl database — GRCh38 coordinates. NF1 has 58 exons spanning about 290 kilobases on chromosome 17."

**What you see:** Gene structure visualization with exon boxes, plus a summary showing exon count, gene span, and coordinates.

### Step 2c: Set cell type and nuclease
- Under "Cell type," select **iPSC**
- Under "Nuclease," select **enFnCas9**

**What to say:**
> "I'm selecting iPSC as the cell type — this activates p53-aware safety scoring and the delivery advisor's iPSC-specific warnings. And enFnCas9 as the primary nuclease — our lab's engineered nuclease with the broader NRG PAM."

### Step 2d: Add mutations
- Mutation 1: Exon **20**, position **50**, Ref **G**, Alt **A** → Click **"Add Mutation"**
- Mutation 2: Exon **50**, Ref **C**, Alt **T** → Click **"Add Mutation"**

**What to say:**
> "I'm defining two mutations — a G>A transition in exon 20 and a C>T transition in exon 50. These represent compound heterozygous pathogenic variants."

**What you see:** Gene visualization with red markers at mutation positions. Mutation table with classification (transition, base-editable: Yes).

---

## Page 3: MOSAIC Strategy Optimizer (2 minutes)

**What to click:** Sidebar → **"Strategy Optimizer (MOSAIC)"**

### Step 3a: Run the analysis
- Click **"Run MOSAIC Analysis"** → wait ~2 seconds

**What you see:**
1. **Mutation Classification Table** — type, base-editable status, prime-editable status
2. **Inter-Site Distance Analysis** — genomic distance between exon 20 and 50
3. **Strategy Ranking Table** — 6-8 strategies ranked by score

**What to say:**
> "Here's the strategy ranking. The tool enumerated every feasible combination — dual base editing, dual prime editing, sequential HDR, hybrid approaches. Dual base editing ranks first because both mutations are transitions and base editing doesn't require DSBs — critical in iPSCs where p53 kills most DSB-bearing cells."

### Step 3b: Expand a strategy (if time)
- Click any strategy row for details: description, reagents, screening effort

---

## Page 4: ConversionSim Monte Carlo (1 minute — show if asked)

**What to click:** Sidebar → **"Conversion Tract Simulator"**

### Step 4a: Configure
- Cut type: **Staggered 5'** (for enFnCas9)
- Overhang: **3 bp**
- Donor: **Circular ssDNA**
- Homology arms: **300 bp**
- Simulations: **10,000**

### Step 4b: Run
- Click **"Run Simulation"** → ~1 second

**What you see:**
1. Tract length histogram (right-skewed, median ~450 bp)
2. Survival curve — P(tract >= distance) vs distance
3. Summary statistics with mean, median, SEs, 95% CIs
4. Conversion probability at key distances

**What to say:**
> "10,000 virtual cells undergoing HDR via SDSA. Each gets random resection, RAD51 filament coverage, and a geometric synthesis tract. Every number has a standard error and 95% CI. Note: this models SDSA only — not ssODN repair, which uses a different pathway."

**If asked about p=0.002 sensitivity:**
> "We swept p from 0.001 to 0.005 — a 10-fold range — and HDR recommendations are robust across all values. The figure is in the presentation slide 25."

---

## Page 5: Strategy Analysis (THE KEY DEMO — 3-4 minutes)

**What to click:** Sidebar → **"Strategy Analysis"**

This is the multi-nuclease, TOPSIS-scored pipeline with delivery-aware recommendations.

### Step 5a: Enter a variant
In the sidebar:
- Gene: **COL7A1**
- Chromosome: **3**
- Position: **48580586**
- Ref: **C**
- Alt: **T**
- Cell type: **iPSC**
- Nuclease: **SpCas9**

### Step 5b: Run
- Click **"Analyze Variant"**
- Takes ~10-15 seconds (Ensembl API call + multi-nuclease evaluation)

**What to say while loading:**
> "The pipeline is now fetching the COL7A1 transcript from Ensembl, normalizing the variant, running PAM scanning across all 5 nucleases, evaluating base editing with 9 editor profiles, prime editing, and HDR feasibility, then scoring with 6-dimensional TOPSIS and running 10,000 sensitivity permutations."

### Step 5c: Walk through results

**What you see (top to bottom):**

1. **Summary metrics** — strategies ranked, rejected, scoring method (TOPSIS), sensitivity runs (10,000)

2. **Variant annotation card:**
   - HGVS c./p. notation
   - Consequence (e.g., stop_gained/nonsense)
   - Reference validation: PASS
   - Splice proximity if relevant

   **What to say:**
   > "The variant is annotated automatically — consequence, HGVS notation, and reference allele validation against GRCh38."

3. **Strategy ranking (TOPSIS 6D):**
   - Each strategy shows: rank, TOPSIS score, Pareto status (optimal/dominated), rank stability bar, 6D dimension scores (Safety, Feasibility, Complexity, Risk, Confidence), evidence tier badge

   **What to say:**
   > "Strategies are ranked by TOPSIS — a formal multi-criteria decision method. The rank stability bar shows how robust this ranking is across 10,000 random weight permutations. Green stability means the recommendation is solid; orange means it's sensitive to weight choices."

4. **Delivery Recommendations (NEW):**
   - For each strategy: deliverability status, delivery method, complexity badge (1-5), donor format recommendation (ssODN/cssDNA/etc.), cell-type warnings, viability enhancer tips

   **What to say:**
   > "This is new — delivery-aware recommendations based on 87 verified references. The tool annotates each strategy with the recommended delivery method, optimal donor format, and cell-type-specific viability tips. For iPSCs, it warns about dsDNA toxicity via p53 and recommends ROCK inhibitor, BCL-XL, and cold shock for HDR. These are post-ranking annotations — they don't change the TOPSIS scores but add practical laboratory guidance."

5. **Feasibility breakdown:**
   - Base editing results for ALL nuclease-editor combinations
   - Prime editing: PBS/RT template lengths
   - HDR: donor type, cut-to-edit distance

   **What to say:**
   > "The feasibility section shows results for every nuclease-editor combination tested — not just SpCas9 but also enFnCas9, SpCas9-NG, and SpRY. This is how the multi-nuclease engine rescued base editing from 0% to 20% of top-ranked strategies."

6. **Rejected strategies** (collapsible) — with rejection reasons

### Alternative demo variant: HBB sickle cell
If someone asks about sickle cell disease:
- Gene: **HBB**
- Chromosome: **11**
- Position: **5227002**
- Ref: **T**
- Alt: **A**
- Cell type: **CD34_HSC** (new in v3.1!)
- Nuclease: **SpCas9**

**What to say:**
> "For the sickle cell E6V mutation in CD34+ HSCs, the delivery advisor recommends nucleofection with mRNA, suggests pre-stimulation with SCF/TPO/FLT3L, and warns about quiescence — HSCs need to be in S/G2 phase for HDR."

---

## Demo Flow Recommendations

### Short demo (5 min): Home → Setup (NF1) → MOSAIC → show ranking
### Standard demo (8 min): Add Strategy Analysis (COL7A1)
### Full demo (12 min): Add ConversionSim + HBB sickle cell in CD34_HSC

### When to show which page

| If someone asks... | Show this page |
|---|---|
| "What does the tool do?" | Home (15 sec) |
| "How do I use it?" | Gene & Mutation Setup (2 min) |
| "What strategies does it recommend?" | MOSAIC Strategy Optimizer (2 min) |
| "How does the Monte Carlo work?" | ConversionSim (1 min) |
| "Show me multi-nuclease / TOPSIS" | v3 Strategy Analysis (3 min) |
| "What about delivery?" | v3 Strategy Analysis → Delivery section |
| "Can you run it on a specific gene?" | Setup → Ensembl fetch, or v3 → enter coordinates |
| "What about sickle cell / HSCs?" | v3 → HBB variant with CD34_HSC |

---

## Troubleshooting During the Demo

| Problem | Solution |
|---|---|
| App won't start | `pip install streamlit plotly pandas requests` then retry |
| Ensembl fetch times out | "API is slow today. Let me use the demo gene." → Select "Generate demo gene" |
| No strategies generated | Probably need 2 mutations (MOSAIC page). Add a second mutation. |
| v3 page shows import error | Check that `core/` is on Python path. Run from project root. |
| Page looks weird | Resize browser. Toggle theme in Streamlit Settings (top-right menu). |
| Delivery section missing | Pipeline ran before delivery advisor was integrated. Re-click "Analyze Variant". |

---

## Pro Tips

1. **Follow the flow:** Home → Setup → MOSAIC → v3 Analysis. Don't skip around.
2. **Talk while things load.** Ensembl calls take 5-10 seconds. Explain the pipeline stages.
3. **Use NF1 for MOSAIC, COL7A1 for v3.** These are well-characterized, multi-exon genes.
4. **Show the delivery section.** It's new and demonstrates practical lab guidance — the PI will appreciate this.
5. **If the PI suggests a gene, try it!** Most HGNC symbols work with Ensembl fetch.
6. **End with the ranking + delivery.** That's the money shot — TOPSIS ranking with practical delivery guidance.
7. **Don't apologize for the interface.** It's a research tool. The science matters more than the design.
