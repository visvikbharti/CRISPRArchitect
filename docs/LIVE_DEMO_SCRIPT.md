# CRISPRArchitect v3 — Live Demo Script

## When to Use This

If during the lab meeting, someone says "can you show it working?" or the PI asks "show me," follow this script. It takes about 2 minutes.

## Setup (before the meeting)

Open a terminal and navigate:
```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
```

## Demo 1: COL7A1 — Base Editing Wins (30 seconds)

This is the case where CRISPRArchitect correctly identifies base editing as optimal.

```bash
python cli.py analyze --gene COL7A1 --hgvs "NM_000094.4:c.5047C>T" --cell iPSC
```

**What to say while it runs (~10 seconds):**
> "I'm giving it a COL7A1 nonsense variant — a C>T that creates a stop codon at position 1683. This causes epidermolysis bullosa. Let's see what the tool recommends."

**Expected output highlights:**
- Variant: COL7A1:c.5047C>T
- Consequence: stop_gained (nonsense)
- Top strategy: **Single-step Base Editing** (TOPSIS score ~0.98)
- Pareto: non-dominated
- Rank stability: >90%

**What to say:**
> "The tool recommends base editing — specifically ABE, because the reverse complement of this C>T mutation is a G>A, which ABE can correct. It scores 0.98 out of 1.0, with 90%+ rank stability. This means the recommendation holds under 9,000 out of 10,000 random weight combinations."

## Demo 2: Fetch a Gene (15 seconds)

If someone asks about a gene the lab works with:

```bash
python cli.py fetch --gene RPE65
```

**What to say:**
> "This pulls the full exon structure from Ensembl in real-time. 14 exons, chromosome 1. You can see the coordinates for each exon. This is the first step — the tool needs to know the gene structure to evaluate PAM sites and editing windows."

## Demo 3: ConversionSim (15 seconds)

```bash
python cli.py simulate --cell iPSC --nuclease enFnCas9 --donor cssDNA -n 10000
```

**What to say:**
> "This runs 10,000 Monte Carlo simulations of HDR with enFnCas9 and a cssDNA donor in iPSCs. Notice the confidence intervals — every number has uncertainty quantification. The HDR rate is about 3.5% with a 95% CI of 3.2 to 3.9%."

## Demo 4: HBB Sickle Cell Disease — Our Lab's Target (45 seconds)

This is directly relevant to the lab's SCD work and GMP/clinical trial.

```bash
python cli.py analyze --gene HBB --hgvs "NM_000518.5:c.20A>T" --cell CD34_HSC
```

**What to say while it runs:**
> "This is the classic sickle cell variant — HBB c.20A>T, the E6V transversion. Our lab is actively working on this. Let's see what CRISPRArchitect recommends."

**Expected output highlights:**
- Variant: HBB:c.20A>T (p.Glu7Val)
- Consequence: missense
- Base editing: **REJECTED** — A>T is a transversion, neither ABE nor CBE can correct it
- Top strategy: **Prime Editing** (PE can correct T back to A)
- HDR: feasible but penalized for DSB in HSCs

**What to say:**
> "The tool correctly rejects base editing — this is a transversion, not a transition. PE is recommended. Notice the delivery advisor now gives HSC-specific guidance: pre-stimulate with SCF/TPO/FLT3L, use nucleofection with mRNA, minimize culture to preserve stemness. This directly applies to our upcoming clinical work."

**Why this demo is powerful:** It connects CRISPRArchitect to the lab's own research program.

---

## Demo 5: Webapp (3 minutes, if time allows)

If the audience wants to see the visual interface:

### Launch
```bash
streamlit run webapp/app.py
```
Opens http://localhost:8501

### Quick webapp walkthrough
1. Sidebar → **Strategy Analysis**
2. Enter: Gene **COL7A1**, Chromosome **3**, Position **48580586**, Ref **C**, Alt **T**
3. Cell type: **iPSC**, Nuclease: **SpCas9**
4. Click **"Analyze Variant"** → wait ~10 seconds

**What to show:**
- **Variant annotation card** — consequence, HGVS, ref validation
- **TOPSIS ranking** — 6D scores, rank stability bars, Pareto badges
- **Delivery recommendations** — donor format, delivery method, viability tips
- **Feasibility breakdown** — all nuclease-editor combinations tested

**What to say about delivery:**
> "Every ranked strategy gets delivery annotations. For this iPSC case, it recommends ssODN for the base editing strategy, warns about p53 toxicity if dsDNA were used, and suggests ROCK inhibitor and p53DD co-delivery for viability. These come from a survey of 87 verified references."

### For sickle cell in webapp
- Gene: **HBB**, Chromosome: **11**, Position: **5227002**, Ref: **T**, Alt: **A**
- Cell type: **CD34_HSC**
- Shows HSC-specific delivery guidance (pre-stimulation, HiFi Cas9, minimize culture)

See `APP_DEMO_WALKTHROUGH.md` for the full step-by-step webapp guide.

---

## Tips

- **Don't apologize for the terminal interface.** The CLI is a feature, not a limitation. It shows the tool is a real computational pipeline, not a pretty GUI wrapper.
- **The webapp and CLI use the same pipeline.** Results are identical — the webapp just adds visualization.
- **If something fails** (Ensembl timeout, etc.): "The Ensembl API can be slow sometimes. Let me show you the pre-computed results instead." Then Alt+Tab to the figures.
- **If someone asks about a gene and it doesn't work:** It might be because Ensembl's gene symbol doesn't match. Try the HGNC symbol. Common issues: GBA (renamed to GBA1 in Ensembl).
- **Show the delivery section in the webapp.** It's the newest feature and demonstrates practical lab guidance.
