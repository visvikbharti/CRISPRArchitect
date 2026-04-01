# CRISPRArchitect — Live Demo Script

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
- Transcript: ENST00000681320 (COL7A1) — 119 exons, reverse strand
- Parsed: NM_000094.4:c.5047C>T → chr3:48580586
- Consequence: nonsense
- HGVS c.: c.5047C>T, HGVS p.: p.Arg1683Ter
- Ref validation: PASS
- Top strategy: **Single-step Base Editing** (score ~0.983)
- Evidence tier: A

**What to say:**
> "The tool recommends base editing — specifically ABE, because the reverse complement of this C>T mutation is a G>A, which ABE can correct. It scores 0.983 out of 1.0 with evidence tier A. Prime editing is ranked second as a backup."

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
python cli.py analyze --gene HBB --hgvs "NM_000518.5:c.20A>T" --cell HSC
```

> **Note:** Use `--cell HSC` on the CLI. The webapp dropdown shows `CD34_HSC`.

**What to say while it runs:**
> "This is the classic sickle cell variant — HBB c.20A>T, the E6V transversion. Our lab is actively working on this. Let's see what CRISPRArchitect recommends."

**Expected output highlights:**
- Transcript: ENST00000335295 (HBB) — 3 exons, reverse strand
- Parsed: NM_000518.5:c.20A>T → chr11:5227002
- Consequence: missense
- HGVS c.: c.20A>T, HGVS p.: p.Glu7Val
- Base editing: **NOT offered** — A>T is a transversion, neither ABE nor CBE can correct it
- Top strategy: **Single-step Prime Editing** (score 1.000, evidence tier A)
- HDR: ranked second, penalized for DSB requirement in HSCs

**What to say:**
> "The tool correctly excludes base editing — this is a transversion, not a transition. Prime editing is recommended with a perfect score. HDR is feasible but penalized because DSBs in HSCs carry higher toxicity risk. This directly applies to our upcoming clinical work."

**Why this demo is powerful:** It connects CRISPRArchitect to the lab's own research program.

---

## Bonus Demos: If Someone Asks About a Specific Disease

These are pre-tested and ready. Use them if the PI or audience asks "what about DMD?" or "does it work for thalassemia?"

### DMD — Duchenne Muscular Dystrophy (30 seconds)

```bash
python cli.py analyze --gene DMD --hgvs "NM_004006.3:c.10108C>T" --cell iPSC
```

**What to say:**
> "This is a DMD nonsense variant — C>T creating a premature stop codon. The tool recommends base editing with a score of 0.996. DMD is the largest human gene — 79 exons, 2.4 megabases — and CRISPRArchitect handles it correctly, resolving the CDS position to chrX."

**Expected output:**
- Transcript: ENST00000357033 (DMD) — 79 exons, chrX, reverse strand
- Consequence: nonsense, p.Arg3370Ter
- Top strategy: **Single-step Base Editing** (score 0.996, tier A)

### Beta-Thalassemia — Codon 39 Nonsense (30 seconds)

This is one of the most common beta-thal mutations worldwide (Mediterranean, Middle Eastern populations).

```bash
python cli.py analyze --gene HBB --hgvs "NM_000518.5:c.118C>T" --cell HSC
```

**What to say:**
> "This is beta-thalassemia codon 39 — a C>T nonsense mutation, one of the most common worldwide. Since it's a C>T transition, base editing can directly correct it. The tool scores it 0.996 with tier A evidence."

**Expected output:**
- Consequence: nonsense, p.Gln40Ter
- Top strategy: **Single-step Base Editing** (score 0.996, tier A)

**Contrast with SCD demo:** Both are HBB mutations in HSCs, but the tool gives different recommendations — BE for beta-thal (transition) vs PE for SCD (transversion). This shows the tool is genuinely analyzing the variant, not just pattern-matching the gene.

### Beta-Thalassemia — IVS-I-5 Splice Variant (30 seconds)

The most common beta-thal mutation in the Indian subcontinent. Good to show if the PI asks about splice variants.

```bash
python cli.py analyze --gene HBB --hgvs "NM_000518.5:c.92+5G>C" --cell HSC
```

**What to say:**
> "This is IVS-I-5 — a splice variant 5 bases into intron 1 of HBB. It's the most common beta-thalassemia mutation in India. The tool correctly identifies it as a splice_region variant, maps it to the intronic position chr11:5226925, and recommends prime editing since base editing cannot handle this G>C transversion."

**Expected output:**
- Parsed: c.92+5G>C → chr11:5226925 (intronic)
- Consequence: splice_region
- Top strategy: **Single-step Prime Editing** (score 1.000, tier A)

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
- **TOPSIS ranking** — multi-dimensional scores, evidence tiers
- **Delivery recommendations** — donor format, delivery method, viability tips
- **Feasibility breakdown** — all nuclease-editor combinations tested

**What to say about delivery:**
> "Every ranked strategy gets delivery annotations. For this iPSC case, it recommends ssODN for the base editing strategy, warns about p53 toxicity if dsDNA were used, and suggests ROCK inhibitor and p53DD co-delivery for viability. These come from a survey of 87 verified references."

### For sickle cell in webapp
- Gene: **HBB**, Chromosome: **11**, Position: **5227002**, Ref: **A**, Alt: **T**
- Cell type: **CD34_HSC**
- Shows HSC-specific delivery guidance (pre-stimulation, HiFi Cas9, minimize culture)

See `APP_DEMO_WALKTHROUGH.md` for the full step-by-step webapp guide.

---

## Tips

- **Don't apologize for the terminal interface.** The CLI is a feature, not a limitation. It shows the tool is a real computational pipeline, not a pretty GUI wrapper.
- **The webapp and CLI use the same pipeline.** Results are identical — the webapp just adds visualization.
- **CLI uses `HSC`, webapp uses `CD34_HSC`.** Both map to the same cell type internally.
- **If something fails** (Ensembl timeout, etc.): "The Ensembl API can be slow sometimes. Let me show you the pre-computed results instead." Then Alt+Tab to the figures.
- **If someone asks about a gene and it doesn't work:** It might be because Ensembl's gene symbol doesn't match. Try the HGNC symbol. Common issues: GBA (renamed to GBA1 in Ensembl).
- **Show the delivery section in the webapp.** It's the newest feature and demonstrates practical lab guidance.
