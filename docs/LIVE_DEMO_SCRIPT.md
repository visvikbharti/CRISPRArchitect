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

## Tips

- **Don't apologize for the terminal interface.** The CLI is a feature, not a limitation. It shows the tool is a real computational pipeline, not a pretty GUI wrapper.
- **If something fails** (Ensembl timeout, etc.): "The Ensembl API can be slow sometimes. Let me show you the pre-computed results instead." Then Alt+Tab to the figures.
- **If someone asks about a gene and it doesn't work:** It might be because Ensembl's gene symbol doesn't match. Try the HGNC symbol. Common issues: GBA (renamed to GBA1 in Ensembl).
