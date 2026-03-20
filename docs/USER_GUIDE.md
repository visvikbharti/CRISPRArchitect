# CRISPRArchitect — Step-by-Step User Guide

## Table of Contents

1. [How to Launch the Web App](#1-how-to-launch-the-web-app)
2. [Step-by-Step Walkthrough (Demo Mode)](#2-step-by-step-walkthrough-demo-mode)
3. [Step-by-Step Walkthrough (Real Gene)](#3-step-by-step-walkthrough-real-gene)
4. [Where is the Human Genome Used?](#4-where-is-the-human-genome-used)
5. [Command-Line Usage (No Web App)](#5-command-line-usage-no-web-app)
6. [Troubleshooting](#6-troubleshooting)

---

## 1. How to Launch the Web App

### Prerequisites

```bash
# You need Python 3.9+ and these packages:
pip install numpy scipy matplotlib streamlit plotly pandas
```

### Launch

```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA/crisprarchitect
streamlit run webapp/app.py
```

Your browser will open to **http://localhost:8501**

### Dark Mode

There is a **Dark Mode toggle** at the top of the left sidebar. Switch it on/off anytime.

---

## 2. Step-by-Step Walkthrough (Demo Mode)

This walkthrough uses a computer-generated demo gene to show you every feature.

### STEP 1: Gene & Mutation Setup

1. Click **"Gene & Mutation Setup"** in the left sidebar
2. Under **Gene Input Mode**, select **"Generate demo gene"**
3. Set **Gene name** to anything (e.g., "MyGene")
4. Set **Number of exons** to **60** (use the slider)
5. Click **"Generate Demo Gene"** button
6. You should see a message: "Demo gene created with 60 exons" and a gene visualization (blue bars = exons, gaps = introns)

### STEP 2: Add Mutations

Still on the Gene & Mutation Setup page:

7. Scroll down to **"Add a Mutation"**
8. For **Mutation 1**:
   - Exon number: **20**
   - Reference allele: **G**
   - Alternate allele: **A**
   - Click **"Add Mutation"**
9. For **Mutation 2**:
   - Exon number: **50**
   - Reference allele: **C**
   - Alternate allele: **T**
   - Click **"Add Mutation"**
10. You should see "2 mutations defined" in the sidebar

### STEP 3: Select Cell Type & Nuclease

11. **Cell type**: select **iPSC** (default)
12. **Nuclease**: select **enFnCas9** (default)

> **What these mean biologically:**
> - iPSC = induced pluripotent stem cells (hard to edit, low HDR, p53 concerns)
> - enFnCas9 = engineered FnCas9 with improved HDR and PAM flexibility

### STEP 4: Run Strategy Optimizer (MOSAIC)

13. Click **"Strategy Optimizer (MOSAIC)"** in the sidebar
14. Click the **"Run MOSAIC Analysis"** button
15. You will see:
    - **Mutation Classification**: shows each mutation type (transition/transversion) and whether it's amenable to base editing or prime editing
    - **Distance Analysis**: genomic distance, exonic distance, intronic distance between the two mutation sites
    - **Strategy Ranking Table**: all feasible strategies ranked by overall score
    - **Bar Chart**: comparing strategies across efficiency, safety, time, and cost
    - **Detailed descriptions** of each strategy (click to expand)

> **How to interpret the ranking:**
> - Score closer to 1.0 = better
> - Green (#1 rank) = recommended strategy
> - Safety score is weighted heavily for iPSCs because of the p53 selection risk

### STEP 5: Run Conversion Tract Simulator

16. Click **"Conversion Tract Simulator"** in the sidebar
17. The inputs should be pre-filled based on your nuclease choice:
    - Cut type: staggered_5prime (for enFnCas9)
    - Overhang: 3 bp
    - Donor topology: circular_ssDNA
    - Homology arm length: 300 bp
    - Number of simulations: 10,000 (increase for more precision, slower)
18. Click **"Run Simulation"**
19. You will see:
    - **Summary statistics**: HDR success rate, median tract length
    - **Histogram**: distribution of gene conversion tract lengths
    - **Survival curve**: probability that the tract extends at least X bp from the cut
    - **Key distance probabilities**: e.g., "P(reaching 500 bp) = 52%"
    - **"Can HDR reach?" box**: if you have a gene defined, it shows whether HDR from one mutation site can reach the other (answer: NO for distant sites)

> **How to interpret:**
> - Median tract ~500 bp means half of HDR events copy less than 500 bp of donor
> - The survival curve drops steeply — at 2 kb, only ~7% of events reach that far
> - For sites separated by 100+ kb: probability is effectively 0

### STEP 6: Run 3D Distance & Risk Analysis

20. Click **"3D Distance & Risk (ChromBridge)"** in the sidebar
21. The genomic distance should be auto-filled from your gene
22. Set **Donor size** (e.g., 3000 bp for a 3 kb cssDNA)
23. Set **Donor type** (e.g., Circular ssDNA)
24. Click **"Run ChromBridge Analysis"**
25. You will see:
    - **3D Distance Prediction**: estimated physical distance in nanometers
    - **Bridgeability Analysis**: can your donor physically reach both sites? (Includes donor coil diameter vs inter-locus distance)
    - **Translocation Risk**: if making two DSBs, what is the risk of deletion/inversion/translocation
    - **Safety Grade**: A through F, with detailed breakdown

> **How to interpret:**
> - If bridgeability ratio < 1.0: donor CANNOT bridge the gap
> - Translocation risk: anything > 1% is concerning for clinical applications
> - Safety grade accounts for viability, rearrangement risk, p53 selection

### STEP 7: Donor Quality Check (TopoPred)

26. Click **"Donor Quality Check (TopoPred)"** in the sidebar
27. Paste a DNA sequence in the text area (this would be your cssDNA donor sequence)
   - For demo: you can paste any DNA sequence, e.g., 500 random bases
28. Set the **Left arm region** (start, end) — e.g., 0 to 200
29. Set the **Right arm region** (start, end) — e.g., 300 to 500
30. Click **"Analyze Donor"**
31. You will see:
    - **G-quadruplex motifs**: any G-rich sequences that could form stable G4 structures
    - **Hairpins**: self-complementary regions that fold into stem-loops
    - **Accessibility scores**: what fraction of each homology arm is available for RAD51 binding
    - **Recommendation**: "Good donor design" or "Redesign recommended"

> **How to interpret:**
> - Accessibility > 60% = good (arms are mostly single-stranded and accessible)
> - Accessibility < 40% = problematic (arms are sequestered in secondary structures)
> - G4 motifs in homology arms = critical problem

### STEP 8: Full Analysis Report

32. Click **"Full Analysis Report"** in the sidebar
33. Click **"Generate Full Report"**
34. This combines results from ALL modules you have run
35. Click **"Download Report"** to save as a text file

---

## 3. Step-by-Step Walkthrough (Real Gene)

To use a real human gene (e.g., NF1, DMD, BRCA2), you use the command-line script which fetches real exon coordinates from the Ensembl database (no genome download needed).

### Option A: Use the pre-built script

```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA
python crisprarchitect/examples/fetch_real_genes.py
```

This automatically:
1. Connects to the Ensembl REST API (internet required)
2. Fetches real GRCh38 exon coordinates for NF1, DMD, and BRCA2
3. Runs the full CRISPRArchitect analysis on each
4. Prints a comparison table

### Option B: Analyze any gene by name

```python
import json, urllib.request

# Fetch gene info from Ensembl
gene_symbol = "NF1"  # Change this to your gene
url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json;expand=1"
req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
with urllib.request.urlopen(req) as response:
    gene_info = json.loads(response.read().decode())

# Find canonical transcript and its exons
for t in gene_info['Transcript']:
    if t.get('is_canonical') == 1:
        exons = sorted(t['Exon'], key=lambda e: e['start'])
        print(f"Gene: {gene_symbol}, Exons: {len(exons)}")

        # Convert to CRISPRArchitect format
        exon_list = []
        for i, e in enumerate(exons, 1):
            exon_list.append({
                'number': i,
                'start': e['start'],
                'end': e['end']
            })
        break

# Now use CRISPRArchitect
from crisprarchitect.mosaic.gene_structure import GeneStructure
gene = GeneStructure.from_manual(gene_symbol, exon_list)
print(f"Genomic distance exon 20 to 50: {gene.genomic_distance(20, 50)} bp")
```

### Option C: Enter real coordinates in the web app

1. Go to **"Gene & Mutation Setup"** page
2. Select **"Enter exon coordinates manually"**
3. Paste exon coordinates in the text area, one per line:
   ```
   1,31094977,31095369
   2,31155983,31156126
   3,31159010,31159093
   ...
   ```
   (Format: exon_number, start, end)
4. Or select **"Upload CSV"** and upload a CSV file with columns: number, start, end

---

## 4. Where is the Human Genome Used?

### Short Answer

**You do NOT need to download the human genome.** We fetch only the gene annotation data (exon coordinates) from the Ensembl REST API over the internet.

### Detailed Explanation

Here is exactly what data is used and where it comes from:

### What We USE:

| Data | Source | How Accessed | Needed Locally? |
|------|--------|-------------|----------------|
| **Exon coordinates** (start/end positions for each exon) | Ensembl GRCh38 annotation | REST API call over internet | NO — fetched on-the-fly |
| **Gene location** (chromosome, start, end) | Ensembl GRCh38 annotation | REST API call | NO |
| **Intron sizes** (calculated as gaps between exons) | Derived from exon coordinates | Calculated by our code | NO |

### What We Do NOT Use:

| Data | Why Not |
|------|---------|
| **Actual DNA sequences** of introns/exons | Our models work with distances and sizes, not sequence content |
| **Whole genome FASTA** (GRCh38.fa, ~3 GB) | Not needed for distance calculations or strategy optimization |
| **Hi-C data files** (for ChromBridge) | We use a polymer physics MODEL instead, calibrated to published Hi-C statistics |
| **ChIP-seq data** (for LoopSim CTCF sites) | We use estimated CTCF spacing; real data would improve accuracy |

### The Ensembl REST API

When you run `fetch_real_genes.py`, it makes HTTP requests like:

```
https://rest.ensembl.org/lookup/symbol/homo_sapiens/NF1?content-type=application/json
```

This returns a JSON response with the gene's exon coordinates. No genome download needed. The API is free and requires no authentication.

### For the Web App

The web app currently works with:
1. **Demo genes** (computer-generated with realistic exon/intron sizes)
2. **Manually entered coordinates** (you paste real coordinates from Ensembl/UCSC)
3. **CSV upload** (you provide a CSV of exon coordinates)

The command-line script (`fetch_real_genes.py`) adds the ability to auto-fetch from Ensembl.

### When Would You Need the Actual Genome?

You would need to download genome sequences only if you wanted to:
- Design actual sgRNA sequences (need to know the target DNA sequence)
- Analyze actual donor template secondary structure (need the real sequence for TopoPred)
- Check for off-target sites (need genome-wide sequence search)

For these tasks, you would use:
```bash
# Download just your gene region (not the whole genome):
# Using samtools/tabix to fetch a specific region from UCSC:
# Or use Ensembl REST API to get the sequence:
curl "https://rest.ensembl.org/sequence/region/human/17:31094927..31377687:1?content-type=text/plain"
```

---

## 5. Command-Line Usage (No Web App)

You can use CRISPRArchitect entirely from Python scripts without the web interface:

### Quick analysis

```python
import sys
sys.path.insert(0, '/Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA')

# 1. ConversionSim
from crisprarchitect.conversion_sim import ConversionSimulator

sim = ConversionSimulator(
    cut_type='staggered_5prime',  # enFnCas9
    overhang_length=3,
    donor_topology='circular_ssDNA',
    homology_arm_length=300,
    cell_type='iPSC',
    n_simulations=10000
)
results = sim.run()
sim.summary()  # Prints full report
sim.probability_at_distance(500)  # Returns float

# 2. ChromBridge
from crisprarchitect.chrombridge import ChromatinDistancePredictor

pred = ChromatinDistancePredictor()
dist = pred.predict_3d_distance(500000)  # 500 kb
bridge = pred.can_donor_bridge(500000, 3000, 'circular_ssDNA')
print(bridge.explanation)

# 3. MOSAIC
from crisprarchitect.mosaic.gene_structure import GeneStructure
from crisprarchitect.mosaic.mutation_classifier import Mutation
from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
from crisprarchitect.mosaic.scorer import StrategyScorer

exons = [{'number': i, 'start': i*5000, 'end': i*5000+150} for i in range(1, 61)]
gene = GeneStructure.from_manual('MyGene', exons)
m1 = Mutation(exon_number=20, position=100000, ref_allele='G', alt_allele='A')
m2 = Mutation(exon_number=50, position=250000, ref_allele='C', alt_allele='T')

strategies = StrategyEnumerator().enumerate_strategies(gene, [m1, m2], 'iPSC', 'enFnCas9')
ranked = StrategyScorer().rank_strategies(strategies, 'iPSC')

for r in ranked:
    print(f"#{r.rank} {r.strategy.name}: score={r.overall_score:.3f}")

# 4. TopoPred
from crisprarchitect.topopred import DonorAnalyzer

analyzer = DonorAnalyzer()
report = analyzer.analyze(
    sequence="ATCGATCG" * 100,  # Your donor sequence
    left_arm=(0, 200),
    right_arm=(600, 800)
)
print(report['summary'])

# 5. LoopSim
from crisprarchitect.loopsim import LoopSimulator

sim = LoopSimulator(0, 2000000, resolution_bp=1000)
sim.setup_fiber(ctcf_spacing=100000)
sim.add_dsb(500000)
sim.add_second_dsb(1500000)
results = sim.run(n_simulations=500)
sim.summary()
```

### Run the full scenario analysis

```bash
cd /Users/vishalbharti/Downloads/DSB_REPAIR_MECHANICS_LITERATURE_REVIEW_cssDNA
python crisprarchitect/examples/your_scenario_analysis.py
```

### Run real gene analysis

```bash
python crisprarchitect/examples/fetch_real_genes.py
```

---

## 6. Troubleshooting

### "ModuleNotFoundError: No module named 'crisprarchitect'"

Make sure you run streamlit from the `crisprarchitect/` directory:
```bash
cd /path/to/crisprarchitect
streamlit run webapp/app.py
```

### "TypeError: unsupported operand type(s) for |"

Your Python version is 3.9. This was fixed — make sure you have the latest code.

### Web app shows error on a page

Click a different page and come back. Or restart:
```bash
pkill -f streamlit
streamlit run webapp/app.py
```

### Ensembl API timeout

The Ensembl REST API occasionally has slow responses. Just re-run the script.

### "No exons found for transcript"

The gene might not have a canonical transcript in Ensembl. Try a different gene symbol.
