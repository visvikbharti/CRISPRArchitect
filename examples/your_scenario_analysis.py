#!/usr/bin/env python3
"""
=============================================================================
CRISPRArchitect — Complete Analysis of Your Experimental Scenario
=============================================================================

SCENARIO:
    - Gene with ~60 exons (median exon size ~150 nt)
    - Mutations in exon 40 and exon 60
    - Proposed approach: single cssDNA (~3 kb exonic) + enFnCas9 + two sgRNAs
    - Cell type: human iPSCs

QUESTION: Can a single cssDNA correct both mutations simultaneously?

This script runs ALL four CRISPRArchitect modules to answer this question
with quantitative evidence.

HOW TO RUN:
    cd crisprarchitect/
    python examples/your_scenario_analysis.py
=============================================================================
"""

import sys
import os
import random

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# ============================================================
# STEP 0: Set up the scenario
# ============================================================
print("""
╔══════════════════════════════════════════════════════════════╗
║       CRISPRArchitect — Full Scenario Analysis              ║
║       Gene: 60-exon gene with mutations in exon 40 & 60    ║
║       System: enFnCas9 + cssDNA in iPSCs                   ║
╚══════════════════════════════════════════════════════════════╝
""")

# ============================================================
# STEP 1: Define the gene structure
# ============================================================
print("=" * 70)
print("STEP 1: GENE STRUCTURE ANALYSIS")
print("=" * 70)

from crisprarchitect.mosaic.gene_structure import GeneStructure

# Create a realistic 60-exon gene
# Exon sizes: ~100-250 bp (median ~150, matching your description)
# Intron sizes: variable (some small, some huge — realistic for human genes)
random.seed(42)  # For reproducibility
exons = []
pos = 0
for i in range(1, 61):
    exon_size = random.randint(100, 250)  # Realistic exon sizes
    exons.append({'number': i, 'start': pos, 'end': pos + exon_size})
    pos += exon_size
    if i < 60:
        # Human introns are highly variable in size
        # Some are tiny (~1 kb), some are massive (>100 kb)
        intron_size = random.choice([
            1500, 2000, 3000, 5000, 8000,     # Small introns
            15000, 25000, 50000,               # Medium introns
            80000, 100000, 150000, 200000,     # Large introns
        ])
        pos += intron_size

gene = GeneStructure.from_manual('MyLargeGene', exons)

print(f"\nGene: {gene.gene_name}")
print(f"Total exons: {len(gene.exons)}")
print(f"Total gene span: {gene.genomic_distance(1, 60) / 1e6:.2f} Mb")
print(f"\n--- Between Exon 40 and Exon 60 ---")
genomic_dist = gene.genomic_distance(40, 60)
exonic_dist = gene.exonic_distance(40, 60)
intronic_dist = gene.intronic_distance(40, 60)
print(f"Genomic distance (including introns): {genomic_dist / 1e6:.2f} Mb")
print(f"Exonic distance (just exon bases):    {exonic_dist / 1e3:.1f} kb")
print(f"Intronic distance:                    {intronic_dist / 1e6:.2f} Mb")
print(f"Can span with single 3kb cssDNA?      {gene.can_span_with_single_donor(40, 60, 3000)}")
print(f"Can span with single 10kb cssDNA?     {gene.can_span_with_single_donor(40, 60, 10000)}")
print(f"Can span with single 20kb cssDNA?     {gene.can_span_with_single_donor(40, 60, 20000)}")

# ============================================================
# STEP 2: ChromBridge — Can a cssDNA physically bridge the gap?
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: 3D CHROMATIN DISTANCE ANALYSIS (ChromBridge)")
print("=" * 70)

from crisprarchitect.chrombridge import ChromatinDistancePredictor

predictor = ChromatinDistancePredictor()

# Predict 3D distance between exon 40 and exon 60
dist_result = predictor.predict_3d_distance(genomic_dist)
print(f"\nGenomic distance: {genomic_dist / 1e6:.2f} Mb")
print(f"Predicted 3D distance: {dist_result.mean_3d_distance_nm:.0f} nm "
      f"({dist_result.mean_3d_distance_nm / 1000:.2f} μm)")
print(f"  5th percentile:  {dist_result.p5_nm:.0f} nm (closest 5% of cell configurations)")
print(f"  95th percentile: {dist_result.p95_nm:.0f} nm (most spread configurations)")

# Can a 3 kb cssDNA bridge this gap?
print("\n--- Bridgeability Analysis ---")
bridge = predictor.can_donor_bridge(
    genomic_distance_bp=genomic_dist,
    donor_size_bp=3000,
    donor_type='circular_ssDNA'
)
print(f"\n3 kb cssDNA donor:")
print(f"  Donor coil diameter: {bridge.donor_coil_diameter_nm:.0f} nm")
print(f"  Inter-locus distance: {bridge.inter_locus_distance_nm:.0f} nm")
print(f"  Bridgeability ratio: {bridge.bridgeability_ratio:.4f}")
print(f"  Can bridge? {bridge.feasible}")

# Also try with the maximum GATALYST cssDNA (20 kb)
bridge_20kb = predictor.can_donor_bridge(
    genomic_distance_bp=genomic_dist,
    donor_size_bp=20000,
    donor_type='circular_ssDNA'
)
print(f"\n20 kb cssDNA donor (maximum GATALYST size):")
print(f"  Donor coil diameter: {bridge_20kb.donor_coil_diameter_nm:.0f} nm")
print(f"  Can bridge? {bridge_20kb.feasible}")
print(f"  Ratio: {bridge_20kb.bridgeability_ratio:.4f}")

# Translocation risk
print("\n--- Translocation Risk ---")
from crisprarchitect.chrombridge.translocation import TranslocationRiskPredictor

trans_pred = TranslocationRiskPredictor()
risk = trans_pred.estimate_risk(
    site1_coord=gene.exons[39].start,
    site2_coord=gene.exons[59].start,
    same_chromosome=True
)
print(f"\nDeletion frequency estimate: {risk.deletion_frequency:.4f} ({risk.deletion_frequency*100:.2f}%)")
print(f"Inversion frequency estimate: {risk.inversion_frequency:.4f}")
print(f"Total rearrangement risk: {risk.total_rearrangement_risk:.4f}")
print(f"Risk level: {risk.risk_level}")

# ============================================================
# STEP 3: ConversionSim — How far can HDR reach from a cut?
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: GENE CONVERSION TRACT SIMULATION (ConversionSim)")
print("=" * 70)

from crisprarchitect.conversion_sim import ConversionSimulator

# Simulate with enFnCas9 + cssDNA in iPSCs
sim = ConversionSimulator(
    cut_type='staggered_5prime',
    overhang_length=3,          # enFnCas9 estimated stagger
    donor_topology='circular_ssDNA',
    homology_arm_length=300,    # Optimal for cssDNA (Iyer et al., 2022)
    cell_type='iPSC',
    n_simulations=50000         # More simulations for better statistics
)

results = sim.run()
stats = sim.summary()

# Can the tract reach from exon 40 to exon 60?
print(f"\n--- Can HDR from exon 40 reach exon 60? ---")
print(f"Distance to cover: {genomic_dist / 1e6:.2f} Mb = {genomic_dist:,} bp")
prob_reach = sim.probability_at_distance(genomic_dist)
print(f"Probability of tract reaching exon 60: {prob_reach:.10f}")
print(f"That's effectively ZERO.")
print(f"\nFor comparison:")
print(f"  P(reaching 500 bp):   {sim.probability_at_distance(500):.1%}")
print(f"  P(reaching 1000 bp):  {sim.probability_at_distance(1000):.1%}")
print(f"  P(reaching 2000 bp):  {sim.probability_at_distance(2000):.1%}")
print(f"  P(reaching 5000 bp):  {sim.probability_at_distance(5000):.1%}")
print(f"  P(reaching {genomic_dist:,} bp): {prob_reach}")

# ============================================================
# STEP 4: MOSAIC — What IS the best strategy?
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: STRATEGY OPTIMIZATION (MOSAIC)")
print("=" * 70)

from crisprarchitect.mosaic.mutation_classifier import Mutation, MutationClassifier
from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
from crisprarchitect.mosaic.scorer import StrategyScorer

# Define the two mutations
# Example: both are transition mutations (G>A and C>T)
mutation1 = Mutation(
    exon_number=40,
    position=gene.exons[39].start + 75,  # Middle of exon 40
    ref_allele='G',
    alt_allele='A',
)
mutation2 = Mutation(
    exon_number=60,
    position=gene.exons[59].start + 75,  # Middle of exon 60
    ref_allele='C',
    alt_allele='T',
)

# Classify mutations
classifier = MutationClassifier()
print(f"\nMutation 1 (Exon 40): {mutation1.ref_allele}>{mutation1.alt_allele}")
print(f"  Type: {classifier.classify(mutation1.ref_allele, mutation1.alt_allele)}")
print(f"  Base editing amenable: {classifier.base_editing_amenable(mutation1)}")
print(f"  Prime editing amenable: {classifier.prime_editing_amenable(mutation1)}")

print(f"\nMutation 2 (Exon 60): {mutation2.ref_allele}>{mutation2.alt_allele}")
print(f"  Type: {classifier.classify(mutation2.ref_allele, mutation2.alt_allele)}")
print(f"  Base editing amenable: {classifier.base_editing_amenable(mutation2)}")
print(f"  Prime editing amenable: {classifier.prime_editing_amenable(mutation2)}")

# Enumerate all strategies
enumerator = StrategyEnumerator()
strategies = enumerator.enumerate_strategies(gene, [mutation1, mutation2], 'iPSC', 'enFnCas9')

print(f"\n--- All Feasible Strategies ---")
for s in strategies:
    print(f"\n  Strategy: {s.name}")
    print(f"  DSBs needed: {s.num_dsbs}")
    print(f"  Editing rounds: {s.num_editing_rounds}")
    print(f"  Donors needed: {s.donor_templates_needed}")
    print(f"  Estimated efficiency: {s.estimated_efficiency:.4f} ({s.estimated_efficiency*100:.2f}%)")
    print(f"  Safety score: {s.safety_score:.2f}")
    print(f"  Clones to screen: ~{s.screening_burden}")

# Score and rank
scorer = StrategyScorer()
ranked = scorer.rank_strategies(strategies, 'iPSC')

print(f"\n{'='*70}")
print("FINAL RECOMMENDATIONS (ranked by overall score)")
print(f"{'='*70}")

for r in ranked:
    print(f"\n{'─'*60}")
    print(f"  Rank #{r.rank}: {r.strategy.name}")
    print(f"{'─'*60}")
    print(f"  Overall Score:    {r.overall_score:.3f}")
    print(f"  Efficiency Score: {r.efficiency_score:.3f}")
    print(f"  Safety Score:     {r.safety_score:.3f}")
    print(f"  Time Score:       {r.time_score:.3f}")
    print(f"  Description: {r.strategy.description}")

# ============================================================
# STEP 5: TopoPred — If we design a cssDNA donor, is it good?
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: DONOR TEMPLATE QUALITY CHECK (cssDNA-TopoPred)")
print("=" * 70)

from crisprarchitect.topopred import DonorAnalyzer

# Create a hypothetical cssDNA donor for exon 40 correction
# 300 bp left arm (intronic) + 150 bp exon 40 (corrected) + 300 bp right arm (intronic)
# Using a random but realistic sequence for demo
random.seed(123)
def random_dna(length):
    return ''.join(random.choice('ATCG') for _ in range(length))

left_arm = random_dna(300)
exon_40_corrected = random_dna(150)
right_arm = random_dna(300)
# Add a phagemid backbone (~2200 bp)
backbone = random_dna(2200)

donor_sequence = left_arm + exon_40_corrected + right_arm + backbone

print(f"\nHypothetical cssDNA donor for exon 40 correction:")
print(f"  Left homology arm:  300 bp (intronic sequence)")
print(f"  Corrected exon 40:  150 bp")
print(f"  Right homology arm: 300 bp (intronic sequence)")
print(f"  Phagemid backbone:  2200 bp")
print(f"  Total cssDNA size:  {len(donor_sequence)} nt")

analyzer = DonorAnalyzer()
report = analyzer.analyze(
    sequence=donor_sequence,
    left_arm=(0, 300),
    right_arm=(450, 750),
    coding_regions=[(300, 450)]  # The exon region
)

print(f"\n{report['summary']}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("COMPREHENSIVE ANSWER TO YOUR QUESTION")
print("=" * 70)

print("""
QUESTION: Can a single cssDNA (~3 kb exonic) serve as HDR template for
          both exon 40 and exon 60 cuts with enFnCas9 in iPSCs?

ANSWER: NO — and here is the quantitative evidence from all four analyses:

1. GENOMIC DISTANCE (Gene Structure Analysis):
""")
print(f"   Exons 40-60 are separated by {genomic_dist/1e6:.2f} Mb of genomic sequence")
print(f"   A 3 kb cssDNA cannot span {genomic_dist/1e6:.2f} Mb (it's {genomic_dist/3000:.0f}x too short)")
print(f"   Even a 20 kb cssDNA (max GATALYST size) cannot span this gap")

print("""
2. PHYSICAL DISTANCE (ChromBridge):
""")
print(f"   In 3D nuclear space, the two loci are ~{dist_result.mean_3d_distance_nm:.0f} nm apart")
print(f"   A 3 kb cssDNA has a coil diameter of ~{bridge.donor_coil_diameter_nm:.0f} nm")
print(f"   The donor is {dist_result.mean_3d_distance_nm/bridge.donor_coil_diameter_nm:.0f}x too small to bridge the gap")

print("""
3. GENE CONVERSION TRACT (ConversionSim):
""")
print(f"   Median HDR tract length: ~{stats['tract_median_bp']:.0f} bp")
print(f"   Even the 95th percentile: ~{stats['tract_p95_bp']:.0f} bp")
print(f"   Distance to exon 60: {genomic_dist:,} bp — unreachable by HDR")

print("""
4. BEST STRATEGY (MOSAIC):
""")
best = ranked[0]
print(f"   Recommended: {best.strategy.name}")
print(f"   {best.strategy.description}")

print("""
================================================================================
BOTTOM LINE: Use two separate cssDNA donors with two separate sgRNAs.
For your specific mutations (both transitions), DUAL BASE EDITING is optimal
— no DSBs needed, highest safety in iPSCs.

If HDR is required: use SEQUENTIAL editing (correct one site, validate,
then correct the other) for maximum safety in iPSCs.
================================================================================
""")
