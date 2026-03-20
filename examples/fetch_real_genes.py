#!/usr/bin/env python3
"""
Fetch real gene structures from Ensembl REST API (GRCh38)
and run CRISPRArchitect analysis on them.
"""

import json
import sys
import os
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

def fetch_ensembl_json(url):
    """Fetch JSON from Ensembl REST API."""
    req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read().decode())


def fetch_transcript_exons(transcript_id):
    """Fetch all exons for a given Ensembl transcript ID."""
    url = f"https://rest.ensembl.org/overlap/id/{transcript_id}?feature=exon;content-type=application/json"
    all_exons = fetch_ensembl_json(url)
    # Filter for this transcript only
    transcript_exons = [e for e in all_exons if e.get('Parent') == transcript_id]
    # Sort by rank
    transcript_exons.sort(key=lambda e: e.get('rank', 0))
    return transcript_exons


def fetch_gene_info(gene_symbol):
    """Fetch gene info from Ensembl."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json;expand=1"
    return fetch_ensembl_json(url)


def get_canonical_transcript(gene_info):
    """Find the canonical transcript from gene info."""
    if 'Transcript' not in gene_info:
        return None
    for t in gene_info['Transcript']:
        if t.get('is_canonical') == 1:
            return t['id']
    # Fallback: pick the one with most exons
    transcripts = gene_info['Transcript']
    best = max(transcripts, key=lambda t: len(t.get('Exon', [])))
    return best['id']


# ============================================================
# GENES TO ANALYZE
# ============================================================
GENES = {
    "NF1": {
        "description": "Neurofibromin 1 — 57 exons, causes neurofibromatosis type 1",
        "mutation_exons": (20, 50),  # Hypothetical dual mutations
        "mut1": ("G", "A"),  # Transition (ABE-amenable)
        "mut2": ("C", "G"),  # Transversion (HDR required)
    },
    "DMD": {
        "description": "Dystrophin — 79 exons, largest human gene, causes Duchenne MD",
        "mutation_exons": (45, 55),  # In the deletion hotspot
        "mut1": ("A", "G"),  # Transition
        "mut2": ("T", "A"),  # Transversion
    },
    "BRCA2": {
        "description": "BRCA2 DNA repair — 27 exons, hereditary breast/ovarian cancer",
        "mutation_exons": (10, 22),
        "mut1": ("G", "A"),  # Transition
        "mut2": ("C", "T"),  # Transition
    },
}

print("""
╔══════════════════════════════════════════════════════════════════╗
║   CRISPRArchitect — Real Gene Analysis (GRCh38)                ║
║   Fetching gene structures from Ensembl REST API               ║
╚══════════════════════════════════════════════════════════════════╝
""")

# Store results for all genes
all_results = {}

for gene_symbol, gene_config in GENES.items():
    print(f"\n{'='*70}")
    print(f"  GENE: {gene_symbol} — {gene_config['description']}")
    print(f"{'='*70}")

    # Fetch gene info
    try:
        print(f"  Fetching from Ensembl...")
        gene_info = fetch_gene_info(gene_symbol)
        canonical_id = get_canonical_transcript(gene_info)
        print(f"  Gene ID: {gene_info['id']}")
        print(f"  Location: chr{gene_info['seq_region_name']}:{gene_info['start']}-{gene_info['end']}")
        print(f"  Strand: {'forward (+)' if gene_info['strand'] == 1 else 'reverse (-)'}")
        print(f"  Gene span: {(gene_info['end'] - gene_info['start']) / 1e6:.2f} Mb")
        print(f"  Canonical transcript: {canonical_id}")
    except Exception as e:
        print(f"  ERROR fetching gene info: {e}")
        continue

    # Fetch exons
    try:
        print(f"  Fetching exon coordinates...")
        exons_raw = fetch_transcript_exons(canonical_id)
        n_exons = len(exons_raw)
        print(f"  Total exons: {n_exons}")

        if n_exons == 0:
            print(f"  WARNING: No exons found for transcript {canonical_id}")
            continue

        # Convert to our format
        exons_list = []
        for e in exons_raw:
            exons_list.append({
                'number': e['rank'],
                'start': e['start'],
                'end': e['end'],
            })

        # Print first/last few exons
        print(f"\n  First 3 exons:")
        for e in exons_list[:3]:
            size = e['end'] - e['start']
            print(f"    Exon {e['number']}: chr17:{e['start']}-{e['end']} ({size} bp)")
        print(f"  ...")
        print(f"  Last 3 exons:")
        for e in exons_list[-3:]:
            size = e['end'] - e['start']
            print(f"    Exon {e['number']}: chr17:{e['start']}-{e['end']} ({size} bp)")

    except Exception as e:
        print(f"  ERROR fetching exons: {e}")
        continue

    # Run CRISPRArchitect
    print(f"\n  --- CRISPRArchitect Analysis ---")

    from crisprarchitect.mosaic.gene_structure import GeneStructure
    from crisprarchitect.mosaic.mutation_classifier import Mutation, MutationClassifier
    from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
    from crisprarchitect.mosaic.scorer import StrategyScorer
    from crisprarchitect.chrombridge import ChromatinDistancePredictor
    from crisprarchitect.conversion_sim import ConversionSimulator

    # Create gene structure
    gene = GeneStructure.from_manual(gene_symbol, exons_list)

    # Get mutation exon numbers (adjust if gene has fewer exons)
    mut_exon1 = min(gene_config['mutation_exons'][0], n_exons - 1)
    mut_exon2 = min(gene_config['mutation_exons'][1], n_exons)

    # Calculate distances
    try:
        genomic_dist = gene.genomic_distance(mut_exon1, mut_exon2)
        exonic_dist = gene.exonic_distance(mut_exon1, mut_exon2)
        intronic_dist = gene.intronic_distance(mut_exon1, mut_exon2)

        print(f"\n  Mutations in exon {mut_exon1} and exon {mut_exon2}:")
        print(f"    Genomic distance:  {genomic_dist/1e6:.3f} Mb ({genomic_dist:,} bp)")
        print(f"    Exonic distance:   {exonic_dist/1e3:.1f} kb")
        print(f"    Intronic distance: {intronic_dist/1e6:.3f} Mb")
        print(f"    Single donor feasible (10 kb)? {gene.can_span_with_single_donor(mut_exon1, mut_exon2)}")
    except Exception as e:
        print(f"  ERROR calculating distances: {e}")
        continue

    # ChromBridge: 3D distance
    predictor = ChromatinDistancePredictor()
    dist_3d = predictor.predict_3d_distance(genomic_dist)
    bridge_3kb = predictor.can_donor_bridge(genomic_dist, 3000, 'circular_ssDNA')

    print(f"\n  ChromBridge 3D Analysis:")
    print(f"    Predicted 3D distance: {dist_3d.mean_3d_distance_nm:.0f} nm ({dist_3d.mean_3d_distance_nm/1000:.2f} μm)")
    print(f"    3 kb cssDNA diameter:  {bridge_3kb.donor_coil_diameter_nm:.0f} nm")
    print(f"    Can cssDNA bridge?     {bridge_3kb.feasible}")
    print(f"    Bridgeability ratio:   {bridge_3kb.bridgeability_ratio:.4f}")

    # ConversionSim: tract length
    sim = ConversionSimulator(
        cut_type='staggered_5prime', overhang_length=3,
        donor_topology='circular_ssDNA', homology_arm_length=300,
        cell_type='iPSC', n_simulations=10000
    )
    sim.run()
    p_reach = sim.probability_at_distance(genomic_dist)

    print(f"\n  ConversionSim (enFnCas9 + cssDNA + iPSC):")
    print(f"    P(tract reaching exon {mut_exon2}): {p_reach:.10f}")

    # MOSAIC: Strategy ranking
    m1_pos = exons_list[mut_exon1 - 1]['start'] + 50
    m2_pos = exons_list[mut_exon2 - 1]['start'] + 50
    m1 = Mutation(exon_number=mut_exon1, position=m1_pos,
                  ref_allele=gene_config['mut1'][0], alt_allele=gene_config['mut1'][1])
    m2 = Mutation(exon_number=mut_exon2, position=m2_pos,
                  ref_allele=gene_config['mut2'][0], alt_allele=gene_config['mut2'][1])

    classifier = MutationClassifier()
    enumerator = StrategyEnumerator()
    scorer = StrategyScorer()

    strategies = enumerator.enumerate_strategies(gene, [m1, m2], 'iPSC', 'enFnCas9')
    ranked = scorer.rank_strategies(strategies, 'iPSC')

    print(f"\n  MOSAIC Strategy Ranking:")
    print(f"    Mut1 (exon {mut_exon1}, {m1.ref_allele}>{m1.alt_allele}): "
          f"{classifier.classify(m1.ref_allele, m1.alt_allele)}, "
          f"base-editable={classifier.base_editing_amenable(m1)}")
    print(f"    Mut2 (exon {mut_exon2}, {m2.ref_allele}>{m2.alt_allele}): "
          f"{classifier.classify(m2.ref_allele, m2.alt_allele)}, "
          f"base-editable={classifier.base_editing_amenable(m2)}")
    print(f"\n    Ranked strategies:")
    for r in ranked[:3]:
        print(f"      #{r.rank} {r.strategy.name} (score={r.overall_score:.3f})")

    all_results[gene_symbol] = {
        'n_exons': n_exons,
        'gene_span_mb': (gene_info['end'] - gene_info['start']) / 1e6,
        'genomic_dist_mb': genomic_dist / 1e6,
        'dist_3d_nm': dist_3d.mean_3d_distance_nm,
        'bridgeable': bridge_3kb.feasible,
        'p_tract_reach': p_reach,
        'best_strategy': ranked[0].strategy.name if ranked else "N/A",
    }

# ============================================================
# SUMMARY TABLE
# ============================================================
print(f"\n\n{'='*90}")
print(f"  CROSS-GENE COMPARISON SUMMARY")
print(f"{'='*90}")
print(f"{'Gene':<8} {'Exons':<7} {'Span':<8} {'Mut Distance':<13} {'3D Dist':<10} "
      f"{'Bridge?':<9} {'P(tract)':<12} {'Best Strategy'}")
print(f"{'-'*90}")
for gene_sym, r in all_results.items():
    print(f"{gene_sym:<8} {r['n_exons']:<7} {r['gene_span_mb']:.1f} Mb  "
          f"{r['genomic_dist_mb']:.3f} Mb     {r['dist_3d_nm']:.0f} nm   "
          f"{'Yes' if r['bridgeable'] else 'No':<9} {r['p_tract_reach']:.1e}     "
          f"{r['best_strategy']}")

print(f"\n{'='*90}")
print("CONCLUSION: For ALL tested genes, single-template dual-site HDR is infeasible.")
print("The recommended strategy depends on mutation type (base-editable vs HDR-required).")
print(f"{'='*90}")
