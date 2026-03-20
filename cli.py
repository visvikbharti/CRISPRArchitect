#!/usr/bin/env python3
"""
CRISPRArchitect — Command Line Interface
==========================================

Usage:
    crisprarchitect analyze --gene NF1 --exon1 20 --mut1 G>A --exon2 50 --mut2 C>T
    crisprarchitect fetch --gene DMD
    crisprarchitect simulate --cell iPSC --nuclease enFnCas9 --donor cssDNA
    crisprarchitect webapp
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))


def cmd_analyze(args):
    """Run full CRISPRArchitect analysis on a gene."""
    from crisprarchitect.utils.ensembl import fetch_gene
    from crisprarchitect.mosaic.gene_structure import GeneStructure
    from crisprarchitect.mosaic.mutation_classifier import Mutation, MutationClassifier
    from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
    from crisprarchitect.mosaic.scorer import StrategyScorer
    from crisprarchitect.chrombridge import ChromatinDistancePredictor
    from crisprarchitect.conversion_sim import ConversionSimulator

    print(f"Fetching {args.gene} from Ensembl...")
    gene_info = fetch_gene(args.gene)
    print(gene_info.summary())
    print()

    gene = GeneStructure.from_manual(args.gene, gene_info.exons)

    # Parse mutations
    ref1, alt1 = args.mut1.split(">")
    ref2, alt2 = args.mut2.split(">")
    m1 = Mutation(exon_number=args.exon1,
                  position=gene_info.exons[args.exon1 - 1]["start"] + 50,
                  ref_allele=ref1, alt_allele=alt1)
    m2 = Mutation(exon_number=args.exon2,
                  position=gene_info.exons[args.exon2 - 1]["start"] + 50,
                  ref_allele=ref2, alt_allele=alt2)

    # Distance analysis
    dist = gene.genomic_distance(args.exon1, args.exon2)
    print(f"Genomic distance exon {args.exon1} to {args.exon2}: {dist:,} bp ({dist/1e6:.3f} Mb)")

    # ChromBridge
    pred = ChromatinDistancePredictor()
    d3d = pred.predict_3d_distance(abs(dist))
    bridge = pred.can_donor_bridge(abs(dist), 3000, "circular_ssDNA")
    print(f"3D distance: {d3d.mean_3d_distance_nm:.0f} nm")
    print(f"3 kb cssDNA can bridge: {bridge.feasible}")
    print()

    # ConversionSim
    sim = ConversionSimulator(
        cut_type="staggered_5prime", overhang_length=3,
        donor_topology="circular_ssDNA", homology_arm_length=300,
        cell_type=args.cell, n_simulations=10000
    )
    sim.run()
    sim.summary()
    print()

    # MOSAIC
    classifier = MutationClassifier()
    print(f"Mut1 (exon {args.exon1}, {ref1}>{alt1}): "
          f"{classifier.classify(ref1, alt1)}, "
          f"base-editable={classifier.base_editing_amenable(m1)}")
    print(f"Mut2 (exon {args.exon2}, {ref2}>{alt2}): "
          f"{classifier.classify(ref2, alt2)}, "
          f"base-editable={classifier.base_editing_amenable(m2)}")

    strategies = StrategyEnumerator().enumerate_strategies(
        gene, [m1, m2], args.cell, "SpCas9"
    )
    ranked = StrategyScorer().rank_strategies(strategies, args.cell)

    print(f"\nMOSAIC Strategy Ranking:")
    for r in ranked:
        print(f"  #{r.rank} {r.strategy.name} (score={r.overall_score:.3f})")
    print(f"\nRecommended: {ranked[0].strategy.name}")
    print(f"  {ranked[0].strategy.description}")


def cmd_fetch(args):
    """Fetch gene info from Ensembl."""
    from crisprarchitect.utils.ensembl import fetch_gene
    gene_info = fetch_gene(args.gene)
    print(gene_info.summary())
    print(f"\nExon coordinates:")
    for e in gene_info.exons:
        size = e["end"] - e["start"]
        print(f"  Exon {e['number']:3d}: {e['start']:>12,} - {e['end']:>12,}  ({size:>5,} bp)")


def cmd_simulate(args):
    """Run ConversionSim."""
    from crisprarchitect.conversion_sim import ConversionSimulator
    sim = ConversionSimulator(
        cut_type="staggered_5prime" if args.nuclease == "enFnCas9" else "blunt",
        overhang_length=3 if args.nuclease == "enFnCas9" else 0,
        donor_topology="circular_ssDNA" if args.donor == "cssDNA" else "linear_ssDNA",
        homology_arm_length=args.arms,
        cell_type=args.cell,
        n_simulations=args.n
    )
    sim.run()
    sim.summary()


def cmd_webapp(args):
    """Launch the Streamlit web app."""
    import subprocess
    app_path = os.path.join(os.path.dirname(__file__), "..", "webapp", "app.py")
    if not os.path.exists(app_path):
        app_path = os.path.join(os.path.dirname(__file__), "webapp", "app.py")
    subprocess.run(["streamlit", "run", app_path, "--server.port", str(args.port)])


def main():
    parser = argparse.ArgumentParser(
        prog="crisprarchitect",
        description="CRISPRArchitect — Multi-Site Genome Editing Strategy Optimizer"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # analyze
    p_analyze = subparsers.add_parser("analyze", help="Full analysis of a gene with two mutations")
    p_analyze.add_argument("--gene", required=True, help="Gene symbol (e.g., NF1)")
    p_analyze.add_argument("--exon1", type=int, required=True, help="First mutation exon number")
    p_analyze.add_argument("--mut1", required=True, help="First mutation (e.g., G>A)")
    p_analyze.add_argument("--exon2", type=int, required=True, help="Second mutation exon number")
    p_analyze.add_argument("--mut2", required=True, help="Second mutation (e.g., C>T)")
    p_analyze.add_argument("--cell", default="iPSC", help="Cell type (default: iPSC)")
    p_analyze.set_defaults(func=cmd_analyze)

    # fetch
    p_fetch = subparsers.add_parser("fetch", help="Fetch gene info from Ensembl")
    p_fetch.add_argument("--gene", required=True, help="Gene symbol")
    p_fetch.set_defaults(func=cmd_fetch)

    # simulate
    p_sim = subparsers.add_parser("simulate", help="Run ConversionSim tract simulation")
    p_sim.add_argument("--cell", default="iPSC", help="Cell type")
    p_sim.add_argument("--nuclease", default="enFnCas9", help="Nuclease")
    p_sim.add_argument("--donor", default="cssDNA", help="Donor topology")
    p_sim.add_argument("--arms", type=int, default=300, help="Homology arm length (bp)")
    p_sim.add_argument("-n", type=int, default=10000, help="Number of simulations")
    p_sim.set_defaults(func=cmd_simulate)

    # webapp
    p_web = subparsers.add_parser("webapp", help="Launch the Streamlit web interface")
    p_web.add_argument("--port", type=int, default=8501, help="Port number")
    p_web.set_defaults(func=cmd_webapp)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
