#!/usr/bin/env python3
"""
CRISPRArchitect v3 — Command Line Interface
=============================================

Usage:
    crisprarchitect analyze --gene NF1 --variants "c.910C>T" "c.4537C>T"
    crisprarchitect analyze-legacy --gene NF1 --exon1 20 --mut1 G>A --exon2 50 --mut2 C>T
    crisprarchitect fetch --gene DMD
    crisprarchitect simulate --cell iPSC --nuclease enFnCas9 --donor cssDNA
    crisprarchitect webapp
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))


def cmd_analyze(args):
    """Run the v3 pipeline: HGVS input -> TOPSIS-ranked strategies."""
    from core.pipeline.strategy_stage import StrategyPipeline
    from core.models import GenomicVariantInput

    # If HGVS strings provided, parse them and resolve to genomic coords
    if args.hgvs:
        from core.sequence.hgvs_parser import HGVSParser
        from core.sequence.fetcher import TranscriptFetcher
        parser = HGVSParser()
        fetcher = TranscriptFetcher()
        # Fetch transcript once by gene symbol (works for any HGVS format)
        try:
            transcript = fetcher.fetch_by_gene(args.gene)
            print(f"  Transcript: {transcript.transcript_id} ({transcript.gene_symbol})")
        except Exception as e:
            print(f"  ERROR fetching transcript for {args.gene}: {e}",
                  file=sys.stderr)
            sys.exit(1)
        variants = []
        for hgvs_str in args.hgvs:
            try:
                parsed = parser.parse(hgvs_str)
                # Map CDS position to genomic coordinate
                genomic_pos = parser._cds_to_genomic(
                    transcript, parsed.cds_position, parsed.intron_offset
                )
                v = GenomicVariantInput(
                    chromosome=transcript.chromosome,
                    position=genomic_pos,
                    ref_allele=parsed.ref_allele,
                    alt_allele=parsed.alt_allele,
                    gene_symbol=args.gene,
                    transcript_id=transcript.transcript_id,
                    name=parsed.raw_notation,
                )
                variants.append(v)
                print(f"  Parsed: {hgvs_str} -> chr{v.chromosome}:{v.position}")
            except ValueError as e:
                print(f"  ERROR parsing '{hgvs_str}': {e}", file=sys.stderr)
                sys.exit(1)
    elif args.chr and args.pos and args.ref and args.alt:
        variants = [GenomicVariantInput(
            chromosome=args.chr,
            position=args.pos,
            ref_allele=args.ref,
            alt_allele=args.alt,
            gene_symbol=args.gene,
            name=f"{args.gene}:{args.ref}{args.pos}{args.alt}",
        )]
    else:
        print("ERROR: Provide either --hgvs variants or --chr/--pos/--ref/--alt",
              file=sys.stderr)
        sys.exit(1)

    print(f"\nRunning CRISPRArchitect v3 pipeline for {args.gene}...")
    print(f"  Cell type: {args.cell}")
    print(f"  Primary nuclease: {args.nuclease}")
    print(f"  Variants: {len(variants)}\n")

    pipeline = StrategyPipeline(
        cell_type=args.cell,
        nuclease=args.nuclease,
    )
    result = pipeline.run(variants)

    # Print transcript info
    if result.transcript:
        ti = result.transcript
        print(f"Transcript: {ti.transcript_id} ({ti.gene_symbol})")
        print(f"  Exons: {len(ti.exons)}, Strand: {ti.strand}")
        print()

    # Print variant annotations
    for nv in result.variants:
        print(f"Variant: {nv.input.name}")
        if nv.coding:
            ca = nv.coding
            print(f"  Consequence: {ca.consequence.value}")
            if ca.hgvs_c:
                print(f"  HGVS c.: {ca.hgvs_c}")
            if ca.hgvs_p:
                print(f"  HGVS p.: {ca.hgvs_p}")
        if nv.ref_validation:
            rv = nv.ref_validation
            status = "PASS" if rv.is_valid else "MISMATCH"
            print(f"  Ref validation: {status}")
        print()

    # Print strategy ranking
    print(f"Strategy Ranking (TOPSIS + sensitivity analysis):")
    print(f"{'Rank':<6} {'Strategy':<35} {'Score':<8} {'Evidence':<10} {'Safety':<8}")
    print("-" * 70)
    for s in result.strategies:
        evidence = s.strategy.evidence_tier.value if s.strategy.evidence_tier else "N/A"
        print(f"  #{s.rank:<4} {s.strategy_name:<35} {s.overall_score:<8.3f} "
              f"{evidence:<10} {s.safety_score:<8.2f}")

    if result.strategies:
        top = result.strategies[0]
        print(f"\nRecommendation: {top.strategy_name}")
        print(f"  Evidence tier: {top.strategy.evidence_tier.value if top.strategy.evidence_tier else 'N/A'}")

    # Print rejected strategies
    if result.rejected_strategies:
        print(f"\nRejected strategies ({len(result.rejected_strategies)}):")
        for r in result.rejected_strategies:
            reasons = "; ".join(r.rejection_reasons) if r.rejection_reasons else "N/A"
            print(f"  {r.strategy_name}: {reasons}")


def cmd_analyze_legacy(args):
    """Run v1 analysis (exon-based input, legacy scoring)."""
    from crisprarchitect.utils.ensembl import fetch_gene
    from crisprarchitect.mosaic.gene_structure import GeneStructure
    from crisprarchitect.mosaic.mutation_classifier import Mutation, MutationClassifier
    from crisprarchitect.mosaic.strategy_enumerator import StrategyEnumerator
    from crisprarchitect.mosaic.scorer import StrategyScorer
    from crisprarchitect.conversion_sim import ConversionSimulator

    print(f"[Legacy mode] Fetching {args.gene} from Ensembl...")
    gene_info = fetch_gene(args.gene)
    print(gene_info.summary())
    print()

    gene = GeneStructure.from_manual(args.gene, gene_info.exons)

    ref1, alt1 = args.mut1.split(">")
    ref2, alt2 = args.mut2.split(">")
    m1 = Mutation(exon_number=args.exon1,
                  position=gene_info.exons[args.exon1 - 1]["start"] + 50,
                  ref_allele=ref1, alt_allele=alt1)
    m2 = Mutation(exon_number=args.exon2,
                  position=gene_info.exons[args.exon2 - 1]["start"] + 50,
                  ref_allele=ref2, alt_allele=alt2)

    dist = gene.genomic_distance(args.exon1, args.exon2)
    print(f"Genomic distance exon {args.exon1} to {args.exon2}: {dist:,} bp")

    print()

    sim = ConversionSimulator(
        cut_type="staggered_5prime", overhang_length=3,
        donor_topology="circular_ssDNA", homology_arm_length=300,
        cell_type=args.cell, n_simulations=10000
    )
    sim.run()
    sim.summary()
    print()

    classifier = MutationClassifier()
    strategies = StrategyEnumerator().enumerate_strategies(
        gene, [m1, m2], args.cell, "SpCas9"
    )
    ranked = StrategyScorer().rank_strategies(strategies, args.cell)

    print(f"MOSAIC Strategy Ranking (legacy v1 scorer):")
    for r in ranked:
        print(f"  #{r.rank} {r.strategy.name} (score={r.overall_score:.3f})")
    print(f"\nRecommended: {ranked[0].strategy.name}")


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
    app_path = os.path.join(os.path.dirname(__file__), "webapp", "app.py")
    if not os.path.exists(app_path):
        app_path = os.path.join(os.path.dirname(__file__), "..", "webapp", "app.py")
    subprocess.run(["streamlit", "run", app_path, "--server.port", str(args.port)])


def main():
    parser = argparse.ArgumentParser(
        prog="crisprarchitect",
        description="CRISPRArchitect v3 — Multi-Nuclease Genome Editing Strategy Optimizer"
    )
    parser.add_argument("--version", action="version", version="%(prog)s 3.0.0")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # analyze (v3 pipeline — primary entry point)
    p_analyze = subparsers.add_parser(
        "analyze",
        help="Analyze variants using v3 pipeline (TOPSIS + multi-nuclease)"
    )
    p_analyze.add_argument("--gene", required=True, help="Gene symbol (e.g., NF1)")
    p_analyze.add_argument("--hgvs", nargs="+",
                           help="HGVS variant(s) (e.g., c.910C>T c.4537C>T)")
    p_analyze.add_argument("--chr", help="Chromosome (e.g., 17)")
    p_analyze.add_argument("--pos", type=int, help="Genomic position (1-based)")
    p_analyze.add_argument("--ref", help="Reference allele")
    p_analyze.add_argument("--alt", help="Alternate allele")
    p_analyze.add_argument("--cell", default="iPSC",
                           help="Cell type (default: iPSC)")
    p_analyze.add_argument("--nuclease", default="SpCas9",
                           help="Primary nuclease (default: SpCas9)")
    p_analyze.set_defaults(func=cmd_analyze)

    # analyze-legacy (v1 pipeline — backward compatibility)
    p_legacy = subparsers.add_parser(
        "analyze-legacy",
        help="Legacy v1 analysis (exon-based input, weighted-sum scoring)"
    )
    p_legacy.add_argument("--gene", required=True, help="Gene symbol")
    p_legacy.add_argument("--exon1", type=int, required=True)
    p_legacy.add_argument("--mut1", required=True, help="e.g., G>A")
    p_legacy.add_argument("--exon2", type=int, required=True)
    p_legacy.add_argument("--mut2", required=True, help="e.g., C>T")
    p_legacy.add_argument("--cell", default="iPSC")
    p_legacy.set_defaults(func=cmd_analyze_legacy)

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
