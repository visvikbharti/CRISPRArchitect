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
    print(f"{'Rank':<6} {'Strategy':<35} {'Score':<8} {'Stability':<12} {'Evidence':<10} {'Safety':<8}")
    print("-" * 82)
    for s in result.strategies:
        evidence = s.strategy.evidence_tier.value if s.strategy.evidence_tier else "N/A"
        stability = f"{s.rank_stability:.1%}" if s.rank_stability is not None else "N/A"
        print(f"  #{s.rank:<4} {s.strategy_name:<35} {s.overall_score:<8.3f} "
              f"{stability:<12} {evidence:<10} {s.safety_score:<8.2f}")

    if result.strategies:
        from core.pipeline.strategy_stage import (
            RANK_STABILITY_ROBUST,
            RANK_STABILITY_STABLE,
            assess_stability,
        )

        top = result.strategies[0]
        assessment = assess_stability(result.strategies)

        print(f"\nRecommendation: {top.strategy_name}")
        print(f"  Evidence tier: {top.strategy.evidence_tier.value if top.strategy.evidence_tier else 'N/A'}")

        if top.rank_stability is not None:
            print(
                f"  Rank stability: {top.rank_stability:.1%} — {assessment.human_label} "
                f"(top-ranked in {top.rank_stability:.1%} of 10,000 weight permutations)"
            )
            if assessment.level == "robust":
                print(
                    f"    Above {RANK_STABILITY_ROBUST:.0%} threshold: unconditional recommendation."
                )
            elif assessment.level == "stable":
                print(
                    f"    In [{RANK_STABILITY_STABLE:.0%}, {RANK_STABILITY_ROBUST:.0%}): single "
                    f"recommendation but weight perturbations sometimes promote alternatives."
                )
            elif assessment.is_flip_sensitive and assessment.runner_up is not None:
                ru = assessment.runner_up
                gap = assessment.score_gap or 0.0
                print(
                    f"    Below {RANK_STABILITY_STABLE:.0%}: the choice between top and alternative "
                    f"depends on context. Consider the tradeoff below."
                )
                print(f"\n  Alternative considered: {ru.strategy_name}")
                if ru.rank_stability is not None:
                    print(f"    Alternative stability: {ru.rank_stability:.1%}")
                print(
                    f"    TOPSIS score gap: {gap:+.3f} "
                    f"(top = {top.overall_score:.3f}, alt = {ru.overall_score:.3f})"
                )
                print(f"\n  Why each might be preferred:")
                for line in assessment.preferential_reasoning:
                    print(f"    - {line}")
        else:
            print(f"  Rank stability: not computed (sensitivity analysis disabled).")

    # Print feasibility details (which specific editors were evaluated)
    if result.bundles:
        print(f"\nFeasibility Details:")
        for b in result.bundles:
            # Base editing details
            best_be = b.best_base_editing_result()
            if best_be and best_be.metadata:
                m = best_be.metadata
                editor = m.get('editor', '?')
                nuclease = m.get('nuclease', '?')
                window = m.get('window', '?')
                tier = m.get('evidence_tier', '?')
                print(f"  Best base editor: {editor} + {nuclease} "
                      f"(window {window}, tier {tier})")
                if best_be.best_guide:
                    g = best_be.best_guide
                    print(f"    Guide: {g.sequence_20mer} | PAM: {g.pam_sequence} "
                          f"| GC: {g.gc_content:.0%}")
                    print(f"    Target at window position {best_be.target_position_in_window}")
                if best_be.bystander_count > 0:
                    print(f"    Bystanders: {best_be.bystander_count} "
                          f"at window positions {best_be.bystander_positions}")

            # Count all tested editors
            n_tested = len(b.base_editing_results)
            n_feasible = sum(1 for r in b.base_editing_results
                           if r.label.value == 'feasible')
            n_marginal = sum(1 for r in b.base_editing_results
                           if r.label.value == 'marginal')
            n_rejected = sum(1 for r in b.base_editing_results
                           if r.label.value == 'not_feasible')
            if n_tested > 0:
                print(f"  Editor-nuclease combos tested: {n_tested} "
                      f"({n_feasible} feasible, {n_marginal} marginal, "
                      f"{n_rejected} not feasible)")

            # Prime editing
            pe = b.prime_editing_result
            if pe and pe.label.value != 'not_feasible':
                print(f"  Prime editing: {pe.label.value} (score {pe.score:.3f})")

            # HDR
            hdr = b.hdr_result
            if hdr and hdr.label.value != 'not_feasible':
                print(f"  HDR: {hdr.label.value} (score {hdr.score:.3f})")

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
