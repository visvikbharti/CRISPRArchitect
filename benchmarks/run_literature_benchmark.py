#!/usr/bin/env python3
"""
Run CRISPRArchitect v3 pipeline on the literature benchmark cases.

For each case with a parseable HGVS variant, runs the full pipeline
and compares the tool's top-ranked strategy against the experimentally
used strategy from the published paper.

Output: JSON results file + printed concordance report.
"""

import json
import re
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.pipeline.strategy_stage import StrategyPipeline
from core.models import GenomicVariantInput


def parse_hgvs_from_benchmark(case):
    """Extract transcript, c.notation, gene, ref, alt from benchmark case."""
    hgvs = case.get("variant_hgvs", "")
    gene = case.get("gene", "")

    # Match NM_XXXXXX.X:c.NNNX>Y
    match = re.search(r"(NM_\d+\.\d+):?(c\.(\d+)([ACGT])>([ACGT]))", hgvs)
    if not match:
        return None

    transcript = match.group(1)
    c_notation = match.group(2)
    position = int(match.group(3))
    ref = match.group(4)
    alt = match.group(5)

    return {
        "transcript": transcript,
        "c_notation": c_notation,
        "cds_position": position,
        "ref": ref,
        "alt": alt,
        "gene": gene,
    }


def classify_published_strategy(case):
    """Map published strategy to our modality categories."""
    strategy = case.get("strategy", "").lower()
    editor = case.get("editor", "").lower()

    if "base editing" in strategy or "abe" in editor or "cbe" in editor or "be3" in editor:
        return "BE"
    elif "prime" in strategy or "pe2" in editor or "pe3" in editor:
        return "PE"
    elif "hdr" in strategy:
        return "HDR"
    elif "multi" in strategy or "comparison" in strategy:
        return "COMPARISON"
    else:
        return "UNKNOWN"


def classify_pipeline_strategy(strategy_name):
    """Map pipeline strategy name to our categories."""
    name = strategy_name.lower()
    if "base edit" in name:
        return "BE"
    elif "prime edit" in name:
        return "PE"
    elif "hdr" in name:
        return "HDR"
    else:
        return "UNKNOWN"


def run_benchmark():
    """Run the literature benchmark."""
    # Load cases
    benchmark_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "literature_benchmark_v1.json"
    )
    with open(benchmark_path) as f:
        data = json.load(f)

    results = []
    skipped = []
    errors = []

    # Deduplicate: for the same variant (gene + HGVS), run the pipeline once
    # and compare against multiple published strategies
    variant_cache = {}  # key: (gene, c_notation) -> pipeline result

    print("=" * 75)
    print("  CRISPRArchitect v3 — Literature Benchmark Evaluation")
    print("=" * 75)
    print()

    for case in data["cases"]:
        case_id = case["case_id"]
        parsed = parse_hgvs_from_benchmark(case)

        if parsed is None:
            skipped.append({
                "case_id": case_id,
                "gene": case.get("gene", "?"),
                "reason": "No parseable HGVS c. notation",
            })
            continue

        # Check if non-human
        cell_type = case.get("cell_type", "").lower()
        if "mouse" in cell_type or "macaque" in cell_type:
            skipped.append({
                "case_id": case_id,
                "gene": parsed["gene"],
                "reason": f"Non-human cell type: {case.get('cell_type', '?')}",
            })
            continue

        # Check if it's a disruption/inactivation (not correction)
        hgvs_text = case.get("variant_hgvs", "").lower()
        if "disruption" in hgvs_text or "inactivat" in hgvs_text:
            skipped.append({
                "case_id": case_id,
                "gene": parsed["gene"],
                "reason": "Therapeutic disruption, not correction",
            })
            continue

        cache_key = (parsed["gene"], parsed["c_notation"])
        published_category = classify_published_strategy(case)

        print(f"  [{case_id}] {parsed['gene']} {parsed['c_notation']} "
              f"(published: {published_category})")

        # Run pipeline (use cache for same variant)
        if cache_key in variant_cache:
            pipeline_result = variant_cache[cache_key]
            print(f"    Using cached result for {cache_key[0]}:{cache_key[1]}")
        else:
            try:
                # Determine cell type for pipeline
                pipeline_cell = "iPSC"  # default
                if "hek" in cell_type or "293" in cell_type:
                    pipeline_cell = "HEK293T"
                elif "k562" in cell_type:
                    pipeline_cell = "K562"
                elif "hspc" in cell_type or "cd34" in cell_type or "hsc" in cell_type:
                    pipeline_cell = "HSC"

                pipeline = StrategyPipeline(
                    cell_type=pipeline_cell,
                    nuclease="SpCas9",
                )

                # Create variant input from HGVS
                variant = GenomicVariantInput(
                    chromosome="",  # will be resolved from gene
                    position=0,  # will be resolved
                    ref_allele=parsed["ref"],
                    alt_allele=parsed["alt"],
                    gene_symbol=parsed["gene"],
                    name=f"{parsed['gene']}:{parsed['c_notation']}",
                )

                # Run pipeline
                start_time = time.time()
                pipeline_result = pipeline.run([variant])
                elapsed = time.time() - start_time

                variant_cache[cache_key] = pipeline_result
                print(f"    Pipeline completed in {elapsed:.1f}s")

                # Brief delay for API rate limiting
                time.sleep(0.5)

            except Exception as e:
                error_msg = str(e)[:200]
                print(f"    ERROR: {error_msg}")
                errors.append({
                    "case_id": case_id,
                    "gene": parsed["gene"],
                    "error": error_msg,
                })
                continue

        # Extract top strategy from pipeline
        if pipeline_result.strategies:
            top = pipeline_result.strategies[0]
            top_category = classify_pipeline_strategy(top.strategy_name)
            top_score = top.overall_score

            # Check for Pareto and rank stability
            is_pareto = any("non-dominated" in n for n in top.annotation_notes)
            rank_stability = top.rank_stability if hasattr(top, 'rank_stability') and top.rank_stability else None

            concordant = (top_category == published_category) or (published_category == "COMPARISON")

            result = {
                "case_id": case_id,
                "gene": parsed["gene"],
                "variant": parsed["c_notation"],
                "cell_type": case.get("cell_type", "?"),
                "published_strategy": published_category,
                "published_editor": case.get("editor", "?"),
                "published_efficiency": case.get("editing_efficiency_pct", "?"),
                "pipeline_top1": top.strategy_name,
                "pipeline_top1_category": top_category,
                "pipeline_top1_score": round(top_score, 4),
                "pipeline_pareto": is_pareto,
                "pipeline_rank_stability": round(rank_stability, 3) if rank_stability else None,
                "concordant": concordant,
                "pmid": case.get("pmid", "?"),
            }

            # Also get top-3 categories
            top3_cats = [classify_pipeline_strategy(s.strategy_name)
                         for s in pipeline_result.strategies[:3]]
            result["pipeline_top3_categories"] = top3_cats
            result["top3_concordant"] = published_category in top3_cats or published_category == "COMPARISON"

            results.append(result)

            status = "CONCORDANT" if concordant else "DISCORDANT"
            print(f"    Pipeline recommends: {top.strategy_name} ({top_category}) "
                  f"score={top_score:.3f}")
            print(f"    Published used: {published_category} | {status}")
        else:
            print(f"    No strategies generated")
            errors.append({
                "case_id": case_id,
                "gene": parsed["gene"],
                "error": "No strategies generated",
            })

        print()

    # Summary
    print("=" * 75)
    print("  BENCHMARK SUMMARY")
    print("=" * 75)
    print(f"  Total cases in dataset:      {len(data['cases'])}")
    print(f"  Cases evaluated:             {len(results)}")
    print(f"  Cases skipped:               {len(skipped)}")
    print(f"  Cases with errors:           {len(errors)}")
    print()

    if results:
        n_concordant = sum(1 for r in results if r["concordant"])
        n_top3_concordant = sum(1 for r in results if r["top3_concordant"])
        n_total = len(results)

        print(f"  Top-1 concordance:           {n_concordant}/{n_total} "
              f"({100*n_concordant/n_total:.1f}%)")
        print(f"  Top-3 concordance:           {n_top3_concordant}/{n_total} "
              f"({100*n_top3_concordant/n_total:.1f}%)")

        # Per-category breakdown
        print(f"\n  Per-category concordance:")
        for cat in ["BE", "PE", "HDR", "COMPARISON"]:
            cat_results = [r for r in results if r["published_strategy"] == cat]
            if cat_results:
                cat_concordant = sum(1 for r in cat_results if r["concordant"])
                print(f"    {cat:12s}: {cat_concordant}/{len(cat_results)}")

        # Pipeline recommendation distribution
        print(f"\n  Pipeline recommendation distribution:")
        rec_counts = {}
        for r in results:
            cat = r["pipeline_top1_category"]
            rec_counts[cat] = rec_counts.get(cat, 0) + 1
        for cat, count in sorted(rec_counts.items(), key=lambda x: -x[1]):
            print(f"    {cat:12s}: {count}/{n_total} ({100*count/n_total:.1f}%)")

    # Save results
    output = {
        "benchmark_version": "literature_v1",
        "pipeline_version": "3.0.0",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_total_cases": len(data["cases"]),
        "n_evaluated": len(results),
        "n_skipped": len(skipped),
        "n_errors": len(errors),
        "results": results,
        "skipped": skipped,
        "errors": errors,
    }

    if results:
        n_concordant = sum(1 for r in results if r["concordant"])
        n_top3 = sum(1 for r in results if r["top3_concordant"])
        output["top1_concordance"] = round(n_concordant / len(results), 4)
        output["top3_concordance"] = round(n_top3 / len(results), 4)

    output_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "benchmark_results",
        "literature_benchmark_results.json",
    )
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  Results saved to: {output_path}")
    print("=" * 75)


if __name__ == "__main__":
    run_benchmark()
