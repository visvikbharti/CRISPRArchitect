#!/usr/bin/env python3
"""
CRISPRArchitect v2 Benchmark Runner
=====================================

CLI script that loads a benchmark dataset, runs the evaluator, prints
a summary table, and saves results as JSON and (optionally) plots.

Usage
-----
    python -m benchmarks.run_benchmark \
        --input benchmarks/dataset_v1.json \
        --output-dir results/benchmark_v1 \
        --prefix v1_run1

    # Compare consequence-aware vs naive:
    python -m benchmarks.run_benchmark \
        --input benchmarks/dataset_v1.json \
        --output-dir results/benchmark_v1 \
        --compare

Outputs
-------
    <output-dir>/<prefix>_results.json        Summary metrics + per-case
    <output-dir>/<prefix>_comparison.json      Aware vs naive (if --compare)
    <output-dir>/<prefix>_summary_metrics.pdf  Bar chart (if --compare)
    <output-dir>/<prefix>_ranking_shift.pdf    Shift chart (if --compare)
    <output-dir>/<prefix>_strategy_distribution.pdf
    <output-dir>/<prefix>_feasibility_heatmap.pdf
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

# Ensure project root is on sys.path so that ``core`` and ``benchmarks``
# are importable regardless of how the script is invoked.
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from benchmarks.evaluator import BenchmarkEvaluator, load_cases

logger = logging.getLogger("benchmarks.run_benchmark")


# ---------------------------------------------------------------------------
# Table formatting
# ---------------------------------------------------------------------------

def _format_summary_table(summary_dict: Dict[str, Any]) -> str:
    """Format a BenchmarkSummary dict as an ASCII table.

    Parameters
    ----------
    summary_dict : dict
        Output of ``BenchmarkSummary.to_dict()``.

    Returns
    -------
    str
        Multi-line table string.
    """
    lines = []
    lines.append("")
    lines.append("=" * 68)
    lines.append("  CRISPRArchitect v2 Benchmark Summary")
    lines.append("=" * 68)
    lines.append("")
    lines.append(
        "  Total cases:                 %d" % summary_dict["n_cases"]
    )
    lines.append(
        "  Top-1 accuracy:              %.1f%%"
        % (summary_dict["top1_accuracy"] * 100)
    )
    lines.append(
        "  Top-3 accuracy:              %.1f%%"
        % (summary_dict["top3_accuracy"] * 100)
    )
    lines.append(
        "  Rejection accuracy:          %.1f%%"
        % (summary_dict["rejection_accuracy"] * 100)
    )
    lines.append(
        "  Consequence-shift fraction:  %.1f%%"
        % (summary_dict["consequence_shift_fraction"] * 100)
    )
    lines.append("")

    # Per-case details
    case_results = summary_dict.get("case_results", [])
    if case_results:
        lines.append("-" * 68)
        header = "  %-20s %-6s %-6s %-6s %s" % (
            "Case ID", "Top1", "Top3", "Rej", "Top Strategy"
        )
        lines.append(header)
        lines.append("-" * 68)

        for cr in case_results:
            case_id = cr.get("case_id", "?")
            top1 = "Y" if cr.get("top1_correct") else "N"
            top3 = "Y" if cr.get("top3_correct") else "N"
            rej = "Y" if cr.get("rejected_correctly") else "N"
            top_strat = cr.get("top_strategy", "")
            error = cr.get("error", cr.get("pipeline_error", ""))

            if error:
                top_strat = "[ERROR] %s" % error[:30]

            lines.append(
                "  %-20s %-6s %-6s %-6s %s"
                % (case_id[:20], top1, top3, rej, top_strat[:30])
            )

        lines.append("-" * 68)

    lines.append("")
    return "\n".join(lines)


def _format_comparison_table(
    aware_dict: Dict[str, Any],
    naive_dict: Dict[str, Any],
    changed: List[str],
    unchanged: List[str],
) -> str:
    """Format a comparison table for aware vs naive pipelines.

    Parameters
    ----------
    aware_dict : dict
        Summary dict for the consequence-aware pipeline.
    naive_dict : dict
        Summary dict for the naive (no consequence penalties) pipeline.
    changed : list of str
        Case IDs where the top strategy shifted.
    unchanged : list of str
        Case IDs where the top strategy stayed the same.

    Returns
    -------
    str
        Multi-line table string.
    """
    lines = []
    lines.append("")
    lines.append("=" * 68)
    lines.append("  Consequence-Aware vs Naive Comparison")
    lines.append("=" * 68)
    lines.append("")
    lines.append(
        "  %-30s %-15s %-15s" % ("Metric", "Aware", "Naive")
    )
    lines.append("-" * 68)
    lines.append(
        "  %-30s %-15.1f %-15.1f"
        % (
            "Top-1 accuracy (%)",
            aware_dict["top1_accuracy"] * 100,
            naive_dict["top1_accuracy"] * 100,
        )
    )
    lines.append(
        "  %-30s %-15.1f %-15.1f"
        % (
            "Top-3 accuracy (%)",
            aware_dict["top3_accuracy"] * 100,
            naive_dict["top3_accuracy"] * 100,
        )
    )
    lines.append(
        "  %-30s %-15.1f %-15.1f"
        % (
            "Rejection accuracy (%)",
            aware_dict["rejection_accuracy"] * 100,
            naive_dict["rejection_accuracy"] * 100,
        )
    )
    lines.append("-" * 68)
    lines.append("  Cases with ranking shift:    %d" % len(changed))
    lines.append("  Cases without shift:         %d" % len(unchanged))

    if changed:
        lines.append("")
        lines.append("  Shifted cases: %s" % ", ".join(changed[:10]))
        if len(changed) > 10:
            lines.append("    ... and %d more" % (len(changed) - 10))

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plot generation helpers
# ---------------------------------------------------------------------------

def _generate_standard_plots(
    summary_dict: Dict[str, Any],
    case_results: list,
    output_dir: Path,
    prefix: str,
) -> None:
    """Generate standard benchmark plots (strategy distribution + heatmap).

    Parameters
    ----------
    summary_dict : dict
        Serialized benchmark summary.
    case_results : list
        List of BenchmarkResult objects from the evaluator.
    output_dir : Path
        Directory for output files.
    prefix : str
        Filename prefix.
    """
    try:
        from benchmarks.plotting import (
            plot_top_strategy_distribution,
            plot_feasibility_heatmap,
        )

        dist_path = output_dir / ("%s_strategy_distribution.pdf" % prefix)
        plot_top_strategy_distribution(
            summary_dict.get("case_results", []), str(dist_path)
        )
        logger.info("Saved strategy distribution to %s", dist_path)

        heatmap_path = output_dir / ("%s_feasibility_heatmap.pdf" % prefix)
        plot_feasibility_heatmap(
            summary_dict.get("case_results", []), str(heatmap_path)
        )
        logger.info("Saved feasibility heatmap to %s", heatmap_path)

    except ImportError as e:
        logger.warning("Plotting unavailable (%s). Skipping plots.", e)
    except Exception as e:
        logger.warning("Plot generation failed: %s", e)


def _generate_comparison_plots(
    aware_dict: Dict[str, Any],
    naive_dict: Dict[str, Any],
    changed: List[str],
    unchanged: List[str],
    case_results_dicts: List[Dict[str, Any]],
    output_dir: Path,
    prefix: str,
) -> None:
    """Generate all comparison plots (aware vs naive).

    Parameters
    ----------
    aware_dict : dict
        Serialized aware summary.
    naive_dict : dict
        Serialized naive summary.
    changed : list of str
        Case IDs where ranking shifted.
    unchanged : list of str
        Case IDs where ranking stayed.
    case_results_dicts : list of dict
        Serialized per-case results.
    output_dir : Path
        Directory for output files.
    prefix : str
        Filename prefix.
    """
    try:
        from benchmarks.plotting import (
            plot_summary_metrics,
            plot_ranking_shift,
            plot_top_strategy_distribution,
            plot_feasibility_heatmap,
        )

        # Summary bar chart
        summary_path = output_dir / ("%s_summary_metrics.pdf" % prefix)
        plot_summary_metrics(aware_dict, naive_dict, str(summary_path))
        logger.info("Saved summary metrics to %s", summary_path)

        # Ranking shift chart
        shift_path = output_dir / ("%s_ranking_shift.pdf" % prefix)
        plot_ranking_shift(changed, unchanged, str(shift_path))
        logger.info("Saved ranking shift to %s", shift_path)

        # Strategy distribution
        dist_path = output_dir / ("%s_strategy_distribution.pdf" % prefix)
        plot_top_strategy_distribution(case_results_dicts, str(dist_path))
        logger.info("Saved strategy distribution to %s", dist_path)

        # Feasibility heatmap
        heatmap_path = output_dir / ("%s_feasibility_heatmap.pdf" % prefix)
        plot_feasibility_heatmap(case_results_dicts, str(heatmap_path))
        logger.info("Saved feasibility heatmap to %s", heatmap_path)

    except ImportError as e:
        logger.warning("Plotting unavailable (%s). Skipping plots.", e)
    except Exception as e:
        logger.warning("Plot generation failed: %s", e)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: List[str] = None) -> int:
    """Entry point for the benchmark runner.

    Parameters
    ----------
    argv : list of str or None
        Command-line arguments (defaults to sys.argv[1:]).

    Returns
    -------
    int
        Exit code (0 = success, 1 = error).
    """
    parser = argparse.ArgumentParser(
        prog="run_benchmark",
        description=(
            "CRISPRArchitect v2 Benchmark Runner: evaluate pipeline "
            "performance against curated truth labels."
        ),
    )
    parser.add_argument(
        "--input",
        type=str,
        default=str(_SCRIPT_DIR / "dataset_v1.json"),
        help=(
            "Path to benchmark dataset JSON "
            "(default: benchmarks/dataset_v1.json)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(_PROJECT_ROOT / "results" / "benchmark"),
        help="Output directory for results (default: results/benchmark/).",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="benchmark",
        help="Filename prefix for output files (default: 'benchmark').",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Run aware-vs-naive comparison and generate shift plots.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation (useful for headless environments).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=True,
        help="Print per-case progress (default: True).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-case logging.",
    )

    args = parser.parse_args(argv)

    # Configure logging
    log_level = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    verbose = not args.quiet

    # --- Load dataset ---
    input_path = args.input
    logger.info("Loading benchmark dataset from %s", input_path)
    try:
        cases = load_cases(input_path)
    except (FileNotFoundError, ValueError) as e:
        logger.error("Failed to load dataset: %s", e)
        print("ERROR: %s" % e, file=sys.stderr)
        return 1

    logger.info("Loaded %d benchmark cases", len(cases))
    print("Loaded %d benchmark cases from %s" % (len(cases), input_path))

    # --- Create output directory ---
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = args.prefix
    evaluator = BenchmarkEvaluator(verbose=verbose)

    # --- Run benchmark ---
    t0 = time.time()

    if args.compare:
        print("Running consequence-aware vs naive comparison...")
        comparison = evaluator.compare_aware_vs_naive(cases)

        aware_summary = comparison["aware_summary"]
        naive_summary = comparison["naive_summary"]
        changed = comparison["changed_cases"]
        unchanged = comparison["unchanged_cases"]

        aware_dict = aware_summary.to_dict()
        naive_dict = naive_summary.to_dict()

        # Print tables
        print(_format_summary_table(aware_dict))
        print(_format_comparison_table(
            aware_dict, naive_dict, changed, unchanged
        ))

        # Save JSON
        output_json = output_dir / ("%s_comparison.json" % prefix)
        with open(str(output_json), "w") as f:
            json.dump(
                {
                    "aware": aware_dict,
                    "naive": naive_dict,
                    "changed_cases": changed,
                    "unchanged_cases": unchanged,
                },
                f,
                indent=2,
                default=str,
            )
        print("Saved comparison results to %s" % output_json)

        # Plots
        if not args.no_plots:
            _generate_comparison_plots(
                aware_dict, naive_dict, changed, unchanged,
                aware_dict.get("case_results", []),
                output_dir, prefix,
            )

    else:
        print("Running standard benchmark evaluation...")
        summary = evaluator.evaluate_cases(cases)
        summary_dict = summary.to_dict()

        # Print table
        print(_format_summary_table(summary_dict))

        # Save JSON
        output_json = output_dir / ("%s_results.json" % prefix)
        with open(str(output_json), "w") as f:
            json.dump(summary_dict, f, indent=2, default=str)
        print("Saved results to %s" % output_json)

        # Plots
        if not args.no_plots:
            _generate_standard_plots(
                summary_dict, summary.case_results, output_dir, prefix,
            )

    elapsed = time.time() - t0
    print("Total wall time: %.1f seconds" % elapsed)
    print("Benchmark complete.")

    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
