"""
Publication-quality figure generators for CRISPRArchitect v2 benchmarks.
=========================================================================

All figures follow journal-standard style:
- Clean sans-serif fonts (Arial / Helvetica fallback)
- High DPI (300) for print quality
- Minimal chart junk; no gridlines unless they aid interpretation
- Colorblind-friendly palette (Okabe-Ito inspired)
- Vector-compatible output (PDF preferred, PNG at 300 DPI)

Four public functions
---------------------
1. plot_summary_metrics     -- Grouped bar: aware vs naive accuracy
2. plot_ranking_shift       -- Bar chart: changed vs unchanged counts
3. plot_top_strategy_distribution -- Horizontal bar: strategy type counts
4. plot_feasibility_heatmap -- Heatmap: cases x modalities
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional, Union

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless servers
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


# ---------------------------------------------------------------------------
# Publication style defaults
# ---------------------------------------------------------------------------

# Attempt Helvetica, fall back to Arial, then generic sans-serif
_FONT_FAMILY = "sans-serif"
_FONT_SANS_SERIF = ["Helvetica", "Arial", "DejaVu Sans"]

plt.rcParams.update({
    "font.family": _FONT_FAMILY,
    "font.sans-serif": _FONT_SANS_SERIF,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.2,
    "patch.linewidth": 0.5,
})


# ---------------------------------------------------------------------------
# Colorblind-friendly palette
# ---------------------------------------------------------------------------

COLORS = {
    # Pipeline comparison
    "aware": "#2166AC",       # steel blue
    "naive": "#B2182B",       # brick red
    # Ranking shift
    "changed": "#4DAF4A",     # green
    "unchanged": "#BDBDBD",   # light gray
    # Strategy modalities
    "Base Editing": "#1B9E77", # teal
    "Prime Editing": "#D95F02",# orange
    "HDR": "#7570B3",          # purple
    "Hybrid": "#E7298A",       # pink
    "Exon Deletion": "#66A61E",# olive
    "Other": "#A6761D",        # brown
    "NONE": "#999999",         # gray (pipeline error)
}

# Heatmap colors: white = not feasible, green = feasible, amber = marginal
_HEATMAP_CMAP = mcolors.ListedColormap(["#F5F5F5", "#FFC107", "#4CAF50"])
_HEATMAP_BOUNDS = [0, 0.5, 1.5, 2.5]
_HEATMAP_NORM = mcolors.BoundaryNorm(_HEATMAP_BOUNDS, _HEATMAP_CMAP.N)


# ---------------------------------------------------------------------------
# 1. Summary metrics (aware vs naive)
# ---------------------------------------------------------------------------

def plot_summary_metrics(
    aware_summary: Dict[str, Any],
    naive_summary: Dict[str, Any],
    output_path: str,
) -> None:
    """Grouped bar chart comparing consequence-aware vs naive accuracy.

    Three metric groups: Top-1 accuracy, Top-3 accuracy, Rejection accuracy.
    Each group has two bars (naive = red, aware = blue).

    Parameters
    ----------
    aware_summary : dict
        Must contain keys ``top1_accuracy``, ``top3_accuracy``,
        ``rejection_accuracy`` (float 0-1).
    naive_summary : dict
        Same keys as ``aware_summary``.
    output_path : str
        File path for the saved figure (PDF or PNG).
    """
    labels = ["Top-1\nAccuracy", "Top-3\nAccuracy", "Rejection\nAccuracy"]
    aware_vals = [
        aware_summary.get("top1_accuracy", 0),
        aware_summary.get("top3_accuracy", 0),
        aware_summary.get("rejection_accuracy", 0),
    ]
    naive_vals = [
        naive_summary.get("top1_accuracy", 0),
        naive_summary.get("top3_accuracy", 0),
        naive_summary.get("rejection_accuracy", 0),
    ]

    x = np.arange(len(labels))
    width = 0.32

    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    bars_naive = ax.bar(
        x - width / 2, [v * 100 for v in naive_vals], width,
        label="Sequence-only", color=COLORS["naive"], alpha=0.85,
        edgecolor="white", linewidth=0.5,
    )
    bars_aware = ax.bar(
        x + width / 2, [v * 100 for v in aware_vals], width,
        label="Consequence-aware", color=COLORS["aware"], alpha=0.85,
        edgecolor="white", linewidth=0.5,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 109)
    ax.set_ylabel("Accuracy (%)")
    ax.set_title("Benchmark Performance: Aware vs Naive Scoring")
    ax.legend(loc="upper left")

    # Value labels
    for bars in [bars_naive, bars_aware]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + 1.5,
                    "%.0f%%" % h,
                    ha="center", va="bottom", fontsize=8,
                )

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. Ranking shift
# ---------------------------------------------------------------------------

def plot_ranking_shift(
    changed: Union[List[str], int],
    unchanged: Union[List[str], int],
    output_path: str,
) -> None:
    """Bar chart showing how many cases had their top strategy re-ranked
    by consequence-aware scoring.

    Parameters
    ----------
    changed : list of str or int
        Case IDs (or count) where the top strategy shifted.
    unchanged : list of str or int
        Case IDs (or count) where the top strategy stayed the same.
    output_path : str
        Output file path.
    """
    n_changed = len(changed) if isinstance(changed, list) else int(changed)
    n_unchanged = len(unchanged) if isinstance(unchanged, list) else int(unchanged)

    categories = ["Ranking\nChanged", "Unchanged"]
    values = [n_changed, n_unchanged]
    colors = [COLORS["changed"], COLORS["unchanged"]]

    fig, ax = plt.subplots(figsize=(4.0, 3.5))
    bars = ax.bar(categories, values, color=colors, width=0.55,
                  edgecolor="white", linewidth=0.5)

    ax.set_ylabel("Number of Cases")
    ax.set_title("Effect of Consequence-Aware Scoring")

    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            str(val),
            ha="center", va="bottom", fontsize=11, fontweight="bold",
        )

    total = n_changed + n_unchanged
    if total > 0:
        pct = n_changed / total * 100
        ax.text(
            0.95, 0.92,
            "%.0f%% of cases\nre-ranked" % pct,
            transform=ax.transAxes,
            ha="right", va="top", fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="#FFFDE7", edgecolor="#BDBDBD",
                alpha=0.9,
            ),
        )

    ax.set_ylim(0, max(values) * 1.25 if values else 1)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3. Top strategy distribution
# ---------------------------------------------------------------------------

def plot_top_strategy_distribution(
    case_results: List[Dict[str, Any]],
    output_path: str,
) -> None:
    """Horizontal bar chart showing distribution of top-ranked strategy types.

    Strategies are grouped into categories: Base Editing, Prime Editing,
    HDR, Hybrid, Exon Deletion, Other.

    Parameters
    ----------
    case_results : list of dict
        Each dict must have ``top_strategy`` (str). Typically the
        ``case_results`` list from ``BenchmarkSummary.to_dict()``.
    output_path : str
        Output file path.
    """
    strategies = [
        r.get("top_strategy", "NONE") for r in case_results
    ]

    categories = [_categorize_strategy(s) for s in strategies]
    counter = Counter(categories)

    # Sort by count descending
    sorted_items = sorted(counter.items(), key=lambda x: x[1])
    labels = [item[0] for item in sorted_items]
    values = [item[1] for item in sorted_items]
    bar_colors = [COLORS.get(l, COLORS["Other"]) for l in labels]

    fig, ax = plt.subplots(figsize=(5.5, max(2.5, len(labels) * 0.55)))
    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, values, color=bar_colors, height=0.6,
                   edgecolor="white", linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Number of Cases")
    ax.set_title("Distribution of Top-Ranked Strategies")

    # Value labels
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
            str(val), ha="left", va="center", fontsize=9,
        )

    ax.set_xlim(0, max(values) * 1.2 if values else 1)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def _categorize_strategy(name: str) -> str:
    """Map a strategy name to a display category."""
    n = name.lower()
    if n in ("none", ""):
        return "NONE"
    if "hybrid" in n:
        return "Hybrid"
    if "dual" in n and ("base" in n or "abe" in n or "cbe" in n):
        return "Base Editing"
    if "base" in n or "abe" in n or "cbe" in n:
        return "Base Editing"
    if "dual" in n and "prime" in n:
        return "Prime Editing"
    if "prime" in n:
        return "Prime Editing"
    if "exon" in n:
        return "Exon Deletion"
    if "hdr" in n:
        return "HDR"
    return "Other"


# ---------------------------------------------------------------------------
# 4. Feasibility heatmap
# ---------------------------------------------------------------------------

def plot_feasibility_heatmap(
    case_results: List[Dict[str, Any]],
    output_path: str,
    max_cases: int = 30,
) -> None:
    """Heatmap showing feasibility across editing modalities for each case.

    Rows = benchmark cases (capped at ``max_cases``).
    Columns = BE, PE, HDR modalities.
    Cell values: 0 = not feasible, 1 = marginal, 2 = feasible.

    Feasibility is inferred from which modalities appear in the ranked
    strategy list for each case.

    Parameters
    ----------
    case_results : list of dict
        Each dict must have ``case_id``, ``top_strategy``, ``top3_strategies``.
    output_path : str
        Output file path.
    max_cases : int
        Maximum number of cases to display (default 30).
    """
    n = min(len(case_results), max_cases)
    if n == 0:
        return

    modalities = ["BE", "PE", "HDR"]
    data = np.zeros((n, 3), dtype=float)
    labels_y = []

    for i in range(n):
        r = case_results[i]
        case_id = r.get("case_id", "Case %d" % (i + 1))
        labels_y.append(case_id)

        top = r.get("top_strategy", "").lower()
        top3 = [s.lower() for s in r.get("top3_strategies", [])]
        all_strats = top3  # use top3 as available set

        # BE feasibility
        if any("base" in s or "abe" in s or "cbe" in s for s in all_strats):
            if "base" in top or "abe" in top or "cbe" in top:
                data[i, 0] = 2  # feasible (top-ranked)
            else:
                data[i, 0] = 1  # marginal (present but not top)

        # PE feasibility
        if any("prime" in s for s in all_strats):
            if "prime" in top:
                data[i, 1] = 2
            else:
                data[i, 1] = 1

        # HDR feasibility
        if any("hdr" in s for s in all_strats):
            if "hdr" in top:
                data[i, 2] = 2
            else:
                data[i, 2] = 1

    # Figure dimensions scale with case count
    fig_height = max(3.0, n * 0.28 + 1.0)
    fig, ax = plt.subplots(figsize=(4.5, fig_height))

    im = ax.imshow(
        data, aspect="auto", cmap=_HEATMAP_CMAP, norm=_HEATMAP_NORM,
        interpolation="nearest",
    )

    ax.set_xticks(range(3))
    ax.set_xticklabels(modalities, fontweight="bold")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels_y, fontsize=max(5, 9 - n // 10))
    ax.set_title("Feasibility Across Editing Modalities")

    # Cell text annotations
    _feasibility_labels = {0: "", 1: "?", 2: "+"}
    for row in range(n):
        for col in range(3):
            val = int(data[row, col])
            text = _feasibility_labels.get(val, "")
            if text:
                text_color = "white" if val == 2 else "#333333"
                ax.text(
                    col, row, text,
                    ha="center", va="center",
                    fontsize=8, fontweight="bold",
                    color=text_color,
                )

    # Colorbar with custom labels
    cbar = plt.colorbar(
        im, ax=ax, shrink=0.4, pad=0.02,
        ticks=[0.75, 1.5, 2.25],
    )
    cbar.ax.set_yticklabels(
        ["Not feasible", "Marginal", "Feasible"],
        fontsize=7,
    )

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile
    import os

    print("Testing plotting.py...")

    tmpdir = tempfile.mkdtemp(prefix="crisprarchitect_plot_test_")

    # 1. Summary metrics
    out1 = os.path.join(tmpdir, "test_summary_metrics.pdf")
    plot_summary_metrics(
        {"top1_accuracy": 0.83, "top3_accuracy": 0.93,
         "rejection_accuracy": 0.97},
        {"top1_accuracy": 0.60, "top3_accuracy": 0.77,
         "rejection_accuracy": 0.83},
        out1,
    )
    assert os.path.exists(out1), "Summary metrics plot not created"
    print("  plot_summary_metrics: OK (%s)" % out1)

    # 2. Ranking shift
    out2 = os.path.join(tmpdir, "test_ranking_shift.pdf")
    plot_ranking_shift(
        ["CASE_001", "CASE_005", "CASE_012"],
        ["CASE_002", "CASE_003", "CASE_004", "CASE_006", "CASE_007"],
        out2,
    )
    assert os.path.exists(out2), "Ranking shift plot not created"
    print("  plot_ranking_shift: OK (%s)" % out2)

    # 3. Strategy distribution
    out3 = os.path.join(tmpdir, "test_strategy_distribution.pdf")
    dummy_cases = [
        {"case_id": "C%d" % i, "top_strategy": s,
         "top3_strategies": [s]}
        for i, s in enumerate([
            "Single-step Base Editing", "Single-step Base Editing",
            "Single-step Prime Editing", "Single-step HDR",
            "Dual Base Editing", "Hybrid ABE + HDR",
            "Sequential HDR", "NONE",
        ])
    ]
    plot_top_strategy_distribution(dummy_cases, out3)
    assert os.path.exists(out3), "Strategy distribution plot not created"
    print("  plot_top_strategy_distribution: OK (%s)" % out3)

    # 4. Feasibility heatmap
    out4 = os.path.join(tmpdir, "test_feasibility_heatmap.pdf")
    plot_feasibility_heatmap(dummy_cases, out4)
    assert os.path.exists(out4), "Feasibility heatmap not created"
    print("  plot_feasibility_heatmap: OK (%s)" % out4)

    print("All plotting self-tests: PASS")
    print("Test outputs in: %s" % tmpdir)
