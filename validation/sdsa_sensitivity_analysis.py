#!/usr/bin/env python3
"""
SDSA Displacement Probability Sensitivity Analysis
====================================================

Sweeps the central [ASSUMED] parameter p (SDSA_DISPLACEMENT_PROB_PER_BP)
from 0.001 to 0.005 and quantifies the effect on:
  - Mean and median gene conversion tract lengths
  - Distance-dependent HDR probability
  - Overall HDR rates across cell types and donor topologies

Motivation: p=0.002 is the default, derived from functional evidence
(Stark lab SDSA assay, successful HDR with 300-1000bp HA), but NOT
directly measured. This analysis demonstrates whether CRISPRArchitect's
HDR recommendations are robust to 5-fold variation in this parameter.

Output:
  - validation/figures/Fig_SDSA_Sensitivity.png (publication quality)
  - validation/figures/Fig_SDSA_Sensitivity.pdf
  - Console summary table

Usage:
    cd crisprarchitect/
    python3 validation/sdsa_sensitivity_analysis.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from utils import constants
from conversion_sim import ConversionSimulator
from conversion_sim import synthesis as _synthesis_mod

# Try matplotlib
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not available, skipping figure generation")


# =========================================================================
# Configuration
# =========================================================================
P_VALUES = [0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005]
N_SIMULATIONS = 10_000
SEED = 42

# Configurations to test
CONFIGS = [
    {
        "name": "SpCas9 + dsDNA (iPSC)",
        "cut_type": "blunt",
        "overhang_length": 0,
        "donor_topology": "linear_dsDNA",
        "homology_arm_length": 800,
        "cell_type": "iPSC",
        "color": "#3498DB",
        "linestyle": "-",
    },
    {
        "name": "SpCas9 + cssDNA (iPSC)",
        "cut_type": "blunt",
        "overhang_length": 0,
        "donor_topology": "circular_ssDNA",
        "homology_arm_length": 300,
        "cell_type": "iPSC",
        "color": "#1B9E77",
        "linestyle": "-",
    },
    {
        "name": "enFnCas9 + cssDNA (iPSC)",
        "cut_type": "staggered_5prime",
        "overhang_length": 3,
        "donor_topology": "circular_ssDNA",
        "homology_arm_length": 300,
        "cell_type": "iPSC",
        "color": "#D95F02",
        "linestyle": "-",
    },
    {
        "name": "enFnCas9 + cssDNA (HEK293T)",
        "cut_type": "staggered_5prime",
        "overhang_length": 3,
        "donor_topology": "circular_ssDNA",
        "homology_arm_length": 300,
        "cell_type": "HEK293T",
        "color": "#D95F02",
        "linestyle": "--",
    },
]

# Distance thresholds to check
DISTANCE_THRESHOLDS = [100, 200, 500, 1000, 2000]


# =========================================================================
# Run sweep
# =========================================================================
def run_sensitivity_sweep():
    """Run the full p-value sweep and return results."""
    results = {cfg["name"]: {} for cfg in CONFIGS}

    for p_val in P_VALUES:
        # Monkey-patch the constant in BOTH the constants module
        # and the synthesis module (which imports it by value)
        constants.SDSA_DISPLACEMENT_PROB_PER_BP = p_val
        _synthesis_mod.SDSA_DISPLACEMENT_PROB_PER_BP = p_val

        for cfg in CONFIGS:
            sim = ConversionSimulator(
                cut_type=cfg["cut_type"],
                overhang_length=cfg["overhang_length"],
                donor_topology=cfg["donor_topology"],
                homology_arm_length=cfg["homology_arm_length"],
                cell_type=cfg["cell_type"],
                n_simulations=N_SIMULATIONS,
                seed=SEED,
            )
            sim.run()
            # Redirect stdout to suppress summary output
            import io, contextlib
            f_buf = io.StringIO()
            with contextlib.redirect_stdout(f_buf):
                summary = sim.summary()

            # Extract tract lengths from simulation results
            tract_lengths = sim._results.tract_lengths_bp[sim._results.hdr_success]
            # If no successful tracts, use all tract lengths for stats
            if len(tract_lengths) == 0:
                tract_lengths = sim._results.tract_lengths_bp

            # Compute statistics
            mean_tract = np.mean(tract_lengths)
            median_tract = np.median(tract_lengths)
            std_tract = np.std(tract_lengths, ddof=1)
            se_mean = std_tract / np.sqrt(len(tract_lengths))

            # Distance probabilities
            dist_probs = {}
            for d in DISTANCE_THRESHOLDS:
                frac = np.mean(tract_lengths >= d)
                n = len(tract_lengths)
                se = np.sqrt(frac * (1 - frac) / n)
                dist_probs[d] = {"prob": frac, "se": se}

            # HDR rate
            hdr_rate = summary.get("hdr_success_rate", 0.0)
            hdr_ci = (summary.get("hdr_rate_ci95_low", 0.0),
                      summary.get("hdr_rate_ci95_high", 0.0))

            results[cfg["name"]][p_val] = {
                "mean_tract": mean_tract,
                "median_tract": median_tract,
                "std_tract": std_tract,
                "se_mean": se_mean,
                "dist_probs": dist_probs,
                "hdr_rate": hdr_rate,
                "hdr_ci_95": hdr_ci,
            }

    # Restore default
    constants.SDSA_DISPLACEMENT_PROB_PER_BP = 0.002
    _synthesis_mod.SDSA_DISPLACEMENT_PROB_PER_BP = 0.002

    return results


def print_summary_table(results):
    """Print a formatted summary table."""
    print("\n" + "=" * 100)
    print("SDSA DISPLACEMENT PROBABILITY SENSITIVITY ANALYSIS")
    print("=" * 100)
    print(f"Parameter swept: SDSA_DISPLACEMENT_PROB_PER_BP")
    print(f"Range: {P_VALUES[0]} to {P_VALUES[-1]} ({len(P_VALUES)} values)")
    print(f"Simulations per point: {N_SIMULATIONS:,}")
    print(f"Default value: 0.002 (mean tract = 500 bp)")
    print()

    for cfg_name, p_results in results.items():
        print(f"\n--- {cfg_name} ---")
        print(f"{'p':>8s}  {'Mean(bp)':>10s}  {'Median(bp)':>10s}  "
              f"{'P(≥500)':>8s}  {'P(≥1000)':>8s}  {'HDR rate':>10s}")
        print("-" * 70)

        for p_val in P_VALUES:
            r = p_results[p_val]
            p500 = r["dist_probs"][500]["prob"]
            p1000 = r["dist_probs"][1000]["prob"]
            hdr = r["hdr_rate"]
            ci = r["hdr_ci_95"]
            print(f"{p_val:>8.4f}  {r['mean_tract']:>10.0f}  {r['median_tract']:>10.0f}  "
                  f"{p500:>8.1%}  {p1000:>8.1%}  "
                  f"{hdr:>6.1%} [{ci[0]:.1%}-{ci[1]:.1%}]")

    # Robustness assessment
    print("\n" + "=" * 100)
    print("ROBUSTNESS ASSESSMENT")
    print("=" * 100)

    for cfg_name, p_results in results.items():
        default_r = p_results[0.002]
        min_r = p_results[P_VALUES[0]]
        max_r = p_results[P_VALUES[-1]]

        mean_range = max_r["mean_tract"] - min_r["mean_tract"]
        hdr_range = abs(max_r["hdr_rate"] - min_r["hdr_rate"])

        print(f"\n{cfg_name}:")
        print(f"  Mean tract range: {min_r['mean_tract']:.0f} - "
              f"{max_r['mean_tract']:.0f} bp "
              f"(default: {default_r['mean_tract']:.0f} bp)")
        print(f"  HDR rate range:   {min(min_r['hdr_rate'], max_r['hdr_rate']):.1%} - "
              f"{max(min_r['hdr_rate'], max_r['hdr_rate']):.1%} "
              f"(default: {default_r['hdr_rate']:.1%})")

        # Check if qualitative recommendation changes
        p500_default = default_r["dist_probs"][500]["prob"]
        p500_min = min_r["dist_probs"][500]["prob"]
        p500_max = max_r["dist_probs"][500]["prob"]

        if p500_min > 0.2 and p500_max > 0.2:
            print(f"  P(tract≥500bp): {p500_max:.1%} - {p500_min:.1%} "
                  f"→ ROBUST (always >20%)")
        elif p500_min < 0.2 and p500_max > 0.2:
            print(f"  P(tract≥500bp): {p500_max:.1%} - {p500_min:.1%} "
                  f"→ SENSITIVE (crosses 20% threshold)")
        else:
            print(f"  P(tract≥500bp): {p500_max:.1%} - {p500_min:.1%} "
                  f"→ consistently low")


def generate_figure(results):
    """Generate publication-quality sensitivity analysis figure."""
    if not HAS_MPL:
        return

    fig = plt.figure(figsize=(16, 10), facecolor="#1B2838")
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3,
                  left=0.08, right=0.95, top=0.92, bottom=0.08)

    # Color theme
    TEXT_COLOR = "#E8EEF4"
    GRID_COLOR = "#2A3A4A"
    BG_COLOR = "#1B2838"
    CARD_COLOR = "#233448"
    DEFAULT_LINE_COLOR = "#FF4444"

    def style_ax(ax, title, xlabel, ylabel):
        ax.set_facecolor(CARD_COLOR)
        ax.set_title(title, color=TEXT_COLOR, fontsize=13, fontweight="bold", pad=10)
        ax.set_xlabel(xlabel, color=TEXT_COLOR, fontsize=11)
        ax.set_ylabel(ylabel, color=TEXT_COLOR, fontsize=11)
        ax.tick_params(colors=TEXT_COLOR, which="both")
        ax.grid(True, alpha=0.3, color=GRID_COLOR)
        for spine in ax.spines.values():
            spine.set_color(GRID_COLOR)

    # ── Panel A: Mean tract length vs p ──
    ax_a = fig.add_subplot(gs[0, 0])
    for cfg in CONFIGS:
        name = cfg["name"]
        p_vals = sorted(results[name].keys())
        means = [results[name][p]["mean_tract"] for p in p_vals]
        se_vals = [results[name][p]["se_mean"] for p in p_vals]
        ax_a.errorbar(
            [p * 1000 for p in p_vals], means,
            yerr=[1.96 * se for se in se_vals],
            label=name, color=cfg["color"], linestyle=cfg["linestyle"],
            marker="o", markersize=5, capsize=3, linewidth=2,
        )
    ax_a.axvline(x=2.0, color=DEFAULT_LINE_COLOR, linestyle=":", alpha=0.7,
                 label="Default (p=0.002)")
    style_ax(ax_a, "A. Mean Tract Length vs. p",
             "p (× 10⁻³ per bp)", "Mean tract length (bp)")
    ax_a.legend(fontsize=8, facecolor=CARD_COLOR, edgecolor=GRID_COLOR,
                labelcolor=TEXT_COLOR, loc="upper right")

    # ── Panel B: P(tract ≥ d) vs p for key distances ──
    ax_b = fig.add_subplot(gs[0, 1])
    # Use the enFnCas9 + cssDNA (iPSC) config as reference
    ref_cfg = "enFnCas9 + cssDNA (iPSC)"
    dist_colors = {100: "#2ECC71", 200: "#1B9E77", 500: "#3498DB",
                   1000: "#D95F02", 2000: "#E74C3C"}
    for d in [200, 500, 1000, 2000]:
        p_vals = sorted(results[ref_cfg].keys())
        probs = [results[ref_cfg][p]["dist_probs"][d]["prob"] for p in p_vals]
        ses = [results[ref_cfg][p]["dist_probs"][d]["se"] for p in p_vals]
        ax_b.errorbar(
            [p * 1000 for p in p_vals], [pr * 100 for pr in probs],
            yerr=[1.96 * se * 100 for se in ses],
            label=f"≥{d} bp", color=dist_colors[d],
            marker="s", markersize=4, capsize=3, linewidth=2,
        )
    ax_b.axvline(x=2.0, color=DEFAULT_LINE_COLOR, linestyle=":", alpha=0.7)
    style_ax(ax_b, f"B. Distance Probabilities ({ref_cfg})",
             "p (× 10⁻³ per bp)", "P(tract ≥ d) (%)")
    ax_b.legend(fontsize=9, facecolor=CARD_COLOR, edgecolor=GRID_COLOR,
                labelcolor=TEXT_COLOR, loc="upper right")

    # ── Panel C: HDR rate vs p ──
    ax_c = fig.add_subplot(gs[1, 0])
    for cfg in CONFIGS:
        name = cfg["name"]
        p_vals = sorted(results[name].keys())
        rates = [results[name][p]["hdr_rate"] * 100 for p in p_vals]
        ax_c.plot(
            [p * 1000 for p in p_vals], rates,
            label=name, color=cfg["color"], linestyle=cfg["linestyle"],
            marker="D", markersize=5, linewidth=2,
        )
    ax_c.axvline(x=2.0, color=DEFAULT_LINE_COLOR, linestyle=":", alpha=0.7,
                 label="Default (p=0.002)")
    style_ax(ax_c, "C. HDR Rate vs. p", "p (× 10⁻³ per bp)", "HDR rate (%)")
    ax_c.legend(fontsize=8, facecolor=CARD_COLOR, edgecolor=GRID_COLOR,
                labelcolor=TEXT_COLOR, loc="upper right")

    # ── Panel D: Summary text ──
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_facecolor(CARD_COLOR)
    ax_d.axis("off")

    # Compute fold-change for default config
    ref = results[ref_cfg]
    default_mean = ref[0.002]["mean_tract"]
    min_mean = ref[P_VALUES[-1]]["mean_tract"]
    max_mean = ref[P_VALUES[0]]["mean_tract"]
    default_hdr = ref[0.002]["hdr_rate"]
    min_hdr = ref[P_VALUES[-1]]["hdr_rate"]
    max_hdr = ref[P_VALUES[0]]["hdr_rate"]
    default_p500 = ref[0.002]["dist_probs"][500]["prob"]

    summary_text = (
        f"ROBUSTNESS SUMMARY\n"
        f"{'─' * 35}\n\n"
        f"Parameter: SDSA displacement prob.\n"
        f"Range:  {P_VALUES[0]} – {P_VALUES[-1]}\n"
        f"        ({P_VALUES[-1]/P_VALUES[0]:.0f}-fold range)\n\n"
        f"Default (p=0.002):\n"
        f"  Mean tract: {default_mean:.0f} bp\n"
        f"  HDR rate:   {default_hdr:.1%}\n"
        f"  P(≥500bp):  {default_p500:.1%}\n\n"
        f"Across full range:\n"
        f"  Tract: {min_mean:.0f} – {max_mean:.0f} bp\n"
        f"  HDR:   {min(min_hdr, max_hdr):.1%} – {max(min_hdr, max_hdr):.1%}\n\n"
    )

    # Qualitative conclusion
    if min_hdr > 0.005 and max_hdr > 0.005:
        summary_text += "Conclusion: HDR remains\nFEASIBLE across full range.\nRecommendations are ROBUST."
        conclusion_color = "#2ECC71"
    else:
        summary_text += "Conclusion: HDR feasibility\nis SENSITIVE to p.\nExperimental calibration needed."
        conclusion_color = "#FF9800"

    ax_d.text(0.05, 0.95, summary_text, transform=ax_d.transAxes,
              fontsize=11, color=TEXT_COLOR, fontfamily="monospace",
              verticalalignment="top", linespacing=1.4)

    for spine in ax_d.spines.values():
        spine.set_color(GRID_COLOR)

    # Title
    fig.suptitle(
        "SDSA Displacement Probability Sensitivity Analysis",
        fontsize=16, fontweight="bold", color=TEXT_COLOR, y=0.98,
    )
    fig.text(0.5, 0.94,
             f"p swept from {P_VALUES[0]} to {P_VALUES[-1]} "
             f"({len(P_VALUES)} values × {N_SIMULATIONS:,} simulations × "
             f"{len(CONFIGS)} configurations)",
             fontsize=10, color="#889AAA", ha="center")

    # Save
    fig_dir = os.path.join(os.path.dirname(__file__), "figures")
    os.makedirs(fig_dir, exist_ok=True)

    for fmt in ["png", "pdf"]:
        path = os.path.join(fig_dir, f"Fig_SDSA_Sensitivity.{fmt}")
        fig.savefig(path, dpi=200, facecolor=fig.get_facecolor(),
                    edgecolor="none")
        print(f"Saved: {path}")

    plt.close(fig)


# =========================================================================
# Main
# =========================================================================
if __name__ == "__main__":
    print("Running SDSA sensitivity analysis...")
    print(f"  p values: {P_VALUES}")
    print(f"  Configurations: {len(CONFIGS)}")
    print(f"  Simulations per point: {N_SIMULATIONS:,}")
    print(f"  Total simulations: {len(P_VALUES) * len(CONFIGS) * N_SIMULATIONS:,}")
    print()

    results = run_sensitivity_sweep()
    print_summary_table(results)
    generate_figure(results)

    print("\nDone.")
