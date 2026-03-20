"""
Visualization Functions for Loop Extrusion Simulations
=======================================================

This module provides publication-quality visualization of cohesin loop
extrusion simulation results. Three main plot types are available:

1. **Extrusion Kymograph** — A space-time diagram showing how the loop
   expands over time. This is the most informative single visualization
   of the extrusion process.

2. **Scan Probability Track** — A genomic track showing the per-position
   probability of being scanned by RAD51. This mimics the format of
   ChIP-seq enrichment profiles and is the most directly useful output
   for experimental biologists.

3. **Dual-DSB Comparison** — Side-by-side or overlaid search domains
   from two DSBs, demonstrating (non-)overlap.

All functions use matplotlib and are designed for:
- High-resolution figure export (300 DPI default)
- Colorblind-friendly palettes
- Clear axis labels with genomic coordinates
- Biological annotations (CTCF sites, DSB positions, chromatin states)

Kymograph Biology Note
----------------------
A "kymograph" in biology is a space-time plot. In the context of loop
extrusion, it shows:
- X-axis: genomic position (Mb)
- Y-axis: time (seconds or minutes)
- Colored region: the extruded loop domain

The domain starts as a single point (the DSB) and expands outward as
cohesin reels chromatin from both sides. When an extrusion arm encounters
a CTCF site and stalls, the boundary becomes a flat horizontal line.
The final shape resembles a "V" or funnel opening upward (or downward,
depending on time axis direction).

Multiple overlaid traces (from different Monte Carlo runs) show the
stochasticity: some traces stall early (strong CTCF encounter), others
extend much further (CTCF bypass events).
"""

from __future__ import annotations

from typing import List, Optional, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .chromatin_fiber import ChromatinFiber
    from .cohesin_extruder import ExtrusionResult
    from .simulator import LoopSimResults


def plot_extrusion_kymograph(
    trace_results: List["ExtrusionResult"],
    fiber: "ChromatinFiber",
    figsize: tuple = (12, 6),
    save_path: Optional[str] = None,
    time_unit: str = "minutes",
):
    """Plot a kymograph (space-time diagram) of loop extrusion.

    Each trace shows one stochastic extrusion simulation as a pair of
    lines expanding outward from the DSB anchor point. The left arm
    extends toward lower genomic coordinates; the right arm toward higher.
    CTCF stall events are visible as horizontal plateaus where the
    boundary stops advancing.

    This plot directly illustrates the core mechanism of Marin-Gonzalez
    et al. (2025): cohesin loaded at a DSB creates an expanding chromatin
    loop whose boundaries are determined by CTCF barriers.

    Parameters
    ----------
    trace_results : list of ExtrusionResult
        Individual extrusion results WITH time traces recorded
        (``time_trace_left`` and ``time_trace_right`` must not be None).
    fiber : ChromatinFiber
        The chromatin fiber (for coordinate conversion and CTCF annotation).
    figsize : tuple, optional
        Figure size (width, height) in inches.
    save_path : str, optional
        If provided, save figure to this path.
    time_unit : str, optional
        "seconds" or "minutes". Default "minutes".

    Returns
    -------
    matplotlib.figure.Figure
        The kymograph figure.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Color palette for traces (colorblind-friendly)
    colors = [
        "#2196F3",  # blue
        "#FF9800",  # orange
        "#4CAF50",  # green
        "#9C27B0",  # purple
        "#F44336",  # red
        "#00BCD4",  # cyan
        "#795548",  # brown
        "#607D8B",  # blue-grey
    ]

    time_divisor = 60.0 if time_unit == "minutes" else 1.0
    time_label = "Time (minutes)" if time_unit == "minutes" else "Time (seconds)"

    dsb_position_mb = None

    for i, result in enumerate(trace_results):
        if result.time_trace_left is None or result.time_trace_right is None:
            continue

        color = colors[i % len(colors)]
        n_steps = len(result.time_trace_left)
        time_axis = np.arange(n_steps) / time_divisor

        # Convert bead indices to Mb coordinates
        left_mb = (
            fiber.start_bp + result.time_trace_left * fiber.resolution_bp
        ) / 1e6
        right_mb = (
            fiber.start_bp + result.time_trace_right * fiber.resolution_bp
        ) / 1e6

        dsb_mb = result.dsb_position_bp / 1e6
        dsb_position_mb = dsb_mb

        # Plot left and right boundaries
        ax.plot(left_mb, time_axis, color=color, linewidth=1.2, alpha=0.8)
        ax.plot(right_mb, time_axis, color=color, linewidth=1.2, alpha=0.8)

        # Fill the extruded region with transparent color
        ax.fill_betweenx(
            time_axis, left_mb, right_mb,
            color=color, alpha=0.08,
        )

        # Mark stall events
        if result.left_stalled_at_ctcf:
            stall_time = time_axis[-1]
            stall_pos = left_mb[-1]
            ax.plot(
                stall_pos, stall_time, "v",
                color=color, markersize=6, zorder=5,
            )
        if result.right_stalled_at_ctcf:
            stall_time = time_axis[-1]
            stall_pos = right_mb[-1]
            ax.plot(
                stall_pos, stall_time, "v",
                color=color, markersize=6, zorder=5,
            )

    # Mark DSB position
    if dsb_position_mb is not None:
        ax.axvline(
            dsb_position_mb, color="#C0392B", linewidth=2,
            linestyle="--", zorder=10, label="DSB position",
        )

    # Mark CTCF sites
    for site in fiber.ctcf_sites:
        site_mb = site.position_bp / 1e6
        if site.orientation == "forward":
            marker = ">"
            color_ctcf = "#E74C3C"
        else:
            marker = "<"
            color_ctcf = "#3498DB"
        ax.plot(
            site_mb, 0, marker, color=color_ctcf, markersize=5,
            alpha=0.6, zorder=4,
        )

    ax.set_xlabel("Genomic Position (Mb)", fontsize=12)
    ax.set_ylabel(time_label, fontsize=12)
    ax.set_title(
        "Loop Extrusion Kymograph\n"
        "Cohesin-mediated chromatin scanning from DSB anchor",
        fontsize=13, fontweight="bold",
    )
    ax.legend(fontsize=10, loc="upper right", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Invert y-axis so time flows downward (convention for kymographs)
    ax.invert_yaxis()

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_scan_probability_track(
    results: "LoopSimResults",
    fiber: "ChromatinFiber",
    figsize: tuple = (14, 6),
    save_path: Optional[str] = None,
):
    """Plot the scan probability as a genomic track with CTCF annotations.

    This is the primary visualization output of LoopSim. It shows a
    genome-browser-style track where:

    - **Top panel**: The per-position probability that a given bead is
      within the cohesin-extruded search domain. Values near 1.0 (dark
      blue) indicate positions always scanned; values near 0 (white)
      indicate positions rarely or never reached.

    - **Bottom annotation**: CTCF sites shown as triangles. Forward-
      oriented sites (blocking rightward extrusion) point right (>).
      Reverse-oriented sites (blocking leftward extrusion) point left (<).

    - **DSB position(s)**: Marked with red vertical lines.

    For dual-DSB simulations, both search domains are overlaid in
    different colors, making it easy to see whether they overlap.

    Parameters
    ----------
    results : LoopSimResults
        Simulation results from LoopSimulator.run().
    fiber : ChromatinFiber
        The chromatin fiber model.
    figsize : tuple, optional
        Figure size. Default (14, 6).
    save_path : str, optional
        If provided, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
        The genomic track figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.gridspec as gridspec

    n_dsbs = len(results.dsb_positions)

    # Use GridSpec: top = scan probability, bottom = CTCF annotations
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        2, 1, height_ratios=[4, 1], hspace=0.08, figure=fig,
    )
    ax_main = fig.add_subplot(gs[0])
    ax_ctcf = fig.add_subplot(gs[1], sharex=ax_main)

    # Genomic position axis (Mb)
    positions_mb = (
        fiber.start_bp + np.arange(fiber.n_beads) * fiber.resolution_bp
    ) / 1e6

    # Domain colors (colorblind-friendly)
    domain_colors = ["#2171B5", "#CB181D"]  # blue, red
    domain_fills = ["#C6DBEF", "#FCBBA1"]  # light blue, light red
    domain_labels = []

    for i, dsb_pos in enumerate(results.dsb_positions):
        domain = results.per_dsb_results[dsb_pos]
        color = domain_colors[i % len(domain_colors)]
        fill = domain_fills[i % len(domain_fills)]

        prob = domain.position_scan_probability

        # Fill area under the probability curve
        ax_main.fill_between(
            positions_mb, 0, prob,
            color=fill, alpha=0.6,
        )
        # Draw the probability curve
        ax_main.plot(
            positions_mb, prob,
            color=color, linewidth=1.5, alpha=0.9,
            label=f"DSB at {dsb_pos / 1e6:.2f} Mb "
                  f"(domain ~{domain.mean_domain_size_bp / 1e6:.1f} Mb)",
        )

        # Mark DSB position
        ax_main.axvline(
            dsb_pos / 1e6, color=color, linewidth=2,
            linestyle="--", alpha=0.8,
        )
        # Add DSB label at top
        ax_main.text(
            dsb_pos / 1e6, 1.05, "DSB",
            ha="center", va="bottom", fontsize=9,
            color=color, fontweight="bold",
            transform=ax_main.get_xaxis_transform(),
        )

    ax_main.set_ylabel("Scan Probability", fontsize=12)
    ax_main.set_ylim(0, 1.05)
    ax_main.set_xlim(positions_mb[0], positions_mb[-1])
    ax_main.legend(fontsize=10, loc="upper right", framealpha=0.9)
    ax_main.set_title(
        "RAD51 Homology Search Domain — Cohesin Loop Extrusion Model\n"
        "(Marin-Gonzalez et al., Science, 2025)",
        fontsize=13, fontweight="bold",
    )
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)
    plt.setp(ax_main.get_xticklabels(), visible=False)

    # Mark heterochromatin regions
    chrom_state = fiber.chromatin_state
    closed_beads = np.where(chrom_state == 1)[0]
    if len(closed_beads) > 0:
        # Find contiguous closed regions
        diffs = np.diff(closed_beads)
        breaks = np.where(diffs > 1)[0]
        starts = np.concatenate([[0], breaks + 1])
        ends = np.concatenate([breaks, [len(closed_beads) - 1]])
        for s, e in zip(starts, ends):
            x_start = positions_mb[closed_beads[s]]
            x_end = positions_mb[closed_beads[e]]
            ax_main.axvspan(
                x_start, x_end,
                color="#BDBDBD", alpha=0.2, zorder=0,
            )

    # ---- CTCF annotation panel ----
    ax_ctcf.set_ylim(-0.5, 0.5)
    ax_ctcf.set_yticks([])
    ax_ctcf.axhline(0, color="#CCCCCC", linewidth=0.5)

    for site in fiber.ctcf_sites:
        site_mb = site.position_bp / 1e6
        if site.orientation == "forward":
            # Forward CTCF: triangle pointing right (blocks rightward extrusion)
            ax_ctcf.plot(
                site_mb, 0.15, ">",
                color="#E74C3C", markersize=max(3, site.strength * 8),
                alpha=0.7,
            )
        else:
            # Reverse CTCF: triangle pointing left (blocks leftward extrusion)
            ax_ctcf.plot(
                site_mb, -0.15, "<",
                color="#3498DB", markersize=max(3, site.strength * 8),
                alpha=0.7,
            )

    # DSB markers in CTCF panel
    for dsb_pos in results.dsb_positions:
        ax_ctcf.axvline(
            dsb_pos / 1e6, color="#C0392B", linewidth=2,
            linestyle="--", alpha=0.5,
        )

    ax_ctcf.set_xlabel("Genomic Position (Mb)", fontsize=12)
    ax_ctcf.set_ylabel("CTCF", fontsize=10)
    ax_ctcf.spines["top"].set_visible(False)
    ax_ctcf.spines["right"].set_visible(False)

    # Add CTCF legend
    ax_ctcf.plot([], [], ">", color="#E74C3C", markersize=6, label="Forward CTCF (blocks right)")
    ax_ctcf.plot([], [], "<", color="#3498DB", markersize=6, label="Reverse CTCF (blocks left)")
    ax_ctcf.legend(fontsize=8, loc="lower right", ncol=2, framealpha=0.9)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_dual_dsb_comparison(
    results1,
    results2,
    fiber: "ChromatinFiber",
    distance_bp: Optional[int] = None,
    figsize: tuple = (14, 8),
    save_path: Optional[str] = None,
):
    """Plot overlaid search domains from two DSBs to show (non-)overlap.

    This visualization is specifically designed to demonstrate the key
    prediction from the loop extrusion model: **two DSBs separated by more
    than one TAD have non-overlapping search domains**. This means:

    - Each DSB undergoes independent repair.
    - The intervening chromatin segment is at risk of deletion.
    - Exogenous donors must be targeted to each DSB independently.

    The plot shows three panels:
    1. **Top**: DSB1 search domain (blue)
    2. **Middle**: DSB2 search domain (red)
    3. **Bottom**: Overlay of both domains, highlighting overlap (if any)

    Parameters
    ----------
    results1 : SearchDomainResult
        Search domain result for DSB1.
    results2 : SearchDomainResult
        Search domain result for DSB2.
    fiber : ChromatinFiber
        The chromatin fiber model.
    distance_bp : int, optional
        Distance between the DSBs (for annotation). Computed if not given.
    figsize : tuple, optional
        Figure size. Default (14, 8).
    save_path : str, optional
        If provided, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
        The comparison figure.
    """
    import matplotlib.pyplot as plt

    if distance_bp is None:
        distance_bp = abs(results2.dsb_position_bp - results1.dsb_position_bp)
    distance_mb = distance_bp / 1e6

    positions_mb = (
        fiber.start_bp + np.arange(fiber.n_beads) * fiber.resolution_bp
    ) / 1e6

    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)

    prob1 = results1.position_scan_probability
    prob2 = results2.position_scan_probability

    # Panel 1: DSB1
    axes[0].fill_between(positions_mb, 0, prob1, color="#C6DBEF", alpha=0.7)
    axes[0].plot(positions_mb, prob1, color="#2171B5", linewidth=1.5)
    axes[0].axvline(
        results1.dsb_position_bp / 1e6, color="#2171B5",
        linewidth=2, linestyle="--",
    )
    axes[0].set_ylabel("P(scanned)", fontsize=11)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title(
        f"DSB 1 at {results1.dsb_position_bp / 1e6:.2f} Mb "
        f"(domain ~{results1.mean_domain_size_bp / 1e6:.1f} Mb)",
        fontsize=11, fontweight="bold", color="#2171B5",
    )
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # Panel 2: DSB2
    axes[1].fill_between(positions_mb, 0, prob2, color="#FCBBA1", alpha=0.7)
    axes[1].plot(positions_mb, prob2, color="#CB181D", linewidth=1.5)
    axes[1].axvline(
        results2.dsb_position_bp / 1e6, color="#CB181D",
        linewidth=2, linestyle="--",
    )
    axes[1].set_ylabel("P(scanned)", fontsize=11)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title(
        f"DSB 2 at {results2.dsb_position_bp / 1e6:.2f} Mb "
        f"(domain ~{results2.mean_domain_size_bp / 1e6:.1f} Mb)",
        fontsize=11, fontweight="bold", color="#CB181D",
    )
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    # Panel 3: Overlay
    axes[2].fill_between(positions_mb, 0, prob1, color="#C6DBEF", alpha=0.5, label="DSB1 domain")
    axes[2].fill_between(positions_mb, 0, prob2, color="#FCBBA1", alpha=0.5, label="DSB2 domain")
    axes[2].plot(positions_mb, prob1, color="#2171B5", linewidth=1.2, alpha=0.8)
    axes[2].plot(positions_mb, prob2, color="#CB181D", linewidth=1.2, alpha=0.8)

    # Highlight overlap region
    overlap = np.minimum(prob1, prob2)
    overlap_mask = overlap > 0.01
    if np.any(overlap_mask):
        axes[2].fill_between(
            positions_mb, 0, overlap,
            where=overlap_mask,
            color="#8E44AD", alpha=0.4,
            label="Overlap region",
        )

    axes[2].axvline(
        results1.dsb_position_bp / 1e6, color="#2171B5",
        linewidth=2, linestyle="--", alpha=0.7,
    )
    axes[2].axvline(
        results2.dsb_position_bp / 1e6, color="#CB181D",
        linewidth=2, linestyle="--", alpha=0.7,
    )
    axes[2].set_ylabel("P(scanned)", fontsize=11)
    axes[2].set_ylim(0, 1.05)
    axes[2].set_xlabel("Genomic Position (Mb)", fontsize=12)
    axes[2].set_title(
        f"Overlay — DSBs separated by {distance_mb:.2f} Mb",
        fontsize=11, fontweight="bold",
    )
    axes[2].legend(fontsize=9, loc="upper right", framealpha=0.9)
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)

    # Mark CTCF sites in the overlay panel
    for site in fiber.ctcf_sites:
        site_mb = site.position_bp / 1e6
        marker = ">" if site.orientation == "forward" else "<"
        color_ctcf = "#E74C3C" if site.orientation == "forward" else "#3498DB"
        axes[2].plot(
            site_mb, -0.03, marker, color=color_ctcf,
            markersize=4, alpha=0.5, clip_on=False,
        )

    # Add overall title with overlap assessment
    max_overlap = float(np.max(overlap)) if np.any(overlap_mask) else 0.0
    if max_overlap < 0.01:
        overlap_text = "NO OVERLAP — independent repair domains"
    elif max_overlap < 0.1:
        overlap_text = "MINIMAL OVERLAP — mostly independent"
    else:
        overlap_text = f"OVERLAP detected (max P = {max_overlap:.2f})"

    fig.suptitle(
        f"Dual-DSB Search Domain Comparison\n{overlap_text}",
        fontsize=14, fontweight="bold", y=1.02,
    )

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig
