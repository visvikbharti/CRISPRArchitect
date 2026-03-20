"""
Visualization Utilities for CRISPRArchitect
============================================

Plotting functions shared across all modules. Uses matplotlib and seaborn
for publication-quality figures.

All plots include biological annotations and explanations to help
researchers interpret the results.
"""

import numpy as np
from typing import List, Optional, Tuple, Dict

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.ticker import FuncFormatter
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


def check_matplotlib():
    """Check if matplotlib is available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )


def setup_style():
    """Set up a clean, publication-ready plot style."""
    check_matplotlib()
    if HAS_SEABORN:
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=1.2)
    else:
        plt.style.use('seaborn-v0_8-whitegrid')


def plot_tract_distribution(
    tract_lengths: np.ndarray,
    title: str = "Gene Conversion Tract Length Distribution",
    key_distances: Optional[List[int]] = None,
    save_path: Optional[str] = None,
) -> None:
    """Plot the distribution of gene conversion tract lengths from simulation.

    This visualization shows how far from the DSB cut site the donor template
    sequence is expected to be incorporated during HDR. It is the key output
    of the ConversionSim module.

    Parameters
    ----------
    tract_lengths : np.ndarray
        Array of simulated tract lengths in base pairs
    title : str
        Plot title
    key_distances : list of int, optional
        Specific distances (bp) to annotate on the plot
        (e.g., [200, 500, 1000, 2000] to show incorporation probability at these distances)
    save_path : str, optional
        If provided, save the figure to this path
    """
    check_matplotlib()
    setup_style()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left panel: Histogram ---
    ax1.hist(tract_lengths, bins=50, color='#2196F3', alpha=0.7, edgecolor='white')
    ax1.axvline(np.median(tract_lengths), color='red', linestyle='--', linewidth=2,
                label=f'Median: {np.median(tract_lengths):.0f} bp')
    ax1.axvline(np.mean(tract_lengths), color='orange', linestyle='--', linewidth=2,
                label=f'Mean: {np.mean(tract_lengths):.0f} bp')
    ax1.set_xlabel('Gene Conversion Tract Length (bp)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(title, fontsize=13)
    ax1.legend(fontsize=10)

    # --- Right panel: Survival curve (probability of reaching distance X) ---
    sorted_tracts = np.sort(tract_lengths)
    survival = 1.0 - np.arange(1, len(sorted_tracts) + 1) / len(sorted_tracts)
    ax2.plot(sorted_tracts, survival, color='#4CAF50', linewidth=2)
    ax2.set_xlabel('Distance from DSB (bp)', fontsize=12)
    ax2.set_ylabel('P(tract extends at least this far)', fontsize=12)
    ax2.set_title('Probability of Edit Incorporation vs. Distance', fontsize=13)
    ax2.set_ylim(-0.05, 1.05)

    # Annotate key distances
    if key_distances is None:
        key_distances = [200, 500, 1000, 2000, 5000]

    for dist in key_distances:
        prob = np.mean(tract_lengths >= dist)
        if prob > 0.001:  # Only annotate if probability is meaningful
            ax2.axvline(dist, color='gray', linestyle=':', alpha=0.5)
            ax2.annotate(f'{dist} bp: {prob:.1%}',
                        xy=(dist, prob), xytext=(dist + 100, prob + 0.05),
                        fontsize=9, arrowprops=dict(arrowstyle='->', color='gray'))

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_distance_comparison(
    genomic_distances_bp: List[int],
    physical_distances_nm: List[float],
    donor_diameters_nm: List[float],
    donor_labels: List[str],
    save_path: Optional[str] = None,
) -> None:
    """Plot genomic vs. physical distance with donor size comparison.

    This is the key visualization for ChromBridge, showing why donor templates
    are too small to bridge distant genomic loci.

    Parameters
    ----------
    genomic_distances_bp : list of int
        Genomic separations to plot (e.g., [1e4, 1e5, 1e6, 1e7])
    physical_distances_nm : list of float
        Predicted 3D distances for each genomic separation
    donor_diameters_nm : list of float
        Random coil diameters for different donor sizes
    donor_labels : list of str
        Labels for each donor (e.g., ["3 kb cssDNA", "10 kb cssDNA", "20 kb cssDNA"])
    save_path : str, optional
        Path to save the figure
    """
    check_matplotlib()
    setup_style()

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot physical distance vs genomic distance
    ax.plot(genomic_distances_bp, physical_distances_nm,
            'o-', color='#E91E63', linewidth=2, markersize=8,
            label='3D inter-locus distance')

    # Plot donor diameters as horizontal bands
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']
    for i, (diameter, label) in enumerate(zip(donor_diameters_nm, donor_labels)):
        color = colors[i % len(colors)]
        ax.axhline(diameter, color=color, linestyle='--', linewidth=2,
                   label=f'{label} (diameter: {diameter:.0f} nm)', alpha=0.7)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Genomic Distance (bp)', fontsize=12)
    ax.set_ylabel('Distance (nm)', fontsize=12)
    ax.set_title('3D Inter-Locus Distance vs. Donor Template Size\n'
                 '(Donor must be larger than inter-locus distance to bridge)',
                 fontsize=13)
    ax.legend(fontsize=10, loc='upper left')

    # Format x-axis with readable labels
    def bp_formatter(x, pos):
        if x >= 1e6:
            return f'{x/1e6:.0f} Mb'
        elif x >= 1e3:
            return f'{x/1e3:.0f} kb'
        return f'{x:.0f} bp'

    ax.xaxis.set_major_formatter(FuncFormatter(bp_formatter))

    # Add annotation
    ax.annotate('Donor CANNOT bridge\n(too small)',
               xy=(0.7, 0.3), xycoords='axes fraction',
               fontsize=14, color='red', fontweight='bold',
               ha='center', va='center',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_gene_structure(
    exons: List[Dict],
    mutations: Optional[List[Dict]] = None,
    highlight_region: Optional[Tuple[int, int]] = None,
    title: str = "Gene Structure",
    save_path: Optional[str] = None,
) -> None:
    """Plot a schematic gene structure diagram with exons, introns, and mutations.

    Parameters
    ----------
    exons : list of dict
        Each dict has: 'number', 'start', 'end', 'size'
    mutations : list of dict, optional
        Each dict has: 'exon_number', 'position', 'label'
    highlight_region : tuple, optional
        (start_exon, end_exon) to highlight a region
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    check_matplotlib()
    setup_style()

    fig, ax = plt.subplots(figsize=(16, 3))

    gene_start = exons[0]['start']
    gene_end = exons[-1]['end']
    gene_length = gene_end - gene_start

    # Draw intron line (thin)
    ax.plot([gene_start, gene_end], [0, 0], color='gray', linewidth=1)

    # Draw exons (thick boxes)
    for exon in exons:
        width = exon['end'] - exon['start']
        color = '#2196F3'
        if highlight_region and highlight_region[0] <= exon['number'] <= highlight_region[1]:
            color = '#FF9800'
        rect = mpatches.FancyBboxPatch(
            (exon['start'], -0.3), width, 0.6,
            boxstyle="round,pad=0.01", facecolor=color,
            edgecolor='black', linewidth=0.5, alpha=0.8
        )
        ax.add_patch(rect)

        # Label every 10th exon
        if exon['number'] % 10 == 0 or exon['number'] == 1:
            ax.text(exon['start'] + width / 2, -0.6, f"E{exon['number']}",
                   ha='center', va='top', fontsize=7)

    # Mark mutations
    if mutations:
        for mut in mutations:
            ax.plot(mut['position'], 0.5, 'v', color='red', markersize=12)
            ax.text(mut['position'], 0.8, mut.get('label', ''),
                   ha='center', fontsize=9, color='red', fontweight='bold')

    ax.set_xlim(gene_start - gene_length * 0.02, gene_end + gene_length * 0.02)
    ax.set_ylim(-1.0, 1.5)
    ax.set_xlabel('Genomic Position (bp)', fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.set_yticks([])

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_accessibility_map(
    accessibility_scores: np.ndarray,
    sequence_length: int,
    left_arm: Tuple[int, int],
    right_arm: Tuple[int, int],
    g4_regions: Optional[List[Tuple[int, int]]] = None,
    hairpin_regions: Optional[List[Tuple[int, int]]] = None,
    title: str = "cssDNA Donor Accessibility Map",
    save_path: Optional[str] = None,
) -> None:
    """Plot per-nucleotide accessibility of a cssDNA donor template.

    Shows which regions of the donor are accessible for RAD51 binding
    vs. sequestered in secondary structures.

    Parameters
    ----------
    accessibility_scores : np.ndarray
        Per-nucleotide scores (0 = structured, 1 = accessible)
    sequence_length : int
        Total donor length
    left_arm : tuple
        (start, end) of left homology arm
    right_arm : tuple
        (start, end) of right homology arm
    g4_regions : list of tuple, optional
        G-quadruplex positions
    hairpin_regions : list of tuple, optional
        Hairpin positions
    save_path : str, optional
        Path to save figure
    """
    check_matplotlib()
    setup_style()

    fig, ax = plt.subplots(figsize=(14, 4))

    positions = np.arange(sequence_length)

    # Color by accessibility: green = accessible, red = structured
    colors = np.where(accessibility_scores > 0.5, '#4CAF50', '#F44336')

    ax.bar(positions, accessibility_scores, width=1.0, color=colors, alpha=0.6)

    # Highlight homology arms
    ax.axvspan(left_arm[0], left_arm[1], alpha=0.15, color='blue',
               label='Left homology arm')
    ax.axvspan(right_arm[0], right_arm[1], alpha=0.15, color='purple',
               label='Right homology arm')

    # Mark G4 regions
    if g4_regions:
        for start, end in g4_regions:
            ax.axvspan(start, end, alpha=0.3, color='red', label='G-quadruplex')

    # Mark hairpin regions
    if hairpin_regions:
        for start, end in hairpin_regions:
            ax.axvspan(start, end, alpha=0.2, color='orange', label='Hairpin')

    ax.set_xlabel('Position in Donor (nt)', fontsize=12)
    ax.set_ylabel('Accessibility Score', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.set_ylim(0, 1.1)

    # Remove duplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=10)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def format_bp(value: float) -> str:
    """Format a base pair count in human-readable form.

    Examples: 1500 -> '1.5 kb', 2500000 -> '2.5 Mb'
    """
    if abs(value) >= 1e6:
        return f"{value / 1e6:.1f} Mb"
    elif abs(value) >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{value:.0f} bp"


def format_nm(value: float) -> str:
    """Format a nanometer distance in human-readable form.

    Examples: 500 -> '500 nm', 1500 -> '1.5 μm'
    """
    if abs(value) >= 1000:
        return f"{value / 1000:.1f} μm"
    return f"{value:.0f} nm"
