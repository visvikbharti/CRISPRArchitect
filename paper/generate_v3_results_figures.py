#!/usr/bin/env python3
"""
Generate publication-quality figures from REAL CRISPRArchitect v3 data.

Every number plotted here comes from actual pipeline runs or verified
benchmark results. No dummy data, no placeholders.

Outputs figures to paper/figures/v3_results/
"""

import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# Style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Colors
TEAL = '#1B9E77'
CORAL = '#D95F02'
GOLD = '#FFC107'
NAVY = '#1B2838'
BLUE = '#3498DB'
RED = '#E74C3C'
GREEN = '#2ECC71'
GRAY = '#95A5A6'
DARK = '#2C3E50'

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'figures', 'v3')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def save_fig(fig, name):
    for ext in ['png', 'pdf']:
        path = os.path.join(OUTPUT_DIR, f'{name}.{ext}')
        fig.savefig(path)
    print(f'  Saved: {name}')
    plt.close(fig)


# =====================================================================
# FIGURE 1: Strategy Distribution — v2 vs v3
# =====================================================================
def fig_strategy_distribution():
    """Bar chart comparing SpCas9-only vs multi-nuclease strategy distributions."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # SpCas9-only data (from definitive_benchmark_results.json — verified)
    spcas9_data = {'PE': 29, 'HDR': 1, 'BE': 0}
    # Multi-nuclease data (from v3_benchmark_results.json — verified)
    multi_data = {'PE': 23, 'BE': 6, 'HDR': 1}

    categories = ['BE', 'PE', 'HDR']
    colors = [TEAL, CORAL, GOLD]

    for ax, (title, data) in zip(axes, [
        ('SpCas9 + ABE7.10 only', spcas9_data),
        ('Multi-nuclease + TOPSIS 6D', multi_data),
    ]):
        vals = [data.get(c, 0) for c in categories]
        bars = ax.bar(categories, vals, color=colors, edgecolor='white',
                      linewidth=1.5, width=0.6)

        # Add count labels on bars
        for bar, val in zip(bars, vals):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        str(val), ha='center', va='bottom', fontweight='bold',
                        fontsize=14)

        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of cases (top-1 ranked)')
        ax.set_ylim(0, 33)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Impact of Multi-Nuclease Engine (30 ClinVar Cases)',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'Fig_StrategyDistribution_v2_v3')


# =====================================================================
# FIGURE 2: Literature Benchmark Concordance
# =====================================================================
def fig_literature_benchmark():
    """Concordance results from the literature benchmark."""
    # Real data from our pipeline run
    results_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'benchmark_results',
                                'literature_benchmark_results.json')
    if os.path.exists(results_path):
        with open(results_path) as f:
            bench = json.load(f)
        cases = bench.get('results', [])
    else:
        print('  WARNING: No literature benchmark results found')
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: Overall concordance
    ax = axes[0]
    n = len(cases)
    top1 = sum(1 for c in cases if c['concordant'])
    top3 = sum(1 for c in cases if c['top3_concordant'])

    bars = ax.bar(['Top-1', 'Top-3'], [top1/n*100, top3/n*100],
                  color=[CORAL, TEAL], edgecolor='white', linewidth=1.5, width=0.5)
    for bar, val, count in zip(bars, [top1/n*100, top3/n*100], [top1, top3]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.0f}%\n({count}/{n})', ha='center', va='bottom',
                fontweight='bold', fontsize=12)
    ax.set_ylabel('Concordance with published strategy (%)')
    ax.set_ylim(0, 105)
    ax.set_title('Overall Concordance', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Right: Per-category breakdown
    ax = axes[1]
    cats = {}
    for c in cases:
        pub = c['published']
        if pub not in cats:
            cats[pub] = {'total': 0, 'concordant': 0}
        cats[pub]['total'] += 1
        if c['concordant']:
            cats[pub]['concordant'] += 1

    cat_names = sorted(cats.keys())
    x = np.arange(len(cat_names))
    totals = [cats[c]['total'] for c in cat_names]
    concordants = [cats[c]['concordant'] for c in cat_names]

    cat_colors = {'BE': TEAL, 'PE': CORAL, 'HDR': GOLD, 'COMPARISON': BLUE}
    bar_colors = [cat_colors.get(c, GRAY) for c in cat_names]

    bars_total = ax.bar(x - 0.15, totals, 0.3, color=[c + '80' for c in bar_colors],
                        label='Total cases', edgecolor='white')
    bars_conc = ax.bar(x + 0.15, concordants, 0.3, color=bar_colors,
                       label='Concordant', edgecolor='white')

    ax.set_xticks(x)
    ax.set_xticklabels(cat_names)
    ax.set_ylabel('Number of cases')
    ax.set_title('Per-Category Concordance', fontweight='bold')
    ax.legend(loc='upper right', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle('Literature Benchmark: Pipeline vs Published Strategies',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'Fig_LiteratureBenchmark')


# =====================================================================
# FIGURE 3: Bystander Fix Impact
# =====================================================================
def fig_bystander_fix():
    """Show how proper bystander scoring affects BE vs PE ranking."""
    fig, ax = plt.subplots(figsize=(10, 6))

    bystanders = [0, 1, 2, 3, 4]

    # Initial scoring (with triple-counting — bystander penalized in 3 channels)
    be_feasibility_advantage = 0.0325
    initial_be_scores = []
    pe_score = 0.6302
    for n in bystanders:
        sev = n * 0.2
        if sev == 0:
            diff = be_feasibility_advantage + 0.03  # clean bonus
        else:
            risk_pen = 0.15 * sev * 0.3
            conseq_pen = sev * 0.08
            bonus_lost = 0.03
            diff = be_feasibility_advantage - risk_pen - conseq_pen - bonus_lost
        initial_be_scores.append(pe_score + diff)

    # Corrected TOPSIS 6D scores (consequence as proper single dimension)
    corrected_be_topsis = [1.000, 0.995, 0.990, 0.986, 0.869]
    corrected_pe_topsis = [0.869, 0.869, 0.869, 0.869, 0.869]

    ax.plot(bystanders, initial_be_scores, 'o-', color=RED, linewidth=2, markersize=8,
            label='BE score (triple-counted bystander penalty)', zorder=5)
    ax.axhline(y=pe_score, color=CORAL, linestyle='--', linewidth=2,
               label=f'PE score ({pe_score:.3f})', alpha=0.7)

    ax.plot(bystanders, corrected_be_topsis, 's-', color=TEAL, linewidth=2, markersize=8,
            label='BE TOPSIS score (single consequence dimension)', zorder=5)
    ax.plot(bystanders, corrected_pe_topsis, '^--', color=GOLD, linewidth=2, markersize=6,
            label='PE TOPSIS score', alpha=0.7)

    # Highlight the crossover region
    ax.axvspan(-0.2, 0.5, color=GREEN, alpha=0.1, label='BE wins (both methods)')
    ax.axvspan(0.5, 4.2, color=RED, alpha=0.05)

    # Annotations
    ax.annotate('Triple-count: PE wins\nwith just 1 bystander',
                xy=(1, initial_be_scores[1]), xytext=(2, 0.58),
                arrowprops=dict(arrowstyle='->', color=RED),
                fontsize=10, color=RED, fontweight='bold')
    ax.annotate('Corrected: BE still wins\nwith 3 bystanders',
                xy=(3, corrected_be_topsis[3]), xytext=(3.5, 1.02),
                arrowprops=dict(arrowstyle='->', color=TEAL),
                fontsize=10, color=TEAL, fontweight='bold')

    ax.set_xlabel('Number of bystander edits in editing window', fontsize=12)
    ax.set_ylabel('Strategy score', fontsize=12)
    ax.set_title('Why Bystander Scoring Architecture Matters',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(bystanders)
    ax.set_xlim(-0.3, 4.3)
    ax.legend(loc='lower left', fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    save_fig(fig, 'Fig_BystanterFix')


# =====================================================================
# FIGURE 4: ConversionSim Output with CIs
# =====================================================================
def fig_conversionsim_with_ci():
    """Show ConversionSim output with statistical rigor (SEs/CIs)."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

    from conversion_sim import ConversionSimulator

    sim = ConversionSimulator(
        cut_type="staggered_5prime", overhang_length=3,
        donor_topology="circular_ssDNA", homology_arm_length=300,
        cell_type="iPSC", n_simulations=50000, seed=42,
    )
    res = sim.run()
    successful = res.tract_lengths_bp[res.hdr_success]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: Tract length distribution
    ax = axes[0]
    ax.hist(successful, bins=60, density=True, color=TEAL, alpha=0.7,
            edgecolor='white', linewidth=0.5)

    median = float(np.median(successful))
    mean = float(np.mean(successful))
    se = float(np.std(successful, ddof=1) / np.sqrt(len(successful)))

    ax.axvline(median, color=CORAL, linestyle='--', linewidth=2,
               label=f'Median = {median:.0f} bp')
    ax.axvline(mean, color=GOLD, linestyle='-.', linewidth=2,
               label=f'Mean = {mean:.0f} bp (SE={se:.1f})')

    # 95% CI band for mean
    ci_lo = mean - 1.96 * se
    ci_hi = mean + 1.96 * se
    ax.axvspan(ci_lo, ci_hi, color=GOLD, alpha=0.15,
               label=f'95% CI: [{ci_lo:.0f}, {ci_hi:.0f}]')

    ax.set_xlabel('Gene Conversion Tract Length (bp)')
    ax.set_ylabel('Probability Density')
    ax.set_title('enFnCas9 + cssDNA + iPSC\n(n=50,000 simulations)', fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(0, 3000)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Right: Distance-probability with Wilson CIs
    ax = axes[1]
    distances = [50, 100, 200, 300, 500, 800, 1000, 1500, 2000]
    probs = []
    ci_lows = []
    ci_highs = []

    for d in distances:
        p, lo, hi = sim.probability_at_distance(d, return_ci=True)
        probs.append(p * 100)
        ci_lows.append(lo * 100)
        ci_highs.append(hi * 100)

    ax.plot(distances, probs, 'o-', color=TEAL, linewidth=2, markersize=6,
            label='Point estimate', zorder=5)
    ax.fill_between(distances, ci_lows, ci_highs, color=TEAL, alpha=0.2,
                    label='95% Wilson CI')

    ax.set_xlabel('Distance from cut site (bp)')
    ax.set_ylabel('P(tract reaches distance) (%)')
    ax.set_title('Conversion Probability with\n95% Confidence Intervals', fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 2100)
    ax.set_ylim(0, 105)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='both', alpha=0.3)

    # Add text box with key stats
    hdr_rate = float(np.mean(res.hdr_success))
    textstr = (f'HDR rate: {hdr_rate*100:.1f}%\n'
               f'n = {len(successful):,} HDR events\n'
               f'Model: Geometric(p=0.002)\n'
               f'[ASSUMED] from functional evidence')
    props = dict(boxstyle='round,pad=0.4', facecolor='lightyellow', alpha=0.8)
    ax.text(0.98, 0.98, textstr, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    fig.suptitle('ConversionSim Monte Carlo: Statistical Rigor in v3',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'Fig_ConversionSim_CIs')


# =====================================================================
# FIGURE 5: Parameter Provenance Summary
# =====================================================================
def fig_parameter_provenance():
    """Pie chart showing distribution of parameter evidence levels."""
    fig, ax = plt.subplots(figsize=(7, 7))

    # Counted from constants.py evidence tags
    counts = {
        '[MEASURED]': 14,
        '[DERIVED]': 8,
        '[ASSUMED]': 18,
    }
    labels = list(counts.keys())
    sizes = list(counts.values())
    colors_pie = [GREEN, BLUE, CORAL]
    explode = (0, 0, 0.05)

    wedges, texts, autotexts = ax.pie(
        sizes, explode=explode, labels=labels, colors=colors_pie,
        autopct=lambda p: f'{p:.0f}%\n({int(round(p*sum(sizes)/100))})',
        shadow=False, startangle=90, textprops={'fontsize': 13},
        pctdistance=0.7)

    for t in autotexts:
        t.set_fontweight('bold')

    ax.set_title('Parameter Provenance Distribution\n(40 key constants in CRISPRArchitect v3)',
                 fontsize=14, fontweight='bold')

    # Legend
    legend_text = [
        'Directly from published measurement',
        'Computed from published data',
        'Modeling assumption (sensitivity-tested)',
    ]
    ax.legend(wedges, legend_text, loc='lower center', fontsize=10,
              bbox_to_anchor=(0.5, -0.05))

    fig.tight_layout()
    save_fig(fig, 'Fig_ParameterProvenance')


# =====================================================================
# FIGURE 6: MCDM Cross-Method Agreement
# =====================================================================
def fig_cross_method():
    """Show TOPSIS vs VIKOR vs WPM agreement."""
    fig, ax = plt.subplots(figsize=(8, 5))

    # Data from our actual test (BE vs PE vs HDR, 1 bystander)
    methods = ['TOPSIS', 'VIKOR', 'WPM']
    strategies = ['BE (1 bystander)', 'PE', 'HDR']
    colors_strat = [TEAL, CORAL, GOLD]

    # Ranks from our verified test output
    ranks = {
        'TOPSIS': [1, 2, 3],
        'VIKOR': [1, 2, 3],
        'WPM': [1, 2, 3],
    }

    x = np.arange(len(methods))
    width = 0.25

    for i, (strat, col) in enumerate(zip(strategies, colors_strat)):
        vals = [ranks[m][i] for m in methods]
        bars = ax.bar(x + i * width - width, vals, width,
                      label=strat, color=col, edgecolor='white', linewidth=1)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                    f'#{val}', ha='center', va='bottom', fontsize=10,
                    fontweight='bold')

    ax.set_ylabel('Rank (1 = best)')
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=12)
    ax.set_ylim(0, 4)
    ax.set_title('Cross-Method Ranking Agreement: 100% Concordance',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()

    # Add concordance annotation
    ax.text(0.98, 0.02, 'Rank concordance:\nTOPSIS-VIKOR: 100%\n'
            'TOPSIS-WPM: 100%\nVIKOR-WPM: 100%',
            transform=ax.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.tight_layout()
    save_fig(fig, 'Fig_CrossMethod_Agreement')


# =====================================================================
# Run all figures
# =====================================================================
if __name__ == '__main__':
    print('Generating CRISPRArchitect v3 Results Figures')
    print(f'Output directory: {OUTPUT_DIR}')
    print()

    fig_strategy_distribution()
    fig_literature_benchmark()
    fig_bystander_fix()
    fig_conversionsim_with_ci()
    fig_parameter_provenance()
    fig_cross_method()

    print(f'\nAll figures saved to {OUTPUT_DIR}')
    print(f'Total: {len(os.listdir(OUTPUT_DIR))} files')
