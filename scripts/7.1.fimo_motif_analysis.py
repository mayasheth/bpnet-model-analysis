import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import os
import pickle
import argparse

# Import analysis functions and style config from utils
from fimo_analysis_utils import (
    # Style configuration
    COLOR_P300_POS, COLOR_P300_NEG, COLORS_4_DISCRETE,
    CMAP_SEQUENTIAL, CMAP_DIVERGING, apply_style,
    # Analysis functions
    motif_enrichment_analysis, motif_score_correlation, motif_density_comparison,
    find_motif_pairs, pair_enrichment_analysis, spacing_analysis, orientation_enrichment,
    score_distribution_by_p300, score_percentile_enrichment, score_threshold_sweep
)

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
parser = argparse.ArgumentParser(description='FIMO motif analysis and plotting')
parser.add_argument('--plot-only', action='store_true',
                    help='Skip analysis, load pre-computed results from TSV files and regenerate plots')
args = parser.parse_args()

# Output directory for plots and results
OUTPUT_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/FIMO/elements_v1/analysis_v1"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

# File paths - FILL THESE IN
REGIONS_FILE = "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"
FIMO_RES_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/FIMO/elements_v1/"
MOTIF_NAMES = ["GATA_1", "STAT_0", "E2F_1", "ATF-CEBP_0", "FOS-JUN_0", "ELF_0", "CTCF_0", "ATF-CREB-JUN_0", "NF2L-NFE_0", "RFX_1",
               "NFI_0", "NFI_1", "CREB_0", "TAL_0", "RUNX_0", "ZNF_0", "ZNF_3"]

FIMO_FILES = [os.path.join(FIMO_RES_DIR, f"{motif}.tsv") for motif in MOTIF_NAMES]

# Print mode
print("="*80)
if args.plot_only:
    print("MODE: Plot-only (skipping analysis, loading pre-computed results)")
else:
    print("MODE: Full analysis")
print(f"OUTPUT DIR: {OUTPUT_DIR}")
print("="*80)

# Load regions file
print("\nLoading regions...")
regions_df = pd.read_csv(REGIONS_FILE, sep='\t')

# Create sequence_name to match FIMO format (chr:start-end)
regions_df['sequence_name'] = (regions_df['chrom'].astype(str) + ':' + 
                                regions_df['start'].astype(str) + '-' + 
                                regions_df['end'].astype(str))

print(f"Loaded {len(regions_df)} regions")
print(f"p300+ regions: {regions_df['EP300_peak_overlap'].sum()}")
print(f"p300- regions: {(~regions_df['EP300_peak_overlap'].astype(bool)).sum()}")

# Load all FIMO files and combine
print("\nLoading FIMO results...")
fimo_dfs = []
for fimo_file in glob(FIMO_FILES) if isinstance(FIMO_FILES, str) else FIMO_FILES:
    df = pd.read_csv(fimo_file, sep='\t')
    fimo_dfs.append(df)
    print(f"  Loaded {fimo_file}: {len(df)} hits for {df['motif_name'].unique()}")

motif_hits = pd.concat(fimo_dfs, ignore_index=True)
print(f"\nTotal motif hits: {len(motif_hits)}")
print(f"Unique motifs: {motif_hits['motif_name'].nunique()}")
print(f"Motifs: {motif_hits['motif_name'].unique()}")

# =============================================================================
# ** LOAD PRE-COMPUTED REGION MAPPING **
# =============================================================================

MAPPING_FILE = os.path.join(OUTPUT_DIR, "fimo_to_annotation_mapping.pkl")

print("\n** Loading pre-computed region mapping **")
print(f"Reading mapping from {MAPPING_FILE}...")

with open(MAPPING_FILE, 'rb') as f:
    mapping_data = pickle.load(f)
    mapping_dict = mapping_data['mapping_dict']
    mapping_df = mapping_data['mapping_df']

print(f"Loaded mapping for {len(mapping_dict)} regions")

# Apply mapping to motif_hits
print("\nApplying mapping to motif hits...")
motif_hits['mapped_region'] = motif_hits['sequence_name'].map(mapping_dict)

unmapped = motif_hits['mapped_region'].isna().sum()
print(f"Unmapped motif hits: {unmapped} ({100*unmapped/len(motif_hits):.2f}%)")

motif_hits = motif_hits[motif_hits['mapped_region'].notna()].copy()
motif_hits['sequence_name'] = motif_hits['mapped_region']
motif_hits = motif_hits.drop(columns=['mapped_region'])

print(f"Motif hits after mapping: {len(motif_hits)}")

# Filter regions to only those with FIMO data
valid_regions = set(mapping_dict.values())
regions_df = regions_df[regions_df['sequence_name'].isin(valid_regions)].copy()
print(f"Regions after filtering: {len(regions_df)}")

# Create p300_regions dictionary for functions
p300_regions = {}
for _, row in regions_df.iterrows():
    p300_regions[row['sequence_name']] = {
        'p300_bound': bool(row['EP300_peak_overlap']),
        'p300_signal': row['EP300.RPM'],
        'accessibility': row['DHS.RPM']
    }

# Create p300_signals DataFrame for correlation analyses
p300_signals = regions_df[['sequence_name', 'EP300.RPM', 'EP300_peak_overlap']].copy()
p300_signals.columns = ['sequence_name', 'p300_signal', 'p300_bound']

print(f"\nData loaded successfully!")

# =============================================================================
# 2. RUN ANALYSES OR LOAD PRE-COMPUTED RESULTS
# =============================================================================

if args.plot_only:
    # Load pre-computed results from TSV files
    print("\n" + "="*80)
    print("PLOT-ONLY MODE: Loading pre-computed results")
    print("="*80)

    print(f"\nLoading analysis results from {OUTPUT_DIR}...")
    enrichment_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_enrichment_fishers.tsv'), sep='\t')
    print(f"  Loaded motif_enrichment_fishers.tsv: {len(enrichment_results)} motifs")

    correlation_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_p300_correlation.tsv'), sep='\t')
    print(f"  Loaded motif_p300_correlation.tsv: {len(correlation_results)} entries")

    density_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_density_comparison.tsv'), sep='\t')
    print(f"  Loaded motif_density_comparison.tsv: {len(density_results)} motifs")

    motif_pairs = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_pairs_all.tsv'), sep='\t')
    print(f"  Loaded motif_pairs_all.tsv: {len(motif_pairs)} pairs")

    pair_enrichment_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'pair_enrichment.tsv'), sep='\t')
    print(f"  Loaded pair_enrichment.tsv: {len(pair_enrichment_results)} pair types")

    spacing_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'pair_spacing_analysis.tsv'), sep='\t')
    print(f"  Loaded pair_spacing_analysis.tsv: {len(spacing_results)} entries")

    orientation_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'pair_orientation_enrichment.tsv'), sep='\t')
    print(f"  Loaded pair_orientation_enrichment.tsv: {len(orientation_results)} entries")

    score_distribution_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_score_distribution_by_p300.tsv'), sep='\t')
    print(f"  Loaded motif_score_distribution_by_p300.tsv: {len(score_distribution_results)} motifs")

    score_percentile_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_score_percentile_enrichment.tsv'), sep='\t')
    print(f"  Loaded motif_score_percentile_enrichment.tsv: {len(score_percentile_results)} entries")

    score_threshold_results = pd.read_csv(os.path.join(OUTPUT_DIR, 'motif_score_threshold_sweep.tsv'), sep='\t')
    print(f"  Loaded motif_score_threshold_sweep.tsv: {len(score_threshold_results)} entries")

else:
    # Run full analyses
    print("\n" + "="*80)
    print("ANALYSIS 1: SINGLE MOTIF ENRICHMENT")
    print("="*80)

    # Analysis 1.1: Binary enrichment (Fisher's exact test)
    print("\n1.1 Running Fisher's exact test for motif enrichment...")
    enrichment_results = motif_enrichment_analysis(motif_hits, p300_regions)
    enrichment_results = enrichment_results.sort_values('p_value')
    print(enrichment_results.head(10))
    enrichment_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_enrichment_fishers.tsv'), sep='\t', index=False)

    # Analysis 1.2: Quantitative correlation
    print("\n1.2 Running correlation analysis with p300 signal...")
    correlation_results = motif_score_correlation(motif_hits, p300_signals)
    correlation_results = correlation_results.sort_values('p_value')
    print(correlation_results.head(20))
    correlation_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_p300_correlation.tsv'), sep='\t', index=False)

    # Analysis 1.3: Motif density comparison
    print("\n1.3 Running motif density analysis...")
    density_results = motif_density_comparison(motif_hits, p300_regions)
    density_results = density_results.sort_values('mann_whitney_p')
    print(density_results)
    density_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_density_comparison.tsv'), sep='\t', index=False)

    # =============================================================================
    # 3. RUN MOTIF PAIR ANALYSES
    # =============================================================================

    print("\n" + "="*80)
    print("ANALYSIS 2: MOTIF PAIR ENRICHMENT AND SPACING")
    print("="*80)

    # Analysis 2.1: Find motif pairs
    print("\n2.1 Finding motif pairs (max distance = 100)...")
    motif_pairs = find_motif_pairs(motif_hits, max_distance=100)
    print(f"Found {len(motif_pairs)} motif pairs")
    print(f"Unique pair types: {motif_pairs['pair_name'].nunique()}")
    print(f"Homotypic pairs: {(motif_pairs['pair_type'] == 'homotypic').sum()}")
    print(f"Heterotypic pairs: {(motif_pairs['pair_type'] == 'heterotypic').sum()}")
    motif_pairs.to_csv(os.path.join(OUTPUT_DIR, 'motif_pairs_all.tsv'), sep='\t', index=False)

    # Analysis 2.2: Pair enrichment
    print("\n2.2 Testing pair enrichment in p300 regions...")
    pair_enrichment_results = pair_enrichment_analysis(motif_pairs, p300_regions)
    pair_enrichment_results = pair_enrichment_results.sort_values('p_value')
    print(pair_enrichment_results.head(20))
    pair_enrichment_results.to_csv(os.path.join(OUTPUT_DIR, 'pair_enrichment.tsv'), sep='\t', index=False)

    # Analysis 2.3: Spacing analysis
    print("\n2.3 Analyzing spacing distributions...")
    spacing_results = spacing_analysis(motif_pairs, p300_regions)
    spacing_results = spacing_results.sort_values('mann_whitney_p')
    print(spacing_results.head(20))
    spacing_results.to_csv(os.path.join(OUTPUT_DIR, 'pair_spacing_analysis.tsv'), sep='\t', index=False)

    # Analysis 2.4: Orientation enrichment
    print("\n2.4 Testing orientation preferences...")
    orientation_results = orientation_enrichment(motif_pairs, p300_regions)
    orientation_results = orientation_results.sort_values('p_value')
    print(orientation_results.head(20))
    orientation_results.to_csv(os.path.join(OUTPUT_DIR, 'pair_orientation_enrichment.tsv'), sep='\t', index=False)

    # =============================================================================
    # ANALYSIS 3: MOTIF SCORE QUALITY VS P300 BINDING
    # =============================================================================

    print("\n" + "="*80)
    print("ANALYSIS 3: MOTIF SCORE QUALITY VS P300 BINDING")
    print("="*80)

    # Analysis 3.1: Score distribution comparison
    print("\n3.1 Comparing score distributions between p300+ and p300- regions...")
    score_distribution_results = score_distribution_by_p300(motif_hits, p300_regions)
    score_distribution_results = score_distribution_results.sort_values('p_value')
    print(score_distribution_results)
    score_distribution_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_score_distribution_by_p300.tsv'), sep='\t', index=False)

    # Analysis 3.2: Score percentile enrichment
    print("\n3.2 Testing score percentile enrichment (quartiles)...")
    score_percentile_results = score_percentile_enrichment(motif_hits, p300_regions, n_bins=4)
    score_percentile_results = score_percentile_results.sort_values(['motif', 'score_bin'])
    print(score_percentile_results.head(20))
    score_percentile_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_score_percentile_enrichment.tsv'), sep='\t', index=False)

    # Analysis 3.3: Score threshold sweep
    print("\n3.3 Running score threshold sweep...")
    score_threshold_results = score_threshold_sweep(motif_hits, p300_regions,
                                                     percentile_thresholds=[0, 25, 50, 75, 90])
    score_threshold_results = score_threshold_results.sort_values(['motif', 'percentile_threshold'])
    print(score_threshold_results.head(20))
    score_threshold_results.to_csv(os.path.join(OUTPUT_DIR, 'motif_score_threshold_sweep.tsv'), sep='\t', index=False)

# =============================================================================
# 3.5 CALCULATE OBSERVED/EXPECTED RATIO FOR PAIR ENRICHMENT
# =============================================================================
# This controls for individual motif frequencies when assessing pair enrichment
# Can be computed from pre-saved results (works in plot-only mode)

print("\n" + "="*80)
print("CALCULATING OBSERVED/EXPECTED PAIR ENRICHMENT")
print("="*80)

# Create lookup dict for individual motif frequencies
motif_freq_lookup = enrichment_results.set_index('motif')[['p300_pos_freq', 'p300_neg_freq']].to_dict('index')

# Calculate observed/expected for each pair
obs_exp_data = []
for _, row in pair_enrichment_results.iterrows():
    pair_name = row['pair_name']
    motifs = pair_name.split('-')

    if len(motifs) == 2:
        m1, m2 = motifs[0], motifs[1]

        # Get individual motif frequencies (handle missing motifs gracefully)
        if m1 in motif_freq_lookup and m2 in motif_freq_lookup:
            # Expected frequency under independence: P(m1) * P(m2)
            exp_p300_pos = motif_freq_lookup[m1]['p300_pos_freq'] * motif_freq_lookup[m2]['p300_pos_freq']
            exp_p300_neg = motif_freq_lookup[m1]['p300_neg_freq'] * motif_freq_lookup[m2]['p300_neg_freq']

            # Observed/Expected ratio
            obs_exp_p300_pos = row['p300_pos_freq'] / exp_p300_pos if exp_p300_pos > 0 else np.nan
            obs_exp_p300_neg = row['p300_neg_freq'] / exp_p300_neg if exp_p300_neg > 0 else np.nan

            # Adjusted enrichment: (Obs/Exp in p300+) / (Obs/Exp in p300-)
            # Values > 1 mean pair is more enriched in p300+ than expected from individual frequencies
            adj_enrichment = obs_exp_p300_pos / obs_exp_p300_neg if obs_exp_p300_neg > 0 else np.nan
        else:
            exp_p300_pos = np.nan
            exp_p300_neg = np.nan
            obs_exp_p300_pos = np.nan
            obs_exp_p300_neg = np.nan
            adj_enrichment = np.nan

        obs_exp_data.append({
            'pair_name': pair_name,
            'expected_p300_pos_freq': exp_p300_pos,
            'expected_p300_neg_freq': exp_p300_neg,
            'obs_exp_p300_pos': obs_exp_p300_pos,
            'obs_exp_p300_neg': obs_exp_p300_neg,
            'adjusted_enrichment': adj_enrichment
        })

obs_exp_df = pd.DataFrame(obs_exp_data)

# Merge back into pair_enrichment_results
pair_enrichment_results = pair_enrichment_results.merge(obs_exp_df, on='pair_name', how='left')

print(f"Calculated observed/expected ratios for {len(obs_exp_df)} pairs")
print("\nTop pairs by adjusted enrichment (controlling for individual motif frequencies):")
print(pair_enrichment_results.nlargest(10, 'adjusted_enrichment')[['pair_name', 'pair_type', 'odds_ratio', 'adjusted_enrichment', 'obs_exp_p300_pos']])

# =============================================================================
# 4. KEY VISUALIZATIONS
# =============================================================================

print("\n" + "="*80)
print("GENERATING PLOTS")
print("="*80)

# Plot 1: Motif enrichment volcano plot
print("\nPlot 1: Motif enrichment volcano plot...")
fig, ax = plt.subplots(figsize=(5, 5))
enrichment_plot = enrichment_results.copy()
enrichment_plot['-log10_pval'] = -np.log10(enrichment_plot['p_value'])
enrichment_plot['log2_OR'] = np.log2(enrichment_plot['odds_ratio'])

ax.scatter(enrichment_plot['log2_OR'], enrichment_plot['-log10_pval'],
           alpha=0.6, s=100, color=COLOR_P300_POS)

# Label significant motifs
sig_threshold = 0.05
for _, row in enrichment_plot[enrichment_plot['p_value'] < sig_threshold].iterrows():
    ax.text(row['log2_OR'], row['-log10_pval'], row['motif'],
            fontsize=9, ha='right')

ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
ax.axvline(0, color='black', linestyle='-', alpha=0.3)
ax.set_xlabel('log2(odds ratio)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Motif enrichment in p300+ regions\n(dashed line: p=0.05)')
apply_style(ax)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'plot1_motif_enrichment_volcano.png'), dpi=300)
plt.close()

# Plot 2: Motif frequency comparison (bar plot)
print("\nPlot 2: Motif frequency comparison...")
fig, ax = plt.subplots(figsize=(12, 6))

# Sort by decreasing odds ratio (enrichment)
enrichment_sorted = enrichment_results.sort_values('odds_ratio', ascending=False)

x = np.arange(len(enrichment_sorted))
width = 0.35

ax.bar(x - width/2, enrichment_sorted['p300_pos_freq'], width,
       label='p300+', alpha=0.8, color=COLOR_P300_POS)
ax.bar(x + width/2, enrichment_sorted['p300_neg_freq'], width,
       label='p300-', alpha=0.8, color=COLOR_P300_NEG)

ax.set_xlabel('Motif')
ax.set_ylabel('Frequency')
ax.set_title('Motif frequency: p300+ vs p300- regions')
ax.set_xticks(x)
ax.set_xticklabels(enrichment_sorted['motif'], rotation=45, ha='right')
ax.legend()
apply_style(ax)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'plot2_motif_frequency_comparison.png'), dpi=300)
plt.close()

# Plot 2.1: Motif density comparison (hits per kb)
print("\nPlot 2.1: Motif density comparison...")
fig, ax = plt.subplots(figsize=(12, 6))

# Sort by fold change (descending)
density_sorted = density_results.sort_values('fold_change', ascending=False)

x = np.arange(len(density_sorted))
width = 0.35

ax.bar(x - width/2, density_sorted['p300_pos_mean_density'], width,
       label='p300+', alpha=0.8, color=COLOR_P300_POS)
ax.bar(x + width/2, density_sorted['p300_neg_mean_density'], width,
       label='p300-', alpha=0.8, color=COLOR_P300_NEG)

# Add significance stars
for i, (_, row) in enumerate(density_sorted.iterrows()):
    max_height = max(row['p300_pos_mean_density'], row['p300_neg_mean_density'])
    if row['mann_whitney_p'] < 0.001:
        ax.text(i, max_height + 0.01, '***', ha='center', fontsize=10)
    elif row['mann_whitney_p'] < 0.01:
        ax.text(i, max_height + 0.01, '**', ha='center', fontsize=10)
    elif row['mann_whitney_p'] < 0.05:
        ax.text(i, max_height + 0.01, '*', ha='center', fontsize=10)

ax.set_xlabel('Motif')
ax.set_ylabel('Mean density (hits per kb)')
ax.set_title('Motif density: p300+ vs p300- regions\n(* p<0.05, ** p<0.01, *** p<0.001)')
ax.set_xticks(x)
ax.set_xticklabels(density_sorted['motif'], rotation=45, ha='right')
ax.legend(loc='upper right')
apply_style(ax)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'plot2.1_motif_density_comparison.png'), dpi=300)
plt.close()

# Plot 3: Pair enrichment heatmaps (odds ratio and observed/expected side-by-side)
print("\nPlot 3: Pair enrichment heatmaps (side-by-side)...")

# Parse pair names to get individual motifs and values
pair_data = []
for _, row in pair_enrichment_results.iterrows():
    motifs = row['pair_name'].split('-')
    if len(motifs) == 2:
        pair_data.append({
            'motif1': motifs[0],
            'motif2': motifs[1],
            'odds_ratio': row['odds_ratio'],
            'obs_exp_p300_pos': row.get('obs_exp_p300_pos', np.nan)
        })

pair_df = pd.DataFrame(pair_data)

# Create symmetric matrices for both metrics
symmetric_pairs = []
for _, row in pair_df.iterrows():
    symmetric_pairs.append(row.to_dict())
    if row['motif1'] != row['motif2']:
        symmetric_pairs.append({
            'motif1': row['motif2'],
            'motif2': row['motif1'],
            'odds_ratio': row['odds_ratio'],
            'obs_exp_p300_pos': row['obs_exp_p300_pos']
        })

symmetric_df = pd.DataFrame(symmetric_pairs)

# Pivot to create matrices
odds_matrix = symmetric_df.pivot_table(
    index='motif1', columns='motif2', values='odds_ratio', aggfunc='first'
).fillna(1)

obs_exp_matrix = symmetric_df.pivot_table(
    index='motif1', columns='motif2', values='obs_exp_p300_pos', aggfunc='first'
).fillna(1)

# Create side-by-side heatmaps
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Left: Odds ratio heatmap
sns.heatmap(odds_matrix, cmap=CMAP_DIVERGING, center=1,
            annot=True, fmt='.2f', cbar_kws={'label': 'Odds ratio'},
            square=True, linewidths=0.5, linecolor='gray', ax=ax1)
ax1.set_title('Odds ratio\n(raw co-occurrence enrichment)', fontsize=12)
ax1.tick_params(axis='both', colors='black')

# Right: Observed/Expected heatmap
sns.heatmap(obs_exp_matrix, cmap=CMAP_DIVERGING, center=1,
            annot=True, fmt='.2f', cbar_kws={'label': 'Observed/Expected'},
            square=True, linewidths=0.5, linecolor='gray', ax=ax2)
ax2.set_title('Observed/Expected\n(controlling for individual motif frequencies)', fontsize=12)
ax2.tick_params(axis='both', colors='black')

fig.suptitle('Motif pair enrichment in p300+ regions', fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'plot3_pair_enrichment_heatmaps.png'), dpi=300)
plt.close()

# Plot 4.1 and 4.2: Spacing distributions with enrichment (heterotypic and homotypic)
from matplotlib.gridspec import GridSpec

def plot_spacing_with_enrichment(pair_names, pair_type_label, output_filename):
    """Helper function to create spacing plots with stacked enrichment subplot."""
    fig = plt.figure(figsize=(18, 16))
    # 2 rows of pairs, 3 columns, each pair has 2 stacked subplots (density + enrichment, equal height)
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 1], hspace=0.35, wspace=0.3)

    for idx, pair_name in enumerate(pair_names):
        row = (idx // 3) * 2  # 0 or 2 (for density plots)
        col = idx % 3

        ax_density = fig.add_subplot(gs[row, col])
        ax_enrich = fig.add_subplot(gs[row + 1, col], sharex=ax_density)

        # Filter data (cut off at 95bp)
        pair_subset = motif_pairs[(motif_pairs['pair_name'] == pair_name) &
                                   (motif_pairs['distance'] < 95)]

        p300_pos_dist = []
        p300_neg_dist = []
        for _, r in pair_subset.iterrows():
            if p300_regions[r['sequence_name']]['p300_bound']:
                p300_pos_dist.append(r['distance'])
            else:
                p300_neg_dist.append(r['distance'])

        if len(p300_pos_dist) < 5 or len(p300_neg_dist) < 5:
            ax_density.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax_density.transAxes)
            ax_density.set_title(pair_name, fontsize=10, fontweight='bold')
            continue

        # Create histograms (cap at 95bp)
        max_dist = max(max(p300_pos_dist), max(p300_neg_dist))
        bins = np.arange(0, min(max_dist + 1, 96), 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        p300_pos_counts, _ = np.histogram(p300_pos_dist, bins=bins, density=True)
        p300_neg_counts, _ = np.histogram(p300_neg_dist, bins=bins, density=True)

        # Top plot: Density distributions
        ax_density.plot(bin_centers, p300_pos_counts, color=COLOR_P300_POS, linewidth=2,
                        label=f'p300+ (n={len(p300_pos_dist)})', alpha=0.8)
        ax_density.plot(bin_centers, p300_neg_counts, color=COLOR_P300_NEG, linewidth=2,
                        label=f'p300- (n={len(p300_neg_dist)})', alpha=0.8)
        ax_density.fill_between(bin_centers, p300_pos_counts, alpha=0.2, color=COLOR_P300_POS)
        ax_density.fill_between(bin_centers, p300_neg_counts, alpha=0.2, color=COLOR_P300_NEG)

        # Auto-scale y-axis (don't start at 0)
        all_counts = np.concatenate([p300_pos_counts, p300_neg_counts])
        ymin = max(0, all_counts[all_counts > 0].min() * 0.8) if any(all_counts > 0) else 0
        ymax = all_counts.max() * 1.1
        ax_density.set_ylim(ymin, ymax)

        ax_density.set_ylabel('Density', fontsize=9)
        ax_density.set_title(pair_name, fontsize=10, fontweight='bold')
        ax_density.legend(fontsize=7, loc='upper right')
        apply_style(ax_density)
        ax_density.set_xlim(0, min(max_dist, 95))
        plt.setp(ax_density.get_xticklabels(), visible=False)

        # Bottom plot: Enrichment ratio (p300+ / p300-)
        with np.errstate(divide='ignore', invalid='ignore'):
            enrichment = np.where(p300_neg_counts > 0, p300_pos_counts / p300_neg_counts, np.nan)

        ax_enrich.plot(bin_centers, enrichment, color='black', linewidth=2, alpha=0.8)
        ax_enrich.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
        ax_enrich.fill_between(bin_centers, 1, enrichment, where=enrichment > 1,
                               alpha=0.3, color=COLOR_P300_POS)
        ax_enrich.fill_between(bin_centers, 1, enrichment, where=enrichment < 1,
                               alpha=0.3, color=COLOR_P300_NEG)

        # Auto-scale y-axis for enrichment
        valid_enrich = enrichment[~np.isnan(enrichment)]
        if len(valid_enrich) > 0:
            enrich_min = max(0.5, valid_enrich.min() * 0.9)
            enrich_max = min(3, valid_enrich.max() * 1.1)
            ax_enrich.set_ylim(enrich_min, enrich_max)

        ax_enrich.set_xlabel('Distance (bp)', fontsize=9)
        ax_enrich.set_ylabel('Enrichment', fontsize=9)
        apply_style(ax_enrich)

    plt.suptitle(f'Top {pair_type_label} pairs: spacing distributions with enrichment',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.savefig(os.path.join(OUTPUT_DIR, output_filename), dpi=300, bbox_inches='tight')
    plt.close()

# Plot 4.1: Heterotypic pairs
print("\nPlot 4.1: Spacing distributions (heterotypic pairs)...")
heterotypic_pairs = pair_enrichment_results[pair_enrichment_results['pair_type'] == 'heterotypic']
top_heterotypic = heterotypic_pairs.nsmallest(6, 'p_value')['pair_name'].values
plot_spacing_with_enrichment(top_heterotypic, 'heterotypic', 'plot4.1_spacing_heterotypic.png')

# Plot 4.2: Homotypic pairs
print("\nPlot 4.2: Spacing distributions (homotypic pairs)...")
homotypic_pairs = pair_enrichment_results[pair_enrichment_results['pair_type'] == 'homotypic']
top_homotypic = homotypic_pairs.nsmallest(6, 'p_value')['pair_name'].values
plot_spacing_with_enrichment(top_homotypic, 'homotypic', 'plot4.2_spacing_homotypic.png')

# Plot 5.1 and 5.2: Orientation preferences (heterotypic and homotypic)
def plot_orientation_preferences(pair_names, pair_type_label, output_filename):
    """Helper function to create orientation preference plot for a set of pairs."""
    orient_data = []
    for pair_name in pair_names:
        pair_subset = motif_pairs[motif_pairs['pair_name'] == pair_name]

        for orientation in ['++', '+-', '-+', '--']:
            oriented_subset = pair_subset[pair_subset['orientation'] == orientation]

            p300_pos_count = sum(p300_regions[r]['p300_bound'] for r in oriented_subset['sequence_name'])
            p300_neg_count = len(oriented_subset) - p300_pos_count

            total_p300_pos = sum(p300_regions[r]['p300_bound'] for r in pair_subset['sequence_name'])
            total_p300_neg = len(pair_subset) - total_p300_pos

            p300_pos_frac = p300_pos_count / total_p300_pos if total_p300_pos > 0 else 0
            p300_neg_frac = p300_neg_count / total_p300_neg if total_p300_neg > 0 else 0

            orient_data.append({
                'pair_name': pair_name,
                'orientation': orientation,
                'p300+ fraction': p300_pos_frac,
                'p300- fraction': p300_neg_frac,
                'enrichment': p300_pos_frac / p300_neg_frac if p300_neg_frac > 0 else np.nan
            })

    orient_df = pd.DataFrame(orient_data)

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), height_ratios=[1, 0.6])
    (ax1, ax2), (ax3, ax4) = axes

    pivot_pos = orient_df.pivot(index='pair_name', columns='orientation', values='p300+ fraction')
    pivot_pos.plot(kind='bar', ax=ax1, color=COLORS_4_DISCRETE)
    ax1.set_xlabel('')
    ax1.set_ylabel('Fraction of pairs')
    ax1.set_title('Orientation preferences in p300+ regions')
    ax1.legend(title='Orientation')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    apply_style(ax1)

    pivot_neg = orient_df.pivot(index='pair_name', columns='orientation', values='p300- fraction')
    pivot_neg.plot(kind='bar', ax=ax2, color=COLORS_4_DISCRETE)
    ax2.set_xlabel('')
    ax2.set_ylabel('Fraction of pairs')
    ax2.set_title('Orientation preferences in p300- regions')
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
    apply_style(ax2)

    pivot_enrich = orient_df.pivot(index='pair_name', columns='orientation', values='enrichment')
    pivot_enrich.plot(kind='bar', ax=ax3, color=COLORS_4_DISCRETE)
    ax3.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Motif pair')
    ax3.set_ylabel('Enrichment (p300+ / p300-)')
    ax3.set_title('Orientation enrichment ratio')
    ax3.legend(title='Orientation')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right')
    apply_style(ax3)

    ax4.axis('off')

    plt.suptitle(f'Orientation preferences: top {pair_type_label} pairs', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, output_filename), dpi=300)
    plt.close()

# Plot 5.1: Heterotypic pairs
print("\nPlot 5.1: Orientation preferences (heterotypic pairs)...")
heterotypic_orient = pair_enrichment_results[pair_enrichment_results['pair_type'] == 'heterotypic']
top_heterotypic_orient = heterotypic_orient.nsmallest(5, 'p_value')['pair_name'].values
plot_orientation_preferences(top_heterotypic_orient, 'heterotypic', 'plot5.1_orientation_heterotypic.png')

# Plot 5.2: Homotypic pairs
print("\nPlot 5.2: Orientation preferences (homotypic pairs)...")
homotypic_orient = pair_enrichment_results[pair_enrichment_results['pair_type'] == 'homotypic']
top_homotypic_orient = homotypic_orient.nsmallest(5, 'p_value')['pair_name'].values
plot_orientation_preferences(top_homotypic_orient, 'homotypic', 'plot5.2_orientation_homotypic.png')


# Plot 6: p300 signal vs motif count (violin plots)
# Red line = mean, Black line = median
print("\nPlot 6: p300 signal vs motif count (violin plots)...")
top_motif = enrichment_results.iloc[0]['motif']
motif_counts = motif_hits[motif_hits['motif_name'] == top_motif].groupby('sequence_name').size().reset_index(name='count')
merged_plot = regions_df[['sequence_name', 'EP300.RPM', 'EP300_peak_overlap']].merge(
    motif_counts, on='sequence_name', how='left'
).fillna(0)

# Convert count to integer for cleaner grouping
merged_plot['count'] = merged_plot['count'].astype(int)

# Create violin plot
fig, ax = plt.subplots(figsize=(14, 6))

# Get unique motif counts and sort them
unique_counts = sorted(merged_plot['count'].unique())

# Prepare data for violin plot, filtering to n>20
plot_data = []
positions = []
colors = []

for count in unique_counts:
    subset = merged_plot[merged_plot['count'] == count]
    if len(subset) > 20:  # Only include groups with n>20
        # Add pseudocount of 1 to avoid log(0)
        plot_data.append(subset['EP300.RPM'].values + 1)
        positions.append(count)

        # Color by fraction with p300 peaks
        frac_with_peak = subset['EP300_peak_overlap'].mean()
        colors.append(frac_with_peak)

# Create violin plot
parts = ax.violinplot(plot_data, positions=positions, widths=0.7,
                       showmeans=True, showmedians=True)

# Color violins by p300 peak fraction using sequential colormap
cmap = cm.get_cmap(CMAP_SEQUENTIAL)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(cmap(colors[i]))
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')
    pc.set_linewidth(1)

# Style the mean and median lines
parts['cmedians'].set_color('black')
parts['cmedians'].set_linewidth(2)
parts['cmeans'].set_color('red')
parts['cmeans'].set_linewidth(1.5)

ax.set_xlabel(f'{top_motif} motif count', fontsize=12)
ax.set_ylabel('p300 signal (RPM + 1, log10)', fontsize=12)
ax.set_yscale('log')
ax.set_ylim(bottom=1)  # Start at log10(1) = 0
ax.set_title(f'p300 signal distribution by {top_motif} motif count\n(red line = mean, black line = median)', fontsize=14)

# Set x-axis to use only integer ticks
ax.set_xticks(positions)
ax.set_xticklabels([str(int(p)) for p in positions])

# Add colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Fraction in p300 peak', fontsize=10)

# Add sample size annotations
for i, count in enumerate(positions):
    n = len(plot_data[i])
    ax.text(count, ax.get_ylim()[1] * 0.8, f'n={n}',
            ha='center', va='top', fontsize=8, alpha=0.7)

apply_style(ax)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'plot6_p300_signal_vs_motif_count.png'), dpi=300)
plt.close()

# Plot 6.1: p300 signal vs motif count for top 6 motifs
print("\nPlot 6.1: p300 signal vs motif count (top 6 motifs)...")
top_6_motifs = enrichment_results.nsmallest(6, 'p_value')['motif'].values

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, motif_name in enumerate(top_6_motifs):
    ax = axes[idx]

    # Get motif counts for this motif
    motif_counts = motif_hits[motif_hits['motif_name'] == motif_name].groupby('sequence_name').size().reset_index(name='count')
    merged_plot = regions_df[['sequence_name', 'EP300.RPM', 'EP300_peak_overlap']].merge(
        motif_counts, on='sequence_name', how='left'
    ).fillna(0)
    merged_plot['count'] = merged_plot['count'].astype(int)

    # Get unique counts and prepare data
    unique_counts = sorted(merged_plot['count'].unique())
    plot_data = []
    positions = []
    colors = []

    for count in unique_counts:
        subset = merged_plot[merged_plot['count'] == count]
        if len(subset) > 20:
            # Add pseudocount of 1 to avoid log(0)
            plot_data.append(subset['EP300.RPM'].values + 1)
            positions.append(count)
            colors.append(subset['EP300_peak_overlap'].mean())

    if len(plot_data) == 0:
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(motif_name)
        continue

    # Create violin plot
    parts = ax.violinplot(plot_data, positions=positions, widths=0.7,
                           showmeans=True, showmedians=True)

    # Color violins
    cmap = cm.get_cmap(CMAP_SEQUENTIAL)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cmap(colors[i]))
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)

    parts['cmedians'].set_color('black')
    parts['cmedians'].set_linewidth(2)
    parts['cmeans'].set_color('red')
    parts['cmeans'].set_linewidth(1.5)

    ax.set_xlabel('Motif count', fontsize=10)
    ax.set_ylabel('p300 signal (RPM + 1)', fontsize=10)
    ax.set_yscale('log')
    ax.set_ylim(bottom=1)  # Start at log10(1) = 0
    ax.set_title(motif_name, fontsize=11, fontweight='bold')
    ax.set_xticks(positions)
    ax.set_xticklabels([str(int(p)) for p in positions])
    apply_style(ax)

plt.suptitle('p300 signal by motif count (top 6 enriched motifs)\n(red=mean, black=median)',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'plot6.1_p300_signal_by_count_top6.png'), dpi=300)
plt.close()

# =============================================================================
# PLOTS 7-9: MOTIF SCORE QUALITY ANALYSES
# =============================================================================

# Plot 7: Score distribution comparison (violin plots for top motifs)
print("\nPlot 7: Score distribution by p300 status (violin plots)...")

# Get top 6 motifs by enrichment significance
top_motifs_for_score = score_distribution_results.nsmallest(6, 'p_value')['motif'].values

# Pre-compute p300 status for motif_hits
p300_status_for_hits = motif_hits['sequence_name'].map(
    lambda x: p300_regions[x]['p300_bound'] if x in p300_regions else None
)
motif_hits_plot = motif_hits.copy()
motif_hits_plot['p300_bound'] = p300_status_for_hits

fig, axes = plt.subplots(2, 3, figsize=(9, 6))
axes = axes.flatten()

for idx, motif_name in enumerate(top_motifs_for_score):
    ax = axes[idx]
    motif_data = motif_hits_plot[motif_hits_plot['motif_name'] == motif_name]

    # Prepare data for violin plot
    p300_pos_scores = motif_data.loc[motif_data['p300_bound'] == True, 'score'].values
    p300_neg_scores = motif_data.loc[motif_data['p300_bound'] == False, 'score'].values

    # Create violin plot
    parts = ax.violinplot([p300_pos_scores, p300_neg_scores], positions=[1, 2],
                           showmeans=True, showmedians=True)

    # Color the violins
    parts['bodies'][0].set_facecolor(COLOR_P300_POS)
    parts['bodies'][0].set_alpha(0.7)
    parts['bodies'][1].set_facecolor(COLOR_P300_NEG)
    parts['bodies'][1].set_alpha(0.7)

    # Style mean/median lines
    parts['cmedians'].set_color('black')
    parts['cmedians'].set_linewidth(2)
    parts['cmeans'].set_color('red')
    parts['cmeans'].set_linewidth(1.5)

    # Get stats for title
    stats = score_distribution_results[score_distribution_results['motif'] == motif_name].iloc[0]
    fc = stats['score_fold_change']
    pval = stats['p_value']

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['p300+', 'p300-'], fontsize=8)
    ax.set_ylabel('FIMO score', fontsize=8)
    ax.set_title(f'{motif_name}\nFC={fc:.2f}, p={pval:.2e}', fontsize=8)
    apply_style(ax)

fig.suptitle('Score distributions: p300+ vs p300- (red=mean, black=median)', fontsize=10, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(OUTPUT_DIR, 'plot7_score_distribution_violin.png'), dpi=300)
plt.close()

# Plot 8: Score percentile enrichment (grouped bar plot)
print("\nPlot 8: Score percentile enrichment by quartile...")

# Get top 6 motifs
top_motifs_pct = score_distribution_results.nsmallest(6, 'p_value')['motif'].values

fig, axes = plt.subplots(2, 3, figsize=(9, 6))
axes = axes.flatten()

for idx, motif_name in enumerate(top_motifs_pct):
    ax = axes[idx]
    motif_pct_data = score_percentile_results[score_percentile_results['motif'] == motif_name]

    if len(motif_pct_data) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(motif_name)
        continue

    # Plot odds ratios by quartile
    x = np.arange(len(motif_pct_data))
    bars = ax.bar(x, motif_pct_data['odds_ratio'], color=COLOR_P300_POS, alpha=0.8)

    # Add significance stars
    for i, (_, row) in enumerate(motif_pct_data.iterrows()):
        if row['p_value'] < 0.001:
            ax.text(i, row['odds_ratio'] + 0.02, '***', ha='center', fontsize=8)
        elif row['p_value'] < 0.01:
            ax.text(i, row['odds_ratio'] + 0.02, '**', ha='center', fontsize=8)
        elif row['p_value'] < 0.05:
            ax.text(i, row['odds_ratio'] + 0.02, '*', ha='center', fontsize=8)

    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(motif_pct_data['score_bin'], fontsize=7)
    ax.set_xlabel('Quartile (Q1=low, Q4=high)', fontsize=7)
    ax.set_ylabel('Odds ratio', fontsize=7)
    ax.set_title(motif_name, fontsize=9, fontweight='bold')
    apply_style(ax)

fig.suptitle('Enrichment by score quartile (* p<0.05, ** p<0.01, *** p<0.001)', fontsize=10, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(OUTPUT_DIR, 'plot8_score_percentile_enrichment.png'), dpi=300)
plt.close()

# Plot 9: Score threshold sweep (line plot)
print("\nPlot 9: Score threshold sweep...")

fig, axes = plt.subplots(2, 3, figsize=(9, 6))
axes = axes.flatten()

for idx, motif_name in enumerate(top_motifs_pct):
    ax = axes[idx]
    motif_sweep = score_threshold_results[score_threshold_results['motif'] == motif_name]

    if len(motif_sweep) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(motif_name)
        continue

    # Sort by threshold
    motif_sweep = motif_sweep.sort_values('percentile_threshold')

    # X-axis: percentage of hits kept (inverted from percentile threshold)
    x = 100 - motif_sweep['percentile_threshold'].values
    y = motif_sweep['log2_odds_ratio'].values

    ax.plot(x, y, marker='o', linewidth=2, markersize=6, color=COLOR_P300_POS)
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)

    # Add significance markers
    for i, (_, row) in enumerate(motif_sweep.iterrows()):
        if row['p_value'] < 0.05:
            ax.scatter([100 - row['percentile_threshold']], [row['log2_odds_ratio']],
                      s=100, facecolors='none', edgecolors='red', linewidths=2)

    ax.set_xlabel('% hits kept', fontsize=7)
    ax.set_ylabel('log2(OR)', fontsize=7)
    ax.set_title(motif_name, fontsize=9, fontweight='bold')
    ax.set_xlim(0, 105)
    ax.invert_xaxis()  # Higher stringency on right
    apply_style(ax)

fig.suptitle('Enrichment vs stringency (red circles: p<0.05; right = stricter)', fontsize=10, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(OUTPUT_DIR, 'plot9_score_threshold_sweep.png'), dpi=300)
plt.close()

print("\n" + "="*80)
if args.plot_only:
    print("PLOTTING COMPLETE! (plot-only mode)")
else:
    print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nOutput directory: {OUTPUT_DIR}")

if not args.plot_only:
    print("\nTSV files:")
    print("  - motif_enrichment_fishers.tsv")
    print("  - motif_p300_correlation.tsv")
    print("  - motif_density_comparison.tsv")
    print("  - motif_pairs_all.tsv")
    print("  - pair_enrichment.tsv")
    print("  - pair_spacing_analysis.tsv")
    print("  - pair_orientation_enrichment.tsv")
    print("  - motif_score_distribution_by_p300.tsv")
    print("  - motif_score_percentile_enrichment.tsv")
    print("  - motif_score_threshold_sweep.tsv")

print("\nPlots:")
print("  - plot1_motif_enrichment_volcano.png")
print("  - plot2_motif_frequency_comparison.png (sorted by decreasing enrichment)")
print("  - plot2.1_motif_density_comparison.png (hits per kb, sorted by fold change)")
print("  - plot3_pair_enrichment_heatmaps.png (odds ratio + obs/exp side-by-side)")
print("  - plot4.1_spacing_heterotypic.png (with enrichment subplot)")
print("  - plot4.2_spacing_homotypic.png (with enrichment subplot)")
print("  - plot5.1_orientation_heterotypic.png")
print("  - plot5.2_orientation_homotypic.png")
print("  - plot6_p300_signal_vs_motif_count.png (log10 scale)")
print("  - plot6.1_p300_signal_by_count_top6.png (top 6 motifs)")
print("  - plot7_score_distribution_violin.png (score distributions p300+ vs p300-)")
print("  - plot8_score_percentile_enrichment.png (enrichment by score quartile)")
print("  - plot9_score_threshold_sweep.png (enrichment vs score stringency)")