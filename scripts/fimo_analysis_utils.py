"""
Utility functions for FIMO motif analysis.

Contains:
- Style configuration for plots
- Analysis functions for motif enrichment, correlations, pairs, etc.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu, ks_2samp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# =============================================================================
# STYLE CONFIGURATION (see STYLE_GUIDELINES.md)
# =============================================================================

# Color palettes
COLOR_P300_POS = '#792374'  # purple for p300+
COLOR_P300_NEG = '#49bcbc'  # cyan for p300-
COLORS_4_DISCRETE = ['#e96a00', '#c5373d', '#0096a0', '#429130']  # orange, red, teal, green
CMAP_SEQUENTIAL = 'PuBu'

# Diverging colormap 'managua' - create custom if not available
try:
    plt.colormaps['managua']
    CMAP_DIVERGING = 'managua'
except (KeyError, ValueError):
    # Fallback: create managua-like diverging colormap (teal-white-orange)
    managua_colors = ['#005f73', '#0a9396', '#94d2bd', '#e9d8a6', '#ee9b00', '#ca6702', '#9b2226']
    CMAP_DIVERGING = LinearSegmentedColormap.from_list('managua', managua_colors)


def apply_style(ax):
    """Apply consistent styling to axes."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.tick_params(axis='both', colors='black')
    ax.grid(False)


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def motif_enrichment_analysis(motif_hits, p300_regions):
    """
    Test if motif presence is enriched in P300+ vs P300- regions.
    Uses Fisher's exact test.
    """
    # Pre-compute p300 status sets once
    all_regions = set(p300_regions.keys())
    p300_pos_regions = {r for r, data in p300_regions.items() if data['p300_bound']}
    p300_neg_regions = all_regions - p300_pos_regions
    n_p300_pos_total = len(p300_pos_regions)
    n_p300_neg_total = len(p300_neg_regions)

    results = []

    for motif_name in motif_hits['motif_name'].unique():
        # Get regions with this motif
        regions_with_motif = set(motif_hits[motif_hits['motif_name'] == motif_name]['sequence_name'])

        # Build contingency table using set operations
        p300_pos_with_motif = len(regions_with_motif & p300_pos_regions)
        p300_neg_with_motif = len(regions_with_motif & p300_neg_regions)
        p300_pos_no_motif = n_p300_pos_total - p300_pos_with_motif
        p300_neg_no_motif = n_p300_neg_total - p300_neg_with_motif

        # Fisher's exact test
        odds_ratio, p_value = fisher_exact([[p300_pos_with_motif, p300_pos_no_motif],
                                            [p300_neg_with_motif, p300_neg_no_motif]])

        results.append({
            'motif': motif_name,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'p300_pos_freq': p300_pos_with_motif / n_p300_pos_total if n_p300_pos_total > 0 else 0,
            'p300_neg_freq': p300_neg_with_motif / n_p300_neg_total if n_p300_neg_total > 0 else 0
        })

    return pd.DataFrame(results)


def motif_score_correlation(motif_hits, p300_signals):
    """Correlate motif count/score with continuous P300 signal."""
    results = []

    for motif_name in motif_hits['motif_name'].unique():
        # Aggregate motif information per region
        motif_data = motif_hits[motif_hits['motif_name'] == motif_name].groupby('sequence_name').agg({
            'score': ['count', 'sum', 'max', 'mean']
        }).reset_index()

        # Flatten the multi-level columns
        motif_data.columns = ['sequence_name', 'count', 'sum', 'max', 'mean']

        # Merge with P300 signal
        merged = pd.merge(p300_signals, motif_data, on='sequence_name', how='left').fillna(0)

        # Correlation tests
        for metric in ['count', 'sum', 'max', 'mean']:
            corr, p_val = spearmanr(merged['p300_signal'], merged[metric])
            results.append({
                'motif': motif_name,
                'metric': metric,
                'spearman_r': corr,
                'p_value': p_val
            })

    return pd.DataFrame(results)


def motif_density_comparison(motif_hits, p300_regions):
    """
    Compare motif density (hits per kb) between P300+ and P300- regions.
    """
    # Pre-compute region lengths and P300 status
    region_info = []
    for region_name, region_data in p300_regions.items():
        coords = region_name.split(':')[1].split('-')
        length_kb = (int(coords[1]) - int(coords[0])) / 1000
        region_info.append({
            'sequence_name': region_name,
            'length_kb': length_kb,
            'p300_bound': region_data['p300_bound']
        })
    region_df = pd.DataFrame(region_info)

    results = []

    for motif_name in motif_hits['motif_name'].unique():
        # Count motifs per region using groupby
        motif_counts = motif_hits[motif_hits['motif_name'] == motif_name].groupby(
            'sequence_name'
        ).size().reset_index(name='n_motifs')

        # Merge with region info
        merged = region_df.merge(motif_counts, on='sequence_name', how='left').fillna(0)
        merged['density'] = merged['n_motifs'] / merged['length_kb']

        # Split by P300 status
        p300_pos_density = merged[merged['p300_bound']]['density'].values
        p300_neg_density = merged[~merged['p300_bound']]['density'].values

        u_stat, p_val = mannwhitneyu(p300_pos_density, p300_neg_density, alternative='greater')

        results.append({
            'motif': motif_name,
            'p300_pos_mean_density': np.mean(p300_pos_density),
            'p300_neg_mean_density': np.mean(p300_neg_density),
            'fold_change': np.mean(p300_pos_density) / (np.mean(p300_neg_density) + 1e-10),
            'mann_whitney_p': p_val
        })

    return pd.DataFrame(results)


def find_motif_pairs(motif_hits, max_distance=500):
    """
    Identify all motif pairs within specified distance in each region.
    """
    # Sort once globally
    motif_hits_sorted = motif_hits.sort_values(['sequence_name', 'start']).reset_index(drop=True)

    # Group by region and get indices
    grouped = motif_hits_sorted.groupby('sequence_name')

    pairs = []

    for region_name, group in grouped:
        # Convert to numpy for faster access
        starts = group['start'].values
        ends = group['end'].values
        motif_names = group['motif_name'].values
        strands = group['strand'].values
        scores = group['score'].values
        n_motifs = len(group)

        # Pairwise comparison with vectorized distance check
        for i in range(n_motifs):
            # Only check motifs that could be within max_distance
            for j in range(i + 1, n_motifs):
                distance = starts[j] - ends[i]

                if distance > max_distance:
                    break  # All subsequent motifs are too far

                # Determine pair type
                if motif_names[i] == motif_names[j]:
                    pair_type = 'homotypic'
                    pair_name = f"{motif_names[i]}-{motif_names[i]}"
                else:
                    pair_type = 'heterotypic'
                    pair_name = '-'.join(sorted([motif_names[i], motif_names[j]]))

                pairs.append({
                    'sequence_name': region_name,
                    'pair_name': pair_name,
                    'pair_type': pair_type,
                    'motif1': motif_names[i],
                    'motif2': motif_names[j],
                    'distance': distance,
                    'orientation': f"{strands[i]}{strands[j]}",
                    'motif1_score': scores[i],
                    'motif2_score': scores[j]
                })

    return pd.DataFrame(pairs)


def pair_enrichment_analysis(motif_pairs, p300_regions):
    """
    Test enrichment of motif pairs in P300+ regions.
    """
    # Pre-compute P300 status for all regions
    all_regions = set(p300_regions.keys())
    p300_pos_regions = {r for r, data in p300_regions.items() if data['p300_bound']}
    p300_neg_regions = all_regions - p300_pos_regions

    n_p300_pos_total = len(p300_pos_regions)
    n_p300_neg_total = len(p300_neg_regions)

    results = []

    for pair_name in motif_pairs['pair_name'].unique():
        # Get regions with this pair
        regions_with_pair = set(motif_pairs[motif_pairs['pair_name'] == pair_name]['sequence_name'])

        # Contingency table using set operations
        p300_pos_with_pair = len(regions_with_pair & p300_pos_regions)
        p300_pos_no_pair = n_p300_pos_total - p300_pos_with_pair
        p300_neg_with_pair = len(regions_with_pair & p300_neg_regions)
        p300_neg_no_pair = n_p300_neg_total - p300_neg_with_pair

        odds_ratio, p_value = fisher_exact([[p300_pos_with_pair, p300_pos_no_pair],
                                            [p300_neg_with_pair, p300_neg_no_pair]])

        # Get pair type (just once, not from every row)
        pair_type = motif_pairs.loc[motif_pairs['pair_name'] == pair_name, 'pair_type'].iloc[0]

        results.append({
            'pair_name': pair_name,
            'pair_type': pair_type,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'n_p300_pos': p300_pos_with_pair,
            'n_p300_neg': p300_neg_with_pair,
            'p300_pos_freq': p300_pos_with_pair / n_p300_pos_total if n_p300_pos_total > 0 else 0,
            'p300_neg_freq': p300_neg_with_pair / n_p300_neg_total if n_p300_neg_total > 0 else 0
        })

    return pd.DataFrame(results)


def spacing_analysis(motif_pairs, p300_regions):
    """Analyze spacing distributions for enriched pairs."""
    # Pre-compute p300 status as a vectorized lookup
    p300_status_map = {r: data['p300_bound'] for r, data in p300_regions.items()}
    motif_pairs = motif_pairs.copy()
    motif_pairs['p300_bound'] = motif_pairs['sequence_name'].map(p300_status_map)

    results = []

    for pair_name in motif_pairs['pair_name'].unique():
        pair_subset = motif_pairs[motif_pairs['pair_name'] == pair_name]

        # Vectorized separation by p300 status
        p300_pos_distances = pair_subset.loc[pair_subset['p300_bound'] == True, 'distance'].values
        p300_neg_distances = pair_subset.loc[pair_subset['p300_bound'] == False, 'distance'].values

        if len(p300_pos_distances) > 5 and len(p300_neg_distances) > 5:
            # Statistical tests
            mw_stat, mw_p = mannwhitneyu(p300_pos_distances, p300_neg_distances)
            ks_stat, ks_p = ks_2samp(p300_pos_distances, p300_neg_distances)

            results.append({
                'pair_name': pair_name,
                'p300_pos_median_spacing': np.median(p300_pos_distances),
                'p300_neg_median_spacing': np.median(p300_neg_distances),
                'p300_pos_mean_spacing': np.mean(p300_pos_distances),
                'p300_neg_mean_spacing': np.mean(p300_neg_distances),
                'mann_whitney_p': mw_p,
                'ks_test_p': ks_p,
                'n_p300_pos_pairs': len(p300_pos_distances),
                'n_p300_neg_pairs': len(p300_neg_distances)
            })

    return pd.DataFrame(results)


def orientation_enrichment(motif_pairs, p300_regions):
    """Test if specific orientations (++, +-, -+, --) are enriched."""
    # Pre-compute p300 status sets
    p300_pos_regions = {r for r, data in p300_regions.items() if data['p300_bound']}

    results = []

    for pair_name in motif_pairs['pair_name'].unique():
        pair_subset = motif_pairs[motif_pairs['pair_name'] == pair_name]

        # Pre-compute totals for this pair using set operations
        all_pair_regions = set(pair_subset['sequence_name'])
        all_p300_pos = len(all_pair_regions & p300_pos_regions)
        all_p300_neg = len(all_pair_regions) - all_p300_pos

        for orientation in ['++', '+-', '-+', '--']:
            oriented_subset = pair_subset[pair_subset['orientation'] == orientation]
            regions_with_oriented_pair = set(oriented_subset['sequence_name'])

            p300_pos = len(regions_with_oriented_pair & p300_pos_regions)
            p300_neg = len(regions_with_oriented_pair) - p300_pos

            if p300_pos + p300_neg > 0:
                odds_ratio, p_value = fisher_exact([[p300_pos, all_p300_pos - p300_pos],
                                                    [p300_neg, all_p300_neg - p300_neg]])

                results.append({
                    'pair_name': pair_name,
                    'orientation': orientation,
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'n_p300_pos': p300_pos,
                    'fraction_this_orientation': (p300_pos + p300_neg) / (all_p300_pos + all_p300_neg) if (all_p300_pos + all_p300_neg) > 0 else 0
                })

    return pd.DataFrame(results)


def score_distribution_by_p300(motif_hits, p300_regions):
    """
    Compare FIMO score distributions between p300+ and p300- regions for each motif.
    Tests whether p300+ regions have systematically higher-scoring motif matches.
    """
    # Pre-compute p300 status as array for fast lookup
    p300_status = motif_hits['sequence_name'].map(
        lambda x: p300_regions[x]['p300_bound'] if x in p300_regions else None
    )
    motif_hits_with_status = motif_hits.copy()
    motif_hits_with_status['p300_bound'] = p300_status
    motif_hits_with_status = motif_hits_with_status[motif_hits_with_status['p300_bound'].notna()]

    results = []

    for motif_name in motif_hits_with_status['motif_name'].unique():
        motif_subset = motif_hits_with_status[motif_hits_with_status['motif_name'] == motif_name]

        p300_pos_scores = motif_subset.loc[motif_subset['p300_bound'], 'score'].values
        p300_neg_scores = motif_subset.loc[~motif_subset['p300_bound'], 'score'].values

        if len(p300_pos_scores) >= 5 and len(p300_neg_scores) >= 5:
            # Mann-Whitney U test (are p300+ scores higher?)
            u_stat, p_val = mannwhitneyu(p300_pos_scores, p300_neg_scores, alternative='greater')

            results.append({
                'motif': motif_name,
                'n_p300_pos_hits': len(p300_pos_scores),
                'n_p300_neg_hits': len(p300_neg_scores),
                'p300_pos_median_score': np.median(p300_pos_scores),
                'p300_neg_median_score': np.median(p300_neg_scores),
                'p300_pos_mean_score': np.mean(p300_pos_scores),
                'p300_neg_mean_score': np.mean(p300_neg_scores),
                'score_fold_change': np.median(p300_pos_scores) / np.median(p300_neg_scores),
                'mann_whitney_u': u_stat,
                'p_value': p_val
            })

    return pd.DataFrame(results)


def score_percentile_enrichment(motif_hits, p300_regions, n_bins=4):
    """
    Bin motif hits by within-motif score percentile and calculate p300+ enrichment per bin.
    Tests whether top-scoring hits are disproportionately in p300+ regions.
    """
    # Pre-compute p300 status
    p300_status_map = {r: data['p300_bound'] for r, data in p300_regions.items()}

    # Count total p300+ and p300- regions (for enrichment denominator)
    n_p300_pos_total = sum(1 for v in p300_status_map.values() if v)
    n_p300_neg_total = sum(1 for v in p300_status_map.values() if not v)

    results = []

    for motif_name in motif_hits['motif_name'].unique():
        motif_subset = motif_hits[motif_hits['motif_name'] == motif_name].copy()

        # Compute within-motif score percentiles
        motif_subset['score_percentile'] = motif_subset['score'].rank(pct=True)

        # Bin into quartiles (or n_bins)
        bin_edges = np.linspace(0, 1, n_bins + 1)
        bin_labels = [f'Q{i+1}' for i in range(n_bins)]
        motif_subset['score_bin'] = pd.cut(
            motif_subset['score_percentile'],
            bins=bin_edges,
            labels=bin_labels,
            include_lowest=True
        )

        # Map p300 status
        motif_subset['p300_bound'] = motif_subset['sequence_name'].map(p300_status_map)
        motif_subset = motif_subset[motif_subset['p300_bound'].notna()]

        for bin_label in bin_labels:
            bin_subset = motif_subset[motif_subset['score_bin'] == bin_label]

            if len(bin_subset) < 10:
                continue

            # Get unique regions with hits in this bin
            regions_with_hit = set(bin_subset['sequence_name'])

            # Contingency table
            p300_pos_with_hit = len([r for r in regions_with_hit if p300_status_map.get(r, False)])
            p300_neg_with_hit = len([r for r in regions_with_hit if not p300_status_map.get(r, True)])
            p300_pos_no_hit = n_p300_pos_total - p300_pos_with_hit
            p300_neg_no_hit = n_p300_neg_total - p300_neg_with_hit

            # Fisher's exact test
            odds_ratio, p_value = fisher_exact([
                [p300_pos_with_hit, p300_pos_no_hit],
                [p300_neg_with_hit, p300_neg_no_hit]
            ])

            # Score range for this bin
            bin_scores = bin_subset['score']

            results.append({
                'motif': motif_name,
                'score_bin': bin_label,
                'bin_percentile_range': f'{bin_edges[bin_labels.index(bin_label)]:.0%}-{bin_edges[bin_labels.index(bin_label)+1]:.0%}',
                'n_hits': len(bin_subset),
                'n_regions': len(regions_with_hit),
                'min_score': bin_scores.min(),
                'max_score': bin_scores.max(),
                'median_score': bin_scores.median(),
                'p300_pos_with_hit': p300_pos_with_hit,
                'p300_neg_with_hit': p300_neg_with_hit,
                'odds_ratio': odds_ratio,
                'p_value': p_value
            })

    return pd.DataFrame(results)


def score_threshold_sweep(motif_hits, p300_regions, percentile_thresholds=None):
    """
    Vary score threshold and measure how enrichment changes.
    Tests whether keeping only high-scoring hits increases enrichment.
    """
    if percentile_thresholds is None:
        percentile_thresholds = [0, 25, 50, 75, 90]  # Keep top 100%, 75%, 50%, 25%, 10%

    # Pre-compute p300 status
    p300_status_map = {r: data['p300_bound'] for r, data in p300_regions.items()}
    n_p300_pos_total = sum(1 for v in p300_status_map.values() if v)
    n_p300_neg_total = sum(1 for v in p300_status_map.values() if not v)

    results = []

    for motif_name in motif_hits['motif_name'].unique():
        motif_subset = motif_hits[motif_hits['motif_name'] == motif_name].copy()

        # Pre-compute score percentiles for this motif
        score_percentiles = motif_subset['score'].rank(pct=True).values
        scores = motif_subset['score'].values
        regions = motif_subset['sequence_name'].values

        for pct_threshold in percentile_thresholds:
            # Keep hits above this percentile threshold
            mask = score_percentiles >= (pct_threshold / 100)
            filtered_regions = regions[mask]
            filtered_scores = scores[mask]

            if len(filtered_regions) < 10:
                continue

            # Unique regions passing threshold
            regions_passing = set(filtered_regions)

            # Contingency table
            p300_pos_with_hit = len([r for r in regions_passing if p300_status_map.get(r, False)])
            p300_neg_with_hit = len([r for r in regions_passing if not p300_status_map.get(r, True)])
            p300_pos_no_hit = n_p300_pos_total - p300_pos_with_hit
            p300_neg_no_hit = n_p300_neg_total - p300_neg_with_hit

            # Fisher's exact test
            odds_ratio, p_value = fisher_exact([
                [p300_pos_with_hit, p300_pos_no_hit],
                [p300_neg_with_hit, p300_neg_no_hit]
            ])

            # Score threshold value
            if pct_threshold > 0:
                score_threshold = np.percentile(scores, pct_threshold)
            else:
                score_threshold = scores.min()

            results.append({
                'motif': motif_name,
                'percentile_threshold': pct_threshold,
                'pct_hits_kept': f'top {100-pct_threshold}%',
                'score_threshold': score_threshold,
                'n_hits_passing': len(filtered_regions),
                'n_regions_passing': len(regions_passing),
                'p300_pos_with_hit': p300_pos_with_hit,
                'p300_neg_with_hit': p300_neg_with_hit,
                'odds_ratio': odds_ratio,
                'log2_odds_ratio': np.log2(odds_ratio) if odds_ratio > 0 else np.nan,
                'p_value': p_value
            })

    return pd.DataFrame(results)
