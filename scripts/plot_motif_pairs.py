#!/usr/bin/env python3
"""
Motif pair plotting module with standalone capability.

Can be used as:
1. Import: from plot_motif_pairs import plot_motif_pairs
2. Standalone script: python plot_motif_pairs.py results_dir --out-dir plots
"""

import os, sys, argparse, logging
from typing import Dict, List, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from scipy.stats import gaussian_kde

# ---------------------------
# Data preparation utilities
# ---------------------------

def prepare_pair_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each motif pair, compute summary statistics across all orientations and spacings.
    Returns df with columns for max, std, mean log2fc values.
    """
    summary = df.groupby(['motif1_name', 'motif2_name']).agg({
        'log2_fc_vs_baseline': ['max', 'std', 'mean'],
        'log2_synergy': ['max', 'std', 'mean']
    }).reset_index()
    
    # Flatten column names
    summary.columns = ['motif1_name', 'motif2_name', 
                      'max_log2fc_baseline', 'std_log2fc_baseline', 'mean_log2fc_baseline',
                      'max_synergy', 'std_synergy', 'mean_synergy']
    
    return summary

def prepare_pair_vs_individual_summary(df: pd.DataFrame, motif_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Compute max log2fc vs individual motifs for each pair.
    Returns df with separate columns for vs_motif1 and vs_motif2.
    """
    results = []
    
    for (m1, m2), group in df.groupby(['motif1_name', 'motif2_name']):
        # Find the column names for individual comparisons
        col1 = f'log2_fc_vs_{m1}_alone'
        col2 = f'log2_fc_vs_{m2}_alone'
        
        if col1 in group.columns and col2 in group.columns:
            max_vs_m1 = group[col1].max()
            max_vs_m2 = group[col2].max()
        else:
            max_vs_m1 = np.nan
            max_vs_m2 = np.nan
        
        results.append({
            'motif1_name': m1,
            'motif2_name': m2,
            'max_vs_motif1': max_vs_m1,
            'max_vs_motif2': max_vs_m2
        })
    
    return pd.DataFrame(results)

def get_cluster_order(matrix):
    distance_matrix = pdist(matrix.fillna(0), metric='euclidean')
    linkage_matrix = linkage(distance_matrix, method='ward')
    dendro = dendrogram(linkage_matrix, no_plot=True)
    indices = dendro['leaves']
    motif_names = matrix.index[indices].tolist()

    return motif_names


def get_optimal_configurations(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each motif pair, find the spacing + orientation that yields:
    1. Maximum log2 fold-change vs baseline
    2. Maximum synergy score

    Args:
        df: DataFrame with motif pair results containing columns:
            motif1_name, motif2_name, orientation_pattern, spacing,
            log2_fc_vs_baseline, log2_synergy

    Returns:
        DataFrame with optimal configurations for each pair
    """
    results = []

    for (m1, m2), group in df.groupby(['motif1_name', 'motif2_name']):
        # Average across model folds if present
        if 'model_fold' in group.columns:
            group_avg = group.groupby(['orientation_pattern', 'spacing']).agg({
                'log2_fc_vs_baseline': 'mean',
                'log2_synergy': 'mean'
            }).reset_index()
        else:
            group_avg = group[['orientation_pattern', 'spacing', 'log2_fc_vs_baseline', 'log2_synergy']].copy()

        # Find max log2_fc_vs_baseline
        idx_max_fc = group_avg['log2_fc_vs_baseline'].idxmax()
        max_fc_row = group_avg.loc[idx_max_fc]

        # Find max synergy
        idx_max_syn = group_avg['log2_synergy'].idxmax()
        max_syn_row = group_avg.loc[idx_max_syn]

        results.append({
            'motif1_name': m1,
            'motif2_name': m2,
            # Max log2_fc configuration
            'max_log2fc_orientation': max_fc_row['orientation_pattern'],
            'max_log2fc_spacing': int(max_fc_row['spacing']),
            'max_log2fc_value': max_fc_row['log2_fc_vs_baseline'],
            # Max synergy configuration
            'max_synergy_orientation': max_syn_row['orientation_pattern'],
            'max_synergy_spacing': int(max_syn_row['spacing']),
            'max_synergy_value': max_syn_row['log2_synergy'],
            # Are they the same configuration?
            'same_optimal_config': (
                max_fc_row['orientation_pattern'] == max_syn_row['orientation_pattern'] and
                max_fc_row['spacing'] == max_syn_row['spacing']
            )
        })

    return pd.DataFrame(results)

def create_symmetric_matrix(df: pd.DataFrame, motif_order: List[str], 
                           value_col: str, fill_value: float = np.nan) -> pd.DataFrame:
    """
    Create a symmetric matrix for heatmap plotting.
    
    Args:
        df: DataFrame with motif1_name, motif2_name, and value columns
        motif_order: Ordered list of motif names for matrix axes
        value_col: Column name containing values to plot
        fill_value: Value to use for missing entries
    
    Returns:
        DataFrame with motif names as both index and columns
    """
    matrix = pd.DataFrame(fill_value, index=motif_order, columns=motif_order)
    
    for _, row in df.iterrows():
        m1, m2 = row['motif1_name'], row['motif2_name']
        val = row[value_col]
        
        if m1 in motif_order and m2 in motif_order:
            matrix.loc[m1, m2] = val
            matrix.loc[m2, m1] = val  # Symmetric
    
    return matrix


def create_asymmetric_matrix(df_summary: pd.DataFrame, df_vs_individual: pd.DataFrame,
                            motif_order: List[str]) -> pd.DataFrame:
    """
    Create asymmetric matrix where upper triangle shows max vs motif1 (row)
    and lower triangle shows max vs motif2 (column).
    
    Args:
        df_summary: Summary with baseline comparisons
        df_vs_individual: Summary with individual motif comparisons
        motif_order: Ordered list of motif names
    
    Returns:
        DataFrame for heatmap with asymmetric values
    """
    n = len(motif_order)
    matrix = pd.DataFrame(np.nan, index=motif_order, columns=motif_order)
    
    # Merge the two dataframes
    df_merged = df_summary.merge(df_vs_individual, on=['motif1_name', 'motif2_name'])
    
    for _, row in df_merged.iterrows():
        m1, m2 = row['motif1_name'], row['motif2_name']
        
        if m1 not in motif_order or m2 not in motif_order:
            continue
        
        i1, i2 = motif_order.index(m1), motif_order.index(m2)
        
        # Diagonal: motif1 = motif2
        if i1 == i2:
            matrix.loc[m1, m2] = row['max_vs_motif1']
        # Upper triangle (i1 < i2): pair vs motif1 alone
        elif i1 < i2:
            matrix.loc[m1, m2] = row['max_vs_motif1']
            # Lower triangle: pair vs motif2 alone
            matrix.loc[m2, m1] = row['max_vs_motif2']
        # Lower triangle (i1 > i2): already handled above
        else:
            matrix.loc[m1, m2] = row['max_vs_motif2']
            matrix.loc[m2, m1] = row['max_vs_motif1']

    return matrix

# ---------------------------
# Plotting utilities
# ---------------------------

def plot_individual_motif_insertions(df, motif_order, out_path, color='#006eae'):
    """
    Plot vertically stacked (ridge) density curves for log2_fc_vs_baseline for each motif.
    
    Args:
        df (pd.DataFrame): Full raw results with log2_fc_vs_baseline.
        motif_order (list of str): Motif names in desired plot order (top first).
        out_path (str): Output filename for plot.
        color (str): Color for filled densities.
    """
    for motif in motif_order:
        subset = df.loc[df['motif_name'] == motif]
        n_total = len(subset)
        n_valid = subset['log2_fc_vs_baseline'].notna().sum()
        print(f"{motif}: {n_total} rows, {n_valid} with valid values")

    plt.figure(figsize=(7, 0.5 * len(motif_order)))
    all_values = df['log2_fc_vs_baseline'].dropna()
    xmin = all_values.min() - 0.1
    xmax = all_values.max() + 0.1

    # For vertical ridge stacking, offset curves by idx
    for idx, motif in enumerate(motif_order):
        values = df.loc[df['motif_name'] == motif, 'log2_fc_vs_baseline'].dropna()
        print(f"{motif}: mean={values.mean()}, std={values.std()}")
        offset = idx
        if len(values) < 2:
            v = values.iloc[0]
            plt.plot([v, v], [offset, offset + 1], color=color, alpha=0.7, linewidth=2)
        else:
            kde = gaussian_kde(values)
            x_grid = np.linspace(xmin, xmax, 200)
            density = kde(x_grid)
            density = density / density.max()
            plt.fill_between(x_grid, offset, density + offset, edgecolor=None, color=color, alpha=0.5)

    # X = 0
    plt.axvline(0, color='#435369', linestyle='--', linewidth=1, zorder=10)
    
    # Y ticks and labels (first motif at top)
    plt.yticks(range(len(motif_order)), motif_order)
    plt.ylim(-1, len(motif_order))
    plt.xlabel('log2 fold-change vs baseline', fontsize=12)
    plt.title('Individual motif insertions', fontsize=14, pad=14)

    # Remove grid outline
    ax = plt.gca()
    for i in ['top', 'left', 'right']:
        ax.spines[i].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

def plot_symmetric_heatmap(matrix: pd.DataFrame,
                          title: str,
                          out_path: str,
                          cmap: str = 'RdBu_r',
                          cbar_label: str = 'Value',
                          center_zero: bool = True,
                          annot: bool = False,
                          fmt: str = '.2f',
                          cluster_order: list = None) -> None:
    """
    Plot a symmetric heatmap with proper formatting.
    
    Args:
        matrix: Square DataFrame with values
        title: Plot title
        out_path: Output file path
        cmap: Colormap name
        cbar_label: Colorbar label
        center_zero: Whether to center colormap at 0
        annot: Whether to annotate cells with values
        fmt: Format string for annotations
    """
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Cluster
    if cluster_order:
        matrix = matrix.loc[cluster_order, cluster_order]

    # Create heatmap
    im = ax.imshow(matrix.values, cmap=cmap, aspect='equal')

    # Compute abs_max for colormap centering and annotation colors
    vmin, vmax = np.nanmin(matrix.values), np.nanmax(matrix.values)
    abs_max = max(abs(vmin), abs(vmax)) if np.isfinite([vmin, vmax]).all() else 1.0

    # Center colormap at 0 if requested
    if center_zero and np.isfinite([vmin, vmax]).all():
        im.set_clim(-abs_max, abs_max)

    # Set ticks and labels
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha='right')
    ax.set_yticklabels(matrix.index)

    # Add value annotations if requested
    if annot:
        for i in range(len(matrix.index)):
            for j in range(len(matrix.columns)):
                val = matrix.iloc[i, j]
                if np.isfinite(val):
                    text = ax.text(j, i, f'{val:{fmt}}',
                                 ha="center", va="center",
                                 color="black" if abs(val) < abs_max * 0.5 else "white",
                                 fontsize=8)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, fontsize=12)
    
    # Title and layout
    ax.set_title(title, fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

def plot_asymmetric_heatmap(matrix: pd.DataFrame,
                           title: str,
                           out_path: str,
                           cmap: str = 'RdBu_r',
                           cbar_label: str = 'Value',
                           center_zero: bool = True,
                           cluster_order: list = None) -> None:
    """
    Plot asymmetric heatmap with annotation explaining the asymmetry.
    """
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Cluster
    if cluster_order:
        matrix = matrix.loc[cluster_order, cluster_order]

    # Create heatmap
    im = ax.imshow(matrix.values, cmap=cmap, aspect='equal')
    
    # Center colormap at 0 if requested
    if center_zero:
        vmin, vmax = np.nanmin(matrix.values), np.nanmax(matrix.values)
        if np.isfinite([vmin, vmax]).all():
            abs_max = max(abs(vmin), abs(vmax))
            im.set_clim(-abs_max, abs_max)
    
    # Set ticks and labels
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha='right')
    ax.set_yticklabels(matrix.index)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, fontsize=12)
    
    # Title and layout
    ax.set_title(title, fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

# ---------------------------
# Main orchestrator
# ---------------------------

def plot_motif_pairs(df_pairs: pd.DataFrame,
                    df_singles: pd.DataFrame,
                    motif_dict: Dict[str, str],
                    out_dir: str,
                    heatmap_cmap: str = 'RdBu_r',
                    custom_colors: Optional[List[str]] = None) -> Dict[str, str]:
    """
    Generate all motif pair heatmap plots.
    
    Args:
        df: DataFrame with motif pair results (from process_motif_pair workers)
        motif_dict: Dictionary of {motif_name: motif_seq}
        out_dir: Output directory for plots
        heatmap_cmap: Colormap for heatmaps
        custom_colors: Custom colors (not used currently, for future extension)
    
    Returns:
        Dictionary of output file paths
    """
    # Style
    plt.style.use('default')
    
    # Ensure output directory
    os.makedirs(out_dir, exist_ok=True)

    # Save motif dictionary as table
    motif_table = pd.DataFrame([
        {'motif_name': name, 'motif_seq': seq} 
        for name, seq in motif_dict.items()
    ])
    motif_table_path = os.path.join(out_dir, 'motif_dictionary.tsv')
    motif_table.to_csv(motif_table_path, sep='\t', index=False)
    print(f"Saved motif dictionary to: {motif_table_path}")
    
    # Get ordered list of motif names
    motif_order = sorted(motif_dict.keys())
    
    # Prepare summary data
    print("Preparing summary statistics...")
    df_summary = prepare_pair_summary(df_pairs)
    df_vs_individual = prepare_pair_vs_individual_summary(df_pairs, motif_dict)
    
    # Output paths
    out_paths = {
        'single_motifs': os.path.join(out_dir, 'individual_motif_insertions.pdf'),
        'max_vs_baseline': os.path.join(out_dir, 'pairs_max_log2fc_vs_baseline.pdf'),
        'max_vs_individual': os.path.join(out_dir, 'pairs_max_log2fc_vs_individual.pdf'),
        'max_synergy': os.path.join(out_dir, 'pairs_max_synergy.pdf'),
        'std_vs_baseline': os.path.join(out_dir, 'pairs_std_log2fc_vs_baseline.pdf'),
    }
    
    # 1. Max log2fc vs baseline (symmetric)
    print("Plotting max log2 FC vs baseline...")
    matrix_max_baseline = create_symmetric_matrix(
        df_summary, motif_order, 'max_log2fc_baseline'
    )
    default_order = get_cluster_order(matrix_max_baseline)

    plot_symmetric_heatmap(
        matrix=matrix_max_baseline,
        title='Max log2 fold-change vs baseline\n(all orientations & spacings)',
        out_path=out_paths['max_vs_baseline'],
        cmap=heatmap_cmap,
        cbar_label='Max log2 fold-change vs baseline',
        center_zero=True,
        annot=False,
        cluster_order=default_order
    )

    # 1.5. Individual motif insertions
    print("Plotting individual motif insertions...")
    plot_individual_motif_insertions(
        df=df_singles,
        motif_order=default_order,
        out_path=out_paths['single_motifs']
    )
    
    # 2. Asymmetric: upper = vs motif1, lower = vs motif2
    print("Plotting max log2 FC vs individual motifs (asymmetric)...")
    matrix_vs_individual = create_asymmetric_matrix(
        df_summary, df_vs_individual, motif_order
    )
    plot_asymmetric_heatmap(
        matrix=matrix_vs_individual,
        title='Max log2 fold-change versus individual motifs\n(upper: vs row motif | lower: vs column motif)',
        out_path=out_paths['max_vs_individual'],
        cmap=heatmap_cmap,
        cbar_label='Max log2 fold-change',
        center_zero=True,
        cluster_order=default_order
    )
    
    # 3. Max synergy (symmetric)
    print("Plotting max synergy...")
    matrix_synergy = create_symmetric_matrix(
        df_summary, motif_order, 'max_synergy'
    )
    plot_symmetric_heatmap(
        matrix=matrix_synergy,
        title='Max synergy score\n(pair - additive expectation)',
        out_path=out_paths['max_synergy'],
        cmap=heatmap_cmap,
        cbar_label='Max synergy (log2)',
        center_zero=True,
        annot=False,
        cluster_order=default_order
    )
    
    # 4. Std of log2fc vs baseline (symmetric)
    print("Plotting standard deviation...")
    matrix_std = create_symmetric_matrix(
        df_summary, motif_order, 'std_log2fc_baseline'
    )
    plot_symmetric_heatmap(
        matrix=matrix_std,
        title='Standard deviation log2 fold-change\nacross all orientations & spacings',
        out_path=out_paths['std_vs_baseline'],
        cmap=heatmap_cmap,
        cbar_label='St. dev. log2 fold-change vs. baseline',
        center_zero=True,
        annot=False,
        cluster_order=default_order
    )

    # 5. Optimal configurations report
    print("Computing optimal spacing/orientation configurations...")
    df_optimal = get_optimal_configurations(df_pairs)
    optimal_path = os.path.join(out_dir, 'optimal_configurations.tsv')
    df_optimal.to_csv(optimal_path, sep='\t', index=False)
    out_paths['optimal_configurations'] = optimal_path

    # Print summary
    n_pairs = len(df_optimal)
    n_same_config = df_optimal['same_optimal_config'].sum()
    print(f"\nOptimal Configuration Summary:")
    print(f"  Total motif pairs: {n_pairs}")
    print(f"  Pairs with same optimal config for log2FC and synergy: {n_same_config} ({100*n_same_config/n_pairs:.1f}%)")
    print(f"  Saved optimal configurations to: {optimal_path}")

    # Show top pairs by max log2FC
    print(f"\nTop 5 pairs by max log2 fold-change vs baseline:")
    top_fc = df_optimal.nlargest(5, 'max_log2fc_value')
    for _, row in top_fc.iterrows():
        print(f"  {row['motif1_name']} + {row['motif2_name']}: "
              f"log2FC={row['max_log2fc_value']:.3f} at {row['max_log2fc_orientation']} spacing={row['max_log2fc_spacing']}")

    # Show top pairs by max synergy
    print(f"\nTop 5 pairs by max synergy:")
    top_syn = df_optimal.nlargest(5, 'max_synergy_value')
    for _, row in top_syn.iterrows():
        print(f"  {row['motif1_name']} + {row['motif2_name']}: "
              f"synergy={row['max_synergy_value']:.3f} at {row['max_synergy_orientation']} spacing={row['max_synergy_spacing']}")

    return out_paths

# ============================================================================
# STANDALONE SCRIPT FUNCTIONALITY
# ============================================================================

def parse_standalone_args():
    """Parse command line arguments for standalone mode"""
    parser = argparse.ArgumentParser(
        description="Plot motif pairs results from saved dataframes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("results_dir", type = str,
                       help="path to output of 6.2.motif_pairs.py containing motif_pairs.raw_results.tsv.gz")
    parser.add_argument("--out-dir", default = "plots",
                       help="Name of directory for output plots (will be created in results_dir)")
    
    # Plot customization
    parser.add_argument("--heatmap-cmap", default = "RdBu_r",
                       help="Colormap for heatmap plots")
    parser.add_argument("--custom-colors", nargs="+",
                       help="Custom colors for different motif counts (hex codes or color names)")
    
    return parser.parse_args()

def main_standalone():
    """Main function for standalone script mode"""
    args = parse_standalone_args()

    # Create output directory
    out_dir = os.path.join(args.results_dir, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    
    try:
        # Load dataframe
        pairs_path = os.path.join(args.results_dir, "motif_pairs.raw_results.tsv.gz")
        singles_path = os.path.join(args.results_dir, "individual_motifs.raw_results.tsv.gz")

        df_pairs = pd.read_csv(pairs_path, sep='\t')
        df_singles = pd.read_csv(singles_path, sep='\t')

        # Infer motif_dict
        df_motif = df_singles[['motif_name', 'motif_seq']].drop_duplicates()
        motif_dict = dict(zip(df_motif['motif_name'], df_motif['motif_seq']))
        
        # Parse custom colors if provided
        custom_colors = args.custom_colors if args.custom_colors else None
        
        plot_motif_pairs(
            df_pairs=df_pairs,
            df_singles=df_singles,
            motif_dict=motif_dict,
            out_dir=out_dir,
            heatmap_cmap=args.heatmap_cmap,
            custom_colors=custom_colors
        )
        
    except Exception as e:
        print(f'Error during processing: {e}')
        sys.exit(1)

# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    # This runs when the script is executed directly
    main_standalone()
else:
    # This runs when the module is imported
    pass