#!/usr/bin/env python3
"""
Motif spacing plotting module with standalone capability.

Can be used as:
1. Import: from plot_motif_spacing import plot_motif_spacing
2. Standalone script: python plot_motif_spacing.py input_data.tsv --out-dir plots/
"""

import os, sys, argparse, logging
from typing import Dict, Iterable, List, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from pathlib import Path

# ---------------------------
# Data preparation utilities
# ---------------------------

def prepare_avg_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-(motif_counts, orientation_pattern, spacing) mean/sem of log2_fc_vs_baseline.
    """
    out = (
        df.groupby(['motif_counts', 'orientation_pattern', 'spacing'])
          .agg({'log2_fc_vs_baseline': ['mean', 'sem']})
          .reset_index()
    )
    out.columns = ['motif_counts', 'orientation_pattern', 'spacing', 'log2_fc_mean', 'log2_fc_sem']
    return out


def pivot_all_spacings(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot to (motif_counts, orientation_pattern) x spacing matrix for heatmap of log2_fc_vs_baseline.
    """
    return df.pivot_table(
        values='log2_fc_vs_baseline',
        index=['motif_counts', 'orientation_pattern'],
        columns='spacing',
        aggfunc='mean'
    )

# ---------------------------
# Color utilities
# ---------------------------

def pick_pattern_colors(patterns: Iterable[str],
                        custom_colors: Optional[List[str]] = None) -> Dict[str, Tuple[float, float, float, float]]:
    """
    Return a mapping from orientation pattern -> color.
    If custom_colors supplied and long enough, use it; otherwise pick a small fixed palette
    (for ≤4) or fall back to a colormap.
    """
    patterns = list(dict.fromkeys(patterns))  # stable-unique
    n = len(patterns)

    if custom_colors and len(custom_colors) >= n:
        colors = custom_colors[:n]
    elif n <= 4:
        # small, readable set
        colors = ['#0096a0', '#429130', '#f29742', '#c5373d'][:n]
    else:
        colors = plt.cm.viridis(np.linspace(0, 1, n))

    return dict(zip(patterns, colors))


def model_fold_colors(n_folds: int = 5, cmap_name: str = 'BuGn') -> List:
    """
    Colors for per-fold lines in the variability plot.
    """
    return list(plt.cm.get_cmap(cmap_name)(np.linspace(0, 1, n_folds)))

# ---------------------------
# Tick-label utilities
# ---------------------------

def every_nth_labels(all_values: List, n: int = 10, mode: str = "by_position") -> List[str]:
    """
    Return labels for x-axis with only every nth shown.
    mode:
      - "by_position": label items at indices i where i % n == 0
      - "by_value": label items whose numeric value % n == 0 (works if values are numeric)
    """
    if mode not in {"by_position", "by_value"}:
        mode = "by_position"

    if mode == "by_position":
        return [str(v) if i % n == 0 else '' for i, v in enumerate(all_values)]
    else:  # by_value
        labels = []
        for v in all_values:
            try:
                labels.append(str(v) if (float(v) % n == 0) else '')
            except Exception:
                labels.append('')
        return labels
    
# ---------------------------
# Plotters
# ---------------------------

def plot_heatmap_all_spacings(pivot_data: pd.DataFrame,
                              title: str,
                              out_path: str,
                              heatmap_cmap: str,
                              x_label: str,
                              cbar_label: str,
                              xtick_n: int = 10,
                              xtick_mode: str = "by_value") -> None:
    """
    Heatmap of mean log2_fc_vs_baseline across spacings, indexed by (motif_counts | orientation_pattern).
    """
    fig, ax = plt.subplots(figsize=(16, 8))

    # Image
    im = ax.imshow(pivot_data.values, cmap=heatmap_cmap, aspect='auto')

    # Center colormap at 0 by making clim symmetric
    vmin, vmax = np.nanmin(pivot_data.values), np.nanmax(pivot_data.values)
    abs_max = max(abs(vmin), abs(vmax)) if np.isfinite([vmin, vmax]).all() else 1.0
    im.set_clim(-abs_max, abs_max)

    # Y ticks (row labels)
    y_labels = [f"{idx[0]} | {idx[1]}" for idx in pivot_data.index]
    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels)

    # X ticks/labels
    all_spacings = list(pivot_data.columns)
    ax.set_xticks(range(len(all_spacings)))
    ax.set_xticklabels(every_nth_labels(all_spacings, n=xtick_n, mode=xtick_mode), rotation=0)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    # Titles/labels
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Motif count | Orientation pattern')

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)

def plot_model_variability_by_fold(
    df,
    motif_counts,
    title,
    out_path,
    x_label,
    y_label,
    n_folds: int = 5,
    fold_cmap: str = 'BuGn',
):
    motif_counts = list(motif_counts)

    # Collect patterns per row (motif_count)
    row_patterns = [
        sorted(df[df['motif_counts'] == mc]['orientation_pattern'].dropna().unique())
        for mc in motif_counts
    ]
    n_rows = len(motif_counts)
    max_cols = max((len(p) for p in row_patterns), default=1)

    # Figure sized to the *max* columns, but each row only creates axes for its existing patterns
    fig = plt.figure(figsize=(4 * max(1, max_cols), 3.5 * max(1, n_rows)))
    outer = GridSpec(nrows=n_rows, ncols=1, hspace=0.35, figure=fig)

    # Shared y-limits across all plots
    if df['log2_fc_vs_baseline'].notna().any():
        y_min = float(np.nanmin(df['log2_fc_vs_baseline'].values))
        y_max = float(np.nanmax(df['log2_fc_vs_baseline'].values))
        pad = 0.05 * (y_max - y_min) if np.isfinite(y_max - y_min) and (y_max > y_min) else 1.0
        ylims = (y_min - pad, y_max + pad)
    else:
        ylims = (0, 1)

    # Colors for folds
    colors = model_fold_colors(n_folds, cmap_name=fold_cmap)

    first_ax = None
    for i, (mc, patterns) in enumerate(zip(motif_counts, row_patterns)):
        n_cols = max(1, len(patterns))
        inner = GridSpecFromSubplotSpec(1, n_cols, subplot_spec=outer[i], wspace=0.25)

        if len(patterns) == 0:
            # Graceful fallback if a row has no patterns
            ax = fig.add_subplot(inner[0, 0])
            ax.text(0.5, 0.5, f'No patterns for {mc} motifs', ha='center', va='center')
            ax.set_axis_off()
            continue

        for j, pattern in enumerate(patterns):
            ax = fig.add_subplot(inner[0, j])

            subset = df[(df['motif_counts'] == mc) & (df['orientation_pattern'] == pattern)]
            for fold_idx in range(n_folds):
                model_data = subset[subset['model_fold'] == fold_idx].sort_values('spacing')
                if model_data.empty:
                    continue
                ax.plot(
                    model_data['spacing'],
                    model_data['log2_fc_vs_baseline'],
                    alpha=0.85,
                    lw=1.8,
                    color=colors[fold_idx],
                )

            ax.set_title(f'{mc} motifs, {pattern}')
            ax.set_xlabel(x_label)
            if j == 0:
                ax.set_ylabel(y_label)
            ax.grid(True, alpha=0.3)
            ax.set_ylim(*ylims)

            if first_ax is None:
                first_ax = ax

    # --- single legend for folds, placed near the title ---
    legend_handles = [Line2D([0], [0], color=colors[k], lw=2, label=f'Fold {k}')
                    for k in range(n_folds)]
    
    # Figure title
    fig.suptitle(title, fontsize=16, y=0.94)

    # Legend just under the top edge, centered (right next to the title area)
    fig.legend(handles=legend_handles,
            title=None,
            loc='upper center',
            bbox_to_anchor=(0.5, 0.95),
            fontsize=12,
            ncol=n_folds,
            frameon=False,
            handlelength=2.0)

    fig.suptitle(title, fontsize=16)
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)

def plot_effects_by_motif_count(
    df: pd.DataFrame,
    motif_counts: Iterable[int],
    title: str,
    out_path: str,
    x_label: str,
    metric_and_labels: List[Tuple[str, str]] = (
        ('log2_fc_vs_baseline', 'log2 fold-change vs. background'),
        ('log2_fc_vs_single',   'log2 fold-change vs. single motif')),
    custom_colors: Optional[List[str]] = None
) -> None:
    """
    For each motif_count, plot average metric across spacings with one line per orientation pattern.
    One row per motif count; one column per metric. Y axes are shared within each column.
    Ensures a horizontal y=0 line and that the y-min is never greater than 0.
    """
    motif_counts = list(motif_counts)
    n_cols = len(metric_and_labels)
    n_rows = len(motif_counts)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.5 * n_cols, 3.5 * n_rows), sharey='col')
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    # Pre-compute y-limits per column (metric) and force bottom <= 0
    col_ylims = []
    df_used = df[df['motif_counts'].isin(motif_counts)]
    for metric, _ in metric_and_labels:
        vals = df_used[metric].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            y_min, y_max = 0.0, 1.0
        else:
            y_min, y_max = float(np.min(vals)), float(np.max(vals))
        rng = y_max - y_min
        pad = 0.05 * rng if np.isfinite(rng) and rng > 0 else 0.1
        bottom = min(0.0, y_min - pad)  # never above zero
        top = y_max + pad
        col_ylims.append((bottom, top))

    for row_j, n_motifs in enumerate(motif_counts):
        data_subset = df[df['motif_counts'] == n_motifs]
        color_map = pick_pattern_colors(data_subset['orientation_pattern'].unique(), custom_colors)

        for col_i, (metric, ylabel) in enumerate(metric_and_labels):
            ax = axes[row_j, col_i]

            for pattern in sorted(data_subset['orientation_pattern'].unique()):
                pdat = data_subset[data_subset['orientation_pattern'] == pattern]
                avg = pdat.groupby('spacing')[metric].mean().reset_index()
                if len(avg) == 0:
                    continue
                ax.plot(
                    avg['spacing'],
                    avg[metric],
                    label=pattern,
                    color=color_map.get(pattern, None),
                    linewidth=2
                )

            # Aesthetics
            ax.set_title(f'{n_motifs} motifs')
            ax.set_xlabel(x_label)
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)

        # Legend per row (rightmost panel)
        axes[row_j, -1].legend(title='Orientation',
                               bbox_to_anchor=(1.05, 1),
                               loc='upper left',
                               fontsize=12,
                               frameon=False)

    # Apply y-lims and draw y=0 line for every axis (per column)
    for col_i in range(n_cols):
        bottom, top = col_ylims[col_i]
        for row_j in range(n_rows):
            ax = axes[row_j, col_i]
            ax.set_ylim(bottom, top)
            ax.axhline(0.0, color='0.35', lw=1.2, ls='--', zorder=5)  # zero reference line

    fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)


# ---------------------------
# Orchestrator
# ---------------------------

def plot_motif_spacing(df: pd.DataFrame,
                       motif_name: str,
                       motif_counts: Iterable[int],
                       out_dir: str,
                       heatmap_cmap: str = 'RdBu_r',
                       custom_colors: Optional[List[str]] = None,
                       n_folds: int = 5,
                       xtick_step: int = 10,
                       xtick_mode: str = "by_value") -> Dict[str, str]:
    """
    Orchestrates plotting of:
      1) Heatmap across all spacings
      2) Per-fold model variability traces
      3) Metric comparison lines (baseline vs single)
    Returns dict of output file paths.
    """
    # Style
    plt.style.use('default')

    # Ensure output directory
    os.makedirs(out_dir, exist_ok=True)

    # Output paths
    out_paths = {
        'all_heatmap': os.path.join(out_dir, "all_spacings_heatmap.pdf"),
        'by_fold':     os.path.join(out_dir, "model_variability.pdf"),
        'all_effects': os.path.join(out_dir, "all_spacings_effects.pdf"),
    }

    # Titles and labels
    plot_titles = {
        'all_heatmap': f'{motif_name} motif insertions at all spacings',
        'by_fold':     f'{motif_name} motif insertions by fold',
        'all_effects': f'{motif_name} motif insertions at all spacings',
    }
    label_key = {
        'bp': 'Spacing (bp)',
        'fc_background': 'log2 fold-change vs. background',
        'fc_single': 'log2 fold-change vs. single motif'
    }

    print("Plotting heatmap...")
    pivot_data = pivot_all_spacings(df)
    plot_heatmap_all_spacings(
        pivot_data=pivot_data,
        title=plot_titles['all_heatmap'],
        out_path=out_paths['all_heatmap'],
        heatmap_cmap=heatmap_cmap,
        x_label=label_key['bp'],
        cbar_label=label_key['fc_background'],
        xtick_n=xtick_step,
        xtick_mode=xtick_mode
    )

    print("Plotting per-fold variability...")
    plot_model_variability_by_fold(
        df=df,
        motif_counts=motif_counts,
        title=plot_titles['by_fold'],
        out_path=out_paths['by_fold'],
        x_label=label_key['bp'],
        y_label=label_key['fc_background'],
        n_folds=n_folds,
    )

    print("Plotting metric comparison...")
    plot_effects_by_motif_count(
        df=df,
        motif_counts=motif_counts,
        title=plot_titles['all_effects'],
        out_path=out_paths['all_effects'],
        x_label=label_key['bp'],
        metric_and_labels=[
            ('log2_fc_vs_baseline', label_key['fc_background']),
            ('log2_fc_vs_single',   label_key['fc_single']),
        ],
        custom_colors=custom_colors
    )

    print("All plots saved successfully!")
    return out_paths


# ============================================================================
# STANDALONE SCRIPT FUNCTIONALITY
# ============================================================================

def setup_logging(log_level="INFO"):
    """Setup basic logging to console"""
    level_val = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=level_val,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    return logging.getLogger(__name__)

def validate_dataframe(df, file_path, logger):
    """Validate that dataframe has required columns"""
    required_cols = ['motif_counts', 'orientation_pattern', 'spacing', 
                    'log2_fc_vs_baseline', 'log2_fc_vs_single', 'model_fold']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        logger.error(f"Missing required columns in {file_path}: {missing_cols}")
        return False
    
    if df.empty:
        logger.warning(f"Empty dataframe in {file_path}")
        return False
    
    logger.info(f"Validated {file_path}: {df.shape[0]} rows, {df.shape[1]} columns")
    return True

def load_and_merge_dataframes(input_paths, logger):
    """Load multiple dataframes and merge them"""
    dataframes = []
    
    for path_str in input_paths:
        path = Path(path_str)
        if not path.exists():
            logger.error(f"Input file does not exist: {path}")
            continue
            
        logger.info(f"Loading {path}")
        
        try:
            # Try to determine file format from extension
            if path.suffix.lower() in ['.tsv', '.txt']:
                df = pd.read_csv(path, sep='\t')
            elif path.suffix.lower() == '.csv':
                df = pd.read_csv(path)
            elif path.suffix.lower() in ['.pkl', '.pickle']:
                df = pd.read_pickle(path)
            elif path.suffix.lower() in ['.parquet', '.pq']:
                df = pd.read_parquet(path)
            else:
                # Default to tab-separated
                logger.warning(f"Unknown extension {path.suffix}, trying tab-separated format")
                df = pd.read_csv(path, sep='\t')
                
            # Validate dataframe
            if validate_dataframe(df, path, logger):
                # Add source file info for tracking
                df['source_file'] = str(path.name)
                dataframes.append(df)
            else:
                logger.error(f"Skipping invalid dataframe: {path}")
                
        except Exception as e:
            logger.error(f"Error loading {path}: {e}")
            continue
    
    if not dataframes:
        raise ValueError("No valid dataframes could be loaded")
    
    # Merge all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Check for duplicates (same experimental conditions from multiple files)
    dup_cols = ['motif_counts', 'orientation_pattern', 'spacing', 'model_fold']
    duplicates = merged_df.duplicated(subset=dup_cols, keep=False)
    
    if duplicates.any():
        n_dups = duplicates.sum()
        logger.warning(f"Found {n_dups} potential duplicate rows across input files")
        logger.info("Duplicate detection based on: motif_counts, orientation_pattern, spacing, model_fold")
        
        # Show which files have overlapping conditions
        dup_summary = merged_df[duplicates].groupby('source_file')[dup_cols].nunique()
        logger.info(f"Duplicate summary by source file:\n{dup_summary}")
        
        # Option to remove duplicates (keeping first occurrence)
        logger.info("Keeping first occurrence of duplicate conditions...")
        merged_df = merged_df.drop_duplicates(subset=dup_cols, keep='first')
        logger.info(f"After removing duplicates: {merged_df.shape[0]} rows")
    
    # Log summary of merged data
    logger.info(f"Successfully merged {len(dataframes)} dataframes")
    logger.info(f"Final dataset: {merged_df.shape[0]} rows, {merged_df.shape[1]} columns")
    logger.info(f"Source files: {sorted(merged_df['source_file'].unique())}")
    
    return merged_df

def infer_parameters_from_data(df, logger):
    """Infer reasonable parameters from the data if not provided"""
    motif_counts = sorted(df['motif_counts'].unique())
    spacings = sorted(df['spacing'].unique())
    
    logger.info(f"Inferred parameters:")
    logger.info(f"  Motif counts: {motif_counts}")
    logger.info(f"  Spacing range: {min(spacings)}-{max(spacings)} bp ({len(spacings)} values)")
    
    return motif_counts

def filter_dataframe(df, args, logger):
    """Apply any requested filters to the dataframe"""
    original_rows = len(df)
    
    # Remove +--. which is redundant with ++-
    df = df[df['orientation_pattern'] != "+--"]

    # Filter by spacing range
    if args.min_spacing is not None:
        df = df[df['spacing'] >= args.min_spacing]
        logger.info(f"Filtered to spacing >= {args.min_spacing}: {len(df)} rows")
    
    if args.max_spacing is not None:
        df = df[df['spacing'] <= args.max_spacing]
        logger.info(f"Filtered to spacing <= {args.max_spacing}: {len(df)} rows")
    
    # Filter by orientation patterns
    if args.filter_patterns:
        df = df[df['orientation_pattern'].isin(args.filter_patterns)]
        logger.info(f"Filtered to patterns {args.filter_patterns}: {len(df)} rows")
    
    # Filter by motif counts if specified
    if args.motif_counts:
        df = df[df['motif_counts'].isin(args.motif_counts)]
        logger.info(f"Filtered to motif counts {args.motif_counts}: {len(df)} rows")
    
    if len(df) == 0:
        raise ValueError("No data remaining after applying filters")
    
    if len(df) != original_rows:
        logger.info(f"Applied filters: {original_rows} -> {len(df)} rows ({len(df)/original_rows*100:.1f}%)")
    
    return df

def parse_standalone_args():
    """Parse command line arguments for standalone mode"""
    parser = argparse.ArgumentParser(
        description="Plot motif spacing analysis results from saved dataframes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("input_paths", nargs="+", 
                       help="Path(s) to input dataframe files (TSV, CSV, pickle, or parquet)")
    parser.add_argument("--out-dir", required=True,
                       help="Output directory for plots")
    
    # Optional parameters (will be inferred if not provided)
    parser.add_argument("--motif-name", 
                       help="Name of motif for plot titles")
    parser.add_argument("--motif-counts", nargs="+", type=int,
                       help="Motif counts to include (default: all found in data)")
    
    # Plot customization
    parser.add_argument("--heatmap-cmap", default = "RdBu_r",
                       help="Colormap for heatmap plots")
    parser.add_argument("--custom-colors", nargs="+",
                       help="Custom colors for different motif counts (hex codes or color names)")
    
    # Data filtering/processing
    parser.add_argument("--min-spacing", type=int,
                       help="Minimum spacing value to include in plots")
    parser.add_argument("--max-spacing", type=int,
                       help="Maximum spacing value to include in plots")
    parser.add_argument("--filter-patterns", nargs="+",
                       help="Only include these orientation patterns (e.g., '++' '--' '+-')")
    
    # Output options
    parser.add_argument("--save-merged-data", action="store_true",
                       help="Save the merged dataframe to output directory")
    parser.add_argument("--log-level", default="INFO",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging level")
    
    return parser.parse_args()

def main_standalone():
    """Main function for standalone script mode"""
    args = parse_standalone_args()
    logger = setup_logging(args.log_level)
    
    logger.info("=== Motif Spacing Plotter ===")
    logger.info(f"Input files: {args.input_paths}")
    logger.info(f"Output directory: {args.out_dir}")
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    try:
        # Load and merge dataframes
        df = load_and_merge_dataframes(args.input_paths, logger)
        
        # Apply filters
        df = filter_dataframe(df, args, logger)
        
        # Infer parameters from data if not provided
        motif_counts_data = infer_parameters_from_data(df, logger)
        
        # Use provided parameters or fall back to inferred ones
        motif_counts = args.motif_counts if args.motif_counts else motif_counts_data
        motif_name = args.motif_name if args.motif_name else ""
        
        logger.info(f"Final plotting parameters:")
        logger.info(f"  Motif name: {motif_name}")
        logger.info(f"  Motif counts: {motif_counts}")
        
        # Parse custom colors if provided
        custom_colors = args.custom_colors if args.custom_colors else None
        
        # Save merged data if requested
        if args.save_merged_data:
            merged_file = os.path.join(args.out_dir, "merged_data.tsv")
            df.to_csv(merged_file, sep='\t', index=False)
            logger.info(f"Saved merged dataframe to: {merged_file}")
    
        # Generate plots
        logger.info("Generating plots...")
        plot_motif_spacing(
            df=df,
            motif_name=motif_name,
            motif_counts=motif_counts,
            out_dir=args.out_dir,
            heatmap_cmap=args.heatmap_cmap,
            custom_colors=custom_colors
        )
        
        logger.info(f"Plots saved to: {args.out_dir}")
        logger.info("=== Complete ===")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
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