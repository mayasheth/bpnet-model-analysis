import os, subprocess, argparse, math
import numpy as np
import pandas as pd
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
sns.set_style('whitegrid')

# output files:
    # hits.bed.bgz file with labeled hits
    # hits.tsv with renamed hits
    # peak summary file predicted and true counts, number and list of hits per peak
    # plots:
        # number of hits x true counts, predicted counts
        # something about hit importance

### --- PARSE ARGUMENTS --- ###

parser = argparse.ArgumentParser(description="")
parser.add_argument("--finemo_dir", type=str, help="path to finemo results_directory")
parser.add_argument("--pred", type=str, help="(optional) path to mean prediction output")
parser.add_argument("--overlap_col", type=str, help="column indicating overlap with peaks in pred file (e.g. EP300_peak_overlap)")
parser.add_argument("--motifs", type=str, help="(optional) path to file with motif names and labels")
parser.add_argument("--target_name", type=str, help="target in motif file")
parser.add_argument("--modisco_id", type=str, help="(if motifs argument given) modisco_id in motif file")
parser.add_argument("--out_dir", type=str, help="output directory")
parser.add_argument("--stratify_hits", action='store_true', help="if true, analyze hits stratified by global coeffient")
args = parser.parse_args()


### --- FUNCTIONS --- ###

def df_to_dict(df, key_col, value_col):
    df_uniq = df[[key_col, value_col]].drop_duplicates()
    df_map = dict(zip(df_uniq[key_col], df_uniq[value_col]))

    return df_map
    
def add_label_to_hits(hits, motif_data):
    if motif_data:
        # add index based on original names 
        motif_idx_map = df_to_dict(motif_data, "motif_id", "motif_idx")
        motif_label_map = df_to_dict(motif_data, "motif_id", "motif_name")
        hits["motif_idx"] = [motif_idx_map[m] for m in hits["motif_name"]] # hits "motif_name" = motif_data "motif_id"

         # label: (orig_idx) motif_renamed
        hits["motif_label"] = [f"({motif_idx_map[n]}) {motif_label_map[n]}" for n in hits["motif_name"]] # new name
        hits["motif_rename"] = [motif_label_map[n] for n in hits["motif_name"]]
    else:
        hits["motif_label"] = hits["motif_name"]

    return hits

def save_hits_to_bed(hits, motif_data, out_prefix): 
    tsv_path = out_prefix + ".tsv"
    bed_path = out_prefix + ".bed"

    hits = add_label_to_hits(hits, motif_data)
    hits.to_csv(tsv_path, sep="\t", index=False, header=True)
    
    hits_bed = hits[["chr", "start", "end", "motif_label", "hit_importance", "strand"]]
    bed_path = out_prefix + ".bed"
    hits_bed.to_csv(bed_path, sep="\t", index=False, header=False)
    subprocess.run(["bgzip", "-f", bed_path], check=True)
    subprocess.run(["tabix", "-p", "bed", f"{bed_path}.gz"], check=True)

    print(f"Saved and indexed hits at {bed_path}.gz!")
    return hits

def compile_predictions_and_hits(pred, peaks, hits, motif_data, out_file):
    ## finemo peaks = pred peaks +/- 250 bp

    peaks["chr_start"] = peaks["chr"] + ":" + (peaks["peak_region_start"] - 250).astype(str) 
    peaks = peaks[['peak_id', 'peak_region_start']]

    ## summarize hits per peak
    # read in hits and summarize (on column peak_id)
    hits = add_label_to_hits(hits, motif_data)
    hits_smry = (
        hits.groupby('peak_id').agg(
            n_hits=('motif_name', 'count'),
            motif_list=('motif_label', lambda x: ' | '.join(x)),
            motif_strand_list=('strand', lambda x: ' | '.join(x)),
            sum_hit_coefficient_global=('hit_coefficient_global', 'sum'),
            sum_hit_coefficient=('hit_coefficient', 'sum'),
            sum_hit_importance=('hit_importance', 'sum')).reset_index()
    )
    print(hits_smry.head())

    # merge hits with peaks and fill missing values
    peaks = peaks.merge(hits_smry, on = "peak_id", how = "left")
    print(peaks.head())

    ## add predictions
    if pred is not None:
        # pred = pred.sort_values(by=['chrom', 'start'], ascending=[True, True])
        # pred["peak_id"] = np.arange(len(pred))
        # pred["peak_name"] = pred["chrom"] + ":" + pred["start"].astype(str) + "-" + pred["end"].astype(str)
        peaks = peaks.merge(pred, on = "peak_id")
        print(f"Predictions without hits from finemo: {peaks['n_hits'].isnull().sum()}")

    ## fill NA
    peaks['n_hits'] = peaks['n_hits'].fillna(0).astype(int)
    peaks['sum_hit_coefficient_global'] = peaks['sum_hit_coefficient_global'].fillna(0).astype(float)
    peaks['sum_hit_coefficient'] = peaks['sum_hit_coefficient'].fillna(0).astype(float)
    peaks['sum_hit_importance'] = peaks['sum_hit_importance'].fillna(0).astype(float)
    peaks[['motif_list', 'motif_strand_list']] = peaks[['motif_list', 'motif_strand_list']].fillna("")
    print(peaks.head())

    # save
    peaks.to_csv(out_file, sep="\t", index=False)
    return peaks

def plot_simple_boxplot(x, y, ax,
                   box_col="#b778b3",
                   xlabel="X",
                   ylabel="Y",
                   title=""):
    """
    Boxplot of y values for each unique x.
    """
    sns.set_theme(style="white")

    # Prepare data
    df = pd.DataFrame({"x": x, "y": y})
    sns.boxplot(x="x", y="y", data=df,
                color=box_col,
                ax=ax,
                fliersize=0.5)

    # Clean aesthetics
    ax.grid(False)
    ax.set_xlabel(xlabel, color='black')
    ax.set_ylabel(ylabel, color='black')
    ax.set_title(title, color='black')
    
    # Set all axes elements to black
    ax.tick_params(colors='black')

def plot_scatter(data, xvar, yvar, ax,
                xlabel="X",
                ylabel="Y",
                title="",
                best_fit_col="#792374"):
    """
    Scatter of large n with:
      - fitted regression line (solid)
      - annotated Spearman rho
    """

    sns.set_theme(style="white")
    cmap = mcolors.LinearSegmentedColormap.from_list("MyPurples", ["#ffffff", best_fit_col])
    sns.kdeplot(data, x=xvar, y=yvar, ax=ax, cmap=cmap, fill=True)

    # remove any grid lines
    ax.grid(False)

    # regression
    x = data[xvar]
    y = data[yvar]
    finite = np.isfinite(x) & np.isfinite(y)

    rho, _ = spearmanr(x[finite], y[finite])
    stats_txt = f"Spearman ρ = {rho:.2f}"
    ax.text(0.04, 0.9, stats_txt,
            transform=ax.transAxes,
            va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor="white", alpha=0.7),
            fontsize=12)

    #ax.set_aspect('equal', 'box')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)

def plot_hits_distribution(n_hits, ax, mn = 0, mx = 15, title = "Hits per peak", color = "#b778b3"):
    ## 1) Distribution of hits per peak
    mn = 0
    mx = 15
    bins = np.arange(mn, mx+2)

    sns.histplot(
        n_hits,
        kde=False,
        bins = bins,
        color=color,
        shrink=0.8,      # controls bar width relative to bin size
        discrete = True,
        ax=ax
    )

    # remove any remaining grid
    ax.grid(False)
    ax.set_xlabel('Motif hits per peak', color='black')
    ax.set_ylabel('Number of peaks', color='black')
    ax.set_title(title, color='black')

def plot_motif_cdfs(hits: pd.DataFrame, output_path: str, motif_col="motif_rename", overlap_col=None):
    """
    For each motif (plus one row for 'All motifs'), plot three CDFs:
      - hit_coefficient (log10 scale)
      - hit_coefficient_global (log10 scale)
      - hit_importance (log10 scale)
      - hit_similarity (log10 scale)

    Each motif gets one row of plots (3 columns), with a motif title on the left.

    If overlap_col is provided (0/1):
      Each subplot shows ECDFs for:
        - all hits (dark)
        - hits_plus  (overlap_col == 1; medium)
        - hits_minus (overlap_col == 0; light)

    Vertical dashed lines mark the 20/40/60/80/100% quantiles of the distribution
    for the subset shown in that row.
    """
    sns.set_theme(style="whitegrid")

    # Colors
    color_plus   = "#a64791"  # dark
    color_all    = "#006eae"  # medium
    color_minus  = "#96a0b3"  # light

    motif_names = ["All motifs"] + sorted(hits[motif_col].unique())
    n_motifs = len(motif_names)

    ncols = 4  
    nrows = n_motifs

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=(2.75 * ncols, 2.25 * nrows),
        squeeze=False
    )

    # --- Precompute global x-axis ranges (log10) for each metric ---
    def safe_log10(x):
        x = np.asarray(x)
        return np.log10(x[x > 0]) if np.any(x > 0) else np.array([])

    ranges = {}
    for metric in ["hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"]:
        vals = safe_log10(hits[metric].dropna().values)
        if vals.size > 0:
            ranges[metric] = (vals.min(), vals.max())
        else:
            ranges[metric] = (0, 1)  # fallback

    # Helper to draw ECDFs with quantile lines
    def draw_ecdfs(ax, df_all, series, log=False, xlabel="", base_color=color_all):
        # Filter positive for log plots
        if log:
            df_all = df_all[df_all[series] > 0].copy()
            if df_all.empty:
                ax.set_title(xlabel, color="black")
                ax.text(0.5, 0.5, "No positive values",
                        ha="center", va="center", fontsize=7,
                        transform=ax.transAxes)
                ax.grid(False)
                return
            xcol = f"log10_{series}"
            df_all[xcol] = np.log10(df_all[series])
        else:
            xcol = series

        # Always plot "all"
        if len(df_all) > 0:
            sns.ecdfplot(data=df_all, x=xcol, ax=ax,
                         color=color_all, linewidth=2, label="all")
        else:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=9, transform=ax.transAxes)

        # If overlap splits provided, overlay plus/minus
        if overlap_col is not None and overlap_col in df_all.columns:
            df_plus = df_all[df_all[overlap_col] == 1]
            df_minus = df_all[df_all[overlap_col] != 1]

            if len(df_plus) > 0:
                sns.ecdfplot(data=df_plus, x=xcol, ax=ax,
                             color=color_plus, linewidth=1.8, label="hits_plus")
            if len(df_minus) > 0:
                sns.ecdfplot(data=df_minus, x=xcol, ax=ax,
                             color=color_minus, linewidth=1.8, label="hits_minus")

        ax.set_xlabel(xlabel, color="black")
        ax.set_ylabel("Proportion", color="black")
        ax.grid(False)

        # --- Vertical markers at 20, 40, 60, 80, 100% quantiles (subset-specific) ---
        if len(df_all) > 0:
            qvals = df_all[xcol].quantile([0.2, 0.4, 0.6, 0.8, 1.0])
            for q in qvals:
                ax.axvline(x=q, color="#c5cad7", linestyle="--", linewidth=0.7)

        # --- Apply consistent x-axis range ---
        if log and series in ranges:
            ax.set_xlim(ranges[series])
        elif not log and series in ranges:
            ax.set_xlim(ranges[series])

    # --- Draw plots motif by motif (including All motifs) ---
    for idx, motif in enumerate(motif_names):
        if motif == "All motifs":
            subset = hits.copy()
        else:
            subset = hits[hits[motif_col] == motif].copy()

        # Column 1: coefficient
        ax = axes[idx][0]
        draw_ecdfs(ax, subset, series="hit_coefficient", log=True,
                   xlabel="Local hit coef. (log10)")

        # Column 2: global coefficient
        ax = axes[idx][1]
        draw_ecdfs(ax, subset, series="hit_coefficient_global", log=True,
                   xlabel="Global hit coef. (log10)")

        # Column 3: importance
        ax = axes[idx][2]
        draw_ecdfs(ax, subset, series="hit_importance", log=True,
                   xlabel="Hit importance (log10)")
        
        # Column 4: similarity
        ax = axes[idx][3]
        draw_ecdfs(ax, subset, series="hit_similarity", log=True,
                   xlabel="Hit similarity (log10)")

        # Add motif heading across row (left margin)
        fig.text(0.01, 1 - ((idx + 0.5) / nrows), motif,
                 rotation=90, va='center', ha='center',
                 fontsize=12 if motif != "All motifs" else 14,
                 color='black', weight="bold" if motif == "All motifs" else "normal")

        # Legend only for first row
        if motif == "All motifs":
            handles, labels = axes[idx][0].get_legend_handles_labels()
            if labels:
                axes[idx][0].legend(handles, labels,
                                    loc="upper left", frameon=False, fontsize=8)

    fig.tight_layout(rect=[0.05, 0, 1, 1])
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

def _prepare_metrics(df: pd.DataFrame, metrics, log10_metrics=True) -> pd.DataFrame:
    """
    Return a clean matrix of metrics for correlation.
    If log10_metrics is True, log10-transform strictly-positive values, dropping others.
    """
    metrics = list(metrics)
    M = df[metrics].copy()
    if log10_metrics:
        for m in metrics:
            # keep only strictly-positive for log10; drop rows that fail
            pos = M[m] > 0
            M = M[pos]
            M[m] = np.log10(M[m])
    # drop rows with any NaNs left
    M = M.dropna(axis=0, how="any")
    return M

def _corr_mats_for_subset(df_subset: pd.DataFrame, metrics, log10_metrics=True):
    """
    Compute Pearson and Spearman correlation matrices for the given subset.
    Returns (pearson_df, spearman_df) or (None, None) if insufficient data.
    """
    M = _prepare_metrics(df_subset, metrics, log10_metrics=log10_metrics)
    # need at least 3 rows to get stable correlation
    if M.shape[0] < 3:
        return None, None
    pearson  = M.corr(method="pearson")
    spearman = M.corr(method="spearman")
    return pearson, spearman

def compute_correlations_by_motif(
    hits: pd.DataFrame,
    motif_col: str = "motif_rename",
    metrics=("hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"),
    log10_metrics: bool = True
):
    """
    Build dicts: {subset_name -> corr_df} for Pearson and Spearman.
    Includes 'All motifs' and each individual motif.
    """
    subsets = ["All motifs"] + sorted(hits[motif_col].dropna().unique().tolist())
    pearson_by_subset  = {}
    spearman_by_subset = {}

    for name in subsets:
        sub_df = hits if name == "All motifs" else hits[hits[motif_col] == name]
        p, s = _corr_mats_for_subset(sub_df, metrics, log10_metrics=log10_metrics)
        pearson_by_subset[name]  = p
        spearman_by_subset[name] = s
    return pearson_by_subset, spearman_by_subset

ABBREV = {
    "hit_coefficient": "hc",
    "hit_coefficient_global": "hcg",
    "hit_importance": "hi",
    "hit_similarity": "hs",      # or "hit_correlation": "hcorr"
}

def _combine_upper_pearson_lower_spearman(P: pd.DataFrame, S: pd.DataFrame, abbrev=True):
    """Return a matrix with Pearson in upper triangle, Spearman in lower."""
    if P is None or S is None:
        return None
    # Ensure same order
    metrics = list(P.index)
    M = P.copy()
    # fill lower with Spearman
    for i in range(len(metrics)):
        for j in range(i):
            M.iloc[i, j] = S.iloc[i, j]
    # optional short names
    if abbrev:
        new_index = [ABBREV.get(m, m) for m in M.index]
        new_cols  = [ABBREV.get(m, m) for m in M.columns]
        M.index = new_index
        M.columns = new_cols
    return M

def plot_corr_compact_triangle(
    df: pd.DataFrame,
    motif_col: str,
    metrics,
    log10_metrics: bool,
    output_path: str,
    col_wrap: int = 5,          # how many panels per row
    cmap: str = "mako_r",
    annot: bool = True,
    annot_fmt: str = ".2f",
    labelsize: int = 8,
):
    # compute Pearson/Spearman per subset (uses your existing helper)
    pearson_by_subset, spearman_by_subset = compute_correlations_by_motif(
        df, motif_col, metrics, log10_metrics
    )

    # Build display list (skip subsets with insufficient data)
    items = []
    for name in pearson_by_subset.keys():
        P, S = pearson_by_subset[name], spearman_by_subset[name]
        M = _combine_upper_pearson_lower_spearman(P, S, abbrev=True)
        if M is not None:
            items.append((name, M))

    n = len(items)
    cols = min(col_wrap, n) if n > 0 else 1
    rows = int(np.ceil(n / cols)) if n > 0 else 1

    # Figure sizing: each panel small & square; single shared cbar on the right
    fig_w = 2.6 * cols + 0.5   # add room for shared colorbar
    fig_h = 2.6 * rows
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(rows, cols + 1,
                           width_ratios=[1]*cols + [0.06],
                           wspace=0.25, hspace=0.35)

    # keep your chosen bounds; if you want full range use vmin, vmax = -1, 1
    vmin, vmax = 0, 1.0
    center = 0 if (vmin < 0 < vmax) else None

    cax = fig.add_subplot(gs[:, -1])  # shared colorbar axis
    first_ax = None

    for k, (name, M) in enumerate(items):
        r, c = divmod(k, cols)
        ax = fig.add_subplot(gs[r, c])
        first_ax = first_ax or ax

        # ---- custom annotation matrix with blank diagonal ----
        if annot:
            # format numbers as strings per annot_fmt
            fmt_func = (lambda x: format(x, annot_fmt))
            annot_vals = np.vectorize(fmt_func)(M.values)
            np.fill_diagonal(annot_vals, "")  # remove diagonal labels
            heatmap_annot = annot_vals
            heatmap_fmt = ""                   # using our own strings
        else:
            heatmap_annot = False
            heatmap_fmt = annot_fmt

        sns.heatmap(
            M, ax=ax, vmin=vmin, vmax=vmax, center=center, cmap=cmap,
            square=True, cbar=False, linewidths=0.5, linecolor="white",
            annot=heatmap_annot, fmt=heatmap_fmt,
            annot_kws={"size": labelsize-1} if annot else None
        )
        ax.set_title(name, fontsize=labelsize+1, pad=4)
        ax.tick_params(axis="both", labelsize=labelsize, length=0)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha="center")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va="center")

    # Shared colorbar
    if first_ax is not None:
        mappable = first_ax.collections[0]
        cb = fig.colorbar(mappable, cax=cax)
        cb.set_label("r (Pearson, upper) / ρ (Spearman, lower)", fontsize=labelsize)
        cb.ax.tick_params(labelsize=labelsize)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

def _align_order(A: pd.DataFrame, B: pd.DataFrame):
    """Make sure A and B have identical row/col order (intersection only)."""
    if A is None or B is None:
        return None, None
    common = [m for m in A.index if m in B.index]
    A2 = A.loc[common, common]
    B2 = B.loc[common, common]
    return A2, B2

def _combine_upper_pearson_lower_spearman(Mu: pd.DataFrame, Ml: pd.DataFrame, abbrev=True):
    """Upper triangle from Mu, lower from Ml."""
    if Mu is None or Ml is None:
        return None
    # enforce same order
    metrics = list(Mu.index)
    M = Mu.copy()
    for i in range(len(metrics)):
        for j in range(i):
            M.iloc[i, j] = Ml.iloc[i, j]
    if abbrev:
        M.index  = [ABBREV.get(m, m) for m in M.index]
        M.columns= [ABBREV.get(m, m) for m in M.columns]
    return M

def compute_delta_correlations_by_motif(
    hits_A: pd.DataFrame,
    hits_B: pd.DataFrame,
    motif_col: str = "motif_rename",
    metrics=("hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"),
    log10_metrics: bool = True,
    diff: str = "B_minus_A",   # or "A_minus_B"
):
    """Δ correlation matrices (Pearson & Spearman) per subset present in either A or B."""
    P_A, S_A = compute_correlations_by_motif(hits_A, motif_col, metrics, log10_metrics)
    P_B, S_B = compute_correlations_by_motif(hits_B, motif_col, metrics, log10_metrics)

    # Use union of subset names; only keep those where both sides have matrices
    subset_names = list(dict.fromkeys(list(P_A.keys()) + list(P_B.keys())))
    dP, dS = {}, {}

    sign = +1 if diff == "B_minus_A" else -1
    for name in subset_names:
        A_P, B_P = _align_order(P_A.get(name), P_B.get(name))
        A_S, B_S = _align_order(S_A.get(name), S_B.get(name))

        if A_P is None or B_P is None or A_S is None or B_S is None:
            dP[name], dS[name] = None, None
            continue
        # require some finite values
        if A_P.isna().all().all() or B_P.isna().all().all():
            dP[name], dS[name] = None, None
            continue

        dP[name] = (B_P - A_P) * sign
        dS[name] = (B_S - A_S) * sign

    return dP, dS

def plot_corr_delta_triangle(
    delta_pearson_by_subset: dict,
    delta_spearman_by_subset: dict,
    output_path: str,
    col_wrap: int = 4,
    cmap: str = "icefire",
    annot: bool = True,
    annot_fmt: str = ".2f",
    labelsize: int = 8,
):
    """One panel/subset: upper=ΔPearson, lower=ΔSpearman. Shared symmetric colorbar."""
    # build matrices for plotting
    items = []
    for name in delta_pearson_by_subset.keys():
        dP = delta_pearson_by_subset.get(name)
        dS = delta_spearman_by_subset.get(name)
        if dP is None or dS is None:
            continue
        M = _combine_upper_pearson_lower_spearman(dP, dS, abbrev=True)
        if M is not None:
            items.append((name, M))

    if not items:
        raise ValueError("No subsets with delta correlations to plot.")

    # symmetric range across all panels (comparable colors)
    all_vals = np.concatenate([x[1].values.ravel() for x in items])
    all_vals = all_vals[np.isfinite(all_vals)]
    maxabs = float(np.nanmax(np.abs(all_vals))) if all_vals.size else 1.0
    #vmin, vmax, center = -maxabs, +maxabs, 0.0
    vmin, vmax, center = -0.2, 0.2, 0

    n = len(items)
    cols = min(col_wrap, n)
    rows = int(np.ceil(n / cols))

    fig_w = 2.6 * cols + 0.6
    fig_h = 2.6 * rows
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(rows, cols + 1,
                           width_ratios=[1]*cols + [0.06],
                           wspace=0.25, hspace=0.35)

    cax = fig.add_subplot(gs[:, -1])
    first_ax = None

    for k, (name, M) in enumerate(items):
        r, c = divmod(k, cols)
        ax = fig.add_subplot(gs[r, c])
        first_ax = first_ax or ax

        # custom annotation with blank diagonal (keep color)
        if annot:
            fmt_func = (lambda x: format(x, annot_fmt))
            annot_vals = np.vectorize(fmt_func)(M.values)
            np.fill_diagonal(annot_vals, "")
            hm_annot, hm_fmt = annot_vals, ""
        else:
            hm_annot, hm_fmt = False, annot_fmt

        sns.heatmap(
            M, ax=ax, vmin=vmin, vmax=vmax, center=center, cmap=cmap,
            square=True, cbar=False, linewidths=0.5, linecolor="white",
            annot=hm_annot, fmt=hm_fmt,
            annot_kws={"size": labelsize-1} if annot else None
        )
        ax.set_title(name, fontsize=labelsize+1, pad=4)
        ax.tick_params(axis="both", labelsize=labelsize, length=0)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha="center")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va="center")

    if first_ax is not None:
        mappable = first_ax.collections[0]
        cb = fig.colorbar(mappable, cax=cax)
        cb.set_label("Δ correlation (upper: ΔPearson, lower: ΔSpearman)", fontsize=labelsize)
        cb.ax.tick_params(labelsize=labelsize)

    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig)

    
### --- SET UP --- ###

# load input files
hits = pd.read_table(os.path.join(args.finemo_dir, "hits.tsv"))
peaks = pd.read_table(os.path.join(args.finemo_dir, "peaks_qc.tsv")) # peak_id, chr, peak_region_start

if args.pred:
    pred = pd.read_table(args.pred, compression = "gzip")
else:
    pred = None

if args.motifs:
    motif_data = pd.read_table(args.motifs)
    motif_data = motif_data.loc[(motif_data["target"] == args.target_name) & (motif_data["modisco_id"] == args.modisco_id)]
else:
    motif_data = None

QUANTILE_METRIC = "hit_coefficient_global"
if args.stratify_hits:
    QUANTILES = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00]
else:
    QUANTILES = [0, 1]

# define output files
os.makedirs(args.out_dir, exist_ok = True)
peak_smry_path = os.path.join(args.out_dir, "hits_per_peak.with_predictions.tsv.gz")
hits_renamed_prefix = os.path.join(args.out_dir, "hits_renamed")
per_peak_plots = os.path.join(args.out_dir, "hits_per_peak.pdf")
per_motif_plots = os.path.join(args.out_dir, "motif_scores.pdf")

### --- REFORMAT --- ###
peak_smry = compile_predictions_and_hits(pred, peaks, hits, motif_data, peak_smry_path)
hits_renamed = save_hits_to_bed(hits, motif_data, hits_renamed_prefix)

## save hits per bin
if args.stratify_hits:
    qs = np.quantile(hits[QUANTILE_METRIC].to_numpy(), QUANTILES)
    for i in range(len(QUANTILES) - 1):
        lo, hi = qs[i], qs[i+1]
        mask   = (hits[QUANTILE_METRIC] >= lo) & (hits[QUANTILE_METRIC] <= hi)
        bin_hits = hits.loc[mask].copy()

        # suffix like ".q00_20", ".q20_40", … ".q80_100"
        quantile_id = f'.q{int(QUANTILES[i]*100):02d}_{int(QUANTILES[i+1]*100):02d}'
        save_hits_to_bed(bin_hits, motif_data, hits_renamed_prefix + quantile_id)

### ---  PLOT PER PEAK METRICS --- ###

# Fixed axes helper
def _apply_axes_limits(ax_hits, ax_true, ax_pred, ax_delta):
    # x for all panels
    for ax in (ax_hits, ax_true, ax_pred, ax_delta):
        ax.set_xlim(-0.5, 10.5); ax.set_xticks(range(0, 11))

    # y-lims (counts / logcounts / delta)
    ax_hits.set_ylim(0, 37000)
    ax_true.set_ylim(0, 15); ax_true.set_yticks(range(0, 16, 3))
    ax_pred.set_ylim(0, 15); ax_pred.set_yticks(range(0, 16, 3))
    ax_delta.set_ylim(-8, 8); ax_delta.set_yticks(range(-9, 9, 3))
    ax_delta.axhline(0, ls="--", color="#c5cad7", alpha=1)

if args.pred:
    peak_smry["delta_logcounts"] = peak_smry["mean_pred_logcounts"] - peak_smry["true_logcounts"]

    # top-k% thresholds: top20, top40, top60, top80, top100
    tops = QUANTILES[1:]
    nrows = len(tops)
    fig, axes = plt.subplots(nrows, 4, figsize=(10, 2.5 * nrows), squeeze=False)

    vals = hits[QUANTILE_METRIC].to_numpy()
    for r, top in enumerate(tops):
        thr = np.quantile(vals, 1.0 - top)   # top-k% means keep >= this threshold
        keep = hits[QUANTILE_METRIC] >= thr
        hits_top = hits.loc[keep].copy()

        # recompute per-peak summary using only this subset of hits
        peak_smry_top = compile_predictions_and_hits(pred, peaks, hits_top, motif_data, None)
        peak_smry_top["delta_logcounts"] = (
            peak_smry_top["mean_pred_logcounts"] - peak_smry_top["true_logcounts"]
        )

        # row title
        quantile_label = f"Top {int(top*100)}% (n={hits_top.shape[0]})"
        # 1) Distribution of hits per peak
        plot_hits_distribution(peak_smry_top["n_hits"], axes[r, 0], title = quantile_label)

        # 2) Observed counts x #hits
        plot_simple_boxplot(
            peak_smry_top["n_hits"],
            peak_smry_top["true_logcounts"],
            ax=axes[r, 1],
            xlabel="Motif hits per peak",
            ylabel="Observed log counts",
            title=None,
        )

        # 3) Predicted counts x #hits
        plot_simple_boxplot(
            peak_smry_top["n_hits"],
            peak_smry_top["mean_pred_logcounts"],
            ax=axes[r, 2],
            xlabel="Motif hits per peak",
            ylabel="Predicted log counts",
            title=None,
        )

        # 4) Delta logcounts x #hits
        plot_simple_boxplot(
            peak_smry_top["n_hits"],
            peak_smry_top["delta_logcounts"],
            ax=axes[r, 3],
            xlabel="Motif hits per peak",
            ylabel="Predicted - observed log counts",
            title=None,
        )

        _apply_axes_limits(axes[r, 0], axes[r, 1], axes[r, 2], axes[r, 3])

    fig.subplots_adjust(wspace=0.35, hspace=0.55)
    sns.despine(trim=True)
    plt.savefig(per_peak_plots)

# no predictions available
else:
    tops = QUANTILES[1:]
    ncols = len(tops)
    fig, axes = plt.subplots(1, ncols, figsize=(3.5 * ncols, 4))

    vals = hits[QUANTILE_METRIC].to_numpy()
    for c, top in enumerate(tops):
        thr = np.quantile(vals, 1.0 - top)   # top-k% means keep >= this threshold
        keep = hits[QUANTILE_METRIC] >= thr
        hits_top = hits.loc[keep].copy()

        # recompute per-peak summary using only this subset of hits
        peak_smry_top = compile_predictions_and_hits(pred, peaks, hits_top, motif_data, None)
        peak_smry_top["delta_logcounts"] = (
            peak_smry_top["mean_pred_logcounts"] - peak_smry_top["true_logcounts"]
        )

        # row title
        quantile_label = f"Top {int(top*100)}% (n={hits_top.shape[0]})"

        # distribution of hits per peak
        plot_hits_distribution(peak_smry_top["n_hits"], axes[c], title = quantile_label)
        axes[c].set_xlim(-0.5, 10.5)
        axes[c].set_xticks(range(0, 11))
        axes[c].set_ylim(0, 37000)

    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(per_peak_plots)

### --- PLOT PER MOTIF METRICS, incl correlations --- ###
if args.pred:
    hits_renamed = hits_renamed.merge(pred, on = "peak_id")
    print(hits_renamed.columns)

if args.motifs:
    cols_iterate = ["motif_rename", "motif_label", "motif_family"]
else:
    cols_iterate = ["motif_name"]

for col_name in cols_iterate:
    if col_name in hits_renamed.columns:
        this_out_file = os.path.join(args.out_dir, f"motif_scores.{col_name}.pdf")
        plot_motif_cdfs(hits_renamed, this_out_file, col_name, args.overlap_col)

    plot_corr_compact_triangle(
        df=hits_renamed,
        output_path=os.path.join(args.out_dir, f"correlations.all_elements.{col_name}.pdf"),
        motif_col= col_name,
        metrics=("hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"),
        log10_metrics=True  
    )

    if (args.pred is not None) & (args.overlap_col is not None):
        hits_renamed_plus = hits_renamed.loc[hits_renamed[args.overlap_col] == 1]
        hits_renamed_minus = hits_renamed.loc[hits_renamed[args.overlap_col] != 1]

        plot_corr_compact_triangle(
            df=hits_renamed_plus,
            output_path=os.path.join(args.out_dir, f"correlations.peak_overlap.{col_name}.pdf"),
            motif_col= col_name,
            metrics=("hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"),
            log10_metrics=True  
        )

        plot_corr_compact_triangle(
            df=hits_renamed_minus,
            output_path=os.path.join(args.out_dir, f"correlations.no_peak_overlap.{col_name}.pdf"),
            motif_col= col_name,
            metrics=("hit_coefficient", "hit_coefficient_global", "hit_importance", "hit_similarity"),
            log10_metrics=True  
        )

        # delta correlations
        dP, dS = compute_delta_correlations_by_motif(
            hits_renamed_plus, hits_renamed_minus,
            motif_col=col_name,
            metrics=("hit_coefficient","hit_coefficient_global","hit_importance","hit_similarity"),
            log10_metrics=True,
            diff="A_minus_B"   # or "A_minus_B"
        )
        plot_corr_delta_triangle(dP, dS,
            output_path=os.path.join(args.out_dir, f"correlations.delta_peaks_no_peaks.{col_name}.pdf"),
            col_wrap=5)

