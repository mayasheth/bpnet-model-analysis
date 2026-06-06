import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import subprocess
from scipy.stats import pearsonr, spearmanr


def plot_scatter(df, x_var, y_var,
                ax,
                marker_size=15,
                alpha=0.2,
                point_col="#a64791",
                best_fit_col="#792374",
                one_to_one_col="#000000",
                xlabel="X",
                ylabel="Y",
                title="",
                show_legend=True):
    """
    Scatter of large n with:
      - 1:1 reference line (dashed)
      - fitted regression line (solid)
      - equal aspect so y=x is 45°
      - annotated Pearson r and Spearman rho
      - no grid, frameless legend
    """
    x = df[x_var]
    y = df[y_var]

    sns.set_theme(style="white")
    cmap = mcolors.LinearSegmentedColormap.from_list("MyPurples", ["#ffffff", best_fit_col])
    sns.kdeplot(data=df, x=x_var, y=y_var, ax=ax, cmap=cmap, fill=True)

    # ax.scatter(x, y,
    #            s=marker_size,
    #            alpha=alpha,
    #            linewidth=0,
    #            color=point_col,
    #            rasterized=True)

    # remove any grid lines
    ax.grid(False)

    # compute limits
    finite = np.isfinite(x) & np.isfinite(y)
    lo = 0 
    hi = max(np.nanmax(x[finite]), np.nanmax(y[finite]))
    hi = max(hi, 20)

    # 1:1 ref
    ref_line, = ax.plot([lo, hi], [lo, hi],
                        color=one_to_one_col, linestyle=(0, (5, 5)), lw=1.5, label="y = x")
    # regression
    m, b = np.polyfit(x[finite], y[finite], 1)
    reg_y = m * np.array([lo, hi]) + b
    reg_line, = ax.plot([lo, hi], reg_y,
                        color=best_fit_col, linestyle="-", lw=1.5, label="fit")

    r, _ = pearsonr(x[finite], y[finite])
    rho, _ = spearmanr(x[finite], y[finite])
    stats_txt = f"Pearson r = {r:.2f}\nSpearman ρ = {rho:.2f}"
    ax.text(0.04, 0.9, stats_txt,
            transform=ax.transAxes,
            va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor="white", alpha=0.7),
            fontsize=12)

    ax.set_aspect('equal', 'box')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    ax.legend(handles=[ref_line, reg_line], loc="lower right", frameon=False)

    if not show_legend:
        ax.legend_.remove()

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
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


# file paths
BASE_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model"
CHROM_ANNOT_FILE = "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"

MEAN_OUT_DIR = os.path.join(BASE_DIR, "predictions_mean", "all_folds"); os.makedirs(MEAN_OUT_DIR, exist_ok = True)
CV_OUT_DIR = os.path.join(BASE_DIR, "predictions_cv", "all_folds"); os.makedirs(CV_OUT_DIR, exist_ok = True)

FOLDS = ["fold0", "fold1", "fold2", "fold3", "fold4"]
NARROWPEAK_SCHEMA = ["chrom", "start", "end", "peak_name", "peak_score", "peak_strand", "peak_signal", "peak_pval", "peak_qval", "peak_summit"]

#### MEAN PREDICTIONS ####
mean_pred_file = os.path.join(MEAN_OUT_DIR, "mean_predictions_counts.tsv.gz")
mean_pred = pd.read_table(mean_pred_file, compression = "gzip")
print(mean_pred.head())
print(len(mean_pred))

cols_keep = ["chrom", "start", "end", "EP300.RPM", "EP300_peak_overlap"]
chrom_annot = pd.read_table(CHROM_ANNOT_FILE, usecols = cols_keep)
print(chrom_annot.head())
print(len(chrom_annot))

mean_pred = mean_pred.merge(chrom_annot, on = ["chrom", "start", "end"], how = "left")
print(np.count_nonzero(mean_pred['EP300_peak_overlap'].isna()))

#plot
fig, axes = plt.subplots(1, 2, figsize=(8, 4))

plot_scatter(mean_pred,
            "mean_pred_logcounts", "true_logcounts",
            ax = axes[0],
            xlabel=f"Predicted log counts",
            ylabel="Observed log counts",
            title=f"EP300 counts (mean predictions)")
plot_scatter(mean_pred[mean_pred["EP300_peak_overlap"] == 1],
            "mean_pred_logcounts", "true_logcounts",
            ax = axes[1],
            xlabel=f"Predicted log counts",
            ylabel="Observed log counts",
            title=f"EP300 counts (mean predictions, p300 peaks)")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(MEAN_OUT_DIR, "pred_vs_observed.pdf"))


#### CV PREDICTIONS ####
# goal: combine predictions across folds and plot performance versus observed counts (for CV)
cv_pred_list = []

for this_fold in FOLDS: 
    print(f"Working on {this_fold}...")

    peaks = pd.DataFrame()

    # load true and predicted counts
    for col in ["true_logcounts_0", "true_logcounts_1", "pred_logcounts_0", "pred_logcounts_1"]:
        peaks[col] = np.load(os.path.join(BASE_DIR, "predictions_cv", this_fold, f"{col}.npy"))
    
    peaks["fold"] = this_fold
    peaks["true_logcounts_total"] = peaks["true_logcounts_0"] + peaks["true_logcounts_1"]
    peaks["pred_logcounts_total"] = peaks["pred_logcounts_0"] + peaks["pred_logcounts_1"]
    
    # print percent of 0 observed counts
    n_zero = np.count_nonzero(peaks["true_logcounts_total"] == 0)
    print(f"{n_zero} zeros out of {len(peaks)} elements")

    cv_pred_list.append(peaks)

# combine formatted peaks, add other annotations
cv_pred = pd.concat(cv_pred_list, axis = 0, ignore_index = True)
cv_pred.to_csv(os.path.join(CV_OUT_DIR, "cv_counts.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")

## plot predicted versus observed for each fold, positive and negative strands
# top row = all elements, bottom row = nonzero elements

nonzero_cv_pred = cv_pred[cv_pred["true_logcounts_total"] > 0]

fig, axes = plt.subplots(1, 5, figsize=(20, 4))
for i in range(len(FOLDS)):
    this_fold = FOLDS[i]
    
    plot_scatter(cv_pred[cv_pred["fold"] == this_fold],
            "pred_logcounts_total", "true_logcounts_total",
            ax = axes[i],
            xlabel=f"Predicted log counts",
            ylabel="Observed log counts",
            title=f"EP300 counts ({this_fold})",
            show_legend = False)

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(CV_OUT_DIR, "pred_vs_observed_by_fold.pdf"))

## plot summary across folds
fig, ax = plt.subplots(1, 1, figsize=(5,5))

# all elements
plot_scatter(cv_pred,
            "pred_logcounts_total", "true_logcounts_total",
            ax = ax,
            xlabel=f"Predicted log counts",
            ylabel="Observed log counts",
            title=f"EP300 counts (all folds)")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(CV_OUT_DIR, "pred_vs_observed.all_folds.pdf"))

