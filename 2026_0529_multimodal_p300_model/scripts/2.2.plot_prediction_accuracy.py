#!/usr/bin/env python
"""
Plot prediction accuracy for MultiModalBPNet on candidate elements.

Reads outputs from 2.1.predict_multimodal.py and produces:
  mean_predictions.pdf      scatter plots for mean predictions
  cv_predictions.pdf        scatter plots for CV predictions

Usage:
  python scripts/2.2.plot_prediction_accuracy.py --predictions-dir predictions/atac
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--predictions-dir", required=True,
                   help="Directory containing mean_predictions.tsv.gz and cv_predictions.tsv.gz")
    return p.parse_args()


def plot_scatter(df, x_col, y_col, ax, color="#792374",
                 xlabel="Observed log counts", ylabel="Predicted log counts", title=""):
    x = df[x_col].values
    y = df[y_col].values
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]

    cmap = mcolors.LinearSegmentedColormap.from_list("density", ["#ffffff", color])
    sns.kdeplot(x=x, y=y, ax=ax, cmap=cmap, fill=True)

    lo = 0
    hi = max(x.max(), y.max(), 20)
    ax.plot([lo, hi], [lo, hi], color="black", linestyle=(0, (5, 5)), lw=1.2, label="y = x")
    m, b = np.polyfit(x, y, 1)
    ax.plot([lo, hi], [m * lo + b, m * hi + b], color=color, lw=1.5, label="fit")

    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    ax.text(0.04, 0.96, f"Pearson r = {r:.3f}\nSpearman ρ = {rho:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="white", alpha=0.7))

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", "box")
    ax.set_xlabel(xlabel, color="black")
    ax.set_ylabel(ylabel, color="black")
    ax.set_title(title, fontsize=11)
    ax.tick_params(colors="black")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
    ax.legend(loc="lower right", frameon=False, fontsize=9)


def accuracy_row(df, x_col, y_col, label):
    x = df[x_col].values
    y = df[y_col].values
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    mse = np.mean((y - x) ** 2)
    return {"subset": label, "n": finite.sum(), "pearson_r": r,
            "spearman_rho": rho, "mse": mse}


def main():
    args = parse_args()

    mean_path = os.path.join(args.predictions_dir, "mean_predictions.tsv.gz")
    cv_path = os.path.join(args.predictions_dir, "cv_predictions.tsv.gz")

    mean_df = pd.read_csv(mean_path, sep="\t")
    cv_df = pd.read_csv(cv_path, sep="\t")

    peaks_only = mean_df["EP300_peak_overlap"] == 1
    non_peaks = mean_df["EP300_peak_overlap"] == 0

    # ── Mean predictions plot ─────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    plot_scatter(mean_df, "true_logcounts", "mean_pred_logcounts", axes[0],
                 title=f"Mean predictions — all elements (n={len(mean_df):,})")
    plot_scatter(mean_df[peaks_only], "true_logcounts", "mean_pred_logcounts", axes[1],
                 title=f"Mean predictions — p300+ (n={peaks_only.sum():,})")
    plot_scatter(mean_df[non_peaks], "true_logcounts", "mean_pred_logcounts", axes[2],
                 title=f"Mean predictions — p300− (n={non_peaks.sum():,})")

    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(os.path.join(args.predictions_dir, "mean_predictions.pdf"))
    plt.close()

    # ── CV predictions plot ───────────────────────────────────────────────
    folds = sorted(cv_df["fold"].unique())
    fig, axes = plt.subplots(1, len(folds) + 1, figsize=(5 * (len(folds) + 1), 5))

    plot_scatter(cv_df, "true_logcounts", "pred_logcounts", axes[0],
                 color="#006eae",
                 title=f"CV predictions — all folds (n={len(cv_df):,})")
    for i, fold in enumerate(folds):
        sub = cv_df[cv_df["fold"] == fold]
        plot_scatter(sub, "true_logcounts", "pred_logcounts", axes[i + 1],
                     color="#006eae",
                     title=f"CV fold {fold} (n={len(sub):,})")

    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(os.path.join(args.predictions_dir, "cv_predictions.pdf"))
    plt.close()

    # ── Accuracy table ────────────────────────────────────────────────────
    cv_peaks    = cv_df[cv_df["EP300_peak_overlap"] == 1]
    cv_non_peaks = cv_df[cv_df["EP300_peak_overlap"] == 0]

    rows = [
        accuracy_row(mean_df, "true_logcounts", "mean_pred_logcounts",
                     "mean — all elements"),
        accuracy_row(mean_df[peaks_only], "true_logcounts", "mean_pred_logcounts",
                     "mean — p300+"),
        accuracy_row(mean_df[non_peaks], "true_logcounts", "mean_pred_logcounts",
                     "mean — p300−"),
        accuracy_row(cv_df, "true_logcounts", "pred_logcounts",
                     "CV — all folds"),
        accuracy_row(cv_peaks, "true_logcounts", "pred_logcounts",
                     "CV — p300+"),
        accuracy_row(cv_non_peaks, "true_logcounts", "pred_logcounts",
                     "CV — p300−"),
    ]
    for fold in folds:
        sub = cv_df[cv_df["fold"] == fold]
        rows.append(accuracy_row(sub, "true_logcounts", "pred_logcounts",
                                 f"CV — fold {fold}"))

    acc_df = pd.DataFrame(rows)
    acc_df["pearson_r"] = acc_df["pearson_r"].round(4)
    acc_df["spearman_rho"] = acc_df["spearman_rho"].round(4)
    acc_df["mse"] = acc_df["mse"].round(4)
    print("\nPrediction accuracy:")
    print(acc_df.to_string(index=False))

    acc_df.to_csv(os.path.join(args.predictions_dir, "prediction_accuracy.tsv"),
                  sep="\t", index=False)


if __name__ == "__main__":
    main()
