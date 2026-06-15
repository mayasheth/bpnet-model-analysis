#!/usr/bin/env python
"""
Plot GM12878 cross-cell-type transferability results.

Fig 2d  — bar chart: Pearson r for 3 models on GM12878 elements, plus
          GM12878 inter-replicate ceiling as a 4th bar group.
Fig S2a — 3-panel scatter (all elements): gm12878_bpnet / k562_v1 / k562_multimodal
Fig S2b — 3-panel scatter (p300+ elements only)

Usage:
  pixi run -e ism python scripts/plot_transferability.py \
      --output-dir 2026_0606_GM12878_transferability/figures/
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

OAK = os.environ.get("OAK", "/oak/stanford/groups/engreitz")
PROJECT = f"{OAK}/Users/sheth/EP300_BPNet"
GM_DIR = f"{PROJECT}/2026_0606_GM12878_transferability"

# ── Model definitions for Fig 2d ─────────────────────────────────────────────
MODELS = [
    {
        "label": "GM12878 BPNet\n(in-cell-type ceiling)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_bpnet/mean/all_folds/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_bpnet/mean/all_folds/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#a5a083",  # stone
    },
    {
        "label": "K562 BPNet\n(seq only, cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/k562_bpnet_v1/mean/all_folds/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/k562_bpnet_v1/mean/all_folds/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#792374",  # purple
    },
    {
        "label": "K562 ATAC-only\n(cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/k562_atac_only/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/k562_atac_only/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#0096a0",  # dark teal
    },
    {
        "label": "K562 multimodal\n(+ ATAC, cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/k562_multimodal_atac/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/k562_multimodal_atac/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#006eae",  # blue
    },
    {
        "label": "GM12878 ATAC-only\n(in-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_atac_only/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_atac_only/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#49bcbc",  # teal
    },
    {
        "label": "GM12878 multimodal\n(+ ATAC, in-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_multimodal_atac/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_multimodal_atac/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "color": "#5496ce",  # medium blue
    },
    {
        "label": "GM12878\ninter-replicate ceiling",
        "acc_tsv": f"{GM_DIR}/GM12878_replicate_correlations.tsv",
        "tsv":     None,
        "subset_all":   "all elements",
        "subset_peaks": "p300+ elements",
        "color": "#a5a083",  # stone
    },
]

PRED_COL = "mean_pred_logcounts"
TRUE_COL = "true_logcounts"
PEAK_COL = "EP300_peak_overlap"


def load_pearson(cfg, subset_key):
    path, subset = cfg["acc_tsv"], cfg[subset_key]
    if not os.path.exists(path):
        print(f"  WARNING: not found: {path}")
        return np.nan
    df = pd.read_table(path)
    row = df[df["subset"] == subset]
    if len(row) == 0:
        print(f"  WARNING: subset '{subset}' not found in {path}")
        return np.nan
    return float(row["pearson_r"].iloc[0])


def plot_scatter(x, y, ax, color, title, max_counts=10):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("d", ["#ffffff", color])
    sns.kdeplot(x=x, y=y, ax=ax, cmap=cmap, fill=True)
    lo, hi = 0, max_counts
    ax.plot([lo, hi], [lo, hi], color="black", ls=(0, (5, 5)), lw=1.2)
    m, b = np.polyfit(x, y, 1)
    ax.plot([lo, hi], [m * lo + b, m * hi + b], color=color, lw=1.5)
    r, _   = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    ax.text(0.04, 0.96, f"Pearson r = {r:.3f}\nSpearman ρ = {rho:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="white", alpha=0.7))
    ax.set(xlim=(lo, hi), ylim=(lo, hi), aspect="equal",
           xlabel="Observed log counts", ylabel="Predicted log counts")
    ax.set_title(title, fontsize=10)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")


def _add_value_labels_h(ax, bars, vals):
    for bar, val in zip(bars, vals):
        if np.isfinite(val):
            ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", ha="left", va="center", fontsize=8, color="black")


def _style_barh_ax(ax):
    ax.set_xlabel("Pearson r (log p300 counts)", color="black")
    ax.set_xlim(0, 1.10)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig2d(out_dir):
    all_vals  = [load_pearson(m, "subset_all")   for m in MODELS]
    peak_vals = [load_pearson(m, "subset_peaks")  for m in MODELS]
    labels    = [m["label"] for m in MODELS]
    colors    = [m["color"] for m in MODELS]

    # Sort by all-elements Pearson r ascending
    order = np.argsort([v if np.isfinite(v) else -1 for v in all_vals])
    all_vals  = [all_vals[i]  for i in order]
    peak_vals = [peak_vals[i] for i in order]
    labels    = [labels[i]    for i in order]
    colors    = [colors[i]    for i in order]
    n = len(labels)

    # Layout 1: grouped (all + p300+ per model), horizontal bars
    y = np.arange(n)
    height, gap = 0.35, 0.08

    fig, ax = plt.subplots(figsize=(7, max(4, n * 0.8 + 1)))
    bars_all  = ax.barh(y + height/2 + gap/2, all_vals,  height, color=colors,
                        edgecolor="black", linewidth=0.8)
    bars_peak = ax.barh(y - height/2 - gap/2, peak_vals, height, color=colors,
                        edgecolor="black", linewidth=0.8, alpha=0.55)
    _add_value_labels_h(ax, bars_all,  all_vals)
    _add_value_labels_h(ax, bars_peak, peak_vals)
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
    solid = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, label="All elements")
    alpha = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, alpha=0.55, label="p300+ elements")
    ax.legend(handles=[solid, alpha], frameon=False, fontsize=8, loc="lower right")
    _style_barh_ax(ax)
    plt.tight_layout()
    path = os.path.join(out_dir, "transferability_bar.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

    # Layout 2: two separate figures, one per subset
    height = 0.6
    y = np.arange(n)
    for vals, title, fname in [
        (all_vals,  "All elements",   "transferability_bar_all_elements.pdf"),
        (peak_vals, "p300+ elements", "transferability_bar_p300plus.pdf"),
    ]:
        fig, ax = plt.subplots(figsize=(7, max(3, n * 0.6 + 1)))
        for yi, val, color in zip(y, vals, colors):
            if np.isfinite(val):
                bar = ax.barh(yi, val, height, color=color, edgecolor="black", linewidth=0.8)
                _add_value_labels_h(ax, bar, [val])
        ax.set_yticks(y)
        ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
        ax.set_title(title, fontsize=11, color="black", fontweight="bold")
        _style_barh_ax(ax)
        plt.tight_layout()
        path = os.path.join(out_dir, fname)
        plt.savefig(path, bbox_inches="tight")
        plt.close()
        print(f"Saved: {path}")


def scatter_panels(out_dir, peaks_only, max_counts=10):
    scatter_models = [m for m in MODELS if m["tsv"] is not None]
    fig, axes = plt.subplots(1, len(scatter_models), figsize=(5 * len(scatter_models), 5))
    for ax, model in zip(axes, scatter_models):
        path = model["tsv"]
        if not os.path.exists(path):
            print(f"  WARNING: not found: {path}")
            ax.set_visible(False)
            continue
        df = pd.read_csv(path, sep="\t")
        if peaks_only:
            df = df[df[PEAK_COL] == 1]
        subset_label = "p300+" if peaks_only else "all elements"
        title = f"{model['label'].replace(chr(10), ' ')}\n{subset_label} (n={len(df):,})"
        plot_scatter(df[TRUE_COL].values, df[PRED_COL].values,
                     ax, model["color"], title, max_counts)
    sns.despine(trim=True)
    plt.tight_layout()
    suffix = "peaks" if peaks_only else "all"
    path = os.path.join(out_dir, f"transferability_scatter_{suffix}.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", default=f"{GM_DIR}/figures",
                   help="Directory to write PDFs (default: %(default)s)")
    p.add_argument("--max-counts", type=float, default=10)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    fig2d(args.output_dir)
    scatter_panels(args.output_dir, peaks_only=False, max_counts=args.max_counts)
    scatter_panels(args.output_dir, peaks_only=True,  max_counts=args.max_counts)


if __name__ == "__main__":
    main()
