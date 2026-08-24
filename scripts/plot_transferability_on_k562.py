#!/usr/bin/env python
"""
Plot K562 cross-cell-type transferability results (mirror of plot_transferability.py).

Fig 2e — bar chart: Pearson r for GM12878-trained models evaluated on K562
         elements, plus K562-trained in-cell-type models and K562 inter-replicate
         ceiling.
Fig S2c — 3-panel scatter (all elements): k562_bpnet / gm12878_bpnet_on_k562 /
          gm12878_multimodal_on_k562
Fig S2d — 3-panel scatter (p300+ elements only)

Usage:
  pixi run -e ism python scripts/plot_transferability_on_k562.py \
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
K562_DIR = f"{PROJECT}/2025_0517_official_EP300_K562_model"
MM_DIR = f"{PROJECT}/2026_0529_multimodal_p300_model"

# ── Model definitions for the K562-elements bar chart ────────────────────────
MODELS = [
    {
        "label": "K562 BPNet\n(seq only, in-cell-type ceiling)",
        "acc_tsv": f"{K562_DIR}/predictions_mean/all_folds/prediction_accuracy.tsv",
        "tsv":     f"{K562_DIR}/predictions_mean/all_folds/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "pred_col": "mean_pred_logcounts",
        "color": "#a5a083",  # stone
    },
    {
        "label": "GM12878 BPNet\n(seq only, cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_bpnet_on_k562/mean/all_folds/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_bpnet_on_k562/mean/all_folds/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "pred_col": "mean_pred_logcounts",
        "color": "#792374",  # purple
    },
    {
        "label": "GM12878 ATAC-only\n(cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_atac_only_on_k562/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_atac_only_on_k562/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "pred_col": "mean_pred_logcounts",
        "color": "#0096a0",  # dark teal
    },
    {
        "label": "GM12878 multimodal\n(+ ATAC, cross-cell-type)",
        "acc_tsv": f"{GM_DIR}/predictions/gm12878_multimodal_on_k562/prediction_accuracy.tsv",
        "tsv":     f"{GM_DIR}/predictions/gm12878_multimodal_on_k562/mean_predictions.tsv.gz",
        "subset_all":   "mean — all elements",
        "subset_peaks": "mean — p300+",
        "pred_col": "mean_pred_logcounts",
        "color": "#006eae",  # blue
    },
    {
        "label": "K562 ATAC-only\n(in-cell-type)",
        "acc_tsv": f"{MM_DIR}/predictions/atac_only/prediction_accuracy.tsv",
        "tsv":     f"{MM_DIR}/predictions/atac_only/cv_predictions.tsv.gz",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "pred_col": "pred_logcounts",
        "color": "#49bcbc",  # teal
    },
    {
        "label": "K562 multimodal\n(+ ATAC, in-cell-type)",
        "acc_tsv": f"{MM_DIR}/predictions/atac/prediction_accuracy.tsv",
        "tsv":     f"{MM_DIR}/predictions/atac/cv_predictions.tsv.gz",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "pred_col": "pred_logcounts",
        "color": "#5496ce",  # medium blue
    },
    {
        "label": "K562\ninter-replicate ceiling",
        "acc_tsv": f"{K562_DIR}/replicate_correlations.tsv",
        "tsv":     None,
        "subset_all":   "all elements",
        "subset_peaks": "p300+ elements",
        "pred_col": None,
        "color": "#a5a083",  # stone
    },
]

TRUE_COL = "true_logcounts"
PEAK_COL = "EP300_peak_overlap"


def load_pearson(cfg, subset_key):
    path, subset = cfg["acc_tsv"], cfg[subset_key]
    if not os.path.exists(path):
        print(f"  WARNING: not found: {path}")
        return np.nan, np.nan
    df = pd.read_table(path)
    row = df[df["subset"] == subset]
    if len(row) == 0:
        print(f"  WARNING: subset '{subset}' not found in {path}")
        return np.nan, np.nan
    return float(row["pearson_r"].iloc[0]), int(row["n"].iloc[0])


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


def _add_value_labels_h(ax, bars, vals, ns):
    for bar, val, n in zip(bars, vals, ns):
        if np.isfinite(val):
            label = f"{val:.3f} (n={n:,})" if np.isfinite(n) else f"{val:.3f}"
            ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height() / 2,
                    label, ha="left", va="center", fontsize=8, color="black")


def _style_barh_ax(ax):
    ax.set_xlabel("Pearson r (log p300 counts)", color="black")
    ax.set_xlim(0, 1.35)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig_bar(out_dir):
    all_pairs  = [load_pearson(m, "subset_all")   for m in MODELS]
    peak_pairs = [load_pearson(m, "subset_peaks")  for m in MODELS]
    all_vals,  all_ns  = zip(*all_pairs)
    peak_vals, peak_ns = zip(*peak_pairs)
    labels    = [m["label"] for m in MODELS]
    colors    = [m["color"] for m in MODELS]

    # Sort by all-elements Pearson r ascending
    order = np.argsort([v if np.isfinite(v) else -1 for v in all_vals])
    all_vals  = [all_vals[i]  for i in order]
    all_ns    = [all_ns[i]    for i in order]
    peak_vals = [peak_vals[i] for i in order]
    peak_ns   = [peak_ns[i]   for i in order]
    labels    = [labels[i]    for i in order]
    colors    = [colors[i]    for i in order]
    n = len(labels)

    # Layout 1: grouped (all + p300+ per model), horizontal bars
    y = np.arange(n)
    height, gap = 0.35, 0.08

    fig, ax = plt.subplots(figsize=(7.5, max(4, n * 0.8 + 1)))
    bars_all  = ax.barh(y + height/2 + gap/2, all_vals,  height, color=colors,
                        edgecolor="black", linewidth=0.8)
    bars_peak = ax.barh(y - height/2 - gap/2, peak_vals, height, color=colors,
                        edgecolor="black", linewidth=0.8, alpha=0.55)
    _add_value_labels_h(ax, bars_all,  all_vals,  all_ns)
    _add_value_labels_h(ax, bars_peak, peak_vals, peak_ns)
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
    solid = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, label="All elements")
    alpha = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, alpha=0.55, label="p300+ elements")
    ax.legend(handles=[solid, alpha], frameon=False, fontsize=8, loc="lower right")
    _style_barh_ax(ax)
    plt.tight_layout()
    path = os.path.join(out_dir, "transferability_bar_on_k562.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

    # Layout 2: two separate figures, one per subset
    height = 0.6
    y = np.arange(n)
    for vals, ns, title, fname in [
        (all_vals,  all_ns,  "All elements",   "transferability_bar_on_k562_all_elements.pdf"),
        (peak_vals, peak_ns, "p300+ elements", "transferability_bar_on_k562_p300plus.pdf"),
    ]:
        fig, ax = plt.subplots(figsize=(7.5, max(3, n * 0.6 + 1)))
        for yi, val, nn, color in zip(y, vals, ns, colors):
            if np.isfinite(val):
                bar = ax.barh(yi, val, height, color=color, edgecolor="black", linewidth=0.8)
                _add_value_labels_h(ax, bar, [val], [nn])
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
        plot_scatter(df[TRUE_COL].values, df[model["pred_col"]].values,
                     ax, model["color"], title, max_counts)
    sns.despine(trim=True)
    plt.tight_layout()
    suffix = "peaks" if peaks_only else "all"
    path = os.path.join(out_dir, f"transferability_scatter_on_k562_{suffix}.pdf")
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
    fig_bar(args.output_dir)
    scatter_panels(args.output_dir, peaks_only=False, max_counts=args.max_counts)
    scatter_panels(args.output_dir, peaks_only=True,  max_counts=args.max_counts)


if __name__ == "__main__":
    main()
