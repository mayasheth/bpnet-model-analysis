#!/usr/bin/env python
"""
Supplemental figure: compare p300 BPNet training region strategies (v1, v2, v3).

Produces:
  training_region_comparison.pdf        — bar chart grouped by model (all / p300+)
  training_region_comparison_by_subset.pdf — bar chart grouped by subset

Usage:
  pixi run -e ism python scripts/plot_training_region_comparison.py \
      --output-dir figures/
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


PROJECT = os.environ.get("OAK", "/oak/stanford/groups/engreitz") + \
          "/Users/sheth/EP300_BPNet"

MODEL_CONFIGS = [
    {
        "label": "p300 BPNet v1\n(GC-matched negatives)",
        "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/predictions_cv/all_folds/prediction_accuracy.tsv",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "color": "#792374",  # purple
    },
    {
        "label": "p300 BPNet v2\n(TF-annotation negatives)",
        "acc_tsv": f"{PROJECT}/2025_0703_retrain_p300_model/predictions_cv/all_folds/prediction_accuracy.tsv",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "color": "#006eae",  # blue
    },
    {
        "label": "p300 BPNet v3\n(v2 + extended training)",
        "acc_tsv": f"{PROJECT}/2025_1016_p300_model_v3/predictions_cv/all_folds/prediction_accuracy.tsv",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "color": "#e96a00",  # orange
    },
    {
        "label": "Inter-replicate\nceiling",
        "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/replicate_correlations.tsv",
        "subset_all":   "all elements",
        "subset_peaks": "p300+ elements",
        "color": "#a5a083",  # stone
        "pearson_col": "pearson_r",
    },
]


def load_pearson(cfg, subset_key):
    path = cfg["acc_tsv"]
    subset = cfg.get(subset_key)
    if path is None or subset is None:
        return np.nan
    if not os.path.exists(path):
        print(f"  WARNING: not found: {path}")
        return np.nan
    df = pd.read_table(path)
    col = cfg.get("pearson_col", "pearson_r")
    row = df[df["subset"] == subset]
    if len(row) == 0:
        print(f"  WARNING: subset '{subset}' not found in {path}")
        return np.nan
    return float(row[col].iloc[0])


def _add_value_labels(ax, bars, vals):
    for bar, val in zip(bars, vals):
        if np.isfinite(val):
            ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", ha="left", va="center", fontsize=8, color="black")


def _style_ax(ax):
    ax.set_xlabel("Pearson r (log p300 counts)", color="black")
    ax.set_xlim(0, 1.10)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_by_model(out_dir, labels, colors, all_vals, peak_vals):
    n = len(labels)
    y = np.arange(n)
    height, gap = 0.35, 0.08

    fig, ax = plt.subplots(figsize=(7, max(4, n * 0.8 + 1)))
    bars_all  = ax.barh(y + height/2 + gap/2, all_vals,  height, color=colors,
                        edgecolor="black", linewidth=0.8)
    bars_peak = ax.barh(y - height/2 - gap/2, peak_vals, height, color=colors,
                        edgecolor="black", linewidth=0.8, alpha=0.55)
    _add_value_labels(ax, bars_all,  all_vals)
    _add_value_labels(ax, bars_peak, peak_vals)
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
    solid = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, label="All elements")
    alpha = mpatches.Patch(facecolor="grey", edgecolor="black", lw=0.8, alpha=0.55, label="p300+ elements")
    ax.legend(handles=[solid, alpha], frameon=False, fontsize=9, loc="lower right")
    _style_ax(ax)
    plt.tight_layout()
    path = os.path.join(out_dir, "training_region_comparison.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def plot_by_subset(out_dir, labels, colors, all_vals, peak_vals):
    n = len(labels)
    height = 0.6
    y = np.arange(n)

    for vals, title, suffix in [
        (all_vals,  "All elements",   "training_region_comparison_all_elements.pdf"),
        (peak_vals, "p300+ elements", "training_region_comparison_p300plus.pdf"),
    ]:
        fig, ax = plt.subplots(figsize=(7, max(3, n * 0.6 + 1)))
        for yi, val, color in zip(y, vals, colors):
            if np.isfinite(val):
                bar = ax.barh(yi, val, height, color=color, edgecolor="black", linewidth=0.8)
                _add_value_labels(ax, bar, [val])
        ax.set_yticks(y)
        ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
        ax.set_title(title, fontsize=11, color="black", fontweight="bold")
        _style_ax(ax)
        plt.tight_layout()
        path = os.path.join(out_dir, suffix)
        plt.savefig(path, bbox_inches="tight")
        plt.close()
        print(f"Saved: {path}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", default=f"{PROJECT}/figures",
                   help="Directory to write PDFs (default: %(default)s)")
    p.add_argument("--no-sort", action="store_true",
                   help="Preserve MODEL_CONFIGS order instead of sorting by all-element Pearson r")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    labels    = [c["label"] for c in MODEL_CONFIGS]
    colors    = [c["color"] for c in MODEL_CONFIGS]
    all_vals  = [load_pearson(c, "subset_all")   for c in MODEL_CONFIGS]
    peak_vals = [load_pearson(c, "subset_peaks")  for c in MODEL_CONFIGS]

    print("\nAll elements:")
    for l, v in zip(labels, all_vals):
        print(f"  {l.replace(chr(10), ' '):<40} {'N/A' if not np.isfinite(v) else f'{v:.4f}'}")
    print("\np300+ elements:")
    for l, v in zip(labels, peak_vals):
        print(f"  {l.replace(chr(10), ' '):<40} {'N/A' if not np.isfinite(v) else f'{v:.4f}'}")

    if not args.no_sort:
        order = np.argsort([v if np.isfinite(v) else -1 for v in all_vals])
        labels    = [labels[i]    for i in order]
        colors    = [colors[i]    for i in order]
        all_vals  = [all_vals[i]  for i in order]
        peak_vals = [peak_vals[i] for i in order]

    plot_by_model(args.output_dir, labels, colors, all_vals, peak_vals)
    plot_by_subset(args.output_dir, labels, colors, all_vals, peak_vals)


if __name__ == "__main__":
    main()
