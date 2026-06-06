#!/usr/bin/env python
"""
Plot model comparison bar chart: Pearson r for p300 prediction across models.

Reads per-model prediction_accuracy.tsv files and replicate_correlations.tsv,
then produces a grouped bar chart with groups for "all elements" and "p300+ elements".

Usage:
  conda activate analysis
  python scripts/plot_model_comparison.py --output figures/model_comparison.pdf

Edit MODEL_CONFIGS below to add/remove models or update file paths.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ── Model definitions ─────────────────────────────────────────────────────────
# Each entry: (display label, prediction_accuracy.tsv path, accuracy_subset strings, color)
# accuracy_subset strings must match the "subset" column in prediction_accuracy.tsv.
# Set path to None to skip that model.

PROJECT = os.environ.get("OAK", "/oak/stanford/groups/engreitz") + \
          "/Users/sheth/EP300_BPNet"

MODEL_CONFIGS = [
    {
        "label": "p300 BPNet\n(sequence only)",
        "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/predictions_cv/all_folds/prediction_accuracy.tsv",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "color": "#792374",  # purple
    },
    {
        "label": "p300 multimodal\n(+ ATAC)",
        "acc_tsv": f"{PROJECT}/2026_0529_multimodal_p300_model/predictions/atac/prediction_accuracy.tsv",
        "subset_all":   "CV — all folds",
        "subset_peaks": "CV — p300+",
        "color": "#006eae",  # blue
    },
    {
        "label": "ATAC\nChromBPNet",
        "acc_tsv": None,
        "subset_all":   None,
        "subset_peaks": None,
        "hardcoded_all":   0.70,   # from ChromBPNet manuscript (ATAC Pearson r)
        "hardcoded_peaks": None,
        "color": "#49bcbc",  # teal
    },
    {
        "label": "GATA1 BPNet",
        "acc_tsv": f"{PROJECT}/K562_GATA1_BPNet/predictions_cv/all_folds/prediction_accuracy.tsv",
        "subset_all":   "CV — all elements",
        "subset_peaks": "CV — p300+",
        "color": "#5eb342",  # green
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
    hardcoded_key = "hardcoded_" + subset_key.replace("subset_", "")
    if hardcoded_key in cfg:
        return cfg[hardcoded_key] if cfg[hardcoded_key] is not None else np.nan
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


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--output", default=f"{PROJECT}/figures/model_comparison.pdf")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    labels     = [c["label"] for c in MODEL_CONFIGS]
    colors     = [c["color"] for c in MODEL_CONFIGS]
    all_vals   = [load_pearson(c, "subset_all")   for c in MODEL_CONFIGS]
    peak_vals  = [load_pearson(c, "subset_peaks")  for c in MODEL_CONFIGS]

    print("\nAll elements:")
    for l, v in zip(labels, all_vals):
        print(f"  {l.replace(chr(10), ' '):<30} {v:.4f}" if np.isfinite(v) else f"  {l:<30} N/A")
    print("\np300+ elements:")
    for l, v in zip(labels, peak_vals):
        print(f"  {l.replace(chr(10), ' '):<30} {v:.4f}" if np.isfinite(v) else f"  {l:<30} N/A")

    n = len(MODEL_CONFIGS)
    x = np.arange(n)
    width = 0.35
    gap   = 0.08

    fig, ax = plt.subplots(figsize=(max(7, n * 1.4 + 2), 5))

    bars_all  = ax.bar(x - width/2 - gap/2, all_vals,  width, color=colors,
                       edgecolor="black", linewidth=0.8, label="All elements")
    bars_peak = ax.bar(x + width/2 + gap/2, peak_vals, width, color=colors,
                       edgecolor="black", linewidth=0.8, alpha=0.55, label="p300+ elements")

    # Value labels on bars
    for bar, val in list(zip(bars_all, all_vals)) + list(zip(bars_peak, peak_vals)):
        if np.isfinite(val):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=8, color="black")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10, color="black")
    ax.set_ylabel("Pearson r (log counts)", color="black")
    ax.set_ylim(0, 1.05)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    solid_patch  = mpatches.Patch(facecolor="grey", edgecolor="black", linewidth=0.8,
                                  label="All elements")
    alpha_patch  = mpatches.Patch(facecolor="grey", edgecolor="black", linewidth=0.8,
                                  alpha=0.55, label="p300+ elements")
    ax.legend(handles=[solid_patch, alpha_patch], frameon=False, fontsize=9,
              loc="upper right")

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches="tight")
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
