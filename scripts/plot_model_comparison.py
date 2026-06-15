#!/usr/bin/env python
"""
Plot model comparison bar charts: Pearson r for p300 prediction across models.

Three named figures (--figure):
  k562-default         Fig 1d: GATA1, p300 seq-only, ATAC ChromBPNet, inter-replicate ceiling
  k562-multimodal  Fig 2 K562: ATAC counts, ATAC-only, seq-only, multimodal, ceiling
  (transferability is a separate script: plot_transferability.py)

Two layouts:
  Default: groups by model (each model has two side-by-side bars: all / p300+)
  --group-by-subset: groups by subset (all models together for "all elements",
                     then all models together for "p300+")

Usage:
  python scripts/plot_model_comparison.py --figure k562-default \
      --output figures/model_comparison.pdf
  python scripts/plot_model_comparison.py --figure k562-multimodal \
      --output figures/k562_multimodal_comparison.pdf
  python scripts/plot_model_comparison.py --figure k562-multimodal \
      --group-by-subset --output figures/k562_multimodal_comparison_by_subset.pdf

Edit FIGURE_CONFIGS below to add/remove models or update file paths.
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


# ── Model definitions ─────────────────────────────────────────────────────────
# Each entry: (display label, prediction_accuracy.tsv path, accuracy_subset strings, color)
# accuracy_subset strings must match the "subset" column in prediction_accuracy.tsv.
# Set path to None to skip that model.

PROJECT = os.environ.get("OAK", "/oak/stanford/groups/engreitz") + \
          "/Users/sheth/EP300_BPNet"

FIGURE_CONFIGS = {
    # Fig 1d: main model comparison (K562, sequence-era models)
    "k562-default": [
        {
            "label": "GATA1 BPNet",
            "acc_tsv": f"{PROJECT}/K562_GATA1_BPNet/predictions_cv/all_folds/prediction_accuracy.tsv",
            "subset_all":   "CV — all elements",
            "subset_peaks": "CV — p300+",
            "color": "#5eb342",  # green
        },
        {
            "label": "p300 BPNet\n(sequence only)",
            "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/predictions_cv/all_folds/prediction_accuracy.tsv",
            "subset_all":   "CV — all elements",
            "subset_peaks": "CV — p300+",
            "color": "#792374",  # purple
        },
        {
            "label": "ATAC\nChromBPNet",
            "acc_tsv": None,
            "hardcoded_all":   0.70,   # from ChromBPNet manuscript
            "hardcoded_peaks": None,
            "color": "#49bcbc",  # teal
        },
        {
            "label": "Inter-replicate\nceiling",
            "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/replicate_correlations.tsv",
            "subset_all":   "all elements",
            "subset_peaks": "p300+ elements",
            "color": "#a5a083",  # stone
            "pearson_col": "pearson_r",
        },
    ],

    # Fig 2 K562: multimodal model progression
    "k562-multimodal": [
        {
            "label": "ATAC counts\ncorrelation",
            "acc_tsv": None,
            "hardcoded_all":   0.590,  # Pearson(ATAC logcounts, p300 logcounts), K562
            "hardcoded_peaks": 0.335,
            "color": "#49bcbc",  # teal
        },
        {
            "label": "ATAC-only\nBPNet",
            "acc_tsv": f"{PROJECT}/2026_0529_multimodal_p300_model/predictions/atac_only/prediction_accuracy.tsv",
            "subset_all":   "CV — all elements",
            "subset_peaks": "CV — p300+",
            "color": "#0096a0",  # dark teal
        },
        {
            "label": "p300 BPNet\n(sequence only)",
            "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/predictions_cv/all_folds/prediction_accuracy.tsv",
            "subset_all":   "CV — all elements",
            "subset_peaks": "CV — p300+",
            "color": "#792374",  # purple
        },
        {
            "label": "p300 BPNet\n(multimodal)",
            "acc_tsv": f"{PROJECT}/2026_0529_multimodal_p300_model/predictions/atac/prediction_accuracy.tsv",
            "subset_all":   "CV — all elements",
            "subset_peaks": "CV — p300+",
            "color": "#006eae",  # blue
        },
        {
            "label": "Inter-replicate\nceiling",
            "acc_tsv": f"{PROJECT}/2025_0517_official_EP300_K562_model/replicate_correlations.tsv",
            "subset_all":   "all elements",
            "subset_peaks": "p300+ elements",
            "color": "#a5a083",  # stone
            "pearson_col": "pearson_r",
        },
    ],
}


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
    p.add_argument("--figure", default="k562-default", choices=list(FIGURE_CONFIGS.keys()),
                   help="Which figure config to use (default: k562-default)")
    p.add_argument("--output", default=f"{PROJECT}/figures/model_comparison.pdf")
    p.add_argument("--group-by-subset", action="store_true",
                   help="Group bars by subset (all / p300+) rather than by model")
    p.add_argument("--no-sort", action="store_true",
                   help="Preserve config order instead of sorting by all-element Pearson r")
    return p.parse_args()


def add_value_labels(ax, bars, vals):
    for bar, val in zip(bars, vals):
        if np.isfinite(val):
            ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height() / 2,
                    f"{val:.3f}", ha="left", va="center", fontsize=8, color="black")


def style_ax(ax):
    ax.set_xlabel("Pearson r (log counts)", color="black")
    ax.set_xlim(0, 1.10)
    ax.tick_params(colors="black")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_by_model(ax, labels, colors, all_vals, peak_vals):
    """Default layout: one group per model, two bars per group (all / p300+)."""
    n = len(labels)
    y = np.arange(n)
    height, gap = 0.35, 0.08

    bars_all  = ax.barh(y + height/2 + gap/2, all_vals,  height, color=colors,
                        edgecolor="black", linewidth=0.8)
    bars_peak = ax.barh(y - height/2 - gap/2, peak_vals, height, color=colors,
                        edgecolor="black", linewidth=0.8, alpha=0.55)

    add_value_labels(ax, bars_all,  all_vals)
    add_value_labels(ax, bars_peak, peak_vals)

    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")

    solid_patch = mpatches.Patch(facecolor="grey", edgecolor="black", linewidth=0.8,
                                 label="All elements")
    alpha_patch = mpatches.Patch(facecolor="grey", edgecolor="black", linewidth=0.8,
                                 alpha=0.55, label="p300+ elements")
    ax.legend(handles=[solid_patch, alpha_patch], frameon=False, fontsize=9,
              loc="lower right")


def plot_by_subset(output, labels, colors, all_vals, peak_vals):
    """Two separate figures: one for all elements, one for p300+ elements."""
    base, ext = os.path.splitext(output)
    n = len(labels)
    height = 0.6
    y = np.arange(n)

    for vals, title, suffix in [
        (all_vals,  "All elements",   "_all_elements"),
        (peak_vals, "p300+ elements", "_p300plus"),
    ]:
        fig, ax = plt.subplots(figsize=(7, max(3, n * 0.6 + 1)))
        for yi, val, color in zip(y, vals, colors):
            if np.isfinite(val):
                bar = ax.barh(yi, val, height, color=color, edgecolor="black", linewidth=0.8)
                add_value_labels(ax, bar, [val])
        ax.set_yticks(y)
        ax.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=9, color="black")
        ax.set_title(title, fontsize=11, color="black", fontweight="bold")
        style_ax(ax)
        plt.tight_layout()
        out_path = base + suffix + ext
        plt.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out_path}")


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    configs   = FIGURE_CONFIGS[args.figure]
    labels    = [c["label"] for c in configs]
    colors    = [c["color"] for c in configs]
    all_vals  = [load_pearson(c, "subset_all")   for c in configs]
    peak_vals = [load_pearson(c, "subset_peaks")  for c in configs]

    if not args.no_sort:
        order = np.argsort([v if np.isfinite(v) else -1 for v in all_vals])
        labels    = [labels[i]    for i in order]
        colors    = [colors[i]    for i in order]
        all_vals  = [all_vals[i]  for i in order]
        peak_vals = [peak_vals[i] for i in order]

    print("\nAll elements:")
    for l, v in zip(labels, all_vals):
        print(f"  {l.replace(chr(10), ' '):<30} {v:.4f}" if np.isfinite(v) else f"  {l.replace(chr(10), ' '):<30} N/A")
    print("\np300+ elements:")
    for l, v in zip(labels, peak_vals):
        print(f"  {l.replace(chr(10), ' '):<30} {v:.4f}" if np.isfinite(v) else f"  {l.replace(chr(10), ' '):<30} N/A")

    n = len(labels)
    if args.group_by_subset:
        plot_by_subset(args.output, labels, colors, all_vals, peak_vals)
    else:
        fig, ax = plt.subplots(figsize=(7, max(4, n * 0.8 + 1)))
        plot_by_model(ax, labels, colors, all_vals, peak_vals)
        style_ax(ax)
        plt.tight_layout()
        plt.savefig(args.output, bbox_inches="tight")
        print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
