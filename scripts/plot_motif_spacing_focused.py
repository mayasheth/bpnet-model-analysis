#!/usr/bin/env python3
"""
Focused spacing curve for GATA x E-box (0–50 bp, 1 bp resolution).
Shows log2_fc_vs_baseline vs inter-motif spacing with two orientation lines
(++ and -+; E-box is palindromic so ++ == +- and -+ == --), plus horizontal
reference lines for GATA alone and E-box alone derived from the dataset.

Usage:
  python3 scripts/plot_motif_spacing_focused.py \
      --data    motif_spacing/GATA_Ebox_50bp/raw_results.tsv \
      --out     motif_spacing/GATA_Ebox_50bp/spacing_focused_v1.pdf \
      --palette managua
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc

# E-box (CACGTG) is a perfect palindrome, so flipping its strand is a no-op.
# In this dataset: motif order = [GATA, Ebox], so orientation chars = (GATA, Ebox).
# ++ == +- and -+ == -- throughout.
# We plot two distinct cases: ++ (GATA forward) and -+ (GATA reversed).
ORIENTATIONS = ["++", "-+"]
ORIENT_LABELS = {
    "++": "GATA+, E-box (±)",
    "-+": "GATA−, E-box (±)",
}

# Stone palette for individual motif reference lines
REF_COLORS = {
    "GATA": "#5e5948",   # stone-dark
    "Ebox": "#a5a083",   # stone-mid
}
REF_LABELS = {
    "GATA": "GATA alone",
    "Ebox": "E-box alone",
}


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors="black", length=3, width=0.8)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


def main():
    DATA_DIR = (
        "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
        "2025_0517_official_EP300_K562_model/motif_spacing/GATA_Ebox_50bp"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",    default=f"{DATA_DIR}/raw_results.tsv")
    parser.add_argument("--out",     default=f"{DATA_DIR}/spacing_focused_v1.pdf")
    parser.add_argument("--palette", default="managua", choices=["managua", "romaO"],
                        help="cmcrameri palette for orientation lines")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    raw = pd.read_csv(args.data, sep="\t")
    spacings = sorted(raw["spacing"].unique())

    # Individual motif reference values derived from the paired experiment:
    # log2_fc_vs_single_X = log2_fc_pair - log2_fc_X_alone
    # → X_alone = mean(log2_fc_vs_baseline - log2_fc_vs_single_X)
    gata_alone = (raw["log2_fc_vs_baseline"] - raw["log2_fc_vs_single_GATA"]).mean()
    ebox_alone = (raw["log2_fc_vs_baseline"] - raw["log2_fc_vs_single_Ebox"]).mean()

    cmap = getattr(cmc, args.palette)
    orient_colors = [cmap(0.25), cmap(0.70)]

    fig, ax = plt.subplots(figsize=(5.5, 2.85))

    for orient, color in zip(ORIENTATIONS, orient_colors):
        grp = raw[raw["orientation_pattern"] == orient]
        stats = (
            grp.groupby("spacing")["log2_fc_vs_baseline"]
            .agg(["mean", "sem"])
            .reindex(spacings)
        )
        means = stats["mean"].values
        sems  = stats["sem"].values

        for fold in sorted(grp["model_fold"].unique()):
            fold_vals = (
                grp[grp["model_fold"] == fold]
                .set_index("spacing")["log2_fc_vs_baseline"]
                .reindex(spacings).values
            )
            ax.plot(spacings, fold_vals, color=color, linewidth=0.6,
                    alpha=0.3, zorder=2)

        ax.plot(spacings, means, color=color, linewidth=1.8,
                label=ORIENT_LABELS[orient], zorder=3)

    # Reference lines for individual motif effects
    ax.axhline(gata_alone, color=REF_COLORS["GATA"], linewidth=1.0,
               linestyle="--", alpha=0.8, zorder=1, label=REF_LABELS["GATA"])
    ax.axhline(ebox_alone, color=REF_COLORS["Ebox"], linewidth=1.0,
               linestyle="--", alpha=0.8, zorder=1, label=REF_LABELS["Ebox"])

    # vertical line at 9 bp peak
    ax.axvline(9, color="black", linewidth=0.8, linestyle="--", alpha=0.5, zorder=1)

    # y = 0 baseline
    ax.axhline(0, color="black", linewidth=0.5, linestyle=":", alpha=0.3, zorder=0)

    # x-axis: show ticks every 10 bp, minor ticks every 5
    ax.set_xticks(range(0, 51, 10))
    ax.set_xticks(range(0, 51, 5), minor=True)
    ax.set_xlim(-1, 51)

    ax.tick_params(axis="both", labelsize=8)
    ax.set_xlabel("Inter-motif spacing (bp)", fontsize=9)
    ax.set_ylabel("Log$_2$ fold change vs baseline", fontsize=9)
    ax.set_title("GATA + E-box", fontsize=9, pad=5)

    ax.legend(fontsize=8, frameon=False, loc="upper right")
    style_ax(ax)

    plt.tight_layout(pad=1.2)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
