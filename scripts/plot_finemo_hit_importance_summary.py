#!/usr/bin/env python3
"""
Two summary figures for FiNeMo hit importance:
  1. motif_importance_violin.pdf       - distribution of hit_importance per motif, ordered by median
  2. motif_importance_vs_frequency.pdf - % total hits vs % total importance per motif

Usage:
  python3 scripts/plot_finemo_hit_importance_summary.py \
      --finemo-dir 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

EXCLUDE_MOTIFS = {"FOS_JUN", "GATA_TAL1_8BP"}
NATURE_COLORS = [
    "#dc6464", "#5496ce", "#e9c54e", "#c5c500",
    "#5eb342", "#49bcbc", "#b778b3", "#f29742", "#bc9678",
]


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--finemo-dir",
        default="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
                "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2",
    )
    args = parser.parse_args()
    finemo_dir = args.finemo_dir
    out_dir = finemo_dir

    # ── load data ──────────────────────────────────────────────────────────────
    print("Loading hits ...")
    hits = pd.read_csv(os.path.join(finemo_dir, "hits.tsv"), sep="\t")
    hits = hits[~hits["motif_name"].isin(EXCLUDE_MOTIFS)]

    report = pd.read_table(os.path.join(finemo_dir, "motif_report.tsv"))
    report = report[~report["motif_name"].isin(EXCLUDE_MOTIFS)]
    report = report.sort_values("num_hits_total", ascending=False)
    motif_color = dict(zip(report["motif_name"], NATURE_COLORS))

    # ── shared data for both figures ──────────────────────────────────────────
    print("Plotting violin ...")
    # ascending so highest-median motif is at top of the horizontal plot
    medians = hits.groupby("motif_name")["hit_importance"].median().sort_values(ascending=True)
    motif_order = medians.index.tolist()
    y = np.arange(len(motif_order))
    groups = [hits.loc[hits["motif_name"] == m, "hit_importance"].values for m in motif_order]

    total_hits = hits["motif_name"].value_counts()
    total_importance = hits.groupby("motif_name")["hit_importance"].sum()
    pct_hits = (total_hits / total_hits.sum() * 100).rename("pct_hits")
    pct_imp  = (total_importance / total_importance.sum() * 100).rename("pct_importance")
    summary = pd.concat([pct_hits, pct_imp], axis=1).loc[motif_order]
    n_per_motif = hits.groupby("motif_name").size()

    def _draw_violin(ax):
        parts = ax.violinplot(
            groups, positions=y,
            showmedians=False, showextrema=False, widths=0.75, vert=False,
        )
        for body, m in zip(parts["bodies"], motif_order):
            body.set_facecolor(motif_color.get(m, "#888888"))
            body.set_alpha(0.7)
            body.set_edgecolor("none")
        ax.boxplot(
            groups, positions=y, widths=0.15, patch_artist=False,
            showfliers=False, vert=False,
            medianprops=dict(color="black", linewidth=1.5),
            boxprops=dict(color="black", linewidth=0.8),
            whiskerprops=dict(color="black", linewidth=0.8),
            capprops=dict(color="black", linewidth=0.8),
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{m}  (n={n_per_motif[m]:,})" for m in motif_order], fontsize=8,
        )
        ax.set_xlabel("Hit importance", fontsize=9)
        style_ax(ax)

    # ── figure 1: standalone violin ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5, 3.5))
    _draw_violin(ax)
    plt.tight_layout()
    out = os.path.join(out_dir, "motif_importance_violin.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")

    # ── figure 2: combined violin + bar chart ─────────────────────────────────
    print("Plotting frequency vs importance ...")
    fig, (ax_vio, ax_bar) = plt.subplots(
        1, 2, figsize=(8, 3.5), sharey=True,
        gridspec_kw={"width_ratios": [1, 1]},
    )
    _draw_violin(ax_vio)

    bar_h = 0.35
    ax_bar.barh(y + bar_h / 2, summary["pct_hits"], height=bar_h, color="#c5cad7")
    ax_bar.barh(y - bar_h / 2, summary["pct_importance"], height=bar_h,
                color=[motif_color.get(m, "#888") for m in motif_order])
    ax_bar.set_xlabel("% of total", fontsize=9)
    ax_bar.axvline(0, color="black", linewidth=0.8)
    ax_bar.legend(handles=[
        mpatches.Patch(color="#c5cad7", label="% of total hits"),
        mpatches.Patch(color="#888888", label="% of total importance"),
    ], fontsize=7, frameon=False, loc="lower right")
    style_ax(ax_bar)
    plt.tight_layout()
    out = os.path.join(out_dir, "motif_importance_vs_frequency.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")

    print("Done.")


if __name__ == "__main__":
    main()
