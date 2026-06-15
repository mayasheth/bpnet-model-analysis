#!/usr/bin/env python3
"""
Composite figure: 3-motif panel illustrating how ChIP correlation identifies
binding TFs (single factor vs large family vs co-binding complex).

Columns: CTCF (clean single-TF ID), AP1_1 (large family), GATA (co-binding complex)
Each column: CWM logo + horizontal bar chart of top-N correlated TFs,
             with a marker for TFs with known direct p300 interactions.

Usage:
  python3 scripts/plot_finemo_composite_figure.py

Output:
  tf_correlation/composite_motif_tf_figure.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import logomaker

FINEMO_DIR = "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs"
CORR_FILE  = os.path.join(FINEMO_DIR, "tf_correlation", "correlations.tsv.gz")
OUT_DIR    = os.path.join(FINEMO_DIR, "tf_correlation")

# Motifs to show and display titles
MOTIFS = ["CTCF", "AP1_1", "GATA"]
MOTIF_TITLES = {
    "CTCF":  "CTCF",
    "AP1_1": "AP-1",
    "GATA":  "GATA",
}
MOTIF_COLORS = {
    "CTCF":  "#f29742",
    "AP1_1": "#5496ce",
    "GATA":  "#dc6464",
}

TOP_N = 8   # TFs per panel

# TFs with documented direct interactions with p300/CBP.
# Update this set as you review the literature.
# Sources: PMID literature, BioGRID, STRING co-complex evidence.
P300_INTERACTORS = {
    "TAL1",    # TAL1/SCL recruits p300 at erythroid enhancers
    "GATA1",   # GATA1-p300 interaction well-documented
    "GATA2",   # GATA2-p300 interaction
    "FOSL1",   # AP-1 family recruits p300
    "FOSL2",
    "FOS",
    "FOSB",
    "JUN",
    "JUNB",
    "JUND",
    "NFE2",    # bZIP, AP-1 superfamily
    "MAFG",    # bZIP
    "MAFF",
    "SPI1",    # PU.1 interacts with p300
    "CEBPB",   # C/EBPβ interacts with p300
    "CEBPG",
    "CREM",
    "CREB1",   # CREB recruits p300 via KID domain
    "ATF1",
}
# Note: CBFA2T3, TCF3, ARID1B — unclear direct p300 interaction; leave unmarked.
# Update above set after literature/database review.


def load_cwm(motif_name):
    path = os.path.join(FINEMO_DIR, "CWMs", motif_name, "modisco_fc.txt")
    data = np.loadtxt(path)
    return pd.DataFrame(data.T, columns=["A", "C", "G", "T"])


def trim_cwm(df, threshold=0.01):
    keep = df.abs().max(axis=1)
    keep = keep[keep > threshold].index
    if len(keep) == 0:
        return df
    return df.loc[keep[0]:keep[-1]].reset_index(drop=True)


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


# ── load correlations ──────────────────────────────────────────────────────────

corr = pd.read_csv(CORR_FILE, sep="\t")

# ── figure layout ──────────────────────────────────────────────────────────────
# 3 columns × 2 rows (logo row + bar chart row)

fig = plt.figure(figsize=(10, 5.5))
outer = gridspec.GridSpec(
    2, 3,
    figure=fig,
    height_ratios=[1, 2.2],
    hspace=0.15,
    wspace=0.45,
)

for col_idx, motif in enumerate(MOTIFS):
    color = MOTIF_COLORS[motif]
    sub = corr[corr["motif"] == motif].sort_values("spearman_r", ascending=False)
    top = sub.head(TOP_N).reset_index(drop=True)

    # ── logo ──────────────────────────────────────────────────────────────────
    ax_logo = fig.add_subplot(outer[0, col_idx])
    try:
        cwm = trim_cwm(load_cwm(motif))
        cwm.index = range(len(cwm))
        logomaker.Logo(cwm, ax=ax_logo, color_scheme="classic")
    except Exception:
        ax_logo.text(0.5, 0.5, motif, ha="center", va="center", transform=ax_logo.transAxes)
    ax_logo.set_title(MOTIF_TITLES[motif], fontsize=10, fontweight="bold", pad=4)
    ax_logo.set_xticks([])
    ax_logo.set_yticks([])
    for sp in ax_logo.spines.values():
        sp.set_visible(False)

    # ── bar chart ─────────────────────────────────────────────────────────────
    ax_bar = fig.add_subplot(outer[1, col_idx])

    y = np.arange(len(top))[::-1]
    ax_bar.barh(y, top["spearman_r"], color=color, height=0.65, zorder=2)

    # p300 interaction marker: dot on the right side of bar + bold label
    labels = []
    for i, row in top.iterrows():
        tf = row["tf"]
        is_interactor = tf in P300_INTERACTORS
        label = tf
        labels.append(label)
        if is_interactor:
            ax_bar.plot(
                row["spearman_r"] + 0.008, y[i],
                marker="*", color="#792374", markersize=7,
                clip_on=False, zorder=3,
            )

    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(labels, fontsize=7.5)
    for tick, tf in zip(ax_bar.get_yticklabels(), top["tf"]):
        if tf in P300_INTERACTORS:
            tick.set_fontweight("bold")
    ax_bar.set_xlabel("Spearman r", fontsize=8)
    ax_bar.axvline(0, color="black", linewidth=0.8)
    ax_bar.set_xlim(left=0)
    style_ax(ax_bar)

# legend for p300 indicator
fig.text(
    0.98, 0.04,
    "★ direct p300 interaction",
    ha="right", va="bottom", fontsize=7, color="#792374",
    style="italic",
)

out = os.path.join(OUT_DIR, "composite_motif_tf_figure.pdf")
fig.savefig(out, bbox_inches="tight")
plt.close()
print(f"Saved {out}")
