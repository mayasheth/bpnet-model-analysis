#!/usr/bin/env python3
"""
Bar chart: Spearman r between each motif's per-peak hit activity and EP300 ChIP-seq RPM.

Data source: tf_correlation/correlations.tsv.gz (EP300 rows only).

Usage:
  python3 scripts/plot_motif_ep300_correlation.py

Output: tf_correlation/motif_ep300_correlation.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt

FINEMO_DIR = "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs"
CORR_FILE = os.path.join(FINEMO_DIR, "tf_correlation", "correlations.tsv.gz")
OUT_DIR = os.path.join(FINEMO_DIR, "tf_correlation")

MOTIF_COLORS = {
    "GATA":        "#dc6464",
    "AP1_1":       "#5496ce",
    "REPEAT_G":    "#e9c54e",
    "GATA_TAL1":   "#c5c500",
    "ETS_1":       "#5eb342",
    "STAT_2":      "#49bcbc",
    "CREB_ATF_3":  "#b778b3",
    "CTCF":        "#f29742",
    "CREB_ATF_1":  "#bc9678",
}


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


df = pd.read_csv(CORR_FILE, sep="\t")
ep300 = (
    df[df["tf"] == "EP300"]
    .sort_values("spearman_r", ascending=True)  # ascending for horizontal barh
    .reset_index(drop=True)
)

fig, ax = plt.subplots(figsize=(3.8, 3.2))

colors = [MOTIF_COLORS.get(m, "#888888") for m in ep300["motif"]]
y = np.arange(len(ep300))

ax.barh(y, ep300["spearman_r"], color=colors, height=0.65)
ax.set_yticks(y)
ax.set_yticklabels(ep300["motif"], fontsize=8)
ax.set_xlabel("Spearman r with p300 ChIP RPM", fontsize=9)
ax.axvline(0, color="black", linewidth=0.8)

style_ax(ax)
plt.tight_layout()

out = os.path.join(OUT_DIR, "motif_ep300_correlation.pdf")
fig.savefig(out, bbox_inches="tight")
plt.close()
print(f"Saved {out}")
