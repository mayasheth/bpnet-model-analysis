#!/usr/bin/env python3
"""
Plot forward CWM logos for all FiNeMo motifs in a single horizontal row,
ordered by decreasing total hit count.

Uses motif_cwms.npy + motif_data.tsv for the pre-computed restricted-region
window (same window used for individual logo SVGs).

Usage:
    pixi run -e ism python scripts/plot_finemo_motif_logos_row.py \
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
import logomaker

EXCLUDE_MOTIFS = {"FOS_JUN", "GATA_TAL1_8BP"}
BASES = ["A", "C", "G", "T"]


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--finemo-dir",
        default="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
                "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2",
    )
    args = parser.parse_args()
    finemo_dir = args.finemo_dir

    cwms = np.load(os.path.join(finemo_dir, "motif_cwms.npy"))   # (N, 4, L)
    meta = pd.read_table(os.path.join(finemo_dir, "motif_data.tsv"))
    report = pd.read_table(os.path.join(finemo_dir, "motif_report.tsv"))

    report = report[~report["motif_name"].isin(EXCLUDE_MOTIFS)]
    report = report.sort_values("num_hits_total", ascending=False)
    motifs = report["motif_name"].tolist()

    fwd = meta[meta["strand"] == "+"].set_index("motif_name")

    n = len(motifs)
    fig, axes = plt.subplots(1, n, figsize=(n * 2.2, 1.5))
    if n == 1:
        axes = [axes]

    for ax, motif in zip(axes, motifs):
        row = fwd.loc[motif]
        cwm_slice = cwms[int(row["motif_id"])].T[int(row["motif_start"]):int(row["motif_end"])]
        df = pd.DataFrame(cwm_slice, columns=BASES)

        logomaker.Logo(df, ax=ax, color_scheme="classic", font_name="DejaVu Sans")
        ax.set_title(motif, fontsize=7, pad=3)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xticks([])
        ax.set_ylabel("")
        style_ax(ax)

    plt.tight_layout(pad=0.5)

    out = os.path.join(finemo_dir, "motif_logos_row.pdf")
    fig.savefig(out, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close()


if __name__ == "__main__":
    main()
