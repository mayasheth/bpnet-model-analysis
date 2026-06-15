#!/usr/bin/env python3
"""
Horizontal violin + boxplot of individual motif insertion effects (log2 FC vs baseline).
Motifs grouped into three sections, sorted by decreasing mean effect within each section.

Usage:
  python3 scripts/plot_individual_motif_insertions.py \
      --raw     <motif_pairs_v1/individual_motifs.raw_results.tsv.gz> \
      --summary <motif_pairs_v1/individual_motifs.summary.tsv> \
      --out     <motif_pairs_v1/plots/individual_motif_insertions_v2.pdf>
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

# ── collapse motifs before any grouping ───────────────────────────────────────
# Maps source motif name → canonical name it is merged into
COLLAPSE = {"CREB_ATF_3": "AP1_2"}

# ── motif display names ────────────────────────────────────────────────────────
DISPLAY_NAMES = {
    "ETS":         "ETS (CGGAAG)",
    "GATA":        "GATA (GATAA)",
    "AP1_2":       "AP-1",
    "CREB_ATF_1":  "ATF::C/EBP (TGACGTCA)",
    "STAT":        "STAT (TTCCNGGAA)",
    "Grepeats":    "Poly-G",
    "CTCF":        "CTCF (CCNNNAGGGGGCG)",
    "NRF1":        "NRF1 (GCGCATGCGC)",
    "NFY":         "NF-Y (CCAAT)",
    "NFI":         "NF-I (TGGCNNNNNGCCA)",
    "RFX_1":       "RFX (GTTGCCATGGCAAC)",
    "REST":        "REST (CAGCACCNNGGACAG)",
    "Ebox_CACGTG": "E-box (CACGTG)",
    "E2F":         "E2F (TTTCCCGCCAAA)",
    "TATAbox":     "TATA-box (TATAAAA)",
}

# ── section membership (after collapse) ───────────────────────────────────────
SECTIONS = {
    "p300-important": ["ETS", "AP1_2", "CREB_ATF_1", "STAT", "Grepeats", "GATA"],
    "DNase-important": ["CTCF", "NRF1", "NFY", "NFI", "RFX_1", "REST", "Ebox_CACGTG"],
    "neither":         ["E2F", "TATAbox"],
}

SECTION_COLORS = {
    "p300-important":  "#b778b3",
    "DNase-important": "#5496ce",
    "neither":         "#96a0b3",
}

SECTION_GAP = 0.0  # no extra spacing between sections


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
    MOTIF_DIR = (
        "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
        "2025_0517_official_EP300_K562_model/motif_spacing/motif_pairs_v1"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw",     default=f"{MOTIF_DIR}/individual_motifs.raw_results.tsv.gz")
    parser.add_argument("--summary", default=f"{MOTIF_DIR}/individual_motifs.summary.tsv")
    parser.add_argument("--out",     default=f"{MOTIF_DIR}/plots/individual_motif_insertions_v2.pdf")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    raw     = pd.read_csv(args.raw,     sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")

    # ── collapse motifs ────────────────────────────────────────────────────────
    raw["motif_name"] = raw["motif_name"].replace(COLLAPSE)
    summary["motif_name"] = summary["motif_name"].replace(COLLAPSE)
    summary = summary.groupby("motif_name", as_index=False)["mean_log2fc"].mean()

    mean_fc = summary.set_index("motif_name")["mean_log2fc"].to_dict()

    # ── build ordered motif list: sections in order, within each by desc mean ──
    ordered_motifs = []
    for section, members in SECTIONS.items():
        present = [m for m in members if m in mean_fc]
        present.sort(key=lambda m: mean_fc[m], reverse=True)
        ordered_motifs.append((section, present))

    # ── assign y positions (bottom to top = first to last in list) ─────────────
    # We'll plot top-of-list at highest y so the figure reads top-down.
    # Collect all motifs in display order (top to bottom), then reverse for y.
    display_order = []  # top to bottom
    for section, members in ordered_motifs:
        display_order.extend([(section, m) for m in members])

    # y=0 is bottom; highest y = top. We want display_order[0] at the top.
    n_motifs = len(display_order)
    y_pos = {}
    y = n_motifs - 1
    prev_section = None
    gap_accumulated = 0
    positions = []  # (section, motif, y_coord)
    for i, (section, motif) in enumerate(display_order):
        if prev_section is not None and section != prev_section:
            gap_accumulated += SECTION_GAP
        coord = y - gap_accumulated
        positions.append((section, motif, coord))
        y_pos[motif] = coord
        prev_section = section
        y -= 1

    # ── figure ─────────────────────────────────────────────────────────────────
    fig_height = n_motifs * 0.22 + 0.8
    fig, ax = plt.subplots(figsize=(5.5, fig_height))

    # draw violins and boxplots per motif
    for section, motif, coord in positions:
        color  = SECTION_COLORS[section]
        values = raw.loc[raw["motif_name"] == motif, "log2_fc_vs_baseline"].dropna().values
        if len(values) == 0:
            continue
        parts = ax.violinplot(
            [values], positions=[coord],
            showmedians=False, showextrema=False, widths=0.7, vert=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.7)
            body.set_edgecolor("none")
        ax.boxplot(
            [values], positions=[coord], widths=0.15, patch_artist=False,
            showfliers=False, vert=False,
            medianprops=dict(color="black", linewidth=1.5),
            boxprops=dict(color="black", linewidth=0.8),
            whiskerprops=dict(color="black", linewidth=0.8),
            capprops=dict(color="black", linewidth=0.8),
        )

    # y-tick labels
    ytick_coords  = [coord for _, _, coord in positions]
    ytick_labels  = [DISPLAY_NAMES.get(motif, motif) for _, motif, _ in positions]
    ax.set_yticks(ytick_coords)
    ax.set_yticklabels(ytick_labels, fontsize=8)

    # reference line at 0
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)

    y_coords_all = [coord for _, _, coord in positions]
    ax.set_ylim(min(y_coords_all) - 0.6, max(y_coords_all) + 0.6)
    ax.set_xlim(-0.5, 1.0)
    ax.set_xlabel("Log$_2$ fold change vs baseline", fontsize=9)
    style_ax(ax)

    plt.tight_layout()
    fig.savefig(args.out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
