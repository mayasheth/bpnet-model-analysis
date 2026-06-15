#!/usr/bin/env python3
"""
Publication-ready heatmaps from motif pair insertion analysis.
  1. Max log2 FC vs baseline (sequential colormap)
  2. Max log2 synergy — observed vs additive (diverging colormap)

Motif order matches individual motif insertions plot: sections
(p300-important, DNase-important, neither), within each by decreasing
mean individual effect.

Usage:
  python3 scripts/plot_motif_pair_heatmaps.py \
      --raw      <motif_pairs_v1/motif_pairs.raw_results.tsv.gz> \
      --singles  <motif_pairs_v1/individual_motifs.summary.tsv> \
      --out-dir  <motif_pairs_v1/plots>
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
import matplotlib.colors as mcolors

# ── collapse motifs before any grouping ───────────────────────────────────────
COLLAPSE = {"CREB_ATF_3": "AP1_2"}

# ── display names (no sequences) ──────────────────────────────────────────────
DISPLAY_NAMES = {
    "ETS":         "ETS",
    "GATA":        "GATA",
    "AP1_2":       "AP-1",
    "CREB_ATF_1":  "ATF::C/EBP",
    "STAT":        "STAT",
    "Grepeats":    "Poly-G",
    "CTCF":        "CTCF",
    "NRF1":        "NRF1",
    "NFY":         "NF-Y",
    "NFI":         "NF-I",
    "RFX_1":       "RFX",
    "REST":        "REST",
    "Ebox_CACGTG": "E-box",
    "E2F":         "E2F",
    "TATAbox":     "TATA-box",
}

SECTIONS = {
    "p300-important":  ["ETS", "AP1_2", "CREB_ATF_1", "STAT", "Grepeats", "GATA"],
    "DNase-important": ["CTCF", "NRF1", "NFY", "NFI", "RFX_1", "REST", "Ebox_CACGTG"],
    "neither":         ["E2F", "TATAbox"],
}


def get_motif_order(singles_summary):
    mean_fc = singles_summary.set_index("motif_name")["mean_log2fc"].to_dict()
    order = []
    for members in SECTIONS.values():
        present = [m for m in members if m in mean_fc]
        present.sort(key=lambda m: mean_fc[m], reverse=True)
        order.extend(present)
    return order


def build_symmetric_matrix(df_agg, motif_order, value_col):
    """Average over model folds, take max over orientations/spacings, build symmetric matrix."""
    # Average across folds per (motif1, motif2, orientation, spacing)
    group_mean = (
        df_agg
        .groupby(["motif1_name", "motif2_name", "orientation_pattern", "spacing"])[value_col]
        .mean()
        .reset_index()
    )
    # Max over orientations and spacings per pair
    pair_max = (
        group_mean
        .groupby(["motif1_name", "motif2_name"])[value_col]
        .max()
        .reset_index()
    )

    n = len(motif_order)
    mat = pd.DataFrame(np.nan, index=motif_order, columns=motif_order)
    for _, row in pair_max.iterrows():
        m1, m2, val = row["motif1_name"], row["motif2_name"], row[value_col]
        if m1 in mat.index and m2 in mat.columns:
            mat.loc[m1, m2] = val
            mat.loc[m2, m1] = val
    return mat


def get_cmaps():
    import cmcrameri.cm as cmc
    return cmc.vik_r, cmc.roma  # log2fc (blue=positive), synergy


def plot_heatmap(matrix, out_path, cbar_label, cmap, vmin=None, vmax=None, center=None):
    labels = [DISPLAY_NAMES.get(m, m) for m in matrix.index]
    n = len(labels)
    fig_size = max(4.5, n * 0.32 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.88))

    data = matrix.values.astype(float)

    if center is not None and vmin is None:
        abs_max = np.nanmax(np.abs(data))
        vmin, vmax = -abs_max, abs_max

    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7, colors="black")
    cbar.outline.set_edgecolor("black")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(colors="black", length=2)

    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    MOTIF_DIR = (
        "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
        "2025_0517_official_EP300_K562_model/motif_spacing/motif_pairs_v1"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw",     default=f"{MOTIF_DIR}/motif_pairs.raw_results.tsv.gz")
    parser.add_argument("--singles", default=f"{MOTIF_DIR}/individual_motifs.summary.tsv")
    parser.add_argument("--out-dir", default=f"{MOTIF_DIR}/plots")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading data ...")
    raw     = pd.read_csv(args.raw,     sep="\t")
    singles = pd.read_csv(args.singles, sep="\t")

    # ── collapse motifs ────────────────────────────────────────────────────────
    raw["motif1_name"] = raw["motif1_name"].replace(COLLAPSE)
    raw["motif2_name"] = raw["motif2_name"].replace(COLLAPSE)
    singles["motif_name"] = singles["motif_name"].replace(COLLAPSE)
    singles = singles.groupby("motif_name", as_index=False)["mean_log2fc"].mean()

    motif_order = get_motif_order(singles)
    print(f"Motif order ({len(motif_order)}): {motif_order}")

    cmap_fc, _ = get_cmaps()
    cmap_syn = cmap_fc

    # ── matrix 1: max log2 FC vs baseline ─────────────────────────────────────
    print("Building max log2 FC vs baseline matrix ...")
    mat_fc = build_symmetric_matrix(raw, motif_order, "log2_fc_vs_baseline")
    plot_heatmap(
        mat_fc,
        out_path=os.path.join(args.out_dir, "pairs_max_log2fc_vs_baseline_v2.pdf"),
        cbar_label="Max log₂ FC vs baseline",
        cmap=cmap_fc,
        center=0,
    )

    # ── matrix 2: max log2 synergy (observed vs additive) ─────────────────────
    print("Building max synergy matrix ...")
    mat_syn = build_symmetric_matrix(raw, motif_order, "log2_synergy")
    plot_heatmap(
        mat_syn,
        out_path=os.path.join(args.out_dir, "pairs_max_synergy_v2.pdf"),
        cbar_label="Max log₂ synergy (obs − additive)",
        cmap=cmap_syn,
        center=0,
    )


if __name__ == "__main__":
    main()
