#!/usr/bin/env python3
"""
Compute and plot what fraction of the total SHAP signal is explained by FiNeMo motif hits.

Denominator: L1 norm of SHAP contributions per peak (from contribution_regions.npz).
Numerator: sum(hit_importance) per peak (from hits.tsv).

EP300_peak_overlap loaded from chromatin annotations joined by peak_id.

Outputs (under <finemo-dir>/annotated_motifs/):
  motif_explained_fraction.pdf          - bar chart per motif + unexplained
  motif_explained_fraction.motifs_only.pdf - motifs only (zoomed)

Usage:
  python3 scripts/plot_finemo_explained_fraction.py \
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

CHROMATIN_ANNOT = (
    "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/"
    "2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"
)
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
    out_dir = os.path.join(finemo_dir, "annotated_motifs")
    os.makedirs(out_dir, exist_ok=True)

    # ── motif order and colors ─────────────────────────────────────────────────
    report = pd.read_table(os.path.join(finemo_dir, "motif_report.tsv"))
    report = report[~report["motif_name"].isin(EXCLUDE_MOTIFS)]
    report = report.sort_values("num_hits_total", ascending=False)
    motif_order = report["motif_name"].tolist()
    motif_color = dict(zip(motif_order, NATURE_COLORS))
    motif_color["Unexplained"] = "#96a0b3"

    # ── load contributions ─────────────────────────────────────────────────────
    print("Loading contribution_regions.npz ...", flush=True)
    npz = np.load(os.path.join(finemo_dir, "contribution_regions.npz"))
    contributions = npz["contributions"].astype(np.float32)  # (n_peaks, 4, 500)
    peak_ids_npz = npz["peak_id"].astype(int)

    total_shap = np.abs(contributions).sum(axis=(1, 2))
    shap_df = pd.DataFrame({"peak_id": peak_ids_npz, "total_shap": total_shap.astype(float)})

    # ── load hits ─────────────────────────────────────────────────────────────
    print("Loading hits ...", flush=True)
    hits = pd.read_csv(os.path.join(finemo_dir, "hits.tsv"), sep="\t")
    hits = hits[~hits["motif_name"].isin(EXCLUDE_MOTIFS)]

    peak_hit_sum = (
        hits.groupby("peak_id")["hit_importance"]
        .sum()
        .reset_index()
        .rename(columns={"hit_importance": "sum_hit_importance"})
    )
    motif_hit_sum = (
        hits.groupby("motif_name")["hit_importance"]
        .sum()
        .reset_index()
        .rename(columns={"hit_importance": "sum_hit_importance"})
    )

    # ── load EP300_peak_overlap from chromatin annotations ────────────────────
    print("Loading chromatin annotations ...", flush=True)
    annot = pd.read_csv(CHROMATIN_ANNOT, sep="\t", usecols=["peak_id", "EP300_peak_overlap"])

    # ── merge ──────────────────────────────────────────────────────────────────
    df = (
        shap_df
        .merge(peak_hit_sum, on="peak_id", how="left")
        .merge(annot, on="peak_id", how="left")
    )
    df["sum_hit_importance"] = df["sum_hit_importance"].fillna(0.0)

    # ── overall stats ──────────────────────────────────────────────────────────
    total_signal = df["total_shap"].sum()
    total_explained = df["sum_hit_importance"].sum()
    explained_frac = total_explained / total_signal

    print(f"\nOverall explained fraction: {explained_frac:.3f}  ({explained_frac*100:.1f}%)")
    print(f"  Total |SHAP| signal:     {total_signal:.1f}")
    print(f"  Motif-explained signal:  {total_explained:.1f}")
    print(f"  Unexplained:             {total_signal - total_explained:.1f}  ({(1-explained_frac)*100:.1f}%)")

    for label, mask in [
        ("p300+ peaks", df["EP300_peak_overlap"] == 1),
        ("p300- peaks", df["EP300_peak_overlap"] == 0),
    ]:
        sub = df[mask]
        frac = sub["sum_hit_importance"].sum() / sub["total_shap"].sum()
        print(f"  {label}: {frac:.3f}  ({frac*100:.1f}%)")

    # ── per-motif fractions ────────────────────────────────────────────────────
    motif_hit_sum["pct_of_total"] = motif_hit_sum["sum_hit_importance"] / total_signal * 100
    motif_hit_sum = motif_hit_sum.sort_values("pct_of_total", ascending=False)

    unexplained_pct = (1 - explained_frac) * 100
    plot_df = pd.concat([
        motif_hit_sum,
        pd.DataFrame({"motif_name": ["Unexplained"],
                      "sum_hit_importance": [total_signal - total_explained],
                      "pct_of_total": [unexplained_pct]}),
    ], ignore_index=True)

    print("\nPer-motif fraction of total SHAP:")
    for _, row in plot_df.iterrows():
        print(f"  {row['motif_name']:20s}  {row['pct_of_total']:.1f}%")

    # ── plot ──────────────────────────────────────────────────────────────────
    def make_bar_chart(df, out_path, figsize=(4.5, 3.5)):
        fig, ax = plt.subplots(figsize=figsize)
        colors = [motif_color.get(m, "#888888") for m in df["motif_name"]]
        y = np.arange(len(df))
        ax.barh(y, df["pct_of_total"], color=colors, height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(df["motif_name"], fontsize=8)
        ax.set_xlabel("% of total SHAP signal", fontsize=9)
        ax.invert_yaxis()
        style_ax(ax)
        plt.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close()
        print(f"Saved {out_path}")

    make_bar_chart(plot_df, os.path.join(out_dir, "motif_explained_fraction.pdf"))
    make_bar_chart(motif_hit_sum, os.path.join(out_dir, "motif_explained_fraction.motifs_only.pdf"),
                   figsize=(4.5, 2.8))


if __name__ == "__main__":
    main()
