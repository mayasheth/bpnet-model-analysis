#!/usr/bin/env python3
"""
Figure-ready 3-panel plot filtered to top 50% of hits by hit_coefficient_global:
  Panel 1 — histogram of motif hits per peak
  Panel 2 — predicted log counts distribution by n_hits (violin + boxplot)
  Panel 3 — observed log counts distribution by n_hits  (violin + boxplot)

Produces two PDFs: all elements and p300+ elements only.

Usage:
  python3 scripts/plot_finemo_counts_per_peak_top50pct.py \
      --hits   <annotated_motifs/hits_renamed.tsv> \
      --peaks  <annotated_motifs/hits_per_peak.with_predictions.tsv.gz> \
      --annot  <chromatin_annotations.tsv> \
      --out    <plot_annotated_motifs/counts_per_peak.top50pct.pdf>
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


VIOLIN_COLOR = "#b778b3"
MAX_N_HITS   = 9
POOL_ABOVE   = 7  # 7, 8, 9 → "7+"
CHROMATIN_ANNOT = (
    "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/"
    "2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"
)


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


def violin_boxplot(ax, groups, positions, color):
    """Vertical violin with overlaid thin boxplot (no fliers)."""
    non_empty = [(g, p) for g, p in zip(groups, positions) if len(g) > 0]
    if not non_empty:
        return
    g_vals, g_pos = zip(*non_empty)

    parts = ax.violinplot(
        g_vals, positions=g_pos,
        showmedians=False, showextrema=False, widths=0.7, vert=True,
    )
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_alpha(0.7)
        body.set_edgecolor("none")

    ax.boxplot(
        g_vals, positions=g_pos, widths=0.15, patch_artist=False,
        showfliers=False, vert=True,
        medianprops=dict(color="black", linewidth=1.5),
        boxprops=dict(color="black", linewidth=0.8),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )


def make_figure(peaks_df, out_path):
    # Pool 7, 8, 9+ into a single "7+" bin
    df = peaks_df.copy()
    df["n_hits_binned"] = df["n_hits"].clip(upper=POOL_ABOVE)
    n_vals  = list(range(POOL_ABOVE + 1))
    xlabels = [str(n) for n in range(POOL_ABOVE)] + [f"{POOL_ABOVE}+"]

    pred_groups = [
        df.loc[df["n_hits_binned"] == n, "mean_pred_logcounts"].dropna().values
        for n in n_vals
    ]
    obs_groups = [
        df.loc[df["n_hits_binned"] == n, "true_logcounts"].dropna().values
        for n in n_vals
    ]
    hist_counts = [int((df["n_hits_binned"] == n).sum()) for n in n_vals]

    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))

    # Panel 1: histogram with rotated count labels on bars
    ax = axes[0]
    ax.bar(n_vals, hist_counts, color=VIOLIN_COLOR, edgecolor="#430b4e", linewidth=0.5)
    y_max = max(hist_counts) * 1.15
    for n, c in zip(n_vals, hist_counts):
        if c > 0:
            ax.text(n, c + y_max * 0.01, f"{c:,}", ha="center", va="bottom",
                    fontsize=6, color="black", rotation=90)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Motif hits per peak", fontsize=9)
    ax.set_ylabel("Number of peaks", fontsize=9)
    ax.set_xticks(n_vals)
    ax.set_xticklabels(xlabels)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    style_ax(ax)

    # Panel 2: predicted logcounts
    ax = axes[1]
    violin_boxplot(ax, pred_groups, n_vals, VIOLIN_COLOR)
    ax.set_xlabel("Motif hits per peak", fontsize=9)
    ax.set_ylabel("Predicted log counts", fontsize=9)
    ax.set_xticks(n_vals)
    ax.set_xticklabels(xlabels)
    ax.set_ylim(0, 15)
    ax.set_yticks(range(0, 16, 3))
    style_ax(ax)

    # Panel 3: observed logcounts
    ax = axes[2]
    violin_boxplot(ax, obs_groups, n_vals, VIOLIN_COLOR)
    ax.set_xlabel("Motif hits per peak", fontsize=9)
    ax.set_ylabel("Observed log counts", fontsize=9)
    ax.set_xticks(n_vals)
    ax.set_xticklabels(xlabels)
    ax.set_ylim(0, 15)
    ax.set_yticks(range(0, 16, 3))
    style_ax(ax)

    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    FINEMO = (
        "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
        "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2"
    )
    parser.add_argument("--hits",  default=f"{FINEMO}/annotated_motifs/hits_renamed.tsv")
    parser.add_argument("--peaks", default=f"{FINEMO}/annotated_motifs/hits_per_peak.with_predictions.tsv.gz")
    parser.add_argument("--annot", default=CHROMATIN_ANNOT)
    parser.add_argument("--out",   default=f"{FINEMO}/plot_annotated_motifs/counts_per_peak.top50pct.pdf")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ── load ──────────────────────────────────────────────────────────────────
    print("Loading hits ...")
    hits  = pd.read_csv(args.hits,  sep="\t")
    peaks = pd.read_csv(args.peaks, sep="\t")

    print("Loading chromatin annotations ...")
    annot = pd.read_csv(args.annot, sep="\t", usecols=["peak_id", "EP300_peak_overlap"])
    peaks = peaks.merge(annot, on="peak_id", how="left")

    # ── filter top 50% by hit_coefficient_global ──────────────────────────────
    threshold = hits["hit_coefficient_global"].quantile(0.5)
    hits_top = hits[hits["hit_coefficient_global"] >= threshold].copy()
    print(f"Top 50% threshold: {threshold:.4e}  | hits retained: {len(hits_top):,} / {len(hits):,}")

    # ── recompute n_hits per peak ──────────────────────────────────────────────
    n_hits_per_peak = hits_top.groupby("peak_id").size().reset_index(name="n_hits")
    peaks_sub = (
        peaks.drop(columns=["n_hits"], errors="ignore")
             .merge(n_hits_per_peak, on="peak_id", how="left")
    )
    peaks_sub["n_hits"] = peaks_sub["n_hits"].fillna(0).astype(int)

    # ── all elements ───────────────────────────────────────────────────────────
    make_figure(peaks_sub, args.out)

    # ── p300+ elements only (EP300_peak_overlap flag) ─────────────────────────
    peaks_p300 = peaks_sub[peaks_sub["EP300_peak_overlap"] == 1]
    print(f"p300+ peaks (overlap flag): {len(peaks_p300):,} / {len(peaks_sub):,}")
    base, ext = os.path.splitext(args.out)
    make_figure(peaks_p300, f"{base}.p300plus{ext}")

    # ── top 20% by observed p300 counts + any remaining p300+ elements ────────
    obs_threshold = peaks_sub["true_logcounts"].quantile(0.8)
    in_top20 = peaks_sub["true_logcounts"] >= obs_threshold
    in_p300  = peaks_sub["EP300_peak_overlap"] == 1
    peaks_top_obs = peaks_sub[in_top20 | in_p300]
    n_top20_only = (in_top20 & ~in_p300).sum()
    n_p300_added = (~in_top20 & in_p300).sum()
    print(f"Top 20% obs (threshold={obs_threshold:.3f}) + p300+ union: {len(peaks_top_obs):,} "
          f"({n_top20_only:,} top-20-only + {n_p300_added:,} p300+ added)")
    make_figure(peaks_top_obs, f"{base}.top20pct_obs{ext}")


if __name__ == "__main__":
    main()
