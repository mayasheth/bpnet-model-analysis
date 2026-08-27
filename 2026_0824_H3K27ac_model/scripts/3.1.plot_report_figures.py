#!/usr/bin/env python
"""
Report figures for the H3K27ac modeling analysis, Nature-formatted.

Reads the per-fold evaluation tables written by 2.2.evaluate_stratified.py and plots
the mean across chromosome-holdout folds with a t-based 95% confidence interval, with
the individual fold values overlaid as points.

Why mean +/- CI rather than one pooled correlation: each fold is an independent replicate
of the train-and-evaluate procedure, so the spread across folds is the only uncertainty
estimate available. A pooled correlation is a single number with no error bar, and it can
also be biased when folds differ in mean or scale (a fold whose predictions are shifted
relative to another contributes off-diagonal spread that belongs to neither fold). With
measured run-to-run variance around 0.018, several comparisons in this analysis differ by
less than that, and a pooled number cannot say so.

Colour encodes input modality throughout: blue = ATAC only, red = sequence only,
purple = sequence + ATAC.

Usage:
  python 3.1.plot_report_figures.py --results results --figdir figures
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, figsize, add_panel_label

# Modality colours, as requested: blue ATAC, red sequence, purple both.
COLOR = {"atac": "#2166AC", "sequence": "#B2182B", "multimodal": "#762A83"}
LABEL = {"atac": "ATAC only", "sequence": "Sequence only",
         "multimodal": "Sequence\n+ ATAC"}
MODE_ORDER = ["sequence", "atac", "multimodal"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="results")
    p.add_argument("--figdir", default="figures")
    return p.parse_args()


def ceiling(r):
    """Bound on achievable correlation from an inter-replicate r.

    Spearman-Brown for the two-replicate merged target, then sqrt because a model
    predicts the expected signal rather than a noisy draw of it.
    """
    rel = 2 * r / (1 + r)
    return rel ** 0.5


def load_per_fold(results, prefix):
    path = os.path.join(results, f"{prefix}stratified_per_fold.tsv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def mean_ci(values):
    v = np.asarray(values, dtype=float)
    n = len(v)
    m = v.mean()
    if n < 2:
        return m, np.nan
    sd = v.std(ddof=1)
    return m, float(stats.t.ppf(0.975, n - 1)) * sd / np.sqrt(n)


def bars_with_folds(ax, groups, ceil_val=None, ylabel=None, title=None):
    """One bar per group: mean across folds, 95% CI, individual folds as points.

    groups: list of (x_label, mode, [per-fold values])
    """
    x = np.arange(len(groups))
    for i, (_, mode, vals) in enumerate(groups):
        m, half = mean_ci(vals)
        ax.bar(i, m, width=0.62, color=COLOR[mode], zorder=2,
               edgecolor="none")
        if np.isfinite(half):
            ax.errorbar(i, m, yerr=half, fmt="none", ecolor="black",
                        elinewidth=0.75, capsize=2.5, capthick=0.75, zorder=4)
        # individual folds, jittered so overlapping values stay visible
        rng = np.random.default_rng(0)
        jit = rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jit, vals, s=4.5, color="black",
                   alpha=0.75, linewidths=0, zorder=5)
    if ceil_val is not None:
        ax.axhline(ceil_val, color="black", lw=0.75, ls=(0, (4, 2)), zorder=1)
        ax.text(len(groups) - 0.45, ceil_val, "ceiling", ha="right", va="bottom",
                fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylim(0, 1.0)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


def fig_three_mode(pf, results, figdir):
    """Headline: three modalities, all elements vs top quintile."""
    ceil_path = os.path.join(results, "fiveprime_replicate_ceiling_by_window.tsv")
    c = pd.read_csv(ceil_path, sep="\t")
    row = c[c["half_window"] == 500].iloc[0]
    ceils = {"all": ceiling(row["pearson_all"]),
             "top_quintile": ceiling(row["pearson_top_quintile"])}

    fig, axes = plt.subplots(1, 2, figsize=figsize(columns=2, aspect=0.40),
                             sharey=True)
    for ax, stratum, panel, ttl in zip(
            axes, ["all", "top_quintile"], "ab",
            ["All elements", "Top H3K27ac quintile"]):
        groups = []
        for mode in MODE_ORDER:
            sub = pf[(pf["mode"] == mode) & (pf["stratum"] == stratum)]
            if sub.empty:
                continue
            groups.append((LABEL[mode], mode, sub["pearson"].tolist()))
        bars_with_folds(ax, groups, ceil_val=ceils[stratum],
                        ylabel="Pearson r (log counts)" if panel == "a" else None,
                        title=ttl)
        add_panel_label(ax, panel)
    fig.tight_layout()
    return save_fig(fig, os.path.join(figdir, "fig1_three_mode_comparison.png"))


def fig_p300_vs_h3k27ac(pf, figdir):
    """p300 and H3K27ac side by side, each target explicitly labelled."""
    targets = [("p300", "p300_"), ("H3K27ac", "h3k27ac_")]
    fig, axes = plt.subplots(1, 2, figsize=figsize(columns=2, aspect=0.44),
                             sharey=True)
    for ax, (tname, tprefix), panel in zip(axes, targets, "ab"):
        groups = []
        for mode in MODE_ORDER:
            sub = pf[(pf["mode"] == mode) & (pf["stratum"] == "top_quintile")
                     & pf["config"].str.startswith(tprefix)]
            if sub.empty:
                continue
            groups.append((LABEL[mode], mode, sub["pearson"].tolist()))
        if not groups:
            continue
        bars_with_folds(ax, groups,
                        ylabel="Pearson r (log counts)" if panel == "a" else None,
                        title=f"{tname} target")
        # the quantity being compared: the sequence margin over accessibility
        by = {m: np.mean(v) for _, m, v in groups}
        if "atac" in by and "multimodal" in by:
            delta = by["multimodal"] - by["atac"]
            ax.annotate("", xy=(2, by["multimodal"]), xytext=(2, by["atac"]),
                        arrowprops=dict(arrowstyle="<->", lw=0.75, color="#333"))
            ax.text(2.28, (by["multimodal"] + by["atac"]) / 2,
                    f"+{delta:.3f}", fontsize=6.5, va="center", color="#333")
        ax.set_xlim(-0.6, 2.75)
        add_panel_label(ax, panel)
    fig.suptitle("Sequence margin over accessibility, top signal quintile", y=1.02)
    fig.tight_layout()
    return save_fig(fig, os.path.join(figdir, "fig4_p300_vs_h3k27ac.png"))


def fig_transfer(pf_gm, pf_k562, results, figdir):
    """K562 in-cell-type vs GM12878 transfer, as a fraction of each ceiling."""
    kc = pd.read_csv(os.path.join(results, "replicate_ceiling_by_window.tsv"),
                     sep="\t")
    gc = pd.read_csv(os.path.join(results,
                                  "gm12878_frag250_replicate_ceiling_by_window.tsv"),
                     sep="\t")
    k_ceil = ceiling(kc[kc["half_window"] == 500]["pearson_top_quintile"].iloc[0])
    g_ceil = ceiling(gc[gc["half_window"] == 500]["pearson_top_quintile"].iloc[0])

    fig, ax = plt.subplots(figsize=figsize(columns=1, aspect=0.78))
    width = 0.38
    for j, mode in enumerate(MODE_ORDER):
        k = pf_k562[(pf_k562["mode"] == mode)
                    & (pf_k562["stratum"] == "top_quintile")
                    & pf_k562["config"].str.startswith("h3k27ac_")]
        g = pf_gm[(pf_gm["mode"] == mode) & (pf_gm["stratum"] == "top_quintile")]
        if k.empty or g.empty:
            continue
        kf = 100 * k["pearson"].to_numpy() / k_ceil
        gf = 100 * g["pearson"].to_numpy() / g_ceil
        for off, vals, hatch in ((-width / 2, kf, None), (width / 2, gf, "///")):
            m, half = mean_ci(vals)
            ax.bar(j + off, m, width=width, color=COLOR[mode], hatch=hatch,
                   edgecolor="white", linewidth=0.4, zorder=2)
            if np.isfinite(half):
                ax.errorbar(j + off, m, yerr=half, fmt="none", ecolor="black",
                            elinewidth=0.75, capsize=2, capthick=0.75, zorder=4)
    ax.set_xticks(range(len(MODE_ORDER)))
    ax.set_xticklabels([LABEL[m] for m in MODE_ORDER])
    ax.set_ylabel("Percent of achievable ceiling")
    ax.set_ylim(0, 100)
    ax.set_title("K562 (solid) vs GM12878 transfer (hatched)")
    fig.tight_layout()
    return save_fig(fig, os.path.join(figdir, "fig6_transfer.png"))


def main():
    args = parse_args()
    apply_rcparams()
    os.makedirs(args.figdir, exist_ok=True)
    written = []

    pf5 = load_per_fold(args.results, "fiveprime_")
    if pf5 is not None:
        written.append(fig_three_mode(pf5, args.results, args.figdir))
    else:
        print("note: fiveprime_stratified_per_fold.tsv not found; skipping Fig 1")

    pfc = load_per_fold(args.results, "p300_vs_h3k27ac_")
    if pfc is not None:
        written.append(fig_p300_vs_h3k27ac(pfc, args.figdir))
    else:
        print("note: p300_vs_h3k27ac per-fold table not found; skipping Fig 4")

    pfg = load_per_fold(args.results, "gm12878_transfer_")
    if pfg is not None and pfc is not None:
        written.append(fig_transfer(pfg, pfc, args.results, args.figdir))
    else:
        print("note: transfer per-fold table not found; skipping Fig 6")

    for w in written:
        print(f"Wrote {w}")


if __name__ == "__main__":
    main()
