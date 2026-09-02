#!/usr/bin/env python
"""Figures 10 and 11: the wide receptive field, and the profile-shape ceiling.

Fig 10 is a slopegraph rather than bars because the comparison is PAIRED within fold. Fold
sd on this project is 0.041-0.046, several times the effect, so bars with error bars would
show two overlapping clouds and hide a difference that is consistent fold by fold. The
slopes carry the actual evidence; the point of the figure is that the all-element slopes
rise together while the top-quintile slopes do not.

Fig 11 shows raw replicate agreement and the corrected ceiling side by side, because the raw
number invites the wrong conclusion: a 1 bp agreement of 0.02 is not a 0.02 ceiling. A model
is scored against the POOLED track, whose noise variance is halved, so the reachable value
is sqrt(2r/(1+r)) -- the same correction 0.3 applies to counts.

Colour follows the project convention (red = sequence-only input) in Fig 10. Fig 11 is
data-only with no model, so it uses cell type alone, GM12878 lighter as elsewhere.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, t as tdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
SEQ = "#B2182B"
K562_C, GM_C = "#33475b", "#8FA5BC"
TCRIT = tdist.ppf(0.975, df=4)
CELLS = [("k562", "K562"), ("gm12878", "GM12878")]
STRATA = [("overall_pearson", "All elements"), ("overall_pearson_topq", "Top quintile")]


def paired(df, metric):
    w = df[df["config"] == "sequence_WIDE"].sort_values("fold")[metric].to_numpy(float)
    n = df[df["config"] == "sequence_narrow"].sort_values("fold")[metric].to_numpy(float)
    d = w - n
    half = TCRIT * d.std(ddof=1) / np.sqrt(len(d))
    return n, w, d.mean(), half, ttest_rel(w, n).pvalue


def fig10():
    apply_rcparams()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax, (key, title) in zip(axes, CELLS):
        f = f"{P}/results/wide_seqonly_{key}_per_fold.tsv"
        if not os.path.exists(f):
            ax.text(.5, .5, "not yet scored", ha="center", transform=ax.transAxes)
            continue
        df = pd.read_csv(f, sep="\t")
        for j, (metric, slab) in enumerate(STRATA):
            x0, x1 = 3 * j, 3 * j + 1
            nar, wid, mu, half, p = paired(df, metric)
            for a, b in zip(nar, wid):
                ax.plot([x0, x1], [a, b], color="0.6", lw=0.6, marker="o", ms=2.2,
                        mfc="0.6", mec="none", zorder=1)
            ax.plot([x0, x1], [nar.mean(), wid.mean()], color=SEQ, lw=2.0,
                    marker="o", ms=4, zorder=2)
            star = "*" if p < 0.05 else "n.s."
            ax.text(x0 + .5, max(wid.max(), nar.max()) + .012,
                    f"{mu:+.3f} {star}", ha="center", color=SEQ, fontsize=6)
            ax.text(x0 + .5, -0.19, slab, ha="center", va="top", transform=
                    ax.get_xaxis_transform(), fontsize=7)
        ax.set_xticks([0, 1, 3, 4])
        ax.set_xticklabels(["1.1 kb", "4.2 kb"] * 2, fontsize=6, rotation=0)
        ax.set_xlim(-0.6, 4.6)
        ax.set_ylabel("Pearson $r$ vs observed H3K27ac")
        ax.set_title(title, fontsize=8)
    add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_fig(fig, f"{P}/figures/fig10_wide_receptive_field.png")
    print("wrote fig10")


def fig11():
    apply_rcparams()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    for key, title, col in [("k562", "K562", K562_C), ("gm12878", "GM12878", GM_C)]:
        f = f"{P}/results/profile_ceiling_binsize_{key}.tsv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep="\t")
        for strat, ls, mk in [("topq", "-", "o"), ("all", "--", "s")]:
            g = d[d["stratum"] == strat].sort_values("bin_bp")
            lab = f"{title}, {'top quintile' if strat == 'topq' else 'all elements'}"
            axes[0].plot(g["bin_bp"], g["shape_r_unstranded"], ls, color=col,
                         marker=mk, ms=3, lw=1.2, label=lab)
            axes[1].plot(g["bin_bp"], g["ceiling_unstranded"], ls, color=col,
                         marker=mk, ms=3, lw=1.2, label=lab)
    for ax, ylab, ttl in [
            (axes[0], "Replicate 1 vs replicate 2 $r$", "Raw agreement"),
            (axes[1], "Ceiling on the pooled track", "Corrected: $\\sqrt{2r/(1+r)}$")]:
        ax.set_xscale("log", base=2)
        ax.set_xticks([1, 5, 10, 25, 50, 100, 250])
        ax.set_xticklabels(["1", "5", "10", "25", "50", "100", "250"])
        ax.set_xlabel("Profile bin size (bp)")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, 1)
        ax.set_title(ttl, fontsize=8)
    axes[1].axhline(1.0, color="0.8", lw=0.6, ls=":")
    axes[0].legend(loc="upper left", fontsize=5.5)
    add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")
    fig.tight_layout()
    save_fig(fig, f"{P}/figures/fig11_profile_ceiling.png")
    print("wrote fig11")


if __name__ == "__main__":
    fig10()
    fig11()
