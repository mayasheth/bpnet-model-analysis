#!/usr/bin/env python
"""Figures 10 and 11: the wide receptive field, and the profile-shape ceiling.

Fig 10 is a slopegraph rather than bars because the comparison is PAIRED within fold. Fold
sd on this project is 0.041-0.046, several times the effect, so bars with error bars would
show two overlapping clouds and hide a difference that is consistent fold by fold.

The figure covers BOTH input modes, because they disagree and an earlier sequence-only
version of this figure supported the wrong conclusion. Sequence-only gains on all elements
and nothing on the top quintile; multimodal gains on the top quintile in both cell types.
Plotting one mode invited reading the sequence null as the whole answer.

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
from nature_style import (apply_rcparams, save_fig, add_panel_label, n_label, n_parts,
                          fold_n_range, annotate_n_fig)

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
SEQ = "#B2182B"
K562_C, GM_C = "#33475b", "#8FA5BC"
TCRIT = tdist.ppf(0.975, df=4)
CELLS = [("k562", "K562"), ("gm12878", "GM12878")]
STRATA = [("overall_pearson", "All elements"), ("overall_pearson_topq", "Top quintile")]


MODES = [("sequence", "Sequence only", "#B2182B"),
         ("multimodal", "Sequence + ATAC", "#762A83")]


def paired(df, metric, mode):
    w = df[df["config"] == f"{mode}_WIDE"].sort_values("fold")[metric].to_numpy(float)
    n = df[df["config"] == f"{mode}_narrow"].sort_values("fold")[metric].to_numpy(float)
    d = w - n
    half = TCRIT * d.std(ddof=1) / np.sqrt(len(d))
    return n, w, d.mean(), half, ttest_rel(w, n).pvalue


def fig10():
    apply_rcparams()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    letters = iter("abcd")
    # Collected per cell type rather than left to the loop variable. This annotation used
    # to read `n_label(df)` after the loops had finished, so a four-panel figure covering
    # two cell types with different element sets (K562 9,298-14,392, GM12878 10,463-14,362)
    # was labelled with whichever one the loop happened to read last -- and raised NameError
    # outright if the first tables were missing.
    seen = {}
    for i, (mode, mlab, colour) in enumerate(MODES):
        for j, (key, title) in enumerate(CELLS):
            ax = axes[i][j]
            f = f"{P}/results/wide_{key}_per_fold.tsv"
            if not os.path.exists(f):
                ax.text(.5, .5, "not yet scored", ha="center", transform=ax.transAxes)
                continue
            df = pd.read_csv(f, sep="\t")
            seen[title] = df
            if f"{mode}_WIDE" not in set(df["config"]):
                ax.text(.5, .5, "not yet scored", ha="center", transform=ax.transAxes)
                continue
            for k, (metric, slab) in enumerate(STRATA):
                x0, x1 = 3 * k, 3 * k + 1
                nar, wid, mu, half, pv = paired(df, metric, mode)
                for aa, bb in zip(nar, wid):
                    ax.plot([x0, x1], [aa, bb], color="0.6", lw=0.6, marker="o", ms=2.2,
                            mfc="0.6", mec="none", zorder=1)
                ax.plot([x0, x1], [nar.mean(), wid.mean()], color=colour, lw=2.0,
                        marker="o", ms=4, zorder=2)
                star = "*" if pv < 0.05 else "n.s."
                ax.text(x0 + .5, max(wid.max(), nar.max()) + .010,
                        f"{mu:+.3f} {star}", ha="center", color=colour, fontsize=6)
                ax.text(x0 + .5, -0.17, slab, ha="center", va="top",
                        transform=ax.get_xaxis_transform(), fontsize=7)
            ax.set_xticks([0, 1, 3, 4])
            ax.set_xticklabels(["1.1 kb", "4.2 kb"] * 2, fontsize=6, rotation=0)
            ax.set_xlim(-0.6, 4.6)
            ax.set_ylabel("Pearson $r$")
            ax.set_title(f"{mlab} — {title}", fontsize=8)
            add_panel_label(ax, next(letters))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    # Every panel draws both strata side by side, so both counts are stated for each cell
    # type. n_topq is the table's own column (2.30), not n/5: `overall_pearson_topq` is
    # computed on a VALUE threshold, obs >= quantile(obs, 0.8), which ties push past a fifth.
    if seen:
        annotate_n_fig(fig, n_parts(*[
            part
            for title, d in seen.items()
            for part in ((f"{title} all", fold_n_range(d)),
                         (f"{title} top q.", fold_n_range(d, n_col="n_topq")))]))
    save_fig(fig, f"{P}/figures/fig10_wide_receptive_field.png")
    print("wrote fig10")


def fig11():
    """Profile ceilings for both marks.

    DNase is overlaid on the same axes as H3K27ac because the comparison is the point: the
    two differ by roughly 4x at 1 bp, which is what decides whether a model's profile head
    can produce a usable base-resolution track for that mark. Only the top quintile is drawn
    for DNase, since showing both strata for both marks puts eight lines on one panel.
    """
    apply_rcparams()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    series = [("k562", "K562 H3K27ac", K562_C, ("topq", "all")),
              ("gm12878", "GM12878 H3K27ac", GM_C, ("topq", "all")),
              ("k562_dnase", "K562 DNase", "#2b8cbe", ("topq",)),
              ("gm12878_dnase", "GM12878 DNase", "#7bccc4", ("topq",))]
    for key, title, col, strata in series:
        f = f"{P}/results/profile_ceiling_binsize_{key}.tsv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, sep="\t")
        for strat in strata:
            ls, mk = ("-", "o") if strat == "topq" else ("--", "s")
            g = d[d["stratum"] == strat].sort_values("bin_bp")
            lab = title if strat == "topq" else f"{title}, all elements"
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
    axes[0].legend(loc="upper left", fontsize=5, ncol=1, frameon=False,
               handlelength=1.6, labelspacing=0.25)
    add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")
    fig.tight_layout()
    # n comes from the table's own n_elements column, not from len(d), which is the number of
    # (bin size, stratum) rows and was silently reporting "n = 14 elements".
    _ns = []
    for key, title, _c, _s in series:
        f = f"{P}/results/profile_ceiling_binsize_{key}.tsv"
        if not os.path.exists(f):
            continue
        _g = pd.read_csv(f, sep="\t")
        _g = _g[(_g["stratum"] == "topq") & (_g["bin_bp"] == 1)]
        if len(_g):
            _ns.append(f"{title} {int(_g['n_elements'].iloc[0]):,}")
    annotate_n_fig(fig, "elements with non-flat profiles in both replicates, top quintile: "
                        + "; ".join(_ns))
    save_fig(fig, f"{P}/figures/fig11_profile_ceiling.png")
    print("wrote fig11")


if __name__ == "__main__":
    fig10()
    fig11()
