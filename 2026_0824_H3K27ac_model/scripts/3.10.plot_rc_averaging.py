#!/usr/bin/env python
"""Figure 13: what test-time reverse-complement averaging is worth.

A forest plot of paired within-fold differences, ordered by how much each model relies on
sequence. That ordering is the point: RC averaging corrects residual strand asymmetry in the
SEQUENCE branch, so the gain should fall away as sequence matters less, and the ATAC-only
models -- whose input is unstranded coverage -- act as a built-in negative control.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, t as tdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label, n_label, annotate_n_fig

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
TCRIT = tdist.ppf(0.975, df=4)
# ordered by decreasing reliance on sequence; colour = input modality convention
ROWS = [("sequence5p",        "Sequence only",              "#B2182B"),
        ("residual_sequence", "Sequence, residual",         "#B2182B"),
        ("multimodal5p",      "Sequence + ATAC",            "#762A83"),
        ("residual_multimodal", "Sequence + ATAC, residual", "#762A83"),
        ("atac5p",            "ATAC only",                  "#2166AC"),
        ("residual_atac",     "ATAC only, residual",        "#2166AC")]

plain = pd.read_csv(f"{P}/results/prof_residual_grid_per_fold.tsv", sep="\t")
rc = pd.read_csv(f"{P}/results/prof_rc_residual_grid_per_fold.tsv", sep="\t")

apply_rcparams()
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
METRICS = [("overall_pearson_topq", "Top quintile, counts"),
           ("profile_pearson_topq", "Top quintile, profile shape")]
for ax, (metric, title) in zip(axes, METRICS):
    for i, (cfg, lab, col) in enumerate(ROWS):
        y = len(ROWS) - 1 - i
        a = rc[rc["config"] == cfg].sort_values("fold")[metric].to_numpy(float)
        b = plain[plain["config"] == cfg].sort_values("fold")[metric].to_numpy(float)
        d = a - b
        mu = d.mean()
        half = TCRIT * d.std(ddof=1) / np.sqrt(len(d))
        pv = ttest_rel(a, b).pvalue
        ax.errorbar(mu, y, xerr=half, fmt="o", ms=4, color=col, lw=1.2, capsize=2)
        ax.text(mu + half + 0.0008, y, "*" if pv < 0.05 else "n.s.",
                va="center", fontsize=6, color=col)
        if ax is axes[0]:
            ax.text(-0.0065, y, lab, ha="right", va="center", fontsize=6.5)
        print(f"{metric:<22} {cfg:<20} {mu:+.4f} [{mu-half:+.4f}, {mu+half:+.4f}] p={pv:.4f}")
    ax.axvline(0, color="0.7", lw=0.7, ls="--")
    ax.set_yticks([])
    ax.set_ylim(-0.7, len(ROWS) - 0.3)
    ax.set_xlabel("RC-averaged − single-pass Pearson $r$")
    ax.set_title(title, fontsize=8)
axes[0].set_xlim(-0.006, 0.024)
add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")
fig.tight_layout(rect=(0.14, 0, 1, 1))
annotate_n_fig(fig, n_label(plain))
save_fig(fig, f"{P}/figures/fig13_rc_averaging.png")
print("wrote fig13")
