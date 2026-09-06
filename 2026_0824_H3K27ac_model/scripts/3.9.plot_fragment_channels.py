#!/usr/bin/env python
"""Figure 12: fragment-size-stratified accessibility channels.

Paired slopegraph for the same reason as Fig 10 -- fold sd is 0.041-0.046, several times the
effect, so the evidence is that every fold moves the same way rather than that two clouds
separate.

The input is a strict superset of the flat baseline: the four length bins sum to the flat 5'
track exactly (0 discrepancy in 545.7 million insertions, 0.24.validate_atac_fragment_channels.py),
so a gain is added information and not a changed input definition.
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
MM = "#762A83"
TCRIT = tdist.ppf(0.975, df=4)
STRATA = [("overall_pearson", "All elements"), ("overall_pearson_topq", "Top quintile")]

f = f"{P}/results/fragchan_k562_per_fold.tsv"
if not os.path.exists(f):
    raise SystemExit(f"{f} not found; run 2.17 fragchan first")
df = pd.read_csv(f, sep="\t")

apply_rcparams()
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
for ax, (metric, slab) in zip(axes, STRATA):
    w = df[df["config"] == "multimodal_FRAGCHAN"].sort_values("fold")[metric].to_numpy(float)
    n = df[df["config"] == "multimodal_flat"].sort_values("fold")[metric].to_numpy(float)
    d = w - n
    mu = d.mean()
    half = TCRIT * d.std(ddof=1) / np.sqrt(len(d))
    pv = ttest_rel(w, n).pvalue
    for a, b in zip(n, w):
        ax.plot([0, 1], [a, b], color="0.6", lw=0.6, marker="o", ms=2.4,
                mfc="0.6", mec="none", zorder=1)
    ax.plot([0, 1], [n.mean(), w.mean()], color=MM, lw=2.0, marker="o", ms=4, zorder=2)
    star = "*" if pv < 0.05 else "n.s."
    ax.text(0.5, max(w.max(), n.max()) + 0.002, f"{mu:+.4f} {star}", ha="center",
            color=MM, fontsize=6.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Flat ATAC\n(1 channel)", "Fragment-stratified\n(5 channels)"],
                       fontsize=6.5)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("Pearson $r$ vs observed H3K27ac")
    ax.set_title(slab, fontsize=8)
    print(f"{metric}: {n.mean():.4f} -> {w.mean():.4f}  "
          f"delta {mu:+.4f} [{mu-half:+.4f}, {mu+half:+.4f}] p={pv:.4f}")
add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")
fig.tight_layout()
annotate_n_fig(fig, n_label(df))
save_fig(fig, f"{P}/figures/fig12_fragment_channels.png")
print("wrote fig12")
