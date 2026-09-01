#!/usr/bin/env python
"""Figure: does the ChromBPNet-convention 5' accessibility input change performance?

Paired within fold: each model retrained on `genomecov -bg -5` single-base insertion counts
against its counterpart on full-interval coverage, identical in every other respect.
The finding is the SPLIT -- it helps the model that has only accessibility, and does nothing
for the model that also has sequence.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import t as tdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label, figsize

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
COLOR = {"atac": "#2166AC", "multimodal": "#762A83"}
TC = tdist.ppf(0.975, df=4)
SETS = [("K562", "accs5p_k562"), ("GM12878", "accs5p_gm12878"), ("p300", "accs5p_p300")]


def mean_ci(v):
    v = np.asarray(v, float)
    return v.mean(), TC * v.std(ddof=1) / np.sqrt(len(v))


apply_rcparams()
fig, axes = plt.subplots(1, 2, figsize=figsize(columns=2, aspect=0.40))

for ax, stratum, title in [
        (axes[0], "top_quintile", "Top signal quintile"),
        (axes[1], "all", "All elements")]:
    W = 0.38
    for i, (cell, pre) in enumerate(SETS):
        d = pd.read_csv(f"{P}/results/{pre}_stratified_per_fold.tsv", sep="\t")
        d = d[d["stratum"] == stratum]
        piv = d.pivot(index="fold", columns="config", values="pearson")
        for j, mode in enumerate(["atac", "multimodal"]):
            a, b = f"{mode}_5PRIME", f"{mode}_OLD"
            if a not in piv or b not in piv:
                continue
            diff = (piv[a] - piv[b]).to_numpy()
            m, half = mean_ci(diff)
            x = i + (j - 0.5) * W
            ax.bar(x, m, width=W, color=COLOR[mode], zorder=2,
                   edgecolor="white", linewidth=0.4)
            ax.errorbar(x, m, yerr=half, fmt="none", ecolor="black", elinewidth=0.7,
                        capsize=1.5, capthick=0.7, zorder=4)
    ax.axhline(0, color="black", lw=0.75, zorder=3)
    ax.set_xticks(range(len(SETS)))
    ax.set_xticklabels([c for c, _ in SETS], fontsize=7)
    ax.set_ylabel("Δ Pearson $r$  (5′ − full-interval)")
    ax.set_title(title, fontsize=8)

axes[0].legend(handles=[Patch(facecolor=COLOR["atac"], label="ATAC-only"),
                        Patch(facecolor=COLOR["multimodal"], label="Sequence + ATAC")],
               frameon=False, fontsize=5.5, loc="upper right")
add_panel_label(axes[0], "a")
add_panel_label(axes[1], "b")
fig.tight_layout()
out = save_fig(fig, f"{P}/figures/fig9_atac_5prime_input.png")
print("Wrote", out, "and .pdf")

for cell, pre in SETS:
    d = pd.read_csv(f"{P}/results/{pre}_stratified_per_fold.tsv", sep="\t")
    for stratum in ["top_quintile", "all"]:
        piv = d[d["stratum"] == stratum].pivot(index="fold", columns="config", values="pearson")
        for mode in ["atac", "multimodal"]:
            m, h = mean_ci((piv[f"{mode}_5PRIME"] - piv[f"{mode}_OLD"]).to_numpy())
            print(f"{cell:<8} {stratum:<13} {mode:<11} {m:+.4f} [{m-h:+.4f}, {m+h:+.4f}]")
