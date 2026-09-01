#!/usr/bin/env python
"""Figure 5: the residual objective across two cell types.

Supersedes the K562-only 3.3 figure. The point of the figure is the REPLICATION: the same
2x2 of {input} x {objective} plus an ATAC negative control, measured independently in K562
and GM12878, with panel c showing that the multimodal cost lands within 0.002 of itself
across two cell types with different targets, libraries and element sets.

Colour = input modality (project convention: blue ATAC, red sequence, purple sequence+ATAC).
Hatching = trained on the residual rather than the total signal.
Cell type = bar position within a group (K562 left, GM12878 right, GM12878 lighter).
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
COLOR = {"atac": "#2166AC", "sequence": "#B2182B", "multimodal": "#762A83"}
TCRIT = tdist.ppf(0.975, df=4)

BARS = [
    ("sequence5p",          "Seq\n(total)",       "sequence",   False),
    ("residual_sequence",   "Seq\n(resid.)",      "sequence",   True),
    ("multimodal5p",        "Seq+ATAC\n(total)",  "multimodal", False),
    ("residual_multimodal", "Seq+ATAC\n(resid.)", "multimodal", True),
    ("residual_atac",       "ATAC\n(ctrl)",       "atac",       True),
]
CELLS = [("K562", "residual_grid_per_fold.tsv", "residual_grid_residual_evaluation.tsv", 1.0),
         ("GM12878", "gm12878_residual_grid_per_fold.tsv",
          "gm12878_residual_grid_residual_evaluation.tsv", 0.55)]


def mean_ci(v):
    v = np.asarray(v, float)
    return v.mean(), TCRIT * v.std(ddof=1) / np.sqrt(len(v))


pf = {c: pd.read_csv(f"{P}/results/{f}", sep="\t") for c, f, _, _ in CELLS}
pooled = {c: pd.read_csv(f"{P}/results/{g}", sep="\t") for c, _, g, _ in CELLS}

apply_rcparams()
fig, axes = plt.subplots(1, 3, figsize=figsize(columns=2, aspect=0.40))
rng = np.random.default_rng(0)
W = 0.38

# --- a: residual r, both cell types -----------------------------------------
ax = axes[0]
for i, (key, lab, mode, hatched) in enumerate(BARS):
    for j, (cell, _, _, alpha) in enumerate(CELLS):
        v = pf[cell].loc[pf[cell]["config"] == key, "residual_pearson"].to_numpy()
        if not len(v):
            continue
        m, half = mean_ci(v)
        x = i + (j - 0.5) * W
        ax.bar(x, m, width=W, color=COLOR[mode], alpha=alpha, zorder=2,
               edgecolor="white", linewidth=0.4, hatch="///" if hatched else None)
        ax.errorbar(x, m, yerr=half, fmt="none", ecolor="black", elinewidth=0.7,
                    capsize=1.5, capthick=0.7, zorder=4)
        ax.scatter(np.full(len(v), x) + (rng.random(len(v)) - 0.5) * 0.14, v,
                   s=3, color="black", alpha=0.7, linewidths=0, zorder=5)
ax.axhline(0, color="black", lw=0.6, zorder=3)
ax.set_xticks(range(len(BARS)))
ax.set_xticklabels([b[1].replace("\n", " ") for b in BARS], fontsize=5, rotation=30, ha="right")
ax.set_ylabel("Residual Pearson $r$")
ax.set_title("Residual beyond ATAC", fontsize=8)
ax.legend(handles=[Patch(facecolor="#666666", label="K562"),
                   Patch(facecolor="#666666", alpha=0.55, label="GM12878")],
          frameon=False, fontsize=5, loc="upper left")
add_panel_label(ax, "a")

# --- b: |residual| quintiles ------------------------------------------------
ax = axes[1]
for key, lab, mode, hatched in BARS:
    for cell, _, _, alpha in CELLS:
        s = pooled[cell][(pooled[cell]["config"] == key)
                         & pooled[cell]["stratum"].str.startswith("abs_resid_q")]
        if not len(s):
            continue
        ax.plot(range(1, len(s) + 1), s["residual_pearson"], marker="o", ms=2.5, lw=0.9,
                color=COLOR[mode], alpha=alpha,
                ls=(0, (3, 1.5)) if hatched else "-")
ax.set_xlabel("|true residual| quintile\n(5 = ATAC most wrong)")
ax.set_ylabel("Residual Pearson $r$")
ax.set_xticks(range(1, 6))
ax.set_title("Same ordering, both lines", fontsize=8)
add_panel_label(ax, "b")

# --- c: the replication -- paired effect of the objective -------------------
ax = axes[2]
PAIRS = [("residual_sequence", "sequence5p", "sequence", "Sequence"),
         ("residual_multimodal", "multimodal5p", "multimodal", "Seq + ATAC")]
for i, (a_, b_, mode, lab) in enumerate(PAIRS):
    for j, (cell, _, _, alpha) in enumerate(CELLS):
        d = pf[cell].pivot(index="fold", columns="config", values="residual_pearson")
        diff = (d[a_] - d[b_]).to_numpy()
        m, half = mean_ci(diff)
        x = i + (j - 0.5) * W
        ax.bar(x, m, width=W, color=COLOR[mode], alpha=alpha, zorder=2,
               edgecolor="white", linewidth=0.4)
        ax.errorbar(x, m, yerr=half, fmt="none", ecolor="black", elinewidth=0.7,
                    capsize=1.5, capthick=0.7, zorder=4)
        ax.scatter(np.full(len(diff), x) + (rng.random(len(diff)) - 0.5) * 0.14, diff,
                   s=3, color="black", alpha=0.7, linewidths=0, zorder=5)
        off = 0.055 if m > 0 else -0.055
        ax.text(x, m + half + off if m > 0 else m - half + off, f"{m:+.3f}",
                ha="center", va="bottom" if m > 0 else "top", fontsize=5)
ax.axhline(0, color="black", lw=0.75, zorder=3)
ax.set_xticks(range(len(PAIRS)))
ax.set_xticklabels([p[3] for p in PAIRS], fontsize=6)
ax.set_ylabel("Effect of residual objective\n(paired over folds)")
ax.set_title("Effect of the objective", fontsize=8)
ax.set_ylim(-0.13, 0.45)
add_panel_label(ax, "c")

fig.tight_layout()
out = save_fig(fig, f"{P}/figures/fig5_residual_grid_both.png")
print("Wrote", out, "and .pdf")

print("\n--- legend values ---")
for cell, _, _, _ in CELLS:
    d = pf[cell].pivot(index="fold", columns="config", values="residual_pearson")
    for key, lab, _, _ in BARS:
        if key in d:
            m, h = mean_ci(d[key])
            print(f"{cell:<9} {key:<22} {m:.3f} [{m-h:.3f}, {m+h:.3f}]")
    for a_, b_, _, lab in PAIRS:
        m, h = mean_ci((d[a_] - d[b_]).to_numpy())
        print(f"{cell:<9} DIFF {a_} - {b_}: {m:+.4f} [{m-h:+.4f}, {m+h:+.4f}]")
