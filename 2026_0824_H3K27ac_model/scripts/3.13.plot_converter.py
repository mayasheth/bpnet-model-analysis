#!/usr/bin/env python
"""Figure 1 for report 4: the ATAC->DNase converter, in-cell and transferred.

Panels a and b are the bin-size curve rather than a single number, because the converter
exists to feed a model that reads base resolution. A converter good only at 50 bp would
produce a smooth bump a flat fill could match, and one number would hide which of those we
have. The x axis is log because the bins span 1 to 250.

Three lines and a ceiling per panel, and each is a different kind of object:
  grey    the OBSERVED ATAC track, no model. This is the comparator that decides deployment,
          because it is what currently sits in the downstream accessibility slot.
  blue    the ATAC-only model, project convention for an accessibility-only input. The
          control for whether sequence is doing anything.
  purple  the converter, project convention for sequence + ATAC.
  black   the DNase inter-replicate ceiling, dashed, computed on the same held-out windows.

Panel c is the paired within-fold delta at 1 bp, which is the resolution the downstream
branch reads. Paired because fold sd on this project runs 0.041-0.046.

Top quintile only. The all-elements stratum is in the tables and the report's Methods; a
figure carrying both strata for two cell types across four series is unreadable.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import (apply_rcparams, save_fig, add_panel_label, figsize,
                          n_parts, fold_n_range, annotate_n_fig)

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
ATAC_C, MM_C, OBS_C = "#2166AC", "#762A83", "#8C8C8C"
CELLS = [("k562", "K562 (in-cell)"), ("gm12878", "GM12878 (transferred)")]
STRAT = "topq"

base = pd.read_csv(f"{P}/results/atac_vs_dnase_profile_baseline.tsv", sep="\t")
def load(lab, cell):
    return pd.read_csv(f"{P}/results/converter_profile_{lab}_{cell}.tsv", sep="\t")


def curve(df, col):
    g = df[df["stratum"] == STRAT].groupby("bin_bp")[col].mean()
    return g.index.to_numpy(float), g.to_numpy(float)


apply_rcparams()
fig, axes = plt.subplots(1, 3, figsize=figsize(columns=2, aspect=0.36))

for ax, (cell, title) in zip(axes[:2], CELLS):
    conv, ao = load("converter5f", cell), load("ataconly5f", cell)
    b = base[(base["cell"] == cell) & (base["stratum"] == STRAT)].groupby("bin_bp")[
        "shape_r_observed_atac_vs_dnase"].mean()
    x, y_conv = curve(conv, "shape_r_model")
    _, y_ao = curve(ao, "shape_r_model")
    _, y_ceil = curve(conv, "ceiling")
    ax.plot(x, y_ceil, color="black", lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.text(x[-1], y_ceil[-1], " ceiling", fontsize=5.5, va="center", color="#333333")
    ax.plot(b.index.to_numpy(float), b.to_numpy(float), color=OBS_C, lw=1.6,
            marker="o", ms=3, label="observed ATAC")
    ax.plot(x, y_ao, color=ATAC_C, lw=1.6, marker="o", ms=3, label="ATAC only, model")
    ax.plot(x, y_conv, color=MM_C, lw=1.6, marker="o", ms=3, label="converter (seq + ATAC)")
    ax.set_xscale("log")
    ax.set_xticks([1, 5, 10, 25, 50, 100, 250])
    ax.set_xticklabels(["1", "5", "10", "25", "50", "100", "250"], fontsize=5.5)
    ax.set_xlabel("bin size (bp)")
    ax.set_ylim(0, 1.05)
    ax.set_title(title, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("shape $r$ with observed DNase\n(top quintile)")
axes[0].legend(loc="lower right", frameon=False, fontsize=5.5, handlelength=1.4,
               labelspacing=0.3)
add_panel_label(axes[0], "a"); add_panel_label(axes[1], "b")

# --- c: paired delta at 1 bp, the resolution the downstream branch reads ------
ax = axes[2]
pair = pd.read_csv(f"{P}/results/converter_vs_atac_paired.tsv", sep="\t")
p1 = pair[(pair["bin_bp"] == 1) & (pair["stratum"] == STRAT)]
rows = [("K562", "converter5f", MM_C), ("K562", "ataconly5f", ATAC_C),
        ("GM12878", "converter5f", MM_C), ("GM12878", "ataconly5f", ATAC_C)]
ypos, labs = [], []
for i, (cellname, lab, col) in enumerate(rows):
    cell = "k562" if cellname == "K562" else "gm12878"
    r = p1[(p1["cell"] == cell) & (p1["label"] == lab)]
    if r.empty:
        continue
    r = r.iloc[0]
    y = len(rows) - 1 - i
    ax.errorbar(r["delta"], y, xerr=[[r["delta"] - r["ci_lo"]], [r["ci_hi"] - r["delta"]]],
                fmt="o", ms=4, color=col, lw=1.3, capsize=2)
    ax.text(r["delta"], y + 0.28, f"{r['delta']:+.3f}", ha="center", fontsize=5.5, color=col)
    ypos.append(y)
    labs.append(f"{cellname}, " + ("converter" if lab == "converter5f" else "ATAC only"))
ax.axvline(0, color="0.6", lw=0.7, ls="--")
ax.set_yticks(ypos); ax.set_yticklabels(labs, fontsize=6)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xlabel("gain over the observed ATAC track\n(1 bp shape $r$, paired within fold)")
ax.set_title("1 bp, top quintile", fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
add_panel_label(ax, "c")

fig.tight_layout()
# Counts differ by cell type and come from the model tables, which record the elements
# actually scored after window extraction rather than the size of the candidate set.
_c = load("converter5f", "k562"); _g = load("converter5f", "gm12878")
annotate_n_fig(fig, n_parts(
    ("K562, top quintile", fold_n_range(_c[_c["stratum"] == STRAT], n_col="n_elements")),
    ("GM12878, top quintile", fold_n_range(_g[_g["stratum"] == STRAT], n_col="n_elements"))))
out = save_fig(fig, f"{P}/figures/fig17_converter.png")
print("wrote", out)

print("\n1 bp, top quintile, mean over 5 folds:")
for cell, title in CELLS:
    conv, ao = load("converter5f", cell), load("ataconly5f", cell)
    b = base[(base["cell"] == cell) & (base["stratum"] == STRAT) & (base["bin_bp"] == 1)]
    f = lambda d, c: d[(d["stratum"] == STRAT) & (d["bin_bp"] == 1)][c].mean()
    print(f"  {title:<26} ATAC {b['shape_r_observed_atac_vs_dnase'].mean():.3f}  "
          f"ATAC-only {f(ao,'shape_r_model'):.3f}  converter {f(conv,'shape_r_model'):.3f}  "
          f"ceiling {f(conv,'ceiling'):.3f}")
