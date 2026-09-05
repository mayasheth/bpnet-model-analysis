#!/usr/bin/env python
"""Figure: CRISPR benchmark AUPRC for the eleven predicted-activity arms.

A forest plot rather than bars, because the CIs are the point: they are ~+/-0.05 wide while
the effects under discussion are ~0.02, so any presentation that hides them overstates the
result. The floor (ATAC only) and ceiling (observed H3K27ac) are drawn as vertical lines so
each predicted arm reads as a position between them.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
SRC = ("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/"
       "workflow/results/2026_0904_predicted_activity/performance_summary.txt")

d = pd.read_csv(SRC, sep="\t")
d["arm"] = d["pred_uid"].str.replace(r"\.ABC\.Score$", "", regex=True)

FLOOR, CEIL = "K562_ATAC_only", "K562_ATAC_H3K27ac_element"
# top to bottom: ceiling, K562-trained, GM12878-trained, floor, sequence arms, distance
ORDER = [
    (CEIL, "Observed H3K27ac", "#404040"),
    ("k27only_k562_multimodal", "Predicted alone (K562 multimodal)", "#238b45"),
    ("pred_k562_multimodal", "ATAC x predicted (K562 multimodal)", "#2171b5"),
    ("k27only_k562_atac", "Predicted alone (K562 ATAC)", "#74c476"),
    ("pred_k562_atac", "ATAC x predicted (K562 ATAC)", "#6baed6"),
    (FLOOR, "ATAC only", "#bdbdbd"),
    ("pred_gm12878_atac", "ATAC x predicted (GM12878 ATAC)", "#fcae91"),
    ("pred_gm12878_multimodal", "ATAC x predicted (GM12878 multimodal)", "#cb181d"),
    ("baseline.distToTSS", "Distance to TSS", "#c5cad7"),
    ("pred_k562_sequence", "ATAC x predicted (K562 sequence)", "#9ecae1"),
    ("pred_gm12878_sequence", "ATAC x predicted (GM12878 sequence)", "#fb6a4a"),
    ("k27only_k562_sequence", "Predicted alone (K562 sequence)", "#c7e9c0"),
]
got = {r["arm"]: r for _, r in d.iterrows()}
missing = [a for a, _l, _c in ORDER if a not in got]
assert not missing, f"missing arms in summary: {missing}"

floor_v = float(got[FLOOR]["AUPRC"])
ceil_v = float(got[CEIL]["AUPRC"])

apply_rcparams()
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.axvspan(floor_v, ceil_v, color="0.93", zorder=0)
ax.axvline(floor_v, color="0.55", lw=0.9, ls="--", zorder=1)
ax.axvline(ceil_v, color="0.25", lw=0.9, ls="--", zorder=1)

for i, (arm, lab, colour) in enumerate(ORDER):
    y = len(ORDER) - 1 - i
    r = got[arm]
    v, lo, hi = float(r["AUPRC"]), float(r["AUPRC_lowerCi"]), float(r["AUPRC_upperCi"])
    ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", ms=4.5, color=colour,
                lw=1.3, capsize=2.5, zorder=3)
    ax.text(0.235, y, lab, ha="right", va="center", fontsize=6.5)
    print(f"{lab:<42} {v:.3f} [{lo:.3f}, {hi:.3f}]")

ax.set_yticks([])
ax.set_ylim(-0.8, len(ORDER) - 0.2)
ax.set_xlim(0.235, 0.60)
ax.set_xlabel("AUPRC against the CRISPR benchmark (95% CI)")
ax.text(ceil_v, len(ORDER) - 0.45, " ceiling", fontsize=6, color="0.25", va="center")
ax.text(floor_v, len(ORDER) - 0.45, "floor ", fontsize=6, color="0.55", va="center",
        ha="right")
ax.set_title("Predicted H3K27ac as the ABC activity term", fontsize=8)
fig.tight_layout(rect=(0.30, 0, 1, 1))
save_fig(fig, f"{P}/figures/fig14_crispr_benchmark.png")
print(f"\nfloor {floor_v:.3f}, ceiling {ceil_v:.3f}, headroom {ceil_v - floor_v:+.3f}")
print("wrote fig14")
