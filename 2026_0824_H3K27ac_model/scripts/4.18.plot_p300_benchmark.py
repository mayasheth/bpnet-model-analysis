#!/usr/bin/env python
"""Figure: p300 as the ABC activity term, with the transfer arm.

Companion to 4.8 (the H3K27ac version) and deliberately built to the same geometry, same
floor and same observed-H3K27ac reference line, so the two figures can be read side by side.
Both runs share those two anchor arms and they agree to 1e-4 across the runs, which is what
licenses the comparison.

The panel is a forest plot for the same reason 4.8 is: the per-predictor intervals are ~+/-0.05
wide while the differences are 0.01-0.10, so a bar chart would imply resolution that the
unpaired intervals do not have. The paired-bootstrap deltas that DO resolve them are annotated
directly on the figure rather than left in a table, because the unpaired intervals are the
thing a reader's eye is drawn to and they are the misleading part.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, annotate_n_fig

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
SRC = ("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/"
       "workflow/results/2026_0906_p300_all/performance_summary.txt")

d = pd.read_csv(SRC, sep="\t")
d["arm"] = d["pred_uid"].str.replace(r"\.ABC\.Score$", "", regex=True)

MERGED = os.path.join(os.path.dirname(SRC), "expt_pred_merged_annot.txt.gz")
_m = pd.read_csv(MERGED, sep="\t", usecols=["pred_uid", "Regulated"], low_memory=False)
_one = _m[_m["pred_uid"] == _m["pred_uid"].iloc[0]]
N_LABEL = (f"n = {len(_one):,} element-gene pairs, "
           f"{int(_one['Regulated'].astype(str).str.upper().eq('TRUE').sum()):,} regulated; "
           f"identical pair set for every arm")

FLOOR, K27 = "K562_ATAC_only", "K562_ATAC_H3K27ac_element"
# Paired-bootstrap deltas against the floor, from 4.17. Annotated on the figure because the
# unpaired intervals drawn here cannot resolve any of them.
PAIRED = {
    "p300obs_k562": "+0.099 [+0.079, +0.120]",
    "p300pred_k562_multimodal": "+0.055 [+0.034, +0.074]",
    "p300only_k562_multimodal": "+0.045 [+0.012, +0.077]",
    "p300pred_gm12878_multimodal": "+0.009 [-0.004, +0.022]  n.s.",
    "p300only_gm12878_multimodal": "+0.005 [-0.019, +0.028]  n.s.",
    "p300pred_k562_atac": "+0.004 [-0.012, +0.021]  n.s.",
}
ORDER = [
    ("p300obs_k562", "ATAC x observed p300", "#54278f"),
    (K27, "ATAC x observed H3K27ac", "#404040"),
    ("p300pred_k562_multimodal", "ATAC x predicted p300 (K562 seq+ATAC)", "#807dba"),
    ("p300only_k562_multimodal", "Predicted p300 alone (K562 seq+ATAC)", "#9e9ac8"),
    ("p300pred_gm12878_multimodal", "ATAC x predicted p300 (GM12878 -> K562)", "#08519c"),
    ("p300only_gm12878_multimodal", "Predicted p300 alone (GM12878 -> K562)", "#6baed6"),
    ("p300pred_k562_atac", "ATAC x predicted p300 (K562 ATAC-only)", "#bcbddc"),
    (FLOOR, "ATAC only", "#bdbdbd"),
    ("p300only_k562_atac", "Predicted p300 alone (K562 ATAC-only)", "#d9d9d9"),
    ("baseline.distToTSS", "Distance to TSS", "#c5cad7"),
]
got = {r["arm"]: r for _, r in d.iterrows()}
missing = [a for a, _l, _c in ORDER if a not in got]
assert not missing, f"missing arms in summary: {missing}"

floor_v, k27_v = float(got[FLOOR]["AUPRC"]), float(got[K27]["AUPRC"])

apply_rcparams()
fig, ax = plt.subplots(figsize=(7.2, 3.9))
ax.axvline(floor_v, color="0.55", lw=0.9, ls="--", zorder=1)
ax.axvline(k27_v, color="0.25", lw=0.9, ls=(0, (4, 2)), zorder=1)

XL, XR = 0.235, 0.66
for i, (arm, lab, colour) in enumerate(ORDER):
    y = len(ORDER) - 1 - i
    r = got[arm]
    v, lo, hi = float(r["AUPRC"]), float(r["AUPRC_lowerCi"]), float(r["AUPRC_upperCi"])
    ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", ms=4.5, color=colour,
                lw=1.3, capsize=2.5, zorder=3)
    ax.text(XL, y, lab, ha="right", va="center", fontsize=6.5)
    if arm in PAIRED:
        ax.text(XR, y, PAIRED[arm], ha="right", va="center", fontsize=5.5,
                color="#333333", family="monospace")
    print(f"{lab:<44} {v:.3f} [{lo:.3f}, {hi:.3f}]")

ax.set_yticks([])
ax.set_ylim(-0.8, len(ORDER) - 0.2)
ax.set_xlim(XL, XR)
ax.set_xlabel("AUPRC against the CRISPR benchmark (unpaired 95% CI)")
ax.text(k27_v, len(ORDER) - 0.45, " observed H3K27ac", fontsize=6, color="0.25", va="center")
ax.text(floor_v, len(ORDER) - 0.45, "floor ", fontsize=6, color="0.55", va="center", ha="right")
ax.text(XR, len(ORDER) - 0.45, "paired delta vs floor", fontsize=5.5, color="#333333",
        ha="right", va="center", style="italic")
ax.set_title("p300 as the ABC activity term: the gain is real in K562 and absent on transfer",
             fontsize=8)
fig.tight_layout(rect=(0.30, 0, 1, 1))
annotate_n_fig(fig, N_LABEL)
save_fig(fig, f"{P}/figures/fig16_p300_benchmark")
print(f"\nfloor {floor_v:.3f}, observed H3K27ac {k27_v:.3f}, "
      f"observed p300 {float(got['p300obs_k562']['AUPRC']):.3f}")
print("wrote fig16")
