#!/usr/bin/env python
"""Figure: p300 as the ABC activity term, with the transfer arm.

Companion to 4.8 (the H3K27ac version), same floor and same observed-H3K27ac reference line
so the two can be read side by side. The two runs share both anchor arms and agree on them to
1e-4, which is what licenses the comparison.

Forest plot rather than bars for the same reason as 4.8: the per-predictor intervals are
~+/-0.05 wide while the differences are 0.004-0.099, so bars would imply resolution the
unpaired intervals do not have. The paired-bootstrap deltas that DO resolve them are printed
alongside, because the overlapping bars are what a reader's eye is drawn to and they are the
misleading layer.

LAYOUT. Row labels and the delta column live OUTSIDE the axes, positioned with
`ax.get_yaxis_transform()` (x in axes fraction, y in data coordinates) and `clip_on=False`,
with margins reserved by `subplots_adjust`. An earlier version placed the deltas inside the
axes with ha="right", which drew them leftward straight across the error bars, and put three
header labels at colliding data coordinates. Text that belongs beside a plot should not be
positioned in data space.
"""
import os
import sys
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
# Paired-bootstrap deltas against the floor, from 4.17.
PAIRED = {
    "p300obs_k562":                ("+0.099", "[+0.079, +0.120]", ""),
    "p300pred_k562_multimodal":    ("+0.055", "[+0.034, +0.074]", ""),
    "p300only_k562_multimodal":    ("+0.045", "[+0.012, +0.077]", ""),
    "p300pred_gm12878_multimodal": ("+0.009", "[-0.004, +0.022]", "n.s."),
    "p300only_gm12878_multimodal": ("+0.005", "[-0.019, +0.028]", "n.s."),
    "p300pred_k562_atac":          ("+0.004", "[-0.012, +0.021]", "n.s."),
}
ORDER = [
    ("p300obs_k562", "ATAC × observed p300", "#54278f"),
    (K27, "ATAC × observed H3K27ac", "#404040"),
    ("p300pred_k562_multimodal", "ATAC × predicted p300 (K562 seq+ATAC)", "#807dba"),
    ("p300only_k562_multimodal", "Predicted p300 alone (K562 seq+ATAC)", "#9e9ac8"),
    ("p300pred_gm12878_multimodal", "ATAC × predicted p300 (GM12878→K562)", "#08519c"),
    ("p300only_gm12878_multimodal", "Predicted p300 alone (GM12878→K562)", "#6baed6"),
    ("p300pred_k562_atac", "ATAC × predicted p300 (K562 ATAC-only)", "#bcbddc"),
    (FLOOR, "ATAC only", "#bdbdbd"),
    ("p300only_k562_atac", "Predicted p300 alone (K562 ATAC-only)", "#d9d9d9"),
    ("baseline.distToTSS", "Distance to TSS", "#c5cad7"),
]
got = {r["arm"]: r for _, r in d.iterrows()}
missing = [a for a, _l, _c in ORDER if a not in got]
assert not missing, f"missing arms in summary: {missing}"

floor_v, k27_v = float(got[FLOOR]["AUPRC"]), float(got[K27]["AUPRC"])
N = len(ORDER)

apply_rcparams()
fig, ax = plt.subplots(figsize=(7.8, 3.6))
# Margins reserved for the label column (left) and the delta column (right).
fig.subplots_adjust(left=0.315, right=0.735, top=0.86, bottom=0.16)
trans = ax.get_yaxis_transform()          # x: axes fraction, y: data

ax.axvline(floor_v, color="0.55", lw=0.9, ls="--", zorder=1)
ax.axvline(k27_v, color="0.30", lw=0.9, ls=(0, (4, 2)), zorder=1)

for i, (arm, lab, colour) in enumerate(ORDER):
    y = N - 1 - i
    r = got[arm]
    v, lo, hi = float(r["AUPRC"]), float(r["AUPRC_lowerCi"]), float(r["AUPRC_upperCi"])
    ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", ms=4.5, color=colour,
                lw=1.3, capsize=2.5, zorder=3)
    ax.text(-0.02, y, lab, ha="right", va="center", fontsize=6.5,
            transform=trans, clip_on=False)
    if arm in PAIRED:
        dv, ci, ns = PAIRED[arm]
        weight = "bold" if not ns else "normal"
        col = "#111111" if not ns else "#888888"
        ax.text(1.03, y, dv, ha="left", va="center", fontsize=6, color=col,
                weight=weight, family="monospace", transform=trans, clip_on=False)
        ax.text(1.13, y, ci, ha="left", va="center", fontsize=5.6, color=col,
                family="monospace", transform=trans, clip_on=False)
        if ns:
            ax.text(1.42, y, ns, ha="left", va="center", fontsize=5.6, color="#888888",
                    style="italic", transform=trans, clip_on=False)
    print(f"{lab:<44} {v:.3f} [{lo:.3f}, {hi:.3f}]")

# Column header for the delta block, outside the axes so it cannot collide with the data.
ax.text(1.03, N - 0.35, "paired Δ vs floor (95% CI)", ha="left", va="center", fontsize=5.8,
        color="#333333", style="italic", transform=trans, clip_on=False)

ax.set_yticks([])
ax.set_ylim(-0.6, N - 0.02)
ax.set_xlim(0.38, 0.63)
ax.set_xlabel("AUPRC against the CRISPR benchmark (unpaired 95% CI)")
# Reference-line labels ride above the top row, one extending left of its line and the other
# right, so they clear each other and the tick labels. Below the axis they collided with the
# 0.45 and 0.50 ticks.
ax.text(floor_v, N - 0.42, "floor ", ha="right", va="center", fontsize=5.8, color="0.45")
ax.text(k27_v, N - 0.42, " observed H3K27ac", ha="left", va="center", fontsize=5.8,
        color="0.30")
ax.set_title("p300 as the ABC activity term: real in K562, absent on transfer", fontsize=8)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
annotate_n_fig(fig, N_LABEL)
save_fig(fig, f"{P}/figures/fig16_p300_benchmark")
print(f"\nfloor {floor_v:.3f}, observed H3K27ac {k27_v:.3f}, "
      f"observed p300 {float(got['p300obs_k562']['AUPRC']):.3f}")
print("wrote fig16")
