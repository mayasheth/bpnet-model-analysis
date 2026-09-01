#!/usr/bin/env python
"""Figure 8: model-free ATAC-H3K27ac coupling, compared ACROSS cell types.

The previous version of this figure showed a single K562 hexbin while the surrounding text
made a claim about K562 versus GM12878, so the panel did not support the claim. This plots
each cell type side by side on identical axes and adds the summary that the claim rests on.
"""
import os, sys
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label, figsize

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
TR = f"{D}/2026_0606_GM12878_transferability"
NP = ["chr", "start", "end", "name", "score", "strand",
      "signalValue", "pValue", "qValue", "summit"]
HW = 500

PANELS = [
    ("K562", f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
     [f"{P}/data/h3k27ac_5p_plus.bw", f"{P}/data/h3k27ac_5p_minus.bw"],
     f"{D}/reference/K562_DNase_candidate_elements.narrowPeak"),
    ("GM12878", f"{TR}/data/atac.bw",
     [f"{P}/data/gm12878_h3k27ac_5p_plus.bw", f"{P}/data/gm12878_h3k27ac_5p_minus.bw"],
     f"{TR}/reference/GM12878_candidate_elements.narrowPeak"),
    ("TeloHAEC", f"{P}/data/panel/TeloHAEC_ctrl_atac.bw",
     [f"{P}/data/panel/TeloHAEC_ctrl_h3k27ac_r1_5p_plus.bw",
      f"{P}/data/panel/TeloHAEC_ctrl_h3k27ac_r1_5p_minus.bw"],
     f"{P}/reference/celltype_elements/TeloHAEC_ctrl_ATAC_candidate_elements.narrowPeak"),
]


def load_els(path):
    d = pd.read_csv(path, sep="\t", header=None, names=NP)
    mid = (d["end"] - d["start"]) // 2
    d["center"] = d["start"] + d["summit"].where(d["summit"] >= 0, mid)
    return d


def sums(paths, els):
    bws = [pyBigWig.open(p) for p in paths]
    ch = bws[0].chroms()
    out = np.full(len(els), np.nan)
    for i, r in enumerate(els.itertuples()):
        s, e = int(r.center) - HW, int(r.center) + HW
        if r.chr not in ch or s < 0 or e > ch[r.chr]:
            continue
        t = 0.0
        for bw in bws:
            t += np.nan_to_num(np.array(bw.values(r.chr, s, e), dtype=np.float64)).sum()
        out[i] = t
    for bw in bws:
        bw.close()
    return out


apply_rcparams()
fig, axes = plt.subplots(1, 4, figsize=figsize(columns=2, aspect=0.30))
rows = []
for k, (name, atac, k27, elp) in enumerate(PANELS):
    els = load_els(elp)
    a = np.log1p(sums([atac], els))
    h = np.log1p(sums(k27, els))
    ok = np.isfinite(a) & np.isfinite(h)
    a, h = a[ok], h[ok]
    top = h >= np.quantile(h, 0.8)
    r_all = pearsonr(a, h)[0]
    r_top = pearsonr(a[top], h[top])[0]
    rows.append({"cell_type": name, "n": int(ok.sum()),
                 "pearson_all": round(float(r_all), 4),
                 "pearson_top_quintile": round(float(r_top), 4)})
    ax = axes[k]
    ax.hexbin(a, h, gridsize=45, bins="log", cmap="Blues", linewidths=0)
    ax.set_xlabel("log1p ATAC")
    if k == 0:
        ax.set_ylabel("log1p H3K27ac")
    ax.set_title(f"{name}\ntop-quintile $r$ = {r_top:.3f}", fontsize=7)
    add_panel_label(ax, "abc"[k])
    print(f"{name}: n={ok.sum():,} r_all={r_all:.4f} r_top={r_top:.4f}", flush=True)

# --- panel d: the summary the claim rests on --------------------------------
ax = axes[3]
cp = pd.read_csv(f"{P}/results/atac_vs_h3k27ac_by_celltype.tsv", sep="\t")
cp = cp[cp["stratum"] == "top_quintile"].set_index("label")["pearson"]
ORDER = [("K562", "K562"), ("GM12878", "GM12878"), ("TeloHAEC_ctrl", "TeloHAEC"),
         ("TeloHAEC_IL1b", "+IL1b"), ("TeloHAEC_TNFa", "+TNFa"), ("TeloHAEC_VEGF", "−VEGF")]
vals = [cp.get(k, np.nan) for k, _ in ORDER]
# neutral tones: blue/red/purple are reserved for input modality elsewhere in the
# report, and these bars distinguish cell types, not modalities.
cols = ["#3A3A3A", "#3A3A3A", "#8C8C8C", "#BDBDBD", "#BDBDBD", "#BDBDBD"]
ax.bar(range(len(vals)), vals, color=cols, width=0.7)
for i, v in enumerate(vals):
    ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=4.5)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([lab for _, lab in ORDER], fontsize=5.5, rotation=40, ha="right")
ax.set_ylabel("Coupling, top quintile ($r$)")
ax.set_ylim(0, 0.60)
ax.set_title("Coupling falls outside K562", fontsize=7)
add_panel_label(ax, "d")

fig.tight_layout()
out = save_fig(fig, f"{P}/figures/fig8_coupling_across_celltypes.png")
print("Wrote", out, "and .pdf")
df = pd.DataFrame(rows)
df.to_csv(f"{P}/results/coupling_panel_recomputed.tsv", sep="\t", index=False)
print(df.to_string(index=False))
