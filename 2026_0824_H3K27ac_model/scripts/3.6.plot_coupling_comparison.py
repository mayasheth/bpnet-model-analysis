#!/usr/bin/env python
"""Figure 8: model-free ATAC-H3K27ac coupling, compared ACROSS cell types.

The previous version of this figure showed a single K562 hexbin while the surrounding text
made a claim about K562 versus GM12878, so the panel did not support the claim. This plots
each cell type side by side on identical axes and adds the summary that the claim rests on.
"""
import argparse
import hashlib
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


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--recompute", action="store_true",
                help="Ignore the cache and re-read every bigwig.")
args = ap.parse_args()

CACHE = f"{P}/results/.cache_coupling_panel"
os.makedirs(CACHE, exist_ok=True)


def cache_key(name, atac, k27, elp):
    """Key on the inputs AND their mtimes, so a rebuilt track invalidates the cache
    rather than silently serving stale vectors."""
    parts = [name, atac, *k27, elp, str(HW)]
    for f in [atac, *k27, elp]:
        parts.append(str(int(os.path.getmtime(f))) if os.path.exists(f) else "missing")
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]


def load_or_compute(name, atac, k27, elp):
    key = cache_key(name, atac, k27, elp)
    path = f"{CACHE}/{name}_{key}.npz"
    if os.path.exists(path) and not args.recompute:
        d = np.load(path)
        print(f"{name}: cache hit ({os.path.basename(path)})", flush=True)
        return d["a"], d["h"]
    els = load_els(elp)
    a = np.log1p(sums([atac], els))
    h = np.log1p(sums(k27, els))
    ok = np.isfinite(a) & np.isfinite(h)
    a, h = a[ok], h[ok]
    np.savez_compressed(path, a=a, h=h)
    print(f"{name}: computed and cached -> {os.path.basename(path)}", flush=True)
    return a, h


apply_rcparams()
# row 1: three hexbins (2 grid cols each); row 2: two bar panels (3 grid cols each)
fig = plt.figure(figsize=figsize(columns=2, aspect=0.62))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.85], hspace=0.75, wspace=0.55)
axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]),
        fig.add_subplot(gs[0, 4:6]),
        fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
rows = []
for k, (name, atac, k27, elp) in enumerate(PANELS):
    a, h = load_or_compute(name, atac, k27, elp)
    ok = np.ones(len(a), bool)
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

# --- panels d, e: coupling across cell types, both strata -------------------
# Both strata are shown because the ordering INVERTS between them: on all elements
# TeloHAEC is the highest, on the top quintile it is the lowest. The all-elements
# number is dominated by the dead-versus-active contrast, so the top quintile is the
# informative one — the same stratification caveat that applies throughout the report.
cp = pd.read_csv(f"{P}/results/atac_vs_h3k27ac_by_celltype.tsv", sep="\t")
ORDER = [("K562", "K562"), ("GM12878", "GM12878"), ("TeloHAEC_ctrl", "TeloHAEC"),
         ("TeloHAEC_IL1b", "+IL1b"), ("TeloHAEC_TNFa", "+TNFa"), ("TeloHAEC_VEGF", "−VEGF")]
# neutral tones: blue/red/purple are reserved for input modality elsewhere in the
# report, and these bars distinguish cell types, not modalities.
cols = ["#3A3A3A", "#3A3A3A", "#8C8C8C", "#BDBDBD", "#BDBDBD", "#BDBDBD"]
for ax, stratum, lab, ttl in [
        (axes[3], "all", "All elements", "All elements: TeloHAEC highest"),
        (axes[4], "top_quintile", "Top quintile", "Top quintile: ordering inverts")]:
    ser = cp[cp["stratum"] == stratum].set_index("label")["pearson"]
    vals = [ser.get(k, np.nan) for k, _ in ORDER]
    ax.bar(range(len(vals)), vals, color=cols, width=0.7)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=5)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([l for _, l in ORDER], fontsize=6, rotation=40, ha="right")
    ax.set_ylabel(f"Coupling, {lab.lower()} ($r$)")
    ax.set_ylim(0, 0.85)
    ax.set_title(ttl, fontsize=7)
add_panel_label(axes[3], "d")
add_panel_label(axes[4], "e")

out = save_fig(fig, f"{P}/figures/fig8_coupling_across_celltypes.png")
print("Wrote", out, "and .pdf")
df = pd.DataFrame(rows)
df.to_csv(f"{P}/results/coupling_panel_recomputed.tsv", sep="\t", index=False)
print(df.to_string(index=False))
