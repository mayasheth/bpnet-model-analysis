#!/usr/bin/env python3
"""Explainer figure for partial correlation given ATAC (the 4.35 metric), on real data.

Non-promoter ABC candidate regions only (132,186), values as ABC saw them, converted to
percentile ranks because ABC keeps only ranks. 3,000 regions sampled for display; every rho
printed on the figure is computed on ALL regions.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata

ABC = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OUT = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/figures"
KEY = ["chr", "start", "end"]
OBS = f"{ABC}/2026_0905_p300_activity/p300obs_k562"
MODELS = [
    ("K562 in-cell",  f"{ABC}/2026_0905_p300_activity/p300pred_k562_multimodal", "#2a78d6"),
    ("multi-cell",    f"{ABC}/2026_0922_p300_multicell/p300pred_multicell_holdK562", "#eb6834"),
    ("GM12878 only",  f"{ABC}/2026_0922_p300_multicell/p300pred_gm12878_corrected_tgtnorm",
     "#1baf7a"),
]
H3K27AC = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element"
INK, INK2, GRAY = "#0b0b0b", "#52514e", "#9a9994"

plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "figure.dpi": 170, "savefig.dpi": 170, "savefig.bbox": "tight"})


def load(d):
    x = pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t",
                    usecols=KEY + ["class", "H3K27ac.RPM", "ATAC.RPM"])
    return x[~x["class"].isin(["tss", "promoter"])].set_index(KEY).sort_index()


def pct(v):
    return 100 * (rankdata(v) - 0.5) / len(v)


def resid(y, z):
    Z = np.column_stack([np.ones_like(z), z])
    return y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]


obs = load(OBS)
atac = pct(obs["ATAC.RPM"].values)
real = pct(obs["H3K27ac.RPM"].values)
preds = {}
for name, d, col in MODELS + [("observed H3K27ac", H3K27AC, "#4a3aa7")]:
    p = load(d).reindex(obs.index)
    assert p["H3K27ac.RPM"].notna().all()
    preds[name] = pct(p["H3K27ac.RPM"].values)

rng = np.random.default_rng(0)
show = rng.choice(len(obs), 3000, replace=False)
r = lambda a, b: float(np.corrcoef(a, b)[0, 1])      # Pearson on ranks = Spearman
rr = resid(real, atac)

fig = plt.figure(figsize=(13.6, 8.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.55, wspace=0.62)


def dots(ax, x, y, col, title, rho, xl, yl, lim=(0, 100)):
    ax.scatter(x[show], y[show], s=4, color=col, alpha=0.28, linewidths=0)
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, pad=6)
    ax.text(0.04, 0.93, f"ρ = {rho:.2f}", transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK, va="top")
    ax.set_xlabel(xl); ax.set_ylabel(yl)


# ---- row 1: the problem, on the in-cell model
inc = preds["K562 in-cell"]
ax = fig.add_subplot(gs[0, 0])
dots(ax, atac, real, GRAY, "1  ATAC predicts p300", r(atac, real), "ATAC (rank)", "real p300 (rank)")
ax = fig.add_subplot(gs[0, 1])
dots(ax, inc, real, "#2a78d6", "2  Model predicts p300", r(inc, real),
     "model prediction (rank)", "real p300 (rank)")
ax = fig.add_subplot(gs[0, 2])
dots(ax, atac, inc, "#2a78d6", "3  …but the model mostly\n    follows ATAC", r(atac, inc),
     "ATAC (rank)", "model prediction (rank)")
ax = fig.add_subplot(gs[0, 3]); ax.axis("off")
ax.text(0.0, 0.78, "So is panel 2 just panel 1\nin disguise?", fontsize=12, color=INK,
        fontweight="bold", va="top")
ax.text(0.0, 0.42, "Partial correlation answers it:\nsubtract what ATAC explains\n"
        "from BOTH, then correlate\nwhat is left over.", fontsize=11, color=INK2, va="top")

# ---- row 2: leftovers, three models + summary
for i, (name, _, col) in enumerate(MODELS):
    ax = fig.add_subplot(gs[1, i])
    rp = resid(preds[name], atac)
    ax.scatter(rp[show], rr[show], s=4, color=col, alpha=0.28, linewidths=0)
    ax.axhline(0, color=GRAY, lw=0.8); ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlim(-60, 60); ax.set_ylim(-60, 60)
    ax.set_title(f"{'456'[i]}  {name}: leftovers", loc="left", fontsize=10.5, color=INK, pad=6)
    ax.text(0.04, 0.93, f"ρ = {r(rp, rr):.2f}", transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK, va="top")
    ax.set_xlabel("prediction, ATAC removed")
    if i == 0:
        ax.set_ylabel("real p300, ATAC removed")

ax = fig.add_subplot(gs[1, 3])
order = ["observed H3K27ac", "K562 in-cell", "multi-cell", "GM12878 only"]
cols = {"observed H3K27ac": "#4a3aa7", **{n: c for n, _, c in MODELS}}
vals = [r(resid(preds[n], atac), rr) for n in order]
ax.barh(range(4), vals, color=[cols[n] for n in order], height=0.55)
for k, v in enumerate(vals):
    ax.text(v + 0.012, k, f"{v:.2f}", va="center", fontsize=10, color=INK)
ax.set_yticks(range(4)); ax.set_yticklabels(order, fontsize=9.5)
ax.invert_yaxis(); ax.set_xlim(0, 0.65)
ax.set_title("7  Knowledge beyond ATAC", loc="left", fontsize=10.5, color=INK, pad=6)
ax.set_xlabel("partial ρ")
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)

fig.text(0.01, 0.985, "Partial correlation: what does a prediction know about p300 that ATAC "
         "doesn't?", fontsize=13.5, fontweight="bold", color=INK, va="top")
fig.text(0.01, 0.955, "Non-promoter ABC candidate regions (n = 132,186), K562. 3,000 shown; "
         "every ρ uses all. Ranks, because ABC only uses ranks.",
         fontsize=9.5, color=INK2, va="top")
out = f"{OUT}/fig25_partial_corr_explainer.png"
fig.savefig(out); print("wrote", out)
