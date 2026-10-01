#!/usr/bin/env python3
"""Report 6 fig39: benchmark results for F-041 (RNA TF lists), F-042 (fragment channels) and
F-043 (p300 x H3K27ac, transferred H3K27ac correction), from the paired-bootstrap tables.

    a  seed-averaged AUPRC differences, ABC and ENCODE-rE2G (95% bootstrap CI)
    b  observed activity terms in ABC with and without the ATAC term, and in rE2G
    c  deployable H3K27ac and combination arms against their baselines

Usage: 5.4.plot_benchmark_round2.py [--outdir figures]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
R = f"{P}/results"
ap = argparse.ArgumentParser()
ap.add_argument("--outdir", default=f"{P}/figures")
a = ap.parse_args()
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
})


def load(f, abc):
    t = pd.read_csv(f"{R}/{f}", sep="\t")
    if not abc:
        t = t.rename(columns={"model_a": "arm_a", "model_b": "arm_b"})
    return t.set_index(["arm_a", "arm_b"])


rna, frag = load("2026_1001_rna_auprc_pairs.tsv", True), load("2026_1001_frag_auprc_pairs.tsv", True)
re2g = load("2026_1001_re2g_auprc_pairs.tsv", False)
k27, obs = load("2026_1001_k27frag_auprc_pairs.tsv", True), load("2026_1001_k27frag_obsonly_auprc_pairs.tsv", True)

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.9})

# a: seed-averaged deltas
rows = [("RNA + ATAC TF list - MSE-only", rna, ("rnaatac", "mse"), ("rnaatac", "mse")),
        ("RNA + ATAC TF list - ATAC TF list", rna, ("rnaatac", "ctcf"), ("rnaatac", "ctcf")),
        ("RNA TF list - MSE-only", rna, ("rnafam", "mse"), ("rnafam", "mse")),
        ("RNA TF list - ATAC TF list", rna, ("rnafam", "ctcf"), ("rnafam", "ctcf")),
        ("fragment channels - flat 5' input", frag, ("frag", "flat"), ("frag", "flat")),
        ("flat 5' input - MSE-only", frag, ("flat", "mse"), ("flat", "mse"))]
y = np.arange(len(rows))[::-1]
for off, col, lab in ((0.15, "#185FA5", "ABC"), (-0.15, "#D85A30", "ENCODE-rE2G")):
    for yi, (name, t, kabc, kre) in zip(y, rows):
        r = t.loc[kabc] if lab == "ABC" else re2g.loc[kre]
        ax[0].errorbar(r["delta"], yi + off, xerr=[[r["delta"] - r["ci_lo"]], [r["ci_hi"] - r["delta"]]],
                       fmt="o", ms=4, color=col, lw=1, label=lab if yi == y[0] else None)
ax[0].axvline(0, color="0.6", lw=0.7)
ax[0].set_yticks(y)
ax[0].set_yticklabels([r[0] for r in rows], fontsize=8)
ax[0].set_xlabel("AUPRC difference, mean of 3 runs (95% CI)")
ax[0].set_title("a  RNA TF lists and fragment channels")
ax[0].legend(frameon=False, fontsize=7, loc="lower left")

# b: observed activity terms
lab = ["H3K27ac", "p300", "p300 x\nH3K27ac"]
abc_with = [k27.loc[("k27frag_obs_p300.ABC.Score", "k27frag_obs_k27.ABC.Score")]["auprc_b"],
            k27.loc[("k27frag_obs_p300.ABC.Score", "k27frag_obs_k27.ABC.Score")]["auprc_a"],
            k27.loc[("k27frag_obs_p300k27.ABC.Score", "k27frag_obs_k27.ABC.Score")]["auprc_a"]]
abc_wo = [obs.loc[("k27frag_obsonly_k27.ABC.Score", "k27frag_obs_k27.ABC.Score")]["auprc_a"],
          obs.loc[("k27frag_obsonly_p300.ABC.Score", "k27frag_obs_p300.ABC.Score")]["auprc_a"],
          obs.loc[("k27frag_obsonly_p300k27.ABC.Score", "k27frag_obsonly_p300.ABC.Score")]["auprc_a"]]
re = [re2g.loc[("ATAC_p300obspaint_powerlaw_all", "ATAC_k27obspaint_powerlaw_all")]["auprc_b"],
      re2g.loc[("ATAC_p300obspaint_powerlaw_all", "ATAC_k27obspaint_powerlaw_all")]["auprc_a"],
      re2g.loc[("ATAC_p300k27obs_powerlaw_all", "ATAC_k27obspaint_powerlaw_all")]["auprc_a"]]
x = np.arange(3)
for dx, v, c, n in ((-0.27, abc_with, "#85B7EB", "ABC, ATAC x mark"), (0, abc_wo, "#185FA5", "ABC, mark alone"),
                    (0.27, re, "#D85A30", "ENCODE-rE2G (learned weights)")):
    b = ax[1].bar(x + dx, v, 0.25, color=c, label=n)
    for xi, vi in zip(x + dx, v):
        ax[1].text(xi, vi + 0.003, f"{vi:.3f}", ha="center", fontsize=7)
floor = k27.loc[("p300pred_multicell_holdK562.ABC.Score", "K562_ATAC_only.ABC.Score")]["auprc_b"]
ax[1].axhline(floor, color="0.5", lw=0.8, ls="--")
ax[1].text(-0.45, floor - 0.012, "ABC ATAC-only floor", fontsize=7, ha="left", color="0.4")
ax[1].set_xticks(x)
ax[1].set_xticklabels(lab)
ax[1].set_ylim(0.44, 0.73)
ax[1].set_ylabel("AUPRC, K562 CRISPR")
ax[1].set_title("b  observed activity terms")
ax[1].legend(frameon=False, fontsize=7, loc="upper left")

# c: deployable H3K27ac and combinations
rows = [("H3K27ac recalibrated\n- transferred", ("k27frag_tx_base.ABC.Score", "k27frag_tx_repaint.ABC.Score"),
         ("ATAC_k27txbase_powerlaw_all", "ATAC_k27txrepaint_powerlaw_all")),
        ("+ fragment and 1 kb ATAC\nfeatures - recalibrated", ("k27frag_tx_breadth_frag.ABC.Score", "k27frag_tx_base.ABC.Score"),
         ("ATAC_k27txbreadthfrag_powerlaw_all", "ATAC_k27txbase_powerlaw_all")),
        ("pred. p300 x corrected\nH3K27ac - pred. p300", ("k27frag_pred_p300k27.ABC.Score", "p300pred_multicell_holdK562.ABC.Score"),
         ("ATAC_p300k27pred_powerlaw_all", "mse"))]
y = np.arange(len(rows))[::-1]
for off, col, key in ((0.15, "#185FA5", 1), (-0.15, "#D85A30", 2)):
    for yi, r_ in zip(y, rows):
        r = k27.loc[r_[1]] if key == 1 else re2g.loc[r_[2]]
        ax[2].errorbar(r["delta"], yi + off, xerr=[[r["delta"] - r["ci_lo"]], [r["ci_hi"] - r["delta"]]],
                       fmt="o", ms=4, color=col, lw=1, label=("ABC" if key == 1 else "ENCODE-rE2G") if yi == y[0] else None)
ax[2].axvline(0, color="0.6", lw=0.7)
ax[2].set_yticks(y)
ax[2].set_yticklabels([r[0] for r in rows], fontsize=8)
ax[2].set_xlabel("AUPRC difference (95% CI)")
ax[2].set_title("c  deployable H3K27ac and p300 x H3K27ac")
ax[2].legend(frameon=False, fontsize=7, loc="lower right")
os.makedirs(a.outdir, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{a.outdir}/fig39_benchmark_round2.{ext}")
print(f"wrote {a.outdir}/fig39_benchmark_round2.png")
