#!/usr/bin/env python3
"""Report 7 figures: element-level TF evidence from Xu K562 multiome for scE2G (F-047 to F-053).

    fig40  headline: a  p300 information beyond ATAC and motif presence, per evidence type
                     b  scE2G CRISPR AUPRC change vs baseline, every variant
    fig41  footprint ChIP validation within accessibility deciles (homemade, TOBIAS, PRINT)
    fig42  footprint p300 gate per motif, TOBIAS and PRINT, two sets of controls
    fig43  footprints inside the ABC activity term: a  AUPRC change by weight  b  activity coefficient
    fig44  peak-TF Kendall (feature A): a  each feature beyond ATAC  b  kA_max beyond added controls

Inputs are the A2 summary tables (`2026_1002_sce2g_tf_activity/results/`). Model coefficients come
from the scE2G runs on $SCRATCH and are copied once into results/2026_1004_sce2g_coefficients.tsv.

Usage: 5.5.plot_tf_evidence_report_figures.py [--outdir figures]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity/results"
RUNS = "/scratch/users/shethm/sce2g_runs"
ap = argparse.ArgumentParser()
ap.add_argument("--outdir", default=f"{P}/figures")
a = ap.parse_args()
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "text.color": "black", "axes.labelcolor": "black", "axes.edgecolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
})
BLUE, LBLUE, RED, LRED, YEL, GREY, TEAL = "#00488d", "#5496ce", "#9b241c", "#dc6464", "#ca9b23", "#6e788d", "#0096a0"


def tsv(f):
    return pd.read_csv(f"{A2}/{f}", sep="\t")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {name}")


def ci_dot(ax, x, y, lo, hi, color, label=None, marker="o", ms=5):
    ax.errorbar(x, y, xerr=[[x - lo], [hi - x]], fmt=marker, ms=ms, color=color, lw=1.1,
                capsize=0, label=label)


gate = {(s, z): tsv(f"2026_1003_footprint_p300_gate{'.print' if s == 'print' else ''}"
                    f"{'.z_encode_atac' if z else ''}.beyond_motif.tsv").set_index("motif")
        for s in ("tobias", "print") for z in (False, True)}
kend = {z: tsv(f"2026_1004_featA_peak_tf_kendall{'.z_encode_atac' if z else ''}.beyond_motif.tsv")
        for z in (False, True)}
kend1 = {z: tsv(f"2026_1004_featA_peak_tf_kendall{'.z_encode_atac' if z else ''}.all_elements.tsv").set_index("feature")
         for z in (False, True)}
bt = tsv("2026_1003_featB_tobias_paired_bootstrap.tsv").set_index("comparison")
bp = tsv("2026_1004_featB_print_paired_bootstrap.tsv").set_index("comparison")
ba = tsv("2026_1004b_abcfp_paired_bootstrap.tsv").set_index("comparison")
B = "multiome_powerlaw_n6"
T = "multiome_n6_arcfp"

# ---------------------------------------------------------------- fig40 headline
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.2), gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.95})
rows = [("TOBIAS erythroid footprint", lambda z: gate[("tobias", z)].loc["erythroid"]),
        ("PRINT erythroid footprint", lambda z: gate[("print", z)].loc["erythroid"]),
        ("peak-TF Kendall (feature A)", lambda z: kend[z].iloc[-1]),
        ("TOBIAS, HNF4A (not expressed)", lambda z: gate[("tobias", z)].loc["HNF4A"]),
        ("PRINT, HNF4A (not expressed)", lambda z: gate[("print", z)].loc["HNF4A"])]
y = np.arange(len(rows))[::-1]
for off, z, col, lab in ((0.14, False, BLUE, "controls: Xu ATAC"), (-0.14, True, LBLUE, "controls: Xu + ENCODE ATAC")):
    for yi, (_, f) in zip(y, rows):
        r = f(z)
        est = r["partial_rho_fp_max"] if "partial_rho_fp_max" in r else r["partial_rho_kA_max"]
        ci_dot(ax[0], est, yi + off, r["ci_lo"], r["ci_hi"], col, lab if yi == y[0] else None)
ax[0].axvline(0, color="black", lw=0.7)
ax[0].set_yticks(y)
ax[0].set_yticklabels([r[0] for r in rows])
ax[0].set_xlabel("Partial Spearman with K562 p300 (95% CI)")
ax[0].set_title("a  p300 information beyond ATAC and motif presence", loc="left")
ax[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2)

rows = [("+ erythroid motif count", bt.loc[f"multiome_n6_ery_motif - {B}"], BLUE),
        ("+ motif count + TOBIAS footprint", bt.loc[f"multiome_n6_ery_fp - {B}"], BLUE),
        ("+ motif count + PRINT footprint", bp.loc[f"multiome_n6_ery_fp - {B}"], BLUE),
        ("activity x TOBIAS footprint, weight 1", ba.loc[f"{T}_tobias_l1 - {B}"], RED),
        ("activity x PRINT footprint, weight 1", ba.loc[f"{T}_print_l1 - {B}"], RED),
        ("activity x TOBIAS footprint, weight 4", ba.loc[f"{T}_tobias_l4 - {B}"], RED),
        ("activity x PRINT footprint, weight 4", ba.loc[f"{T}_print_l4 - {B}"], RED)]
y = np.arange(len(rows))[::-1]
for yi, (_, r, col) in zip(y, rows):
    ci_dot(ax[1], r["estimate"], yi, r["ci_lo"], r["ci_hi"], col)
ax[1].axvline(0, color="black", lw=0.7)
ax[1].set_yticks(y)
ax[1].set_yticklabels([r[0] for r in rows])
ax[1].set_xlabel(f"CRISPR AUPRC minus baseline ({bt.loc[B, 'estimate']:.3f}), 95% CI")
ax[1].set_title("b  scE2G on the K562 CRISPR benchmark", loc="left")
ax[1].text(0, -0.17, "blue: extra feature    red: inside the ABC activity term", transform=ax[1].transAxes,
           ha="left", va="top", fontsize=7)
ax[1].xaxis.set_major_locator(MaxNLocator(6))
save(fig, "fig40_tf_evidence_headline")

# ---------------------------------------------------------------- fig41 ChIP validation
home = tsv("2026_1002_footprint_validation.chip.tsv")
tob, prn = tsv("2026_1002_tobias_validation.chip.tsv"), tsv("2026_1002_print_validation.chip.tsv")
mot = ["CTCF", "GATA1", "GATA1::TAL1"]
fig, ax = plt.subplots(figsize=(5.6, 3.4))
x = np.arange(len(mot))
w = 0.26
ax.bar(x - w, home["auroc_fs_within_ins_decile"].values, w, color=GREY, label="homemade (8-mer bias)")
for k, (t, col, lab) in enumerate(((tob, BLUE, "TOBIAS"), (prn, RED, "PRINT"))):
    v = t["auroc_score_within_ins_decile"].values
    ax.bar(x + k * w, v, w, color=col, label=lab,
           yerr=[v - t["ci_lo"].values, t["ci_hi"].values - v], error_kw={"lw": 0.9, "ecolor": "black"})
ax.axhline(0.5, color="black", lw=0.7, ls="--")
ax.set_xticks(x)
ax.set_xticklabels(mot)
ax.set_ylim(0.3, 0.9)
ax.set_ylabel("AUROC, ChIP-bound vs unbound sites,\nwithin accessibility deciles")
ax.legend(frameon=False, fontsize=7, loc="upper right")
save(fig, "fig41_footprint_chip_validation")

# ---------------------------------------------------------------- fig42 gate per motif
order = ["erythroid", "GATA1::TAL1", "MAF::NFE2", "NFE2", "TAL1::TCF3", "GATA1", "HNF4A", "CTCF"]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True, gridspec_kw={"wspace": 0.12})
y = np.arange(len(order))[::-1]
for j, z in enumerate((False, True)):
    for off, s, col in ((0.14, "tobias", BLUE), (-0.14, "print", RED)):
        g = gate[(s, z)].loc[order]
        for yi, (_, r) in zip(y, g.iterrows()):
            ci_dot(ax[j], r["partial_rho_fp_max"], yi + off, r["ci_lo"], r["ci_hi"], col,
                   s.upper() if yi == y[0] else None, ms=4)
    ax[j].axvline(0, color="black", lw=0.7)
    ax[j].xaxis.set_major_locator(MaxNLocator(5))
    ax[j].set_xlabel("Partial Spearman with K562 p300 (95% CI)")
    ax[j].set_title(("a  controls: Xu ATAC, composition, motif count and strength",
                     "b  same + ENCODE K562 ATAC")[j], loc="left", fontsize=9)
ax[0].set_yticks(y)
ax[0].set_yticklabels([f"{m} ({int(gate[('tobias', False)].loc[m, 'n_elements']):,})" for m in order])
ax[0].legend(frameon=False, fontsize=7, loc="lower right")
save(fig, "fig42_footprint_p300_gate")

# ---------------------------------------------------------------- fig43 footprints in ABC activity
cf = f"{A2}/2026_1004_sce2g_coefficients.tsv"
if not os.path.exists(cf):
    rec = []
    for run, models in (("2026_1004b_k562_xu_abcfp", [B] + [f"{T}_{s}_l{l}" for s in ("tobias", "print") for l in (1, 4)]),
                        ("2026_1003_k562_xu_featB_tobias", ["multiome_n6_ery_motif", "multiome_n6_ery_fp"]),
                        ("2026_1004_k562_xu_featB_print", ["multiome_n6_ery_fp"])):
        for m in models:
            t = pd.read_csv(f"{RUNS}/{run}/K562_Xu_pl/{m}/model/model_coefficients.tsv", sep="\t")
            t = t[t["test_chr"] == "none"].assign(run=run, model=m)
            rec.append(t[["run", "model", "feature", "coefficient"]])
    pd.concat(rec).to_csv(cf, sep="\t", index=False)
    print(f"wrote {cf}")
C = pd.read_csv(cf, sep="\t")
ab = C[C["run"] == "2026_1004b_k562_xu_abcfp"]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"wspace": 0.4})
for off, s, col in ((-0.06, "tobias", BLUE), (0.06, "print", RED)):
    xs, est, lo, hi, co = [0], [0.0], [0.0], [0.0], []
    co.append(ab.loc[(ab["model"] == B) & (ab["feature"] == "ARC.E2G.Score"), "coefficient"].item())
    for l in (1, 4):
        r = ba.loc[f"{T}_{s}_l{l} - {B}"]
        xs.append(l); est.append(r["estimate"]); lo.append(r["ci_lo"]); hi.append(r["ci_hi"])
        co.append(ab.loc[(ab["model"] == f"{T}_{s}_l{l}") & ab["feature"].str.startswith("ARCfp"), "coefficient"].item())
    xs = np.array(xs, float) + off
    est, lo, hi = map(np.array, (est, lo, hi))
    ax[0].errorbar(xs, est, yerr=[est - lo, hi - est], fmt="o-", ms=4, lw=1.1, color=col, label=s.upper())
    ax[1].plot(xs - off, co, "o-", ms=4, lw=1.1, color=col, label=s.upper())
ax[0].axhline(0, color="black", lw=0.7)
for j in (0, 1):
    ax[j].set_xticks([0, 1, 4])
    ax[j].set_xlabel("Footprint weight (0 = ATAC only)")
ax[0].set_ylabel("CRISPR AUPRC minus baseline (95% CI)")
ax[0].set_title("a  benchmark", loc="left")
ax[1].set_ylabel("Model coefficient of the activity score")
ax[1].set_title("b  weight the model gives the score", loc="left")
ax[0].legend(frameon=False, fontsize=7, loc="lower left")
save(fig, "fig43_footprint_abc_activity")

# ---------------------------------------------------------------- fig44 feature A
feats = [("kA_max", "own-motif TFs, max"), ("kA_max_perm", "motif-TF labels permuted, max"),
         ("kA_any_max", "any of 348 TFs, max"), ("kA_ery_max", "erythroid TFs, max"),
         ("kA_sumpos", "own-motif TFs, sum of positive"), ("n_tf", "number of TFs with a motif")]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.9})
y = np.arange(len(feats))[::-1]
for off, z, col, lab in ((0.14, False, BLUE, "controls: Xu ATAC"), (-0.14, True, LBLUE, "controls: Xu + ENCODE ATAC")):
    for yi, (f, _) in zip(y, feats):
        r = kend1[z].loc[f]
        ci_dot(ax[0], r["partial_rho"], yi + off, r["ci_lo"], r["ci_hi"], col, lab if yi == y[0] else None, ms=4)
    t = kend[z]
    y2 = np.arange(len(t))[::-1]
    for yi, (_, r) in zip(y2, t.iterrows()):
        ci_dot(ax[1], r["partial_rho_kA_max"], yi + off, r["ci_lo"], r["ci_hi"], col, ms=4)
ax[0].set_yticks(y)
ax[0].set_yticklabels([f[1] for f in feats])
ax[0].set_title("a  each feature, beyond ATAC and composition", loc="left", fontsize=9)
ax[1].set_yticks(np.arange(3)[::-1])
ax[1].set_yticklabels(["+ number of TFs", "+ permuted-label score", "+ any-TF score"])
ax[1].set_title("b  own-motif max, adding controls in turn", loc="left", fontsize=9)
for j in (0, 1):
    ax[j].axvline(0, color="black", lw=0.7)
    ax[j].set_xlabel("Partial Spearman with K562 p300 (95% CI)")
ax[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2)
save(fig, "fig44_peak_tf_kendall")
