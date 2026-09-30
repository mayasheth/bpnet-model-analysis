#!/usr/bin/env python3
"""Figures for report 6: improving the transferred p300 model (F-030 to F-038, 4.49-4.52).

Every panel reads a results file on disk; nothing is hand-typed. The per-run seed AUPRCs are
read from the 4.17 / 4.34 logs, which print each replicate group's members (the TSVs keep only
the pairwise contrasts).

Usage: 5.3.plot_tf_activity_report_figures.py [--outdir figures]
"""
import argparse
import os
import re
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
R = f"{P}/results"
ap = argparse.ArgumentParser()
ap.add_argument("--outdir", default=f"{P}/figures")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
})
COL = {"rank": "#08519c", "nbr": "#969696", "oracle": "#7f2704", "track": "#d94801",
       "seed": "#a50f15", "ref": "#54278f"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {a.outdir}/{name}.png")


def row(tsv, A, B):
    """One comparison from a 4.17 (arm_a/arm_b) or 4.34 (model_a/model_b) table."""
    t = pd.read_csv(f"{R}/{tsv}", sep="\t")
    ca, cb = ("arm_a", "arm_b") if "arm_a" in t else ("model_a", "model_b")
    s = f".ABC.Score" if ca == "arm_a" else "_powerlaw_all"
    fix = lambda x: x if x in ("mse", "oracle", "ctcf") or x.endswith(s) else x + s
    hit = t[(t[ca] == fix(A)) & (t[cb] == fix(B))]
    assert len(hit) == 1, f"{tsv}: {A} vs {B} not found"
    h = hit.iloc[0]
    return float(h["delta"]), float(h["ci_lo"]), float(h["ci_hi"]), h


MC, MCR = "p300pred_multicell_holdK562", "ATAC_p300multicell"

# ======================================================== fig 26: every variant, headline
ROWS = [
    ("Ranking loss added to training", None, None, "rank"),
    ("  ranking weight 2", ("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank2", MC),
     ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank2", MCR), "rank"),
    ("  ranking weight 10", ("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank10", MC),
     ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank10", MCR), "rank"),
    ("  ranking weight 20", ("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank20", MC),
     ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank20", MCR), "rank"),
    ("  ranking weight 50", ("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank50", MC),
     ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank50", MCR), "rank"),
    ("  ranking weight 100", ("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank100", MC),
     ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank100", MCR), "rank"),
    ("Neighbouring elements (re-ranking)*", None, None, "nbr"),
    ("  same-class neighbours", ("2026_0925_nbr_auprc_pairs.tsv", "p300pred_nbr_multicell_boostnbr",
                                 "p300pred_nbr_multicell_boostnonbr"), None, "nbr"),
    ("  neighbours of both classes", ("2026_0926_nbrx_auprc_pairs.tsv",
                                      "p300pred_nbrx_multicell_boostnbr",
                                      "p300pred_nbr_multicell_boostnonbr"), None, "nbr"),
    ("GATA elements given better ranks (oracle)", None, None, "oracle"),
    ("  in-cell model's ranks, GATA elements", ("2026_0927_gataoracle_auprc_pairs.tsv",
                                                "p300pred_gataoracle_gata_incell", MC), None, "oracle"),
    ("  same, matched non-GATA elements", ("2026_0927_gataoracle_auprc_pairs.tsv",
                                           "p300pred_gataoracle_ctrl_incell", MC), None, "oracle"),
    ("  observed-p300 ranks, GATA elements", ("2026_0927_gataoracle_auprc_pairs.tsv",
                                              "p300pred_gataoracle_gata_obs", MC), None, "oracle"),
    ("  same, matched non-GATA elements", ("2026_0927_gataoracle_auprc_pairs.tsv",
                                           "p300pred_gataoracle_ctrl_obs", MC), None, "oracle"),
    ("TF-activity track as an input (one run each)", None, None, "track"),
    ("  true TF list", ("2026_0928_motifact_auprc_pairs.tsv", f"{MC}_motifact_oracle", MC),
     ("2026_0928_re2g_motifact_auprc_pairs.tsv", "ATAC_p300mcmotifact", MCR), "track"),
    ("  ATAC TF list, CTCF removed", ("2026_0928_motifact2_auprc_pairs.tsv",
                                      f"{MC}_motifact_assoc_noctcf", MC),
     ("2026_0928_re2g_motifact2_auprc_pairs.tsv", "ATAC_p300mcmotifctcf", MCR), "track"),
    ("  true list + ranking loss", ("2026_0928_motifact2_auprc_pairs.tsv",
                                    f"{MC}_rank10_motifact_oracle", MC),
     ("2026_0928_re2g_motifact2_auprc_pairs.tsv", "ATAC_p300mcmotifrank", MCR), "track"),
    ("  true list, 'head'", ("2026_0929_teach_auprc_pairs.tsv", f"{MC}_motifact_oracle_head", MC),
     ("2026_0929_re2g_teach_auprc_pairs.tsv", "ATAC_p300mcteachhead", MCR), "track"),
    ("  true list, 'pair'", ("2026_0929_teach_auprc_pairs.tsv", f"{MC}_motifact_oracle_pair", MC),
     ("2026_0929_re2g_teach_auprc_pairs.tsv", "ATAC_p300mcteachpair", MCR), "track"),
    ("  true list, 'headpair'", ("2026_0929_teach_auprc_pairs.tsv",
                                 f"{MC}_motifact_oracle_headpair", MC),
     ("2026_0929_re2g_teach_auprc_pairs.tsv", "ATAC_p300mcteachheadpair", MCR), "track"),
    ("  true list, two-stage", ("2026_0930_seeds_auprc_pairs.tsv",
                                f"{MC}_motifact_oracle_twostage", MC),
     ("2026_0930_re2g_seeds_auprc_pairs.tsv", "ATAC_p300mctwostageoracle", MCR), "track"),
    ("  ATAC list, two-stage", ("2026_0930_seeds_auprc_pairs.tsv",
                                f"{MC}_motifact_assoc_noctcf_twostage", MC),
     ("2026_0930_re2g_seeds_auprc_pairs.tsv", "ATAC_p300mctwostagectcf", MCR), "track"),
    ("TF-activity track, mean of 3 training runs", None, None, "seed"),
    ("  true TF list (3 runs vs 3 runs)", ("2026_0930_seeds_auprc_pairs.tsv", "oracle", "mse"),
     ("2026_0930_re2g_seeds_auprc_pairs.tsv", "oracle", "mse"), "seed"),
    ("  ATAC TF list (3 runs vs 3 runs)", ("2026_0930_seeds_auprc_pairs.tsv", "ctcf", "mse"),
     ("2026_0930_re2g_seeds_auprc_pairs.tsv", "ctcf", "mse"), "seed"),
]
incell_abc = row("2026_0927_gataoracle_auprc_pairs.tsv", "p300pred_k562_multimodal", MC)[0]
floor_abc = -row("2026_0926_ranksweep_auprc_pairs.tsv", MC, "K562_ATAC_only")[0]
h = row("2026_0928_re2g_motifact_auprc_pairs.tsv", "ATAC_p300mcmotifact", "ATAC_p300incell")[3]
h2 = row("2026_0928_re2g_motifact_auprc_pairs.tsv", "ATAC_p300mcmotifact", MCR)[3]
incell_re2g = float(h["auprc_b"]) - float(h2["auprc_b"])
floor_re2g = -row("2026_0926_re2g_ranksweep_auprc_pairs.tsv", MCR, "ATAC")[0]

CORR = pd.read_csv(f"{R}/2026_0930_k562_correlation_by_arm.tsv", sep="\t").set_index("arm")
GROUPS = {"mse": [MC, f"{MC}_seed1", f"{MC}_seed2"],
          "oracle": [f"{MC}_motifact_oracle", f"{MC}_motifact_oracle_seed1", f"{MC}_motifact_oracle_seed2"],
          "ctcf": [f"{MC}_motifact_assoc_noctcf", f"{MC}_motifact_assoc_noctcf_seed1",
                   f"{MC}_motifact_assoc_noctcf_seed2"]}


def corr_delta(spec, col):
    """Correlation of the row's arm minus its comparator's, on the same K562 candidates."""
    _, A_, B_ = spec
    val = lambda x: CORR.loc[GROUPS[x], col].mean() if x in GROUPS else CORR.loc[x, col]
    return val(A_) - val(B_)


fig, axes = plt.subplots(1, 3, figsize=(14.5, 8.2), sharey=True,
                         gridspec_kw={"width_ratios": [0.8, 1, 1]})
n = len(ROWS)
ax = axes[0]
for i, (name, abc, re2g, grp) in enumerate(ROWS):
    if abc is None:
        continue
    y = n - 1 - i
    ax.plot(corr_delta(abc, "partial"), y, "o", color=COL[grp], ms=6,
            mec="k" if grp == "seed" else COL[grp])
    ax.plot(corr_delta(abc, "spearman"), y, "D", color="none", mec=COL[grp], ms=5)
ax.axvline(0, color="k", lw=0.8)
ax.plot([], [], "o", color="#555", label="beyond ATAC (partial)")
ax.plot([], [], "D", color="none", mec="#555", label="overall")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
ax.set_xlabel("Spearman with observed K562 p300,\nminus the unmodified multi-cell model")
ax.set_title("a  Correlation over 153,545 elements", loc="left")
for j, (ax, lab) in enumerate(zip(axes[1:], ("ABC", "ENCODE-rE2G"))):
    for i, (name, abc, re2g, grp) in enumerate(ROWS):
        y = n - 1 - i
        spec = abc if j == 0 else re2g
        if spec is None:
            if name.startswith("  "):
                ax.text(0.0, y, "not tested", fontsize=7, color="#999", va="center", ha="center")
            continue
        d, lo, hi, _ = row(*spec)
        ax.plot([lo, hi], [y, y], color=COL[grp], lw=1.6)
        ax.plot(d, y, "o", color=COL[grp], ms=5.5 if grp != "seed" else 7,
                mec="k" if grp == "seed" else COL[grp])
    ax.axvline(0, color="k", lw=0.8)
    ref, flo = (incell_abc, floor_abc) if j == 0 else (incell_re2g, floor_re2g)
    ax.axvline(ref, color=COL["ref"], ls="--", lw=1)
    ax.axvline(flo, color="#bdbdbd", ls="--", lw=1)
    ax.text(ref, n - 0.2, "in-cell\nmodel", color=COL["ref"], fontsize=7, ha="center", va="bottom")
    ax.text(flo, n - 0.2, "ATAC\nonly", color="#888", fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("AUPRC minus the unmodified multi-cell model")
    ax.set_title(f"{'bc'[j]}  {lab} AUPRC, K562 CRISPR benchmark", loc="left")
axes[0].set_yticks(range(n))
axes[0].set_yticklabels([r[0] for r in ROWS][::-1], fontsize=8)
for t, r_ in zip(axes[0].get_yticklabels(), ROWS[::-1]):
    if not r_[0].startswith("  "):
        t.set_fontweight("bold")
        t.set_color(COL[r_[3]])
axes[0].set_ylim(-0.7, n + 0.9)
fig.text(0.01, -0.01, "* compared with the same re-ranking fit without neighbour terms",
         fontsize=7, color="#666")
save(fig, "fig26_variant_forest")

# ======================================================== fig 27: what each variant changes
fig, ax = plt.subplots(figsize=(10.5, 4.6))
ax.set_xlim(0, 100)
ax.set_ylim(0, 50)
ax.axis("off")


def box(x, y, w, h, text, fc="#f0f0f0", ec="#555", fs=8, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                fc=fc, ec=ec, lw=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal")


def arrow(x0, y0, x1, y1, c="#555"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1))


box(1, 36, 13, 8, "DNA sequence\n(2 kb window)", bold=True)
box(1, 24, 13, 8, "ATAC signal\n(target cell)", bold=True)
box(1, 10, 13, 9, "TF-activity track\n(motif hits x how\nactive each TF is)", fc="#fee6ce",
    ec=COL["track"], bold=True)
box(24, 24, 18, 18, "Neural network\n(BPNet-style,\ntrained on 4 other\ncell types)", fc="#deebf7",
    ec=COL["rank"], bold=True)
box(52, 29, 14, 8, "predicted p300\nper element", bold=True)
box(76, 29, 22, 8, "ABC / ENCODE-rE2G\nscored on K562 CRISPR", bold=True)
arrow(14, 40, 24, 36)
arrow(14, 28, 24, 30)
arrow(14, 15, 24, 26, COL["track"])
arrow(42, 33, 52, 33)
arrow(66, 33, 76, 33)
box(21, 45, 24, 4.5, "ranking loss: also score the ORDER", fc="#fff", ec=COL["rank"], fs=7.5)
box(52, 41, 14, 7, "neighbour re-ranking:\nuse elements within\n1 Mb", fc="#fff", ec=COL["nbr"],
    fs=7.2)
box(52, 18, 14, 8, "GATA oracle:\ngive GATA elements\nbetter ranks", fc="#fff", ec=COL["oracle"],
    fs=7.2)
arrow(59, 26, 59, 29, COL["oracle"])
arrow(59, 41, 59, 37, COL["nbr"])
box(22, 1, 17, 12, "'head': track skips\nthe network, adds to\nthe count through\n4 weights", fc="#fff",
    ec=COL["track"], fs=7)
box(41, 1, 17, 12, "'pair': also score the\np300 DIFFERENCE for\nthe same DNA in two\ntraining cells", fc="#fff",
    ec=COL["track"], fs=7)
box(60, 1, 18, 12, "'two-stage': freeze the\nnetwork, train a small\nDNA-free correction\nfrom track + ATAC", fc="#fff",
    ec=COL["track"], fs=7)
box(80, 1, 18, 12, "TF list 'true': from\nthe cell's own p300\n(upper bound)\nTF list 'ATAC': from\nits ATAC only", fc="#fff",
    ec=COL["track"], fs=7)
ax.set_title("Where each variant changes the pipeline", loc="left")
save(fig, "fig27_variants_schematic")

# ======================================================== fig 28: rank sweep and neighbours
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.4))
ax = axes[0]
W = [2, 10, 20, 50, 100]
for tsv, s, lab, c, off in (("2026_0926_ranksweep_auprc_pairs.tsv", f"{MC}_rank", "ABC", COL["rank"], -0.04),
                            ("2026_0926_re2g_ranksweep_auprc_pairs.tsv", "ATAC_p300mcrank", "ENCODE-rE2G", "#6baed6", 0.04)):
    base_ = MC if "re2g" not in tsv else MCR
    vals = [row(tsv, f"{s}{w}", base_) for w in W]
    x = np.log10(W) + off
    ax.errorbar(x, [v[0] for v in vals], yerr=[[v[0] - v[1] for v in vals], [v[2] - v[0] for v in vals]],
                fmt="o-", color=c, capsize=2.5, label=lab, ms=5)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(np.log10(W))
ax.set_xticklabels(W)
ax.set_xlabel("weight on the ranking term (MSE term kept)")
ax.set_ylabel("AUPRC minus MSE-only model")
ax.set_title("a  Ranking-loss weight sweep (one run each)", loc="left")
ax.legend(frameon=False, fontsize=8)
ax = axes[1]
nb = pd.read_csv(f"{R}/2026_0924_neighbour_test_crossclass.tsv", sep="\t")
nb = nb[nb.cls == "nonpromoter"]
arms = [("multi-cell, K562 held out", "p300pred_nbrx_multicell_boostnbr", "p300pred_nbr_multicell_boostnonbr", "multi-cell"),
        ("K562 in-cell", "p300pred_nbrx_incell_boostnbr", "p300pred_nbr_incell_boostnonbr", "in-cell"),
        ("ATAC only (pred = ATAC)", "p300pred_nbrx_atac_boostnbr", "K562_ATAC_only", "ATAC only")]
for arm, A_, B_, lab in arms:
    s = nb[nb.arm == arm].set_index("model")["partial_rho"]
    dr = s["M3 boosting, + neighbour terms"] - s["M3a boosting, no neighbours"]
    d, lo, hi, _ = row("2026_0926_nbrx_auprc_pairs.tsv", A_, B_)
    ax.errorbar(dr, d, yerr=[[d - lo], [hi - d]], fmt="o", color="#555", capsize=3, ms=6)
    ax.annotate(lab, (dr, d), xytext=(6, 4), textcoords="offset points", fontsize=8)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlim(0, 0.26)
ax.set_xlabel("gain in partial correlation with p300 (given ATAC)")
ax.set_ylabel("ABC AUPRC change")
ax.set_title("b  Neighbour terms: better correlation, no ABC gain", loc="left")
save(fig, "fig28_rank_and_neighbour")

# ======================================================== fig 29: error anatomy + GATA oracle
ea = pd.read_csv(f"{R}/2026_0925_p300_error_anatomy.features.tsv", sep="\t")
ea = ea[(ea.cls == "nonpromoter") & ea.feature.str.startswith("motif:")].copy()
ea["name"] = ea.feature.str.split(" ", n=1).str[1]
top = pd.concat([ea.nlargest(8, "coef_r_transfer"), ea.nsmallest(6, "coef_r_transfer")[::-1]])
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axes[0]
yy = np.arange(len(top))[::-1]
ax.barh(yy, top["coef_r_transfer"], color=np.where(top["coef_r_transfer"] > 0, "#cb181d", "#2171b5"),
        height=0.7)
ax.set_yticks(yy)
ax.set_yticklabels(top["name"], fontsize=8)
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("transfer-specific ranking error (percentile points)\n+ = transferred model ranks these elements too low")
ax.set_title("a  Which motifs the transferred model mis-ranks", loc="left")
ax = axes[1]
g = [("in-cell\nmodel ranks", "gata_incell", "ctrl_incell"), ("observed\np300 ranks", "gata_obs", "ctrl_obs")]
for k, (lab, ga, ct) in enumerate(g):
    for m, (arm, c, nm) in enumerate(((ga, "#7f2704", "GATA elements"), (ct, "#fdae6b", "matched non-GATA"))):
        d, lo, hi, _ = row("2026_0927_gataoracle_auprc_pairs.tsv", f"p300pred_gataoracle_{arm}", MC)
        x = k + (m - 0.5) * 0.36
        ax.bar(x, d, width=0.34, color=c, label=nm if k == 0 else None)
        ax.errorbar(x, d, yerr=[[d - lo], [hi - d]], color="k", capsize=3, lw=0.9)
ax.axhline(incell_abc, color=COL["ref"], ls="--", lw=1)
ax.text(1.45, incell_abc, " in-cell\n model", color=COL["ref"], fontsize=7, va="center")
ax.set_xticks([0, 1])
ax.set_xticklabels([x[0] for x in g])
ax.set_ylabel("ABC AUPRC minus multi-cell model")
ax.set_title("b  Re-ranking GATA elements recovers half the gap", loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
save(fig, "fig29_error_anatomy_gata")

# ======================================================== fig 30: TF-list estimators
lv = pd.read_csv(f"{R}/2026_0928_learned_importance_variants.scores.tsv", sep="\t")
cf = pd.read_csv(f"{R}/2026_0929_corrected_tf_list.scores.tsv", sep="\t")
EST = [("ATAC association (raw)", cf, "raw"), ("  CTCF removed", cf, "noctcf"),
       ("  flagged motifs removed", cf, "flagged"), ("  GC regressed out", cf, "gc"),
       ("  flagged + GC", cf, "gc+flagged"),
       ("learned (gradient boosting)", lv, "learned (4.47)"),
       ("learned + ATAC association", lv, "stack: learned + assoc"),
       ("average + predicted departure", lv, "departure: prior + model"),
       ("gap-corrected (panel-relative)", cf, "gap")]
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), sharey=True)
for j, (col, lab) in enumerate((("spearman", "whole ranking vs true importance"),
                                ("dep_spearman", "cell-specific part vs true"))):
    ax = axes[j]
    v = [t[t.estimate == e][col].mean() for _, t, e in EST]
    c = ["#d94801" if nm.startswith("  ") or nm.startswith("ATAC") else
         ("#969696" if "gap" in nm else "#6a51a3") for nm, _, _ in EST]
    yy = np.arange(len(EST))[::-1]
    ax.barh(yy, v, color=c, height=0.7)
    for y_, vv in zip(yy, v):
        ax.text(vv + 0.01, y_, f"{vv:.2f}", va="center", fontsize=7.5)
    ax.set_xlim(0, 0.8)
    ax.set_xlabel("Spearman, mean over 5 held-out cells")
    ax.set_title(f"{'ab'[j]}  {lab}", loc="left")
axes[0].set_yticks(np.arange(len(EST))[::-1])
axes[0].set_yticklabels([e[0] for e in EST], fontsize=8)
save(fig, "fig30_tf_list_estimators")

# ======================================================== fig 31: residuals by model
tt = pd.read_csv(f"{R}/2026_0928_p300_error_anatomy_teach.features.tsv", sep="\t")
tt = tt[(tt.cls == "nonpromoter") & tt.feature.str.startswith("motif:")].copy()
tt["name"] = tt.feature.str.split(" ", n=1).str[1]
FAM = {"GATA family": tt.name.str.upper().str.contains("GATA") | tt.name.str.upper().str.contains("TRPS1"),
       "CEBP family": tt.name.str.upper().str.startswith("CEBP"),
       "ESR1/2": tt.name.str.upper().str.match(r"^ESR[12]$")}
MODELS = [("in-cell\n(K562-trained)", "incell"), ("multi-cell,\nno changes", "multicell"),
          ("+ ranking\nloss w=10", "rank10"), ("+ true\nTF list", "oracle"), ("+ ATAC\nTF list", "ctcf"),
          ("true list,\n'head'", "head"), ("true list,\n'pair'", "pair"), ("true list,\n'headpair'", "headpair")]
Mx = np.array([[tt.loc[m, f"coef_r_{k}"].mean() for _, k in MODELS] for m in FAM.values()]
              + [[tt[f"coef_r_{k}"].abs().mean() for _, k in MODELS]])
fig, ax = plt.subplots(figsize=(10.5, 2.9))
vmax = 6
im = ax.imshow(Mx[:3], cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
for i in range(3):
    for j in range(len(MODELS)):
        ax.text(j, i, f"{Mx[i, j]:+.1f}", ha="center", va="center", fontsize=8.5,
                color="w" if abs(Mx[i, j]) > 4 else "k")
ax.set_xticks(range(len(MODELS)))
ax.set_xticklabels([m[0] for m in MODELS], fontsize=8)
ax.set_yticks(range(3))
ax.set_yticklabels(list(FAM), fontsize=8.5)
ax.set_title("Ranking error at elements carrying each motif family (percentile points; "
             "+ = ranked too low, - = too high)", loc="left", fontsize=9.5)
for j in range(len(MODELS)):
    ax.text(j, 3.05, f"{Mx[3, j]:.2f}", ha="center", va="top", fontsize=7.5, color="#444")
ax.text(-0.6, 3.05, "mean |error|,\nall motifs", ha="right", va="top", fontsize=7.5, color="#444")
fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
save(fig, "fig31_channel_residuals")

# ======================================================== fig 32: seed replicates


def groups(log):
    out = {}
    for line in open(log):
        m = re.match(r"group (\w+): ([0-9., ]+)\s+->", line)
        if m:
            out[m.group(1)] = [float(x) for x in m.group(2).split(",")]
    assert set(out) == {"mse", "oracle", "ctcf"}, f"{log}: groups {sorted(out)}"
    return out


GA = groups(f"{P}/log/seeds_boot.46073940.txt")
GR = groups(f"{P}/log/re2g_seeds_boot.46073942.txt")
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.3))
LAB = [("mse", "multi-cell,\nno TF track", "#238b45"), ("oracle", "+ true\nTF list", "#7f2704"),
       ("ctcf", "+ ATAC\nTF list", "#d94801")]
for j, (Gd, lab) in enumerate(((GA, "ABC"), (GR, "ENCODE-rE2G"))):
    ax = axes[j]
    for k, (key, nm, c) in enumerate(LAB):
        v = Gd[key]
        ax.scatter([k] * 3 + np.array([-0.08, 0, 0.08]), v, color=c, s=34, zorder=3)
        ax.plot([k - 0.22, k + 0.22], [np.mean(v)] * 2, color="k", lw=1.6)
        ax.text(k + 0.26, np.mean(v), f"{np.mean(v):.4f}\nsd {np.std(v, ddof=1):.4f}", fontsize=7,
                va="center")
    ax.set_xticks(range(3))
    ax.set_xticklabels([x[1] for x in LAB])
    ax.set_xlim(-0.5, 2.9)
    ax.set_ylabel("AUPRC")
    ax.set_title(f"{'ab'[j]}  {lab}: three training runs per model", loc="left")
save(fig, "fig32_seed_replicates")

# ======================================================== fig 33: panel diversity heatmap
src = f"{R}/2026_0929_panel_tf_diversity_loose.heatmap.png"
shutil.copy(src, f"{a.outdir}/fig33_panel_tf_diversity.png")
print(f"copied {src} -> {a.outdir}/fig33_panel_tf_diversity.png")

# ======================================================== fig 34: nucleosome signature
PRn = pd.read_csv(f"{R}/2026_0930_nucleosome_signature.partial_rho.tsv", sep="\t", index_col=0)
TFn = pd.read_csv(f"{R}/2026_0930_nucleosome_signature.tf_list_scores.tsv", sep="\t")
NAMES = {"nfr": "short fragments at the centre", "short_frac": "fraction of short fragments",
         "mono_centre": "nucleosome over the centre", "centre_occ": "centre nucleosome / open",
         "mono_flank": "flanking nucleosomes", "flank_ratio": "flanking nucleosomes / open",
         "di_flank": "di-nucleosome fragments, flanks", "nfr_width": "open-region width",
         "positioning": "flanking nucleosome positioning"}
cells = ["K562", "GM12878", "HepG2", "A549", "MCF-7"]
PRn = PRn.loc[list(NAMES), cells]
fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.8), gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.55})
ax = axes[0]
im = ax.imshow(PRn.values, cmap="RdBu_r", vmin=-0.2, vmax=0.2, aspect="auto")
for i in range(PRn.shape[0]):
    for j in range(PRn.shape[1]):
        ax.text(j, i, f"{PRn.values[i, j]:+.2f}", ha="center", va="center", fontsize=7.5)
ax.set_xticks(range(len(cells)))
ax.set_xticklabels(cells)
ax.set_yticks(range(len(NAMES)))
ax.set_yticklabels(list(NAMES.values()), fontsize=8)
ax.set_title("a  Partial Spearman with p300, beyond ATAC", loc="left")
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
ax = axes[1]
EN = [("raw ATAC association", "ATAC"), ("nucleosome score association", "nucleosome\nsignature"),
      ("ATAC + nucleosome (z sum)", "ATAC +\nnucleosome")]
x = np.arange(len(EN))
for k, (col, lab, c) in enumerate((("spearman", "whole ranking", "#6a51a3"),
                                   ("dep_spearman", "cell-specific part", "#d94801"))):
    v = [TFn[TFn.estimate == e][col].mean() for e, _ in EN]
    ax.bar(x + (k - 0.5) * 0.36, v, width=0.34, color=c, label=lab)
    for xx, vv in zip(x + (k - 0.5) * 0.36, v):
        ax.text(xx, vv + 0.01, f"{vv:.2f}", ha="center", fontsize=7.5)
ax.set_xticks(x)
ax.set_xticklabels([e[1] for e in EN])
ax.set_ylabel("Spearman vs true TF importance,\nmean over 5 cells")
ax.set_ylim(0, 0.7)
ax.legend(frameon=False, fontsize=8)
ax.set_title("b  As a TF-list estimator", loc="left")
save(fig, "fig34_nucleosome_signature")
