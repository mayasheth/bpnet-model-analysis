#!/usr/bin/env python3
"""Figures for report 5: multi-cell-type p300 transfer, upstream and downstream.

Every panel reads a file on disk; nothing is hand-typed except the curated cross-experiment
table in figure 5, whose every entry carries a finding ID so it can be checked against
`.living/findings/`.

Usage: 5.1.plot_multicell_report_figures.py [--outdir figures]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
CC = ("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/"
      "workflow/results/2026_0922_p300_multicell")

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", default=f"{P}/figures")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 160,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
})

C = {
    "floor":     "#bdbdbd",
    "k27ceil":   "#404040",
    "obs":       "#54278f",
    "incell":    "#807dba",
    "own":       "#cb181d",
    "tgtnorm":   "#fd8d3c",
    "multicell": "#238b45",
}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {a.outdir}/{name}.png")


# --------------------------------------------------------------- fig 1: the panel
# Training peak counts come from F-026's ledger rather than from wc -l on the peak files:
# the trainer drops peaks whose window runs off a contig or into the blacklist, so the file
# line count (21,068 for GM12878) overstates what the model saw (14,602) by a constant ~31%.
TRAIN_PEAKS = {"GM12878": 14602, "MCF-7": 8279, "A549": 6130, "HepG2": 5658}
HELD_OUT = ("K562", 28532)

depth = pd.read_csv(f"{P}/results/panel_depth_audit.tsv", sep="\t")
atac = (depth[depth.assay == "ATAC"].groupby("cell")["reads"].sum() / 1e6)
# The filtered BAMs double-count: prefer the unfiltered row per accession where both exist.
accstats = json.load(open(
    f"{P}/models/p300_multicell_holdK562_hw500_clw10/fold0/acc_normalization_stats.json"
))["per_cell_type"]

cells = ["GM12878", "MCF-7", "A549", "HepG2"]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))

ax = axes[0]
vals = [TRAIN_PEAKS[c] for c in cells]
ax.bar(cells, vals, color=C["multicell"], width=0.68)
ax.bar([""], [HELD_OUT[1]], color=C["floor"], width=0.68)
ax.text(4, HELD_OUT[1] * 1.02, f"K562\nheld out\n{HELD_OUT[1]:,}", ha="center", va="bottom",
        fontsize=7.5, color="#555")
ax.set_ylabel("training peaks")
ax.set_title(f"a  Pooled training set, {sum(vals):,} peaks", loc="left")
ax.set_ylim(0, HELD_OUT[1] * 1.45)
for i, v in enumerate(vals):
    ax.text(i, v + 400, f"{v:,}", ha="center", fontsize=7.5)
ax.tick_params(axis="x", labelrotation=30)

ax = axes[1]
d = [atac.get(c, np.nan) for c in cells]
ax.bar(cells, d, color=["#9ecae1" if c != "MCF-7" else "#cb181d" for c in cells], width=0.68)
ax.set_ylabel("ATAC reads (millions)")
ax.set_title("b  ATAC library depth", loc="left")
for i, v in enumerate(d):
    ax.text(i, v + 12, f"{v:,.0f}M", ha="center", fontsize=7.5)
ax.set_ylim(0, max(d) * 1.2)
ax.tick_params(axis="x", labelrotation=30)

ax = axes[2]
m = [accstats[c]["acc_mean"] for c in cells]
s = [accstats[c]["acc_std"] for c in cells]
ax.errorbar(range(len(cells)), m, yerr=s, fmt="o", color="#333", capsize=3, ms=5)
ax.scatter([cells.index("MCF-7")], [accstats["MCF-7"]["acc_mean"]], color="#cb181d", zorder=5,
           s=42)
ax.set_xticks(range(len(cells)))
ax.set_xticklabels(cells, rotation=30)
ax.set_ylabel("log1p accessibility, mean $\\pm$ sd")
ax.set_title("c  Accessibility mean per cell type", loc="left")
ax.axhline(np.mean([accstats[c]["acc_mean"] for c in cells if c != "MCF-7"]), ls=":",
           color="#999", lw=1)
ax.annotate("MCF-7 is 2 log units lower;\na single pooled mean would\nleave that offset in the input",
            xy=(1, accstats["MCF-7"]["acc_mean"]), xytext=(1.35, 3.2), fontsize=7,
            color="#cb181d",
            arrowprops=dict(arrowstyle="->", color="#cb181d", lw=0.9))
save(fig, "fig18_multicell_panel")


# ------------------------------------------------- fig 2: upstream, per fold
pf = pd.read_csv(f"{P}/results/p300_multicell_k562_per_fold.tsv", sep="\t")
ORDER = [("K562_incell", "K562 in-cell", C["incell"]),
         ("multicell_holdK562", "multi-cell, K562 held out", C["multicell"]),
         ("GM_single_corrected_tgtnorm", "single-source, target norm", C["tgtnorm"]),
         ("GM_single_corrected", "single-source, own stats", C["own"]),
         ("K562_atac_only_floor", "accessibility-only floor", C["floor"])]

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.4))
for ax, metric, lab in zip(axes, ["overall_pearson", "overall_pearson_topq"],
                           ["overall Pearson", "top-quintile Pearson"]):
    for i, (key, name, col) in enumerate(ORDER):
        v = pf[pf.config == key][metric].values
        ax.scatter(v, [i] * len(v), color=col, s=26, alpha=0.75, zorder=3)
        ax.scatter([v.mean()], [i], color=col, s=95, marker="|", zorder=4, lw=2.2)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([n for _, n, _ in ORDER])
    ax.invert_yaxis()
    ax.set_xlabel(lab)
    ax.grid(axis="x", color="#eee", zorder=0)
axes[0].set_title("a  Counts correlation with observed K562 p300, 5 folds", loc="left")
axes[1].set_title("b  Top signal quintile", loc="left")
save(fig, "fig19_multicell_upstream")


# ------------------------------------------------- fig 3: PRC curves
pr = pd.read_csv(f"{CC}/pr_table.txt.gz", sep="\t")
CURVES = [("p300obs_k562.ABC.Score", "observed p300 (ceiling)", C["obs"]),
          ("p300pred_k562_multimodal.ABC.Score", "predicted p300, K562 in-cell", C["incell"]),
          ("p300pred_gm12878_corrected_tgtnorm.ABC.Score", "single-source, target norm",
           C["tgtnorm"]),
          ("p300pred_multicell_holdK562.ABC.Score", "multi-cell, K562 held out",
           C["multicell"]),
          ("K562_ATAC_only.ABC.Score", "ATAC only (floor)", C["floor"])]
fig, ax = plt.subplots(figsize=(5.4, 4.3))
for uid, name, col in CURVES:
    d = pr[pr.pred_uid == uid].dropna(subset=["precision", "recall"]).sort_values("recall")
    lw = 2.1 if "multicell" in uid else 1.5
    ax.plot(d.recall, d.precision, color=col, label=name, lw=lw)
ax.set_xlabel("recall")
ax.set_ylabel("precision")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(fontsize=7.6, frameon=False, loc="upper right")
ax.set_title("Precision-recall, K562 CRISPR, shared pair set", loc="left")
save(fig, "fig20_multicell_prc")


# ------------------------------------------------- fig 4: the forest plot
pairs = pd.read_csv(f"{P}/results/2026_0922_multicell_auprc_pairs.tsv", sep="\t")


def row(arm_a, arm_b):
    r = pairs[(pairs.arm_a == arm_a) & (pairs.arm_b == arm_b)]
    assert len(r) == 1, f"{arm_a} - {arm_b}"
    return r.iloc[0]


FOREST = [
    ("multi-cell $-$ single-source\n(policy-matched)", "p300pred_multicell_holdK562.ABC.Score",
     "p300pred_gm12878_corrected_tgtnorm.ABC.Score", C["multicell"], True),
    ("multi-cell $-$ floor", "p300pred_multicell_holdK562.ABC.Score",
     "K562_ATAC_only.ABC.Score", C["multicell"], False),
    ("single-source $-$ floor", "p300pred_gm12878_corrected_tgtnorm.ABC.Score",
     "K562_ATAC_only.ABC.Score", C["tgtnorm"], False),
    ("target norm $-$ own stats", "p300pred_gm12878_corrected_tgtnorm.ABC.Score",
     "p300pred_gm12878_corrected.ABC.Score", C["own"], False),
    ("K562 in-cell $-$ floor", "p300pred_k562_multimodal.ABC.Score",
     "K562_ATAC_only.ABC.Score", C["incell"], False),
    ("observed p300 $-$ floor", "p300obs_k562.ABC.Score", "K562_ATAC_only.ABC.Score",
     C["obs"], False),
]
fig, ax = plt.subplots(figsize=(7.4, 3.6))
for i, (lab, aa, bb, col, bold) in enumerate(FOREST):
    r = row(aa, bb)
    ax.plot([r.ci_lo, r.ci_hi], [i, i], color=col, lw=2.6 if bold else 1.8,
            solid_capstyle="round")
    ax.scatter([r.delta], [i], color=col, s=64 if bold else 40, zorder=4)
    ax.text(0.127, i, f"{r.delta:+.4f}   {r.sign_kept_pct:.1f}%", va="center", fontsize=7.8,
            color="#333", fontweight="bold" if bold else "normal")
ax.axvline(0, color="#444", lw=1)
ax.set_yticks(range(len(FOREST)))
ax.set_yticklabels([f[0] for f in FOREST], fontsize=8.2)
ax.invert_yaxis()
ax.set_xlim(-0.03, 0.178)
ax.set_xlabel("$\\Delta$ AUPRC, paired bootstrap (2,000 resamples)")
ax.text(0.127, -0.85, "delta   sign kept", fontsize=7.6, color="#666", fontweight="bold")
ax.set_title("Paired AUPRC differences, K562 CRISPR benchmark", loc="left")
save(fig, "fig21_multicell_auprc_forest")


# ------------------------------------------------- fig 5: upstream against downstream
# EVERY POINT IS COMPUTED, not read off a finding. Each row names an upstream per-fold table
# and a downstream paired-bootstrap table produced by 4.17, and the same intervention is
# looked up in both. Upstream uses top-quintile Pearson throughout, so the x axis is one
# metric rather than a mixture; the interval is a 95% t interval on the 5 paired fold
# differences. Downstream is delta AUPRC with its bootstrap interval.
from scipy import stats as _st

CROSS = [
    # label, upstream table, cfg_a, cfg_b, downstream table, arm_a, arm_b
    ("DNase instead of ATAC\nas model input",
     "dnaseswap_k562_per_fold.tsv", "DNase_input", "ATAC_input",
     "dnase_input_crispr_paired.tsv",
     "pred_dnase_atacacc.ABC.Score", "pred_k562_multimodal.ABC.Score"),
    ("Training to 100 epochs\ninstead of early stopping",
     "multitask_k562_per_fold.tsv", "epoch_matched_baseline", "atac_input",
     "multitask_crispr_paired_controlled.tsv",
     "pred_ep100_atacacc.ABC.Score", "pred_k562_multimodal.ABC.Score"),
    ("Correcting the p300\ntraining target",
     "p300_multicell_k562_per_fold.tsv", "GM_single_corrected", "GM_single_corrupted",
     "p300_corrected_paired.tsv",
     "p300pred_gm12878_corrected.ABC.Score", "p300pred_gm12878_multimodal.ABC.Score"),
    ("Single-source transfer\nagainst the floor",
     "p300_multicell_k562_per_fold.tsv", "GM_single_corrected_tgtnorm",
     "K562_atac_only_floor",
     "2026_0922_multicell_auprc_pairs.tsv",
     "p300pred_gm12878_corrected_tgtnorm.ABC.Score", "K562_ATAC_only.ABC.Score"),
    ("Four source cell types\ninstead of one",
     "p300_multicell_k562_per_fold.tsv", "multicell_holdK562",
     "GM_single_corrected_tgtnorm",
     "2026_0922_multicell_auprc_pairs.tsv",
     "p300pred_multicell_holdK562.ABC.Score",
     "p300pred_gm12878_corrected_tgtnorm.ABC.Score"),
]


def upstream_delta(table, cfg_a, cfg_b, metric="overall_pearson_topq"):
    t = pd.read_csv(f"{P}/results/{table}", sep="\t")
    t = t[t.fold != "config"]
    t["fold"] = t["fold"].astype(int)
    va = t[t.config == cfg_a].sort_values("fold")[metric].astype(float).values
    vb = t[t.config == cfg_b].sort_values("fold")[metric].astype(float).values
    assert len(va) == len(vb) == 5, f"{table} {cfg_a} {cfg_b}: {len(va)} vs {len(vb)}"
    d = va - vb
    se = d.std(ddof=1) / np.sqrt(len(d))
    h = _st.t.ppf(0.975, len(d) - 1) * se
    return d.mean(), d.mean() - h, d.mean() + h


def downstream_delta(table, arm_a, arm_b):
    t = pd.read_csv(f"{P}/results/{table}", sep="\t")
    r = t[(t.arm_a == arm_a) & (t.arm_b == arm_b)]
    assert len(r) == 1, f"{table}: {arm_a} - {arm_b} matched {len(r)}"
    r = r.iloc[0]
    return float(r.delta), float(r.ci_lo), float(r.ci_hi)


fig, ax = plt.subplots(figsize=(6.8, 5.0))
print("\ncross-experiment comparison (fig 5)")
rows = []
for lab, ut, ca, cb, dt, aa, bb in CROSS:
    u, ulo, uhi = upstream_delta(ut, ca, cb)
    d, dlo, dhi = downstream_delta(dt, aa, bb)
    # A row is marked as a mismatch when one evaluation's interval excludes zero and the
    # other's does not, or when the point estimates have opposite signs.
    u_sig = (ulo > 0) or (uhi < 0)
    d_sig = (dlo > 0) or (dhi < 0)
    mismatch = (u_sig != d_sig) or (np.sign(u) != np.sign(d))
    rows.append((lab, u, d, mismatch))
    print(f"  {lab.replace(chr(10), ' '):48s} up {u:+.4f} [{ulo:+.4f},{uhi:+.4f}]  "
          f"down {d:+.4f} [{dlo:+.4f},{dhi:+.4f}]  mismatch={mismatch}")
    col = "#cb181d" if mismatch else "#2171b5"
    ax.errorbar(u, d, yerr=[[d - dlo], [dhi - d]], xerr=[[u - ulo], [uhi - u]],
                fmt="o", color=col, ms=7, capsize=2.5, lw=1.2, zorder=3)

# Hand-placed so the two points near (0.01, 0.015) do not overprint each other.
OFFSETS = {0: (8, 6), 1: (-118, -22), 2: (10, -26), 3: (8, 10), 4: (8, 8)}
for i, (lab, u, d, mismatch) in enumerate(rows):
    ax.annotate(lab, xy=(u, d), xytext=OFFSETS[i], textcoords="offset points", fontsize=7.4,
                color="#333")

ax.axhline(0, color="#888", lw=0.9)
ax.axvline(0, color="#888", lw=0.9)
xs = [r[1] for r in rows]
ys = [r[2] for r in rows]
ax.set_xlim(min(xs) - 0.035, max(xs) + 0.045)
ax.set_ylim(min(ys) - 0.012, max(ys) + 0.016)
ax.set_xlabel("upstream: $\\Delta$ top-quintile Pearson (95% t on 5 folds)")
ax.set_ylabel("downstream: $\\Delta$ AUPRC (95% bootstrap)")
ax.set_title("Upstream correlation difference against downstream AUPRC difference",
             loc="left", fontsize=9)
ax.text(0.02, 0.965, "red: the two evaluations differ in sign, or one excludes zero "
                     "and the other does not",
        transform=ax.transAxes, fontsize=7.2, color="#cb181d", va="top")
r_p = np.corrcoef(xs, ys)[0, 1]
ax.text(0.02, 0.915, f"Pearson r = {r_p:+.2f} across the five interventions (n=5)",
        transform=ax.transAxes, fontsize=7.2, color="#666", va="top")
save(fig, "fig22_upstream_vs_downstream")


# ------------------------------------------------- fig 6: what normalisation bought
fig, ax = plt.subplots(figsize=(5.8, 2.6))
up = pf[pf.config == "GM_single_corrected_tgtnorm"]["overall_pearson_topq"].values - \
     pf[pf.config == "GM_single_corrected"]["overall_pearson_topq"].values
r = row("p300pred_gm12878_corrected_tgtnorm.ABC.Score", "p300pred_gm12878_corrected.ABC.Score")
ax.errorbar([up.mean()], [1], xerr=[[up.mean() - up.min()], [up.max() - up.mean()]],
            fmt="o", color=C["tgtnorm"], ms=7, capsize=3)
ax.errorbar([r.delta], [0], xerr=[[r.delta - r.ci_lo], [r.ci_hi - r.delta]], fmt="o",
            color=C["own"], ms=7, capsize=3)
ax.axvline(0, color="#444", lw=1)
ax.set_yticks([0, 1])
ax.set_yticklabels(["downstream\n$\\Delta$ AUPRC", "upstream\n$\\Delta$ top-q Pearson"],
                   fontsize=8)
ax.set_ylim(-0.6, 1.6)
ax.set_xlabel("target-derived normalisation $-$ own training statistics")
ax.set_title("Effect of the accessibility normalisation policy", loc="left")
save(fig, "fig23_normalisation_policy")


# ------------------------------------------------- fig 7: the ladder
LADDER = [("observed p300", "p300obs_k562.ABC.Score", C["obs"]),
          ("predicted p300, K562 in-cell", "p300pred_k562_multimodal.ABC.Score", C["incell"]),
          ("single-source transfer", "p300pred_gm12878_corrected_tgtnorm.ABC.Score",
           C["tgtnorm"]),
          ("multi-cell transfer", "p300pred_multicell_holdK562.ABC.Score", C["multicell"]),
          ("ATAC only (floor)", "K562_ATAC_only.ABC.Score", C["floor"])]
auprc = {}
for _, r_ in pairs.iterrows():
    auprc[r_.arm_a] = r_.auprc_a
    auprc[r_.arm_b] = r_.auprc_b
floor = auprc["K562_ATAC_only.ABC.Score"]

fig, ax = plt.subplots(figsize=(6.6, 3.0))
for i, (lab, uid, col) in enumerate(LADDER):
    ax.barh(i, auprc[uid] - floor, left=floor, color=col, height=0.6)
    ax.text(auprc[uid] + 0.0015, i, f"{auprc[uid]:.4f}", va="center", fontsize=7.8)
ax.axvline(floor, color="#888", lw=1, ls="--")
ax.set_yticks(range(len(LADDER)))
ax.set_yticklabels([l for l, _, _ in LADDER], fontsize=8.2)
ax.invert_yaxis()
ax.set_xlim(floor - 0.004, auprc["p300obs_k562.ABC.Score"] + 0.012)
ax.set_xlabel("AUPRC on the shared pair set")
ax.set_title(f"AUPRC from the ATAC-only floor ({floor:.4f}) to observed p300 "
             f"({auprc['p300obs_k562.ABC.Score']:.4f})", loc="left", fontsize=9)
save(fig, "fig24_multicell_ladder")

print("done")
