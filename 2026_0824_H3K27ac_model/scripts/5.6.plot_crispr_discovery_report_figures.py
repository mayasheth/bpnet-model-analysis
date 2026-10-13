#!/usr/bin/env python3
"""Report 7, Part C figures: features selected against CRISPR for scE2G (F-057 to F-061, 8.10).

    fig46  the p300 route: predicted p300 (F-057) and features chosen against p300 (F-058) in scE2G
    fig47  nested CRISPR feature discovery (F-059): a  held-out gain per procedure vs its null
                                                    b  picks per outer fold (forward selection)
                                                    c  official scE2G run
    fig48  CGI distance vs activity beyond ATAC in five cells (F-060): a  partial rho  b  deciles
    fig49  how the CGI-distance gain arises (F-061, 8.10): a  encodings  b  gene-density control
    fig50  eQTL transfer test and positive controls (F-063 to F-065): a  enrichment-recall  b, c  readouts vs n6

Inputs are the A2 summary tables (`2026_1002_sce2g_tf_activity/results/`).

Usage: 5.6.plot_crispr_discovery_report_figures.py [--outdir figures]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity/results"
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
BLUE, LBLUE, RED, GREY, GREEN, OCHRE, PURPLE, LPURPLE = (
    "#00488d", "#5496ce", "#9b241c", "#6e788d", "#429130", "#9b740a", "#792374", "#b778b3")
CELLC = {"K562": RED, "GM12878": BLUE, "HepG2": GREEN, "A549": OCHRE, "MCF-7": PURPLE}
B = "multiome_powerlaw_n6"


def tsv(f):
    return pd.read_csv(f"{A2}/{f}", sep="\t")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {name}")


def ci_dot(ax, x, y, lo, hi, color, label=None, marker="o", ms=5):
    ax.errorbar(x, y, xerr=[[x - lo], [hi - x]], fmt=marker, ms=ms, color=color, lw=1.1, capsize=0, label=label)


def forest(ax, rows, xlabel, title):
    y = np.arange(len(rows))[::-1]
    for yi, (_, est, lo, hi, col) in zip(y, rows):
        ci_dot(ax, est, yi, lo, hi, col)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")


# ---------------------------------------------------------------- fig46 the p300 route
pp = tsv("2026_1005_predp300_paired_bootstrap.tsv").set_index("comparison")
sc = tsv("2026_1005_screen_paired_bootstrap.tsv").set_index("comparison")
ac = tsv("2026_1004_activity_paired_bootstrap.tsv").set_index("comparison")
r = lambda t, m: (t.loc[f"{m} - {B}", "estimate"], t.loc[f"{m} - {B}", "ci_lo"], t.loc[f"{m} - {B}", "ci_hi"])  # noqa: E731
rows = [("+ predicted p300, multi-cell (deployable)", *r(pp, "multiome_n6_plus_predp300mc"), BLUE),
        ("+ predicted p300, in-cell K562 model", *r(pp, "multiome_n6_plus_predp300ic"), LBLUE),
        ("+ ENCODE K562 ATAC (control)", *r(pp, "multiome_n6_plus_atacenc"), GREY),
        ("ARC-E2G with sqrt(ATAC x predicted p300)", *r(pp, "multiome_n6_arcpredp300geo"), BLUE),
        ("+ GC (F-055 reference)", *r(sc, "multiome_n6_gc"), GREEN),
        ("+ GC + short-fragment share", *r(sc, "multiome_n6_gc_fragshort"), GREEN),
        ("+ GC + p300 composite", *r(sc, "multiome_n6_gc_composite"), GREEN),
        ("+ observed p300 (ceiling)", *r(ac, "multiome_n6_plus_p300"), OCHRE)]
fig, ax = plt.subplots(figsize=(6.4, 3.8))
forest(ax, rows, f"CRISPR AUPRC minus n6 ({pp.loc[B, 'estimate']:.3f}), 95% CI", "scE2G arms built from p300")
ax.text(0, -0.2, "blue: predicted p300 (F-057)   green: features chosen against p300 (F-058)\n"
        "grey: second-ATAC control   ochre: observed p300, not deployable (F-056)",
        transform=ax.transAxes, ha="left", va="top", fontsize=7)
save(fig, "fig46_p300_route")

# ---------------------------------------------------------------- fig47 nested discovery
nl = tsv("2026_1006_nested_null.tsv").set_index("model")
ch = tsv("2026_1006_nested_real.chosen.tsv")
db = tsv("2026_1007_discovery_paired_bootstrap.tsv").set_index("comparison")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.2, 0.9, 1.1], "wspace": 0.9})
order = [("forward selection (<= 5)", "forward"), ("L1 over all candidates", "l1"),
         ("family: CpG island", "family:cpg_island"), ("family: composition", "family:composition"),
         ("family: fragments", "family:fragments"), ("family: shape", "family:shape"),
         ("family: repeats", "family:repeats"), ("family: footprints", "family:footprints"),
         ("family: p300 composite", "family:p300_composite"), ("family: conservation", "family:conservation"),
         ("family: cells", "family:cells"), ("family: accessibility", "family:accessibility"),
         ("family: neighbourhood", "family:neighbourhood"), ("family: all motifs", "family:motifs")]
y = np.arange(len(order))[::-1]
for yi, (_, k) in zip(y, order):
    q = nl.loc[k]
    ci_dot(ax[0], q["real_gain"], yi, q["real_lo"], q["real_hi"], PURPLE, "held-out gain" if yi == y[0] else None)
    ax[0].plot([q["null_q95"]], [yi], "|", color=GREY, ms=9, mew=1.5, label="null 95th percentile" if yi == y[0] else None)
    ax[0].plot([q["null_max"]], [yi], "x", color=GREY, ms=5, label="null max (100 runs)" if yi == y[0] else None)
ax[0].axvline(0, color="black", lw=0.7)
ax[0].set_yticks(y)
ax[0].set_yticklabels([o[0] for o in order])
ax[0].set_xlabel("AUPRC gain over the split base (0.627), 95% CI")
ax[0].set_title("a  Nested selection, held-out chromosomes", loc="left")
ax[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.13), ncol=3)

fw = ch[ch["procedure"] == "forward"]["chosen"].value_counts()
lab = {"cgi_log_dist": "log distance to nearest CGI", "motif:GATA6_MA1104.3+2": "GATA cluster (GATA6/Gata3/TRPS1)",
       "motif:IKZF2_MA2326.1": "IKZF2 motif", "ins_centre_frac": "insertion centre fraction"}
fw = fw[fw >= 1].sort_values()
ax[1].barh(range(len(fw)), fw.values, color=[PURPLE if v >= 4 else LPURPLE for v in fw.values])
ax[1].set_yticks(range(len(fw)))
ax[1].set_yticklabels([lab.get(i, i.replace("motif:", "")) for i in fw.index], fontsize=7)
ax[1].axvline(4, color=GREY, lw=0.8, ls="--")
ax[1].set_xlim(0, 5.3)
ax[1].set_xlabel("Outer folds that picked it (of 5)")
ax[1].set_title("b  Forward-selection picks", loc="left")

rows = [("+ CGI distance", "multiome_n6_cgi", B, LPURPLE), ("+ GC (F-055)", "multiome_n6_gc", B, GREEN),
        ("+ the three (primary, pre-registered)", "multiome_n6_disc3", B, PURPLE),
        ("+ GC + the three", "multiome_n6_gc_disc3", B, PURPLE),
        ("the three, vs n6 + GC", "multiome_n6_disc3", "multiome_n6_gc", GREY),
        ("GC + the three, vs n6 + GC", "multiome_n6_gc_disc3", "multiome_n6_gc", GREY)]
rows = [(l_, *db.loc[f"{m} - {ref}", ["estimate", "ci_lo", "ci_hi"]].values, c) for l_, m, ref, c in rows]
forest(ax[2], rows, f"CRISPR AUPRC minus n6 ({db.loc[B, 'estimate']:.3f}), 95% CI", "c  Official scE2G run")
ax[2].axvline(0.021, color=PURPLE, lw=0.8, ls=":")
ax[2].text(0.021, -0.6, " nested prediction", color=PURPLE, fontsize=7, va="center")
save(fig, "fig47_crispr_discovery")

# ---------------------------------------------------------------- fig48 CGI distance x activity, five cells
pc = tsv("2026_1009_cgi_activity_panel.partial.tsv")
de = tsv("2026_1009_cgi_activity_panel.deciles.tsv")
cells = ["K562", "GM12878", "HepG2", "A549", "MCF-7"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"wspace": 0.45})
p3 = pc[pc["mark"] == "p300"].set_index("cell")
k27 = pc[(pc["mark"] == "H3K27ac") & (pc["window"] == "500bp")].set_index("cell")
y = np.arange(len(cells))[::-1]
for yi, c in zip(y, cells):
    q = p3.loc[c]
    ci_dot(ax[0], q["partial_Z1"], yi + 0.2, q["lo_Z1"], q["hi_Z1"], BLUE, "p300, ATAC + TSS distance removed" if c == "K562" else None)
    ci_dot(ax[0], q["partial_Z2"], yi, q["lo_Z2"], q["hi_Z2"], LBLUE, "p300, also GC + CpG o/e removed" if c == "K562" else None)
    if c in k27.index:
        q = k27.loc[c]
        ci_dot(ax[0], q["partial_Z1"], yi - 0.2, q["lo_Z1"], q["hi_Z1"], RED, "H3K27ac, ATAC + TSS distance removed" if c == "K562" else None, marker="s", ms=4)
ax[0].axvline(0, color="black", lw=0.7)
ax[0].axvline(p3.loc["K562", "partial_Z1"] / 2, color=GREY, lw=0.8, ls="--")
ax[0].set_yticks(y)
ax[0].set_yticklabels(cells)
ax[0].set_xlabel("Partial Spearman of log CGI distance with the mark (95% CI)")
ax[0].set_title("a  Farther from a CGI, more activity than ATAC predicts", loc="left")
ax[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.15), ncol=2)
d = de[de["mark"] == "p300"].set_index("cell")
xs = np.arange(1, 11)
for c in cells:
    ax[1].plot(xs, d.loc[c, [f"d{i}" for i in xs]].values.astype(float), "-o", ms=3, color=CELLC[c], lw=1.2, label=c)
ax[1].axhline(0, color="black", lw=0.7)
ax[1].set_xticks(xs)
ax[1].set_xlabel("CGI-distance decile (1 = nearest)")
ax[1].set_ylabel("Mean p300 residual given ATAC + TSS distance\n(percentile units)")
ax[1].set_title("b  Where along the distance range", loc="left")
ax[1].legend(frameon=False, fontsize=7, loc="lower right")
save(fig, "fig48_cgi_activity_panel")

# ---------------------------------------------------------------- fig49 encodings and gene density
en = tsv("2026_1009_near_cgi_encoding.tsv").set_index("model")
gd = tsv("2026_1009_gene_density_control.tsv").set_index(["model", "ref"])
fig, ax = plt.subplots(1, 2, figsize=(13, 3.9), gridspec_kw={"wspace": 1.0})
rows = [("percentile (F-059 encoding)", "n6 + pct", PURPLE), ("nearest-10% flag (primary)", "n6 + near10", BLUE),
        ("nearest-20% flag", "n6 + near20", LBLUE), ("graded within nearest 10%", "n6 + cap10", LBLUE),
        ("graded within nearest 20%", "n6 + cap20", LBLUE)]
forest(ax[0], [(l_, en.loc[m, "vs_n6"], en.loc[m, "lo_n6"], en.loc[m, "hi_n6"], c) for l_, m, c in rows],
       "CRISPR AUPRC minus n6 (replica), 95% CI", "a  Encoding CGI distance")
ax[0].axvline(0.75 * en.loc["n6 + pct", "vs_n6"], color=GREY, lw=0.8, ls="--")
rows = [("+ gene density G", ("n6 + G", "n6"), GREY), ("+ CGI distance", ("n6 + pct", "n6"), PURPLE),
        ("+ CGI distance, on top of G", ("n6 + G + pct", "n6 + G"), PURPLE),
        ("+ the three", ("n6 + 3", "n6"), PURPLE), ("+ the three, on top of G", ("n6 + G + 3", "n6 + G"), PURPLE),
        ("+ G, on top of the three", ("n6 + G + 3", "n6 + 3"), GREY)]
forest(ax[1], [(l_, gd.loc[k, "diff"], gd.loc[k, "lo"], gd.loc[k, "hi"], c) for l_, k, c in rows],
       "CRISPR AUPRC difference (replica), 95% CI", "b  Gene-density control")
save(fig, "fig49_cgi_encoding_gene_density")

# ---------------------------------------------------------------- fig50 eQTL transfer test and positive controls
ro = tsv("2026_1013_eqtl_activity_control.readouts.tsv")
cv = tsv("2026_1013_eqtl_activity_control.primary_curves.tsv")
PRIM = "primary (all matches except K562)"
ARMC = {"n6": ("n6 (baseline)", "black"), "n6_three": ("+ the three (F-059)", PURPLE), "n6_cgi": ("+ CGI distance", LPURPLE),
        "n6_abc": ("ARC-E2G -> ABC (activity-only control)", OCHRE), "scatac_published": ("scATAC-only scE2G (control)", GREY)}
fig, ax = plt.subplots(1, 3, figsize=(17, 4.4), gridspec_kw={"width_ratios": [1, 0.8, 0.8], "wspace": 0.95})
for k, (lab_, col) in ARMC.items():
    r_, e_ = cv[f"{k}_recall"].values, cv[f"{k}_enrichment"].values
    ok = (r_ >= 0.01) & np.isfinite(e_)
    ax[0].plot(r_[ok], e_[ok], color=col, lw=1.3, label=lab_)
ax[0].axvline(0.10, color=GREY, lw=0.7, ls=":")
ax[0].axvspan(0.02, 0.20, color=GREY, alpha=0.08, lw=0)
ax[0].set_xlabel("Linking recall (eQTL-eGene pairs recovered)")
ax[0].set_ylabel("eQTL enrichment over background SNPs")
ax[0].set_ylim(0, 50)
ax[0].set_title("a  Enrichment-recall, 62 matched immune pairs", loc="left")
ax[0].legend(frameon=False, fontsize=7, loc="upper right")
p_ = ro[ro["group"] == PRIM].set_index("arm")
order = ["n6_three", "n6_cgi", "n6_abc", "scatac_published"]
for j, (rn, tt) in enumerate([("enr_at_r10", "b  Enrichment at recall 0.10, vs n6"), ("enr_mean", "c  Mean enrichment, recall 0.02-0.20, vs n6")]):
    rows = [(ARMC[k][0], p_.loc[k, f"{rn}_vs_n6"], p_.loc[k, f"{rn}_vs_n6_lo"], p_.loc[k, f"{rn}_vs_n6_hi"], ARMC[k][1]) for k in order]
    forest(ax[1 + j], rows, f"Difference from n6 ({p_.loc['n6', rn]:.1f}), 95% CI", tt)
    if j:
        ax[2].set_yticklabels([])
        ax[2].set_ylabel("")
save(fig, "fig50_eqtl_transfer")
