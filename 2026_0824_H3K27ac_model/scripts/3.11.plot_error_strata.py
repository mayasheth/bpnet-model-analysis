#!/usr/bin/env python
"""Figure: what the mis-predicted elements are.

Three panels, because the finding has three parts and the tables buried all of them:

  a  WHERE the two error tails sit in (accessibility, H3K27ac) space. This is the whole
     phenotype in one panel: over-predicted elements are open and unacetylated, under-
     predicted ones are acetylated and comparatively closed. Everything else follows.
  b  WHICH marks are actually enriched, as fold change over the typical stratum, with IgG
     drawn as the background control. Log axis because the range spans 50x.
  c  WHY peak overlap and quantitative signal disagree on the over-predicted tail. Peaks
     say it is enriched for everything; signal says CTCF only. The disagreement is a
     finding, so it gets a panel rather than a footnote.
  d  WHETHER CTCF explains the over-predicted tail. It does not: CTCF-high elements are a
     quarter of the tail and the CTCF-low remainder is over-predicted more, not less. This
     is the panel that refutes the reading panel b invites, which is exactly why both are
     shown -- panel b is a POOLED per-stratum statistic and a minority of elements sets it.

All three read the SAME per-element stratum assignment written by 4.10 (--emit-elements),
so no panel can drift from the tables.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import (apply_rcparams, save_fig, add_panel_label, figsize,
                          n_label, annotate_n)

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
EL = f"{P}/results/error_strata_elements.tsv"
CTCF_EL = f"{P}/results/error_strata_ctcf_elements.tsv"
RPKM = f"{P}/results/error_strata_rpkm.tsv"
PEAK = f"{P}/results/error_strata_peak_overlap.tsv"

OVER, UNDER, TYP = "over-predicted 1%", "under-predicted 1%", "typical (middle 50%)"
# Panels b-d use the 5% tails, matching the tables in the report and panel d's split;
# panel a uses the 1% tails, where fewer points make the phenotype legible.
OVER5, UNDER5 = "over-predicted 5%", "under-predicted 5%"
C_OVER, C_UNDER, C_BG = "#d94801", "#2171b5", "#d9d9d9"
MARKS = ["CTCF", "EP300", "H3K4me1", "H3K27me3", "H3K27ac", "IgG"]

apply_rcparams()
d = pd.read_csv(EL, sep="\t")
rp = pd.read_csv(RPKM, sep="\t").set_index("stratum")
pk = pd.read_csv(PEAK, sep="\t").set_index("stratum")

# The BED-derived stratum labels use underscores and a 'pct' suffix; map to the element
# table's labels so a rename in one place cannot silently mismatch the other.
def row(tbl, lab):
    key = {OVER: "over-predicted_1pct", UNDER: "under-predicted_1pct",
           OVER5: "over-predicted_5pct", UNDER5: "under-predicted_5pct",
           TYP: "typical_middle_50pct"}[lab]
    if key not in tbl.index:
        raise SystemExit(f"{lab!r} -> {key!r} not in table; have {list(tbl.index)}")
    return tbl.loc[key]

ce = pd.read_csv(CTCF_EL, sep="\t")

fig, axes = plt.subplots(2, 2, figsize=figsize(columns=2, aspect=0.78))
axes = axes.ravel()
ax = axes[0]

# --- a: the phenotype, in (ATAC, H3K27ac) space -------------------------------
x, y = np.log10(d["ATAC.RPM"] + 0.1), np.log10(d["k27_obs"] + 0.1)
ax.hexbin(x, y, gridsize=55, bins="log", cmap="Greys", mincnt=1, linewidths=0)
for lab, col in ((UNDER, C_UNDER), (OVER, C_OVER)):
    m = d["stratum"] == lab
    ax.scatter(x[m], y[m], s=1.4, c=col, alpha=0.55, linewidths=0, rasterized=True,
               label=f"{lab.replace('-predicted', '-pred.')} (n = {m.sum():,})")
ax.set_xlabel("observed ATAC (log$_{10}$ RPM)")
ax.set_ylabel("observed H3K27ac (log$_{10}$ RPM)")
ax.legend(loc="upper left", frameon=False, fontsize=5, handletextpad=0.3,
          borderpad=0.1, labelspacing=0.25, markerscale=3)
annotate_n(ax, n_label(n=len(d), unit="candidate regions"), loc="lower right")
add_panel_label(ax, "a")

# --- b: fold enrichment over the typical stratum ------------------------------
ax = axes[1]
t = row(rp, TYP)
fo = [row(rp, OVER5)[m] / t[m] for m in MARKS]
fu = [row(rp, UNDER5)[m] / t[m] for m in MARKS]
ypos = np.arange(len(MARKS))[::-1]
for yy, a, b in zip(ypos, fo, fu):
    ax.plot([min(a, b), max(a, b)], [yy, yy], color="#cccccc", lw=0.8, zorder=1)
ax.scatter(fo, ypos, s=16, c=C_OVER, zorder=3, label="over-predicted 5%")
ax.scatter(fu, ypos, s=16, c=C_UNDER, zorder=3, label="under-predicted 5%")
ax.axvline(1.0, color="#444444", lw=0.7, ls=(0, (3, 2)), zorder=0)
ax.set_xscale("log")
ax.set_yticks(ypos)
ax.set_yticklabels([m + "  (control)" if m == "IgG" else m for m in MARKS])
for lbl in ax.get_yticklabels():
    if "control" in lbl.get_text():
        lbl.set_color("#888888")
ax.set_xlabel("RPKM, fold vs typical elements")
ax.legend(loc="lower right", frameon=False, fontsize=5, handletextpad=0.3, borderpad=0.1)
annotate_n(ax, "n = 7,678 per tail; 76,771 typical", loc="upper right")
add_panel_label(ax, "b")

# --- c: peaks say everything, signal says CTCF --------------------------------
ax = axes[2]
pm = [m for m in MARKS if m != "IgG"]        # no IgG peak calls exist
pe = [row(pk, OVER5)[m] / row(pk, TYP)[m] for m in pm]
se = [row(rp, OVER5)[m] / row(rp, TYP)[m] for m in pm]
xx = np.arange(len(pm))
w = 0.36
ax.bar(xx - w / 2, pe, w, color="#bdbdbd", edgecolor="none", label="peak overlap")
ax.bar(xx + w / 2, se, w, color=C_OVER, edgecolor="none", label="quantitative signal")
ax.axhline(1.0, color="#444444", lw=0.7, ls=(0, (3, 2)))
ax.set_xticks(xx)
ax.set_xticklabels(pm, rotation=35, ha="right")
ax.set_ylabel("over-predicted 5%, fold vs typical")
ax.legend(loc="upper right", frameon=False, fontsize=5, handletextpad=0.4, borderpad=0.1)
annotate_n(ax, "n = 7,678 vs 76,771", loc="upper left")
add_panel_label(ax, "c")

# --- d: does CTCF explain the tail? Split it and see. ------------------------
ax = axes[3]
thr = ce.loc[ce["stratum"] == TYP, "CTCF_rpkm"].quantile(0.90)
over5 = ce[ce["over5"] == 1]
hi = over5[over5["CTCF_rpkm"] > thr]
lo = over5[over5["CTCF_rpkm"] <= thr]
typ_e = ce[ce["stratum"] == TYP]
groups = [("typical", typ_e, C_BG), ("over-pred.,\nCTCF-low", lo, C_OVER),
          ("over-pred.,\nCTCF-high", hi, "#8c2d04")]
bp = ax.boxplot([g[1]["err"].to_numpy() for g in groups],
                positions=[0, 1, 2], widths=0.55, showfliers=False,
                medianprops=dict(color="#111111", lw=1.1), patch_artist=True)
for patch, g in zip(bp["boxes"], groups):
    patch.set_facecolor(g[2]); patch.set_edgecolor("#555555"); patch.set_linewidth(0.6)
ax.axhline(0, color="#444444", lw=0.7, ls=(0, (3, 2)))
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([g[0] for g in groups])
ax.set_ylabel("prediction error\n(residualised, >0 = over-predicts)")
for i, g in enumerate(groups[1:], start=1):
    ax.annotate(f"{100 * len(g[1]) / len(over5):.0f}% of tail\nmedian {g[1]['err'].median():.2f}",
                xy=(i, ax.get_ylim()[1]), xytext=(0, -2), textcoords="offset points",
                ha="center", va="top", fontsize=5, color="#333333")
annotate_n(ax, f"CTCF-high cut = {thr:.1f} RPKM (90th pct of typical); "
               f"n = {len(lo):,} / {len(hi):,} / {len(typ_e):,}", loc="lower left")
add_panel_label(ax, "d")

for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = save_fig(fig, f"{P}/figures/fig15_error_strata")
print(f"wrote {out}")

print("\npanel d: does removing CTCF-high rescue the tail?")
print(f"   threshold {thr:.2f} RPKM; CTCF-high {len(hi):,} ({100*len(hi)/len(over5):.1f}%) "
      f"median err {hi['err'].median():.3f}")
print(f"   CTCF-low {len(lo):,} ({100*len(lo)/len(over5):.1f}%) "
      f"median err {lo['err'].median():.3f}")
print(f"   typical {len(typ_e):,} median err {typ_e['err'].median():.3f}")

print("\nvalues plotted (5% tails, fold vs typical):")
print(f"{'mark':<10}{'over signal':>13}{'over peaks':>12}{'under signal':>14}")
for m in MARKS:
    ps = row(pk, OVER5)[m] / row(pk, TYP)[m] if m in pk.columns else np.nan
    print(f"{m:<10}{row(rp, OVER5)[m] / t[m]:>13.2f}{ps:>12.2f}"
          f"{row(rp, UNDER5)[m] / t[m]:>14.2f}")
