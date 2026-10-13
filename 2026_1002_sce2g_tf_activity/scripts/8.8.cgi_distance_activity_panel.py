#!/usr/bin/env python3
"""Does distance to the nearest CpG island track activity beyond ATAC, in K562 and in other cells?

WHY. F-059: log distance to the nearest CpG island is the cell-independent part of the three
features that raise scE2G K562 CRISPR AUPRC; alone it adds +0.024 [+0.014, +0.033], and its
coefficient is positive (+3.18 on top of ARC-E2G in the official run): elements FAR from CpG
islands are regulated more often than ATAC-based ARC predicts. F-007 and F-055 point the same way
(ATAC over-ranks GC- and CpG-rich elements). There is no CRISPR set in a second cell, so this asks
the label-free version (Maya, 2026-10-08): if the CRISPR gain comes from ATAC over-ranking
CGI-proximal elements relative to their activity, CGI distance should correlate positively with
p300 and H3K27ac beyond ATAC in K562, and a similar correlation in the other panel cells would
mean the correction is not K562-specific.

REGIONS AND SIGNALS. 4.60's union set and caches (`2026_0824_H3K27ac_model`): ABC candidates from
each cell's own ATAC, merged (356,587 regions); a region is in cell c if it overlaps c's own
candidates; promoters (GENCODE TSS +/- 250 bp) dropped. ATAC and p300 (+/- 250 bp) from 4.43's
cache; H3K27ac (K562, GM12878) and H3K4me1 (all five) at +/- 250 bp and +/- 1 kb from 4.60's cache.
CGI distance exactly as 7.1 / 8.2: log10(1 + bp from the region midpoint to the nearest UCSC
cpgIslandExt midpoint on the same chromosome).

MEASURE. Partial Spearman of CGI distance with each mark, per cell over its own non-promoter
regions: both percentile-ranked, residualised on Z, then correlated.
    Z1 (primary): ATAC percentile, its square, log distance to TSS. GC is left out on purpose:
       it is part of the same signal (F-059: GC adds nothing on top of CGI distance).
    Z2 (descriptive): Z1 + GC + CpG o/e (+/- 250 bp), 4.60's controls; how much of the
       correlation is local composition rather than distance to an island.
Intervals: 1 Mb block bootstrap (1,000 resamples) over the residual pairs; residualisation is
not refitted per resample, so intervals are slightly narrow.
Shape: mean Z1 residual of each mark (in percentile units) by CGI-distance decile.

READING, fixed before the run:
- K562 check: p300 under Z1 should be positive. If its interval includes 0 or it is below 0, the
  CRISPR gain does not run through p300 beyond ATAC, and this proxy cannot speak to other cells
  (inconclusive, not negative). H3K27ac in K562 is read the same way.
- Transfer: if the K562 check passes, the effect counts as shared when p300 under Z1 has its
  interval above 0 in >= 3 of the 4 other cells (GM12878, HepG2, A549, MCF-7), each at least half
  the K562 point estimate. Otherwise call it K562-specific.
- H3K27ac in GM12878 and H3K4me1 in all five are secondary, read the same way.
- Even a shared effect shows only that the mechanism is shared; F-058 showed p300 patterns need
  not carry to CRISPR. It supports, but does not replace, a CRISPR test in a second cell.

Usage: 8.8.cgi_distance_activity_panel.py [--out-prefix results/2026_1009_cgi_activity_panel]
"""
import argparse
import os

import numpy as np
import pandas as pd
import pyfaidx
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
H = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CANDS = {"K562": f"{ABC}/2026_0905_p300_activity/p300obs_k562",
         **{c: f"{ABC}/2026_0924_p300_panel_obs/p300obs_{c.replace('-', '_')}"
            for c in ("GM12878", "HepG2", "A549", "MCF-7")}}
CACHE = f"{H}/data/tf_activity_recovery/union_signals.tsv.gz"
K4CACHE = f"{H}/data/h3k4me1_panel/union_h3k4me1.tsv.gz"
CGI = f"{os.environ['SCRATCH']}/ucsc_hg38/cpgIslandExt.txt.gz"
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
KEY = ["chr", "start", "end"]
CELLS = list(CANDS)
NBOOT = 1000

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_1009_cgi_activity_panel")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)


def overlaps(reg, iv):
    hit = np.zeros(len(reg), bool)
    for c, g in reg.groupby("chr"):
        e = iv[iv["chr"] == c].sort_values("start")
        if e.empty:
            continue
        es, ee = e["start"].values, np.maximum.accumulate(e["end"].values)
        k = np.searchsorted(es, g["end"].values, side="left") - 1
        ok = k >= 0
        hit[g.index.values[ok]] = ee[k[ok]] > g["start"].values[ok]
    return hit


def nearest(reg_mid, chrs, pts):
    """Distance from each midpoint to the nearest point on the same chromosome (1e7 if none)."""
    d = np.full(len(reg_mid), 1e7)
    for c in np.unique(chrs):
        t = np.sort(pts.get(c, np.array([])))
        i = np.where(chrs == c)[0]
        if len(t) > 1:
            k = np.clip(np.searchsorted(t, reg_mid[i]), 1, len(t) - 1)
            d[i] = np.minimum(np.abs(reg_mid[i] - t[k - 1]), np.abs(reg_mid[i] - t[k]))
    return d


# ---------------------------------------------------------------- union set, as 4.43 / 4.60
cand = {c: pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t", usecols=KEY)
        for c, d in CANDS.items()}
allr = pd.concat(cand.values()).sort_values(KEY).reset_index(drop=True)
merged = []
for c, g in allr.groupby("chr", sort=True):
    s, e = g["start"].values, g["end"].values
    cs, ce = s[0], e[0]
    for x, y in zip(s[1:], e[1:]):
        if x <= ce:
            ce = max(ce, y)
        else:
            merged.append((c, cs, ce)); cs, ce = x, y
    merged.append((c, cs, ce))
U = pd.DataFrame(merged, columns=KEY)
U["mid"] = (U["start"] + U["end"]) // 2
S = pd.read_csv(CACHE, sep="\t")
K = pd.read_csv(K4CACHE, sep="\t")
for name, T in (("4.43 cache", S), ("4.60 cache", K)):
    assert len(T) == len(U) and (T[KEY].values == U[KEY].values).all(), f"union set differs from {name}"
for c in CELLS:
    U[f"in_{c}"] = overlaps(U, cand[c])
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
U["prom"] = overlaps(U, tss)
tss_m = {c: ((g["start"] + g["end"]) // 2).values for c, g in tss.groupby("chr")}
U["log_dist_tss"] = np.log10(nearest(U["mid"].values, U["chr"].values, tss_m) + 1)
G = pd.read_csv(CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])
cgi_m = {c: ((g["s"] + g["e"]) // 2).values for c, g in G.groupby("chr")}
U["cgi_log_dist"] = np.log10(nearest(U["mid"].values, U["chr"].values, cgi_m) + 1)
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
gc, cpg = np.empty(len(U)), np.empty(len(U))
for i, (ch, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = str(fa[ch][max(0, mid - 250):mid + 250])
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
U["gc"], U["cpg_oe"] = gc, cpg
U["block"] = U["chr"] + ":" + (U["mid"] // 1_000_000).astype(str)
assert U["cgi_log_dist"].between(0, 7.01).all()
print(f"union: {len(U):,} regions, {int((~U['prom']).sum()):,} non-promoter; "
      f"median CGI distance {10 ** U['cgi_log_dist'].median():,.0f} bp", flush=True)

# ---------------------------------------------------------------- measures
pct = lambda v: (rankdata(v) - 0.5) / len(v)  # noqa: E731
cor = lambda x, y: float(np.corrcoef(x, y)[0, 1])  # noqa: E731


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


def Zmat(sel, c, comp):
    at = pct(S.loc[sel, f"atac_{c}"].values)
    g = U.loc[sel]
    cols = [np.ones(len(at)), at, at ** 2, g["log_dist_tss"].values]
    if comp:
        cols += [g["gc"].values, g["cpg_oe"].values]
    return np.column_stack(cols)


def block_boot(x, y, blocks):
    """Percentile interval of corr(x, y) resampling 1 Mb blocks, from per-block sums."""
    codes, inv = np.unique(blocks, return_inverse=True)
    sums = np.stack([np.bincount(inv, w, len(codes))
                     for w in (np.ones_like(x), x, y, x * x, y * y, x * y)], 1)
    W = rng.multinomial(len(codes), np.full(len(codes), 1 / len(codes)), size=NBOOT)
    n, sx, sy, sxx, syy, sxy = (W @ sums).T
    r = (sxy - sx * sy / n) / np.sqrt((sxx - sx ** 2 / n) * (syy - sy ** 2 / n))
    return np.quantile(r, [0.025, 0.975])


nonprom = ~U["prom"].values
rows, shape = [], []
for c in CELLS:
    sel = U[f"in_{c}"].values & nonprom
    x = pct(U.loc[sel, "cgi_log_dist"].values)
    blocks = U.loc[sel, "block"].values
    dec = np.minimum((x * 10).astype(int), 9)
    marks = [("ATAC", None, S.loc[sel, f"atac_{c}"].values), ("p300", "500bp", S.loc[sel, f"p300_{c}"].values)]
    for w in ("500bp", "2kb"):
        if f"k27_{c}_{w}" in K:
            marks.append(("H3K27ac", w, K.loc[sel, f"k27_{c}_{w}"].values))
        marks.append(("H3K4me1", w, K.loc[sel, f"k4_{c}_{w}"].values))
    for mark, w, v in marks:
        y = pct(v)
        row = dict(cell=c, mark=mark, window=w or "500bp", n=int(sel.sum()), rho_raw=cor(x, y))
        if mark != "ATAC":
            for zname, comp in (("Z1", False), ("Z2", True)):
                Z = Zmat(sel, c, comp)
                rx, ry = resid(x, Z), resid(y, Z)
                lo, hi = block_boot(rx, ry, blocks)
                row.update({f"partial_{zname}": cor(rx, ry), f"lo_{zname}": lo, f"hi_{zname}": hi})
                if zname == "Z1":
                    m = pd.Series(ry).groupby(dec).mean()
                    shape.append(dict(cell=c, mark=mark, window=w, **{f"d{k + 1}": m.get(k, np.nan)
                                                                     for k in range(10)}))
        rows.append(row)
    print(f"{c}: {int(sel.sum()):,} non-promoter regions", flush=True)

PC = pd.DataFrame(rows)
SH = pd.DataFrame(shape)
pd.set_option("display.width", 220)
print("\n=== partial Spearman of log CGI distance with each mark (Z1 = ATAC + ATAC^2 + log TSS dist;"
      " Z2 = Z1 + GC + CpG o/e); 1 Mb block-bootstrap 95% intervals ===")
print(PC.round(3).to_string(index=False))
print("\n=== mean Z1 residual (percentile units) by CGI-distance decile, d1 = nearest ===")
print(SH.round(3).to_string(index=False))

# pre-registered reading, printed so the verdict is in the log
p = PC[(PC["mark"] == "p300")].set_index("cell")
k = p.loc["K562"]
print(f"\nK562 check, p300 Z1: {k.partial_Z1:+.3f} [{k.lo_Z1:+.3f}, {k.hi_Z1:+.3f}] -> "
      + ("PASS" if k.lo_Z1 > 0 else "FAIL (inconclusive for other cells)"))
if k.lo_Z1 > 0:
    ok = [(c, (p.loc[c, "lo_Z1"] > 0) and (p.loc[c, "partial_Z1"] >= k.partial_Z1 / 2))
          for c in CELLS if c != "K562"]
    n_ok = sum(v for _, v in ok)
    print("transfer: " + ", ".join(f"{c} {'yes' if v else 'no'}" for c, v in ok)
          + f" -> {n_ok}/4: " + ("SHARED" if n_ok >= 3 else "K562-SPECIFIC"))

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
PC.to_csv(f"{a.out_prefix}.partial.tsv", sep="\t", index=False)
SH.to_csv(f"{a.out_prefix}.deciles.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.partial.tsv, .deciles.tsv")
