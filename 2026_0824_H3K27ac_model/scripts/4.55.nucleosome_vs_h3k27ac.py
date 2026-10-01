#!/usr/bin/env python3
"""Does the nucleosome signature predict H3K27ac beyond ATAC the way it predicts p300?

WHY. 4.51 found ATAC fragment-length features add partial rho 0.06-0.16 for p300 beyond ATAC,
but the signs flip in MCF-7, so the signal looked library-dependent. H3K27ac sits on the
nucleosomes that flank the open centre, so a real nucleosome effect should be at least as strong
for H3K27ac as for p300 and should have the same sign in each cell. K562 and GM12878 are the two
cells with H3K27ac ChIP-seq in this project.

DATA. Fragment features: 4.51's per-cell tables ($SCRATCH/nucsig). Targets, both 5' ChIP-seq
bigwigs (plus + minus strands), counted two ways:
    p300      over the candidate region (as 4.43/4.51)
    k27_reg   H3K27ac over the candidate region
    k27_1kb   H3K27ac within +/- 1 kb of the region centre (covers the flanking nucleosomes)
Same regions, controls and residualisation as 4.51 (non-promoter, >= 10 fragments; ATAC rank and
rank^2, GC, CpG o/e, log distance to TSS). Joint score: OLS of the target residual on all
features, fitted in the other cell (deployable) and in-cell (optimistic).

Usage: 4.55.nucleosome_vs_h3k27ac.py [--out-prefix results/2026_0930_nucleosome_vs_h3k27ac]
"""
import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
CACHE = f"{A}/data/tf_activity_recovery"
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
OUT = f"{os.environ.get('SCRATCH', '/tmp')}/nucsig"
K27 = {"K562": f"{A}/data/h3k27ac_5p_{{}}.bw", "GM12878": f"{A}/data/gm12878_h3k27ac_5p_{{}}.bw"}
CELLS = list(K27)
KEY = ["chr", "start", "end"]
FEATS = ["nfr", "mono_flank", "mono_centre", "di_flank", "nfr_width", "positioning",
         "flank_ratio", "centre_occ", "short_frac"]
TARGETS = ["p300", "k27_reg", "k27_1kb"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_0930_nucleosome_vs_h3k27ac")
a = ap.parse_args()


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


def bw_sum(path, ch, st, en):
    bw = pyBigWig.open(path)
    L = bw.chroms()
    out = np.zeros(len(ch))
    for i, (c, s, e) in enumerate(zip(ch, st, en)):
        if c in L:
            v = bw.values(c, max(0, int(s)), min(int(e), L[c]), numpy=True)
            out[i] = np.nansum(np.abs(v))
    bw.close()
    return out


pct = lambda v: (rankdata(v) - 0.5) / len(v)


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


S = pd.read_csv(f"{CACHE}/union_signals.tsv.gz", sep="\t")
U = S[KEY].copy()
U["mid"] = (U["start"] + U["end"]) // 2
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
U["prom"] = overlaps(U, tss)
tss["m"] = (tss["start"] + tss["end"]) // 2
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)

DATA = {}
rows = []
RES = []
for c in CELLS:
    f = pd.read_csv(f"{OUT}/{c}.tsv.gz", sep="\t")
    g = U.loc[f["u"].values].reset_index(drop=True)
    keep = ~g["prom"].values & (f["n_total"].values >= 10)
    f, g = f[keep].reset_index(drop=True), g[keep].reset_index(drop=True)
    dist = np.full(len(g), 1e7)
    for ch, gg in g.groupby("chr"):
        t = np.sort(tss.loc[tss["chr"] == ch, "m"].values)
        if len(t) > 1:
            k = np.clip(np.searchsorted(t, gg["mid"].values), 1, len(t) - 1)
            dist[gg.index.values] = np.minimum(np.abs(gg["mid"].values - t[k - 1]),
                                               np.abs(gg["mid"].values - t[k]))
    gc, cpg = np.empty(len(g)), np.empty(len(g))
    for i, (ch, mid) in enumerate(zip(g["chr"].values, g["mid"].values)):
        s = str(fa[ch][max(0, mid - 250):mid + 250])
        n = max(sum(s.count(b) for b in "ACGT"), 1)
        nc, ng = s.count("C"), s.count("G")
        gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
    atac = pct(S.loc[f["u"].values, f"atac_{c}"].values)
    Z = np.column_stack([np.ones(len(g)), atac, atac ** 2, gc, cpg, np.log10(dist + 1)])
    Y = {"p300": S.loc[f["u"].values, f"p300_{c}"].values}
    Y["k27_reg"] = sum(bw_sum(K27[c].format(s), g["chr"], g["start"], g["end"])
                       for s in ("plus", "minus"))
    Y["k27_1kb"] = sum(bw_sum(K27[c].format(s), g["chr"], g["mid"] - 1000, g["mid"] + 1000)
                       for s in ("plus", "minus"))
    X = []
    for k in FEATS:
        v = f[k].values.astype(float)
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        X.append(resid(pct(v), Z))
    X = np.column_stack(X)
    YR = {t: resid(pct(y), Z) for t, y in Y.items()}
    raw = {t: np.corrcoef(pct(y), atac)[0, 1] for t, y in Y.items()}
    for t in TARGETS:
        for j, k in enumerate(FEATS):
            rows.append(dict(cell=c, target=t, feature=k, partial_rho=np.corrcoef(X[:, j], YR[t])[0, 1]))
    DATA[c] = dict(X=X, YR=YR)
    RES.append(pd.DataFrame({"cell": c, "u": f["u"].values, "atac_pct": atac,
                             **{f"resid_{t}": YR[t] for t in TARGETS}}))
    print(f"{c}: {len(g):,} non-promoter regions; rho with ATAC: "
          + ", ".join(f"{t} {raw[t]:.3f}" for t in TARGETS), flush=True)
    print(f"  rho(p300 resid, k27_1kb resid) beyond ATAC+composition: "
          f"{np.corrcoef(YR['p300'], YR['k27_1kb'])[0, 1]:+.3f}", flush=True)

PR = pd.DataFrame(rows).pivot_table(index="feature", columns=["target", "cell"],
                                    values="partial_rho").reindex(FEATS)
pd.set_option("display.width", 200)
print("\n=== partial Spearman with each target beyond ATAC and composition ===")
print(PR.round(3).to_string())
print("\nsign agreement K562 vs GM12878 (of 9 features):",
      ", ".join(f"{t} {int((np.sign(PR[(t, 'K562')]) == np.sign(PR[(t, 'GM12878')])).sum())}"
                for t in TARGETS))

print("\n=== joint nucleosome score, partial rho with target ===")
J = []
for t in TARGETS:
    for c in CELLS:
        o = [d for d in CELLS if d != c][0]
        w = np.linalg.lstsq(np.column_stack([np.ones(len(DATA[o]["X"])), DATA[o]["X"]]),
                            DATA[o]["YR"][t], rcond=None)[0]
        xfer = np.corrcoef(DATA[c]["X"] @ w[1:], DATA[c]["YR"][t])[0, 1]
        win = np.linalg.lstsq(DATA[c]["X"], DATA[c]["YR"][t], rcond=None)[0]
        incell = np.corrcoef(DATA[c]["X"] @ win, DATA[c]["YR"][t])[0, 1]
        J.append(dict(target=t, cell=c, other_cell_weights=xfer, in_cell_weights=incell))
J = pd.DataFrame(J)
print(J.round(3).to_string(index=False))

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
PR.to_csv(f"{a.out_prefix}.partial_rho.tsv", sep="\t")
pd.concat(RES).to_csv(f"{a.out_prefix}.residuals.tsv.gz", sep="\t", index=False)
J.to_csv(f"{a.out_prefix}.joint.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.partial_rho.tsv, .joint.tsv, .residuals.tsv.gz")
