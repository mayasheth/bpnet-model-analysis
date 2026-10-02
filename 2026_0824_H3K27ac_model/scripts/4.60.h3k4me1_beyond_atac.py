#!/usr/bin/env python3
"""Does H3K4me1 carry information beyond ATAC, reproducibly and consistently across cells?

WHY. F-043 found the ATAC term in ABC's activity score costs observed p300 and H3K27ac AUPRC, and
H3K4me1 is the one enhancer mark with ChIP-seq in all five panel cells. Training an H3K4me1 model
is only worth it if (1) H3K4me1 is not just ATAC again, (2) what is left after ATAC is
reproducible between replicates rather than noise, (3) it tracks p300 beyond ATAC the way H3K27ac
does, and (4) it behaves the same way across cells, so a model trained on four cells has
something to transfer to the fifth. This measures all four before any training.

REGIONS. 4.43's union set (ABC candidates from each cell's own ATAC, merged; 356,587 regions),
rebuilt here for the per-cell membership flags and checked against its cached signal table. A
region is "in" cell c if it overlaps one of c's own candidates. Promoters (GENCODE TSS +/- 250 bp)
are dropped. ATAC and p300 come from 4.43's cache (+/- 250 bp of the region midpoint, counted
from the bigwigs the models train on). H3K4me1 is counted from 0.47's stranded 5' tracks
(plus + minus) in two windows: +/- 250 bp (same as ATAC and p300) and +/- 1 kb (H3K4me1 sits on
the flanking nucleosomes, as H3K27ac does in 4.55). H3K27ac from 0.14 in K562 and GM12878 is
counted the same way, as the reference mark.

CONTROLS Z, as in 4.43/4.55: ATAC percentile and its square, GC, CpG o/e, log distance to TSS.
All signals are percentile-ranked before residualising, so every number is a partial Spearman.

MEASURES, per cell over its own non-promoter regions:
    rho_atac      Spearman of the mark with ATAC
    r2_Z          variance of the mark's rank explained by Z; 1 - r2_Z is the part "beyond ATAC"
    rep_resid     mean over replicate pairs of corr(resid rep_i, resid rep_j): how much of the
                  beyond-ATAC part is reproducible. Upper bound for anything a model can learn
    p300_partial  corr(resid mark, resid p300): does the beyond-ATAC part track p300
Across cells:
    xcell_resid   for each pair (c, d), over regions in both, corr(resid mark in c, resid mark in
                  d), each residualised on its own cell's Z. Reported beside the same for p300
                  and for raw ATAC rank, which bracket it
    loo_p300      in cell c, partial corr of c's p300 residual with the MEAN H3K4me1 residual of
                  the other four cells (each residualised on its own Z, then on c's Z). A
                  deployable proxy: if H3K4me1 beyond ATAC is shared across cells, the other
                  cells' values should already predict c's p300. p300's own LOO is the reference

Usage: 4.60.h3k4me1_beyond_atac.py [--out-prefix results/2026_1001_h3k4me1_beyond_atac]
"""
import argparse
import itertools
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CANDS = {"K562": f"{ABC}/2026_0905_p300_activity/p300obs_k562",
         **{c: f"{ABC}/2026_0924_p300_panel_obs/p300obs_{c.replace('-', '_')}"
            for c in ("GM12878", "HepG2", "A549", "MCF-7")}}
CACHE = f"{A}/data/tf_activity_recovery/union_signals.tsv.gz"
K4DIR = f"{A}/data/h3k4me1_panel"
K4CACHE = f"{K4DIR}/union_h3k4me1.tsv.gz"
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
KEY = ["chr", "start", "end"]
CELLS = list(CANDS)
NREP = {"K562": 3, "GM12878": 2, "HepG2": 2, "A549": 3, "MCF-7": 2}
K27 = {"K562": f"{A}/data/h3k27ac_5p_{{}}.bw", "GM12878": f"{A}/data/gm12878_h3k27ac_5p_{{}}.bw"}
WIN = {"500bp": 250, "2kb": 1000}

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_1001_h3k4me1_beyond_atac")
ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
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


# ---------------------------------------------------------------- union set, as 4.43
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
assert len(S) == len(U) and (S[KEY].values == U[KEY].values).all(), "union set differs from 4.43"
for c in CELLS:
    U[f"in_{c}"] = overlaps(U, cand[c])
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
U["prom"] = overlaps(U, tss)
tss["m"] = (tss["start"] + tss["end"]) // 2
dist = np.full(len(U), 1e7)
for c, g in U.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == c, "m"].values)
    if len(t) > 1:
        k = np.clip(np.searchsorted(t, g["mid"].values), 1, len(t) - 1)
        dist[g.index.values] = np.minimum(np.abs(g["mid"].values - t[k - 1]),
                                          np.abs(g["mid"].values - t[k]))
U["log_dist_tss"] = np.log10(dist + 1)
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
gc, cpg = np.empty(len(U)), np.empty(len(U))
for i, (ch, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = str(fa[ch][max(0, mid - 250):mid + 250])
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
U["gc"], U["cpg_oe"] = gc, cpg
print(f"union: {len(U):,} regions, {int((~U['prom']).sum()):,} non-promoter", flush=True)


# ---------------------------------------------------------------- mark counts
def count(job):
    """Sum both strands of one 5' track pair in each window around every region midpoint."""
    name, pattern = job
    out = {w: np.zeros(len(U)) for w in WIN}
    big = max(WIN.values())
    for strand in ("plus", "minus"):
        bw = pyBigWig.open(pattern.format(strand))
        L = bw.chroms()
        for i, (c, m) in enumerate(zip(U["chr"].values, U["mid"].values)):
            if c not in L:
                continue
            s0 = max(0, m - big)
            v = np.nan_to_num(bw.values(c, s0, min(m + big, L[c]), numpy=True))
            for w, h in WIN.items():
                out[w][i] += np.abs(v[max(0, m - h - s0):m + h - s0]).sum()
        bw.close()
    return name, out


jobs = []
for c in CELLS:
    safe = c.replace("-", "_").lower()
    jobs.append((f"k4_{c}", f"{K4DIR}/{safe}_h3k4me1_5p_{{}}.bw"))
    jobs += [(f"k4rep{r}_{c}", f"{K4DIR}/{safe}_h3k4me1_rep{r}_5p_{{}}.bw")
             for r in range(1, NREP[c] + 1)]
jobs += [(f"k27_{c}", p) for c, p in K27.items()]
for _, p in jobs:
    for s in ("plus", "minus"):
        assert os.path.getsize(p.format(s)) > 0, p.format(s)

if os.path.exists(K4CACHE):
    K = pd.read_csv(K4CACHE, sep="\t")
    assert len(K) == len(U)
else:
    K = U[KEY].copy()
    with Pool(a.procs) as pool:
        for name, out in pool.imap_unordered(count, jobs):
            for w, v in out.items():
                K[f"{name}_{w}"] = v
            print(f"  counted {name}", flush=True)
    K.to_csv(K4CACHE, sep="\t", index=False)
    print(f"wrote {K4CACHE}", flush=True)

# ---------------------------------------------------------------- measures
pct = lambda v: (rankdata(v) - 0.5) / len(v)
cor = lambda x, y: float(np.corrcoef(x, y)[0, 1])


def Zmat(sel, c):
    at = pct(S.loc[sel, f"atac_{c}"].values)
    g = U.loc[sel]
    return np.column_stack([np.ones(len(at)), at, at ** 2, g["gc"].values, g["cpg_oe"].values,
                            g["log_dist_tss"].values])


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


def r2(y, Z):
    r = resid(y, Z)
    return 1 - r.var() / y.var()


nonprom = ~U["prom"].values
rows = []
for c in CELLS:
    sel = U[f"in_{c}"].values & nonprom
    Z = Zmat(sel, c)
    at = pct(S.loc[sel, f"atac_{c}"].values)
    p3 = resid(pct(S.loc[sel, f"p300_{c}"].values), Z)
    marks = [("H3K4me1", "k4")] + ([("H3K27ac", "k27")] if c in K27 else [])
    for mark, pre in marks:
        for w in WIN:
            y = pct(K.loc[sel, f"{pre}_{c}_{w}"].values)
            row = dict(cell=c, mark=mark, window=w, n=int(sel.sum()), rho_atac=cor(y, at),
                       r2_Z=r2(y, Z), p300_partial=cor(resid(y, Z), p3))
            if pre == "k4":
                R = [resid(pct(K.loc[sel, f"k4rep{r}_{c}_{w}"].values), Z)
                     for r in range(1, NREP[c] + 1)]
                row["rep_resid"] = np.mean([cor(x, z) for x, z in itertools.combinations(R, 2)])
            rows.append(row)
    print(f"{c}: {int(sel.sum()):,} non-promoter regions", flush=True)
PC = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== per cell: the mark against ATAC, its reproducible residual, and p300 beyond ATAC ===")
print(PC.round(3).to_string(index=False))

# cross-cell agreement of residuals, regions in both cells
rows = []
for c, d in itertools.combinations(CELLS, 2):
    sel = U[f"in_{c}"].values & U[f"in_{d}"].values & nonprom
    Zc, Zd = Zmat(sel, c), Zmat(sel, d)
    row = dict(pair=f"{c}-{d}", n=int(sel.sum()),
               atac_raw=cor(pct(S.loc[sel, f"atac_{c}"].values), pct(S.loc[sel, f"atac_{d}"].values)),
               p300_resid=cor(resid(pct(S.loc[sel, f"p300_{c}"].values), Zc),
                              resid(pct(S.loc[sel, f"p300_{d}"].values), Zd)))
    for w in WIN:
        row[f"k4_resid_{w}"] = cor(resid(pct(K.loc[sel, f"k4_{c}_{w}"].values), Zc),
                                   resid(pct(K.loc[sel, f"k4_{d}_{w}"].values), Zd))
    rows.append(row)
XC = pd.DataFrame(rows)
print("\n=== cross-cell correlation of residuals (each on its own cell's ATAC + composition) ===")
print(XC.round(3).to_string(index=False))

# leave-one-cell-out: other cells' mean residual against this cell's p300 residual
rows = []
for c in CELLS:
    sel = U[f"in_{c}"].values & nonprom
    Zc = Zmat(sel, c)
    p3 = resid(pct(S.loc[sel, f"p300_{c}"].values), Zc)
    others = [d for d in CELLS if d != c]
    row = dict(cell=c, n=int(sel.sum()))
    proxy = np.mean([resid(pct(S.loc[sel, f"p300_{d}"].values), Zmat(sel, d)) for d in others], 0)
    row["loo_p300_from_p300"] = cor(resid(proxy, Zc), p3)
    for w in WIN:
        proxy = np.mean([resid(pct(K.loc[sel, f"k4_{d}_{w}"].values), Zmat(sel, d))
                         for d in others], 0)
        row[f"loo_p300_from_k4_{w}"] = cor(resid(proxy, Zc), p3)
        own = resid(pct(K.loc[sel, f"k4_{c}_{w}"].values), Zc)
        row[f"own_k4_vs_others_k4_{w}"] = cor(resid(proxy, Zc), own)
    rows.append(row)
LO = pd.DataFrame(rows)
print("\n=== leave-one-cell-out: do the other four cells' residuals predict this cell's? ===")
print(LO.round(3).to_string(index=False))

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
PC.to_csv(f"{a.out_prefix}.per_cell.tsv", sep="\t", index=False)
XC.to_csv(f"{a.out_prefix}.cross_cell.tsv", sep="\t", index=False)
LO.to_csv(f"{a.out_prefix}.loo.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.per_cell.tsv, .cross_cell.tsv, .loo.tsv")
