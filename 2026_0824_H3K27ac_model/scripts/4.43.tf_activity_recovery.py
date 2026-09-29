#!/usr/bin/env python3
"""Test 2: can the motifs that matter for p300 in a cell be recovered from its ATAC alone?

WHY. 4.41 found the transferred p300 model mis-ranks elements by which TFs are active in K562
(GATA under-predicted, CEBP/ESR over-predicted). A model can only be told that at deployment if
the important motifs can be estimated from the target cell's ATAC and sequence. This measures
how well they can, in all five p300-panel cell types, before any model is built.

REGIONS. One union set: every ABC candidate region called from each cell's own ATAC (K562's July
run, the four others from 4.39/4.40), merged. Signals are counted from the bigwigs the models
train on (config/p300_panel.json), so every cell is counted the same way. Promoter = overlaps a
GENCODE TSS +/- 250 bp. A region counts as "in" cell c if it overlaps one of c's own candidates.

TRUTH, per cell, from that cell's p300 (not deployable):
    T[c, m] = coefficient of p300 rank on motif-m presence, over c's own non-promoter regions,
              partialling out ATAC rank (+ square), GC, CpG o/e and log distance to TSS.
    i.e. "does motif m predict p300 BEYOND accessibility in c". Also T_plain without ATAC.

ESTIMATES, per cell, from ATAC + sequence only (deployable):
    E_within[c, m] = coefficient of ATAC rank on motif presence over c's regions, partialling out
                     GC, CpG o/e and distance: does the motif mark MORE accessible sites in c.
    E_chromvar[c, m] = chromVAR deviation z-score over the union set: reads in motif-m regions in
                     cell c relative to the expectation from all five cells, against background
                     sets matched on GC and mean accessibility (50 permutations within bins).
    E_spec[c, m]   = E_within[c, m] minus its mean over the other four cells: the motif's
                     accessibility association that is SPECIFIC to c. Ubiquitous motifs (CTCF)
                     cancel without being removed by hand.
    E_discrim[c, m] = coefficient, over all union non-promoter regions, of c's accessibility
                     SPECIFICITY (c's ATAC rank minus the mean rank in the other cells) on motif
                     presence, partialling out mean ATAC rank, GC, CpG o/e and distance: which
                     motifs mark sites that are more open in c than elsewhere.

SCORES. Per cell, Spearman across motifs between each estimate and the truth, and top-20
overlap. Then the contrast that matters for transfer to K562: (K562 - panel mean) for the
estimates against the same for the truth, and against 4.41's transfer-specific residual
coefficients.

Usage: 4.43.tf_activity_recovery.py [--out-prefix results/2026_0926_tf_activity_recovery]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import pyBigWig
from scipy import sparse
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CANDS = {"K562": f"{ABC}/2026_0905_p300_activity/p300obs_k562",
         **{c: f"{ABC}/2026_0924_p300_panel_obs/p300obs_{c.replace('-', '_')}"
            for c in ("GM12878", "HepG2", "A549", "MCF-7")}}
PANEL = json.load(open(f"{A}/config/p300_panel.json"))
GENOME = f"{D}/hg38_resources/hg38.fa"
JASPAR = f"{P}/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
ANAT = f"{A}/results/2026_0925_p300_error_anatomy.features.tsv"
KEY = ["chr", "start", "end"]
HALF = 250
CELLS = list(CANDS)

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_0926_tf_activity_recovery")
ap.add_argument("--cache-dir", default=f"{A}/data/tf_activity_recovery")
ap.add_argument("--n-bg", type=int, default=50)
a = ap.parse_args()
os.makedirs(a.cache_dir, exist_ok=True)


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


# ---------------------------------------------------------------- union region set
cand = {}
for c, d in CANDS.items():
    cand[c] = pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t", usecols=KEY)
    print(f"{c}: {len(cand[c]):,} candidate regions", flush=True)
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
print(f"union: {len(U):,} regions", flush=True)
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


# ---------------------------------------------------------------- signals from bigwigs
def bw_sum(path, reg, absval=False):
    bw = pyBigWig.open(path)
    ch = bw.chroms()
    out = np.zeros(len(reg))
    for i, (c, s, e) in enumerate(zip(reg["chr"].values, reg["mid"].values - HALF,
                                      reg["mid"].values + HALF)):
        if c not in ch:
            continue
        v = bw.values(c, max(0, int(s)), min(int(e), ch[c]), numpy=True)
        v = np.nan_to_num(v)
        out[i] = np.abs(v).sum() if absval else v.sum()
    bw.close()
    return out


sig_cache = f"{a.cache_dir}/union_signals.tsv.gz"
if os.path.exists(sig_cache) and len(pd.read_csv(sig_cache, sep="\t", usecols=["chr"])) == len(U):
    S = pd.read_csv(sig_cache, sep="\t")
else:
    S = U[KEY].copy()
    for c in CELLS:
        sp = PANEL[c]
        S[f"atac_{c}"] = bw_sum(sp["accessibility_bw"], U)
        S[f"p300_{c}"] = bw_sum(sp["signal_plus_bw"], U, True) + \
            bw_sum(sp["signal_minus_bw"], U, True)
        print(f"signals: {c}", flush=True)
    S.to_csv(sig_cache, sep="\t", index=False)

# ---------------------------------------------------------------- sequence, GC, motifs
import pyfaidx  # noqa: E402

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
LUT = np.full(256, -1, np.int8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
idx = np.empty((len(U), 2 * HALF), np.int8)
for i, (c, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = str(fa[c][max(0, mid - HALF):mid + HALF]).encode().ljust(2 * HALF, b"N")
    idx[i] = LUT[np.frombuffer(s, np.uint8)]
valid = idx >= 0
U["gc"] = ((idx == 1) | (idx == 2)).sum(1) / np.maximum(valid.sum(1), 1)
cg = ((idx[:, :-1] == 1) & (idx[:, 1:] == 2)).sum(1)
U["cpg_oe"] = cg * valid.sum(1) / np.maximum((idx == 1).sum(1) * (idx == 2).sum(1), 1)

mcache = f"{a.cache_dir}/union_fimo_presence.npz"
if os.path.exists(mcache):
    z = np.load(mcache, allow_pickle=True)
    M, mnames = z["presence"], list(z["names"])
    assert M.shape[0] == len(U)
else:
    from memelite import fimo
    from memelite.io import read_meme
    motifs = read_meme(JASPAR)
    mnames = list(motifs)
    M = np.zeros((len(U), len(mnames)), bool)
    B = 20000
    for s in range(0, len(U), B):
        n = min(B, len(U) - s)
        X = np.zeros((n, 4, 2 * HALF), np.float32)
        ii, jj = np.nonzero(valid[s:s + n])
        X[ii, idx[s:s + n][ii, jj], jj] = 1.0
        for j, h in enumerate(fimo(motifs, X, threshold=1e-4)):
            if len(h):
                M[s + h["sequence_name"].astype(int).values, j] = True
        print(f"  fimo {s + n:,}/{len(U):,}", flush=True)
    np.savez_compressed(mcache, presence=M, names=np.array(mnames, dtype=object))
keep = (M.sum(0) >= 200) & (M.sum(0) <= len(U) - 200)
M, mnames = M[:, keep], [n for n, k in zip(mnames, keep) if k]
print(f"{len(mnames)} motifs with >= 200 hits", flush=True)


# ---------------------------------------------------------------- truth and within-cell estimate
def fwl_coefs(y, covs, F):
    """OLS coefficient of y on each column of F, partialling out covs (FWL)."""
    Q, _ = np.linalg.qr(np.column_stack([np.ones(len(y))] + covs))
    yr = y - Q @ (Q.T @ y)
    Fr = F - Q @ (Q.T @ F)
    return (Fr * yr[:, None]).sum(0) / (Fr * Fr).sum(0)


pct = lambda v: 100 * (rankdata(v) - 0.5) / len(v)
T, T_plain, E_within = {}, {}, {}
for c in CELLS:
    sel = U[f"in_{c}"].values & ~U["prom"].values
    g = U[sel]
    atac, p3 = pct(S.loc[sel, f"atac_{c}"].values), pct(S.loc[sel, f"p300_{c}"].values)
    comp = [g["gc"].values, g["cpg_oe"].values, g["log_dist_tss"].values]
    F = M[sel].astype(float)
    T[c] = fwl_coefs(p3, [atac, atac ** 2] + comp, F)
    T_plain[c] = fwl_coefs(p3, comp, F)
    E_within[c] = fwl_coefs(atac, comp, F)
    print(f"truth/within: {c} ({int(sel.sum()):,} non-promoter regions)", flush=True)

E_spec = {c: E_within[c] - np.mean([E_within[o] for o in CELLS if o != c], 0) for c in CELLS}
selu = ~U["prom"].values
PR = {c: pct(S.loc[selu, f"atac_{c}"].values) for c in CELLS}
meanpr = np.mean([PR[c] for c in CELLS], 0)
gu = U[selu]
compu = [gu["gc"].values, gu["cpg_oe"].values, gu["log_dist_tss"].values]
Fu = M[selu].astype(float)
E_discrim = {}
for c in CELLS:
    spec = PR[c] - np.mean([PR[o] for o in CELLS if o != c], 0)
    E_discrim[c] = fwl_coefs(spec, [meanpr] + compu, Fu)
print("specificity and discriminative estimates done", flush=True)

# ---------------------------------------------------------------- chromVAR-style deviations
X = np.column_stack([S[f"atac_{c}"].values for c in CELLS])
X = X / X.sum(0, keepdims=True) * 1e6                     # CPM per cell
e = X.sum(1) / X.sum()                                      # expected fraction per region
Ms = sparse.csr_matrix(M.astype(np.float32))


def deviations(Xm):
    obs = Ms.T @ Xm                                         # motifs x cells
    exp = (Ms.T @ e)[:, None] * Xm.sum(0, keepdims=True)
    return (obs - exp) / exp


dev = deviations(X)
mean_acc = np.log1p(X.mean(1))
gbin = np.clip((rankdata(U["gc"]) / len(U) * 10).astype(int), 0, 9)
abin = np.clip((rankdata(mean_acc) / len(U) * 10).astype(int), 0, 9)
bins = gbin * 10 + abin
rng = np.random.default_rng(0)
bg = np.empty((a.n_bg,) + dev.shape)
order = [np.nonzero(bins == b)[0] for b in range(100)]
for it in range(a.n_bg):
    perm = np.arange(len(U))
    for ix in order:                                        # permute regions within bins
        perm[ix] = rng.permutation(ix)
    bg[it] = deviations(X[perm])
E_cv = (dev - bg.mean(0)) / bg.std(0)
E_chromvar = {c: E_cv[:, i] for i, c in enumerate(CELLS)}
print("chromVAR deviations done", flush=True)

# ---------------------------------------------------------------- scores
rows = []
top = lambda v, k=20: set(np.argsort(-v)[:k])
for c in CELLS:
    for est_name, E in (("within-cell ATAC association", E_within), ("cell-specific ATAC association", E_spec),
                        ("discriminative (open here, not elsewhere)", E_discrim),
                        ("chromVAR deviation", E_chromvar)):
        for tname, TT in (("beyond ATAC", T), ("plain", T_plain)):
            rows.append(dict(cell=c, estimate=est_name, truth=tname,
                             spearman=spearmanr(E[c], TT[c])[0],
                             top20_overlap=len(top(E[c]) & top(TT[c]))))
R = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== per cell: estimate vs truth across motifs ===")
print(R.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

panel = [c for c in CELLS if c != "K562"]
dT = T["K562"] - np.mean([T[c] for c in panel], 0)
contrast = []
EST = dict((("within-cell ATAC association", E_within), ("cell-specific ATAC association", E_spec),
                        ("discriminative (open here, not elsewhere)", E_discrim),
                        ("chromVAR deviation", E_chromvar)))
for est_name, E in EST.items():
    dE = E["K562"] - np.mean([E[c] for c in panel], 0)
    contrast.append(dict(estimate=est_name, spearman_vs_truth_contrast=spearmanr(dE, dT)[0],
                         top20_overlap=len(top(dE) & top(dT))))
an = pd.read_csv(ANAT, sep="\t")
an = an[(an.cls == "nonpromoter") & an.feature.str.startswith("motif:")].copy()
an["motif"] = an.feature.str[len("motif:"):]
rt = an.set_index("motif")["coef_r_transfer"].reindex(mnames)
ok = rt.notna().values
print("\n=== K562 minus panel: does the estimated contrast match the true one? ===")
for r in contrast:
    E = EST[r["estimate"]]
    dE = E["K562"] - np.mean([E[c] for c in panel], 0)
    r["spearman_vs_4.41_transfer_residual"] = spearmanr(dE[ok], rt.values[ok])[0]
    print(r)
print(f"truth contrast vs 4.41 transfer residual: "
      f"{spearmanr(dT[ok], rt.values[ok])[0]:.3f}")

out = pd.DataFrame({"motif": mnames,
                    **{f"truth_{c}": T[c] for c in CELLS},
                    **{f"truth_plain_{c}": T_plain[c] for c in CELLS},
                    **{f"within_{c}": E_within[c] for c in CELLS},
                    **{f"chromvar_{c}": E_chromvar[c] for c in CELLS},
                    **{f"spec_{c}": E_spec[c] for c in CELLS},
                    **{f"discrim_{c}": E_discrim[c] for c in CELLS}})
out["truth_contrast_K562"] = dT
os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
out.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
R.to_csv(f"{a.out_prefix}.scores.tsv", sep="\t", index=False)
gi = [i for i, n in enumerate(mnames) if "GATA" in n.upper() or "TRPS1" in n.upper()]
ci = [i for i, n in enumerate(mnames) if n.split(" ", 1)[1].upper().startswith(("CEBP", "ESR"))]
for lab, ii in (("GATA family", gi), ("CEBP/ESR", ci)):
    print(f"\n{lab}: truth contrast (K562 - panel) {dT[ii].mean():+.2f}")
    for nm, E in EST.items():
        print(f"    {nm:<44} " + ", ".join(f"{c} {E[c][ii].mean():+.2f}" for c in CELLS))
ct = [i for i, n in enumerate(mnames) if n.split(" ", 1)[1].upper().startswith("CTCF")]
print("\nwhere CTCF motifs rank (of %d; 1 = top) in K562:" % len(mnames))
for nm, E in list(EST.items()) + [("TRUTH beyond ATAC", T)]:
    rk = rankdata(-E["K562"])
    print(f"    {nm:<44} " + ", ".join(str(int(rk[i])) for i in ct))
print(f"\nwrote {a.out_prefix}.motifs.tsv and .scores.tsv")
