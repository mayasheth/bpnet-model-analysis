#!/usr/bin/env python3
"""Feature A: per-element Kendall between peak accessibility and expression of the TFs whose
motif is in the peak (Xu K562 multiome), and its gate against observed K562 p300.

KENDALL. scE2G's Kendall peak x cell ATAC matrix (5.0 export; 7,821 cells), binarised
(accessible / not), against log1p(CP10k) RNA of each TF detected in >= --min-detect of the same
cells. Kendall tau-b for a binary x reduces to a Mann-Whitney U on y's midranks, so every
peak x TF tau comes from one sparse product (peaks x cells) @ (cells x TFs); checked against
scipy.stats.kendalltau on random pairs.

MOTIFS. 4.44 genome-wide FIMO hits (JASPAR 2024 CORE vertebrates, p < 1e-4); a hit belongs to
the peak containing its centre. A motif names one TF, or two for a dimer (GATA1::TAL1).

ELEMENT FEATURES (scE2G candidate regions, max over overlapping Kendall peaks):
  kA_max      max tau over expressed TFs with a motif in the peak
  kA_sumpos   sum of positive tau over those TFs
  n_tf        number of expressed TFs with a motif in the peak
  kA_ery_max  max tau over GATA1, TAL1, NFE2, KLF1 when their motif is present
  *_perm      same with the motif -> TF assignment permuted among expressed TFs: keeps each
              peak's number of TFs, breaks the motif-TF match. The null for "a TF's own motif".
  kA_any_max  max tau over all expressed TFs, motif or not: generic peak covariation.

GATE (as 2.4; partial Spearman with p300 over non-promoter elements, Z = Xu ATAC + square, GC,
CpG o/e, log distance to TSS; --z-encode-atac adds ENCODE K562 ATAC). TEST 1: each feature over
all elements with a Kendall peak. TEST 2: kA_max beyond Z + n_tf + kA_max_perm + kA_any_max,
elements with >= 1 expressed-TF motif. Bar from F-050/F-051: ~0.1 partial rho did not move
CRISPR, so only a much larger value justifies a scE2G run.

Run in the tobias pixi env. Usage: 5.1.peak_tf_kendall.py [--z-encode-atac] [--out-prefix ...]
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import kendalltau, rankdata

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--atac-prefix", default=f"{SC}/featA/atac_matrix")
ap.add_argument("--rna", default="/oak/stanford/groups/engreitz/Projects/scE2G/data/K562_260122/K562_Xu/rna_count_matrix.csv.gz")
ap.add_argument("--motif-dir", default=f"{SC}/motif_hits")
ap.add_argument("--elements", default=f"{A2}/results/2026_1003_footprint_p300_gate.element_features.tsv.gz")
ap.add_argument("--min-detect", type=float, default=0.05)
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--z-encode-atac", action="store_true")
ap.add_argument("--out-prefix", default=f"{A2}/results/2026_1004_featA_peak_tf_kendall")
a = ap.parse_args()
if a.z_encode_atac:
    a.out_prefix += ".z_encode_atac"
KEY = ["chr", "start", "end"]
ERY = {"GATA1", "TAL1", "NFE2", "KLF1"}
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- ATAC
pre = a.atac_prefix
nr, nc = (int(x) for x in open(f"{pre}.dim.txt").read().split())
A = sp.csc_matrix((np.fromfile(f"{pre}.x.float64", np.float64),
                   np.fromfile(f"{pre}.i.int32", np.int32),
                   np.fromfile(f"{pre}.p.int32", np.int32)), shape=(nr, nc))
peaks = pd.Series(open(f"{pre}.rows.txt").read().split("\n")[:nr])
cells = open(f"{pre}.cols.txt").read().split("\n")[:nc]
B = A.tocsr()
B.data[:] = 1.0
n1 = np.asarray(B.sum(1)).ravel()
P = peaks.str.extract(r"^(chr[^-]+)-(\d+)-(\d+)$")
assert P.notna().all().all(), "unparsed peak names"
P.columns = KEY
P[["start", "end"]] = P[["start", "end"]].astype(int)
print(f"ATAC {nr:,} peaks x {nc:,} cells; accessible in median {np.median(n1):.0f} cells", flush=True)

# ---------------------------------------------------------------- motifs -> TF symbols
names = list(np.load(f"{a.motif_dir}/chr1.npz", allow_pickle=True)["names"])
motif_tfs = []
for n in names:
    sym = re.sub(r"\(VAR\.\d+\)", "", n.split(" ", 1)[1].upper())
    motif_tfs.append(sym.split("::"))
tf_all = sorted({t for ts in motif_tfs for t in ts})

# ---------------------------------------------------------------- RNA (stream; keep motif TFs)
cix = None
lib = None
rows = {}
with open(a.rna) as fh:  # plain CSV despite the .gz suffix
    head = fh.readline().rstrip("\n").split(",")[1:]
    pos = {c: i for i, c in enumerate(head)}
    missing = [c for c in cells if c not in pos]
    assert not missing, f"{len(missing)} ATAC cells not in RNA"
    cix = np.array([pos[c] for c in cells])
    lib = np.zeros(len(cells))
    want = set(tf_all)
    for line in fh:
        g, rest = line.split(",", 1)
        v = np.fromstring(rest, sep=",")[cix]
        lib += v
        if g in want:
            rows[g] = v
print(f"RNA: {len(head):,} cells in matrix, {len(rows):,} of {len(tf_all):,} motif TFs present", flush=True)
tfs = sorted(g for g, v in rows.items() if (v > 0).mean() >= a.min_detect)
Y = np.column_stack([np.log1p(rows[g] / lib * 1e4) for g in tfs])
print(f"{len(tfs)} TFs detected in >= {a.min_detect:.0%} of cells", flush=True)

# ---------------------------------------------------------------- tau-b, peak x TF
N = nc
npairs = N * (N - 1) / 2
R = np.column_stack([rankdata(Y[:, j]) for j in range(Y.shape[1])])
Ty = np.array([((c := np.unique(Y[:, j], return_counts=True)[1]) * (c - 1) / 2).sum() for j in range(Y.shape[1])])
R1 = np.asarray(B @ R)  # peaks x TFs, sum of midranks over accessible cells
n0 = N - n1
U1 = R1 - (n1 * (n1 + 1) / 2)[:, None]
S = 2 * U1 - (n0 * n1)[:, None]
Tx = n0 * (n0 - 1) / 2 + n1 * (n1 - 1) / 2
den = np.sqrt((npairs - Tx)[:, None] * (npairs - Ty)[None, :])
with np.errstate(invalid="ignore", divide="ignore"):
    TAU = np.where(den > 0, S / den, 0.0).astype(np.float32)
del R1, U1, S, den
chk = []
for _ in range(30):
    i, j = rng.integers(nr), rng.integers(len(tfs))
    x = B[i].toarray().ravel()
    chk.append(abs(kendalltau(x, Y[:, j]).statistic - TAU[i, j]) if x.any() and not x.all() else 0.0)
print(f"tau check vs scipy, 30 random pairs: max |diff| {max(chk):.2e}", flush=True)
assert max(chk) < 1e-5, "rank-sum tau does not match scipy"
print("tau distribution:", pd.Series(TAU.ravel()[rng.integers(0, TAU.size, 1_000_000)]).describe().round(4).to_dict(), flush=True)

# ---------------------------------------------------------------- motif hits -> (peak, TF)
tfi = {t: k for k, t in enumerate(tfs)}
m2t = [[tfi[t] for t in ts if t in tfi] for ts in motif_tfs]
pairs = []
for f in sorted(glob.glob(f"{a.motif_dir}/chr*.npz")):
    c = os.path.basename(f)[:-4]
    pc = P.index[P["chr"] == c].values
    if len(pc) == 0:
        continue
    o = np.argsort(P.loc[pc, "start"].values)
    pc = pc[o]
    ps, pe = P.loc[pc, "start"].values, P.loc[pc, "end"].values
    assert (ps[1:] >= pe[:-1]).all(), f"{c}: overlapping Kendall peaks"
    z = np.load(f)
    mid = (z["start"].astype(np.int64) + z["end"]) // 2
    k = np.searchsorted(ps, mid, side="right") - 1
    ok = (k >= 0) & (mid < pe[np.clip(k, 0, None)])
    pm = np.unique(np.stack([pc[k[ok]], z["motif"][ok].astype(np.int64)], 1), axis=0)
    pairs.append(pm)
    print(f"{c}: {ok.sum():,} hits in peaks, {len(pm):,} peak-motif pairs", flush=True)
PM = np.concatenate(pairs)
MT = pd.DataFrame([(m, t) for m, ts in enumerate(m2t) for t in ts], columns=["m", "t"])
PT = pd.DataFrame(PM, columns=["p", "m"]).merge(MT, on="m")[["p", "t"]].drop_duplicates().to_numpy(np.int64)
print(f"{len(PT):,} peak-TF pairs (expressed TFs with a motif in the peak)", flush=True)


def aggregate(pt, tag):
    v = TAU[pt[:, 0], pt[:, 1]]
    mx = np.full(nr, np.nan)
    np.fmax.at(mx, pt[:, 0], v)
    sp_ = np.zeros(nr)
    np.add.at(sp_, pt[:, 0], np.clip(v, 0, None))
    n = np.bincount(pt[:, 0], minlength=nr)
    return {f"kA_max{tag}": mx, f"kA_sumpos{tag}": sp_, f"n_tf{tag}": n}


F = pd.DataFrame(aggregate(PT, ""))
perm = rng.permutation(len(tfs))
PTp = np.unique(np.stack([PT[:, 0], perm[PT[:, 1]]], 1), axis=0)
F = F.join(pd.DataFrame(aggregate(PTp, "_perm")).drop(columns="n_tf_perm"))
ery = np.array([tfs[t] in ERY for t in PT[:, 1]])
F["kA_ery_max"] = aggregate(PT[ery], "")["kA_max"] if ery.any() else np.nan
F["kA_any_max"] = TAU.max(1)
F["n_access"] = n1
F = pd.concat([P, F], axis=1)

# ---------------------------------------------------------------- peaks -> scE2G elements
E = pd.read_csv(a.elements, sep="\t", usecols=KEY + ["prom", "log_dist_tss", "gc", "cpg_oe", "atac_xu", "atac_encode", "p300"])
FEAT = ["kA_max", "kA_sumpos", "n_tf", "kA_max_perm", "kA_sumpos_perm", "kA_ery_max", "kA_any_max"]
for col in FEAT + ["n_peaks"]:
    E[col] = np.nan
for c, g in E.groupby("chr"):
    fp = F[F["chr"] == c].sort_values("start")
    if fp.empty:
        continue
    lo = np.searchsorted(fp["end"].values, g["start"].values, side="right")
    hi = np.searchsorted(fp["start"].values, g["end"].values, side="left")
    E.loc[g.index, "n_peaks"] = hi - lo
    has = hi > lo
    for col in FEAT:
        v = fp[col].values
        E.loc[g.index[has], col] = [np.nanmax(v[l:h]) if np.isfinite(v[l:h]).any() else np.nan
                                    for l, h in zip(lo[has], hi[has])]
print(f"{len(E):,} elements; {(E['n_peaks'] > 0).mean():.3f} overlap a Kendall peak; "
      f"{(E['n_tf'] > 0).mean():.3f} have an expressed-TF motif", flush=True)

# ---------------------------------------------------------------- gate
pct = lambda v: (rankdata(v) - 0.5) / len(v)


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


def partial(x, y, Z):
    rx, ry = resid(x, Z), resid(y, Z)
    est = float(np.corrcoef(rx, ry)[0, 1])
    bs = []
    for _ in range(a.n_boot):
        i = rng.integers(0, len(rx), len(rx))
        bs.append(np.corrcoef(rx[i], ry[i])[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return est, lo, hi


def Zmat(g, extra=()):
    at = pct(g["atac_xu"].values)
    cols = [np.ones(len(g)), at, at ** 2, g["gc"].values, g["cpg_oe"].values, g["log_dist_tss"].values]
    if a.z_encode_atac:
        ae = pct(g["atac_encode"].values)
        cols += [ae, ae ** 2]
    return np.column_stack(cols + [pct(g[c].values) for c in extra])


G = E[(~E["prom"]) & (E["n_peaks"] > 0)].reset_index(drop=True)
G[FEAT] = G[FEAT].fillna(G[FEAT].min() - 1)  # no expressed-TF motif -> lowest rank
Z = Zmat(G)
p3 = pct(G["p300"].values)
rows = []
for f in FEAT + (["atac_encode"] if not a.z_encode_atac else []):
    est, lo, hi = partial(pct(G[f].values), p3, Z)
    rows.append(dict(feature=f, n=len(G), partial_rho=est, ci_lo=lo, ci_hi=hi))
T1 = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== test 1: partial Spearman with p300 beyond Z, non-promoter elements with a Kendall peak ===")
print(T1.round(3).to_string(index=False))

g = G[G["n_tf"] > 0].reset_index(drop=True)
rows = []
for label, extra in [("Z + n_tf", ("n_tf",)), ("Z + n_tf + perm", ("n_tf", "kA_max_perm")),
                     ("Z + n_tf + perm + any", ("n_tf", "kA_max_perm", "kA_any_max"))]:
    est, lo, hi = partial(pct(g["kA_max"].values), pct(g["p300"].values), Zmat(g, extra))
    rows.append(dict(controls=label, n_elements=len(g), partial_rho_kA_max=est, ci_lo=lo, ci_hi=hi))
T2 = pd.DataFrame(rows)
print("\n=== test 2: kA_max beyond motif presence and the permuted / any-TF nulls ===")
print(T2.round(3).to_string(index=False))
strict = T2.iloc[-1]
print(f"\nGATE (kA_max beyond Z + n_tf + perm + any, interval excludes 0): "
      f"{'PASS' if strict['ci_lo'] > 0 else 'FAIL'}  {strict['partial_rho_kA_max']:.3f} "
      f"[{strict['ci_lo']:.3f}, {strict['ci_hi']:.3f}]")

T1.to_csv(f"{a.out_prefix}.all_elements.tsv", sep="\t", index=False)
T2.to_csv(f"{a.out_prefix}.beyond_motif.tsv", sep="\t", index=False)
E.to_csv(f"{a.out_prefix}.element_features.tsv.gz", sep="\t", index=False)
F.to_csv(f"{SC}/featA/peak_features.tsv.gz", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.{{all_elements,beyond_motif}}.tsv, .element_features.tsv.gz")
