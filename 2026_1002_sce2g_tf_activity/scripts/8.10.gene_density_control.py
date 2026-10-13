#!/usr/bin/env python3
"""8.10 Is the CpG-island distance gain local gene density? A control on the 8.0 replica.

WHY. F-061: the K562 CRISPR gain from CGI distance (F-059) comes from the graded far range, not the
near-CGI depletion found in every cell (F-060). Most promoters sit on CpG islands, so being far from
any CGI largely means being in a gene-poor region. scE2G's n6 sees gene density only through
numTSSEnhGene / numCandidateEnhGene, and 8.3's distance control used the target gene's own TSS.
If the gain is gene density, the feature is a neighbourhood effect (fewer competing genes), not a
correction to activity, and F-060's activity result would not explain it.

GENE-DENSITY COVARIATES G, label-free, each 1 + percentile over all 155,134 Xu elements (as 8.2):
  tss_any_logdist  log distance from the element midpoint to the nearest TSS of any gene (7.1's
                   `log_dist_tss`, GENCODE v43 collapsed gene bounds)
  genes_1mb        number of gene TSSs within +/- 500 kb of the element midpoint
Models (8.0 replica, 2026_1007 discovery run tables, LOCO + paired bootstrap 2,000):
  n6, n6 + G, n6 + pct, n6 + G + pct, n6 + 3 (Phase 4), n6 + G + 3

READING, fixed before the run:
- pct "survives gene density" if its gain on top of n6 + G is >= 50% of its gain over n6 (+0.023,
  F-061) with the interval above 0. Then CGI distance carries something beyond gene density.
- Otherwise call the far-range gain largely gene density; the external test should then also carry
  G, and the F-060 activity link does not explain the CRISPR gain.
- n6 + G vs n6 is descriptive (whether gene density alone helps scE2G).

Usage: 8.10.gene_density_control.py [--out results/2026_1009_gene_density_control.tsv]
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("replica", os.path.join(HERE, "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
TSS = ("/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
KEY = ["chr", "start", "end", "TargetGene", "name"]
EK = ["chr", "start", "end"]
HALF = 500_000
ap = argparse.ArgumentParser()
ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1007_k562_xu_discovery/K562_Xu_pl")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out", default="results/2026_1009_gene_density_control.tsv")
a = ap.parse_args()

pct1 = lambda v: 1 + rankdata(np.nan_to_num(np.asarray(v, float), nan=-np.inf)) / len(v)  # noqa: E731

# ---------------------------------------------------------------- element-level covariates
E = pd.read_csv(f"{SC}/feature_screen/library.elements.tsv.gz", sep="\t",
                usecols=EK + ["mid", "log_dist_tss", "cgi_log_dist", "atac_xu"])
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "s", "e", "gene"])
tss["m"] = (tss["s"] + tss["e"]) // 2
tss = tss.drop_duplicates("gene")
E["genes_1mb"] = 0
for c, g in E.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == c, "m"].values)
    E.loc[g.index, "genes_1mb"] = (np.searchsorted(t, g["mid"].values + HALF, side="right")
                                   - np.searchsorted(t, g["mid"].values - HALF, side="left"))
E["tss_any_logdist"] = pct1(E["log_dist_tss"])
E["genes_1mb_p"] = pct1(E["genes_1mb"])
E["cgi_p"] = pct1(E["cgi_log_dist"])
print(f"{len(E):,} Xu elements; {len(tss):,} genes; genes within +/- 500 kb median "
      f"{E['genes_1mb'].median():.0f} (IQR {E['genes_1mb'].quantile(.25):.0f}-{E['genes_1mb'].quantile(.75):.0f})", flush=True)
print("Spearman over elements: CGI distance vs nearest-TSS distance "
      f"{spearmanr(E['cgi_log_dist'], E['log_dist_tss'], nan_policy='omit')[0]:.3f}, "
      f"vs genes within 500 kb {spearmanr(E['cgi_log_dist'], E['genes_1mb'], nan_policy='omit')[0]:.3f}", flush=True)

# ---------------------------------------------------------------- pairs
df, miss, params = rep.load_run(a.run, "multiome_powerlaw_n6")
cand = pd.read_csv(f"{SC}/crispr_discovery/candidates.tsv.gz", sep="\t").drop(columns=["Regulated"]).drop_duplicates(KEY)
MOT = ["motif:GATA6_MA1104.3+2", "motif:IKZF2_MA2326.1"]
df = df.merge(cand[KEY + ["cgi_log_dist"] + MOT], on=KEY, how="left", validate="1:1")
# map each CRISPR element to the overlapping Xu element with the highest Xu ATAC, no overlap -> 1.0
# (lowest), exactly as 8.2 built the candidates
el = df[EK].drop_duplicates().reset_index(drop=True)
best = np.full(len(el), -1)
for c, g in el.groupby("chr"):
    u_ = E.index[E["chr"] == c].values
    if len(u_) == 0:
        continue
    u_ = u_[np.argsort(E.loc[u_, "start"].values)]
    us, ue, at = E.loc[u_, "start"].values, E.loc[u_, "end"].values, E.loc[u_, "atac_xu"].values
    lo = np.searchsorted(ue, g["start"].values, side="right")
    hi = np.searchsorted(us, g["end"].values, side="left")
    for i, l, h in zip(g.index.values, lo, hi):
        if h > l:
            best[i] = u_[l + np.argmax(at[l:h])]
COV = ["tss_any_logdist", "genes_1mb_p", "cgi_p"]
el[COV] = np.where(best[:, None] >= 0, E[COV].values[np.maximum(best, 0)], 1.0)
print(f"{len(el):,} CRISPR elements, {(best >= 0).mean():.4f} overlap a Xu element", flush=True)
n0 = len(df)
df = df.merge(el, on=EK, how="left", validate="m:1")
assert len(df) == n0 and df[["cgi_log_dist"] + COV + MOT].notna().all().all()
bad = ~np.isclose(df["cgi_log_dist"], df["cgi_p"])
assert not bad.any(), f"{bad.sum()} pairs: mapped CGI percentile differs from 8.2's candidates (mapping differs)"
df = df.rename(columns={"cgi_log_dist": "pct", "genes_1mb_p": "genes_1mb"})
y = df["Regulated"].astype(bool).values
print(f"{len(df):,} pairs, {y.sum()} regulated (+ {len(miss):,} missing, {miss.sum()} regulated)", flush=True)

N6 = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene", "ARC.E2G.Score"]
G = ["tss_any_logdist", "genes_1mb"]
THREE = ["pct"] + MOT
MODELS = {"n6": N6, "n6 + G": N6 + G, "n6 + pct": N6 + ["pct"], "n6 + G + pct": N6 + G + ["pct"],
          "n6 + 3": N6 + THREE, "n6 + G + 3": N6 + G + THREE}
S = {m: rep.loco(df, f, params) for m, f in MODELS.items()}
ys = np.r_[y, miss]
Sx = {m: np.r_[s, np.zeros(len(miss))] for m, s in S.items()}
rng = np.random.default_rng(0)
nb = len(ys)
PAIRS = [("n6 + G", "n6"), ("n6 + pct", "n6"), ("n6 + G + pct", "n6 + G"), ("n6 + G + pct", "n6"),
         ("n6 + 3", "n6"), ("n6 + G + 3", "n6 + G"), ("n6 + G + 3", "n6 + 3")]
boot = {p: [] for p in PAIRS}
for _ in range(a.n_boot):
    ii = rng.integers(0, nb, nb)
    au = {m: rep.auprc(ys[ii], Sx[m][ii]) for m in MODELS}
    for m, r in PAIRS:
        boot[(m, r)].append(au[m] - au[r])
A = {m: rep.auprc(y, S[m], miss.sum(), len(miss)) for m in MODELS}
R = pd.DataFrame([dict(model=m, ref=r, auprc=A[m], diff=A[m] - A[r], lo=np.quantile(boot[(m, r)], 0.025),
                       hi=np.quantile(boot[(m, r)], 0.975)) for m, r in PAIRS])
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))

g0 = A["n6 + pct"] - A["n6"]
gG = A["n6 + G + pct"] - A["n6 + G"]
lo = R.set_index(["model", "ref"]).loc[("n6 + G + pct", "n6 + G"), "lo"]
share = gG / g0
print(f"\npct on top of G keeps {share:.0%} of its gain over n6 ({gG:+.4f} of {g0:+.4f}, lower bound {lo:+.4f}) -> "
      + ("survives gene density" if (share >= 0.5 and lo > 0) else "largely gene density"))
R.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
print(f"wrote {a.out}")
