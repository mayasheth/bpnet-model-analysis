#!/usr/bin/env python3
"""8.7 Promoter-side CpG features for scE2G (Maya's suggestion, 2026-10-07), on the 8.0 replica.

HYPOTHESIS, stated before the run: CpG-island (housekeeping-like) promoters respond less to distal
enhancers, so promoter CpG status should lower a pair's score beyond ubiqExpressed. A pre-specified test of
three gene-level features, not a search, so LOCO + paired bootstrap on the replica (no nested selection).

PROMOTER FEATURES (target gene TSS +/- 500 bp from the training table's chrTSS / startTSS):
  prom_gc        GC fraction
  prom_cpg_oe    CpG observed / expected
  prom_cgi_frac  share of the window inside a UCSC CpG island
each encoded 1 + percentile over the benchmark's unique TSSs (label-free), like 8.2.
ELEMENT FEATURES: the three Phase 4 features from 8.2's candidates (cgi_log_dist, GATA cluster, IKZF2).

Usage: 8.7.promoter_cpg.py [--out results/2026_1007_promoter_cpg.tsv]
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd
import pysam
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("replica", os.path.join(HERE, "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
GENOME = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
CGI = f"{SC}/ucsc_hg38/cpgIslandExt.txt.gz"
KEY = ["chr", "start", "end", "TargetGene", "name"]
HALF = 500
ap = argparse.ArgumentParser()
ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1005_k562_xu_screen/K562_Xu_pl")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out", default="results/2026_1007_promoter_cpg.tsv")
a = ap.parse_args()

df, miss, params = rep.load_run(a.run, "multiome_n6_gc")
cand = pd.read_csv(f"{SC}/crispr_discovery/candidates.tsv.gz", sep="\t").drop(columns=["Regulated"]).drop_duplicates(KEY)
ELEM = ["cgi_log_dist", "motif:GATA6_MA1104.3+2", "motif:IKZF2_MA2326.1"]
df = df.merge(cand[KEY + ELEM], on=KEY, how="left", validate="1:1")
assert df[ELEM].notna().all().all()

T = df[["chrTSS", "startTSS"]].drop_duplicates().dropna().astype({"startTSS": int}).reset_index(drop=True)
fa = pysam.FastaFile(GENOME)
G = pd.read_csv(CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])
vals = np.zeros((len(T), 3))
for i, (c, t) in enumerate(zip(T["chrTSS"], T["startTSS"])):
    s = fa.fetch(c, max(0, t - HALF), t + HALF).upper()
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    g = G[(G["chr"] == c) & (G["e"] > t - HALF) & (G["s"] < t + HALF)]
    cov = sum(min(e, t + HALF) - max(b, t - HALF) for b, e in zip(g["s"], g["e"]))
    vals[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1), min(cov / (2 * HALF), 1.0)
PROM = ["prom_gc", "prom_cpg_oe", "prom_cgi_frac"]
for k, f in enumerate(PROM):
    T[f] = 1 + rankdata(vals[:, k]) / len(T)
print(f"{len(T):,} unique TSSs; CpG-island promoter (>= half the window in a CGI) at {(vals[:, 2] >= 0.5).mean():.3f}", flush=True)
df = df.merge(T, on=["chrTSS", "startTSS"], how="left")
n_na = df[PROM].isna().any(axis=1).sum()
df[PROM] = df[PROM].fillna(1.5)
print(f"{n_na} pairs without a TSS -> median percentile", flush=True)
cgi_prom = df["prom_cgi_frac"] > df["prom_cgi_frac"].quantile(0.5)
print("regulated rate, non-CGI vs CGI-overlapping promoter half:",
      round(df.loc[~cgi_prom, "Regulated"].mean(), 4), round(df.loc[cgi_prom, "Regulated"].mean(), 4), flush=True)

N6 = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene", "ARC.E2G.Score"]
MODELS = {
    "n6": N6,
    "n6 + prom_cgi_frac": N6 + ["prom_cgi_frac"],
    "n6 + prom_cpg_oe": N6 + ["prom_cpg_oe"],
    "n6 + prom_gc": N6 + ["prom_gc"],
    "n6 + 3 prom": N6 + PROM,
    "n6 + 3 element (Phase 4)": N6 + ELEM,
    "n6 + 3 element + 3 prom": N6 + ELEM + PROM,
}
y = df["Regulated"].astype(bool).values
S = {m: rep.loco(df, f, params) for m, f in MODELS.items()}
ys = np.r_[y, miss]
Sx = {m: np.r_[s, np.zeros(len(miss))] for m, s in S.items()}
rng = np.random.default_rng(0)
nb = len(ys)
REFS = {"n6": "n6", "n6 + 3 element + 3 prom": "n6 + 3 element (Phase 4)"}
boot = {(m, r): [] for m in MODELS for r in REFS}
for _ in range(a.n_boot):
    ii = rng.integers(0, nb, nb)
    au = {m: rep.auprc(ys[ii], Sx[m][ii]) for m in MODELS}
    for m in MODELS:
        for r, ref in REFS.items():
            boot[(m, r)].append(au[m] - au[ref])
rows = []
for m in MODELS:
    au = rep.auprc(y, S[m], miss.sum(), len(miss))
    row = dict(model=m, auprc=au, vs_n6=au - rep.auprc(y, S["n6"], miss.sum(), len(miss)),
               lo=np.quantile(boot[(m, "n6")], 0.025), hi=np.quantile(boot[(m, "n6")], 0.975))
    rows.append(row)
R = pd.DataFrame(rows)
e3, ep = "n6 + 3 element (Phase 4)", "n6 + 3 element + 3 prom"
d = np.array(boot[(ep, "n6 + 3 element + 3 prom")])
print(R.round(4).to_string(index=False))
print(f"\n3 prom on top of 3 element: {R.set_index('model').loc[ep, 'auprc'] - R.set_index('model').loc[e3, 'auprc']:+.4f} "
      f"[{np.quantile(d, 0.025):+.4f}, {np.quantile(d, 0.975):+.4f}]")
X = np.log(np.abs(df[N6 + PROM].values.astype(float)) + rep.EPS)
c = LogisticRegression(**params).fit(X, y.astype(int)).coef_[0]
print("coefficients, n6 + 3 prom:", {f: round(v, 3) for f, v in zip(N6 + PROM, c)})
R.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
print(f"wrote {a.out}")
