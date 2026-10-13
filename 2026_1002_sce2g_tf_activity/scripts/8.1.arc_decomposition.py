#!/usr/bin/env python3
"""8.1 ARC.E2G.Score vs its parts (ABC.Score, Kendall) as separate scE2G features, on the 8.0 replica.

WHY. On a deep sample like Xu, ARC-E2G performs like ABC or Kendall alone; it combines them to stay
calibrated across depth. Splitting it lets the logistic regression weight each part, take products
(interaction terms), and lets later phases test features against each part separately.

KENDALL SIGN. Kendall tau can be negative and scE2G's log(|x| + 0.01) folds -t onto +t, so four
encodings: raw (as the pipeline would treat it), clipped at 0, 1 + tau, and exp(tau), which matches
how ARC itself uses tau (log ARC = log ABC + w * tau, compute_arc_e2g.R).
Interaction terms are products of the log-transformed features, added after the transform.

Paired bootstrap over pairs (missing CRISPR pairs at score 0, as 4.2), 2,000 resamples, vs n6.

Usage: 8.1.arc_decomposition.py [--run DIR] [--out results/2026_1006_arc_decomposition.tsv]
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

# the replica's file name starts with a digit, so load it by path
spec = importlib.util.spec_from_file_location("replica", os.path.join(os.path.dirname(os.path.abspath(__file__)), "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1005_k562_xu_screen/K562_Xu_pl")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out", default="results/2026_1006_arc_decomposition.tsv")
a = ap.parse_args()
rng = np.random.default_rng(0)

df, miss, params = rep.load_run(a.run, "multiome_n6_gc")
print(f"{len(df):,} pairs, {int(df['Regulated'].sum())} positives; + {len(miss):,} missing ({int(miss.sum())} positives)", flush=True)
print(f"Kendall < 0 at {(df['Kendall'] < 0).mean():.3f} of pairs; spearman ARC vs ABC "
      f"{df['ARC.E2G.Score'].rank().corr(df['ABC.Score'].rank()):.3f}, ARC vs Kendall "
      f"{df['ARC.E2G.Score'].rank().corr(df['Kendall'].rank()):.3f}, ABC vs Kendall {df['ABC.Score'].rank().corr(df['Kendall'].rank()):.3f}", flush=True)
df["Kendall_clip0"] = df["Kendall"].clip(lower=0)
df["Kendall_plus1"] = 1 + df["Kendall"]
# ARC (compute_arc_e2g.R) is log ABC + w * tau with tau signed and linear; exp(tau) makes the
# log(|x| + 0.01) transform return ~tau, so the regression can fit w freely
df["Kendall_exp"] = np.exp(df["Kendall"])
BASE = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene"]


def fit(feats, inter=()):
    """LOCO as 8.0, with optional products of log-transformed features."""
    X = np.log(np.abs(df[feats].values.astype(float)) + rep.EPS)
    col = {f: i for i, f in enumerate(feats)}
    if inter:
        X = np.column_stack([X] + [X[:, col[u]] * X[:, col[v]] for u, v in inter])
    y = df["Regulated"].values.astype(np.int64)
    out = np.full(len(df), np.nan)
    for c in np.unique(df["chr"]):
        te = (df["chr"] == c).values
        out[te] = LogisticRegression(**params).fit(X[~te], y[~te]).predict_proba(X[te])[:, 1]
    return out


K = "Kendall_clip0"
MODELS = {
    "n6 (ARC)": (BASE + ["ARC.E2G.Score"], ()),
    "ABC only": (BASE + ["ABC.Score"], ()),
    "Kendall raw only": (BASE + ["Kendall"], ()),
    "Kendall clip0 only": (BASE + [K], ()),
    "ABC + Kendall raw": (BASE + ["ABC.Score", "Kendall"], ()),
    "ABC + Kendall clip0": (BASE + ["ABC.Score", K], ()),
    "ABC + Kendall plus1": (BASE + ["ABC.Score", "Kendall_plus1"], ()),
    "ABC + Kendall exp": (BASE + ["ABC.Score", "Kendall_exp"], ()),
    "ABC + Kendall exp + GC": (BASE + ["ABC.Score", "Kendall_exp", "elementGC"], ()),
    "ABC + Kendall clip0 + ABCxKendall": (BASE + ["ABC.Score", K], (("ABC.Score", K),)),
    "ARC + ABC + Kendall clip0": (BASE + ["ARC.E2G.Score", "ABC.Score", K], ()),
    "ABC split: ATAC_enh + ABC + Kendall clip0": (BASE + ["normalizedATAC_enh", "ABC.Score", K], ()),
    "n6 + GC": (BASE + ["ARC.E2G.Score", "elementGC"], ()),
    "ABC + Kendall clip0 + GC": (BASE + ["ABC.Score", K, "elementGC"], ()),
}
y = df["Regulated"].astype(bool).values
S = {m: fit(f, i) for m, (f, i) in MODELS.items()}
ys = np.r_[y, miss]
Sx = {m: np.r_[s, np.zeros(len(miss))] for m, s in S.items()}
n = len(ys)
rows = []
for m in MODELS:
    rows.append(dict(model=m, auprc=rep.auprc(y, S[m], miss.sum(), len(miss))))
R = pd.DataFrame(rows).set_index("model")
ref = "n6 (ARC)"
boot = {m: [] for m in MODELS}
for _ in range(a.n_boot):
    ii = rng.integers(0, n, n)
    a0 = rep.auprc(ys[ii], Sx[ref][ii])
    for m in MODELS:
        boot[m].append(rep.auprc(ys[ii], Sx[m][ii]) - a0)
R["vs_n6"] = R["auprc"] - R.loc[ref, "auprc"]
R["lo"] = [np.quantile(boot[m], 0.025) for m in MODELS]
R["hi"] = [np.quantile(boot[m], 0.975) for m in MODELS]
coef_rows = []
for m in ("ABC + Kendall clip0", "ABC + Kendall exp", "ABC + Kendall clip0 + ABCxKendall", "ARC + ABC + Kendall clip0"):
    f, inter = MODELS[m]
    X = np.log(np.abs(df[f].values.astype(float)) + rep.EPS)
    col = {g: i for i, g in enumerate(f)}
    names = list(f) + [f"{u} x {v}" for u, v in inter]
    if inter:
        X = np.column_stack([X] + [X[:, col[u]] * X[:, col[v]] for u, v in inter])
    c = LogisticRegression(**params).fit(X, df["Regulated"].values.astype(np.int64)).coef_[0]
    coef_rows += [dict(model=m, feature=g, coefficient=v) for g, v in zip(names, c)]
print(R.round(4).to_string(), flush=True)
print(pd.DataFrame(coef_rows).round(3).to_string(index=False))
R.reset_index().to_csv(a.out, sep="\t", index=False, float_format="%.5g")
print(f"wrote {a.out}")
