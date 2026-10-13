#!/usr/bin/env python3
"""8.3 Phase 2 of plans/8_crispr_feature_discovery.md: feature selection inside nested chromosome CV.

BASE      the five non-ARC n6 features + ABC.Score + exp(Kendall) (8.1: within 0.004 of n6 with ARC), so
          candidates can be read against ABC and Kendall separately.
OUTER     5 chromosome groups balanced for positives; each group is predicted by a model whose features
          were chosen on the other 4 only. Pooled held-out AUPRC (missing CRISPR pairs at score 0, as 4.2)
          vs the base under the same folds, paired bootstrap over pairs.
INNER     4-fold by chromosome group within the outer training set; criterion = pooled inner AUPRC.
FEATURES  8.2 candidates (1 + percentile), log(|x| + 0.01) like every scE2G feature.

PROCEDURES (fixed before any run)
  forward   add the candidate with the best inner AUPRC, up to 5, stop when the gain < 0.002
  l1        L1 logistic over base + candidates + candidate x log ABC + candidate x tau (all
            standardised), C from an inner grid
  family    each family added as a block with L2 penalty (C from an inner grid); every family's nested
            gain is reported, and the procedure keeps the family with the best inner AUPRC if > +0.002
  trees     gradient-boosted trees on base + candidates (depth 3), iterations by inner early stopping

Usage: 8.3.nested_selection.py [--procedures forward,l1,family,trees] [--feature-perm-seed S] [--out-prefix P]
  --feature-perm-seed (Phase 3 null): labels and base untouched; whole candidate vectors are shuffled among
  pairs within 10 TSS-distance bins, so selection searches features with the library's structure but no
  link to the labels. Its gain is the optimism of searching the library.
  --base-logdist adds log TSS distance to the base (control: does a candidate only proxy distance?). The
  73 pairs without a TSS position are left unshuffled by the feature null (none is a positive).
  --permute-seed (within-gene label shuffle) is kept only for the record: INVALID as a null (2026-10-06):
  218 of 429 positives sit in single-pair or all-positive genes the shuffle cannot change, and it moves
  positives distal (median 28 -> 45 kb).
"""
import argparse
import importlib.util
import os
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, precision_recall_curve

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("replica", os.path.join(HERE, "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1005_k562_xu_screen/K562_Xu_pl")
ap.add_argument("--cand-dir", default=f"{SC}/crispr_discovery")
ap.add_argument("--procedures", default="forward,l1,family,trees")
ap.add_argument("--permute-seed", type=int, default=None, help="INVALID null, kept for the record")
ap.add_argument("--feature-perm-seed", type=int, default=None)
ap.add_argument("--base-logdist", action="store_true",
                help="control: add log10 TSS distance to the base, so a candidate cannot gain by proxying distance")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--jobs", type=int, default=5)
ap.add_argument("--out-prefix", default=None)
a = ap.parse_args()
assert a.permute_seed is None or a.feature_perm_seed is None
tag = ("real" if a.permute_seed is None and a.feature_perm_seed is None
       else f"perm{a.permute_seed}" if a.permute_seed is not None else f"fperm{a.feature_perm_seed}")
tag += "_basedist" if a.base_logdist else ""
a.out_prefix = a.out_prefix or f"{a.cand_dir}/nested_{tag}"
PROCS = a.procedures.split(",")
KEY = ["chr", "start", "end", "TargetGene", "name"]
BASE5 = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene"]
MIN_GAIN = 0.002
C_GRID = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0]

df, miss, params = rep.load_run(a.run, "multiome_n6_gc")
cand = pd.read_csv(f"{a.cand_dir}/candidates.tsv.gz", sep="\t")
fam = pd.read_csv(f"{a.cand_dir}/families.tsv", sep="\t")
cand = cand.drop(columns=["Regulated"]).drop_duplicates(KEY)
df = df.merge(cand, on=KEY, how="left", validate="1:1")
CANDS = fam["feature"].tolist()
assert df[CANDS].notna().all().all(), "candidates missing for some pairs"
df["Kendall_exp"] = np.exp(df["Kendall"])
BASE = BASE5 + ["ABC.Score", "Kendall_exp"]
if a.base_logdist:
    d_ = (((df["start"] + df["end"]) // 2) - df["startTSS"]).abs()
    # 10 ** log10(d + 1): the log transform below returns ~log10-scaled distance; 73 pairs lack a TSS -> median
    df["dist_tss"] = (d_.fillna(d_.median()) + 1).values
    BASE.append("dist_tss")
y = df["Regulated"].astype(bool).values
if a.permute_seed is not None:
    r = np.random.default_rng(a.permute_seed)
    y = y.copy()
    for _, ix in df.groupby("TargetGene").indices.items():
        y[ix] = y[ix][r.permutation(len(ix))]
n_pos_miss, n_miss = int(miss.sum()), len(miss)

L = lambda cols: np.log(np.abs(df[cols].values.astype(float)) + rep.EPS)  # noqa: E731
XB, XC, XA = L(BASE), L(CANDS), L(["ARC.E2G.Score"])
if a.feature_perm_seed is not None:
    dist = (((df["start"] + df["end"]) // 2) - df["startTSS"]).abs().values
    dbin = pd.qcut(np.log10(dist + 1), 10, labels=False, duplicates="drop")
    r = np.random.default_rng(a.feature_perm_seed)
    XC0 = XC.copy()
    for b in np.unique(dbin):
        ix = np.flatnonzero(dbin == b)
        XC[ix] = XC0[ix[r.permutation(len(ix))]]
    print(f"feature null: candidate vectors shuffled within {len(np.unique(dbin))} distance bins; "
          f"{(XC != XC0).any(1).mean():.3f} of pairs changed", flush=True)
lab, ltau = XB[:, BASE.index("ABC.Score")], XB[:, BASE.index("Kendall_exp")]
XI = np.column_stack([(XC - XC.mean(0)) * (lab - lab.mean())[:, None], (XC - XC.mean(0)) * (ltau - ltau.mean())[:, None]])
INAMES = CANDS + [f"{c} x logABC" for c in CANDS] + [f"{c} x tau" for c in CANDS]

# 5 chromosome groups balanced for positives (greedy, deterministic)
pos = pd.Series(y, index=df["chr"]).groupby(level=0).sum().sort_values(ascending=False)
grp, load = {}, np.zeros(5)
for c, n in pos.items():
    k = int(np.argmin(load))
    grp[c], load[k] = k, load[k] + n
G = df["chr"].map(grp).values
print(f"{len(df):,} pairs, {y.sum()} positives ({tag}); groups: positives {load.astype(int).tolist()}; "
      f"{len(CANDS)} candidates", flush=True)


def lr(X, yy, **kw):
    return LogisticRegression(**{**params, **kw}).fit(X, yy) if kw else LogisticRegression(**params).fit(X, yy)


def pr_auc(yy, s):
    p, r_, _ = precision_recall_curve(yy, s)
    return auc(r_, p)


def inner_cv(fit_predict, tr):
    """Pooled AUPRC over the 4 inner chromosome groups of an outer training set."""
    s = np.full(len(df), np.nan)
    for g in np.unique(G[tr]):
        itr, ite = tr & (G != g), tr & (G == g)
        s[ite] = fit_predict(itr, ite)
    return pr_auc(y[tr], s[tr])


def std_fit(X, tr):
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1
    return (X - mu) / sd


def run_fold(k, proc):
    tr, te = G != k, G == k
    chosen = []
    if proc == "base":
        return k, lr(XB[tr], y[tr]).predict_proba(XB[te])[:, 1], chosen
    if proc == "arc":
        X = np.column_stack([XB[:, :5], XA])
        return k, lr(X[tr], y[tr]).predict_proba(X[te])[:, 1], chosen
    if proc == "forward":
        sel, best = [], inner_cv(lambda i, j: lr(XB[i], y[i]).predict_proba(XB[j])[:, 1], tr)
        for _ in range(5):
            scores = {}
            for c in range(XC.shape[1]):
                if c in sel:
                    continue
                X = np.column_stack([XB, XC[:, sel + [c]]])
                scores[c] = inner_cv(lambda i, j: lr(X[i], y[i]).predict_proba(X[j])[:, 1], tr)
            c, s = max(scores.items(), key=lambda kv: kv[1])
            if s - best < MIN_GAIN:
                break
            sel.append(c)
            best = s
        X = np.column_stack([XB, XC[:, sel]])
        return k, lr(X[tr], y[tr]).predict_proba(X[te])[:, 1], [CANDS[c] for c in sel]
    if proc == "l1":
        X = std_fit(np.column_stack([XB, XC, XI]), tr)
        fit = lambda C, i: lr(X[i], y[i], penalty="l1", solver="liblinear", C=C, max_iter=5000)  # noqa: E731
        Cb = max(C_GRID, key=lambda C: inner_cv(lambda i, j: fit(C, i).predict_proba(X[j])[:, 1], tr))
        m = fit(Cb, tr)
        nz = np.flatnonzero(m.coef_[0][len(BASE):])
        return k, m.predict_proba(X[te])[:, 1], [INAMES[i] for i in nz] + [f"C={Cb}"]
    if proc == "family":
        res = {}
        for f, cs in fam.groupby("family")["feature"]:
            idx = [CANDS.index(c) for c in cs]
            X = std_fit(np.column_stack([XB, XC[:, idx]]), tr)
            fit = lambda C, i: lr(X[i], y[i], penalty="l2", C=C, solver="lbfgs", max_iter=5000)  # noqa: E731
            sc = {C: inner_cv(lambda i, j: fit(C, i).predict_proba(X[j])[:, 1], tr) for C in C_GRID}
            Cb = max(sc, key=sc.get)
            res[f] = (sc[Cb], fit(Cb, tr).predict_proba(X[te])[:, 1], Cb)
        base_in = inner_cv(lambda i, j: lr(XB[i], y[i]).predict_proba(XB[j])[:, 1], tr)
        bf = max(res, key=lambda f: res[f][0])
        keep = res[bf][0] - base_in >= MIN_GAIN
        s = res[bf][1] if keep else lr(XB[tr], y[tr]).predict_proba(XB[te])[:, 1]
        return k, s, ([bf] if keep else []) + [f"per_family:{f}:{res[f][0] - base_in:+.4f}" for f in res], \
            {f: res[f][1] for f in res}
    if proc == "trees":
        X = np.column_stack([XB, XC])
        m = HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=500, early_stopping=True,
                                           validation_fraction=0.2, random_state=0).fit(X[tr], y[tr])
        return k, m.predict_proba(X[te])[:, 1], [f"iters={m.n_iter_}"]
    raise ValueError(proc)


def total_auprc(s):
    return rep.auprc(y, s, n_pos_miss, n_miss)


S, CH, FAMS = {}, {}, {}
for proc in ["base", "arc"] + PROCS:
    out = Parallel(n_jobs=a.jobs)(delayed(run_fold)(k, proc) for k in range(5))
    s = np.full(len(df), np.nan)
    for res in out:
        k, sk, ch = res[:3]
        s[G == k] = sk
        CH.setdefault(proc, {})[k] = ch
        if len(res) > 3:
            for f, v in res[3].items():
                FAMS.setdefault(f, np.full(len(df), np.nan))[G == k] = v
    S[proc] = s
    print(f"{proc}: nested AUPRC {total_auprc(s):.4f}", flush=True)
for f, v in FAMS.items():
    S[f"family:{f}"] = v

ys = np.r_[y, np.ones(n_pos_miss, bool), np.zeros(n_miss - n_pos_miss, bool)]
Sx = {m: np.r_[s, np.zeros(n_miss)] for m, s in S.items()}
rng = np.random.default_rng(0)
nb = len(ys)
boot = {m: np.empty(a.n_boot) for m in S}
for b in range(a.n_boot):
    ii = rng.integers(0, nb, nb)
    a0 = pr_auc(ys[ii], Sx["base"][ii])
    for m in S:
        boot[m][b] = pr_auc(ys[ii], Sx[m][ii]) - a0
q = lambda v, x: np.quantile(v, x) if len(v) else np.nan  # noqa: E731  (--n-boot 0 for permutations)
R = pd.DataFrame([dict(model=m, auprc=total_auprc(S[m]), vs_base=total_auprc(S[m]) - total_auprc(S["base"]),
                       lo=q(boot[m], 0.025), hi=q(boot[m], 0.975)) for m in S])
R["permute_seed"] = a.permute_seed
R.to_csv(f"{a.out_prefix}.auprc.tsv", sep="\t", index=False, float_format="%.5g")
rows = [dict(procedure=p, outer_group=k, chosen=c) for p, d in CH.items() for k, cs in d.items() for c in cs]
pd.DataFrame(rows).to_csv(f"{a.out_prefix}.chosen.tsv", sep="\t", index=False)
print(R.round(4).to_string(index=False))
for p in PROCS:
    if p == "family":
        continue
    cnt = pd.Series([c for cs in CH[p].values() for c in cs if not c.startswith(("C=", "iters="))]).value_counts()
    print(f"\n{p}: chosen in >= 2 of 5 outer folds\n{cnt[cnt >= 2].to_string() if len(cnt) else 'none'}")
print(f"\nwrote {a.out_prefix}.auprc.tsv, .chosen.tsv")
