#!/usr/bin/env python3
"""7.2 Screen the 7.1 feature library against observed K562 p300 beyond Xu ATAC (plan in
plans/7_feature_screen.md). No CRISPR label is read.

TARGET  rank(observed p300) residualised on Z = Xu ATAC percentile and its square, log distance
        to TSS, promoter flag, mappability (technical: short-read ChIP), and ranked GC and CpG o/e; every feature is ranked and residualised
        on Z the same way, within each chromosome half, so the univariate statistic is a partial
        Spearman. GC and CpG o/e are in Z because they dominate otherwise (smoke test chr21-22: GC
        -0.55 beyond ATAC alone) and because the 2.4 gate that set the 0.15 bar controlled for them.
        rho_*_base columns repeat the univariate step without them (--no-z-composition: drop them).
HALVES  discovery = odd autosomes; replication = even autosomes + chrX.

1. univariate  partial Spearman per feature on both halves; BH q on discovery; replication interval
               from a 1 Mb block bootstrap (features with discovery q < 0.05 only).
2. models      gradient-boosted trees and an elastic net on the residual, fitted on one half and
               scored on the other, both directions: held-out R2 = share of the p300-beyond-ATAC
               rank variance recovered. Out-of-fold tree predictions = one composite per element.
3. stability   lasso (alpha from LassoCV on discovery) on --n-stab half-samples of discovery 1 Mb
               blocks; selection frequency.
SHORTLIST      the rule fixed in the plan: same sign, replication interval excluding 0,
               |rho_rep| >= 0.15, selection >= 0.8, |rho| < 0.7 with any feature already listed; at most 5.

Run in the tobias pixi env. Usage: 7.2.screen_features.py [--lib-dir DIR] [--out-prefix P] [--n-stab 50]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as tdist
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import ElasticNetCV, Lasso, LassoCV
from sklearn.model_selection import GroupKFold

SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--lib-dir", default=f"{SC}/feature_screen")
ap.add_argument("--out-prefix", default="results/2026_1005_feature_screen")
ap.add_argument("--n-stab", type=int, default=50)
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--min-present", type=float, default=0.005, help="drop motif features in fewer elements")
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--no-z-composition", action="store_true", help="leave GC and CpG o/e out of Z")
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
KEY = ["chr", "start", "end"]
COMP = [] if a.no_z_composition else ["gc", "cpg_oe"]
# mappability is technical: p300 ChIP reads are short and drop in repeats (smoke test: top hit, 0.32)
TECH = ["mappability"]
CONTROL = ["atac_xu", "log_dist_tss", "prom"] + TECH + COMP
NOT_FEATURES = set(KEY) | {"mid", "p300"} | set(CONTROL)
CALIB = ["motif_n:CTCF_MA0139.2", "ref_ery_fp_print", "ref_ery_fp_tobias", "ref_ctcf_fp_print",
         "ref_fp_any_tobias"]

E = pd.read_csv(f"{a.lib_dir}/library.elements.tsv.gz", sep="\t")
M = np.load(f"{a.lib_dir}/library.motifs.npz", allow_pickle=True)
assert (M["chr"] == E["chr"].values).all() and (M["start"] == E["start"].values).all(), "library row order differs"
names = [str(x) for x in M["names"]]
present = (M["n"] > 0).mean(0)
keepm = present >= a.min_present
X = pd.concat([E[[c for c in E.columns if c not in NOT_FEATURES]].reset_index(drop=True),
               pd.DataFrame(M["n"][:, keepm].astype(np.float32), columns=[f"motif_n:{n}" for n, k in zip(names, keepm) if k]),
               pd.DataFrame(M["score_max"][:, keepm], columns=[f"motif_max:{n}" for n, k in zip(names, keepm) if k])], axis=1)
eryth = [f"motif_n:{n}" for n in names if n.split("_")[0] in {"GATA1", "GATA1::TAL1", "NFE2", "MAF::NFE2", "TAL1::TCF3"}]
X["motif_n:erythroid_sum"] = X[[c for c in eryth if c in X]].sum(1)
CALIB.append("motif_n:erythroid_sum")
fam = np.array(["motif" if c.startswith("motif") else "reference" if c.startswith("ref_") else "library" for c in X.columns])
for c in ("width", "gc", "atac_xu"):
    r = pd.Series(X[[k for k in X if k.startswith("motif_n:")]].sum(1)).rank().corr(E[c].rank())
    print(f"  total motif sites vs {c}: spearman {r:.3f}", flush=True)
print(f"{len(E):,} elements; {X.shape[1]:,} features ({keepm.sum()} of {len(names)} motifs in >= "
      f"{a.min_present:.1%} of elements)", flush=True)

chrom = E["chr"].str.replace("chr", "")
odd = chrom.map(lambda c: c.isdigit() and int(c) % 2 == 1).values
halves = {"discovery": odd, "replication": ~odd}
assert all(v.sum() > 1000 for v in halves.values()), {h: int(v.sum()) for h, v in halves.items()}
block = (E["chr"] + ":" + (E["mid"] // 1_000_000).astype(str)).values


def zmat(idx, comp):
    at = rankdata(E.loc[idx, "atac_xu"]) / idx.sum()
    Z = [np.ones(idx.sum()), at, at ** 2, E.loc[idx, "log_dist_tss"], E.loc[idx, "prom"].astype(float)]
    Z += [rankdata(E.loc[idx, c]) / idx.sum() for c in TECH + comp]
    return np.linalg.qr(np.column_stack(Z))[0]


def resid_ranks(V, Q):
    R = rankdata(np.asarray(V, dtype=float), axis=0)
    return R - Q @ (Q.T @ R)


Y, XR = {}, {}
for h, idx in halves.items():
    Q = zmat(idx, COMP)
    Y[h] = resid_ranks(E.loc[idx, "p300"].values, Q)
    XR[h] = resid_ranks(np.nan_to_num(X.loc[idx].values.astype(float), nan=-1.0), Q)


def corr_cols(A, y):
    A = A - A.mean(0)
    y = y - y.mean()
    den = np.sqrt((A ** 2).sum(0) * (y ** 2).sum())
    return np.where(den > 0, (A * y[:, None]).sum(0) / np.where(den > 0, den, 1), 0.0)


# ---------------------------------------------------------------- 1. univariate
U = pd.DataFrame({"feature": X.columns, "family": fam, "present_frac": (np.nan_to_num(X.values) != 0).mean(0)})
for h, idx in halves.items():
    U[f"rho_{h[:4]}"] = corr_cols(XR[h], Y[h])
base = {}
if COMP:
    for h, idx in halves.items():
        Q0 = zmat(idx, [])
        V = np.column_stack([np.nan_to_num(X.loc[idx].values.astype(float), nan=-1.0), E.loc[idx, COMP].values])
        base[h] = corr_cols(resid_ranks(V, Q0), resid_ranks(E.loc[idx, "p300"].values, Q0))
        U[f"rho_{h[:4]}_base"] = base[h][:X.shape[1]]
    print("composition in Z; beyond ATAC alone: " + ", ".join(
        f"{c} {base['discovery'][X.shape[1] + k]:.3f} / {base['replication'][X.shape[1] + k]:.3f}"
        for k, c in enumerate(COMP)) + " (discovery / replication)", flush=True)
n_d = halves["discovery"].sum()
df = n_d - (5 + len(TECH) + len(COMP)) - 1
tt = U["rho_disc"] * np.sqrt(df / np.maximum(1 - U["rho_disc"] ** 2, 1e-12))
U["p_disc"] = 2 * tdist.sf(np.abs(tt), df)
o = np.argsort(U["p_disc"].values)
q = U["p_disc"].values[o] * len(U) / np.arange(1, len(U) + 1)
U.loc[o, "q_disc"] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
sel = np.where(U["q_disc"] < 0.05)[0]
rep = halves["replication"]
bl = block[rep]
ub = np.unique(bl)
bidx = {b: np.where(bl == b)[0] for b in ub}
boots = np.empty((a.n_boot, len(sel)))
for k in range(a.n_boot):
    ii = np.concatenate([bidx[b] for b in rng.choice(ub, len(ub))])
    boots[k] = corr_cols(XR["replication"][ii][:, sel], Y["replication"][ii])
U["rep_lo"], U["rep_hi"] = np.nan, np.nan
U.loc[sel, "rep_lo"], U.loc[sel, "rep_hi"] = np.quantile(boots, 0.025, 0), np.quantile(boots, 0.975, 0)
print(f"univariate: {len(sel):,} features at discovery q < 0.05; {len(ub)} replication blocks", flush=True)

# ---------------------------------------------------------------- 2. models, both directions
Xraw = X.values.astype(np.float32)
oof = np.full(len(E), np.nan)
rows = []
for tr, te in (("discovery", "replication"), ("replication", "discovery")):
    itr, ite = halves[tr], halves[te]
    gb = HistGradientBoostingRegressor(max_iter=600, learning_rate=0.05, max_leaf_nodes=31, early_stopping=True,
                                       validation_fraction=0.1, random_state=a.seed).fit(Xraw[itr], Y[tr])
    pred = gb.predict(Xraw[ite])
    oof[ite] = pred
    r2 = 1 - ((Y[te] - pred) ** 2).sum() / ((Y[te] - Y[te].mean()) ** 2).sum()
    sd = XR[tr].std(0)
    sd[sd == 0] = 1
    folds = list(GroupKFold(3).split(XR[tr], Y[tr], groups=block[itr]))
    en = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=folds, max_iter=2000, n_jobs=-1).fit(XR[tr] / sd, Y[tr])
    sd_te = XR[te].std(0)
    sd_te[sd_te == 0] = 1
    pe = en.predict(XR[te] / sd_te)
    r2e = 1 - ((Y[te] - pe) ** 2).sum() / ((Y[te] - Y[te].mean()) ** 2).sum()
    rows.append(dict(train=tr, test=te, r2_trees=r2, r2_elasticnet=r2e, enet_nonzero=int((en.coef_ != 0).sum())))
    print(f"  train {tr} -> test {te}: held-out R2 trees {r2:.4f}, elastic net {r2e:.4f} "
          f"({rows[-1]['enet_nonzero']} nonzero)", flush=True)
pd.DataFrame(rows).to_csv(f"{a.out_prefix}.models.tsv", sep="\t", index=False)
comp = E[KEY].copy()
comp["screen_composite"] = oof
comp.to_csv(f"{a.lib_dir}/screen_composite.tsv.gz", sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- 3. stability selection (discovery)
Xd = XR["discovery"]
sd = Xd.std(0)
sd[sd == 0] = 1
Xd = Xd / sd
folds = list(GroupKFold(3).split(Xd, Y["discovery"], groups=block[halves["discovery"]]))
alpha = LassoCV(alphas=20, cv=folds, max_iter=2000, n_jobs=-1).fit(Xd, Y["discovery"]).alpha_
bd = block[halves["discovery"]]
ubd = np.unique(bd)
hits = np.zeros(X.shape[1])
for k in range(a.n_stab):
    take = np.isin(bd, rng.choice(ubd, len(ubd) // 2, replace=False))
    hits += Lasso(alpha=alpha, max_iter=2000).fit(Xd[take], Y["discovery"][take]).coef_ != 0
U["stability"] = hits / a.n_stab
print(f"stability: lasso alpha {alpha:.3g}, {(U['stability'] >= 0.8).sum()} features selected in >= 80% of "
      f"{a.n_stab} half-samples", flush=True)

# ---------------------------------------------------------------- shortlist (rule fixed in the plan)
cand = U[(np.sign(U["rho_disc"]) == np.sign(U["rho_repl"])) & ((U["rep_lo"] > 0) | (U["rep_hi"] < 0))
         & (U["rho_repl"].abs() >= 0.15) & (U["stability"] >= 0.8) & (U["family"] != "reference")]
cand = cand.reindex(cand["rho_repl"].abs().sort_values(ascending=False).index)
short = []
for i in cand.index:
    if all(abs(np.corrcoef(XR["replication"][:, i], XR["replication"][:, j])[0, 1]) < 0.7 for j in short):
        short.append(i)
    if len(short) == 5:
        break
U["shortlist"] = U.index.isin(short)
U.sort_values("p_disc").to_csv(f"{a.out_prefix}.univariate.tsv", sep="\t", index=False, float_format="%.4g")
cols = ["feature", "rho_disc", "rho_repl", "rep_lo", "rep_hi", "q_disc", "stability"] + (["rho_repl_base"] if COMP else [])
print("\n=== calibration (known CRISPR outcomes) ===")
print(U[U["feature"].isin(CALIB)][cols].round(4).to_string(index=False))
print("\n=== top 25 by discovery p ===")
print(U.sort_values("p_disc").head(25)[cols].round(4).to_string(index=False))
print("\n=== shortlist ===")
print(U[U["shortlist"]][cols].round(4).to_string(index=False) if short else "none passes the rule")
print(f"\nwrote {a.out_prefix}.univariate.tsv, .models.tsv, {a.lib_dir}/screen_composite.tsv.gz")
