#!/usr/bin/env python3
"""Error analysis of the scE2G TF-evidence models (motif count, TOBIAS / PRINT footprints as features
and inside the activity term, and observed p300 when its run exists) on the K562 CRISPR pairs.

Why do element-level TF features that track p300 (F-049, F-051) not move CRISPR AUPRC (F-050, F-051,
F-053)? Four checks, over the 9,267 CRISPR pairs that have a candidate element (held-out-chromosome
E2G.Score from each model; pairs without a candidate are identical across models and left out):

  1 coverage   how many CRISPR elements and regulated pairs carry an erythroid site at all
  2 gate       footprint vs observed p300 beyond ATAC, recomputed on the CRISPR-tested elements only
               (the screen ran over all 135,965 candidates; the benchmark sees ~3,400 elements)
  3 residual   does a feature carry CRISPR information the baseline lacks? Logistic regression of
               Regulated on logit(baseline E2G.Score) plus the feature's percentile; coefficient and
               95% bootstrap interval (2,000 resamples of pairs). Run for footprints, motif count,
               observed p300, ATAC, GC, and p300 beyond ATAC (residual of p300 rank on ATAC rank).
  4 shifts     how each model moves pairs relative to the baseline: mean change in within-model
               percentile rank of regulated vs not-regulated pairs, split by erythroid site presence
               and by distance; a model only gains AUPRC if it moves regulated pairs up more.
Plus AUPRC by distance bin for every model.

Element features come from the 2.4 element tables (max over candidates overlapping the CRISPR
element). Run in the scprinter pixi env (sklearn). Usage: 6.1.sce2g_tf_error_analysis.py [--out-prefix ...]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, precision_recall_curve

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
RUNS = "/scratch/users/shethm/sce2g_runs"
ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default=f"{A2}/results/2026_1004_sce2g_tf_error_analysis")
ap.add_argument("--n-boot", type=int, default=2000)
a = ap.parse_args()
rng = np.random.default_rng(0)
KEY = ["chr", "start", "end", "TargetGene", "name"]

MODELS = {
    "baseline": ("2026_1003_k562_xu_featB_tobias", "multiome_powerlaw_n6"),
    "+ motif count": ("2026_1003_k562_xu_featB_tobias", "multiome_n6_ery_motif"),
    "+ TOBIAS footprint": ("2026_1003_k562_xu_featB_tobias", "multiome_n6_ery_fp"),
    "+ PRINT footprint": ("2026_1004_k562_xu_featB_print", "multiome_n6_ery_fp"),
    "activity x TOBIAS, w1": ("2026_1004b_k562_xu_abcfp", "multiome_n6_arcfp_tobias_l1"),
    "activity x PRINT, w1": ("2026_1004b_k562_xu_abcfp", "multiome_n6_arcfp_print_l1"),
    "activity x TOBIAS, w4": ("2026_1004b_k562_xu_abcfp", "multiome_n6_arcfp_tobias_l4"),
    "activity x PRINT, w4": ("2026_1004b_k562_xu_abcfp", "multiome_n6_arcfp_print_l4"),
    "activity x TOBIAS, w1 (exact)": ("2026_1004_k562_xu_activity", "multiome_n6_arcfp_tobias_l1"),
    "activity x PRINT, w1 (exact)": ("2026_1004_k562_xu_activity", "multiome_n6_arcfp_print_l1"),
    "activity x TOBIAS, w4 (exact)": ("2026_1004_k562_xu_activity", "multiome_n6_arcfp_tobias_l4"),
    "activity x PRINT, w4 (exact)": ("2026_1004_k562_xu_activity", "multiome_n6_arcfp_print_l4"),
    "+ observed p300": ("2026_1004_k562_xu_activity", "multiome_n6_plus_p300"),
    "activity sqrt(ATAC x p300)": ("2026_1004_k562_xu_activity", "multiome_n6_arcp300geo"),
    "activity p300": ("2026_1004_k562_xu_activity", "multiome_n6_arcp300only"),
    "+ GC": ("2026_1004_k562_xu_corr", "multiome_n6_plus_gc"),
    "+ CTCF footprint": ("2026_1004_k562_xu_corr", "multiome_n6_plus_ctcf_fp"),
    "+ CTCF motif count": ("2026_1004_k562_xu_corr", "multiome_n6_plus_ctcf_n"),
    "+ footprint pct": ("2026_1004_k562_xu_corr", "multiome_n6_plus_fp_pct"),
    "+ GC + CTCF + footprint pct": ("2026_1004_k562_xu_corr", "multiome_n6_corr_all"),
}
D = None
for lab, (run, m) in MODELS.items():
    f = f"{RUNS}/{run}/K562_Xu_pl/{m}/model/training_predictions.tsv"
    if not os.path.exists(f):
        print(f"skip {lab}: {f} missing", flush=True)
        continue
    t = pd.read_csv(f, sep="\t", usecols=KEY + ["Regulated", "startTSS", "E2G.Score"]).rename(columns={"E2G.Score": lab})
    D = t if D is None else D.merge(t.drop(columns="startTSS"), on=KEY + ["Regulated"], how="inner", validate="1:1")
labs = [c for c in MODELS if c in D.columns]
assert len(D) > 9000, len(D)
D["y"] = D["Regulated"].astype(bool)
D["dist"] = (np.abs((D["start"] + D["end"]) / 2 - D["startTSS"])).astype(float)
print(f"{len(D):,} pairs, {D['y'].sum()} regulated, {D.loc[D['y'], 'TargetGene'].nunique()} genes with a regulated pair; "
      f"models: {labs}", flush=True)

# ---------------------------------------------------------------- element features
FEAT = ["p300", "atac_xu", "atac_encode", "gc", "prom", "erythroid_n", "erythroid_fp_sum", "CTCF_n"]
E = {}
for src, f in (("tobias", "2026_1003_footprint_p300_gate.element_features.tsv.gz"),
               ("print", "2026_1003_footprint_p300_gate.print.element_features.tsv.gz")):
    E[src] = pd.read_csv(f"{A2}/results/{f}", sep="\t", usecols=["chr", "start", "end"] + FEAT)
U = E["tobias"].copy()
U["fp_print"] = E["print"]["erythroid_fp_sum"].values
assert (E["print"][["chr", "start", "end"]].values == U[["chr", "start", "end"]].values).all()
U = U.rename(columns={"erythroid_fp_sum": "fp_tobias"})
cols = FEAT[:-3] + ["erythroid_n", "fp_tobias", "fp_print", "CTCF_n"]
U["prom"] = U["prom"].astype(str).str.lower().eq("true").astype(float)
el = D[["chr", "start", "end"]].drop_duplicates().reset_index(drop=True)
vals = np.full((len(el), len(cols)), np.nan)
for c, g in el.groupby("chr"):
    u = U[U["chr"] == c].sort_values("start")
    if u.empty:
        continue
    us, ue = u["start"].values, u["end"].values
    assert (us[1:] >= ue[:-1]).all(), f"{c}: overlapping candidates"
    lo = np.searchsorted(ue, g["start"].values, side="right")
    hi = np.searchsorted(us, g["end"].values, side="left")
    X = u[cols].values
    for i, l, h in zip(g.index.values, lo, hi):
        if h > l:
            vals[i] = X[l:h].max(0)
el[cols] = vals
print(f"{len(el):,} CRISPR elements; {np.isfinite(vals[:, 0]).mean():.3f} overlap a candidate", flush=True)
D = D.merge(el, on=["chr", "start", "end"], how="left", validate="m:1")
pct = lambda v: (rankdata(v) - 0.5) / len(v)

# ---------------------------------------------------------------- 1 coverage
el_reg = D.groupby(["chr", "start", "end"])["y"].any()
cov = pd.DataFrame({
    "set": ["CRISPR elements", "elements with a regulated pair", "regulated pairs", "not-regulated pairs"],
    "n": [len(el), int(el_reg.sum()), int(D["y"].sum()), int((~D["y"]).sum())],
    "frac_with_erythroid_site": [
        float((el["erythroid_n"] > 0).mean()),
        float((el.set_index(["chr", "start", "end"]).loc[el_reg[el_reg].index, "erythroid_n"] > 0).mean()),
        float((D.loc[D["y"], "erythroid_n"] > 0).mean()), float((D.loc[~D["y"], "erythroid_n"] > 0).mean())],
})
print("\n=== 1 coverage ===\n" + cov.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 2 gate on CRISPR-tested elements
def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


def partial(x, y, Z, n=200):
    rx, ry = resid(x, Z), resid(y, Z)
    est = float(np.corrcoef(rx, ry)[0, 1])
    bs = [np.corrcoef(rx[i], ry[i])[0, 1] for i in (rng.integers(0, len(rx), len(rx)) for _ in range(n))]
    return est, *np.percentile(bs, [2.5, 97.5])


rows = []
g = el.dropna(subset=["p300"]).reset_index(drop=True)
g = g[g["erythroid_n"] > 0].reset_index(drop=True)
at, ae = pct(g["atac_xu"].values), pct(g["atac_encode"].values)
Z = np.column_stack([np.ones(len(g)), at, at ** 2, ae, ae ** 2, g["gc"].values, pct(g["erythroid_n"].values)])
for f in ("fp_tobias", "fp_print"):
    est, lo, hi = partial(pct(g[f].values), pct(g["p300"].values), Z)
    rows.append(dict(test=f"{f} vs p300, CRISPR elements with a site", n=len(g), partial_rho=est, ci_lo=lo, ci_hi=hi))
G2 = pd.DataFrame(rows)
print("\n=== 2 gate on CRISPR-tested elements (controls Xu + ENCODE ATAC, GC, motif count) ===\n"
      + G2.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 3 residual information beyond baseline
logit = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
base = logit(D["baseline"].values)
y = D["y"].values.astype(int)
ok = D["p300"].notna().values
p300_beyond = np.full(len(D), np.nan)
p300_beyond[ok] = resid(pct(D.loc[ok, "p300"].values),
                        np.column_stack([np.ones(ok.sum()), pct(D.loc[ok, "atac_xu"].values), pct(D.loc[ok, "atac_xu"].values) ** 2]))
D["p300_beyond_atac"] = p300_beyond
rows = []
for f in ("p300", "p300_beyond_atac", "fp_tobias", "fp_print", "erythroid_n", "atac_xu", "atac_encode", "gc", "CTCF_n"):
    m = D[f].notna().values
    x = pct(D.loc[m, f].values)
    X = np.column_stack([base[m], x])
    fit = lambda idx: LogisticRegression(C=1e6, max_iter=2000).fit(X[idx], y[m][idx]).coef_[0][1]
    est = fit(np.arange(m.sum()))
    bs = [fit(rng.integers(0, m.sum(), m.sum())) for _ in range(a.n_boot // 4)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append(dict(feature=f, n_pairs=int(m.sum()), coef_pct_beyond_baseline=est, ci_lo=lo, ci_hi=hi))
R3 = pd.DataFrame(rows)
print("\n=== 3 logistic coefficient of the feature's percentile beyond logit(baseline score) ===\n"
      + R3.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 4 rank shifts and AUPRC by distance
D["dbin"] = pd.cut(D["dist"], [0, 1e4, 5e4, 1e5, 1e7], labels=["<10 kb", "10-50 kb", "50-100 kb", ">100 kb"], right=False)
D["site"] = np.where(D["erythroid_n"] > 0, "erythroid site", "no site")
rk = {m: pct(D[m].values) for m in labs}
rows = []
for m in labs[1:]:
    d = rk[m] - rk["baseline"]
    for strat, grp in [("all", np.ones(len(D), bool))] + [(s, (D["site"] == s).values) for s in ("erythroid site", "no site")] \
            + [(str(b), (D["dbin"] == b).values) for b in D["dbin"].cat.categories]:
        rows.append(dict(model=m, stratum=strat, n_reg=int((grp & D["y"]).sum()),
                         shift_regulated=float(d[grp & D["y"].values].mean()),
                         shift_not_regulated=float(d[grp & ~D["y"].values].mean())))
S4 = pd.DataFrame(rows)
S4["net"] = S4["shift_regulated"] - S4["shift_not_regulated"]
print("\n=== 4 mean change in percentile rank vs baseline (x100) ===")
print((S4.set_index(["model", "stratum"])[["n_reg", "shift_regulated", "shift_not_regulated", "net"]]
       .assign(shift_regulated=lambda t: 100 * t.shift_regulated, shift_not_regulated=lambda t: 100 * t.shift_not_regulated,
               net=lambda t: 100 * t.net)).round(2).to_string(), flush=True)


def auprc(s, yy):
    p, r, _ = precision_recall_curve(yy, s)
    return auc(r, p)


rows = []
for b in D["dbin"].cat.categories:
    g = D[D["dbin"] == b]
    rows.append(dict(distance=str(b), n=len(g), n_reg=int(g["y"].sum()), **{m: auprc(g[m].values, g["y"].values) for m in labs}))
for cls, sel in (("promoter element", D["prom"] == 1), ("non-promoter element", D["prom"] == 0)):
    g = D[sel]
    rows.append(dict(distance=cls, n=len(g), n_reg=int(g["y"].sum()), **{m: auprc(g[m].values, g["y"].values) for m in labs}))
A5 = pd.DataFrame(rows)
print("\n=== AUPRC by distance bin and element class ===\n" + A5.round(3).to_string(index=False), flush=True)

for name, t in (("coverage", cov), ("gate_crispr_elements", G2), ("residual", R3), ("rank_shift", S4), ("auprc_by_distance", A5)):
    t.to_csv(f"{a.out_prefix}.{name}.tsv", sep="\t", index=False)
D.drop(columns=["dbin"]).to_csv(f"{a.out_prefix}.pairs.tsv.gz", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.*")
