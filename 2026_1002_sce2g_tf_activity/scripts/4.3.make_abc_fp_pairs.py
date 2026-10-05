#!/usr/bin/env python3
"""Footprints inside the ABC activity term: footprint-weighted ABC and ARC-E2G scores per E-G pair.

Activity of element e becomes activity_e * w_e with w_e = 1 + lam * q_e, where q_e is the
percentile rank of the element's erythroid footprint sum (2.4 element table, TOBIAS or PRINT)
among elements with any erythroid site, and 0 for elements without one.

ABC_e = num_e / D_g with D_g the gene's denominator. Scaling D_g by sum(num w) / sum(num) over
the gene's candidate pairs (which keeps any extra terms scE2G puts in D_g) gives
    ABC_fp_e = ABC_e * w_e * sum_g(ABC) / sum_g(ABC * w)
so the numerators, which scE2G writes with 2-3 significant digits, are not needed. ARC-E2G is
ABC times a per-pair Kendall factor (checked below: log(ARC/ABC) is linear in Kendall), so the
same per-gene factor applies: ARC_fp_e = ARC_e * w_e * sum_g(ABC) / sum_g(ABC * w).

Writes one pairs file for scE2G's external-feature join (overlap on chr:TargetGene + coords).

Usage: 4.3.make_abc_fp_pairs.py --out pairs.tsv.gz [--lams 1,4]
"""
import argparse

import numpy as np
import pandas as pd

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
ap = argparse.ArgumentParser()
ap.add_argument("--arc", default="/scratch/users/shethm/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz")
ap.add_argument("--tobias", default=f"{A2}/results/2026_1003_footprint_p300_gate.element_features.tsv.gz")
ap.add_argument("--print", dest="print_", default=f"{A2}/results/2026_1003_footprint_p300_gate.print.element_features.tsv.gz")
ap.add_argument("--lams", default="1,4")
ap.add_argument("--out", required=True)
a = ap.parse_args()
lams = [float(x) for x in a.lams.split(",")]
KEY = ["chr", "start", "end"]

P = pd.read_csv(a.arc, sep="\t", usecols=KEY + ["TargetGene", "ABC.Score", "ARC.E2G.Score", "Kendall"])
print(f"{len(P):,} pairs, {P['TargetGene'].nunique():,} genes; "
      f"{(P['ABC.Score'] == 0).sum():,} with ABC = 0, {(P['ARC.E2G.Score'] == 0).sum():,} with ARC = 0", flush=True)
assert not P.duplicated(KEY + ["TargetGene"]).any()
assert (P["ABC.Score"] >= 0).all() and (P["ARC.E2G.Score"] >= 0).all()

k = P["Kendall"].notna() & (P["ABC.Score"] > 1e-3) & (P["ARC.E2G.Score"] > 0)
lr = np.log(P.loc[k, "ARC.E2G.Score"] / P.loc[k, "ABC.Score"])
fit = np.polyfit(P.loc[k, "Kendall"], lr, 1)
res = lr - np.polyval(fit, P.loc[k, "Kendall"])
print(f"log(ARC/ABC) ~ Kendall over {k.sum():,} pairs with ABC > 1e-3: slope {fit[0]:.3f}, "
      f"intercept {fit[1]:.3f}, R^2 {1 - res.var() / lr.var():.4f}", flush=True)

gene = P["TargetGene"].values
abc_sum = P.groupby("TargetGene")["ABC.Score"].transform("sum").values
out = P[KEY + ["TargetGene"]].copy()
for src, path in [("tobias", a.tobias), ("print", a.print_)]:
    E = pd.read_csv(path, sep="\t", usecols=KEY + ["erythroid_n", "erythroid_fp_sum"])
    assert not E.duplicated(KEY).any()
    has = E["erythroid_n"] > 0
    q = np.zeros(len(E))
    q[has.values] = E.loc[has, "erythroid_fp_sum"].rank(pct=True).values
    E["q"] = q
    M = P[KEY].merge(E[KEY + ["q"]], on=KEY, how="left", validate="m:1")
    print(f"{src}: {has.sum():,} of {len(E):,} elements with an erythroid site; "
          f"{M['q'].notna().mean():.4f} of pairs matched", flush=True)
    assert M["q"].notna().mean() > 0.95
    qq = M["q"].fillna(0).values
    for lam in lams:
        w = 1 + lam * qq
        aw = P["ABC.Score"].values * w
        den = pd.Series(aw).groupby(gene).transform("sum").values
        g = np.divide(abc_sum, den, out=np.ones_like(den), where=den > 0)  # all-zero genes stay 0
        tag = f"{src}_l{lam:g}"
        out[f"ABCfp_{tag}"] = aw * g
        out[f"ARCfp_{tag}"] = P["ARC.E2G.Score"].values * w * g
        s_new = out.groupby("TargetGene")[f"ABCfp_{tag}"].sum()
        s_old = P.groupby("TargetGene")["ABC.Score"].sum()
        assert np.allclose(s_new.values, s_old.loc[s_new.index].values, rtol=1e-6), f"{tag}: per-gene ABC sum changed"
        nz = P["ARC.E2G.Score"].values > 0
        r = np.corrcoef(np.log(out.loc[nz, f"ARCfp_{tag}"]), np.log(P.loc[nz, "ARC.E2G.Score"]))[0, 1]
        print(f"  {tag}: log ARC_fp vs log ARC pearson {r:.4f}; ARC_fp / ARC "
              f"{pd.Series(out.loc[nz, f'ARCfp_{tag}'] / P.loc[nz, 'ARC.E2G.Score']).describe().round(3).to_dict()}", flush=True)

out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}")
