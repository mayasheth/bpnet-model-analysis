#!/usr/bin/env python3
"""Per E-G pair activity-term variants for scE2G, with the exact ABC reweighting (abc_reweight.py).
Supersedes 4.3 (whose per-gene factor treated the self-promoter's ABC = 1 as a share) and 4.5.

Columns (each an ARC-E2G score recomputed with a new element activity, plus raw element p300):
  ARCfp_{tobias,print}_l{1,4}  activity x (1 + lam x erythroid footprint percentile), as F-053
  ARCp300geo                   activity sqrt(ATAC x p300'): w = sqrt(p300' / A)
  ARCp300only                  activity p300': w = p300' / A
  p300_obs                     element p300 (ENCSR000EGE, +/- 250 bp), for the extra-feature model
A = activity_base (= the ABC numerator / powerlaw_contact, checked in abc_columns 2026-10-04);
p300' = p300 + the 1st percentile of positive element p300.

Usage: 4.9.make_activity_pairs.py --out pairs.tsv.gz
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_reweight import check, reweight  # noqa: E402

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
ap = argparse.ArgumentParser()
ap.add_argument("--arc", default="/scratch/users/shethm/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

P = pd.read_csv(a.arc, sep="\t", usecols=KEY + ["TargetGene", "isSelfPromoter", "activity_base",
                                                 "ABC.Score", "ARC.E2G.Score"])
P["isSelfPromoter"] = P["isSelfPromoter"].astype(str).str.lower() == "true"
assert not P.duplicated(KEY + ["TargetGene"]).any()
gene, self_, abc, arc = P["TargetGene"].values, P["isSelfPromoter"].values, P["ABC.Score"].values, P["ARC.E2G.Score"].values
s_ns = check(gene, self_, abc, arc)
print(f"{len(P):,} pairs, {P['TargetGene'].nunique():,} genes, {self_.sum():,} self-promoter pairs (ABC = 1); "
      f"non-self ABC sum per gene: median {s_ns.median():.3f}, min {s_ns.min():.3f}", flush=True)

out = P[KEY + ["TargetGene"]].copy()


def add(tag, w):
    a_new, r_new = reweight(gene, self_, abc, arc, w)
    out[f"ARC{tag}"] = r_new
    nz = (arc > 0) & ~self_
    r = np.corrcoef(np.log(r_new[nz]), np.log(arc[nz]))[0, 1]
    print(f"  ARC{tag}: log vs original pearson {r:.4f}; ratio quantiles "
          f"{np.round(np.quantile(r_new[nz] / arc[nz], [0.05, 0.5, 0.95]), 3).tolist()}", flush=True)


for src, f in (("tobias", "2026_1003_footprint_p300_gate.element_features.tsv.gz"),
               ("print", "2026_1003_footprint_p300_gate.print.element_features.tsv.gz")):
    E = pd.read_csv(f"{A2}/results/{f}", sep="\t", usecols=KEY + ["erythroid_n", "erythroid_fp_sum", "p300"])
    has = E["erythroid_n"] > 0
    q = np.zeros(len(E))
    q[has.values] = E.loc[has, "erythroid_fp_sum"].rank(pct=True).values
    E["q"] = q
    M = P[KEY].merge(E, on=KEY, how="left", validate="m:1")
    assert M["q"].notna().mean() > 0.95
    for lam in (1, 4):
        add(f"fp_{src}_l{lam}", 1 + lam * M["q"].fillna(0).values)
    if src == "tobias":
        p = M["p300"].fillna(0).values
        pc = float(E.loc[E["p300"] > 0, "p300"].quantile(0.01))
        A = P["activity_base"].values
        safe = np.where(A > 0, A, 1.0)
        add("p300geo", np.where(A > 0, np.sqrt((p + pc) / safe), 1.0))
        add("p300only", np.where(A > 0, (p + pc) / safe, 1.0))
        out["p300_obs"] = p
        print(f"  p300 pseudocount {pc:.3g}; element p300 quantiles {E['p300'].quantile([0.1, 0.5, 0.9]).round(1).tolist()}", flush=True)
out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}")
