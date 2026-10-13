#!/usr/bin/env python3
"""Predicted p300 as scE2G evidence, per E-G pair: the deployable version of the F-056 ceiling.

WHY. Observed K562 p300 adds +0.056 to +0.060 CRISPR AUPRC to scE2G (F-056), mostly by demoting
accessible elements that are not active. This asks how much of that a PREDICTED p300 track keeps.

TRACKS (2026_0824_H3K27ac_model/data/abc_predicted, painted one value per July ENCODE candidate,
input ENCODE K562 ATAC, so not the Xu multiome the scE2G model sees):
  predp300mc    multi-cell model, K562 held out (deployable), mean of seeds 0, 1, 2 (F-038)
  predp300ic    in-cell model trained on K562 p300 (not deployable; upper bound for a prediction)
  atac_encode   ENCODE K562 ATAC at E (+/- 250 bp, from the 2.4 table): CONTROL. The predicted
                tracks see ENCODE ATAC, a deeper sample than Xu, so a gain must beat this arm to
                count as coming from the model rather than the second ATAC sample.
  ARCpredp300geo  ARC-E2G recomputed with activity sqrt(ATAC x predp300mc'), exact reweighting as
                4.9 / F-053 (abc_reweight.py); p300mc' = predp300mc + 1st percentile of positives
An element's predicted value is the mean of the painted track over its span; elements no July
candidate overlaps get 0 and are counted.

Usage: 4.10.make_predp300_pairs.py --out pairs.tsv.gz
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyBigWig

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_reweight import check, reweight  # noqa: E402

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
PRED = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/data/abc_predicted"
TRACKS = {"mc0": "predp300mc_multicell_holdK562.bw", "mc1": "predp300mc_multicell_holdK562_seed1.bw",
          "mc2": "predp300mc_multicell_holdK562_seed2.bw", "ic": "predp300_multimodal.bw"}
ap = argparse.ArgumentParser()
ap.add_argument("--arc", default="/scratch/users/shethm/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "ARC/EnhancerPredictionsAllPutative_ARC.tsv.gz")
ap.add_argument("--elements", default=f"{A2}/results/2026_1003_footprint_p300_gate.print.element_features.tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

E = pd.read_csv(a.elements, sep="\t", usecols=KEY + ["atac_encode", "p300"])
assert not E.duplicated(KEY).any()
for k, f in TRACKS.items():
    bw = pyBigWig.open(f"{PRED}/{f}")
    ch = bw.chroms()
    v = np.full(len(E), np.nan)
    for i, (c, s, e) in enumerate(zip(E["chr"].values, E["start"].values, E["end"].values)):
        if c in ch:
            x = bw.stats(c, int(s), int(min(e, ch[c])), type="mean", exact=True)[0]
            v[i] = np.nan if x is None else x
    bw.close()
    print(f"  {k}: {np.isfinite(v).mean():.4f} of {len(E):,} elements covered", flush=True)
    E[k] = np.nan_to_num(v)
E["predp300mc"] = E[["mc0", "mc1", "mc2"]].mean(1)
E["predp300ic"] = E["ic"]
for c in ("predp300mc", "predp300ic", "atac_encode"):
    r = pd.Series(E[c]).rank().corr(pd.Series(E["p300"]).rank())
    print(f"  {c}: spearman vs observed p300 over all elements {r:.3f}", flush=True)
print(f"  seed agreement (spearman mc0 vs mc1) {E['mc0'].rank().corr(E['mc1'].rank()):.3f}", flush=True)

P = pd.read_csv(a.arc, sep="\t", usecols=KEY + ["TargetGene", "isSelfPromoter", "activity_base",
                                                 "ABC.Score", "ARC.E2G.Score"])
P["isSelfPromoter"] = P["isSelfPromoter"].astype(str).str.lower() == "true"
assert not P.duplicated(KEY + ["TargetGene"]).any()
M = P[KEY].merge(E[KEY + ["predp300mc", "predp300ic", "atac_encode"]], on=KEY, how="left", validate="m:1")
assert M["predp300mc"].notna().mean() > 0.95, "element coordinates do not match the ARC table"
gene, self_, abc, arc = P["TargetGene"].values, P["isSelfPromoter"].values, P["ABC.Score"].values, P["ARC.E2G.Score"].values
check(gene, self_, abc, arc)

out = P[KEY + ["TargetGene"]].copy()
for c in ("predp300mc", "predp300ic", "atac_encode"):
    out[c] = M[c].fillna(0).values
p = out["predp300mc"].values
pc = float(np.quantile(p[p > 0], 0.01))
A = P["activity_base"].values
safe = np.where(A > 0, A, 1.0)
_, r_new = reweight(gene, self_, abc, arc, np.where(A > 0, np.sqrt((p + pc) / safe), 1.0))
out["ARCpredp300geo"] = r_new
nz = (arc > 0) & ~self_
print(f"  ARCpredp300geo: pseudocount {pc:.3g}; log vs original pearson "
      f"{np.corrcoef(np.log(r_new[nz]), np.log(arc[nz]))[0, 1]:.4f}", flush=True)
out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}: {len(out):,} pairs")
