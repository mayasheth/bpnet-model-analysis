#!/usr/bin/env python3
"""7.3 E-G pairs for the scE2G confirmation of the 7.2 feature screen (plans/7_feature_screen.md;
pass rule in .living/decisions.md, 2026-10-05).

  gc               element GC fraction (+/- 250 bp), the F-055 feature: every arm carries it
  frag_short_frac  share of Xu insertions from fragments < 100 bp (+/- 250 bp); the one shortlisted
                   feature; 0 where an element has no insertions
  screen_comp1     1 + percentile of the 7.2 out-of-fold composite (trees fitted on the other
                   chromosome half). The composite is on a signed rank-residual scale, and scE2G takes
                   log(|x| + 0.01), which would fold negative values onto positive ones.

Usage: 7.3.make_screen_pairs.py --out pairs.tsv.gz
"""
import argparse
import os

import pandas as pd

SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--lib-dir", default=f"{SC}/feature_screen")
ap.add_argument("--predictions", default=f"{SC}/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "Predictions/EnhancerPredictionsAllPutative.tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

E = pd.read_csv(f"{a.lib_dir}/library.elements.tsv.gz", sep="\t", usecols=KEY + ["gc", "frag_short_frac"])
C = pd.read_csv(f"{a.lib_dir}/screen_composite.tsv.gz", sep="\t")
assert C["screen_composite"].notna().all(), "composite missing for some elements"
E = E.merge(C, on=KEY, how="inner", validate="1:1")
assert len(E) == len(C)
E["frag_short_frac"] = E["frag_short_frac"].fillna(0)
E["screen_comp1"] = 1 + E["screen_composite"].rank(pct=True)
cols = ["gc", "frag_short_frac", "screen_comp1"]
print(E[cols].describe().round(3).to_string(), flush=True)
print(E[cols].corr("spearman").round(3).to_string(), flush=True)
P = pd.read_csv(a.predictions, sep="\t", usecols=KEY + ["TargetGene"])
M = P.merge(E[KEY + cols], on=KEY, how="left")
assert len(M) == len(P)
assert M["gc"].notna().mean() > 0.95, "element coordinates do not match the predictions"
M.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}: {len(M):,} pairs, {M['gc'].notna().mean():.4f} with features")
