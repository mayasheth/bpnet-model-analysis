#!/usr/bin/env python3
"""Element features for the systematic-correction scE2G run (from the 6.1 error analysis), expanded
to E-G pairs like 4.1.

  ctcf_fp_print  PRINT CTCF footprint max in the element (0 without a CTCF site); 6.1: CTCF sites
                 carry negative CRISPR information beyond the baseline score
  ctcf_n         CTCF motif count, to separate motif presence from footprint
  gc             element GC fraction (+/- 250 bp); 6.1: negative beyond the baseline
  fp_pct1        1 + percentile of the PRINT erythroid footprint sum among elements with a site,
                 1 without one. scE2G takes log(|x| + 0.01) of every feature, so a raw sum that is 0
                 at 39% of elements becomes a site / no-site step; 1 + percentile removes the step.

Usage: 4.7.make_correction_pairs.py --out pairs.tsv.gz
"""
import argparse

import numpy as np
import pandas as pd

A2 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity"
ap = argparse.ArgumentParser()
ap.add_argument("--elements", default=f"{A2}/results/2026_1003_footprint_p300_gate.print.element_features.tsv.gz")
ap.add_argument("--predictions", default="/scratch/users/shethm/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "Predictions/EnhancerPredictionsAllPutative.tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

E = pd.read_csv(a.elements, sep="\t", usecols=KEY + ["gc", "CTCF_n", "CTCF_fp_max", "erythroid_n", "erythroid_fp_sum"])
assert not E.duplicated(KEY).any()
has = E["erythroid_n"] > 0
q = np.zeros(len(E))
q[has.values] = E.loc[has, "erythroid_fp_sum"].rank(pct=True).values
E = E.assign(fp_pct1=1 + q, ctcf_fp_print=E["CTCF_fp_max"], ctcf_n=E["CTCF_n"])
cols = ["ctcf_fp_print", "ctcf_n", "gc", "fp_pct1"]
print(E[cols].describe().round(3).to_string(), flush=True)
print(f"CTCF site in {(E['ctcf_n'] > 0).mean():.3f} of elements; erythroid site in {has.mean():.3f}", flush=True)
P = pd.read_csv(a.predictions, sep="\t", usecols=KEY + ["TargetGene"])
M = P.merge(E[KEY + cols], on=KEY, how="left")
assert len(M) == len(P)
assert M[cols[0]].notna().mean() > 0.95, "element coordinates do not match the predictions"
M.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}: {len(M):,} pairs")
