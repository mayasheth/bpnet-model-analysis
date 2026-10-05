#!/usr/bin/env python3
"""Expand per-element footprint features (2.4 element table) to scE2G E-G pairs.

scE2G's merge_external_features.R joins external features to E-G pairs by overlap on
chr:TargetGene + coordinates, so an element feature has to be written once per pair: the pairs
come from ABC's EnhancerPredictionsAllPutative for the Xu run, joined to the 2.4 element table on
exact coordinates (both are scE2G's candidate regions). Elements not in the table (promoter-only
chromosomes TOBIAS did not score) keep NA and are filled with the feature table's fill value.

Usage: 4.1.make_element_feature_pairs.py --elements E.tsv.gz --out pairs.tsv.gz [--cols a,b,...]
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--elements", required=True)
ap.add_argument("--predictions", default="/scratch/users/shethm/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "Predictions/EnhancerPredictionsAllPutative.tsv.gz")
ap.add_argument("--cols", default="erythroid_n,erythroid_fp_sum,erythroid_fp_max,CTCF_n,CTCF_fp_max")
ap.add_argument("--out", required=True)
a = ap.parse_args()
cols = a.cols.split(",")
KEY = ["chr", "start", "end"]

E = pd.read_csv(a.elements, sep="\t", usecols=KEY + cols)
assert not E.duplicated(KEY).any()
P = pd.read_csv(a.predictions, sep="\t", usecols=KEY + ["TargetGene"])
M = P.merge(E, on=KEY, how="left")
assert len(M) == len(P)
hit = M[cols[0]].notna()
print(f"{len(P):,} pairs, {P[KEY].drop_duplicates().shape[0]:,} elements; "
      f"{hit.mean():.4f} of pairs matched an element row", flush=True)
assert hit.mean() > 0.95, "element coordinates do not match the predictions"
M.to_csv(a.out, sep="\t", index=False)
print(M[cols].describe().round(3).to_string())
print(f"wrote {a.out}")
