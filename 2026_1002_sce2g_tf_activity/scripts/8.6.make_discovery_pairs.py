#!/usr/bin/env python3
"""8.6 Phase 4 E-G pairs (plans/8_crispr_feature_discovery.md; pre-registered in .living/decisions.md
[2026-10-07]): the three features picked in >= 4 of 5 outer folds by 8.3 forward selection, which beat the
8.4 feature-shuffle null and survived log TSS distance in the base.

  cgi_logdist_p1  1 + percentile of log10 distance from the element midpoint to the nearest CpG island
  gata_n_p1       1 + percentile of the GATA cluster count (max of GATA6 MA1104.3, Gata3 MA0037.5, TRPS1
                  MA1970.2 sites within +/- 250 bp; 8.2 motif cluster motif:GATA6_MA1104.3+2)
  ikzf2_n_p1      1 + percentile of IKZF2 MA2326.1 sites within +/- 250 bp
  gc              raw GC fraction (+/- 250 bp), the F-055 feature, for the GC arm
Percentiles over all 155,134 Xu elements, exactly as 8.2 encoded them (NaN -> lowest).

Usage: 8.6.make_discovery_pairs.py --out pairs.tsv.gz
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SC = os.environ["SCRATCH"]
LIB = f"{SC}/feature_screen"
ap = argparse.ArgumentParser()
ap.add_argument("--predictions", default=f"{SC}/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/"
                "Predictions/EnhancerPredictionsAllPutative.tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

E = pd.read_csv(f"{LIB}/library.elements.tsv.gz", sep="\t", usecols=KEY + ["cgi_log_dist", "gc"])
M = np.load(f"{LIB}/library.motifs.npz", allow_pickle=True)
assert (M["chr"] == E["chr"].values).all() and (M["start"] == E["start"].values).all()
names = [str(x) for x in M["names"]]
col = lambda n: M["n"][:, names.index(n)].astype(float)  # noqa: E731
gata = np.max([col(n) for n in ("GATA6_MA1104.3", "Gata3_MA0037.5", "TRPS1_MA1970.2")], axis=0)
pct = lambda v: 1 + rankdata(np.nan_to_num(np.asarray(v, float), nan=-np.inf)) / len(v)  # noqa: E731
E["cgi_logdist_p1"] = pct(E["cgi_log_dist"])
E["gata_n_p1"] = pct(gata)
E["ikzf2_n_p1"] = pct(col("IKZF2_MA2326.1"))
cols = ["cgi_logdist_p1", "gata_n_p1", "ikzf2_n_p1", "gc"]
print(f"GATA cluster site in {(gata > 0).mean():.3f} of elements; IKZF2 in {(col('IKZF2_MA2326.1') > 0).mean():.3f}", flush=True)
print(E[cols].describe().round(3).to_string(), flush=True)
P = pd.read_csv(a.predictions, sep="\t", usecols=KEY + ["TargetGene"])
out = P.merge(E[KEY + cols], on=KEY, how="left")
assert len(out) == len(P) and out["gc"].notna().mean() > 0.95
out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out}: {len(out):,} pairs, {out['gc'].notna().mean():.4f} with features")
