#!/usr/bin/env python3
"""Repaint a predicted activity bigwig so its within-class ORDER follows a new score.

WHAT THIS IS FOR. 4.36 produces per-region scores (e.g. a held-out boosting fit that adds
neighbour terms) that are in percentile-rank units, not counts. ABC's rank qnorm keeps only
each region's within-class rank (F-022), so to run such a score through ABC it is enough to
give every region the predicted-track value sitting at its score's within-class quantile:
the track's value distribution is unchanged and only the ordering moves. This is 4.31's rank
transfer applied to every region instead of the tested ones.

Classes follow ABC: promoter = class in {tss, promoter}, nonpromoter = everything else.

Usage:
  4.37.build_reranked_activity.py --pred-bw PRED.bw \\
      --pred-enhancers PRED_ARM/Neighborhoods/EnhancerList.txt \\
      --scores nbr_tracks_multicell.tsv.gz --score-col boost_nbr --out-bw OUT.bw
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--pred-bw", required=True)
ap.add_argument("--pred-enhancers", required=True)
ap.add_argument("--scores", required=True)
ap.add_argument("--score-col", required=True)
ap.add_argument("--out-bw", required=True)
a = ap.parse_args()
KEY = ["chr", "start", "end"]

d = pd.read_csv(a.pred_enhancers, sep="\t")
rc = [c for c in d.columns if c.startswith("H3K27ac.") and c.endswith(".readCount")]
assert len(rc) == 1, f"expected one H3K27ac readCount column, got {rc}"
P = d[KEY + ["class", rc[0]]].rename(columns={rc[0]: "readCount"})
S = pd.read_csv(a.scores, sep="\t", usecols=KEY + [a.score_col])
m = P.merge(S, on=KEY, how="inner", validate="one_to_one")
assert len(m) == len(P) == len(S), f"region sets differ: {len(P):,} / {len(S):,} / {len(m):,}"
m["is_prom"] = m["class"].isin(["tss", "promoter"])

m["target"] = np.nan
for prom in (True, False):
    cls = m.is_prom == prom
    q = m.loc[cls, a.score_col].rank(method="first") / cls.sum()
    vals = np.sort(m.loc[cls, "readCount"].values)
    pos = np.clip(np.ceil(q * len(vals)).astype(int) - 1, 0, len(vals) - 1)
    m.loc[cls, "target"] = vals[pos]
assert m["target"].notna().all()

new = {(r.chr, int(r.start), int(r.end)): r.target / (r.end - r.start) for r in m.itertuples()}
src = pyBigWig.open(a.pred_bw)
out = pyBigWig.open(a.out_bw, "w")
chroms = src.chroms()
out.addHeader([(c, chroms[c]) for c in sorted(chroms)])
hit = 0
for c in sorted(chroms):
    iv = src.intervals(c) or []
    if not iv:
        continue
    v = []
    for st, en, val in iv:
        k = (c, st, en)
        if k in new:
            v.append(float(new[k]))
            hit += 1
        else:
            v.append(float(val))
    out.addEntries([c] * len(iv), [x[0] for x in iv], ends=[x[1] for x in iv], values=v)
out.close()
src.close()
assert hit == len(m), f"repainted {hit:,} of {len(m):,} regions: intervals do not line up"
print(f"repainted {hit:,} regions -> {a.out_bw}")
for prom in (True, False):
    cls = m.is_prom == prom
    print(f"  {'promoter' if prom else 'nonpromoter':<12} rho(new track, score) "
          f"{spearmanr(m.loc[cls, 'target'], m.loc[cls, a.score_col])[0]:.4f}   "
          f"rho(new, old) {spearmanr(m.loc[cls, 'target'], m.loc[cls, 'readCount'])[0]:.4f}")
