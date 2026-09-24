#!/usr/bin/env python3
"""Oracle activity track: observed-p300 RANKS at CRISPR-tested elements, predicted values
everywhere else.

WHAT THIS IS FOR. The multi-cell arm scores 0.4810 AUPRC and observed p300 scores 0.5673. An
ABC score is activity x contact over the SUM of activity x contact across every candidate near
the gene, so the gap can come from two places: the tested elements themselves (the numerator)
or the ~97% of candidates that were never tested (the denominator, and the reference
distribution the rank qnorm uses). This track fixes only the first. Scored alongside the
all-predicted and all-observed arms it splits the gap:

    B - A  = AUPRC from ranking the tested elements correctly
    C - B  = AUPRC from ranking everything else correctly

It is an oracle and is not deployable; it probes where the loss is.

WHY RANKS AND NOT VALUES. ABC's rank qnorm (neighborhoods.run_qnorm) maps each element's
within-class rank onto a reference distribution and discards the input's scale (F-022).
Pasting observed p300 values into a predicted track would give the tested elements ranks set
by the SCALE MISMATCH between the two tracks, not by their observed ordering. So each tested
region gets the predicted-track value that sits at its observed within-class quantile. Classes
follow ABC exactly: promoter = class in {tss, promoter}, nonpromoter = everything else.

HOW IT PAINTS. 4.1 paints each candidate region with a constant predicted_counts / width, so
ABC's readCount over the region is the painted value times the width. Tested regions are
repainted with target_readCount / width; every other interval is copied unchanged.

Usage:
  4.31.build_oracle_activity.py --pred-bw PRED.bw \\
      --pred-enhancers PRED_ARM/Neighborhoods/EnhancerList.txt \\
      --obs-enhancers  OBS_ARM/Neighborhoods/EnhancerList.txt \\
      --crispr-merged expt_pred_merged_annot.txt.gz --out-bw ORACLE.bw [--report-tsv R.tsv]
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--pred-bw", required=True)
ap.add_argument("--pred-enhancers", required=True)
ap.add_argument("--obs-enhancers", required=True)
ap.add_argument("--crispr-merged", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--report-tsv", default=None)
a = ap.parse_args()

KEY = ["chr", "start", "end"]


def load(path, need_count):
    # Only the predicted arm needs a readCount, since that is what gets repainted. The
    # observed arm contributes ranks only, taken from ABC's combined H3K27ac.RPM, which is
    # what ABC itself ranks when an arm has several replicate files.
    d = pd.read_csv(path, sep="\t")
    cols = KEY + ["class", "H3K27ac.RPM"]
    if not need_count:
        return d[cols]
    rc = [c for c in d.columns if c.startswith("H3K27ac.") and c.endswith(".readCount")]
    assert len(rc) == 1, f"{path}: expected one H3K27ac readCount column, got {rc}"
    return d[cols + [rc[0]]].rename(columns={rc[0]: "readCount"})


P = load(a.pred_enhancers, need_count=True).rename(columns={"readCount": "readCount_p"})
O = load(a.obs_enhancers, need_count=False)
m = P.merge(O, on=KEY, suffixes=("_p", "_o"), how="inner", validate="one_to_one")
assert len(m) == len(P) == len(O), (
    f"region sets differ: pred {len(P):,}, obs {len(O):,}, shared {len(m):,}. "
    "The arms must share one candidate region set.")
assert (m["class_p"] == m["class_o"]).all(), "class labels disagree between arms"
m["is_prom"] = m["class_p"].isin(["tss", "promoter"])
print(f"{len(m):,} candidate regions, {int(m.is_prom.sum()):,} promoter-class")

# ---------------------------------------------------------------- tested regions
cr = pd.read_csv(a.crispr_merged, sep="\t", usecols=["chrom", "chromStart", "chromEnd"])
cr = cr.drop_duplicates()
m["tested"] = False
for c, g in m.groupby("chr"):
    e = cr[cr.chrom == c]
    if e.empty:
        continue
    s = g["start"].values
    en = g["end"].values
    order = np.argsort(s)
    s, en, idx = s[order], en[order], g.index.values[order]
    for es, ee in zip(e.chromStart.values, e.chromEnd.values):
        # candidate regions do not overlap each other, so the overlapping set is contiguous
        lo = np.searchsorted(en, es, side="right")
        hi = np.searchsorted(s, ee, side="left")
        if hi > lo:
            m.loc[idx[lo:hi], "tested"] = True
n_el = len(cr)
n_reg = int(m.tested.sum())
print(f"{n_el:,} unique tested elements overlap {n_reg:,} candidate regions "
      f"({100 * n_reg / len(m):.2f}% of all regions)")

# ---------------------------------------------------------------- rank transfer
m["target"] = m["readCount_p"].astype(float)
for prom in (True, False):
    cls = m.is_prom == prom
    q_obs = m.loc[cls, "H3K27ac.RPM_o"].rank() / cls.sum()     # ABC's own quantile rule
    pred_sorted = np.sort(m.loc[cls, "readCount_p"].values)
    t = cls & m.tested
    # value in the predicted within-class distribution at the observed within-class quantile
    pos = np.clip(np.ceil(q_obs[t] * len(pred_sorted)).astype(int) - 1, 0,
                  len(pred_sorted) - 1)
    m.loc[t, "target"] = pred_sorted[pos]
    print(f"  {'promoter' if prom else 'nonpromoter':<12} tested {int(t.sum()):,} of "
          f"{int(cls.sum()):,}")

# ---------------------------------------------------------------- repaint
new = {(r.chr, int(r.start), int(r.end)): r.target / (r.end - r.start)
       for r in m[m.tested].itertuples()}
src = pyBigWig.open(a.pred_bw)
out = pyBigWig.open(a.out_bw, "w")
chroms = src.chroms()
out.addHeader([(c, chroms[c]) for c in sorted(chroms)])
hit = 0
for c in sorted(chroms):
    iv = src.intervals(c) or []
    if not iv:
        continue
    s = [x[0] for x in iv]
    e = [x[1] for x in iv]
    v = []
    for st, en, val in iv:
        k = (c, st, en)
        if k in new:
            v.append(float(new[k]))
            hit += 1
        else:
            v.append(float(val))
    out.addEntries([c] * len(s), s, ends=e, values=v)
out.close()
src.close()
assert hit == n_reg, (
    f"repainted {hit:,} intervals but {n_reg:,} tested regions: the bigwig's intervals do not "
    "line up with the EnhancerList regions, so the repaint is incomplete")
print(f"repainted {hit:,} intervals -> {a.out_bw}")

# ---------------------------------------------------------------- self-check
# Within class, the tested regions should now carry the OBSERVED ordering and the untested
# ones the PREDICTED ordering. Checked in rank space because that is all qnorm keeps.
rows = []
for prom in (True, False):
    cls = m.is_prom == prom
    t, u = cls & m.tested, cls & ~m.tested
    rows.append({
        "class": "promoter" if prom else "nonpromoter",
        "n_tested": int(t.sum()),
        "rho_tested_vs_obs": spearmanr(m.loc[t, "target"], m.loc[t, "H3K27ac.RPM_o"])[0],
        "rho_tested_pred_vs_obs_before": spearmanr(m.loc[t, "readCount_p"],
                                                   m.loc[t, "H3K27ac.RPM_o"])[0],
        "rho_untested_unchanged": spearmanr(m.loc[u, "target"], m.loc[u, "readCount_p"])[0],
    })
rep = pd.DataFrame(rows)
print(rep.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
if a.report_tsv:
    rep.assign(n_elements=n_el, n_regions_touched=n_reg).to_csv(a.report_tsv, sep="\t",
                                                               index=False)
    print("wrote", a.report_tsv)
