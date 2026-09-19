#!/usr/bin/env python3
"""Predicted-activity profile by accessibility decile, straight from the bigwigs.

WHY THIS DOES NOT NEED ABC. F-021's finding is a property of the PREDICTED TRACKS, not of the
benchmark: the ratio of a transferred model's predicted activity to an in-cell model's rises
monotonically with accessibility, 0.64x at ATAC decile 1 to 2.16x at decile 10. Both arms read
the same ATAC track, so the ratio isolates what transfer changed. None of that requires ABC
scores, element-gene pairing or CRISPR labels, so the diagnostic can run the moment a predicted
bigwig exists and does not wait on an hours-long ABC chain.

WHAT IT IS FOR HERE. Testing whether standardizing the accessibility input on the prediction
regions narrows that swing. If the ratio profile flattens, the normalization mismatch was
driving the slope error; if it does not, the slope is learned behaviour and the input scale was
a red herring. Either answer is worth having before the benchmark returns, because the
benchmark can only say whether AUPRC moved, not why.

Deciles are on the accessibility track itself, over the same regions, so they are identical
for every arm by construction.

Usage:
  4.28.slope_profile.py --regions candidateRegions.bed --accessibility-bw atac.bw \\
      --track in_cell=predp300_multimodal.bw --track transferred=predp300gm_....bw \\
      --observed observed=p300_obs.bw --ratio transferred/in_cell
"""
import argparse
import numpy as np
import pandas as pd
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--regions", required=True)
ap.add_argument("--accessibility-bw", required=True)
ap.add_argument("--track", action="append", default=[], metavar="LABEL=PATH", required=True)
ap.add_argument("--ratio", action="append", default=[], metavar="A/B",
                help="report median(A)/median(B) per decile. Repeatable.")
ap.add_argument("--deciles", type=int, default=10)
ap.add_argument("--accessibility-sparse", action="store_true",
                help="the accessibility track holds one interval per region rather than dense "
                     "per-base coverage, so it can be read the fast way")
ap.add_argument("--sample", type=int, default=0, help="0 = all regions")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

TRACKS = dict(kv.split("=", 1) for kv in a.track)
reg = pd.read_csv(a.regions, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["chr", "start", "end"])


def per_region(path, dense=False):
    """Mean track value over each region. Regions absent from the bigwig read as NaN.

    TWO CODE PATHS, AND THE CHOICE IS NOT AN OPTIMISATION. A predicted track from 4.1 holds
    exactly one interval per candidate region, so reading a whole chromosome with intervals()
    is cheap and gives the stored value exactly. A raw accessibility track is dense per-base
    coverage: intervals() on one chromosome of it returns tens of millions of tuples and gets
    the process OOM-killed, which is how this function failed the first time it ran. For those,
    query per region with stats() instead, which touches only the bases asked for.
    """
    b = pyBigWig.open(path)
    ch = b.chroms()
    out = np.full(len(reg), np.nan)
    if dense:
        for i, c, s, e in zip(reg.index, reg["chr"], reg["start"].to_numpy(),
                              reg["end"].to_numpy()):
            if c not in ch or int(e) > ch[c]:
                continue
            v = b.stats(c, int(s), int(e), type="mean", exact=True)[0]
            if v is not None:
                out[i] = v
    else:
        for c, g in reg.groupby("chr"):
            if c not in ch:
                continue
            iv = b.intervals(c)
            if not iv:
                continue
            d = {(s, e): v for s, e, v in iv}
            for i, s, e in zip(g.index, g["start"].to_numpy(), g["end"].to_numpy()):
                v = d.get((int(s), int(e)))
                if v is not None:
                    out[i] = v
    b.close()
    return out


cols = {"atac": per_region(a.accessibility_bw, dense=not a.accessibility_sparse)}
for label, path in TRACKS.items():
    cols[label] = per_region(path)
df = pd.DataFrame(cols).dropna()
if a.sample and len(df) > a.sample:
    df = df.sample(a.sample, random_state=a.seed)
print(f"{len(df):,} regions with a value in every track\n")

df["dec"] = pd.qcut(df["atac"].rank(method="first"), a.deciles, labels=False) + 1
labels = list(TRACKS)
rows = []
for dec, g in df.groupby("dec"):
    r = {"dec": int(dec), "n": len(g), "med_atac": g["atac"].median()}
    for l in labels:
        r[f"med_{l}"] = g[l].median()
    for spec in a.ratio:
        x, y = spec.split("/")
        r[spec] = g[x].median() / max(g[y].median(), 1e-12)
    rows.append(r)
t = pd.DataFrame(rows)
print(t.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))

for spec in a.ratio:
    v = t[spec].to_numpy()
    print(f"\n{spec}: decile 1 {v[0]:.3f} -> decile {a.deciles} {v[-1]:.3f}, "
          f"swing x{v[-1] / max(v[0], 1e-12):.2f}")
    print(f"  {'monotone increasing' if np.all(np.diff(v) > -0.02) else 'not monotone'}; "
          f"spread max/min = x{v.max() / max(v.min(), 1e-12):.2f}")
