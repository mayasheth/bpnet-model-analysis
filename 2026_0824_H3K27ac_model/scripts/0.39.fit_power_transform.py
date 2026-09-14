#!/usr/bin/env python
"""Fit a power transform that un-compresses the painted converter track's dynamic range.

WHY NOT QUANTILE MAPPING. That was tried (0.37) and made things worse: it fixed the marginal
distribution exactly, as designed, and dynamic range fell from 18x to 9x against real DNase's
39x while 1 bp top-quintile shape dropped 0.779 to 0.659. Quantile mapping can only fix a
track's marginal; it cannot fix its ORDERING, and the painted track's rank correlation with
real DNase is 0.827. Forcing the marginal to match while the ordering is wrong hands reference
mass to over-ranked quiet bases, and that mass comes out of the strong elements.

WHAT THIS DOES INSTEAD. `v -> scale * v^gamma` with gamma > 1. Big values grow faster than
small ones, so strong and weak elements are pulled apart. It is strictly monotone, so the
ordering, which is the part F-011 established and the part this cannot repair anyway, is
untouched.

FITTED IN GM12878, APPLIED TO K562, WHICH IS WHAT MAKES IT DEPLOYMENT-LEGAL. GM12878 has both
ATAC and DNase, so gamma can be chosen there by comparing a painted track against real DNase.
Choosing gamma on K562's own DNase would leak the assay the converter exists to avoid needing.
Fitting on one cell type and applying to another is exactly the deployment story: you have a
converter and someone else's paired data, and no DNase in your own cell type.

THE OBJECTIVE IS DYNAMIC RANGE, NOT CORRELATION. Specifically the ratio of mean signal in the
top observed-DNase quintile to that in the bottom, which is the quantity F-013's failure points
at. Shape correlation is reported alongside as a guard rail: a gamma that expands range while
destroying shape is not a fix, and the printout shows both so the trade is visible rather than
assumed.

Usage:
  0.39.fit_power_transform.py --painted-glob 'parts/gm12878_*.bw' \\
      --observed-bw gm12878_dnase_5p.bw --elements ELEMENTS.narrowPeak
"""
import argparse, glob
import numpy as np
import pandas as pd
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--painted-glob", required=True)
ap.add_argument("--observed-bw", required=True)
ap.add_argument("--elements", required=True)
ap.add_argument("--chrom", default="chr8", help="one chromosome is enough to fit one scalar")
ap.add_argument("--gammas", type=float, nargs="+",
                default=[1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0])
ap.add_argument("--half", type=int, default=500)
a = ap.parse_args()

part = None
for f in sorted(glob.glob(a.painted_glob)):
    b = pyBigWig.open(f)
    ch = list(b.chroms())
    b.close()
    if ch == [a.chrom]:
        part = f
        break
if part is None:
    raise SystemExit(f"no painted part for {a.chrom} matched {a.painted_glob}")
print(f"painted: {part}\nobserved: {a.observed_bw}")

els = pd.read_csv(a.elements, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["chr", "start", "end"])
els = els[els["chr"] == a.chrom].reset_index(drop=True)
bp, bo = pyBigWig.open(part), pyBigWig.open(a.observed_bw)
size = bp.chroms()[a.chrom]
Pv, Ov = [], []
for _, r in els.iterrows():
    c = (int(r["start"]) + int(r["end"])) // 2
    s, e = c - a.half, c + a.half
    if s < 0 or e > size:
        continue
    pv = bp.values(a.chrom, s, e, numpy=True)
    ov = bo.values(a.chrom, s, e, numpy=True)
    if pv is None or ov is None or len(pv) != 2 * a.half or len(ov) != 2 * a.half:
        continue
    Pv.append(np.nan_to_num(pv, nan=0.0))
    Ov.append(np.nan_to_num(ov, nan=0.0))
bp.close(); bo.close()
Pv, Ov = np.stack(Pv), np.stack(Ov)
print(f"{len(Pv):,} {a.chrom} elements\n")


def shape_corr(x, y, b=1):
    n, L = x.shape
    xb = x.reshape(n, L // b, b).sum(2)
    yb = y.reshape(n, L // b, b).sum(2)
    xc = xb - xb.mean(1, keepdims=True)
    yc = yb - yb.mean(1, keepdims=True)
    num = (xc * yc).sum(1)
    den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1))
    good = den > 0
    return float((num[good] / den[good]).mean())


tot_o = Ov.sum(1)
q = np.quantile(tot_o, [0.2, 0.8])
lo, hi = tot_o < q[0], tot_o >= q[1]
target_dr = Ov[hi].sum(1).mean() / max(Ov[lo].sum(1).mean(), 1e-9)
print(f"target dynamic range (observed DNase, Q5/Q1) = {target_dr:.1f}x\n")
print(f"{'gamma':>6}{'Q5/Q1':>9}{'% of target':>13}{'1bp topq shape':>16}{'1bp all':>10}")

best = None
for g in a.gammas:
    X = Pv ** g
    # Rescale so the element-level total matches the observed, which is what training's own
    # normalisation would do anyway; it keeps the numbers readable and cannot affect the ratio.
    sc = Ov.sum() / max(X.sum(), 1e-9)
    X = X * sc
    dr = X[hi].sum(1).mean() / max(X[lo].sum(1).mean(), 1e-9)
    sh_t = shape_corr(Ov[hi], X[hi])
    sh_a = shape_corr(Ov, X)
    print(f"{g:>6.1f}{dr:>9.1f}{100 * dr / target_dr:>12.0f}%{sh_t:>16.4f}{sh_a:>10.4f}")
    # Closest to the target range, but refuse a gamma that costs more than 0.02 of the
    # top-quintile shape correlation: shape is what the converter is FOR.
    if sh_t >= shape_corr(Ov[hi], Pv[hi]) - 0.02:
        score = abs(np.log(dr / target_dr))
        if best is None or score < best[1]:
            best = (g, score, dr, sh_t)

print()
if best is None:
    print("NO gamma keeps top-quintile shape within 0.02 of untransformed. The power transform\n"
          "cannot expand range here without paying more in shape than the converter can afford.")
else:
    g, _, dr, sh = best
    print(f"CHOSEN gamma = {g:.1f}: dynamic range {dr:.1f}x against the observed {target_dr:.1f}x "
          f"({100 * dr / target_dr:.0f}%), top-quintile shape {sh:.4f}")
    print(f"Apply with 0.40 to the K562 painted parts. gamma was fitted HERE, in GM12878, so no\n"
          f"K562 DNase enters the choice.")
