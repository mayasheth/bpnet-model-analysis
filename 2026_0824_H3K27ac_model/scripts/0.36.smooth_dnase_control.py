#!/usr/bin/env python
"""The control track: real DNase with its base-resolution structure erased, magnitude kept.

WHAT THIS IS FOR. F-010 showed a DNase accessibility input beats an ATAC one for predicting
H3K27ac, in-cell and on transfer. Nobody has measured WHICH property of DNase wins. Two
candidates are tangled together in that result:

  magnitude   how much DNase signal sits in and around the element
  shape       the base-resolution cut-site structure inside it

The converter is excellent at shape (95% of the inter-replicate ceiling at 1 bp) and
mediocre at magnitude (its counts bought nothing in the ABC benchmark). So the two
candidates have OPPOSITE implications for whether the converter can ever work, and a
converter arm on its own cannot distinguish them: a negative result would be equally
consistent with "shape is irrelevant" and "this converter's shape is not good enough".

This track breaks the tie. A box filter of width W replaces each base with the local mean,
so the total signal over any region much larger than W is preserved while all structure
finer than W is destroyed. Train the H3K27ac model on it and:

  advantage over ATAC SURVIVES   -> shape never mattered; the converter's strength is worth
                                    nothing downstream and better counts is the only path.
  advantage over ATAC COLLAPSES  -> shape is what DNase was winning on, which is precisely
                                    what the converter is good at.

WHY W = 250 AND NOT AN ARBITRARY NUMBER. 250 bp is the scale at which ATAC and DNase already
agree about shape: the observed ATAC-vs-DNase within-element correlation is 0.915 at 250 bp
top quintile, against 0.292 at 1 bp (0.35). So smoothing DNase at 250 bp removes close to
exactly the information that ATAC does not already carry, which is the information the
converter claims to supply. A wider filter would also blur across element boundaries and
start destroying magnitude, which is the thing being held fixed.

THE FILTER PRESERVES THE TOTAL. A box filter is mean-preserving away from the edges, so the
genome-wide sum is conserved to within the filter width at each contig end. That is asserted
at the end rather than assumed: if the sum moved, the control would differ from real DNase in
magnitude too and would no longer isolate shape.

Usage: 0.36.smooth_dnase_control.py --in-bw IN.bw --out-bw OUT.bw --chrom-sizes SIZES [-w 250]
"""
import argparse
import numpy as np
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--in-bw", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("-w", "--width", type=int, default=250)
ap.add_argument("--chroms", default=None)
a = ap.parse_args()

sizes = {}
for line in open(a.chrom_sizes):
    p = line.split()
    if len(p) >= 2:
        sizes[p[0]] = int(p[1])

src = pyBigWig.open(a.in_bw)
have = src.chroms()
chroms = [c for c in sizes if c in have]
if a.chroms:
    keep = set(a.chroms.split(","))
    chroms = [c for c in chroms if c in keep]
chroms.sort()

out = pyBigWig.open(a.out_bw, "w")
out.addHeader([(c, sizes[c]) for c in chroms])

W = a.width
sum_in = sum_out = 0.0
for chrom in chroms:
    n = sizes[chrom]
    x = np.nan_to_num(src.values(chrom, 0, n, numpy=True), nan=0.0).astype(np.float64)
    sum_in += x.sum()
    # Box filter by cumulative sum, which is exact and O(n) rather than an FFT approximation.
    c = np.concatenate([[0.0], np.cumsum(x)])
    half = W // 2
    lo = np.clip(np.arange(n) - half, 0, n)
    hi = np.clip(np.arange(n) - half + W, 0, n)
    y = (c[hi] - c[lo]) / np.maximum(hi - lo, 1)
    sum_out += y.sum()
    del x, c, lo, hi

    nz = y > 0
    if nz.any():
        idx = np.flatnonzero(nz)
        v = y[idx].astype(np.float32)
        brk = np.flatnonzero(np.diff(idx) != 1) + 1
        seg = np.concatenate([[0], brk, [len(idx)]])
        cs, ce, cv = [], [], []
        for b0, b1 in zip(seg[:-1], seg[1:]):
            ii, vv = idx[b0:b1], v[b0:b1]
            ch = np.flatnonzero(np.diff(vv) != 0) + 1
            rb = np.concatenate([[0], ch, [len(vv)]])
            for r0, r1 in zip(rb[:-1], rb[1:]):
                cs.append(int(ii[r0])); ce.append(int(ii[r1 - 1]) + 1)
                cv.append(float(vv[r0]))
        out.addEntries([chrom] * len(cs), cs, ends=ce, values=cv)
    print(f"  {chrom}: n={n:,} covered={int(nz.sum()):,}", flush=True)
    del y, nz

src.close()
out.close()
r = sum_out / sum_in if sum_in else float("nan")
print(f"\nsum in  {sum_in:,.0f}\nsum out {sum_out:,.0f}\nratio   {r:.5f}")
if not (0.995 <= r <= 1.005):
    raise SystemExit(
        f"ERROR: the box filter moved the total by {100 * (r - 1):+.2f}%. It is supposed to "
        f"be mean-preserving, so the control would now differ from real DNase in MAGNITUDE "
        f"as well as shape and could not isolate shape. Not usable.")
print(f"wrote {a.out_bw} (width {W} bp, total preserved)")
