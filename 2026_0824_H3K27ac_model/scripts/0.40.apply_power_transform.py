#!/usr/bin/env python
"""Apply the fitted power transform to the painted parts and write one track.

`v -> scale * v^gamma`, gamma from 0.39 (fitted in GM12878, so no target-cell DNase enters
the choice). Strictly monotone, so the ordering is untouched; only the spacing between strong
and weak changes, which is the defect F-013 pointed at.

GAMMA FITTED ON ONE CELL TYPE WILL NOT LAND PERFECTLY ON ANOTHER, AND THAT IS THE HONEST
DEPLOYMENT SITUATION. Dynamic range goes roughly as DR^gamma: GM12878's painted 10.2x became
17.3x against its 16.0x target, and K562's painted 18x should reach about 32x against its 39x
target, i.e. ~82%. Tuning gamma on K562's own DNase would close that gap and would also be
the leak this whole design avoids.

RESCALE THEN THRESHOLD, in that order. The power changes the value scale, so the track is
rescaled back to the original painted genome-wide total before the sub-one-read floor from
D-22 is applied; otherwise "less than one read" would mean something different after every
gamma. The floor is a property of what an integer count track can represent, not a quantity
borrowed from the target cell type, which is what keeps it legal.

Usage:
  0.40.apply_power_transform.py --in-glob 'parts/*.bw' --out-bw OUT.bw \\
      --chrom-sizes SIZES --gamma 1.2 [--min-value 0.5]
"""
import argparse, glob
import numpy as np
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--in-glob", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("--gamma", type=float, required=True)
ap.add_argument("--min-value", type=float, default=0.5)
ap.add_argument("--chunk", type=int, default=10_000_000)
a = ap.parse_args()

sizes = {}
for line in open(a.chrom_sizes):
    p = line.split()
    if len(p) >= 2:
        sizes[p[0]] = int(p[1])

parts = {}
for f in sorted(glob.glob(a.in_glob)):
    b = pyBigWig.open(f)
    ch = [c for c, n in b.chroms().items() if n > 0]
    b.close()
    if len(ch) == 1 and ch[0] in sizes:
        parts[ch[0]] = f
chroms = sorted(parts)
if not chroms:
    raise SystemExit(f"no parts matched {a.in_glob}")
print(f"{len(chroms)} parts, gamma = {a.gamma}")

# Pass 1: totals before and after the power, so the rescale factor is exact rather than
# estimated from a sample.
print("pass 1: totals")
tot_raw = tot_pow = 0.0
for c in chroms:
    src = pyBigWig.open(parts[c])
    for s0 in range(0, sizes[c], a.chunk):
        v = src.values(c, s0, min(s0 + a.chunk, sizes[c]), numpy=True)
        if v is None:
            continue
        v = np.nan_to_num(v, nan=0.0).astype(np.float64)
        tot_raw += v.sum()
        tot_pow += (v ** a.gamma).sum()
        del v
    src.close()
scale = tot_raw / tot_pow if tot_pow else 1.0
print(f"  raw total {tot_raw:,.0f}; after power {tot_pow:,.3e}; rescale x{scale:.6g}")

print("pass 2: write")
bw = pyBigWig.open(a.out_bw, "w")
bw.addHeader([(c, sizes[c]) for c in chroms])
total = 0.0
n_iv = 0
for c in chroms:
    src = pyBigWig.open(parts[c])
    csum = 0.0
    niv = 0
    for s0 in range(0, sizes[c], a.chunk):
        e0 = min(s0 + a.chunk, sizes[c])
        v = src.values(c, s0, e0, numpy=True)
        if v is None:
            continue
        v = np.nan_to_num(v, nan=0.0).astype(np.float64)
        out = (v ** a.gamma) * scale
        out[out < a.min_value] = 0.0
        out = out.astype(np.float32)
        nz = out > 0
        if nz.any():
            idx = np.flatnonzero(nz)
            vv = out[idx]
            brk = np.flatnonzero((np.diff(idx) != 1) | (np.diff(vv) != 0)) + 1
            rb = np.concatenate([[0], brk, [len(idx)]])
            st = (idx[rb[:-1]] + s0).astype(np.int64)
            en = (idx[rb[1:] - 1] + 1 + s0).astype(np.int64)
            va = vv[rb[:-1]].astype(np.float64)
            bw.addEntries([c] * len(st), st.tolist(), ends=en.tolist(), values=va.tolist())
            csum += float(((en - st) * va).sum())
            niv += len(st)
        del v, out, nz
    src.close()
    total += csum
    n_iv += niv
    print(f"  {c}: {niv:,} intervals  sum {csum:,.0f}", flush=True)
bw.close()
print(f"\n{len(chroms)} chromosomes, {n_iv:,} intervals, genome-wide sum {total:,.0f}")
print(f"wrote {a.out_bw}")
