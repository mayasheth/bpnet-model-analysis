#!/usr/bin/env python
"""Thin a 5' insertion-count bigwig to a target library depth by binomial subsampling.

WHY BINOMIAL AND NOT A SCALE FACTOR. A 5' track's value at a base IS a count of reads whose
5' end fell there. Drawing Binomial(n, p) per base is therefore exactly equivalent to
subsampling reads at rate p, and it reproduces what a shallower library would actually look
like: bases with one read mostly drop to zero, the zero fraction rises, and the sparsity the
model sees changes. Multiplying every value by p keeps the exact same non-zero structure and
would leave the zero fraction untouched, so it would test nothing.

WHY MATCH THE GENOME-WIDE TOTAL RATHER THAN THE PER-WINDOW MEAN. Per-window means are
computed over each cell type's own DNase-derived elements, and those element calls are
themselves depth-dependent, so matching on them would be partly circular. The genome-wide
total is the library property, and it is the thing that would have been controlled at the
sequencing stage. Per-window means are reported after thinning so the residual gap is
visible rather than assumed away.

WHAT THIS IS FOR. F-016: no H3K27ac model transfers across cell types once the floor uses
the same assay as the model, in any of six directions. The named candidate mechanism is that
per-window DNase depth differs sharply (K562 894.6, THP-1 281.4, GM12878 126.3), so a model
trained at one depth meets a different input distribution at another, which is the
fragment-channel lesson of F-005. The depth RATIOS do not order with transfer quality
(K562->GM12878 is a 7.1x mismatch at -0.045 while K562->THP-1 is 3.2x at -0.388), so this is
a candidate rather than a likely cause, and it is the only one that has not been tested.

Usage:
  0.41.thin_5prime_track.py --in-bw IN.bw --out-bw OUT.bw --target-total 50299042
  0.41.thin_5prime_track.py --in-bw IN.bw --out-bw OUT.bw -p 0.1715
"""
import argparse
import numpy as np
import pyBigWig

MAIN = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

ap = argparse.ArgumentParser()
ap.add_argument("--in-bw", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--target-total", type=float, default=None,
                help="thin so the main-chromosome total lands here")
ap.add_argument("-p", "--prob", type=float, default=None,
                help="explicit retention probability; overrides --target-total")
ap.add_argument("--chunk", type=int, default=10_000_000)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--tol", type=float, default=0.02,
                help="fail if the achieved total misses the target by more than this")
a = ap.parse_args()

src = pyBigWig.open(a.in_bw)
chroms = {c: n for c, n in src.chroms().items() if c in MAIN and n > 0}
if not chroms:
    raise SystemExit(f"no main chromosomes in {a.in_bw}")

if a.prob is None:
    if a.target_total is None:
        raise SystemExit("give --target-total or -p")
    total = 0.0
    for c, n in chroms.items():
        total += float(src.stats(c, 0, n, type="sum", exact=True)[0] or 0.0)
    p = a.target_total / total
    print(f"input total {total:,.0f} -> target {a.target_total:,.0f}: p = {p:.4f}")
    if p > 1.0:
        raise SystemExit(f"p = {p:.4f} > 1; cannot thin UP to a deeper library")
else:
    p = a.prob
    print(f"using p = {p:.4f}")

rng = np.random.default_rng(a.seed)
out = pyBigWig.open(a.out_bw, "w")
out.addHeader([(c, chroms[c]) for c in MAIN if c in chroms])
kept = 0.0
seen = 0.0
for c in MAIN:
    if c not in chroms:
        continue
    n = chroms[c]
    ck = 0.0
    for s0 in range(0, n, a.chunk):
        e0 = min(s0 + a.chunk, n)
        v = src.values(c, s0, e0, numpy=True)
        if v is None:
            continue
        v = np.nan_to_num(v, nan=0.0)
        seen += float(v.sum())
        nz = v > 0
        if not nz.any():
            continue
        # Counts are integers in a 5' track; round defensively in case of float storage.
        counts = np.rint(v[nz]).astype(np.int64)
        thinned = rng.binomial(counts, p).astype(np.float64)
        w = np.zeros_like(v)
        w[nz] = thinned
        keep = w > 0
        if keep.any():
            idx = np.flatnonzero(keep)
            vv = w[idx]
            brk = np.flatnonzero((np.diff(idx) != 1) | (np.diff(vv) != 0)) + 1
            rb = np.concatenate([[0], brk, [len(idx)]])
            st = (idx[rb[:-1]] + s0).astype(np.int64)
            en = (idx[rb[1:] - 1] + 1 + s0).astype(np.int64)
            va = vv[rb[:-1]]
            out.addEntries([c] * len(st), st.tolist(), ends=en.tolist(),
                           values=va.tolist())
            ck += float(((en - st) * va).sum())
        del v, w, nz, keep
    kept += ck
    print(f"  {c}: {ck:,.0f}", flush=True)
src.close()
out.close()

print(f"\ninput {seen:,.0f} -> output {kept:,.0f}  (ratio {kept / max(seen, 1):.4f}, "
      f"requested {p:.4f})")
if a.target_total is not None:
    miss = abs(kept - a.target_total) / a.target_total
    print(f"target {a.target_total:,.0f}, achieved {kept:,.0f}, off by {100 * miss:.2f}%")
    if miss > a.tol:
        raise SystemExit(f"ERROR: missed the target depth by {100 * miss:.2f}%, over the "
                         f"{100 * a.tol:.0f}% tolerance; the thinned track is not "
                         f"depth-matched and must not be trained on")
print(f"wrote {a.out_bw}")
