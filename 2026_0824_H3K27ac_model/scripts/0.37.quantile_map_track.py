#!/usr/bin/env python
"""Quantile-map the painted converter track onto a real accessibility value distribution.

WHY. F-013: the painted track reproduces the DNase profile at 95% of its inter-replicate
ceiling and still makes H3K27ac predictions resolvably WORSE than raw ATAC (-0.0099
[-0.0175, -0.0023]). The control in the same run showed the advantage DNase carries is
sub-250 bp structure, which is exactly what the painted track is good at, so the defect is
not shape. It is magnitude: a softmax profile times a predicted count inflates the lowest
observed-signal quintile 2.7x (1.54x even after thresholding, D-22) and lands the
genome-wide total at 34% of real. Thresholding fixed the zero fraction and nothing else.

WHAT THIS DOES. Replaces each painted value with the value at the same rank in a real
accessibility track, so the output's marginal distribution equals the reference's while the
painted track's ORDERING is untouched. Ordering is where the shape information lives, so
this discards the part F-013 indicted and keeps the part F-011 established.

THE REFERENCE IS ANOTHER CELL TYPE'S DNase, AND THAT IS THE POINT. Mapping K562's painted
track onto K562's own DNase distribution would leak the very assay the converter exists to
avoid needing, and the resulting number would not describe anything deployable. GM12878 DNase
is used instead: the deployment claim is "I have ATAC here, a converter, and a DNase
distribution from some other cell type", and that is what this reproduces. No K562 DNase
enters the mapping.

GM12878 DNase is the shallower library (53.8M reads against K562's 301.1M), so the output
carries GM12878-like absolute depth. Training computes its own accessibility normalisation, so
a global scale factor is absorbed; it is the SHAPE of the value distribution that matters here
and that is what is being transferred.

BUILT FROM THE UNTHRESHOLDED PARTS. Quantile mapping subsumes thresholding, since matching
the reference's zero fraction sets the floor, so starting from the thresholded merge would
discard information twice.

Usage:
  0.37.quantile_map_track.py --in-glob 'parts/*.bw' --reference-bw REF.bw \\
      --out-bw OUT.bw --chrom-sizes SIZES
"""
import argparse, glob
import numpy as np
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--in-glob", required=True, help="per-chromosome painted parts")
ap.add_argument("--reference-bw", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("--n-sample-windows", type=int, default=600)
ap.add_argument("--sample-window", type=int, default=200_000)
ap.add_argument("--n-quantiles", type=int, default=4096)
ap.add_argument("--chunk", type=int, default=10_000_000)
ap.add_argument("--seed", type=int, default=0)
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
if not parts:
    raise SystemExit(f"no usable parts matched {a.in_glob}")
chroms = sorted(parts)
print(f"{len(chroms)} painted parts")


def sample(path_or_map, chroms_, label):
    """Random windows, weighted by chromosome length, so the sample is representative."""
    rng = np.random.default_rng(a.seed)
    w = np.array([sizes[c] for c in chroms_], dtype=float)
    w /= w.sum()
    out = []
    single = isinstance(path_or_map, str)
    bw = pyBigWig.open(path_or_map) if single else None
    for _ in range(a.n_sample_windows):
        c = chroms_[rng.choice(len(chroms_), p=w)]
        n = sizes[c]
        if n <= a.sample_window:
            s = 0
            e = n
        else:
            s = int(rng.integers(0, n - a.sample_window))
            e = s + a.sample_window
        src = bw if single else pyBigWig.open(path_or_map[c])
        try:
            if c not in src.chroms():
                continue
            v = src.values(c, s, e, numpy=True)
        finally:
            if not single:
                src.close()
        if v is None:
            continue
        out.append(np.nan_to_num(v, nan=0.0).astype(np.float32))
    if single:
        bw.close()
    v = np.concatenate(out)
    print(f"  {label}: sampled {len(v):,} bp, {100 * (v == 0).mean():.1f}% zero, "
          f"max {v.max():.2f}")
    return v


print("sampling value distributions")
ref = sample(a.reference_bw, [c for c in chroms], "reference")
pnt = sample(parts, chroms, "painted")

# Match the zero fraction: the painted floor is set by where the reference stops being empty.
z_ref = float((ref == 0).mean())
cut = float(np.quantile(pnt, z_ref)) if z_ref > 0 else 0.0
print(f"\nreference is {100 * z_ref:.1f}% zero -> painted floor at {cut:.4f}")

probs = np.linspace(0.0, 1.0, a.n_quantiles)
pnt_hi = pnt[pnt > cut]
ref_hi = ref[ref > 0]
if len(pnt_hi) < 1000 or len(ref_hi) < 1000:
    raise SystemExit("too few non-floor values to build a stable mapping")
pnt_q = np.quantile(pnt_hi, probs)
ref_q = np.quantile(ref_hi, probs)
# np.interp needs a strictly increasing x; ties in a discrete count track are common.
keep = np.concatenate([[True], np.diff(pnt_q) > 0])
pnt_q, ref_q = pnt_q[keep], ref_q[keep]
print(f"mapping built on {len(pnt_q):,} knots; painted "
      f"[{pnt_q[0]:.3f}, {pnt_q[-1]:.1f}] -> reference [{ref_q[0]:.3f}, {ref_q[-1]:.1f}]")
del ref, pnt, pnt_hi, ref_hi

bw = pyBigWig.open(a.out_bw, "w")
bw.addHeader([(c, sizes[c]) for c in chroms])
total = 0.0
n_iv = 0
for c in chroms:
    src = pyBigWig.open(parts[c])
    n = sizes[c]
    csum = 0.0
    niv = 0
    for s0 in range(0, n, a.chunk):
        e0 = min(s0 + a.chunk, n)
        v = src.values(c, s0, e0, numpy=True)
        if v is None:
            continue
        v = np.nan_to_num(v, nan=0.0).astype(np.float32)
        out = np.zeros_like(v)
        m = v > cut
        if m.any():
            out[m] = np.interp(v[m], pnt_q, ref_q).astype(np.float32)
        nz = out > 0
        if nz.any():
            idx = np.flatnonzero(nz)
            vv = out[idx]
            brk = np.flatnonzero((np.diff(idx) != 1) | (np.diff(vv) != 0)) + 1
            rb = np.concatenate([[0], brk, [len(idx)]])
            st = (idx[rb[:-1]] + s0).astype(np.int64)
            en = (idx[rb[1:] - 1] + 1 + s0).astype(np.int64)
            va = vv[rb[:-1]].astype(np.float64)
            bw.addEntries([c] * len(st), st.tolist(), ends=en.tolist(),
                          values=va.tolist())
            csum += float(((en - st) * va).sum())
            niv += len(st)
        del v, out, m, nz
    src.close()
    total += csum
    n_iv += niv
    print(f"  {c}: {niv:,} intervals  sum {csum:,.0f}", flush=True)
bw.close()

print(f"\n{len(chroms)} chromosomes, {n_iv:,} intervals, genome-wide sum {total:,.0f}")

# The output distribution must now match the reference's, or the mapping did not do its job.
print("\nverifying the output distribution against the reference")
got = sample(a.out_bw, chroms, "output")
ref2 = sample(a.reference_bw, chroms, "reference")
zg, zr = float((got == 0).mean()), float((ref2 == 0).mean())
print(f"  zero fraction  output {100 * zg:.1f}%  reference {100 * zr:.1f}%")
for q in (0.9, 0.99, 0.999, 0.9999):
    print(f"  q{q:<7} output {np.quantile(got, q):>10.3f}  "
          f"reference {np.quantile(ref2, q):>10.3f}")
if abs(zg - zr) > 0.02:
    raise SystemExit(f"ERROR: zero fraction off by {100 * abs(zg - zr):.1f} points; the "
                     f"floor did not transfer and the magnitude distortion is not fixed")
print(f"\nwrote {a.out_bw}")
