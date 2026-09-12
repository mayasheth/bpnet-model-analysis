#!/usr/bin/env python
"""Concatenate per-chromosome bigwigs into one, streaming.

4.20 paints one chromosome per job so the genome-wide track can be built in parallel.
Merging is a concatenation, not a real merge, because chromosomes are disjoint.

STREAMS IN CHUNKS. The first version called pyBigWig.intervals(chrom), which materialises
every interval as a Python tuple. The painted track is DENSE -- a softmax profile never
emits zero, so essentially every base carries a value, unlike a real 5-prime read track
which is mostly empty -- so chr1 alone is ~250M tuples and the process was SIGKILLed on the
login node. Values are read and written in windows instead, which bounds memory regardless
of chromosome size.

--MIN-VALUE EXISTS BECAUSE OF THAT SAME DENSITY. Real 5-prime tracks are integer read
counts: a base either has reads or has none. A painted base carrying 0.02 of a read is
below the resolution the real track can represent, and in aggregate that floor is what
makes the converted track 22 GB against real DNase's 828 MB and inflates the lowest
observed-signal quintile 2.7x. Zeroing sub-threshold bases puts the two tracks in the same
space. It defaults to 0 (off) because it changes the data, and the cost to the profile
correlation should be measured on one chromosome before it is switched on.

REFUSES TO WRITE A PARTIAL TRACK unless --allow-missing. A track silently missing a
chromosome trains a model that reads zero accessibility there as closed rather than unknown.

Usage: 4.22.merge_chrom_bigwigs.py --in-glob 'PREFIX_*.bw' --out-bw OUT.bw \
           --chrom-sizes SIZES [--min-value 0.5] [--allow-missing]
"""
import argparse, glob
import numpy as np
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--in-glob", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("--min-value", type=float, default=0.0,
                help="Zero any base below this. 0.5 means 'less than one read', which is "
                     "what a real count track would record as empty.")
ap.add_argument("--chunk", type=int, default=10_000_000)
ap.add_argument("--allow-missing", action="store_true")
a = ap.parse_args()

sizes = {}
for line in open(a.chrom_sizes):
    p = line.split()
    if len(p) >= 2:
        sizes[p[0]] = int(p[1])

files = sorted(glob.glob(a.in_glob))
if not files:
    raise SystemExit(f"no files matched {a.in_glob}")

by_chrom = {}
for f in files:
    b = pyBigWig.open(f)
    ch = [c for c, n in b.chroms().items() if n > 0]
    b.close()
    if len(ch) != 1:
        raise SystemExit(f"{f} carries {len(ch)} chromosomes ({ch[:4]}); expected exactly 1")
    if ch[0] in by_chrom:
        raise SystemExit(f"two inputs both carry {ch[0]}")
    by_chrom[ch[0]] = f

missing = [c for c in sizes if c not in by_chrom]
if missing and not a.allow_missing:
    raise SystemExit(
        f"no input for {sorted(missing)}\nWriting anyway would give the model zero "
        f"accessibility there, which it reads as closed chromatin rather than as missing "
        f"data. Pass --allow-missing only if that is genuinely what you want.")
if missing:
    print(f"WARNING: proceeding without {sorted(missing)}")

out_chroms = [c for c in sorted(by_chrom) if c in sizes]
bw = pyBigWig.open(a.out_bw, "w")
bw.addHeader([(c, sizes[c]) for c in out_chroms])

total = 0.0
n_written = 0
for c in out_chroms:
    src = pyBigWig.open(by_chrom[c])
    n = sizes[c]
    csum = 0.0
    niv = 0
    for s0 in range(0, n, a.chunk):
        e0 = min(s0 + a.chunk, n)
        v = src.values(c, s0, e0, numpy=True)
        if v is None:
            continue
        v = np.nan_to_num(v, nan=0.0).astype(np.float32)
        if a.min_value > 0:
            v[v < a.min_value] = 0.0
        nz = v > 0
        if not nz.any():
            del v, nz
            continue
        idx = np.flatnonzero(nz)
        vv = v[idx]
        # Run-length encode: split where positions are not contiguous OR the value changes.
        cut = np.flatnonzero((np.diff(idx) != 1) | (np.diff(vv) != 0)) + 1
        rb = np.concatenate([[0], cut, [len(idx)]])
        starts = (idx[rb[:-1]] + s0).astype(np.int64)
        ends = (idx[rb[1:] - 1] + 1 + s0).astype(np.int64)
        vals = vv[rb[:-1]].astype(np.float64)
        bw.addEntries([c] * len(starts), starts.tolist(), ends=ends.tolist(),
                      values=vals.tolist())
        csum += float(((ends - starts) * vals).sum())
        niv += len(starts)
        del v, nz, idx, vv, starts, ends, vals
    src.close()
    total += csum
    n_written += niv
    print(f"  {c}: {niv:,} intervals  sum {csum:,.0f}", flush=True)

bw.close()
print(f"\n{len(out_chroms)} chromosomes, {n_written:,} intervals, "
      f"genome-wide sum {total:,.0f}")
print(f"min-value threshold: {a.min_value}")
print(f"wrote {a.out_bw}")
