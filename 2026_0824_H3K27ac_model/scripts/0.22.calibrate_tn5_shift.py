#!/usr/bin/env python3
"""Score candidate Tn5 shifts by correlating a BAM-derived 5'-end track against atac_5p.bw.

Input: a two-column file of `P|M <0-based 5' position>` for one chromosome.
Output: a table of Pearson r for every (plus_delta, minus_delta) in -8..+8, best first.
"""
import sys
import numpy as np
import pyBigWig

ends_path, chrom = sys.argv[1], sys.argv[2]
# Reference bigwig is an argument so the same calibration runs for any cell type; the K562
# path is kept as the default for backward compatibility with the original invocation.
BW = sys.argv[3] if len(sys.argv) > 3 else (
    "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
    "2026_0529_multimodal_p300_model/data/atac_5p.bw")

bw = pyBigWig.open(BW)
n = bw.chroms()[chrom]
ref = np.nan_to_num(bw.values(chrom, 0, n, numpy=True), nan=0.0).astype(np.float64)
bw.close()

plus, minus = [], []
with open(ends_path) as f:
    for line in f:
        s, p = line.split()
        (plus if s == "P" else minus).append(int(p))
plus = np.array(plus, dtype=np.int64)
minus = np.array(minus, dtype=np.int64)
print(f"{chrom}: {len(plus)} plus-strand ends, {len(minus)} minus-strand ends, "
      f"bigwig sum {ref.sum():.3g}", file=sys.stderr)

# Bin to counts once; a shift is then a cheap np.roll.
cp = np.bincount(plus, minlength=n).astype(np.float64)
cm = np.bincount(minus, minlength=n).astype(np.float64)

# Correlate on a dense interior slice to avoid roll wraparound at the edges.
lo, hi = 1_000_000, n - 1_000_000
r_ref = ref[lo:hi]

rows = []
for dp in range(-8, 9):
    sp = np.roll(cp, dp)[lo:hi]
    for dm in range(-8, 9):
        sm = np.roll(cm, dm)[lo:hi]
        rows.append((np.corrcoef(sp + sm, r_ref)[0, 1], dp, dm))

rows.sort(reverse=True)
print(f"{'rank':>4} {'plus_delta':>11} {'minus_delta':>12} {'pearson_r':>10}")
for i, (r, dp, dm) in enumerate(rows[:10]):
    print(f"{i+1:>4} {dp:>11} {dm:>12} {r:>10.4f}")
zero = [r for r, dp, dm in rows if dp == 0 and dm == 0][0]
print(f"\nr at (0,0) = {zero:.4f}")
print(f"best       = {rows[0][0]:.4f} at (plus={rows[0][1]}, minus={rows[0][2]})")
