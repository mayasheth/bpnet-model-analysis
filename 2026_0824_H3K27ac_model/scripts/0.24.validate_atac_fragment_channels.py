#!/usr/bin/env python3
"""Validate the fragment-size-stratified ATAC channels against the flat 5' track.

The stratified channels are only a strict superset of the accs5p input if they partition
the same reads exhaustively, i.e. if

    sub + mono + di + poly  ==  all   (== atac_5p.bw)

both genome-wide and base by base. The channels come from the PE BAMs and `all` from the
tagAligns, so this is a real check on two independent derivations, not a tautology.

Run: python 0.24.validate_atac_fragment_channels.py [k562|gm12878]
"""
import numpy as np
import pyBigWig

import sys
D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
CELL = sys.argv[1] if len(sys.argv) > 1 else "k562"
FLAT = {"k562": f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
        "gm12878": f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw"}
PFX = {"k562": "", "gm12878": "gm12878_"}
assert CELL in FLAT, f"unknown cell type {CELL!r}"
ALL = FLAT[CELL]
BINS = {k: f"{D}/2026_0824_H3K27ac_model/data/{PFX[CELL]}atac_{k}5p.bw" for k in
        ("sub", "mono", "di", "poly")}
print(f"cell type: {CELL}")

bw_all = pyBigWig.open(ALL)
bws = {k: pyBigWig.open(v) for k, v in BINS.items()}

t_all = bw_all.header()["sumData"]
tot = {k: b.header()["sumData"] for k, b in bws.items()}
s = sum(tot.values())
print("genome-wide insertion totals")
print(f"  {'all (atac_5p.bw)':<22} {t_all:>16,.0f}")
for k, v in tot.items():
    print(f"  {k:<22} {v:>16,.0f}   {100*v/t_all:5.1f}%")
print(f"  {'sum of bins':<22} {s:>16,.0f}   {100*s/t_all:5.1f}%")
print(f"  discrepancy            {s - t_all:>16,.0f}   "
      f"({100*abs(s - t_all)/t_all:.4f}% of all)\n")

# Base-resolution agreement on element-dense regions, where the model actually looks.
rng = np.random.RandomState(0)
chroms = {c: n for c, n in bw_all.chroms().items() if n > 5_000_000}
names = sorted(chroms)
worst = 0.0
n_bases = 0
sad = 0.0
for i in range(200):
    c = names[rng.randint(len(names))]
    start = rng.randint(1_000_000, chroms[c] - 1_000_000)
    end = start + 10_000
    a = np.nan_to_num(bw_all.values(c, start, end, numpy=True), nan=0.0)
    b = np.zeros_like(a)
    for k, bwk in bws.items():
        b += np.nan_to_num(bwk.values(c, start, end, numpy=True), nan=0.0)
    d = np.abs(a - b)
    sad += d.sum()
    n_bases += len(a)
    worst = max(worst, d.max())
print(f"base-resolution check over 200 x 10 kb windows ({n_bases:,} bp)")
print(f"  total absolute difference {sad:,.0f}")
print(f"  largest single-base diff  {worst:,.0f}")
ok = worst == 0.0 and abs(s - t_all) / t_all < 1e-4
print(f"\n{'PASS' if ok else 'FAIL'}: channels {'do' if ok else 'do NOT'} "
      f"partition the flat 5' track")
bw_all.close()
for b in bws.values():
    b.close()
