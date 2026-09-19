#!/usr/bin/env python3
"""log1p mean and std of an accessibility track over the regions a model will predict on.

WHY THIS IS NOT A DETAIL. `4.1` standardizes the accessibility input with statistics saved at
TRAINING time. In-cell that is right. For a TRANSFERRED model it is wrong twice over: the mean
was computed on the source cell type's library, and on the source's training windows rather
than on the ABC candidate regions. Measured on 2026-09-17 for the p300 arms, against K562 ATAC
on the ABC regions (log1p mean 3.7232, std 1.1428):

    K562 p300 model,    stored 4.155 / 1.313   ->  centre off by +0.38 sd
    GM12878 p300 model, stored 4.552 / 1.360   ->  centre off by +0.73 sd

So the transferred arm was evaluated on inputs shifted about twice as far from centre as the
in-cell arm's, which is a difference between the arms that has nothing to do with transfer.
Feed the numbers this prints to `4.1 --acc-mean/--acc-std` to remove it.

Sampling is on by default because the statistic is a mean over ~10^8 bases and 20,000 windows
already pins it to three decimals; use --all for the exact value.

Usage:
  4.26.acc_norm_stats.py --regions candidateRegions.bed --accessibility-bw atac.bw
"""
import argparse
import numpy as np
import pandas as pd
import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--regions", required=True)
ap.add_argument("--accessibility-bw", required=True)
ap.add_argument("--in-window", type=int, default=2114,
                help="must match the model's in_window, since the statistic is taken over "
                     "exactly the windows the model will see")
ap.add_argument("--sample", type=int, default=20000)
ap.add_argument("--all", action="store_true")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

half = a.in_window // 2
reg = pd.read_csv(a.regions, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["chr", "start", "end"])
b = pyBigWig.open(a.accessibility_bw)
chroms = b.chroms()
reg = reg[reg["chr"].isin(chroms)].reset_index(drop=True)

idx = np.arange(len(reg))
if not a.all and len(reg) > a.sample:
    idx = np.random.default_rng(a.seed).choice(len(reg), size=a.sample, replace=False)

# Streaming moments: the full matrix would be 153,545 x 2,114 floats, and only the mean and
# variance are wanted.
n = 0
tot = 0.0
totsq = 0.0
skipped = 0
for i in idx:
    row = reg.iloc[i]
    c = (int(row["start"]) + int(row["end"])) // 2
    s, e = c - half, c + half
    if s < 0 or e > chroms[row["chr"]]:
        skipped += 1
        continue
    v = b.values(row["chr"], s, e, numpy=True)
    if v is None:
        skipped += 1
        continue
    l = np.log1p(np.clip(np.nan_to_num(v, nan=0.0), 0, None))
    n += l.size
    tot += float(l.sum())
    totsq += float((l ** 2).sum())
b.close()

mean = tot / n
std = float(np.sqrt(max(totsq / n - mean ** 2, 0.0)))
print(f"regions      {a.regions}")
print(f"track        {a.accessibility_bw}")
print(f"windows      {len(idx) - skipped:,} of {len(reg):,} "
      f"({'exhaustive' if a.all else 'sampled'}, {skipped:,} out of bounds)")
print(f"in_window    {a.in_window}")
print()
print(f"--acc-mean {mean:.4f} --acc-std {std:.4f}")
