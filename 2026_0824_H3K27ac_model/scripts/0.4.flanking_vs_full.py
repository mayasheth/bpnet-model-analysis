#!/usr/bin/env python
"""
Is a flanking-only H3K27ac count a better target than a full-window count?

H3K27ac sits on the nucleosomes flanking an element, and the meta-profile shows a dip
over the nucleosome-free center. That suggests excluding the center might sharpen the
target. But the dip is shallow (~9% below the shoulders), partly because the coverage
track extends 36 bp reads to a fixed 250 bp fragment, which smears signal across the
NFR. Excluding the center therefore discards a lot of real signal.

The decisive question is not "is the center depleted" but "does excluding it raise the
inter-replicate ceiling", since a target nobody can measure reproducibly cannot be
predicted. Fewer reads means more Poisson noise, so there is a real cost to set against
any gain in specificity.

For each outer half-window and inner exclusion radius, reports:
  - inter-replicate Pearson on log1p counts (all elements and top quintile)
  - correlation between the flanking-only and full-window counts, i.e. how different
    the two targets even are

Usage:
  python 0.4.flanking_vs_full.py --elements <narrowPeak> \
      --rep1-bw data/h3k27ac_rep1.bw --rep2-bw data/h3k27ac_rep2.bw \
      --outdir results
"""

import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import pearsonr

SEED = 0
COLS = ["chr", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", required=True)
    p.add_argument("--rep1-bw", required=True)
    p.add_argument("--rep2-bw", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--outer", type=int, nargs="+", default=[500, 1000])
    p.add_argument("--inner", type=int, nargs="+", default=[0, 125, 250, 375])
    p.add_argument("--n-sample", type=int, default=0, help="0 = all elements")
    return p.parse_args()


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    mid = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, mid)
    return df


def collect(df, bw_path, max_outer):
    """Per-element coverage vector over +/-max_outer; NaN row if out of bounds."""
    bw = pyBigWig.open(bw_path)
    sizes = bw.chroms()
    out = np.full((len(df), 2 * max_outer), np.nan, dtype=np.float32)
    chroms = df["chr"].to_numpy()
    centers = df["center"].to_numpy()
    for i in range(len(df)):
        c, ct = chroms[i], int(centers[i])
        if c not in sizes or ct - max_outer < 0 or ct + max_outer > sizes[c]:
            continue
        v = bw.values(c, ct - max_outer, ct + max_outer, numpy=True)
        if v is None or len(v) != 2 * max_outer:
            continue
        out[i] = np.nan_to_num(v, nan=0.0)
    bw.close()
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    df = load_elements(args.elements)
    if args.n_sample and args.n_sample < len(df):
        df = df.sample(n=args.n_sample, random_state=SEED).reset_index(drop=True)
    max_outer = max(args.outer)
    print(f"Elements: {len(df):,}; reading +/-{max_outer} bp per replicate")

    m1 = collect(df, args.rep1_bw, max_outer)
    m2 = collect(df, args.rep2_bw, max_outer)
    ok = ~(np.isnan(m1).any(axis=1) | np.isnan(m2).any(axis=1))
    m1, m2 = m1[ok], m2[ok]
    print(f"Usable: {len(m1):,}")

    rows = []
    for outer in args.outer:
        lo, hi = max_outer - outer, max_outer + outer
        for inner in args.inner:
            if inner >= outer:
                continue
            a, b = m1[:, lo:hi], m2[:, lo:hi]
            if inner == 0:
                s1, s2 = a.sum(axis=1), b.sum(axis=1)
            else:
                # keep only the two flanks, excluding +/-inner around the center
                keep = np.ones(2 * outer, dtype=bool)
                keep[outer - inner:outer + inner] = False
                s1, s2 = a[:, keep].sum(axis=1), b[:, keep].sum(axis=1)
            x, y = np.log1p(s1), np.log1p(s2)
            mean_sig = (x + y) / 2
            top = mean_sig >= np.quantile(mean_sig, 0.8)
            full = np.log1p(a.sum(axis=1) + b.sum(axis=1))
            comb = np.log1p(s1 + s2)
            rows.append({
                "outer_half_window": outer,
                "inner_exclusion": inner,
                "bp_counted": 2 * (outer - inner),
                "frac_reads_kept": round(float(comb.sum() / full.sum()), 4),
                "ceiling_all": round(float(pearsonr(x, y)[0]), 4),
                "ceiling_top_quintile": round(float(pearsonr(x[top], y[top])[0]), 4),
                "r_vs_full_window": round(float(pearsonr(comb, full)[0]), 4),
            })
    res = pd.DataFrame(rows)
    path = os.path.join(args.outdir, "flanking_vs_full_window.tsv")
    res.to_csv(path, sep="\t", index=False)
    print(f"Wrote {path}\n")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
