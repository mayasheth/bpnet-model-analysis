#!/usr/bin/env python
"""
Compute inter-replicate p300 ChIP-seq correlations, stratified by peak overlap.

Reads chromatin annotations TSV (must have EP300_R1.RPM, EP300_R2.RPM, and an
overlap indicator column), computes Pearson/Spearman on log1p(RPM) values for
all elements, peak-overlapping elements, and non-overlapping elements.

Usage:
  python scripts/compute_replicate_correlations.py \
    --chromatin-annot /path/to/chromatin_annotations.tsv \
    --rep1-col EP300_R1.RPM \
    --rep2-col EP300_R2.RPM \
    --overlap-col EP300_peak_overlap \
    --output results/replicate_correlations.tsv
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chromatin-annot", required=True,
                   help="TSV with replicate RPM columns and peak overlap indicator")
    p.add_argument("--rep1-col", default="EP300_R1.RPM")
    p.add_argument("--rep2-col", default="EP300_R2.RPM")
    p.add_argument("--overlap-col", default="EP300_peak_overlap")
    p.add_argument("--output", required=True, help="Output TSV path")
    return p.parse_args()


def corr_row(x_log, y_log, label):
    ok = np.isfinite(x_log) & np.isfinite(y_log)
    x, y = x_log[ok], y_log[ok]
    r, _   = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    return {"subset": label, "n": int(ok.sum()),
            "pearson_r": round(r, 4), "spearman_rho": round(rho, 4)}


def main():
    args = parse_args()

    df = pd.read_table(args.chromatin_annot,
                       usecols=lambda c: c in
                       ["chrom", "chr", "start", "end",
                        args.rep1_col, args.rep2_col, args.overlap_col])

    r1 = np.log1p(df[args.rep1_col].values.astype(float))
    r2 = np.log1p(df[args.rep2_col].values.astype(float))
    ov = df[args.overlap_col].values

    rows = [
        corr_row(r1, r2, "all elements"),
        corr_row(r1[ov == 1], r2[ov == 1], "p300+ elements"),
        corr_row(r1[ov == 0], r2[ov == 0], "p300- elements"),
    ]

    result = pd.DataFrame(rows)
    print(result.to_string(index=False))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
