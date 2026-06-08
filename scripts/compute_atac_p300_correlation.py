#!/usr/bin/env python
"""
Compute Pearson/Spearman correlation between observed ATAC counts and observed p300 counts
across K562 candidate elements, stratified by p300 peak status.

Inputs:
  - Candidate elements narrowPeak (for coordinates)
  - ATAC BigWig (observed accessibility signal)
  - CV predictions TSV (for true p300 logcounts + peak overlap annotation)

Output:
  - TSV with per-element ATAC logcounts + p300 logcounts
  - Printed correlation summary
"""

import argparse
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import pearsonr, spearmanr


NARROWPEAK_COLS = ["chrom", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--elements", required=True,
                   help="Candidate elements narrowPeak (same used for predictions)")
    p.add_argument("--atac-bw", required=True,
                   help="ATAC BigWig (observed accessibility)")
    p.add_argument("--cv-predictions", required=True,
                   help="cv_predictions.tsv.gz from compute_prediction_performance.py")
    p.add_argument("--in-window", type=int, default=2114,
                   help="Window size to extract ATAC signal over (centered on element)")
    p.add_argument("--output", required=True, help="Output TSV path")
    return p.parse_args()


def extract_atac_counts(elements_df, atac_bw_path, in_window):
    """Extract mean ATAC signal over each element window from BigWig."""
    bw = pyBigWig.open(atac_bw_path)
    chrom_sizes = dict(bw.chroms())
    half = in_window // 2
    atac_counts = []

    for _, row in elements_df.iterrows():
        center = (row["start"] + row["end"]) // 2
        w_start = max(0, center - half)
        w_end = min(chrom_sizes.get(row["chrom"], 0), center + half)
        if w_end <= w_start:
            atac_counts.append(np.nan)
            continue
        vals = bw.values(row["chrom"], w_start, w_end, numpy=True)
        vals = np.nan_to_num(vals, nan=0.0)
        atac_counts.append(float(np.sum(vals)))

    bw.close()
    return np.array(atac_counts)


def main():
    args = parse_args()

    # Load elements
    elements = pd.read_csv(args.elements, sep="\t", header=None,
                           names=NARROWPEAK_COLS, usecols=range(10))
    elements["region_name"] = (elements["chrom"] + ":" +
                               (elements["start"] + 75).astype(str) + "-" +
                               (elements["end"] - 75).astype(str))

    # Load CV predictions (true p300 logcounts + peak overlap)
    preds = pd.read_csv(args.cv_predictions, sep="\t")

    # Extract ATAC counts
    print(f"Extracting ATAC signal for {len(elements):,} elements...")
    atac_raw = extract_atac_counts(elements, args.atac_bw, args.in_window)
    elements["atac_logcounts"] = np.log1p(atac_raw)

    # Merge on region_name — prediction TSV uses 1000bp output windows so
    # chrom/start/end don't match the narrowPeak input coordinates directly
    merged = preds.merge(
        elements[["name", "atac_logcounts"]].rename(columns={"name": "region_name"}),
        on="region_name",
        how="inner"
    )
    print(f"Merged {len(merged):,} elements")

    merged.to_csv(args.output, sep="\t", index=False)

    # Compute correlations
    def report(label, df):
        x = df["atac_logcounts"].values
        y = df["true_logcounts"].values
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        pr, _ = pearsonr(x, y)
        sr, _ = spearmanr(x, y)
        print(f"  {label:20s}  n={len(x):7,}  Pearson={pr:.3f}  Spearman={sr:.3f}")

    print("\nATAC counts vs. observed p300 logcounts:")
    report("All elements", merged)
    report("p300+", merged[merged["EP300_peak_overlap"] == 1])
    report("p300-", merged[merged["EP300_peak_overlap"] == 0])


if __name__ == "__main__":
    main()
