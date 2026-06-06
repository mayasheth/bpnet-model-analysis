#!/usr/bin/env python
"""
Compute mean SHAP scores across ChromBPNet folds.
Adapted from bpnet/utils/mean_shap_plus_peaks.py.

Input:  per-fold shap_counts_merged.h5 files (one per fold)
Output: counts_mean_shap_scores.h5 and counts_peaks_valid_scores.bed
"""

import argparse
import os
import gc
import numpy as np
import pandas as pd
import h5py


def parse_args():
    parser = argparse.ArgumentParser(description="Mean SHAP scores across ChromBPNet folds")
    parser.add_argument("--shap-h5s", required=True,
                        help="Comma-separated list of per-fold shap_counts_merged.h5 paths")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()

    shap_files = [f.strip() for f in args.shap_h5s.split(",")]
    print(f"Averaging {len(shap_files)} fold SHAP files")

    hyp_scores_lst, input_seqs_lst = [], []
    chrom_lst, start_lst, end_lst = [], [], []

    for path in shap_files:
        print(f"  Loading: {path}")
        with h5py.File(path, "r") as f:
            hyp_scores_lst.append(f["hyp_scores"][()])
            input_seqs_lst.append(f["input_seqs"][()])
            chrom_lst.append(f["coords_chrom"][()])
            start_lst.append(f["coords_start"][()])
            end_lst.append(f["coords_end"][()])

    # Verify all folds have the same regions
    aligned = (
        all(np.array_equal(x, chrom_lst[0]) for x in chrom_lst) and
        all(np.array_equal(x, start_lst[0]) for x in start_lst) and
        all(np.array_equal(x, end_lst[0]) for x in end_lst)
    )
    if not aligned:
        raise ValueError("Region coordinates do not match across folds")
    print("All folds aligned. Computing mean...")

    hyp_scores_mean = np.nanmean(np.array(hyp_scores_lst), axis=0)
    print(f"Mean hyp_scores shape: {hyp_scores_mean.shape}")

    os.makedirs(args.output_dir, exist_ok=True)

    h5_path = os.path.join(args.output_dir, "counts_mean_shap_scores.h5")
    with h5py.File(h5_path, "w") as f:
        f.create_dataset("coords_chrom", data=chrom_lst[0], compression="gzip")
        f.create_dataset("coords_start", data=start_lst[0], compression="gzip")
        f.create_dataset("coords_end", data=end_lst[0], compression="gzip")
        f.create_dataset("hyp_scores", data=hyp_scores_mean, compression="gzip")
        f.create_dataset("input_seqs", data=input_seqs_lst[0], compression="gzip")
    print(f"Saved mean SHAP h5: {h5_path}")

    # Write BED file of regions used
    chroms = chrom_lst[0].astype("U")
    starts = start_lst[0]
    ends = end_lst[0]
    n = len(chroms)
    bed_df = pd.DataFrame({
        "chrom": chroms, "start": starts, "end": ends,
        "name": ["."] * n, "score": [0] * n, "strand": ["."] * n,
        "signalValue": [0.0] * n, "p": [0.0] * n, "q": [0.0] * n,
        "summit": [(e - s) // 2 for s, e in zip(starts, ends)]
    })
    bed_path = os.path.join(args.output_dir, "counts_peaks_valid_scores.bed")
    bed_df.to_csv(bed_path, sep="\t", header=False, index=False)
    print(f"Saved regions BED: {bed_path}")

    del hyp_scores_lst, input_seqs_lst, hyp_scores_mean
    gc.collect()


if __name__ == "__main__":
    main()
