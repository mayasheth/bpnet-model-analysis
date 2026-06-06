#!/usr/bin/env python
"""
Compute mean predictions across ChromBPNet folds.

Reads per-fold h5 files (output of predict_chrombpnet_regions.py) and writes:
  - TSV: mean pred_log_counts per region
  - h5:  mean logcounts (N,) and mean profile (N, output_len) across folds
"""

import argparse
import gc
import glob
import os
import numpy as np
import pandas as pd
import h5py


def parse_args():
    parser = argparse.ArgumentParser(description="Mean ChromBPNet predictions across folds")
    parser.add_argument("--pred-dir", required=True,
                        help="Directory containing fold*_predictions.h5 files")
    parser.add_argument("--output", required=True,
                        help="Output TSV path (mean log-counts)")
    parser.add_argument("--output-h5", required=True,
                        help="Output h5 path (mean log-counts + mean profile)")
    return parser.parse_args()


def main():
    args = parse_args()

    h5_files = sorted(glob.glob(os.path.join(args.pred_dir, "fold*_predictions.h5")))
    if not h5_files:
        raise FileNotFoundError(f"No fold*_predictions.h5 files found in {args.pred_dir}")
    print(f"Found {len(h5_files)} fold h5s: {[os.path.basename(f) for f in h5_files]}")

    logcounts_lst = []
    profile_lst = []
    chrom_lst = []
    center_lst = []
    start_lst = []
    end_lst = []

    for path in h5_files:
        print(f"  Loading: {path}")
        with h5py.File(path, "r") as f:
            logcounts_lst.append(f["predictions/logcounts"][()])
            profile_lst.append(f["predictions/profile"][()])
            chrom_lst.append(f["coords/coords_chrom"][()].astype("U"))
            center_lst.append(f["coords/coords_center"][()])
            start_lst.append(f["coords/coords_start"][()])
            end_lst.append(f["coords/coords_end"][()])

    # Verify all folds have the same regions
    for i in range(1, len(h5_files)):
        if not (np.array_equal(chrom_lst[i], chrom_lst[0]) and
                np.array_equal(center_lst[i], center_lst[0])):
            raise ValueError(f"Region mismatch between fold 0 and fold {i}")
    print("All folds aligned. Computing means...")

    mean_logcounts = np.mean(np.array(logcounts_lst), axis=0)
    mean_profile = np.mean(np.array(profile_lst), axis=0).astype(np.float32)
    del logcounts_lst, profile_lst
    gc.collect()

    chroms = chrom_lst[0]
    centers = center_lst[0]
    starts = start_lst[0]
    ends = end_lst[0]

    # Save TSV
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    mean_df = pd.DataFrame({
        "chrom": chroms,
        "start": starts,
        "end": ends,
        "center": centers,
        "pred_log_counts": mean_logcounts,
    })
    mean_df.to_csv(args.output, sep="\t", index=False)
    print(f"Saved mean TSV ({len(mean_df)} regions): {args.output}")

    # Save h5
    os.makedirs(os.path.dirname(args.output_h5) or ".", exist_ok=True)
    with h5py.File(args.output_h5, "w") as hf:
        hf.create_dataset("coords_chrom", data=chroms.astype("S"), compression="gzip")
        hf.create_dataset("coords_center", data=centers, compression="gzip")
        hf.create_dataset("coords_start", data=starts, compression="gzip")
        hf.create_dataset("coords_end", data=ends, compression="gzip")
        hf.create_dataset("pred_logcounts", data=mean_logcounts.astype(np.float32), compression="gzip")
        hf.create_dataset("pred_prof", data=mean_profile, compression="gzip")
    print(f"Saved mean h5: {args.output_h5}")
    print(f"  pred_logcounts shape: {mean_logcounts.shape}")
    print(f"  pred_prof shape:      {mean_profile.shape}")


if __name__ == "__main__":
    main()
