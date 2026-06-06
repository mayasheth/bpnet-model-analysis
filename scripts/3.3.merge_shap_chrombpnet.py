#!/usr/bin/env python
"""
Merge per-chromosome SHAP h5 files into a single per-fold h5.
Adapted from bpnet/utils/merge_shap_across_chrom.py.

Input:  <input_dir>/chr*/shap.counts_scores.h5
Output: <output_file>
"""

import argparse
import glob
import os
import numpy as np
import h5py


def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-chromosome SHAP h5 files for one fold")
    parser.add_argument("--input-dir", required=True,
                        help="Per-fold SHAP directory (contains chr*/ subdirs)")
    parser.add_argument("--h5-filename", default="shap.counts_scores.h5",
                        help="Filename within each chromosome directory")
    parser.add_argument("--output-file", required=True,
                        help="Output merged h5 path")
    return parser.parse_args()


def main():
    args = parse_args()

    pattern = os.path.join(args.input_dir, "chr*/", args.h5_filename)
    h5_files = sorted(glob.glob(pattern))
    if not h5_files:
        raise FileNotFoundError(f"No h5 files found matching: {pattern}")
    print(f"Merging {len(h5_files)} chromosome files")

    hyp_scores_lst, input_seqs_lst = [], []
    chrom_lst, start_lst, end_lst = [], [], []

    for path in h5_files:
        with h5py.File(path, "r") as f:
            hyp_scores_lst.append(f["hyp_scores"][()])
            input_seqs_lst.append(f["input_seqs"][()])
            chrom_lst.append(f["coords_chrom"][()])
            start_lst.append(f["coords_start"][()])
            end_lst.append(f["coords_end"][()])
        chrom = os.path.basename(os.path.dirname(path))
        print(f"  {chrom}: {len(chrom_lst[-1])} regions")

    hyp_scores = np.concatenate(hyp_scores_lst)
    input_seqs = np.concatenate(input_seqs_lst)
    chroms = np.concatenate(chrom_lst)
    starts = np.concatenate(start_lst)
    ends = np.concatenate(end_lst)
    print(f"Total: {len(hyp_scores)} regions, shape: {hyp_scores.shape}")

    os.makedirs(os.path.dirname(args.output_file) or ".", exist_ok=True)
    with h5py.File(args.output_file, "w") as f:
        f.create_dataset("coords_chrom", data=chroms, compression="gzip")
        f.create_dataset("coords_start", data=starts, compression="gzip")
        f.create_dataset("coords_end", data=ends, compression="gzip")
        f.create_dataset("hyp_scores", data=hyp_scores, compression="gzip")
        f.create_dataset("input_seqs", data=input_seqs, compression="gzip")

    print(f"Saved: {args.output_file}")


if __name__ == "__main__":
    main()
