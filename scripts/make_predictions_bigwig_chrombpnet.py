#!/usr/bin/env python
"""
Convert ChromBPNet mean profile predictions (h5) to a BigWig.

Reads pred_prof (N, output_len) from the mean h5 and writes per-base predicted
counts to a BigWig. Regions are centered on coords_center; profile covers
[center - output_len//2, center + output_len//2]. Overlapping regions are split
at their midpoints.

Input h5 (output of 2.5.mean_predictions_chrombpnet.py):
  coords_chrom (N,), coords_center (N,), pred_prof (N, output_len)
"""

import argparse
import numpy as np
import h5py
import pyBigWig


def parse_args():
    parser = argparse.ArgumentParser(description="Profile h5 to BigWig for ChromBPNet")
    parser.add_argument("--mean-h5", required=True,
                        help="Mean predictions h5 (output of 2.5.mean_predictions_chrombpnet.py)")
    parser.add_argument("--chrom-sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--output-bw", required=True, help="Output BigWig path")
    return parser.parse_args()


def main():
    args = parse_args()

    chrom_sizes = []
    with open(args.chrom_sizes) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom_sizes.append((parts[0], int(parts[1])))
    chrom_order = {c: i for i, (c, _) in enumerate(chrom_sizes)}

    with h5py.File(args.mean_h5, "r") as hf:
        chroms = hf["coords_chrom"][()].astype("U")
        centers = hf["coords_center"][()]
        profiles = hf["pred_prof"][()]

    output_len = profiles.shape[1]
    half_out = output_len // 2
    print(f"Regions: {len(chroms)}, profile output_len: {output_len}")

    # Filter to chromosomes present in chrom_sizes, then sort
    mask = np.array([c in chrom_order for c in chroms])
    chroms = chroms[mask]
    centers = centers[mask]
    profiles = profiles[mask]

    sort_idx = np.lexsort((centers, [chrom_order[c] for c in chroms]))
    chroms = chroms[sort_idx]
    centers = centers[sort_idx]
    profiles = profiles[sort_idx]

    bw = pyBigWig.open(args.output_bw, "w")
    bw.addHeader(chrom_sizes)

    cur_chr = ""
    cur_end = 0

    for i in range(len(chroms)):
        chrom = chroms[i]
        center = int(centers[i])
        reg_start = center - half_out
        reg_end = center + half_out

        if chrom != cur_chr:
            cur_chr = chrom
            cur_end = 0

        # find midpoint with next region to resolve overlaps
        if i + 1 < len(chroms) and chroms[i + 1] == chrom:
            next_center = int(centers[i + 1])
            write_end = min(reg_end, (center + next_center) // 2)
        else:
            write_end = reg_end

        write_start = max(reg_start, cur_end)

        if write_start >= write_end:
            cur_end = write_end
            continue

        # slice profile values to the region we're writing
        prof_start_idx = write_start - reg_start
        prof_end_idx = write_end - reg_start
        vals = profiles[i, prof_start_idx:prof_end_idx].tolist()

        n = write_end - write_start
        bw.addEntries(
            [chrom] * n,
            list(range(write_start, write_end)),
            ends=list(range(write_start + 1, write_end + 1)),
            values=vals,
        )
        cur_end = write_end

    bw.close()
    print(f"Saved: {args.output_bw}")


if __name__ == "__main__":
    main()
