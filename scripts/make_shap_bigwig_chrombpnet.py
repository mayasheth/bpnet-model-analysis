#!/usr/bin/env python
"""
Convert ChromBPNet mean SHAP scores (h5) to a BigWig.

Wraps the existing bpnet importance_hdf5_to_bigwig utility, which expects:
  h5: hyp_scores (N, L, 4), input_seqs (N, L, 4)
  bed: 10-column narrowPeak where col 9 is summit offset from col 1 (start)

These match the outputs of 3.4.mean_shap_chrombpnet.py.
"""

import argparse
import sys

BPNET_DIR = "/oak/stanford/groups/engreitz/Users/sheth/bpnet-refactor"
sys.path.insert(0, BPNET_DIR)
from bpnet.utils.importance_hdf5_to_bigwig import importance_hdf5_to_bigwig


def parse_args():
    parser = argparse.ArgumentParser(description="SHAP h5 to BigWig for ChromBPNet")
    parser.add_argument("--shap-h5", required=True,
                        help="counts_mean_shap_scores.h5 (output of 3.4.mean_shap_chrombpnet.py)")
    parser.add_argument("--regions-bed", required=True,
                        help="counts_peaks_valid_scores.bed (output of 3.4.mean_shap_chrombpnet.py)")
    parser.add_argument("--chrom-sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--output-bw", required=True, help="Output BigWig path")
    parser.add_argument("--output-stats", required=True, help="Output stats text path")
    return parser.parse_args()


def main():
    args = parse_args()
    importance_hdf5_to_bigwig(
        args.shap_h5,
        args.regions_bed,
        args.output_bw,
        args.output_stats,
        args.chrom_sizes,
    )
    print(f"Saved: {args.output_bw}")
    print(f"Stats: {args.output_stats}")


if __name__ == "__main__":
    main()
