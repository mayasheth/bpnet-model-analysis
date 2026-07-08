#!/usr/bin/env python
import argparse
import os
import sys
import traceback
from bpnet.cli.shap_scores import shap_scores_main as _main

def parse_args():
    p = argparse.ArgumentParser(description="Compute SHAP for a single chromosome")
    p.add_argument("--reference-genome", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--bed-file", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--output-directory", required=True)
    p.add_argument("--input-seq-len", type=int, required=True)
    p.add_argument("--control-len", type=int, required=True)
    p.add_argument("--task-id", type=int, required=True)
    p.add_argument("--input-data", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--counts-only", action="store_true")
    p.add_argument("--generate-shap-bigWigs", action="store_true")
    return p.parse_args()

def main():
    args = parse_args()

    sys.argv = [
        "bpnet-shap",
        "--reference-genome", args.reference_genome,
        "--model", args.model,
        "--bed-file", args.bed_file,
        "--chroms", args.chrom,
        "--output-directory", args.output_directory,
        "--input-seq-len", str(args.input_seq_len),
        "--control-len", str(args.control_len),
        "--task-id", str(args.task_id),
        "--input-data", args.input_data,
        "--chrom-sizes", args.chrom_sizes,
    ]
    if args.counts_only:
        sys.argv.append("--counts-only")
    if args.generate_shap_bigWigs:
        sys.argv.append("--generate-shap-bigWigs")

    try:
        _main()
        with open(f"{args.output_directory}/DONE.txt", "w") as f:
            f.write("SHAP computation completed successfully.\n")

    except Exception as e:
        print(f"[ERROR] SHAP computation failed: {e}")
        traceback.print_exc()

        # Clean up potentially incomplete outputs
        for fname in ["counts_scores.h5", "profile_scores.h5", "peaks_valid_scores.bed"]:
            fpath = os.path.join(args.output_directory, fname)
            if os.path.exists(fpath):
                os.remove(fpath)
                print(f"[CLEANUP] Removed incomplete: {fpath}")
        sys.exit(1)

if __name__ == "__main__":
    main()
