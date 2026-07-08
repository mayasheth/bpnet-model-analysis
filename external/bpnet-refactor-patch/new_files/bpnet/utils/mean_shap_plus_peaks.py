import os
import pandas as pd
import h5py
import hdf5plugin
import numpy as np
import argparse
import gc

parser = argparse.ArgumentParser(description="Calculate mean SHAP over the given HDF5s.")
parser.add_argument("--counts_shaps", type=str, default=None,
                    help="Comma-separated list of counts SHAP .h5 files")
parser.add_argument("--profile_shaps", type=str, default=None,
                    help="Comma-separated list of profile SHAP .h5 files")
parser.add_argument("--output_dir", type=str, default='/cromwell_root/',
                    help="Output directory for the mean SHAP scores")
args = parser.parse_args()

def mean_shap(shaps_list, output_prefix):
    hyp_scores_lst = []
    chrom_lst = []
    start_lst = []
    end_lst = []
    input_seqs_lst = []

    for shap_h5 in shaps_list.split(','):
        try:
            with h5py.File(shap_h5, 'r') as f:
                hyp_scores_lst.append(f['hyp_scores'][()])
                chrom_lst.append(f['coords_chrom'][()])
                start_lst.append(f['coords_start'][()])
                end_lst.append(f['coords_end'][()])
                input_seqs_lst.append(f['input_seqs'][()])
        except Exception as e:
            print(f"[WARNING] Failed to load {shap_h5}: {e}")

    [print(len(x) for x in chrom_lst)]

    print(all(np.array_equal(x, chrom_lst[0]) for x in chrom_lst))
    print(all(np.array_equal(x, start_lst[0]) for x in start_lst))
    print(all(np.array_equal(x, end_lst[0]) for x in end_lst))

    if sum([
        all(np.array_equal(x, chrom_lst[0]) for x in chrom_lst),
        all(np.array_equal(x, start_lst[0]) for x in start_lst),
        all(np.array_equal(x, end_lst[0]) for x in end_lst)
    ]) == 3:
        print(f"[INFO] All input SHAP files aligned. Proceeding with mean calculation.")

        hyp_scores_mean = np.nanmean(np.array(hyp_scores_lst), axis=0)
        print(hyp_scores_mean.shape)
        num_examples = hyp_scores_mean.shape[0]

        os.makedirs(args.output_dir, exist_ok=True)

        h5_path = os.path.join(args.output_dir, f"{output_prefix}_mean_shap_scores.h5")
        with h5py.File(h5_path, "w") as f:
            f.create_dataset("coords_chrom", data=chrom_lst[0].astype("S"), **hdf5plugin.Blosc())
            f.create_dataset("coords_start", data=start_lst[0], **hdf5plugin.Blosc())
            f.create_dataset("coords_end", data=end_lst[0], **hdf5plugin.Blosc())
            f.create_dataset("hyp_scores", data=hyp_scores_mean, **hdf5plugin.Blosc())
            f.create_dataset("input_seqs", data=input_seqs_lst[0], **hdf5plugin.Blosc())

        # Save consensus peaks_valid_scores.bed
        chroms = chrom_lst[0].astype("U")
        starts = start_lst[0]
        ends = end_lst[0]
        bed_df = pd.DataFrame({
            "chrom": chroms,
            "start": starts,
            "stop": ends,
            "name": ["."] * num_examples,
            "score": [0] * num_examples,
            "strand": ["."] * num_examples,
            "signalValue": [0.0] * num_examples,
            "p": [0.0] * num_examples,
            "q": [0.0] * num_examples,
            "summit": [(e - s) // 2 for s, e in zip(starts, ends)]
        })
        bed_path = os.path.join(args.output_dir, f"{output_prefix}_peaks_valid_scores.bed")
        bed_df.to_csv(bed_path, sep="\t", header=False, index=False)

        print(f"[INFO] Wrote mean SHAP HDF5 to: {h5_path}")
        print(f"[INFO] Wrote BED to: {bed_path}")

        del input_seqs_lst, hyp_scores_mean, chrom_lst, start_lst, end_lst, hyp_scores_lst
        gc.collect()

if args.counts_shaps:
    mean_shap(args.counts_shaps, "counts")

if args.profile_shaps:
    mean_shap(args.profile_shaps, "profile")
