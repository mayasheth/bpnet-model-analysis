#!/usr/bin/env python
"""
Compute DeepSHAP contribution scores (counts head) for a ChromBPNet model.

Processes one chromosome at a time. Outputs per-region hypothetical contribution
scores in HDF5 format compatible with MoDISCo and the rest of this pipeline.

Output h5 keys:
    hyp_scores  (N, L, 4) - hypothetical contribution scores
    input_seqs  (N, L, 4) - one-hot encoded sequences
    coords_chrom (N,)
    coords_start (N,)
    coords_end   (N,)
"""

import argparse
import os
import numpy as np
import pandas as pd
import pyfaidx
import h5py

import tensorflow as tf

# Must disable eager execution before importing shap or loading model
tf.compat.v1.disable_eager_execution()

# Keras 3 removed get_session; patch it back for shap.TFDeepExplainer compatibility
if not hasattr(tf.compat.v1.keras.backend, "get_session"):
    def _get_session():
        sess = tf.compat.v1.get_default_session()
        if sess is None:
            sess = tf.compat.v1.Session()
            sess.__enter__()
        return sess
    tf.compat.v1.keras.backend.get_session = _get_session

import shap
from deeplift.dinuc_shuffle import dinuc_shuffle


NARROWPEAK_SCHEMA = ["chr", "start", "end", "name", "score", "strand",
                     "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    parser = argparse.ArgumentParser(description="Compute SHAP scores (counts head) for ChromBPNet")
    parser.add_argument("-m", "--model-h5", required=True, help="Path to ChromBPNet nobias h5 model")
    parser.add_argument("-r", "--regions", required=True, help="narrowPeak file of regions")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix (dir/prefix)")
    parser.add_argument("--chrom", required=True, help="Chromosome to process (e.g. chr1)")
    parser.add_argument("--num-shuffles", type=int, default=20, help="Dinucleotide shuffles for background")
    return parser.parse_args()


def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    enc = np.zeros((len(seq), 4), dtype=np.float32)
    for i, c in enumerate(seq.upper()):
        if c in mapping:
            enc[i, mapping[c]] = 1.0
    return enc


def shuffle_several_times(s, num_shuffles=20):
    return [np.array([dinuc_shuffle(s[0]) for _ in range(num_shuffles)])]


def combine_mult_and_diffref(mult, orig_inp, bg_data):
    projected = np.zeros_like(bg_data[0]).astype("float")
    for i in range(orig_inp[0].shape[-1]):
        hyp_inp = np.zeros_like(orig_inp[0]).astype("float")
        hyp_inp[:, i] = 1.0
        hyp_diff = hyp_inp[None, :, :] - bg_data[0]
        hyp_contrib = hyp_diff * mult[0]
        projected[:, :, i] = np.sum(hyp_contrib, axis=-1)
    return [np.mean(projected, axis=0)]


def get_sequences(regions_df, genome, input_len):
    half = input_len // 2
    seqs, starts, ends = [], [], []
    valid = np.zeros(len(regions_df), dtype=bool)

    for idx, (_, row) in enumerate(regions_df.iterrows()):
        center = int(row["start"]) + int(row["summit"])
        s, e = center - half, center + half
        chrom = row["chr"]
        if s < 0 or e > len(genome[chrom]):
            continue
        seq = str(genome[chrom][s:e])
        if len(seq) != input_len:
            continue
        seqs.append(one_hot_encode(seq))
        starts.append(s)
        ends.append(e)
        valid[idx] = True

    return np.array(seqs, dtype=np.float32), np.array(starts), np.array(ends), valid


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.output_prefix) or ".", exist_ok=True)

    # Load model (after disable_eager_execution)
    model = tf.keras.models.load_model(args.model_h5, compile=False)
    input_len = model.input_shape[1]
    print(f"Model input length: {input_len}")
    print(f"Model output shapes: {[o.shape for o in model.outputs]}")

    # Load and filter regions for this chromosome
    regions_df = pd.read_csv(args.regions, sep="\t", header=None,
                             names=NARROWPEAK_SCHEMA, usecols=range(10))
    regions_df = regions_df[regions_df["chr"] == args.chrom].reset_index(drop=True)
    print(f"Regions on {args.chrom}: {len(regions_df)}")

    if len(regions_df) == 0:
        print(f"No regions found for {args.chrom}, exiting.")
        return

    # Extract sequences
    genome = pyfaidx.Fasta(args.genome)
    seqs, starts, ends, valid = get_sequences(regions_df, genome, input_len)
    genome.close()

    n_filtered = (~valid).sum()
    if n_filtered > 0:
        print(f"Filtered {n_filtered} regions near chromosome boundaries")
    print(f"Computing SHAP for {len(seqs)} sequences")

    chroms_out = np.array([args.chrom] * len(seqs), dtype="S")

    # Build DeepSHAP explainer for counts head
    # model.outputs[1] is the counts head; reduce_sum over the last axis
    counts_output = tf.keras.layers.Lambda(lambda x: tf.reduce_sum(x, axis=-1))(model.outputs[1])
    explainer = shap.explainers.deep.TFDeepExplainer(
        (model.input, counts_output),
        lambda s: shuffle_several_times(s, num_shuffles=args.num_shuffles),
        combine_mult_and_diffref=combine_mult_and_diffref,
    )

    print("Running SHAP...")
    hyp_scores = explainer.shap_values(seqs, progress_message=100)
    # hyp_scores shape: (N, L, 4)

    # Save output
    out_h5 = f"{args.output_prefix}.counts_scores.h5"
    with h5py.File(out_h5, "w") as hf:
        hf.create_dataset("hyp_scores", data=hyp_scores.astype(np.float16), compression="gzip")
        hf.create_dataset("input_seqs", data=seqs.astype(np.int8), compression="gzip")
        hf.create_dataset("coords_chrom", data=chroms_out, compression="gzip")
        hf.create_dataset("coords_start", data=starts, compression="gzip")
        hf.create_dataset("coords_end", data=ends, compression="gzip")

    print(f"Saved: {out_h5}")

    # Write DONE marker
    with open(f"{args.output_prefix}.DONE.txt", "w") as f:
        f.write("done\n")


if __name__ == "__main__":
    main()
