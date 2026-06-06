#!/usr/bin/env python
"""
Predict ChromBPNet log-counts and per-base profile for a set of genomic regions.
Outputs per-region predicted log-counts as TSV and h5 (with profile).

Profile is saved as per-base predicted counts (not log-softmax):
  counts_per_base = softmax(profile_head) * exp(counts_head)
averaged over forward and reverse-complement passes.
"""

import argparse
import numpy as np
import pandas as pd
import pyfaidx
import h5py
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Predict ChromBPNet log-counts on genomic regions")
    parser.add_argument("-m", "--model-h5", required=True, help="Path to ChromBPNet h5 model")
    parser.add_argument("-r", "--regions", required=True, help="narrowPeak or BED file of regions")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-op", "--output-prefix", required=True, help="Output file prefix")
    parser.add_argument("-bs", "--batch-size", type=int, default=64)
    return parser.parse_args()


def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    enc = np.zeros((len(seq), 4), dtype=np.float32)
    for i, c in enumerate(seq.upper()):
        if c in mapping:
            enc[i, mapping[c]] = 1
    return enc


def rc_onehot(ohe):
    return ohe[::-1, ::-1]


def main():
    args = parse_args()

    import tensorflow as tf

    # compile=False skips loading the loss function (not needed for inference)
    model = tf.keras.models.load_model(args.model_h5, compile=False)

    input_len = model.input_shape[1]
    half = input_len // 2
    output_len = model.output_shape[0][1]
    print(f"Model input length: {input_len}, profile output length: {output_len}")

    # Read regions — narrowPeak has summit offset in col 9 (0-indexed)
    col_names = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit"]
    regions = pd.read_csv(args.regions, sep="\t", header=None, names=col_names[:10])
    regions = regions[["chrom", "start", "end"] + (["summit"] if regions.shape[1] >= 10 else [])].copy()

    if "summit" in regions.columns and regions["summit"].notna().all() and (regions["summit"] != -1).all():
        regions["center"] = regions["start"] + regions["summit"].astype(int)
    else:
        regions["center"] = ((regions["start"] + regions["end"]) // 2).astype(int)

    genome = pyfaidx.Fasta(args.genome)

    # Filter regions too close to chromosome boundaries
    chrom_sizes = {name: len(genome[name]) for name in genome.keys()}
    valid = regions.apply(
        lambda r: r["center"] - half >= 0 and r["center"] + half <= chrom_sizes.get(r["chrom"], 0),
        axis=1,
    )
    n_filtered = (~valid).sum()
    if n_filtered > 0:
        print(f"Filtering {n_filtered} regions near chromosome boundaries")
    regions = regions[valid].reset_index(drop=True)

    log_counts = []
    profiles = []
    for i in tqdm(range(0, len(regions), args.batch_size), desc="Predicting"):
        batch = regions.iloc[i : i + args.batch_size]
        seqs_fwd, seqs_rev = [], []

        for _, row in batch.iterrows():
            s = int(row["center"]) - half
            e = int(row["center"]) + half
            seq = str(genome[row["chrom"]][s:e])
            ohe = one_hot_encode(seq)
            seqs_fwd.append(ohe)
            seqs_rev.append(rc_onehot(ohe))

        X_fwd = np.array(seqs_fwd)
        X_rev = np.array(seqs_rev)

        pred_fwd = model(X_fwd, training=False)
        pred_rev = model(X_rev, training=False)

        # Average forward and reverse complement log-counts
        lc = (pred_fwd[1].numpy()[:, 0] + pred_rev[1].numpy()[:, 0]) / 2
        log_counts.extend(lc.tolist())

        # Per-base predicted counts = softmax(profile) * exp(logcounts)
        # RC profile is reversed to align with forward-strand coordinates
        prof_fwd = np.exp(pred_fwd[0].numpy() + pred_fwd[1].numpy())
        prof_rev = np.exp(pred_rev[0].numpy()[:, ::-1] + pred_rev[1].numpy())
        profiles.append((prof_fwd + prof_rev) / 2)

    regions["pred_log_counts"] = log_counts
    profiles = np.concatenate(profiles, axis=0).astype(np.float32)

    # Save TSV
    out_tsv = args.output_prefix + "_predictions.tsv"
    regions[["chrom", "start", "end", "center", "pred_log_counts"]].to_csv(
        out_tsv, sep="\t", index=False
    )
    print(f"Saved TSV: {out_tsv}")

    # Save h5
    out_h5 = args.output_prefix + "_predictions.h5"
    with h5py.File(out_h5, "w") as hf:
        coords = hf.create_group("coords")
        coords.create_dataset("coords_chrom", data=regions["chrom"].values.astype("S"), compression="gzip")
        coords.create_dataset("coords_center", data=regions["center"].values, compression="gzip")
        coords.create_dataset("coords_start", data=regions["start"].values, compression="gzip")
        coords.create_dataset("coords_end", data=regions["end"].values, compression="gzip")
        preds_g = hf.create_group("predictions")
        preds_g.create_dataset(
            "logcounts", data=np.array(log_counts, dtype=np.float32), compression="gzip"
        )
        preds_g.create_dataset(
            "profile", data=profiles, compression="gzip"
        )
    print(f"Saved h5: {out_h5}")


if __name__ == "__main__":
    main()
