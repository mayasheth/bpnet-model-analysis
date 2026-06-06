#!/usr/bin/env python
"""
Compute DeepLIFT/SHAP attribution scores for a MultiModalBPNet model.

Input X has shape (N, 5, L):
    Channels 0-3: one-hot DNA sequence
    Channel 4:    accessibility signal (normalized)

Attribution scores are split by modality:
    seq_scores  (N, 4, L) - sequence importance
    acc_scores  (N, 1, L) - accessibility importance

Reference sequences for DeepLIFT:
    Channels 0-3: dinucleotide shuffle (standard BPNet approach)
    Channel 4:    random permutation of accessibility values within each sequence

Output HDF5 keys:
    seq_hyp_scores  (N, L, 4) - hypothetical sequence importance
    acc_hyp_scores  (N, L, 1) - hypothetical accessibility importance
    input_seqs      (N, L, 4) - one-hot sequences
    input_acc       (N, L, 1) - accessibility input (normalized)
    coords_chrom    (N,)
    coords_start    (N,)
    coords_end      (N,)
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyfaidx
import pyBigWig
import h5py
import torch

from tangermeme.ersatz import dinucleotide_shuffle
from bpnetlite.attribute import deep_lift_shap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from train_multimodal_bpnet import one_hot_encode, normalize_accessibility

NARROWPEAK_COLS = ["chr", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-m", "--model", required=True, help="Path to .torch model file")
    p.add_argument("-r", "--regions", required=True, help="narrowPeak of regions")
    p.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    p.add_argument("--accessibility-bw", required=True,
                   help="Accessibility BigWig")
    p.add_argument("--acc-stats", required=True,
                   help="JSON with acc_mean and acc_std (from training output)")
    p.add_argument("-o", "--output-prefix", required=True)
    p.add_argument("--chrom", required=True, help="Chromosome to process")
    p.add_argument("--num-shuffles", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=8,
                   help="Sequence-reference pairs per batch (reduce if OOM)")
    p.add_argument("--target", type=int, default=0,
                   help="Model output index to attribute (0=profile, 1=counts)")
    p.add_argument("--device", default="cuda")
    return p.parse_args()


def build_references(X, n_shuffles=20, random_state=None):
    """Build reference tensor for 5-channel multimodal input.

    Channels 0-3 (sequence): dinucleotide shuffle
    Channel 4 (accessibility): random permutation of values

    Parameters
    ----------
    X: torch.Tensor, shape (N, 5, L)

    Returns
    -------
    refs: torch.Tensor, shape (N, n_shuffles, 5, L)
    """
    rng = np.random.RandomState(random_state)
    N, C, L = X.shape
    refs = torch.zeros(N, n_shuffles, C, L, dtype=X.dtype)

    for i in range(N):
        seq_refs = dinucleotide_shuffle(
            X[i:i+1, :4, :], n=n_shuffles,
            random_state=rng.randint(2**31)
        )  # shape: (1, n_shuffles, 4, L) or (n_shuffles, 4, L)
        if seq_refs.dim() == 4:
            seq_refs = seq_refs.squeeze(0)  # (n_shuffles, 4, L)
        refs[i, :, :4, :] = seq_refs

        acc = X[i, 4, :].numpy()
        for j in range(n_shuffles):
            refs[i, j, 4, :] = torch.from_numpy(rng.permutation(acc))

    return refs


def main():
    args = parse_args()

    with open(args.acc_stats) as f:
        stats = json.load(f)
    acc_mean = stats["acc_mean"]
    acc_std = stats["acc_std"]

    model = torch.load(args.model, map_location="cpu")
    model.eval()
    in_window = 2114

    regions_df = pd.read_csv(
        args.regions, sep="\t", header=None,
        names=NARROWPEAK_COLS, usecols=range(10)
    )
    regions_df = regions_df[regions_df["chr"] == args.chrom].reset_index(drop=True)
    print(f"Regions on {args.chrom}: {len(regions_df)}")

    if len(regions_df) == 0:
        print("No regions, exiting.")
        return

    half = in_window // 2
    genome = pyfaidx.Fasta(args.genome)
    acc_bw = pyBigWig.open(args.accessibility_bw)

    seqs_list, accs_list, starts_list, ends_list = [], [], [], []

    for _, row in regions_df.iterrows():
        center = int(row["start"]) + int(row["summit"])
        s, e = center - half, center + half
        chrom = row["chr"]

        if s < 0 or e > len(genome[chrom]):
            continue

        seq_str = str(genome[chrom][s:e])
        if len(seq_str) != in_window:
            continue

        acc = acc_bw.values(chrom, s, e, numpy=True)
        if acc is None or len(acc) != in_window:
            continue
        acc = np.nan_to_num(acc, nan=0.0).astype(np.float32)

        seqs_list.append(one_hot_encode(seq_str))
        accs_list.append(acc)
        starts_list.append(s)
        ends_list.append(e)

    genome.close()
    acc_bw.close()

    if len(seqs_list) == 0:
        print("No valid sequences extracted, exiting.")
        return

    seqs = np.stack(seqs_list)    # (N, 4, L)
    accs = np.stack(accs_list)    # (N, L)

    # Normalize accessibility
    accs = np.log1p(np.clip(accs, 0, None))
    accs = (accs - acc_mean) / (acc_std if acc_std > 1e-8 else 1.0)
    accs = accs[:, np.newaxis, :]  # (N, 1, L)

    X = torch.from_numpy(np.concatenate([seqs, accs], axis=1))  # (N, 5, L)
    print(f"Building references for {len(X)} sequences ({args.num_shuffles} shuffles)...")
    references = build_references(X, n_shuffles=args.num_shuffles, random_state=0)
    # shape: (N, n_shuffles, 5, L)

    print(f"Computing attributions for {len(X)} sequences...")
    attr = deep_lift_shap(
        model, X,
        target=args.target,
        batch_size=args.batch_size,
        references=references,
        hypothetical=True,
        device=args.device,
        verbose=True
    )
    # attr shape: (N, 5, L)

    seq_attr = attr[:, :4, :].numpy()  # (N, 4, L)
    acc_attr = attr[:, 4:, :].numpy()  # (N, 1, L)

    out_h5 = f"{args.output_prefix}.shap_scores.h5"
    os.makedirs(os.path.dirname(out_h5) or ".", exist_ok=True)
    with h5py.File(out_h5, "w") as hf:
        # Transpose to (N, L, C) to match chrombpnet/MoDISCo convention
        hf.create_dataset("seq_hyp_scores",
                          data=seq_attr.transpose(0, 2, 1).astype(np.float16),
                          compression="gzip")
        hf.create_dataset("acc_hyp_scores",
                          data=acc_attr.transpose(0, 2, 1).astype(np.float16),
                          compression="gzip")
        hf.create_dataset("input_seqs",
                          data=seqs.transpose(0, 2, 1).astype(np.int8),
                          compression="gzip")
        hf.create_dataset("input_acc",
                          data=accs.transpose(0, 2, 1).astype(np.float16),
                          compression="gzip")
        hf.create_dataset("coords_chrom",
                          data=np.array([args.chrom] * len(seqs_list), dtype="S"),
                          compression="gzip")
        hf.create_dataset("coords_start",
                          data=np.array(starts_list), compression="gzip")
        hf.create_dataset("coords_end",
                          data=np.array(ends_list), compression="gzip")

    print(f"Saved: {out_h5}")

    with open(f"{args.output_prefix}.DONE.txt", "w") as f:
        f.write("done\n")


if __name__ == "__main__":
    main()
