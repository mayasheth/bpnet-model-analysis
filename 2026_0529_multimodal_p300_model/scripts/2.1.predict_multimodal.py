#!/usr/bin/env python
"""
Predict p300 log counts from MultiModalBPNet models on candidate elements.

Outputs in --output-dir:
  mean_predictions.tsv.gz   mean predictions across all 5 folds (all elements)
  cv_predictions.tsv.gz     held-out test chromosome predictions per fold

Usage:
  pixi run -e multimodal python scripts/2.1.predict_multimodal.py \
    --elements reference/K562_DNase_candidate_elements.narrowPeak \
    --genome /path/to/hg38.fa \
    --signal-plus-bw /path/to/chip_plus.bw \
    --signal-minus-bw /path/to/chip_minus.bw \
    --accessibility-bw data/atac.bw \
    --model-dir models/atac \
    --fold-json reference/hg38_five_folds.json \
    --peaks reference/ENCSR000EGE_peaks_inliers.narrowPeak \
    --output-dir predictions/atac
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyfaidx
import pyBigWig
import torch

# multimodal_bpnet.py lives in the top-level scripts/ dir, two levels up
_PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, os.path.join(_PROJECT_DIR, "scripts"))
from tangermeme.predict import predict as tg_predict

NARROWPEAK_COLS = ["chrom", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--elements", required=True, help="Candidate elements narrowPeak")
    p.add_argument("--genome", required=True)
    p.add_argument("--signal-plus-bw", required=True)
    p.add_argument("--signal-minus-bw", required=True)
    p.add_argument("--accessibility-bw", required=True)
    p.add_argument("--model-dir", required=True,
                   help="Dir containing fold0..fold4/ subdirs with multimodal_bpnet.torch")
    p.add_argument("--fold-json", required=True)
    p.add_argument("--peaks", required=True,
                   help="p300 peaks narrowPeak for overlap annotation")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--in-window", type=int, default=2114)
    p.add_argument("--out-window", type=int, default=1000)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--device", default="cpu")
    return p.parse_args()


def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    enc = np.zeros((4, len(seq)), dtype=np.float32)
    for i, c in enumerate(seq.upper()):
        if c in mapping:
            enc[mapping[c], i] = 1.0
    return enc


def extract_windows(elements_df, genome_fa, signal_plus_bw, signal_minus_bw,
                    accessibility_bw, in_window, out_window):
    """Extract seq, accessibility, and signal windows centered on each element.

    Returns
    -------
    seqs:    (N, 4, in_window) float32
    accs:    (N, 1, in_window) float32  (raw, not normalized)
    signals: (N, 2, out_window) float32  (plus, minus strands)
    coords:  DataFrame(chrom, out_start, out_end, region_name), RangeIndex 0..N-1
    """
    half_in = in_window // 2
    half_out = out_window // 2

    genome = pyfaidx.Fasta(genome_fa)
    plus_bw = pyBigWig.open(signal_plus_bw)
    minus_bw = pyBigWig.open(signal_minus_bw)
    acc_bw = pyBigWig.open(accessibility_bw)

    # Use the intersection of chrom sizes across all sources to avoid invalid intervals
    genome_chroms = set(genome.keys())
    chrom_sizes = {}
    for chrom in genome_chroms:
        sizes = [
            len(genome[chrom]),
            plus_bw.chroms().get(chrom, 0),
            minus_bw.chroms().get(chrom, 0),
            acc_bw.chroms().get(chrom, 0),
        ]
        chrom_sizes[chrom] = min(sizes)

    seqs, accs, signals, coords = [], [], [], []
    n_skipped = 0

    for _, row in elements_df.iterrows():
        chrom = row["chrom"]
        center = (int(row["start"]) + int(row["end"])) // 2

        s_in, e_in = center - half_in, center + half_in
        s_out, e_out = center - half_out, center + half_out

        chrom_len = chrom_sizes.get(chrom, 0)
        if chrom_len == 0 or s_in < 0 or e_in > chrom_len or s_out < 0 or e_out > chrom_len:
            n_skipped += 1
            continue

        seq_str = str(genome[chrom][s_in:e_in])
        if len(seq_str) != in_window:
            n_skipped += 1
            continue

        sig_plus = plus_bw.values(chrom, s_out, e_out, numpy=True)
        sig_minus = minus_bw.values(chrom, s_out, e_out, numpy=True)
        acc = acc_bw.values(chrom, s_in, e_in, numpy=True)

        if (sig_plus is None or len(sig_plus) != out_window
                or sig_minus is None or len(sig_minus) != out_window
                or acc is None or len(acc) != in_window):
            n_skipped += 1
            continue

        seqs.append(one_hot_encode(seq_str))
        accs.append(np.nan_to_num(acc, nan=0.0).astype(np.float32)[np.newaxis, :])
        signals.append(np.stack([
            np.nan_to_num(sig_plus, nan=0.0).astype(np.float32),
            np.nan_to_num(sig_minus, nan=0.0).astype(np.float32),
        ], axis=0))
        coords.append({"chrom": chrom, "out_start": s_out, "out_end": e_out,
                       "region_name": row["name"]})

    genome.close(); plus_bw.close(); minus_bw.close(); acc_bw.close()

    print(f"  Extracted {len(seqs)} windows ({n_skipped} skipped)")
    return (np.stack(seqs), np.stack(accs), np.stack(signals),
            pd.DataFrame(coords))


def compute_peak_overlap(coords_df, peaks_path):
    """Return int array: 1 if region overlaps a p300 peak, 0 otherwise."""
    peaks = pd.read_csv(peaks_path, sep="\t", header=None, usecols=[0, 1, 2],
                        names=["chrom", "start", "end"])
    overlap = np.zeros(len(coords_df), dtype=np.int8)

    for chrom, chrom_peaks in peaks.groupby("chrom"):
        r_mask = coords_df["chrom"].values == chrom
        if r_mask.sum() == 0:
            continue
        r_idx = np.where(r_mask)[0]
        r_starts = coords_df["out_start"].values[r_idx]
        r_ends = coords_df["out_end"].values[r_idx]
        p_starts = chrom_peaks["start"].values
        p_ends = chrom_peaks["end"].values
        # (N_regions, N_peaks) overlap matrix
        has_overlap = ((r_starts[:, None] < p_ends[None, :]) &
                       (r_ends[:, None] > p_starts[None, :])).any(axis=1)
        overlap[r_idx] = has_overlap.astype(np.int8)

    return overlap


def normalize_accessibility(accs, mean, std):
    accs = np.log1p(np.clip(accs, 0, None))
    return (accs - mean) / std


def predict_fold(model_path, seqs, accs_raw, acc_stats, batch_size, device):
    """Load model, normalize accessibility, return predicted log counts (N,)."""
    model = torch.load(model_path, map_location="cpu", weights_only=False)
    accs_norm = normalize_accessibility(
        accs_raw.copy(), acc_stats["acc_mean"], acc_stats["acc_std"]
    )
    X = torch.from_numpy(np.concatenate([seqs, accs_norm], axis=1))  # (N, 5, L)
    preds = tg_predict(model, X, func=lambda out: out[1],
                       batch_size=batch_size, device=device)  # (N, 1)
    return preds.squeeze(1).numpy()  # (N,)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.fold_json) as f:
        fold_json = json.load(f)
    # Map chrom -> fold index for CV (using "test" split)
    chrom_to_cv_fold = {}
    for k, v in fold_json.items():
        for c in v["test"]:
            chrom_to_cv_fold[c] = int(k)

    print("Loading candidate elements...")
    elements = pd.read_csv(args.elements, sep="\t", header=None,
                           usecols=range(10), names=NARROWPEAK_COLS)
    print(f"  {len(elements)} elements")

    print("Extracting windows (seq + accessibility + signal)...")
    seqs, accs_raw, signals, coords = extract_windows(
        elements, args.genome,
        args.signal_plus_bw, args.signal_minus_bw, args.accessibility_bw,
        args.in_window, args.out_window,
    )

    true_logcounts = np.log1p(signals.sum(axis=(1, 2)))  # log1p(total counts)

    print("Computing p300 peak overlap...")
    peak_overlap = compute_peak_overlap(coords, args.peaks)

    # Per-fold predictions
    fold_preds = np.zeros((5, len(seqs)), dtype=np.float32)
    for fold_idx in range(5):
        model_path = os.path.join(args.model_dir, f"fold{fold_idx}",
                                  "multimodal_bpnet.torch")
        stats_path = os.path.join(args.model_dir, f"fold{fold_idx}",
                                  "acc_normalization_stats.json")
        with open(stats_path) as f:
            acc_stats = json.load(f)
        print(f"Predicting fold {fold_idx} "
              f"(acc_mean={acc_stats['acc_mean']:.4f}, "
              f"acc_std={acc_stats['acc_std']:.4f})...")
        fold_preds[fold_idx] = predict_fold(
            model_path, seqs, accs_raw, acc_stats, args.batch_size, args.device
        )

    # Assemble mean predictions TSV
    mean_df = coords.rename(columns={"out_start": "start", "out_end": "end"}).copy()
    mean_df["true_logcounts"] = true_logcounts
    mean_df["mean_pred_logcounts"] = fold_preds.mean(axis=0)
    mean_df["EP300_peak_overlap"] = peak_overlap
    for i in range(5):
        mean_df[f"pred_fold{i}"] = fold_preds[i]

    out_cols = ["chrom", "start", "end", "region_name", "true_logcounts",
                "mean_pred_logcounts", "EP300_peak_overlap"] + \
               [f"pred_fold{i}" for i in range(5)]
    mean_df[out_cols].to_csv(
        os.path.join(args.output_dir, "mean_predictions.tsv.gz"),
        sep="\t", index=False, compression="gzip",
    )
    print(f"Saved mean_predictions.tsv.gz ({len(mean_df)} elements)")

    # Assemble CV predictions TSV
    cv_fold_per_region = np.array(
        [chrom_to_cv_fold.get(c, -1) for c in coords["chrom"].values]
    )
    cv_mask = cv_fold_per_region >= 0
    cv_idx = np.where(cv_mask)[0]
    cv_folds = cv_fold_per_region[cv_idx]

    cv_df = coords.iloc[cv_idx].rename(
        columns={"out_start": "start", "out_end": "end"}
    ).copy()
    cv_df["true_logcounts"] = true_logcounts[cv_idx]
    cv_df["pred_logcounts"] = fold_preds[cv_folds, cv_idx]
    cv_df["EP300_peak_overlap"] = peak_overlap[cv_idx]
    cv_df["fold"] = cv_folds

    cv_df[["chrom", "start", "end", "region_name", "true_logcounts",
           "pred_logcounts", "EP300_peak_overlap", "fold"]].to_csv(
        os.path.join(args.output_dir, "cv_predictions.tsv.gz"),
        sep="\t", index=False, compression="gzip",
    )
    print(f"Saved cv_predictions.tsv.gz ({len(cv_df)} elements)")


if __name__ == "__main__":
    main()
