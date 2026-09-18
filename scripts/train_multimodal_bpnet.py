#!/usr/bin/env python
"""
Train a MultiModalBPNet model on p300 ChIP-seq peaks + GC-matched negatives.

Supports three input modes via --mode:
  multimodal  5-channel input: one-hot DNA (4ch) + accessibility profile (1ch)
  sequence    4-channel input: one-hot DNA only  (equivalent to standard BPNet)
  atac        1-channel input: base-pair accessibility profile only

Accessibility is log1p-normalized and standardized using the mean and standard
deviation computed from training peaks. Pass --acc-mean / --acc-std to reuse
pre-computed statistics (required when running multiple folds for consistency).
--genome is not required in atac mode.

Stranded and unstranded targets are both supported. Pass --signal-minus-bw for
stranded data (n_outputs=2, the p300 setup); omit it for unstranded data such as
histone-mark ChIP-seq (n_outputs=1). --accessibility-bw may be omitted in
sequence mode, where the model never reads it.

Usage:
    pixi run -e multimodal python scripts/train_multimodal_bpnet.py \\
        --mode [multimodal|sequence|atac] \\
        --peaks reference/ENCSR000EGE_peaks_inliers.narrowPeak \\
        --negatives reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed \\
        --genome /path/to/hg38.fa \\        # not required for --mode atac
        --signal-plus-bw /path/to/chip_plus.bw \\
        --signal-minus-bw /path/to/chip_minus.bw \\   # omit for unstranded
        --accessibility-bw /path/to/atac.bw \\        # not required for --mode sequence
        --fold reference/hg38_five_folds.json \\
        --fold-key 0 \\
        --output-dir models/atac/fold0
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multimodal_bpnet import MultiModalBPNet

NARROWPEAK_COLS = ["chr", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", default="multimodal",
                   choices=["multimodal", "sequence", "atac"],
                   help="Input mode: multimodal (DNA+ATAC), sequence (DNA only), "
                        "or atac (ATAC profile only). Default: multimodal")
    p.add_argument("--peaks", required=True, help="narrowPeak of training peaks")
    p.add_argument("--negatives", required=True,
                   help="BED of GC-matched negative regions (col 1-3 used)")
    p.add_argument("--genome", default=None,
                   help="Reference genome FASTA (not required for --mode atac)")
    p.add_argument("--signal-plus-bw", required=True,
                   help="ChIP-seq signal BigWig, plus strand (target)")
    p.add_argument("--signal-minus-bw", default=None,
                   help="ChIP-seq signal BigWig, minus strand (target)")
    p.add_argument("--profile-loss-weight", type=float, default=1.0,
                   help="Multiplies the profile term in the loss. Default 1.0 reproduces "
                        "every earlier run. Set it when --profile-target-plus-bw points "
                        "the profile head at a track of different read depth: MNLL scales "
                        "with depth, so an unweighted swap to a 17.8x deeper target "
                        "divides the effective --count-loss-weight by about 18.")
    p.add_argument("--profile-target-plus-bw", default=None,
                   help="Separate target for the PROFILE head, plus strand. The counts "
                        "head keeps --signal-plus-bw. Use to give the profile head a "
                        "learnable auxiliary task: H3K27ac's 1 bp inter-replicate ceiling "
                        "is 0.21 so that head trains on noise, while K562 DNase's is 0.848. "
                        "Needs no second model and nothing extra at inference.")
    p.add_argument("--profile-target-minus-bw", default=None,
                   help="Minus strand of --profile-target-plus-bw. Omit for an unstranded "
                        "profile target; n_outputs then follows the PROFILE target.")
    p.add_argument("--accessibility-bw", default=None,
                   # comma-separate for fragment-size-stratified channels, e.g.
                   #   --accessibility-bw all.bw,sub.bw,mono.bw,di.bw
                   # Channel ORDER is part of the model contract: prediction must pass
                   # the same tracks in the same order, or the learned filters are
                   # applied to the wrong inputs. The order used is recorded in
                   # training_target.json.
                   help="Accessibility BigWig (ATAC or DNase)")
    p.add_argument("--cell-types-json", default=None,
                   help="Train on SEVERAL cell types. JSON mapping a cell-type name to an "
                        "object with peaks, signal_plus_bw, and optionally signal_minus_bw, "
                        "accessibility_bw and negatives; anything omitted falls back to the "
                        "matching single-cell-type flag. When this is given, --peaks and "
                        "--signal-plus-bw are ignored. Windows from every cell type are "
                        "pooled into one training set and the chromosome-fold split is "
                        "applied within each, so a held-out chromosome is held out "
                        "EVERYWHERE and no cell type can leak it.")
    p.add_argument("--holdout-cell", default=None,
                   help="Name of a cell type in --cell-types-json to EXCLUDE from training, "
                        "for leave-one-out transfer experiments. Named rather than indexed "
                        "so the run is self-describing in training_target.json.")
    p.add_argument("--fold", required=True, help="Fold JSON file")
    p.add_argument("--fold-key", default="0", help="Key within fold JSON (default: 0)")
    p.add_argument("--output-dir", required=True, help="Output directory for model")
    p.add_argument("--in-window", type=int, default=2114)
    p.add_argument("--out-window", type=int, default=1000)
    p.add_argument("--max-jitter", type=int, default=50)
    p.add_argument("--n-filters", type=int, default=64)
    p.add_argument("--n-acc-filters", type=int, default=8)
    p.add_argument("--n-layers", type=int, default=8)
    p.add_argument("--count-loss-weight", type=float, default=1.0)
    p.add_argument("--gate-accessibility", action="store_true",
                   help="Let sequence gate the accessibility branch: "
                        "X_acc *= sigmoid(conv(X_seq)). Initialised open, so the model "
                        "starts identical to the ungated one and can only learn to suppress "
                        "accessibility where that helps. Multimodal mode only.")
    p.add_argument("--overprediction-weight", type=float, default=1.0,
                   help="Multiply squared log-count error by this factor where the "
                        "prediction EXCEEDS the truth. 1.0 reproduces bpnetlite exactly. "
                        "Pair it with --gate-accessibility: the gate has no gradient "
                        "pressure to close while the loss is symmetric.")
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-epochs", type=int, default=100)
    p.add_argument("--early-stopping", type=int, default=10)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--negative-ratio", type=float, default=0.1)
    p.add_argument("--count-offset-model", default=None,
                   help="Directory of a trained ATAC-only model (containing fold*/). Its "
                        "held-out predicted logcounts are used as a per-region OFFSET on "
                        "the count loss, so this run learns the RESIDUAL beyond "
                        "accessibility: observed - atac_pred. Final prediction at "
                        "inference is model_output + atac_pred. Requires "
                        "--accessibility-bw even in sequence mode, since the offset "
                        "model reads it. Omit for normal training.")
    p.add_argument("--max-negatives", type=int, default=None,
                   help="Cap on distinct negative windows held in memory. Default "
                        "None = 10x the training peak count (historical behaviour). "
                        "Set this when training on a large peak set: 10x of ~120k "
                        "elements is ~1.2M windows and will OOM. The sampler only "
                        "draws --negative-ratio of each batch from negatives, so a "
                        "much smaller pool is usually sufficient.")
    p.add_argument("--acc-mean", type=float, default=None,
                   help="Pre-computed log1p accessibility mean for standardization")
    p.add_argument("--acc-std", type=float, default=None,
                   help="Pre-computed log1p accessibility std for standardization")
    p.add_argument("--device", default="cuda")
    p.add_argument("--n-workers", type=int, default=4)
    return p.parse_args()


def load_peaks(path, chroms=None):
    df = pd.read_csv(path, sep="\t", header=None, usecols=range(10),
                     names=NARROWPEAK_COLS)
    if chroms is not None:
        df = df[df["chr"].isin(chroms)].reset_index(drop=True)
    return df


def load_negatives(path, chroms=None):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chr", "start", "end"])
    if chroms is not None:
        df = df[df["chr"].isin(chroms)].reset_index(drop=True)
    return df


def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    enc = np.zeros((4, len(seq)), dtype=np.float32)
    for i, c in enumerate(seq.upper()):
        if c in mapping:
            enc[mapping[c], i] = 1.0
    return enc


def extract_windows(regions_df, genome_fa, signal_plus_bw, signal_minus_bw,
                    accessibility_bw, in_window, out_window, max_jitter,
                    is_peak=True, profile_plus_bw=None, profile_minus_bw=None):
    """Extract sequence, signal, and accessibility windows around each region center.

    genome_fa may be None for atac mode (seq extraction is skipped; seqs returned
    as empty array). signal_minus_bw may be None for unstranded targets, giving a
    1-channel signal array. accessibility_bw may be None in sequence mode, giving a
    zero accs array the model never reads.

    profile_plus_bw OPTIONALLY EXTRACTS A SECOND TARGET, in this same loop and over the
    same window bounds, and the return becomes a 5-tuple. Doing it here rather than in a
    second call is the whole point: a region that any one track cannot serve is skipped,
    and a second call would build its own `valid` mask, so a track with one short contig
    would return a different number of rows and silently pair every row after that point
    with the wrong region. Ten existing callers pass nothing and keep the 4-tuple.

    Returns
    -------
    seqs:    np.ndarray, (N, 4, in_window + 2*max_jitter)  — zeros if genome_fa is None
    signals: np.ndarray, (N, C, out_window + 2*max_jitter) — C=2 stranded, 1 unstranded
    accs:    np.ndarray, (N, C, in_window + 2*max_jitter)  — C = number of comma-
             separated accessibility tracks; a single zero channel if none given
    valid:   bool array
    profiles: np.ndarray, same shape convention as signals — ONLY when profile_plus_bw
             is given, in which case this is a 5-tuple
    """
    jitter = max_jitter if is_peak else 0
    half_in = (in_window + 2 * jitter) // 2
    half_out = (out_window + 2 * jitter) // 2

    use_seq = genome_fa is not None
    genome = pyfaidx.Fasta(genome_fa) if use_seq else None
    plus_bw = pyBigWig.open(signal_plus_bw)
    stranded = signal_minus_bw is not None
    minus_bw = pyBigWig.open(signal_minus_bw) if stranded else None
    use_prof = profile_plus_bw is not None
    prof_plus_bw = pyBigWig.open(profile_plus_bw) if use_prof else None
    prof_stranded = use_prof and profile_minus_bw is not None
    prof_minus_bw = pyBigWig.open(profile_minus_bw) if prof_stranded else None
    # accessibility_bw may be a comma-separated list -> one input channel per track
    acc_paths = ([p for p in accessibility_bw.split(",") if p]
                 if accessibility_bw else [])
    use_acc = len(acc_paths) > 0
    acc_bws = [pyBigWig.open(p) for p in acc_paths]
    n_acc = max(1, len(acc_bws))
    # chrom sizes come from the accessibility track only when there is no genome
    acc_chrom_sizes = dict(acc_bws[0].chroms()) if use_acc else dict(plus_bw.chroms())

    seqs, signals, accs, profiles = [], [], [], []
    valid = np.zeros(len(regions_df), dtype=bool)

    for idx, row in regions_df.iterrows():
        chrom = row["chr"]

        if is_peak:
            center = int(row["start"]) + int(row["summit"])
        else:
            center = (int(row["start"]) + int(row["end"])) // 2

        s_in = center - half_in
        e_in = center + half_in
        s_out = center - half_out
        e_out = center + half_out

        # Bounds check against genome or BigWig chrom sizes
        if use_seq:
            chrom_len = len(genome[chrom])
        else:
            chrom_len = acc_chrom_sizes.get(chrom, 0)
        if s_in < 0 or e_in > chrom_len or s_out < 0 or e_out > chrom_len:
            continue

        if use_seq:
            seq_str = str(genome[chrom][s_in:e_in])
            if len(seq_str) != in_window + 2 * jitter:
                continue

        sig_plus = plus_bw.values(chrom, s_out, e_out, numpy=True)
        if sig_plus is None or len(sig_plus) != out_window + 2 * jitter:
            continue
        sig_plus = np.nan_to_num(sig_plus, nan=0.0).astype(np.float32)

        if stranded:
            sig_minus = minus_bw.values(chrom, s_out, e_out, numpy=True)
            if sig_minus is None or len(sig_minus) != out_window + 2 * jitter:
                continue
            sig_minus = np.nan_to_num(sig_minus, nan=0.0).astype(np.float32)

        if use_prof:
            pp = prof_plus_bw.values(chrom, s_out, e_out, numpy=True)
            if pp is None or len(pp) != out_window + 2 * jitter:
                continue
            pp = np.nan_to_num(pp, nan=0.0).astype(np.float32)
            if prof_stranded:
                pm = prof_minus_bw.values(chrom, s_out, e_out, numpy=True)
                if pm is None or len(pm) != out_window + 2 * jitter:
                    continue
                pm = np.nan_to_num(pm, nan=0.0).astype(np.float32)

        if use_acc:
            chans, bad_acc = [], False
            for bw in acc_bws:
                a = bw.values(chrom, s_in, e_in, numpy=True)
                if a is None or len(a) != in_window + 2 * jitter:
                    bad_acc = True
                    break
                chans.append(np.nan_to_num(a, nan=0.0).astype(np.float32))
            if bad_acc:
                continue
            acc = np.stack(chans, axis=0)                       # (C, L)
        else:
            acc = np.zeros((1, in_window + 2 * jitter), dtype=np.float32)

        if use_seq:
            seqs.append(one_hot_encode(seq_str))
        signals.append(np.stack([sig_plus, sig_minus], axis=0) if stranded
                       else sig_plus[np.newaxis, :])
        accs.append(acc)
        if use_prof:
            profiles.append(np.stack([pp, pm], axis=0) if prof_stranded
                            else pp[np.newaxis, :])
        valid[idx] = True

    if use_seq:
        genome.close()
    plus_bw.close()
    if stranded:
        minus_bw.close()
    if use_prof:
        prof_plus_bw.close()
        if prof_stranded:
            prof_minus_bw.close()
    for bw in acc_bws:
        bw.close()

    n = len(accs)
    L_in = in_window + 2 * jitter
    L_out = out_window + 2 * jitter
    if n == 0:
        empty = (np.zeros((0, 4, L_in), dtype=np.float32),
                 np.zeros((0, 2 if signal_minus_bw is not None else 1, L_out),
                          dtype=np.float32),
                 np.zeros((0, n_acc, L_in), dtype=np.float32),
                 valid)
        if use_prof:
            return empty + (np.zeros((0, 2 if prof_stranded else 1, L_out),
                                     dtype=np.float32),)
        return empty

    seqs_arr = np.stack(seqs) if use_seq else np.zeros((n, 4, L_in), dtype=np.float32)
    out = (seqs_arr, np.stack(signals), np.stack(accs), valid)
    if use_prof:
        return out + (np.stack(profiles),)
    return out


def normalize_accessibility(acc, mean=None, std=None):
    """Apply log1p then standardize accessibility signal."""
    acc = np.log1p(np.clip(acc, 0, None))
    if mean is None:
        mean = acc.mean()
    if std is None:
        std = acc.std()
    if std < 1e-8:
        std = 1.0
    return (acc - mean) / std, mean, std


class MultiModalPeakNegativeSampler(torch.utils.data.Dataset):
    """Dataset yielding model inputs shaped according to mode.

    mode='multimodal': Xi shape (5, L) — cat([seq, acc])
    mode='sequence':   Xi shape (4, L) — seq only
    mode='atac':       Xi shape (1, L) — acc only

    Reverse complement augmentation:
    - Sequence channels: flip both channel and position dims
    - Accessibility channel: flip position dim only
    - Signal (target): swap plus/minus and flip positions
    """

    def __init__(self, peak_seqs, peak_accs, peak_signals,
                 neg_seqs, neg_accs, neg_signals,
                 negative_ratio=0.1, in_window=2114, out_window=1000,
                 max_jitter=0, reverse_complement=False, random_state=None,
                 mode='multimodal', peak_offsets=None, neg_offsets=None,
                 peak_profiles=None, neg_profiles=None):
        self.peak_seqs = peak_seqs
        self.peak_accs = peak_accs
        self.peak_signals = peak_signals
        self.neg_seqs = neg_seqs
        self.neg_accs = neg_accs
        self.neg_signals = neg_signals

        # Separate target for the profile head; None = both heads read peak_signals.
        self.peak_profiles = peak_profiles
        self.neg_profiles = neg_profiles
        self.use_profile_target = peak_profiles is not None
        if self.use_profile_target:
            if neg_profiles is None:
                raise ValueError("peak_profiles given without neg_profiles")
            # The profile term is masked to labels == 1 so negatives never reach it, but
            # they still have to collate to the right shape, and an all-zero stand-in
            # would be a silent lie if that masking ever changed.
            if len(peak_profiles) != len(peak_seqs) or len(neg_profiles) != len(neg_seqs):
                raise ValueError(
                    f"profile target rows do not match: peaks {len(peak_profiles)} vs "
                    f"{len(peak_seqs)}, negatives {len(neg_profiles)} vs {len(neg_seqs)}")

        self.n_peaks = len(peak_seqs)
        self.n_negatives = len(neg_seqs)
        self.negative_ratio = negative_ratio
        self.negative_likelihood = 1 / (1 + 1 / negative_ratio)

        self.in_window = in_window
        self.out_window = out_window
        self.max_jitter = max_jitter
        self.reverse_complement = reverse_complement
        self.mode = mode

        # Per-region count offsets for residual training; None = normal training.
        self.peak_offsets = peak_offsets
        self.neg_offsets = neg_offsets
        self.use_offsets = peak_offsets is not None

        self.rng = np.random.RandomState(random_state)
        self.n_peaks_seen = 0
        self.peak_ordering = self.rng.permutation(self.n_peaks)

    def __len__(self):
        return self.n_peaks + int(self.n_peaks * self.negative_ratio)

    def __getitem__(self, idx):
        if idx == 0:
            self.peak_ordering = np.arange(self.n_peaks)
            self.rng.shuffle(self.peak_ordering)

        if self.rng.uniform() >= self.negative_likelihood:
            i = self.peak_ordering[self.n_peaks_seen % self.n_peaks]
            jitter = self.rng.randint(self.max_jitter * 2) if self.max_jitter > 0 else 0
            label = 1
            seqs, accs, signals = self.peak_seqs, self.peak_accs, self.peak_signals
            offsets = self.peak_offsets
            prof = self.peak_profiles
            self.n_peaks_seen += 1
        else:
            i = self.rng.randint(self.n_negatives)
            jitter = 0
            label = 0
            seqs, accs, signals = self.neg_seqs, self.neg_accs, self.neg_signals
            offsets = self.neg_offsets
            prof = self.neg_profiles

        Xi_seq = torch.from_numpy(seqs[i][:, jitter:jitter + self.in_window])
        Xi_acc = torch.from_numpy(accs[i][:, jitter:jitter + self.in_window])
        yi = torch.from_numpy(signals[i][:, jitter:jitter + self.out_window])
        # Same `jitter` scalar as the sequence and accessibility crops above, which is
        # what keeps the second target aligned with the first.
        yp = (torch.from_numpy(prof[i][:, jitter:jitter + self.out_window])
              if self.use_profile_target else None)

        if self.reverse_complement and self.rng.randint(2) == 1:
            Xi_seq = torch.flip(Xi_seq, [0, 1])
            Xi_acc = torch.flip(Xi_acc, [1])
            yi = torch.flip(yi, [0, 1])
            # Flipped as its own block. Concatenating the two targets into one tensor and
            # flipping dim 0 once would reverse the block order too, handing the profile
            # head the counts target.
            if yp is not None:
                yp = torch.flip(yp, [0, 1])

        if self.mode == 'multimodal':
            Xi = torch.cat([Xi_seq, Xi_acc], dim=0)  # (5, in_window)
        elif self.mode == 'sequence':
            Xi = Xi_seq                               # (4, in_window)
        elif self.mode == 'atac':
            Xi = Xi_acc                               # (1, in_window)

        # Layouts, all read by fit() from an EXPLICIT profile_target flag rather than
        # from arity, since (X, offset, y, label) and (X, y, label, y_profile) are both
        # length 4 and sniffing would read the label column as a target:
        #   (X, y, label)                      plain
        #   (X, offset, y, label)              residual training
        #   (X, y, label, y_profile)           separate profile target
        #   (X, offset, y, label, y_profile)   both
        if self.use_offsets:
            off = torch.tensor(offsets[i], dtype=torch.float32)
            if yp is not None:
                return Xi, off, yi, label, yp
            return Xi, off, yi, label
        if yp is not None:
            return Xi, yi, label, yp
        return Xi, yi, label


def compute_count_offsets(offset_model_dir, fold_key, accs_list, batch_size, device):
    """Per-region logcount predictions from a trained ATAC-only model.

    Used as the offset for residual training. Accessibility is normalized with the
    OFFSET MODEL's own saved statistics, not this run's, or the offset model would be
    fed inputs on a different scale than it was trained on.

    accs_list: list of RAW (unnormalized) accessibility arrays, one per region set.
    Returns a list of 1-D offset arrays in the same order.
    """
    model_path = os.path.join(offset_model_dir, f"fold{fold_key}", "multimodal_bpnet.torch")
    stats_path = os.path.join(offset_model_dir, f"fold{fold_key}",
                              "acc_normalization_stats.json")
    for pth in (model_path, stats_path):
        if not os.path.exists(pth):
            raise SystemExit(f"error: --count-offset-model missing {pth}")
    with open(stats_path) as f:
        st = json.load(f)
    model = torch.load(model_path, map_location="cpu", weights_only=False)
    if not hasattr(model, "mode"):
        model.mode = "atac"
    if model.mode != "atac":
        raise SystemExit(f"error: offset model mode is '{model.mode}', expected 'atac'. "
                         "The offset must come from an accessibility-only model, "
                         "otherwise the residual is not 'beyond accessibility'.")
    model = model.to(device).eval()

    out = []
    for accs in accs_list:
        normed = normalize_accessibility(accs.copy(), mean=st["acc_mean"],
                                        std=st["acc_std"])[0]
        preds = []
        with torch.no_grad():
            for i in range(0, len(normed), batch_size):
                xb = torch.from_numpy(normed[i:i + batch_size]).to(device).float()
                _, lc = model(xb)
                preds.append(lc.squeeze(-1).cpu().numpy())
        out.append(np.concatenate(preds) if preds else np.zeros(0, dtype=np.float32))
    model.to("cpu")
    return out


def main():
    args = parse_args()

    if args.mode != 'atac' and args.genome is None:
        raise ValueError("--genome is required for --mode sequence and --mode multimodal")

    if args.mode in ('multimodal', 'atac') and args.accessibility_bw is None:
        raise ValueError("--accessibility-bw is required for --mode multimodal and "
                         "--mode atac")

    if args.profile_target_minus_bw is not None and args.profile_target_plus_bw is None:
        raise ValueError("--profile-target-minus-bw needs --profile-target-plus-bw")

    if args.count_offset_model is not None and args.accessibility_bw is None:
        raise ValueError("--count-offset-model requires --accessibility-bw, since the "
                         "offset model reads accessibility even when this model does not")

    # Record what this model was trained against. Without this a later run cannot tell
    # whether an offset model was fit to the same target, which silently produced an
    # uninterpretable residual result on 2026-08-26 (see .living/learnings.md).
    target_record = {
        "signal_plus_bw": args.signal_plus_bw,
        "signal_minus_bw": args.signal_minus_bw,
        "profile_target_plus_bw": args.profile_target_plus_bw,
        "profile_target_minus_bw": args.profile_target_minus_bw,
        "accessibility_bw": args.accessibility_bw,
        "accessibility_channels": (
            [p for p in args.accessibility_bw.split(",") if p]
            if args.accessibility_bw else []),
        "peaks": args.peaks,
        "mode": args.mode,
        "in_window": args.in_window,
        "out_window": args.out_window,
        "n_layers": args.n_layers,
        "count_loss_weight": args.count_loss_weight,
        "profile_loss_weight": args.profile_loss_weight,
        "count_offset_model": args.count_offset_model,
    }

    if args.count_offset_model is not None:
        # The offset must be an accessibility prediction OF THIS TARGET. A model fit to a
        # differently-processed signal (e.g. fragment-extended vs 5'-end, a ~250x scale
        # difference) trains without error because correlation is shift-invariant, but the
        # residual it produces is not "beyond accessibility" and is not comparable to
        # anything. Fail loudly instead.
        rec_path = os.path.join(args.count_offset_model, f"fold{args.fold_key}",
                                "training_target.json")
        if not os.path.exists(rec_path):
            raise SystemExit(
                f"error: offset model has no training_target.json at {rec_path}, so its "
                "target cannot be verified against this run's. Retrain the offset model "
                "with the current trainer, or backfill the record with "
                "2026_0824_H3K27ac_model/scripts/1.5.backfill_training_target.py.")
        with open(rec_path) as f:
            off_rec = json.load(f)
        for key in ("signal_plus_bw", "signal_minus_bw", "out_window"):
            mine, theirs = target_record[key], off_rec.get(key)
            if mine != theirs:
                raise SystemExit(
                    f"error: offset model target mismatch on '{key}'.\n"
                    f"  this run:     {mine}\n"
                    f"  offset model: {theirs}\n"
                    "The offset model must be trained on the same signal and window. "
                    "Train an ATAC-only model on this target first.")
        if off_rec.get("mode") != "atac":
            raise SystemExit(
                f"error: offset model mode is '{off_rec.get('mode')}', expected 'atac'. "
                "The offset must come from an accessibility-only model.")
        print(f"Offset model target verified against {rec_path}")

    # Stranded targets give 2 output tracks (the p300 setup); unstranded give 1.
    # n_outputs sizes the PROFILE head, so with a separate profile target it follows that
    # target's strandedness, not the counts target's. The counts head is a single scalar
    # either way.
    use_prof_target = args.profile_target_plus_bw is not None
    if use_prof_target:
        n_outputs = 2 if args.profile_target_minus_bw is not None else 1
        print(f"Profile head target: {args.profile_target_plus_bw}"
              + (f" + {args.profile_target_minus_bw}"
                 if args.profile_target_minus_bw else "")
              + f" -> n_outputs={n_outputs}")
        print(f"Counts head target:  {args.signal_plus_bw} (unchanged)")
    else:
        n_outputs = 2 if args.signal_minus_bw is not None else 1
        print(f"Target is {'stranded' if n_outputs == 2 else 'unstranded'} "
              f"-> n_outputs={n_outputs}")

    with open(args.fold) as f:
        fold_data = json.load(f)[str(args.fold_key)]
    train_chroms = fold_data["train"]
    val_chroms = fold_data["val"]

    os.makedirs(args.output_dir, exist_ok=True)
    model_prefix = os.path.join(args.output_dir, "multimodal_bpnet")

    # ---- multi-cell-type panel -----------------------------------------------------
    # PANEL is an ordered mapping name -> spec. The single-cell-type path is expressed as
    # a one-entry panel so there is exactly one extraction code path, and the arrays it
    # produces for a one-cell panel are identical to what the old code produced.
    if args.cell_types_json:
        with open(args.cell_types_json) as _f:
            _panel = json.load(_f)
        if args.holdout_cell:
            if args.holdout_cell not in _panel:
                raise SystemExit(f"error: --holdout-cell {args.holdout_cell!r} is not in "
                                 f"{args.cell_types_json}; have {sorted(_panel)}")
            _panel = {k: v for k, v in _panel.items() if k != args.holdout_cell}
            if not _panel:
                raise SystemExit("error: holding that cell out leaves nothing to train on")
        PANEL = _panel
    else:
        PANEL = {"_single": {}}

    def _spec(spec, key, fallback):
        return spec.get(key, fallback)

    print(f"Panel: {len(PANEL)} cell type(s): {', '.join(PANEL)}"
          + (f"  (holding out {args.holdout_cell})" if args.holdout_cell else ""))

    # Provenance is written AFTER the panel is resolved, so the file names every cell type
    # the model actually trained on and every track it read. The single-cell-type shape is
    # unchanged, so existing readers (1.5's backfill, 4.1's mode sniffing) still work.
    if args.cell_types_json:
        target_record["cell_types_json"] = args.cell_types_json
        target_record["holdout_cell"] = args.holdout_cell
        target_record["trained_on"] = sorted(PANEL)
        target_record["panel"] = {
            cell: {k: _spec(spec, k, getattr(args, k, None))
                   for k in ("peaks", "signal_plus_bw", "signal_minus_bw",
                             "accessibility_bw", "negatives")}
            for cell, spec in PANEL.items()}

    with open(os.path.join(args.output_dir, "training_target.json"), "w") as f:
        json.dump(target_record, f, indent=2)

    if args.cell_types_json and args.count_offset_model:
        raise SystemExit("error: --count-offset-model with --cell-types-json is not "
                         "implemented. The offset model is trained on one cell type's "
                         "accessibility, so pooling would mix offset scales silently.")

    genome = args.genome  # None for atac mode -- extract_windows handles this
    prof_kw = dict(profile_plus_bw=args.profile_target_plus_bw,
                   profile_minus_bw=args.profile_target_minus_bw)

    def _extract(spec, regions, jitter, is_peak, label):
        out = extract_windows(
            regions, genome,
            _spec(spec, "signal_plus_bw", args.signal_plus_bw),
            _spec(spec, "signal_minus_bw", args.signal_minus_bw),
            _spec(spec, "accessibility_bw", args.accessibility_bw),
            args.in_window, args.out_window, jitter,
            is_peak=is_peak, **prof_kw
        )
        seqs, sigs, accs, _ = out[:4]
        prof = out[4] if use_prof_target else None
        print(f"  Extracted {len(accs)} {label}")
        if prof is not None:
            # Both come out of one loop over one `valid` mask, so a mismatch here means
            # the extractor itself is broken rather than a track being short.
            assert len(prof) == len(sigs), (
                f"{label}: {len(prof)} profile-target rows against {len(sigs)} counts rows")
        return seqs, sigs, accs, prof

    bins = {k: [] for k in ("tr_seqs", "tr_sigs", "tr_accs", "tr_prof",
                            "neg_seqs", "neg_sigs", "neg_accs", "neg_prof",
                            "val_seqs", "val_sigs", "val_accs", "val_prof")}
    acc_stats, per_cell_counts = {}, {}

    for cell, spec in PANEL.items():
        tag = "" if cell == "_single" else f"[{cell}] "
        pk_path = _spec(spec, "peaks", args.peaks)
        ng_path = _spec(spec, "negatives", args.negatives)
        acc_bw = _spec(spec, "accessibility_bw", args.accessibility_bw)

        print(f"{tag}Loading peaks...")
        c_train_peaks = load_peaks(pk_path, train_chroms)
        c_val_peaks = load_peaks(pk_path, val_chroms)
        print(f"  Train peaks: {len(c_train_peaks)}, Val peaks: {len(c_val_peaks)}")

        print(f"{tag}Loading negatives...")
        c_negs = load_negatives(ng_path, train_chroms)
        # Cap the negative pool to bound memory when loading genome windows. The cap is
        # PER CELL TYPE, so a panel of N cell types draws N times this many; that is
        # intended, since each cell type needs its own negatives against its own tracks.
        max_negs = (args.max_negatives if args.max_negatives is not None
                    else len(c_train_peaks) * 10)
        if len(c_negs) > max_negs:
            c_negs = c_negs.sample(max_negs, random_state=42).reset_index(drop=True)
        print(f"  Train negatives: {len(c_negs)} (capped at {max_negs})")
        if args.max_negatives is None and len(c_negs) > 200_000:
            est_gb = (len(c_negs) * args.in_window * 6 * 4) / 1e9
            print(f"  WARNING: --max-negatives was not set, so the pool defaulted to "
                  f"10x peaks = {len(c_negs):,} windows (~{est_gb:.0f} GB of extracted "
                  f"arrays). Set --max-negatives explicitly; the sampler only draws "
                  f"{args.negative_ratio:.0%} of each batch from negatives.")

        print(f"{tag}Extracting windows...")
        a_seqs, a_sigs, a_accs, a_prof = _extract(spec, c_train_peaks, args.max_jitter,
                                                  True, "training peaks")
        n_seqs, n_sigs, n_accs, n_prof = _extract(spec, c_negs, 0, False,
                                                  "training negatives")
        v_seqs, v_sigs, v_accs, v_prof = _extract(spec, c_val_peaks, 0, True,
                                                  "validation peaks")

        if use_prof_target:
            for nm, a, b in ((f"{tag}training peaks", a_prof, a_sigs),
                             (f"{tag}validation peaks", v_prof, v_sigs)):
                if float(a.sum()) <= 0:
                    raise SystemExit(f"error: profile target is all zero over {nm}; check "
                                     f"--profile-target-plus-bw covers these regions")
                print(f"  {nm}: profile-target total {a.sum():,.0f} against counts-target "
                      f"total {b.sum():,.0f}")

        # ACCESSIBILITY IS NORMALIZED PER CELL TYPE, BEFORE POOLING, AND THAT IS THE WHOLE
        # POINT. Each cell type's ATAC has its own library depth and its own signal
        # distribution. Standardizing the pooled array with one mean and std would leave
        # those differences in the input, and the model would read library rather than
        # biology -- which is exactly the mechanism F-017 caught driving a transfer
        # collapse. Per-cell-type standardization puts every cell type on mean 0, sd 1.
        # CONSEQUENCE AT INFERENCE: there is no single training statistic to reuse, so a
        # multi-cell-type model MUST be predicted with statistics derived from the target
        # cell type's own track. That is `4.1 --acc-mean/--acc-std`, fed by `4.26`.
        if args.count_offset_model is not None:
            # Residual training needs RAW accessibility, and normalization below replaces
            # these arrays. Copy before that happens rather than extracting a second time.
            raw_accs_for_offset = (a_accs.copy(), n_accs.copy(), v_accs.copy())
        if acc_bw is not None:
            a_accs, m, sd = normalize_accessibility(a_accs, mean=args.acc_mean,
                                                    std=args.acc_std)
            n_accs = normalize_accessibility(n_accs, mean=m, std=sd)[0]
            v_accs = normalize_accessibility(v_accs, mean=m, std=sd)[0]
            acc_stats[cell] = {"acc_mean": float(m), "acc_std": float(sd)}
            print(f"  {tag}accessibility normalization: mean={m:.4f}, std={sd:.4f}")

        for k, v in (("tr_seqs", a_seqs), ("tr_sigs", a_sigs), ("tr_accs", a_accs),
                     ("tr_prof", a_prof), ("neg_seqs", n_seqs), ("neg_sigs", n_sigs),
                     ("neg_accs", n_accs), ("neg_prof", n_prof), ("val_seqs", v_seqs),
                     ("val_sigs", v_sigs), ("val_accs", v_accs), ("val_prof", v_prof)):
            bins[k].append(v)
        per_cell_counts[cell] = {"train_peaks": int(len(a_sigs)),
                                 "train_negatives": int(len(n_sigs)),
                                 "val_peaks": int(len(v_sigs))}

    def _cat(key):
        vals = [v for v in bins[key] if v is not None]
        if not vals:
            return None
        return vals[0] if len(vals) == 1 else np.concatenate(vals, axis=0)

    tr_seqs, tr_sigs, tr_accs, tr_prof = (_cat("tr_seqs"), _cat("tr_sigs"),
                                          _cat("tr_accs"), _cat("tr_prof"))
    neg_seqs, neg_sigs, neg_accs, neg_prof = (_cat("neg_seqs"), _cat("neg_sigs"),
                                              _cat("neg_accs"), _cat("neg_prof"))
    val_seqs, val_sigs, val_accs, val_prof = (_cat("val_seqs"), _cat("val_sigs"),
                                              _cat("val_accs"), _cat("val_prof"))

    if len(PANEL) > 1:
        print("\nPooled across the panel:")
        for cell, c in per_cell_counts.items():
            print(f"  {cell:<12} train={c['train_peaks']:>7,} "
                  f"neg={c['train_negatives']:>7,} val={c['val_peaks']:>7,}")
        print(f"  {'TOTAL':<12} train={len(tr_sigs):>7,} neg={len(neg_sigs):>7,} "
              f"val={len(val_sigs):>7,}")

    # Residual training. Rejected above for a multi-cell panel, so this is the
    # single-cell path only, using the raw arrays copied inside the loop.
    tr_off = neg_off = val_off = None
    if args.count_offset_model is not None:
        print(f"Computing count offsets from {args.count_offset_model} "
              f"(fold {args.fold_key})...")
        tr_off, neg_off, val_off = compute_count_offsets(
            args.count_offset_model, args.fold_key,
            list(raw_accs_for_offset), args.batch_size, args.device)
        print(f"  offsets: train {tr_off.shape} mean {tr_off.mean():.3f}, "
              f"neg {neg_off.shape} mean {neg_off.mean():.3f}, "
              f"val {val_off.shape} mean {val_off.mean():.3f}")
        print("  this run learns the RESIDUAL; final prediction is "
              "model_output + offset")

    if acc_stats:
        stats_path = os.path.join(args.output_dir, "acc_normalization_stats.json")
        with open(stats_path, "w") as f:
            if len(acc_stats) == 1:
                # ONE cell type means one training statistic genuinely applies, so write the
                # flat {"acc_mean":..,"acc_std":..} shape that 4.1 and compute_count_offsets
                # already read. Keyed on the number of cell types, not on the placeholder
                # name: a panel naming a single cell type is still a single-statistic model,
                # and writing the nested shape there would make its own models unreadable by
                # 4.1's default path for no benefit.
                json.dump(next(iter(acc_stats.values())), f)
            else:
                json.dump({"per_cell_type": acc_stats,
                           "note": "Multi-cell-type model: no single training statistic "
                                   "applies. Predict with statistics derived from the "
                                   "TARGET cell type's own accessibility track, via "
                                   "4.1 --acc-mean/--acc-std fed by 4.26."}, f, indent=2)
    else:
        print("No --accessibility-bw given; skipping accessibility normalization")

    # Build datasets
    train_dataset = MultiModalPeakNegativeSampler(
        tr_seqs, tr_accs, tr_sigs,
        neg_seqs, neg_accs, neg_sigs,
        negative_ratio=args.negative_ratio,
        in_window=args.in_window,
        out_window=args.out_window,
        max_jitter=args.max_jitter,
        reverse_complement=True,
        random_state=42,
        mode=args.mode,
        peak_offsets=tr_off,
        neg_offsets=neg_off,
        peak_profiles=tr_prof,
        neg_profiles=neg_prof
    )
    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.batch_size,
        num_workers=0, pin_memory=True
    )

    # Validation tensors — shape depends on mode
    if args.mode == 'multimodal':
        X_valid = torch.cat([torch.from_numpy(val_seqs),
                             torch.from_numpy(val_accs)], dim=1)  # (N, 5, L)
    elif args.mode == 'sequence':
        X_valid = torch.from_numpy(val_seqs)   # (N, 4, L)
    elif args.mode == 'atac':
        X_valid = torch.from_numpy(val_accs)   # (N, 1, L)
    y_valid = torch.from_numpy(val_sigs)

    # Number of accessibility input channels, from the extracted arrays rather than
    # re-parsing the argument, so the model can never disagree with the data.
    n_acc_channels = tr_accs.shape[1]
    if args.mode in ("multimodal", "atac"):
        print(f"Accessibility input channels: {n_acc_channels}")

    # Model
    model = MultiModalBPNet(
        n_filters=args.n_filters,
        n_acc_filters=args.n_acc_filters,
        n_acc_channels=n_acc_channels,
        n_layers=args.n_layers,
        n_outputs=n_outputs,
        mode=args.mode,
        count_loss_weight=args.count_loss_weight,
        profile_loss_weight=args.profile_loss_weight,
        gate_accessibility=args.gate_accessibility,
        overprediction_weight=args.overprediction_weight,
        name=model_prefix,
        verbose=True
    )
    if args.gate_accessibility:
        assert args.mode == "multimodal", "--gate-accessibility needs multimodal mode"
        print("accessibility gating: ON (initialised open)")
    if args.overprediction_weight != 1.0:
        print(f"over-prediction weight: {args.overprediction_weight}")
    if args.profile_loss_weight != 1.0:
        print(f"profile loss weight: {args.profile_loss_weight} "
              f"(profile term is logged UNWEIGHTED, only its gradient is scaled)")
    if use_prof_target and args.profile_loss_weight == 1.0:
        print("WARNING: a separate profile target at --profile-loss-weight 1.0. If the two "
              "targets differ in read depth this silently reweights the counts term; see "
              "_asymmetric_mixture_loss.")
    print(f"Model trimming: {model.trimming} (output window: "
          f"{args.in_window - 2*model.trimming})")
    assert args.in_window - 2 * model.trimming == args.out_window, (
        f"Trimming mismatch: in_window={args.in_window}, "
        f"trimming={model.trimming}, out_window should be {args.out_window}"
    )

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5
    )

    model.fit(
        training_data=train_loader,
        optimizer=optimizer,
        scheduler=scheduler,
        offset_valid=None if val_off is None else torch.from_numpy(val_off),
        X_valid=X_valid,
        y_valid=y_valid,
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        device=args.device,
        early_stopping=args.early_stopping,
        profile_target=use_prof_target,
        y_profile_valid=(torch.from_numpy(val_prof) if use_prof_target else None)
    )

    # Completion marker, written LAST. A preempted job (the `owners` partition is
    # preemptible) leaves a best-so-far checkpoint that is indistinguishable from a
    # finished one, so evaluators must be able to tell them apart. Absence of this file
    # means the run did not finish.
    with open(os.path.join(args.output_dir, "training_complete.json"), "w") as f:
        json.dump({"epochs_requested": args.max_epochs,
                   "early_stopping": args.early_stopping,
                   "model": f"{model_prefix}.torch"}, f, indent=2)
    print(f"Training complete. Model saved to {model_prefix}.torch")


if __name__ == "__main__":
    main()
