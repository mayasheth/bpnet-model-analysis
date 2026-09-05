#!/usr/bin/env python
"""Predict H3K27ac over ABC's candidate regions and write it as a bigWig ABC can consume.

WHY A BIGWIG. `neighborhoods.py:count_bigwig` sums bigwig values over each candidate region
and `count_bigwig_total` sums genome-wide for the RPM denominator, so a painted bigwig is
counted exactly the way a real track is. No fork of ABC is needed -- the file goes in the
`H3K27ac` column of the biosample table.

WHY THE SCALE DOES NOT MATTER. `run_qnorm` is RANK-BASED: it maps each region's
within-sample quantile onto the K562 reference distribution, so `activity_base` depends only
on the rank order of the injected values. That is what makes it legitimate to feed a
model trained on 5' end counts in a +/-500 bp window into a pipeline built around read
counts over the element, and it is why use_qnorm must stay True.

WHY ONE MODEL PER CHROMOSOME. ABC scores the whole genome, but each fold model trained on
four fifths of it. Every chromosome is held out (val or test) in at least one fold, so the
genome-wide track is assembled from the fold that never saw each chromosome. Without this
the CRISPR benchmark is contaminated and will look good for the wrong reason. The same rule
is applied to models from a DIFFERENT cell type, where leakage is not a concern, so that a
cross-cell-type arm cannot gain an advantage from ensembling that the same-cell-type arm is
forbidden.

PAINTING. Each region is painted with predicted_counts / width, so ABC's sum over the region
recovers the model's prediction exactly. The alternative -- painting the prediction itself,
so the sum scales with region width as real read counts do -- is the `--paint density` flag.
Regions the model cannot score (window off the chromosome end) are left empty and reported;
they read as 0 and therefore rank lowest, so the count needs to stay negligible.

Usage:
  4.1.predict_h3k27ac_for_abc.py --regions REG.bed --model-dir DIR --mode MODE \
      --accessibility-bw ATAC.bw --chrom-sizes SIZES --out-bw OUT.bw
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import pyBigWig
import torch

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, normalize_accessibility

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
OUT_W = 1000

ap = argparse.ArgumentParser()
ap.add_argument("--regions", required=True, help="ABC candidateRegions.bed (chr start end)")
ap.add_argument("--model-dir", required=True, help="contains fold0..fold4")
ap.add_argument("--mode", required=True, choices=["sequence", "atac", "multimodal"])
ap.add_argument("--accessibility-bw", default=None)
ap.add_argument("--chrom-sizes", required=True, help="ABC's chrom sizes, for the output bigwig")
ap.add_argument("--out-bw", required=True)
ap.add_argument("--signal-bw", default=None,
                help="Any valid bigwig on this genome. extract_windows requires a signal "
                     "track and uses it for the bounds check; its VALUES ARE DISCARDED "
                     "here. Defaults to the K562 H3K27ac 5-prime plus track so the kept-"
                     "region set matches what training and evaluation saw.")
ap.add_argument("--chroms", default=None,
                help="Comma-separated subset, for a fast end-to-end test (e.g. chr22).")
ap.add_argument("--paint", default="counts", choices=["counts", "density"],
                help="counts: sum over region equals the prediction (default). "
                     "density: sum scales with region width, as real read counts do.")
ap.add_argument("--batch", type=int, default=256)
ap.add_argument("--no-rc-average", dest="rc_average", action="store_false", default=True,
                help="Disable test-time reverse-complement averaging. RC averaging is ON by "
                     "default as of 2026-09-05, matching the 2.15 evaluator: the number we "
                     "report and the track we ship should be produced the same way, or the "
                     "quoted metric is not the pipeline's actual output. Measured gain on "
                     "top-quintile counts: +0.0162 sequence, +0.0081 multimodal, +0.0016 "
                     "and non-significant for ATAC only.")
a = ap.parse_args()

SIGNAL_DEFAULT = ("/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
                  "2026_0824_H3K27ac_model/data/h3k27ac_5p_plus.bw")
if a.signal_bw is None:
    a.signal_bw = SIGNAL_DEFAULT
if not os.path.exists(a.signal_bw):
    raise SystemExit(f"error: signal bigwig {a.signal_bw} not found")
if a.mode in ("multimodal", "atac") and not a.accessibility_bw:
    raise SystemExit(f"error: mode {a.mode} needs --accessibility-bw")

dev = "cuda" if torch.cuda.is_available() else "cpu"
folds = json.load(open(FOLDS))

# chromosome -> the lowest-numbered fold that did not train on it
chrom_fold = {}
for k in sorted(folds):
    for c in folds[k]["val"] + folds[k]["test"]:
        chrom_fold.setdefault(c, k)

sizes = {}
for line in open(a.chrom_sizes):
    p = line.split()
    if len(p) >= 2:
        sizes[p[0]] = int(p[1])

reg = pd.read_csv(a.regions, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["chr", "start", "end"])
reg = reg[reg["chr"].isin(sizes)].reset_index(drop=True)
if a.chroms:
    keep = set(a.chroms.split(","))
    reg = reg[reg["chr"].isin(keep)].reset_index(drop=True)
    print(f"restricted to {sorted(keep)}")
print(f"{len(reg):,} regions on {reg['chr'].nunique()} chromosomes", flush=True)
unmapped = sorted(set(reg["chr"]) - set(chrom_fold))
if unmapped:
    raise SystemExit(f"error: no held-out fold for {unmapped}; cannot predict without leakage")

# One bigwig interval per region, in sorted order, so pyBigWig accepts them.
records = {}
n_pred = n_drop = 0
for chrom, g in reg.groupby("chr", sort=True):
    fold = chrom_fold[chrom]
    md = f"{a.model_dir}/fold{fold}"
    if not os.path.exists(f"{md}/training_complete.json"):
        raise SystemExit(f"error: {md} has no training_complete.json")
    m = torch.load(f"{md}/multimodal_bpnet.torch", map_location="cpu", weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = a.mode
    in_w = OUT_W + 2 * m.trimming
    m = m.to(dev).eval()

    df = g.copy().reset_index(drop=True)
    df["summit"] = ((df["end"] - df["start"]) // 2).astype(int)
    seqs, _, accs, valid = extract_windows(
        df, GEN if a.mode != "atac" else None,
        a.signal_bw, None,
        a.accessibility_bw, in_w, OUT_W, 0, is_peak=True)
    kept = df[valid].reset_index(drop=True)
    n_drop += int((~valid).sum())

    x = accs
    if a.mode in ("multimodal", "atac"):
        st = json.load(open(f"{md}/acc_normalization_stats.json"))
        x = normalize_accessibility(accs, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, x], axis=1) if a.mode == "multimodal"
         else seqs if a.mode == "sequence" else x).astype(np.float32)

    # Reverse complement: one-hot channels are ACGT, so flipping the channel axis maps
    # A<->T and C<->G, and flipping the length axis completes it. Accessibility channels are
    # strand-agnostic coverage and are reversed along length only. Only the counts head is
    # used here, so unlike 2.15 there is no profile to un-flip.
    n_seq = 4 if a.mode in ("multimodal", "sequence") else 0

    def rc(xb):
        if n_seq:
            seq_part = torch.flip(xb[:, :n_seq], dims=[1, 2])
            if xb.shape[1] > n_seq:
                return torch.cat([seq_part, torch.flip(xb[:, n_seq:], dims=[2])], dim=1)
            return seq_part
        return torch.flip(xb, dims=[2])

    out = []
    with torch.no_grad():
        for i in range(0, len(X), a.batch):
            xb = torch.from_numpy(X[i:i + a.batch]).to(dev)
            _, lc = m(xb)
            if a.rc_average:
                _, lc_rc = m(rc(xb))
                lc = (lc + lc_rc) / 2
            out.append(lc.squeeze(-1).cpu().numpy())
    m.to("cpu"); del X, seqs, accs, x
    lc = np.concatenate(out) if out else np.zeros(0)

    counts = np.expm1(lc).clip(min=0.0)          # model emits log1p(counts)
    width = (kept["end"] - kept["start"]).to_numpy(float)
    vals = counts / width if a.paint == "counts" else counts / OUT_W
    records[chrom] = (kept["start"].to_numpy(int), kept["end"].to_numpy(int),
                      vals.astype(np.float32))
    n_pred += len(kept)
    print(f"  {chrom}: fold{fold} in_window={in_w} n={len(kept):,} "
          f"dropped={int((~valid).sum())} median_counts={np.median(counts):.1f} "
          f"rc={'on' if a.rc_average else 'off'}", flush=True)

print(f"\npredicted {n_pred:,} regions; {n_drop:,} dropped "
      f"({100.0*n_drop/max(1,n_pred+n_drop):.3f}% -- these read as 0 in ABC)", flush=True)

bw = pyBigWig.open(a.out_bw, "w")
bw.addHeader([(c, sizes[c]) for c in sorted(records)])
for c in sorted(records):
    s, e, v = records[c]
    o = np.argsort(s)
    bw.addEntries([c] * len(s), s[o].tolist(), ends=e[o].tolist(), values=v[o].tolist())
bw.close()
print("wrote", a.out_bw)

chk = pyBigWig.open(a.out_bw)
tot = sum(l * (chk.stats(c, 0, l, "mean", exact=True)[0] or 0) for c, l in chk.chroms().items())
print(f"genome-wide sum (ABC's RPM denominator) = {tot:,.1f}")
assert abs(tot) > 0, "empty bigwig; ABC asserts a nonzero total"
chk.close()
