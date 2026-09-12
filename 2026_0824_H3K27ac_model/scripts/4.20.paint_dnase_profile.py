#!/usr/bin/env python
"""Paint converted DNase GENOME-WIDE at base resolution: profile x counts.

WHY 4.1 WILL NOT DO. 4.1 paints `predicted_counts / width`, constant across each region,
because ABC only ever sums the track over a region and the profile is discarded. Feeding
that to the downstream model's accessibility branch would hand it exactly the flat track
the converter exists to avoid. This writes `softmax(profile) * exp(counts)` per base, which
is the "profile x counts at base resolution" 0.33 specified as the converter's deliverable.

WHY GENOME-WIDE AND NOT JUST AT CANDIDATE ELEMENTS. Three reasons, and each one alone would
force it:
  * the downstream model reads a 2,114 bp accessibility window per element, wider than this
    model's 1,000 bp output window, so element-only painting leaves the flanks empty;
  * training draws 50,000 genome-wide GC-matched negatives, and the ATAC and real-DNase arms
    see real accessibility there. A converted track that is zero off-element would give the
    negatives no signal at all, which is a confound in the arm's favour or against it
    depending on the model's use of negatives, not a detail;
  * a track that only exists where someone already called peaks is not deployable, which was
    the point of building a converter.

WINDOWS OVERLAP AND ARE AVERAGED, to remove seams. At stride == OUT_W the windows abut and
each base is predicted exactly once, which leaves a hard discontinuity every 1,000 bp; the
downstream accessibility branch has a ~1.1 kb receptive field and reads across those
constantly. At stride 250 each base is predicted from four offsets and the mean is written.

DO NOT EXPECT THIS TO FIX THE COUNTS. It does not, and an earlier version of this comment
claimed otherwise. Both strides were trialled on chr22: per-element log1p count correlation
against observed DNase went 0.682 -> 0.687 and the painted/observed total ratio 0.767 ->
0.770. Essentially unchanged, so the seams were NOT the cause of the count shortfall.

THE ACTUAL CAUSE of that shortfall is that these windows sit at fixed genomic offsets
(multiples of the stride from position 0) and never coincide with an element's own centre,
whereas 2.15 predicts with the window centred exactly on the element and reports
overall_pearson 0.835. Summing a tiled track over an element therefore blends several
off-centre predictions. That is inherent to producing a genome-wide track rather than a
per-element number, not a bug to be fixed, but it does mean the painted track's per-element
magnitudes are weaker than 2.15's headline count metric suggests. The profile is unaffected:
1 bp top-quintile shape against real DNase is 0.781 here versus 2.31's fold-2 value of 0.770.

The ~77% scale factor is absorbed downstream, because training computes its own
accessibility normalisation from whatever track it is given. A SIGNAL-DEPENDENT shortfall
would not be absorbed and would compress dynamic range, which is the diagnosed failure mode
of predicted H3K27ac -- worth checking before reading any negative result as "the profile
does not help".

LEAKAGE. Each chromosome is painted by the lowest-numbered fold that did not train on it,
the same rule 4.1 uses. Without it the downstream benchmark is contaminated and would look
good for the wrong reason.

MEMORY. A whole chromosome of windows does not fit: chr1 at stride 1,000 is ~248,000 windows
and the sequence array alone would be ~8 GB. Windows are therefore extracted and predicted
in chunks, accumulating into one float32 array per chromosome (chr1 ~1 GB).

Usage:
  4.20.paint_dnase_profile.py --model-dir DIR --mode MODE --accessibility-bw ATAC.bw \
      --chrom-sizes SIZES --out-bw OUT.bw [--stride 1000] [--chroms chr22]
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
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
OUT_W = 1000

ap = argparse.ArgumentParser()
ap.add_argument("--model-dir", required=True)
ap.add_argument("--mode", default="multimodal", choices=["sequence", "atac", "multimodal"])
ap.add_argument("--accessibility-bw", default=None)
ap.add_argument("--chrom-sizes", required=True)
ap.add_argument("--out-bw", required=True)
ap.add_argument("--signal-bw", default=f"{P}/data/k562_dnase_5p_plus.bw",
                help="extract_windows requires a signal track for its bounds check; the "
                     "VALUES ARE DISCARDED. Defaults to the K562 DNase plus track.")
ap.add_argument("--stride", type=int, default=250,
                help="Window spacing. Default 250 against a 1,000 bp output window, so "
                     "every base is predicted from four offsets and the mean is written. "
                     "Set equal to OUT_W for abutting windows, which is faster and leaves "
                     "seams; see the docstring for what that cost.")
ap.add_argument("--chunk", type=int, default=8000, help="windows extracted per chunk")
ap.add_argument("--batch", type=int, default=256)
ap.add_argument("--chroms", default=None, help="comma-separated subset, for a smoke test")
ap.add_argument("--no-rc-average", dest="rc_average", action="store_false", default=True)
a = ap.parse_args()

if a.stride > OUT_W:
    raise SystemExit(f"stride {a.stride} exceeds the {OUT_W} bp output window; that would "
                     f"leave unpainted gaps between windows")
if OUT_W % a.stride:
    raise SystemExit(f"stride {a.stride} does not divide the {OUT_W} bp output window, so "
                     f"coverage would be uneven across positions")
if a.mode in ("multimodal", "atac") and not a.accessibility_bw:
    raise SystemExit(f"error: mode {a.mode} needs --accessibility-bw")

dev = "cuda" if torch.cuda.is_available() else "cpu"
folds = json.load(open(FOLDS))
chrom_fold = {}
for k in sorted(folds):
    for c in folds[k]["val"] + folds[k]["test"]:
        chrom_fold.setdefault(c, k)

sizes = {}
for line in open(a.chrom_sizes):
    p = line.split()
    if len(p) >= 2:
        sizes[p[0]] = int(p[1])
chroms = sorted(c for c in sizes if c in chrom_fold)
if a.chroms:
    keep = set(a.chroms.split(","))
    chroms = [c for c in chroms if c in keep]
    print(f"restricted to {chroms}")
if not chroms:
    raise SystemExit("no chromosomes with a held-out fold to paint")

bw = pyBigWig.open(a.out_bw, "w")
bw.addHeader([(c, sizes[c]) for c in chroms])

total_pred = total_drop = 0
grand_sum = 0.0
for chrom in chroms:
    fold = chrom_fold[chrom]
    md = f"{a.model_dir}/fold{fold}"
    if not os.path.exists(f"{md}/training_complete.json"):
        raise SystemExit(f"error: {md} has no training_complete.json")
    m = torch.load(f"{md}/multimodal_bpnet.torch", map_location="cpu", weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = a.mode
    in_w = OUT_W + 2 * m.trimming
    m = m.to(dev).eval()
    st = (json.load(open(f"{md}/acc_normalization_stats.json"))
          if a.mode in ("multimodal", "atac") else None)
    n_seq = 4 if a.mode in ("multimodal", "sequence") else 0

    def rc_input(xb):
        if n_seq:
            seq_part = torch.flip(xb[:, :n_seq], dims=[1, 2])
            if xb.shape[1] > n_seq:
                return torch.cat([seq_part, torch.flip(xb[:, n_seq:], dims=[2])], dim=1)
            return seq_part
        return torch.flip(xb, dims=[2])

    # Two accumulators: summed prediction and how many windows covered each base. The
    # divisor is counted rather than assumed to be OUT_W/stride, because coverage tapers at
    # the chromosome ends and wherever a window was dropped by the bounds check.
    track = np.zeros(sizes[chrom], dtype=np.float32)
    cover = np.zeros(sizes[chrom], dtype=np.uint8)
    starts = np.arange(0, sizes[chrom] - OUT_W + 1, a.stride, dtype=np.int64)
    n_c = n_d = 0
    for i0 in range(0, len(starts), a.chunk):
        s = starts[i0:i0 + a.chunk]
        df = pd.DataFrame({"chr": chrom, "start": s, "end": s + OUT_W,
                           "summit": OUT_W // 2})
        seqs, _, accs, valid = extract_windows(
            df, GEN if a.mode != "atac" else None, a.signal_bw, None,
            a.accessibility_bw, in_w, OUT_W, 0, is_peak=True)
        n_d += int((~valid).sum())
        if valid.sum() == 0:
            continue
        kept = s[valid]
        x = accs
        if st is not None:
            x = normalize_accessibility(accs, mean=st["acc_mean"], std=st["acc_std"])[0]
        X = (np.concatenate([seqs, x], axis=1) if a.mode == "multimodal"
             else seqs if a.mode == "sequence" else x).astype(np.float32)
        del seqs, accs, x

        with torch.no_grad():
            for j in range(0, len(X), a.batch):
                xb = torch.from_numpy(X[j:j + a.batch]).to(dev)
                pr, lc = m(xb)
                sh = pr.shape
                p = torch.softmax(pr.reshape(sh[0], -1), dim=-1)
                if a.rc_average:
                    pr_rc, lc_rc = m(rc_input(xb))
                    lc = (lc + lc_rc) / 2
                    # Undo the RC transform before averaging: reverse positions AND swap
                    # strands. Averaging without the un-flip silently blends a profile with
                    # its mirror image, which looks like a smoother, better-behaved track.
                    p_rc = torch.softmax(torch.flip(pr_rc, dims=[1, 2]).reshape(sh[0], -1),
                                         dim=-1)
                    p = (p + p_rc) / 2
                # Sum the two strands: the downstream accessibility branch takes a single
                # unstranded channel, and the profile head's multinomial is defined over
                # (strand, position) jointly, so the per-base probability is their sum.
                p = p.reshape(sh).sum(dim=1)
                # Counts head is trained on log1p of the window total, so invert it.
                tot = torch.expm1(lc.squeeze(-1)).clamp_min(0.0)
                vals = (p * tot.unsqueeze(-1)).cpu().numpy()
                for k, st0 in enumerate(kept[j:j + a.batch]):
                    track[st0:st0 + OUT_W] += vals[k]
                    cover[st0:st0 + OUT_W] += 1
                n_c += len(vals)
        del X

    m.to("cpu")
    np.divide(track, np.maximum(cover, 1), out=track)
    uncovered = int((cover == 0).sum())
    nz = track > 0
    if nz.any():
        # Runs of equal value: a 1,000 bp window is one interval only where the profile is
        # flat, so write per-base intervals compressed by run-length.
        idx = np.flatnonzero(nz)
        v = track[idx]
        brk = np.flatnonzero(np.diff(idx) != 1) + 1
        seg_bounds = np.concatenate([[0], brk, [len(idx)]])
        cs, ce, cv = [], [], []
        for b0, b1 in zip(seg_bounds[:-1], seg_bounds[1:]):
            ii, vv = idx[b0:b1], v[b0:b1]
            ch = np.flatnonzero(np.diff(vv) != 0) + 1
            rb = np.concatenate([[0], ch, [len(vv)]])
            for r0, r1 in zip(rb[:-1], rb[1:]):
                cs.append(int(ii[r0])); ce.append(int(ii[r1 - 1]) + 1); cv.append(float(vv[r0]))
        bw.addEntries([chrom] * len(cs), cs, ends=ce, values=cv)
    grand_sum += float(track.sum())
    total_pred += n_c; total_drop += n_d
    print(f"  {chrom}: fold{fold} windows={n_c:,} dropped={n_d} uncovered_bp={uncovered:,} "
          f"sum={track.sum():,.0f} rc={'on' if a.rc_average else 'off'}", flush=True)
    del track, cover

bw.close()
print(f"\npainted {total_pred:,} windows; {total_drop:,} dropped at chromosome edges")
print(f"genome-wide sum = {grand_sum:,.1f}")
print(f"wrote {a.out_bw}")
