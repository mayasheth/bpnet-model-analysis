#!/usr/bin/env python
"""
Overlay H3K27ac and ATAC meta-profiles around candidate elements.

Tests the premise behind the whole window design: that ATAC is centered on the element
while H3K27ac is displaced onto the flanking nucleosomes. If true, the sequence
determinant and the measured signal are spatially separated, which has direct
architectural consequences (the model must propagate information outward from the
element, and its receptive field must cover the displacement).

All tracks are profiled on the SAME elements, stratified by H3K27ac signal, so the
comparison is not confounded by each track picking its own strong elements. Shapes are
normalized per track (divided by each track's own maximum) so the comparison is about
LOCATION, not magnitude.

Also profiles the new 5'-end H3K27ac track alongside the fragment-extended one, to show
how much of the shallow central dip was an artifact of the 250 bp extension.

Usage:
  python 0.6.compare_profiles.py --elements <narrowPeak> --outdir results --figdir figures \
      --track "ATAC:/path/atac.bw" \
      --track "H3K27ac (250bp frag):/path/coverage.bw" \
      --track "H3K27ac (5' ends):/path/h3k27ac_5p_plus.bw,/path/h3k27ac_5p_minus.bw" \
      --stratify-by "H3K27ac (250bp frag)"
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams  # Nature-spec fonts, spines, dpi


SEED = 0
COLS = ["chr", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"]
# ATAC blue, H3K27ac purple-family, matching the modality palette used in the
# model figures (blue = accessibility, red = sequence, purple = both).
COLORS = ["#2166AC", "#762A83", "#B2182B", "#5AAE61", "#E08214"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", required=True)
    p.add_argument("--track", action="append", required=True,
                   help='"Label:path[,path2]" — multiple paths are summed (e.g. 5\' '
                        "plus and minus strands). Repeatable.")
    p.add_argument("--stratify-by", required=True,
                   help="Label of the track whose signal defines the quintiles")
    p.add_argument("--outdir", required=True)
    p.add_argument("--figdir", required=True)
    p.add_argument("--flank", type=int, default=2000)
    p.add_argument("--bin-size", type=int, default=25)
    p.add_argument("--n-sample", type=int, default=30000)
    return p.parse_args()


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    mid = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, mid)
    return df


def extract(df, paths, flank, bin_size):
    """Binned profile matrix, summing across paths (for stranded pairs)."""
    nbin = (2 * flank) // bin_size
    bws = [pyBigWig.open(p) for p in paths]
    sizes = bws[0].chroms()
    rows, keep = [], []
    chroms, centers = df["chr"].to_numpy(), df["center"].to_numpy()
    for i in range(len(df)):
        c, ct = chroms[i], int(centers[i])
        if c not in sizes or ct - flank < 0 or ct + flank > sizes[c]:
            continue
        acc = np.zeros(2 * flank, dtype=np.float64)
        bad = False
        for bw in bws:
            v = bw.values(c, ct - flank, ct + flank, numpy=True)
            if v is None or len(v) != 2 * flank:
                bad = True
                break
            acc += np.nan_to_num(v, nan=0.0)
        if bad:
            continue
        rows.append(acc.reshape(nbin, bin_size).mean(axis=1).astype(np.float32))
        keep.append(i)
    for bw in bws:
        bw.close()
    return np.vstack(rows), np.array(keep)


def main():
    args = parse_args()
    apply_rcparams()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    specs = []
    for t in args.track:
        label, paths = t.split(":", 1)
        specs.append((label, [p for p in paths.split(",") if p]))
    labels = [s[0] for s in specs]
    if args.stratify_by not in labels:
        raise SystemExit(f"--stratify-by '{args.stratify_by}' not among {labels}")

    df = load_elements(args.elements)
    if args.n_sample and args.n_sample < len(df):
        df = df.sample(n=args.n_sample, random_state=SEED).reset_index(drop=True)
    print(f"Elements: {len(df):,}")

    mats, keeps = {}, {}
    for label, paths in specs:
        m, k = extract(df, paths, args.flank, args.bin_size)
        mats[label], keeps[label] = m, k
        print(f"  {label}: {m.shape[0]:,} elements, {len(paths)} track(s)")

    # restrict to elements usable in EVERY track so all profiles use the same set
    common = keeps[labels[0]]
    for lb in labels[1:]:
        common = np.intersect1d(common, keeps[lb])
    for lb in labels:
        idx = np.searchsorted(keeps[lb], common)
        mats[lb] = mats[lb][idx]
    print(f"Common to all tracks: {len(common):,} elements")

    nbin = mats[labels[0]].shape[1]
    pos = np.arange(nbin) * args.bin_size - args.flank + args.bin_size / 2

    strat = mats[args.stratify_by].sum(axis=1)
    quint = pd.qcut(strat.argsort().argsort(), 5, labels=False)
    top = quint == 4
    print(f"Top quintile: {top.sum():,} elements "
          f"(stratified by {args.stratify_by})")

    out, summary = {"position": pos}, []
    for lb in labels:
        prof = mats[lb][top].mean(axis=0)
        out[lb] = prof
        norm = prof / prof.max()
        out[f"{lb} (normalized)"] = norm
        peak_pos = float(pos[prof.argmax()])
        center = float(prof[np.abs(pos) < args.bin_size].mean())
        # mean |offset| of signal mass, a shape-free measure of how flanking a track is
        w = np.clip(prof - prof.min(), 0, None)
        spread = float((np.abs(pos) * w).sum() / w.sum())
        summary.append({
            "track": lb,
            "peak_offset_bp": peak_pos,
            "center_over_peak": round(center / prof.max(), 3),
            "mean_abs_offset_bp": round(spread, 1),
        })
    pd.DataFrame(out).to_csv(
        os.path.join(args.outdir, "profile_comparison.tsv"), sep="\t", index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(args.outdir, "profile_comparison_summary.tsv"),
              sep="\t", index=False)
    print()
    print(sm.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, normalize in zip(axes, (True, False)):
        for color, lb in zip(COLORS, labels):
            prof = mats[lb][top].mean(axis=0)
            y = prof / prof.max() if normalize else prof
            ax.plot(pos, y, color=color, lw=1.5, label=lb)
        ax.axvline(0, color="black", lw=0.7, ls="--")
        ax.set_xlabel("Distance from element center (bp)")
        ax.set_ylabel("Normalized signal (fraction of own maximum)"
                      if normalize else "Mean signal")
        ax.set_title("Shape comparison" if normalize else "Raw magnitudes",
                     fontsize=9)
        if not normalize:
            ax.set_yscale("log")
        ax.legend(frameon=False, fontsize=7)
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("black")
        ax.tick_params(colors="black")
    fig.suptitle("Top H3K27ac quintile: is H3K27ac displaced relative to ATAC?",
                 fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.figdir, f"profile_comparison.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
