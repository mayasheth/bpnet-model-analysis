#!/usr/bin/env python
"""
Inter-replicate H3K27ac ceiling as a function of counting window.

Two replicates measure the same biology, so their agreement bounds how well any model
can predict counts. Because a wider window sums more reads, Poisson noise falls and the
ceiling rises with window size — which is the other half of the window trade-off in
`0.1.profile_signal_vs_distance.py` (that script measures signal captured and neighbour
contamination; this one measures how much signal is even reproducible).

Reports Pearson and Spearman on log1p counts per window, over all elements and over
the top signal quintile (the elements that actually carry H3K27ac).

Requires per-replicate BigWigs from 0.2.make_replicate_bigwigs.sh.

Usage:
  python 0.3.replicate_ceiling_by_window.py \
      --elements reference/K562_DNase_candidate_elements.narrowPeak \
      --rep1-bw data/h3k27ac_rep1.bw --rep2-bw data/h3k27ac_rep2.bw \
      --outdir results --figdir figures
"""

import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 0
NARROWPEAK_COLS = ["chr", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", required=True)
    p.add_argument("--rep1-bw", required=True)
    p.add_argument("--rep2-bw", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--figdir", required=True)
    p.add_argument("--half-windows", type=int, nargs="+",
                   default=[250, 500, 750, 1000, 1500, 2000])
    p.add_argument("--n-sample", type=int, default=0,
                   help="Elements to sample; 0 = all (default)")
    return p.parse_args()


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=NARROWPEAK_COLS)
    midpoint = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, midpoint)
    return df


def window_sums(df, bw_path, half_windows):
    """Per-element summed coverage for each half-window. NaN where out of bounds."""
    bw = pyBigWig.open(bw_path)
    sizes = bw.chroms()
    max_hw = max(half_windows)
    out = {hw: np.full(len(df), np.nan) for hw in half_windows}
    centers = df["center"].to_numpy()
    chroms = df["chr"].to_numpy()
    for i in range(len(df)):
        chrom, center = chroms[i], int(centers[i])
        if chrom not in sizes or center - max_hw < 0 or center + max_hw > sizes[chrom]:
            continue
        v = bw.values(chrom, center - max_hw, center + max_hw, numpy=True)
        if v is None or len(v) != 2 * max_hw:
            continue
        v = np.nan_to_num(v, nan=0.0)
        for hw in half_windows:
            out[hw][i] = v[max_hw - hw:max_hw + hw].sum()
    bw.close()
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    df = load_elements(args.elements)
    if args.n_sample and args.n_sample < len(df):
        df = df.sample(n=args.n_sample, random_state=SEED).reset_index(drop=True)
    print(f"Elements: {len(df):,}")

    print("Summing rep1...")
    s1 = window_sums(df, args.rep1_bw, args.half_windows)
    print("Summing rep2...")
    s2 = window_sums(df, args.rep2_bw, args.half_windows)

    rows = []
    for hw in args.half_windows:
        a, b = s1[hw], s2[hw]
        ok = np.isfinite(a) & np.isfinite(b)
        x, y = np.log1p(a[ok]), np.log1p(b[ok])
        # top quintile by mean of the two replicates — the elements carrying signal
        mean_sig = (x + y) / 2
        thr = np.quantile(mean_sig, 0.8)
        top = mean_sig >= thr
        rows.append({
            "half_window": hw,
            "out_window": 2 * hw,
            "n_elements": int(ok.sum()),
            "pearson_all": round(float(pearsonr(x, y)[0]), 4),
            "spearman_all": round(float(spearmanr(x, y)[0]), 4),
            "pearson_top_quintile": round(float(pearsonr(x[top], y[top])[0]), 4),
            "spearman_top_quintile": round(float(spearmanr(x[top], y[top])[0]), 4),
            "mean_counts_rep1": round(float(np.mean(a[ok])), 1),
        })
    res = pd.DataFrame(rows)
    path = os.path.join(args.outdir, "replicate_ceiling_by_window.tsv")
    res.to_csv(path, sep="\t", index=False)
    print(f"Wrote {path}")
    print(res.to_string(index=False))

    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.plot(res["half_window"], res["pearson_all"], marker="o", ms=4,
            color="#0096a0", label="All elements")
    ax.plot(res["half_window"], res["pearson_top_quintile"], marker="s", ms=4,
            color="#792374", label="Top signal quintile")
    ax.set_xlabel("Window half-width (bp)")
    ax.set_ylabel("Inter-replicate Pearson r (log1p counts)")
    ax.set_title("Achievable ceiling rises with window", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("black")
    ax.tick_params(colors="black")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.figdir, f"h3k27ac_replicate_ceiling.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
