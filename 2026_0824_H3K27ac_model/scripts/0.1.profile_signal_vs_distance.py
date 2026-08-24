#!/usr/bin/env python
"""
Choose the H3K27ac counting window empirically.

H3K27ac is deposited on the nucleosomes flanking a regulatory element, not on the
nucleosome-free element itself, so a counts window sized to the element (as used for
p300) systematically misses most of the signal. This measures where the signal
actually sits, so `--out-window` is set from data rather than guessed.

Outputs:
  1. Meta-profile — mean H3K27ac coverage vs distance from element center, stratified
     by element signal decile. Expected shape is bimodal: a dip over the
     nucleosome-free element with flanking shoulders.
  2. Window trade-off — for each candidate half-width w, how far the profile has
     decayed toward its distal plateau at w, against the fraction of elements with a
     *different* candidate element inside w. Widening buys signal and buys neighbour
     contamination; this trade-off picks w.

     (An earlier version reported "fraction of +/-3 kb signal captured". That is
     normalised to an arbitrary outer bound, so it rises near-linearly and never
     saturates - it cannot identify a window. Replaced.)

--flank must be wide enough to reach the distal plateau; the profile is binned during
extraction so memory stays flat in --flank.

Usage:
  python 0.1.profile_signal_vs_distance.py \
      --elements  reference/K562_DNase_candidate_elements.narrowPeak \
      --signal-bw /path/to/ENCSR000AKP_coverage.bw \
      --outdir results --figdir figures
"""

import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 0
NARROWPEAK_COLS = ["chr", "start", "end", "name", "score", "strand",
                   "signalValue", "pValue", "qValue", "summit"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", required=True,
                   help="narrowPeak of candidate elements to center windows on")
    p.add_argument("--signal-bw", required=True, help="H3K27ac coverage BigWig")
    p.add_argument("--outdir", required=True, help="Directory for TSV outputs")
    p.add_argument("--figdir", required=True, help="Directory for figures")
    p.add_argument("--flank", type=int, default=10000,
                   help="Half-width profiled, bp; must reach the plateau (default: 10000)")
    p.add_argument("--bin-size", type=int, default=50,
                   help="Meta-profile bin size in bp (default: 50)")
    p.add_argument("--n-sample", type=int, default=30000,
                   help="Elements to sample (default: 30000; 0 = all)")
    p.add_argument("--half-windows", type=int, nargs="+",
                   default=[250, 500, 750, 1000, 1500, 2000],
                   help="Candidate half-widths for the trade-off curve")
    p.add_argument("--n-layers", type=int, default=8,
                   help="Model n_layers, to report the required in-window (default: 8)")
    return p.parse_args()


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=NARROWPEAK_COLS)
    # summit is an offset from start; in this file it equals width//2, i.e. the midpoint
    midpoint = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, midpoint)
    return df


def extract_binned(df, bw_path, flank, bin_size):
    """Coverage matrix (n_kept, 2*flank/bin_size), binned during extraction."""
    if (2 * flank) % bin_size:
        raise SystemExit(f"error: 2*flank ({2*flank}) must be divisible by "
                         f"--bin-size ({bin_size})")
    nbin = (2 * flank) // bin_size
    bw = pyBigWig.open(bw_path)
    chrom_sizes = bw.chroms()
    rows = []
    for _, row in df.iterrows():
        chrom, center = row["chr"], int(row["center"])
        s, e = center - flank, center + flank
        if chrom not in chrom_sizes or s < 0 or e > chrom_sizes[chrom]:
            continue
        v = bw.values(chrom, s, e, numpy=True)
        if v is None or len(v) != 2 * flank:
            continue
        v = np.nan_to_num(v, nan=0.0).reshape(nbin, bin_size).mean(axis=1)
        rows.append(v.astype(np.float32))
    bw.close()
    if not rows:
        raise SystemExit("error: no elements yielded usable signal windows")
    return np.vstack(rows)


def neighbour_contamination(elements, half_windows):
    """Fraction of elements with another element's center within +/-hw."""
    el = elements.sort_values(["chr", "center"])
    out = {}
    for hw in half_windows:
        n_near = n_tot = 0
        for _, g in el.groupby("chr", sort=False):
            ce = g["center"].to_numpy()
            n_tot += len(ce)
            if len(ce) < 2:
                continue
            d = np.diff(ce)
            near = np.zeros(len(ce), dtype=bool)
            near[:-1] |= d <= hw
            near[1:] |= d <= hw
            n_near += int(near.sum())
        out[hw] = n_near / n_tot
    return out


def style(ax):
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("black")
    ax.tick_params(colors="black")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    all_elements = load_elements(args.elements)
    print(f"Loaded {len(all_elements):,} candidate elements from {args.elements}")

    df = all_elements
    if args.n_sample and args.n_sample < len(df):
        df = df.sample(n=args.n_sample, random_state=SEED)
        print(f"Sampled {len(df):,} elements (seed={SEED})")

    mat = extract_binned(df, args.signal_bw, args.flank, args.bin_size)
    nbin = mat.shape[1]
    positions = np.arange(nbin) * args.bin_size - args.flank + args.bin_size / 2
    print(f"Extracted binned signal for {mat.shape[0]:,} elements "
          f"({nbin} bins of {args.bin_size} bp)")

    total = mat.sum(axis=1)
    decile = pd.qcut(total.argsort().argsort(), 10, labels=False)  # rank-based, tie-safe

    # --- 1. Meta-profile ----------------------------------------------------
    profile = pd.DataFrame({"position": positions,
                            "mean_all": mat.mean(axis=0),
                            "median_all": np.median(mat, axis=0)})
    for d in (3, 5, 7, 9):
        profile[f"mean_decile{d + 1}"] = mat[decile == d].mean(axis=0)
    profile_path = os.path.join(args.outdir, "signal_vs_distance.tsv")
    profile.to_csv(profile_path, sep="\t", index=False)
    print(f"Wrote {profile_path}")

    # --- 2. Window trade-off ------------------------------------------------
    top = mat[decile == 9].mean(axis=0)
    plateau = float(top[np.abs(positions) > 0.8 * args.flank].mean())
    peak = float(top.max())
    peak_pos = float(positions[top.argmax()])
    contam = neighbour_contamination(all_elements, args.half_windows)
    trimming = 47 + sum(2 ** i for i in range(1, args.n_layers + 1))

    rows = []
    for hw in args.half_windows:
        edge = float((top[np.abs(positions - hw).argmin()] +
                      top[np.abs(positions + hw).argmin()]) / 2)
        rows.append({
            "half_window": hw,
            "out_window": 2 * hw,
            "profile_at_edge": round(edge, 2),
            # 1.0 = still at the peak, 0.0 = fully decayed to the distal plateau
            "frac_above_plateau": round((edge - plateau) / (peak - plateau), 3),
            "frac_elements_contaminated": round(contam[hw], 3),
            "required_in_window": 2 * hw + 2 * trimming,
        })
    trade = pd.DataFrame(rows)
    trade_path = os.path.join(args.outdir, "window_tradeoff.tsv")
    trade.to_csv(trade_path, sep="\t", index=False)
    print(f"Wrote {trade_path}")
    print(f"peak {peak:.2f} at {peak_pos:+.0f} bp; distal plateau {plateau:.2f} "
          f"({peak / plateau:.2f}x enrichment); trimming={trimming}")
    print(trade.to_string(index=False))

    # --- Figure -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))

    ax = axes[0]
    shades = plt.get_cmap("PuBu")(np.linspace(0.45, 0.95, 4))
    for color, d in zip(shades, (3, 5, 7, 9)):
        ax.plot(positions, mat[decile == d].mean(axis=0), color=color, lw=1.4,
                label=f"Decile {d + 1}")
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.axhline(plateau, color="black", lw=0.7, ls=":")
    ax.set_xlim(-3000, 3000)
    ax.set_xlabel("Distance from element center (bp)")
    ax.set_ylabel("Mean H3K27ac coverage")
    ax.set_title("Signal sits on the flanking nucleosomes", fontsize=9)
    ax.legend(frameon=False, fontsize=7, title="Signal decile", title_fontsize=7)
    style(ax)

    ax = axes[1]
    ax.plot(trade["half_window"], trade["frac_above_plateau"], marker="o", ms=4,
            color="#0096a0", label="Signal remaining at window edge")
    ax.plot(trade["half_window"], trade["frac_elements_contaminated"], marker="s",
            ms=4, color="#c5373d", label="Elements with a neighbour inside")
    ax.set_xlabel("Window half-width (bp)")
    ax.set_ylabel("Fraction")
    ax.set_title("Wider window: more signal, more neighbours", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    style(ax)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(args.figdir, f"h3k27ac_signal_vs_distance.{ext}")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
