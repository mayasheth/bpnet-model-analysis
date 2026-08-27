#!/usr/bin/env python
"""
How much of H3K27ac does accessibility explain, with no model involved?

The ATAC-only model reaches 0.551 on K562 top-quintile elements, but that mixes what
accessibility can explain in principle with how well a network learns it. This measures
the former directly: the correlation between summed ATAC and summed H3K27ac over the same
element windows. It is the model-free reference for the ATAC-only baseline, and it needs
no GPU.

Run per cell type and the results accumulate in one TSV, so the panel can be compared.
If accessibility explains less of H3K27ac in other cell types, the ATAC-only model's
advantage in K562 is partly a property of K562 rather than of the assay pair.

Both tracks must be processed the same way across cell types or the comparison measures
processing rather than biology.

Usage:
  python 0.12.atac_vs_h3k27ac.py --label K562 \
      --elements <narrowPeak> --atac-bw <atac.bw> \
      --h3k27ac-bw <5p_plus.bw>,<5p_minus.bw> \
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

COLS = ["chr", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"]
SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--label", required=True, help="Cell type / sample label")
    p.add_argument("--elements", required=True)
    p.add_argument("--atac-bw", required=True,
                   help="Accessibility BigWig; comma-separate to sum tracks")
    p.add_argument("--h3k27ac-bw", required=True,
                   help="H3K27ac BigWig; comma-separate to sum (e.g. 5' plus,minus)")
    p.add_argument("--half-window", type=int, default=500)
    p.add_argument("--n-sample", type=int, default=0, help="0 = all elements")
    p.add_argument("--outdir", default="results")
    p.add_argument("--figdir", default="figures")
    return p.parse_args()


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    mid = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, mid)
    return df


def window_sums(df, spec, hw):
    """Summed signal in +/-hw around each element centre; NaN if out of bounds."""
    paths = [p for p in spec.split(",") if p]
    bws = [pyBigWig.open(p) for p in paths]
    sizes = bws[0].chroms()
    out = np.full(len(df), np.nan)
    chroms, centers = df["chr"].to_numpy(), df["center"].to_numpy()
    for i in range(len(df)):
        c, ct = chroms[i], int(centers[i])
        if c not in sizes or ct - hw < 0 or ct + hw > sizes[c]:
            continue
        tot, bad = 0.0, False
        for bw in bws:
            v = bw.values(c, ct - hw, ct + hw, numpy=True)
            if v is None or len(v) != 2 * hw:
                bad = True
                break
            tot += float(np.nan_to_num(v, nan=0.0).sum())
        if not bad:
            out[i] = tot
    for bw in bws:
        bw.close()
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    df = load_elements(args.elements)
    if args.n_sample and args.n_sample < len(df):
        df = df.sample(n=args.n_sample, random_state=SEED).reset_index(drop=True)
    print(f"{args.label}: {len(df):,} elements, +/-{args.half_window} bp windows")

    atac = window_sums(df, args.atac_bw, args.half_window)
    k27 = window_sums(df, args.h3k27ac_bw, args.half_window)
    ok = np.isfinite(atac) & np.isfinite(k27)
    x, y = np.log1p(atac[ok]), np.log1p(k27[ok])
    print(f"  usable: {ok.sum():,}")

    # stratify by H3K27ac, matching how every other evaluation in this project reports
    q = pd.qcut(y.argsort().argsort(), 5, labels=False)
    rows = []
    for name, mask in [("all", np.ones(len(y), bool)),
                       ("top_quintile", q == 4),
                       ("top_two_quintiles", q >= 3)]:
        if mask.sum() < 10:
            continue
        rows.append({"label": args.label, "stratum": name, "n": int(mask.sum()),
                     "pearson": round(float(pearsonr(x[mask], y[mask])[0]), 4),
                     "spearman": round(float(spearmanr(x[mask], y[mask])[0]), 4),
                     "mean_log_atac": round(float(x[mask].mean()), 3),
                     "mean_log_h3k27ac": round(float(y[mask].mean()), 3),
                     "half_window": args.half_window})
    new = pd.DataFrame(rows)
    path = os.path.join(args.outdir, "atac_vs_h3k27ac_by_celltype.tsv")
    if os.path.exists(path):
        old = pd.read_csv(path, sep="\t")
        old = old[old["label"] != args.label]          # replace this label's rows
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(path, sep="\t", index=False)
    print(f"Wrote {path}")
    print(new.to_string(index=False))

    # hexbin scatter for this cell type
    fig, ax = plt.subplots(figsize=(3.5, 3.1))
    hb = ax.hexbin(x, y, gridsize=55, bins="log", cmap="Blues", mincnt=1, linewidths=0)
    r = pearsonr(x, y)[0]
    ax.set_xlabel("log1p ATAC counts")
    ax.set_ylabel("log1p H3K27ac counts")
    ax.set_title(f"{args.label}: r = {r:.3f} (all elements)", fontsize=9)
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.colorbar(hb, ax=ax, label="elements", shrink=0.85)
    fig.tight_layout()
    safe = args.label.replace(" ", "_").replace("+", "plus")
    for ext in ("pdf", "png"):
        fp = os.path.join(args.figdir, f"atac_vs_h3k27ac_{safe}.{ext}")
        fig.savefig(fp, dpi=300, bbox_inches="tight")
        print(f"Wrote {fp}")


if __name__ == "__main__":
    main()
