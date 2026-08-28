#!/usr/bin/env python
"""
ATAC fragment-length distribution, with the channel boundaries marked.

Justifies (or refutes) the fragment-size-stratified channel split. A usable ATAC library
shows two populations: sub-nucleosomal fragments from open chromatin, and a
mono-nucleosomal bump from Tn5 cutting either side of a single nucleosome. If the
mono-nucleosomal population is absent, the "mono" channel carries no nucleosome
positioning information and the whole idea is moot.

Reads TLEN from properly-paired reads in the paired-end BAMs (the tagAlign files are
per-read and lose fragment length entirely, which is why the BAMs had to be re-downloaded).

Usage:
  python 0.10.plot_fragment_lengths.py --bams a.bam b.bam --outdir results --figdir figures
"""

import argparse
import os
import sys
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams  # Nature-spec fonts, spines, dpi


SAMTOOLS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/samtools"
SUB = (1, 99)         # channel: sub-nucleosomal
MONO = (180, 247)     # channel: mono-nucleosomal


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bams", nargs="+", required=True)
    p.add_argument("--region", default="chr1")
    p.add_argument("--n-reads", type=int, default=3000000)
    p.add_argument("--max-len", type=int, default=800)
    p.add_argument("--outdir", required=True)
    p.add_argument("--figdir", required=True)
    return p.parse_args()


def tlens(bam, region, n_reads, max_len):
    """Positive TLEN values from properly-paired reads (each pair counted once)."""
    cmd = f"{SAMTOOLS} view -f 0x2 {bam} {region} | head -{n_reads} | cut -f9"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    v = np.fromiter((int(x) for x in out.split() if x), dtype=np.int32)
    return v[(v > 0) & (v <= max_len)]


def main():
    args = parse_args()
    apply_rcparams()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    rows = []
    colors = ["#0096a0", "#792374", "#e96a00"]
    bins = np.arange(0, args.max_len + 5, 5)
    for color, bam in zip(colors, args.bams):
        v = tlens(bam, args.region, args.n_reads, args.max_len)
        name = os.path.basename(bam).replace(".pe.bam", "")
        if len(v) == 0:
            print(f"  {name}: no paired reads found, skipping")
            continue
        dens, edges = np.histogram(v, bins=bins, density=True)
        centres = (edges[:-1] + edges[1:]) / 2
        ax.plot(centres, dens, color=color, lw=1.4, label=name)
        sub_f = float(((v >= SUB[0]) & (v <= SUB[1])).mean())
        mono_f = float(((v >= MONO[0]) & (v <= MONO[1])).mean())
        rows.append({"replicate": name, "n_fragments": int(len(v)),
                     "median_length": int(np.median(v)),
                     "mode_length_5bp_bin": int(centres[dens.argmax()]),
                     "frac_sub_1_99": round(sub_f, 4),
                     "frac_mono_180_247": round(mono_f, 4),
                     "frac_neither": round(1 - sub_f - mono_f, 4)})
        print(f"  {name}: n={len(v):,} median={np.median(v):.0f} "
              f"sub={sub_f:.3f} mono={mono_f:.3f}")

    ax.axvspan(SUB[0], SUB[1], color="#0096a0", alpha=0.10)
    ax.axvspan(MONO[0], MONO[1], color="#792374", alpha=0.10)
    ymax = ax.get_ylim()[1]
    ax.text(50, ymax * 0.96, "sub", ha="center", fontsize=8, color="#0f6e56")
    ax.text(213, ymax * 0.96, "mono", ha="center", fontsize=8, color="#3c3489")
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("ATAC fragment lengths and the two channel windows", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("black")
    ax.tick_params(colors="black")
    fig.tight_layout()

    df = pd.DataFrame(rows)
    p = os.path.join(args.outdir, "atac_fragment_length_summary.tsv")
    df.to_csv(p, sep="\t", index=False)
    print(f"\nWrote {p}")
    print(df.to_string(index=False))
    for ext in ("pdf", "png"):
        fp = os.path.join(args.figdir, f"atac_fragment_lengths.{ext}")
        fig.savefig(fp, dpi=300, bbox_inches="tight")
        print(f"Wrote {fp}")


if __name__ == "__main__":
    main()
