#!/usr/bin/env python3
"""Average ATAC fragment profiles around elements with more, or less, p300 or H3K27ac than
their accessibility predicts.

WHY. 4.55 (F-039) found fragment-length features predict H3K27ac beyond ATAC with the
opposite sign to p300. This shows the raw fragments behind that: for each target, the top and
bottom quintile of the residual (target after removing ATAC rank, rank^2, GC, CpG o/e and log
TSS distance, from 4.55's residuals table), so the two groups have matched accessibility.

For 3,000 sampled non-promoter regions per group, paired-end fragments (MAPQ >= 30, proper
pairs) with midpoint within +/- 1 kb of the region centre, from the cell's ATAC BAMs (as 4.51):
    a  fragment length distribution, midpoints within +/- 450 bp, per region
    b  midpoints of short fragments (< 120 bp, open DNA), 20 bp bins, per region
    c  midpoints of mono-nucleosome fragments (150-250 bp), 20 bp bins, per region

Usage: 4.56.nucleosome_profiles.py [--cells K562,GM12878] [--n 3000]
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pysam  # noqa: E402

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
CACHE = f"{A}/data/tf_activity_recovery"
SCR = os.environ.get("SCRATCH", "/tmp")
BAMS = {"K562": sorted(glob.glob(f"{SCR}/atac_pe/ENCFF*.pe.bam")),
        "GM12878": [f"{D}/Data/ENCODE/GM12878/ATAC/{x}.bam"
                    for x in ("ENCFF440GRZ", "ENCFF962FMH", "ENCFF981FXV")]}
TGT = {"p300": ("p300", "#534AB7"), "k27_reg": ("H3K27ac", "#D85A30")}
LBINS = np.arange(20, 705, 5)
PBINS = np.arange(-1000, 1001, 20)

ap = argparse.ArgumentParser()
ap.add_argument("--cells", default="K562,GM12878")
ap.add_argument("--n", type=int, default=3000)
ap.add_argument("--residuals", default=f"{A}/results/2026_0930_nucleosome_vs_h3k27ac.residuals.tsv.gz")
ap.add_argument("--out-prefix", default=f"{A}/results/2026_0930_nucleosome_profiles")
ap.add_argument("--outdir", default=f"{A}/figures")
a = ap.parse_args()

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
})

S = pd.read_csv(f"{CACHE}/union_signals.tsv.gz", sep="\t", usecols=["chr", "start", "end"])
mid_all = ((S["start"] + S["end"]) // 2).values
chr_all = S["chr"].values
RS = pd.read_csv(a.residuals, sep="\t")
rng = np.random.default_rng(0)
rows = []


def scan(bams, us):
    fl_h = np.zeros(len(LBINS) - 1)
    nfr_h = np.zeros(len(PBINS) - 1)
    mono_h = np.zeros(len(PBINS) - 1)
    for u in us:
        ch, mid = chr_all[u], mid_all[u]
        rel, fl = [], []
        for bam in bams:
            for r in bam.fetch(ch, max(0, mid - 1400), mid + 1400):
                if (r.is_secondary or r.is_supplementary or not r.is_proper_pair
                        or r.mapping_quality < 30 or r.template_length <= 0):
                    continue
                t = r.template_length
                rel.append(r.reference_start + t // 2 - mid)
                fl.append(t)
        rel, fl = np.asarray(rel), np.asarray(fl)
        fl_h += np.histogram(fl[np.abs(rel) <= 450], LBINS)[0]
        nfr_h += np.histogram(rel[fl < 120], PBINS)[0]
        mono_h += np.histogram(rel[(fl >= 150) & (fl <= 250)], PBINS)[0]
    n = len(us)
    return fl_h / n, nfr_h / n, mono_h / n


for c in a.cells.split(","):
    idx_dir = f"{SCR}/bam_index"
    bams = []
    for b in BAMS[c]:
        bai = b + ".bai" if os.path.exists(b + ".bai") else f"{idx_dir}/{os.path.basename(b)}.bai"
        bams.append(pysam.AlignmentFile(b, index_filename=bai))
    R = RS[RS["cell"] == c].reset_index(drop=True)
    prof = {}
    for t in TGT:
        q = R[f"resid_{t}"].rank(pct=True).values
        for grp, sel in (("top", q >= 0.8), ("bottom", q <= 0.2)):
            pick = rng.choice(np.flatnonzero(sel), min(a.n, int(sel.sum())), replace=False)
            fl_h, nfr_h, mono_h = scan(bams, R["u"].values[pick])
            prof[(t, grp)] = (fl_h, nfr_h, mono_h)
            atac = R["atac_pct"].values[pick].mean()
            print(f"{c} {t:<8} {grp:<6} n={len(pick):,} mean ATAC pct {atac:.3f}; "
                  f"per region: short {nfr_h.sum():.1f}, mono {mono_h.sum():.1f}", flush=True)
            for kind, h, bins in (("length", fl_h, LBINS), ("short_mid", nfr_h, PBINS),
                                  ("mono_mid", mono_h, PBINS)):
                ctr = (bins[:-1] + bins[1:]) / 2
                rows += [dict(cell=c, target=t, group=grp, kind=kind, x=x, per_region=v,
                              mean_atac_pct=atac) for x, v in zip(ctr, h)]

    fig, ax = plt.subplots(2, 3, figsize=(10, 5.6), sharex="col")
    for i, (t, (lab, col)) in enumerate(TGT.items()):
        for grp, ls, alpha in (("top", "-", 1.0), ("bottom", "--", 0.6)):
            fl_h, nfr_h, mono_h = prof[(t, grp)]
            name = f"{'top' if grp == 'top' else 'bottom'} 20% {lab} beyond ATAC"
            ax[i, 0].plot((LBINS[:-1] + LBINS[1:]) / 2, fl_h / 5, ls, color=col, alpha=alpha, label=name)
            ax[i, 1].plot((PBINS[:-1] + PBINS[1:]) / 2, nfr_h / 20, ls, color=col, alpha=alpha)
            ax[i, 2].plot((PBINS[:-1] + PBINS[1:]) / 2, mono_h / 20, ls, color=col, alpha=alpha)
        ax[i, 0].legend(frameon=False, fontsize=7, loc="upper right")
        ax[i, 0].set_ylabel(f"split by {lab}\nfragments per region per bp")
        for j in (1, 2):
            ax[i, j].axvline(0, color="0.8", lw=0.6)
        for x0 in (147, 294):
            ax[i, 0].axvline(x0, color="0.85", lw=0.6)
    ax[0, 0].set_title("a  fragment lengths (centre +/- 450 bp)")
    ax[0, 1].set_title("b  short fragments (< 120 bp): open DNA")
    ax[0, 2].set_title("c  1-nucleosome fragments (150-250 bp)")
    ax[1, 0].set_xlabel("fragment length (bp); grey = 1 and 2 nucleosomes")
    ax[1, 1].set_xlabel("fragment midpoint vs element centre (bp)")
    ax[1, 2].set_xlabel("fragment midpoint vs element centre (bp)")
    fig.suptitle(f"{c}: ATAC fragments at elements with more (solid) or less (dashed) signal "
                 f"than accessibility predicts; groups matched on ATAC", fontsize=10)
    fig.tight_layout()
    name = f"fig35_nucleosome_profiles_{c.lower()}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {a.outdir}/{name}.png", flush=True)

pd.DataFrame(rows).to_csv(f"{a.out_prefix}.tsv.gz", sep="\t", index=False)
print(f"wrote {a.out_prefix}.tsv.gz")
