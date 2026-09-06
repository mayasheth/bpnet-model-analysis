#!/usr/bin/env python3
"""Verify that the GC-matched negatives are actually GC-matched.

`bpnet-gc-background` does not carry the GC column into its output, so the match cannot be
read off the file it writes. This recovers GC by joining the output coordinates back to the
tiling's own annotation -- the same estimator used on both sides of the match, which is the
only way the comparison means anything.

Reports, for each element set: the elements' GC distribution (the matching target), the
matched pool's, and the unmatched pool's, plus the largest per-bin discrepancy. If the
matched column does not sit on the elements column, the retrain that uses these negatives
tests nothing and should not be read.

Positives' GC comes from the same tiling by nearest-bin lookup rather than from the genome
fasta, so that positives and negatives are measured identically.
"""
import os
import numpy as np
import pandas as pd

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
TILING = f"{D}/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
SETS = [("k562_h3k27ac", f"{D}/reference/K562_DNase_candidate_elements.narrowPeak"),
        ("gm12878_h3k27ac",
         f"{D}/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak")]
BIN = 0.02

til = pd.read_csv(TILING, sep="\t", header=None, names=["chr", "start", "end", "gc"])
til["gc"] = pd.to_numeric(til["gc"], errors="coerce")
exact = {(c, s): g for c, s, g in zip(til["chr"], til["start"], til["gc"])}

# Nearest-bin lookup, per chromosome, for regions that are not on the tiling grid.
by_chrom = {}
for c, sub in til.groupby("chr"):
    mid = ((sub["start"] + sub["end"]) // 2).to_numpy()
    o = np.argsort(mid)
    by_chrom[c] = (mid[o], sub["gc"].to_numpy()[o])


def gc_of(df, c_col=0, s_col=1, e_col=2):
    out = np.full(len(df), np.nan)
    ch = df.iloc[:, c_col].to_numpy()
    mid = ((df.iloc[:, s_col] + df.iloc[:, e_col]) // 2).to_numpy()
    st = df.iloc[:, s_col].to_numpy()
    for i in range(len(df)):
        hit = exact.get((ch[i], st[i]))
        if hit is not None and np.isfinite(hit):
            out[i] = hit
            continue
        arr = by_chrom.get(ch[i])
        if arr is None:
            continue
        m, g = arr
        j = np.clip(np.searchsorted(m, mid[i]), 1, len(m) - 1)
        out[i] = g[j - 1] if abs(mid[i] - m[j - 1]) <= abs(m[j] - mid[i]) else g[j]
    return out[np.isfinite(out)]


def line(label, v, ref=None):
    edges = np.arange(0.0, 1.0 + BIN, BIN)
    extra = ""
    if ref is not None and len(ref) and len(v):
        a = np.histogram(ref, bins=edges, density=True)[0] * BIN
        b = np.histogram(v, bins=edges, density=True)[0] * BIN
        extra = f"{np.abs(a - b).max():>12.4f}"
    print(f"{label:<34}{len(v):>10,}{np.mean(v):>10.4f}"
          f"{np.percentile(v, 10):>8.3f}{np.percentile(v, 90):>8.3f}{extra}")


unmatched = til["gc"].dropna().to_numpy()
print(f"{'set':<34}{'n':>10}{'mean GC':>10}{'p10':>8}{'p90':>8}{'max bin diff':>12}")
ok = True
for label, peaks in SETS:
    f = f"{P}/data/gc_negatives_{label}.bed"
    print(f"\n--- {label}")
    el = pd.read_csv(peaks, sep="\t", header=None)
    el_gc = gc_of(el)
    line("candidate elements (target)", el_gc)
    line("unmatched tiling (what we used)", unmatched, el_gc)
    if not os.path.exists(f):
        print(f"  MISSING: {f}")
        ok = False
        continue
    neg = pd.read_csv(f, sep="\t", header=None)
    neg_gc = gc_of(neg)
    line("GC-matched negatives", neg_gc, el_gc)
    d_match = abs(np.mean(neg_gc) - np.mean(el_gc))
    d_un = abs(np.mean(unmatched) - np.mean(el_gc))
    print(f"    mean-GC gap to elements: matched {d_match:.4f}, unmatched {d_un:.4f} "
          f"({d_un / d_match:.1f}x closer)" if d_match > 1e-9 else "")
    if d_match > 0.5 * d_un:
        print("    WARNING: matching barely improved the mean-GC gap; do not read the retrain")
        ok = False

print("\nVERDICT:", "matched pools are usable" if ok
      else "at least one pool is missing or poorly matched")
