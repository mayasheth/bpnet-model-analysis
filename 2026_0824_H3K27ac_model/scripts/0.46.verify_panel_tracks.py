#!/usr/bin/env python3
"""Are the new panel tracks comparable to the existing ones, and built the same way?

TWO DIFFERENT QUESTIONS, AND ONLY ONE OF THEM IS ABOUT THE NEW CELL TYPES.

  COMPARABILITY. Signal-to-noise at each cell type's own EP300 peaks. A new cell type whose
  enrichment is far below K562's contributes mostly noise to a pooled training set, and the
  leave-one-out result for it would be uninterpretable. This is descriptive: cell types
  genuinely differ, so there is no threshold, only an ordering to look at.

  CONSTRUCTION. Whether 0.45's code path reproduces how the existing tracks were made. This
  is the one that can silently poison the panel, because a phase or scale error in the new
  tracks would look like a biological transfer failure. 0.31 made the same argument when it
  built a full-depth control by the same path as its subsampled track.

THE CONSTRUCTION CONTROL IS GM12878 AND ONLY GM12878. Rebuilding a cell type with 0.45 and
correlating against its existing track requires that cell type's source BAMs. GM12878 has
both ATAC and EP300 BAMs on disk. K562's ATAC exists only as tagAlign (the PE BAMs 0.22 used
were on $SCRATCH and are gone), so K562 cannot be rebuilt from source here and is not a
control. One control is enough to catch a systematic construction error, which is what this
is for; it is not enough to catch a cell-type-specific one, and nothing here claims otherwise.

Correlation is computed on chr8 at 1 bp, over that cell type's peaks, because a genome-wide
correlation is dominated by the empty 97% and would read ~1.0 for almost any pair of tracks.

Usage:
  0.46.verify_panel_tracks.py --panel-dir data/p300_panel [--rebuilt-prefix GM12878]
"""
import argparse
import glob
import os

import numpy as np
import pyBigWig

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
EXISTING = {
    "K562": {
        "ep300_plus": f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig",
        "atac": f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
        "peaks": f"{D}/reference/ENCSR000EGE_peaks_inliers.narrowPeak",
    },
    "GM12878": {
        # ENCFF557UDP, the OBSERVED plus strand. NOT the 2026_0606 model's
        # ENCFF960OFK_plus.bw, which F-025 identified as a BPNet model's PREDICTED plus
        # strand. Comparing a rebuild against that file measures how well someone else's
        # model predicts the data, which is not a construction control.
        "ep300_plus": f"{D}/2026_0824_H3K27ac_model/data/p300_panel/ref/ENCFF557UDP.bigWig",
        "atac": f"{D}/2026_0606_GM12878_transferability/data/atac.bw",
        "peaks": "/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878/EP300/"
                 "ENCFF926AKK.bed.gz",
    },
}

ap = argparse.ArgumentParser()
ap.add_argument("--panel-dir", default="data/p300_panel")
ap.add_argument("--rebuilt-prefix", default=None,
                help="cell type also rebuilt by 0.45, for the construction control")
ap.add_argument("--chrom", default="chr8")
ap.add_argument("--half", type=int, default=500)
a = ap.parse_args()


def load_peaks(path, chrom, limit=4000):
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 3 and p[0] == chrom:
                out.append((int(p[1]) + int(p[2])) // 2)
    return np.array(sorted(out)[:limit])


def windows(bw_path, centres, half):
    b = pyBigWig.open(bw_path)
    if a.chrom not in b.chroms():
        b.close()
        return None
    size = b.chroms()[a.chrom]
    rows = []
    for c in centres:
        s, e = int(c) - half, int(c) + half
        if s < 0 or e > size:
            continue
        v = b.values(a.chrom, s, e, numpy=True)
        if v is not None:
            rows.append(np.nan_to_num(v, nan=0.0))
    b.close()
    return np.stack(rows) if rows else None


def bg_mean(bw_path, centres, half, seed=0):
    """Background from RANDOM windows read exactly like the peak windows.

    The first version called pyBigWig's stats(..., "mean"), which averages over the bases
    the bigwig actually stores, while the peak windows are read with values() and count
    absent bases as 0. On a sparse 5'-end track that is two different denominators, and it
    produced enrichments of 0.1x, i.e. peaks apparently quieter than background. Sampling
    random windows through the same `windows()` path makes the two numbers comparable.

    Random positions are drawn on the same chromosome and matched in count to the peak set,
    so the comparison is like for like apart from location.
    """
    b = pyBigWig.open(bw_path)
    if a.chrom not in b.chroms():
        b.close()
        return None
    size = b.chroms()[a.chrom]
    b.close()
    rng = np.random.default_rng(seed)
    rand = rng.integers(half, max(size - half, half + 1), size=max(len(centres), 1))
    W = windows(bw_path, rand, half)
    return float(W.mean()) if W is not None else None


print("=" * 78)
print("1. COMPARABILITY: EP300 enrichment at each cell type's own peaks")
print("=" * 78)
print(f"{'cell':<10}{'source':<10}{'peaks':>8}{'peak mean':>12}{'rand mean':>12}"
      f"{'enrichment':>12}")
rows = []
for cell in ("K562", "GM12878", "A549", "HepG2", "MCF-7"):
    safe = cell.replace("-", "_")
    new_bw = os.path.join(a.panel_dir, f"{safe}_ep300_5p_plus.bw")
    new_pk = os.path.join(a.panel_dir, f"{safe}_ep300_peaks.narrowPeak")
    if os.path.exists(new_bw) and os.path.exists(new_pk):
        bw, pk, src = new_bw, new_pk, "0.45"
    elif cell in EXISTING and os.path.exists(EXISTING[cell]["ep300_plus"]):
        bw, pk, src = EXISTING[cell]["ep300_plus"], EXISTING[cell]["peaks"], "existing"
    else:
        print(f"{cell:<10}{'-':<10}{'(no track yet)':>30}")
        continue
    centres = load_peaks(pk, a.chrom)
    W = windows(bw, centres, a.half)
    if W is None or not len(centres):
        print(f"{cell:<10}{src:<10}{'(no data on ' + a.chrom + ')':>30}")
        continue
    pm, cm = float(W.mean()), bg_mean(bw, centres, a.half)
    enr = pm / max(cm, 1e-12)
    rows.append((cell, enr))
    print(f"{cell:<10}{src:<10}{len(centres):>8,}{pm:>12.4f}{cm:>12.6f}{enr:>12.1f}x")
if len(rows) > 1:
    e = [r[1] for r in rows]
    print(f"\n  enrichment spread {max(e)/max(min(e),1e-9):.2f}x "
          f"({max(rows, key=lambda r: r[1])[0]} highest, "
          f"{min(rows, key=lambda r: r[1])[0]} lowest)")
    print("  Cell types genuinely differ, so read this as an ordering, not a threshold.")
    print("  A cell type far below the rest contributes mostly noise to a pooled set.")

if a.rebuilt_prefix:
    cell = a.rebuilt_prefix
    safe = cell.replace("-", "_")
    print("\n" + "=" * 78)
    print(f"2. CONSTRUCTION CONTROL: 0.45's rebuild of {cell} against its existing tracks")
    print("=" * 78)
    ex = EXISTING.get(cell)
    if not ex:
        raise SystemExit(f"no existing tracks recorded for {cell}")
    centres = load_peaks(ex["peaks"], a.chrom)
    for label, new_name, old_path in (
            ("EP300 5' plus", f"{safe}_ep300_5p_plus.bw", ex["ep300_plus"]),
            ("ATAC coverage", f"{safe}_atac.bw", ex["atac"])):
        new_path = os.path.join(a.panel_dir, new_name)
        if not os.path.exists(new_path):
            print(f"  {label:<16} rebuilt track absent ({new_name}); control not run")
            continue
        A_, B_ = windows(new_path, centres, a.half), windows(old_path, centres, a.half)
        if A_ is None or B_ is None or A_.shape != B_.shape:
            print(f"  {label:<16} shape mismatch {None if A_ is None else A_.shape} "
                  f"vs {None if B_ is None else B_.shape}")
            continue
        x, y = A_.ravel(), B_.ravel()
        r = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else float("nan")
        scale = (y.sum() / x.sum()) if x.sum() else float("nan")
        # A phase error shows up as a correlation that improves when one track is shifted.
        best_shift, best_r = 0, r
        for sh in range(-6, 7):
            if sh == 0:
                continue
            xs = np.roll(A_, sh, axis=1)[:, 6:-6].ravel()
            ys = B_[:, 6:-6].ravel()
            if xs.std() > 0 and ys.std() > 0:
                rr = float(np.corrcoef(xs, ys)[0, 1])
                if rr > best_r:
                    best_shift, best_r = sh, rr
        flag = "" if best_shift == 0 else f"  <-- SHIFTED, best r {best_r:.4f} at {best_shift:+d} bp"
        print(f"  {label:<16} r={r:.4f}  old/new scale={scale:.3f}{flag}")
    print("\n  r near 1.0 with no better shift means 0.45 reproduces the existing path.")
    print("  A scale factor far from 1.0 is usually a depth or replicate-set difference,")
    print("  not a construction error, and is harmless once training normalises.")
