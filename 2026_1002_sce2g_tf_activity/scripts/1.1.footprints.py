#!/usr/bin/env python3
"""Bias-corrected TF footprint scores at motif hits from pseudobulk multiome ATAC (feature B).

WHY. The transferred p300 model mis-ranks elements by which TFs are bound in the target cell
(F-032), and every input tried so far told it which TFs MATTER in the cell, not which ELEMENTS
they are bound at (F-035 to F-038). A footprint is direct, element-level evidence of binding
from the target cell's own ATAC: Tn5 cannot cut where a TF sits, so a bound motif shows fewer
insertions at its core than at its flanks. This is a TOBIAS/HINT-style score with no new
dependencies; scPrinter/PRINT are the better-modelled alternatives if this gate passes.

INSERTIONS. Pseudobulk over all cells in a 10x fragments file: both fragment ends, start and
end - 1 (0-based), which 10x already Tn5-adjusts.

TN5 BIAS. Tn5 prefers some sequences, which by itself makes dips at motifs. The bias is an
8-mer around each insertion (offsets -4..+3), canonical (min of k-mer and its reverse
complement, since both strands' cuts are pooled), estimated as insertions per occurrence over
accessible positions (200 bp bins with >= --min-bin insertions). Pass `bias` writes per-chrom
counts; pass `score` pools them and scores hits.

SCORE, per motif hit (4.44 genome-wide FIMO hits, JASPAR 2024 non-redundant, p < 1e-4):
    core = the motif span, flanks = 30 bp either side
    o = insertions, e = bias weight, both summed per segment; e is rescaled so the window's
        expected total equals its observed total
    FS = log2((o_flank + 1) / (e_flank + 1)) - log2((o_core + 1) / (e_core + 1))
Higher = more protected. Only hits with >= --min-ins insertions in core + flanks are kept, i.e.
accessible ones. Also written: the hit's insertion total (its accessibility).

Usage: 1.1.footprints.py bias  --chrom chr1
       1.1.footprints.py score --chrom chr1
"""
import argparse
import glob
import os
import subprocess

import numpy as np
import pandas as pd
import pyfaidx

FRAG = "/oak/stanford/groups/engreitz/Projects/scE2G/data/K562_260122/K562_Xu/atac_fragments.tsv.gz"
GENOME = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
HITS = os.path.join(os.environ["SCRATCH"], "motif_hits")
OUT = os.path.join(os.environ["SCRATCH"], "footprints_k562_xu")
K, OFF = 8, 4
FLANK = 30

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=("bias", "score"))
ap.add_argument("--chrom", required=True)
ap.add_argument("--frag", default=FRAG)
ap.add_argument("--min-bin", type=int, default=20)
ap.add_argument("--min-ins", type=int, default=10)
a = ap.parse_args()
os.makedirs(OUT, exist_ok=True)

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
L = len(fa[a.chrom])
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
seq = LUT[np.frombuffer(str(fa[a.chrom][:]).encode(), np.uint8)]

# insertion counts per base
proc = subprocess.Popen(["tabix", a.frag, a.chrom], stdout=subprocess.PIPE)
ins = np.zeros(L, np.int32)
n_frag = 0
for chunk in pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[1, 2], names=["s", "e"],
                         dtype={"s": np.int64, "e": np.int64}, chunksize=5_000_000):
    s, e = chunk["s"].values, chunk["e"].values - 1
    for x in (s, e):
        ins += np.bincount(x[(x >= 0) & (x < L)], minlength=L).astype(np.int32)
    n_frag += len(chunk)
proc.wait()
print(f"{a.chrom}: {n_frag:,} fragments, {int(ins.sum()):,} insertions", flush=True)

# canonical 8-mer code at every position (insertion at p uses bases p-4 .. p+3); -1 if any N
code = np.zeros(L, np.int64)
valid = np.ones(L, bool)
for j in range(K):
    sh = np.full(L, 255, np.uint8)
    lo, hi = max(0, OFF - j), min(L, L + OFF - j)
    sh[lo:hi] = seq[lo - OFF + j:hi - OFF + j]
    valid &= sh != 255
    code = code * 4 + np.where(sh == 255, 0, sh)
rc_tab = np.zeros(4 ** K, np.int64)
for c in range(4 ** K):
    x, r = c, 0
    for _ in range(K):
        r = r * 4 + (3 - (x & 3))
        x >>= 2
    rc_tab[c] = r
code = np.minimum(code, rc_tab[code])
code[~valid] = -1

if a.mode == "bias":
    nb = L // 200
    binsum = ins[:nb * 200].reshape(nb, 200).sum(1)
    acc = np.repeat(binsum >= a.min_bin, 200)
    acc = np.concatenate([acc, np.zeros(L - len(acc), bool)]) & (code >= 0)
    occ = np.bincount(code[acc], minlength=4 ** K)
    hit = np.bincount(code[acc], weights=ins[acc], minlength=4 ** K)
    np.savez_compressed(f"{OUT}/bias_{a.chrom}.npz", occ=occ, ins=hit)
    print(f"{a.chrom}: {int(acc.sum()):,} accessible positions, {int(hit.sum()):,} insertions there")
    raise SystemExit(0)

# score: pooled bias
occ = sum(np.load(f)["occ"] for f in glob.glob(f"{OUT}/bias_chr*.npz"))
hit = sum(np.load(f)["ins"] for f in glob.glob(f"{OUT}/bias_chr*.npz"))
rate = (hit + 1) / (occ + 1)
w = rate / np.average(rate, weights=occ)
print(f"bias pooled over {len(glob.glob(f'{OUT}/bias_chr*.npz'))} chromosomes; "
      f"8-mer weight 1st-99th pct {np.percentile(w, 1):.2f}-{np.percentile(w, 99):.2f}", flush=True)
bw = np.where(code >= 0, w[np.maximum(code, 0)], 1.0)
co = np.concatenate([[0], np.cumsum(ins, dtype=np.int64)])
ce = np.concatenate([[0], np.cumsum(bw)])

z = np.load(f"{HITS}/{a.chrom}.npz", allow_pickle=True)
st, en, mo, pv = z["start"].astype(np.int64), z["end"].astype(np.int64), z["motif"], z["neglog10p"]
ok = (st - FLANK >= 0) & (en + FLANK <= L)
st, en, mo, pv = st[ok], en[ok], mo[ok], pv[ok]
seg = lambda c, s, e: c[e] - c[s]
oc, ec = seg(co, st, en), seg(ce, st, en)
of = seg(co, st - FLANK, st) + seg(co, en, en + FLANK)
ef = seg(ce, st - FLANK, st) + seg(ce, en, en + FLANK)
tot = oc + of
keep = tot >= a.min_ins
scale = tot[keep] / (ec[keep] + ef[keep])
fs = (np.log2((of[keep] + 1) / (ef[keep] * scale + 1))
      - np.log2((oc[keep] + 1) / (ec[keep] * scale + 1)))
np.savez_compressed(f"{OUT}/fs_{a.chrom}.npz", start=st[keep].astype(np.int32),
                    end=en[keep].astype(np.int32), motif=mo[keep], neglog10p=pv[keep],
                    n_ins=tot[keep].astype(np.int32), fs=fs.astype(np.float32))
pd.Series(list(z["names"])).to_csv(f"{OUT}/motif_names.tsv", sep="\t", header=False)
print(f"{a.chrom}: {len(st):,} hits, {int(keep.sum()):,} with >= {a.min_ins} insertions; "
      f"FS median {np.median(fs):.3f}")
