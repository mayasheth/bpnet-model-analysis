#!/usr/bin/env python3
"""Genome-wide FIMO scan of the JASPAR 2024 CORE vertebrate motifs, one chromosome per call.

WHY. Test 1 gives the p300 model a TF-activity input: a per-base track, per cell type, equal to
the sum over motifs of (motif hit at this base) x (that motif's activity in the cell). The
motif hits are the same for every cell type and every weighting, so they are scanned ONCE and
the per-cell tracks are built from them afterwards (4.45). The scan has to be genome-wide
because training negatives are drawn from a genome-wide pool.

OUTPUT. One compressed npz per chromosome under $SCRATCH (job output belongs there; Sherlock
/etc/agents/AGENTS.md): motif index (int16), hit start (int32), hit end (int32), -log10 p
(float16), plus the motif names. p < 1e-4, both strands, as in 4.41 and 4.43.

N BASES ARE EXCLUDED. They are encoded as all-zero columns, and a zero PWM score passes
p < 1e-4 for many motifs, so the first chr21 run returned 155M hits of which 94% touched an N
(2026-09-26). A hit is kept only if every base it covers is A, C, G or T; the clean hit rate is
then ~200 per Mb per motif, what p < 1e-4 on two strands predicts.

The chromosome is scanned in chunks with a motif-length overlap so hits spanning a chunk
boundary are found once: a hit is kept only if it STARTS inside its chunk's own span.

Usage: 4.44.scan_genome_motifs.py --chrom chr1 [--out-dir $SCRATCH/motif_hits]
"""
import argparse
import os
import time

import numpy as np

D = "/oak/stanford/groups/engreitz/Users/sheth"
GENOME = f"{D}/hg38_resources/hg38.fa"
JASPAR = f"{D}/EP300_BPNet/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"

ap = argparse.ArgumentParser()
ap.add_argument("--chrom", required=True)
ap.add_argument("--out-dir", default=os.path.join(os.environ.get("SCRATCH", "/tmp"), "motif_hits"))
ap.add_argument("--chunk", type=int, default=2_000_000)
ap.add_argument("--threshold", type=float, default=1e-4)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
out = f"{a.out_dir}/{a.chrom}.npz"
if os.path.exists(out):
    raise SystemExit(f"exists: {out}")

import pyfaidx  # noqa: E402
from memelite import fimo  # noqa: E402
from memelite.io import read_meme  # noqa: E402

motifs = read_meme(JASPAR)
names = list(motifs)
maxw = max(m.shape[1] for m in motifs.values())
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
L = len(fa[a.chrom])
LUT = np.full(256, -1, np.int8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i

mi, st, en, lp = [], [], [], []
t0 = time.time()
for s in range(0, L, a.chunk):
    e = min(L, s + a.chunk + maxw)
    idx = LUT[np.frombuffer(str(fa[a.chrom][s:e]).encode(), np.uint8)]
    ok = idx >= 0
    if ok.sum() < maxw:
        continue
    X = np.zeros((1, 4, e - s), np.float32)
    X[0, idx[ok], np.nonzero(ok)[0]] = 1.0
    ncum = np.concatenate([[0], np.cumsum(~ok)])            # N count up to each position
    for j, h in enumerate(fimo(motifs, X, threshold=a.threshold)):
        if not len(h):
            continue
        hs = h["start"].values.astype(np.int64)
        he_ = h["end"].values.astype(np.int64)
        keep = (hs < a.chunk) & (ncum[np.minimum(he_, len(ok))] - ncum[hs] == 0)
        if not keep.any():
            continue
        hs, he = hs[keep] + s, h["end"].values[keep].astype(np.int64) + s
        mi.append(np.full(len(hs), j, np.int16))
        st.append(hs.astype(np.int32))
        en.append(he.astype(np.int32))
        lp.append((-np.log10(np.maximum(h["p-value"].values[keep], 1e-300))).astype(np.float16))
    print(f"{a.chrom}: {min(s + a.chunk, L) / 1e6:.0f}/{L / 1e6:.0f} Mb "
          f"({time.time() - t0:.0f} s)", flush=True)

np.savez_compressed(out, motif=np.concatenate(mi), start=np.concatenate(st),
                    end=np.concatenate(en), neglog10p=np.concatenate(lp),
                    names=np.array(names, dtype=object), chrom_len=L)
print(f"wrote {out}: {sum(len(x) for x in mi):,} hits")
