#!/usr/bin/env python3
"""Per-cell TF-activity bigwig: at each base, the sum over motifs of (hit covers this base) x
(the motif's activity weight in this cell type).

INPUTS. Genome-wide motif hits from 4.44 (one npz per chromosome) and a weight per motif per
cell type. Weightings:
    oracle   truth_<cell> from 4.43: the motif's coefficient for p300 BEYOND ATAC in that
             cell, from that cell's own p300. Not deployable; it is Test 1's upper bound.
    assoc_noctcf  within_<cell> from 4.43: the motif's ATAC association within the cell, from
             the cell's own ATAC and sequence only, with the CTCF-family motifs set to 0.
             Deployable. 4.48 (2026-09-28): its departure from the panel mean tracks the
             true departure best of every estimator (K562 Spearman 0.77); CTCF, ranked
             first by raw association, is its one large error in K562.
Weights are standardised within each cell type (divided by their sd across motifs), so the
track has the same scale in every cell type and the model cannot read a cell's identity from
the scale. Motifs 4.43 dropped (fewer than 200 hits) get weight 0.

THE TRACK IS RESCALED TO GENOME-WIDE SD 1 (over every base of the scanned chromosomes, zeros
included). Standardising the weights alone left the channel's spread cell-dependent (sd 8.6 in
A549 against 4.3 in MCF-7 over training windows, 2026-09-27), because it depends on how many
large-weight motifs a cell has, and the model could read cell identity off that scale. The
rescaling uses only the cell's own track, not any choice of regions.

The track is written as a run-length bigwig over the bases some weighted hit covers; every
other base is absent, which the extractor reads as 0.

Usage: 4.45.build_activity_track.py --cell K562 --weighting oracle
       [--hits-dir $SCRATCH/motif_hits] [--out-dir $SCRATCH/motif_activity]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyBigWig

A = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
REC = f"{A}/results/2026_0926_tf_activity_recovery.motifs.tsv"
RNAW = f"{A}/results/2026_0930_rna_tf_weights.motifs.tsv"
WEIGHTS = {"oracle": (REC, "truth_{cell}", False), "assoc_noctcf": (REC, "within_{cell}", True),
           # 4.54 (2026-09-30): TF expression, family maximum; and ATAC association + that
           "rnafam": (RNAW, "rnafam_{cell}", False), "rnaatac": (RNAW, "rnaatac_{cell}", False)}
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SCR = os.environ.get("SCRATCH", "/tmp")

ap = argparse.ArgumentParser()
ap.add_argument("--cell", required=True)
ap.add_argument("--weighting", required=True, choices=sorted(WEIGHTS))
ap.add_argument("--hits-dir", default=f"{SCR}/motif_hits")
ap.add_argument("--out-dir", default=f"{SCR}/motif_activity")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

path, col, drop_ctcf = WEIGHTS[a.weighting]
W = pd.read_csv(path, sep="\t")
w = W.set_index("motif")[col.format(cell=a.cell)].fillna(0.0)   # unmapped motifs: no weight
if drop_ctcf:
    is_ctcf = w.index.str.split(" ", n=1).str[1].str.upper().str.startswith("CTCF")
    print(f"CTCF-family motifs set to 0: {', '.join(w.index[is_ctcf])}", flush=True)
    w[is_ctcf] = 0.0
w = w / w.std()
files = {os.path.basename(f)[:-4]: f for f in glob.glob(f"{a.hits_dir}/chr*.npz")}
missing = [c for c in CHROMS if c not in files]
assert not missing, f"motif hits missing for {missing}"

out = f"{a.out_dir}/motifact_{a.weighting}_{a.cell.replace('-', '_')}.bw"


def chrom_track(c, wv):
    z = np.load(files[c], allow_pickle=True)
    L = int(z["chrom_len"])
    m, s, e = z["motif"].astype(int), z["start"].astype(np.int64), z["end"].astype(np.int64)
    keep = wv[m] != 0
    diff = np.zeros(L + 1, np.float64)
    np.add.at(diff, s[keep], wv[m[keep]])
    np.add.at(diff, np.minimum(e[keep], L), -wv[m[keep]])
    val = np.cumsum(diff[:-1])
    val[np.abs(val) < 1e-9] = 0.0
    return val


names0 = list(np.load(files[CHROMS[0]], allow_pickle=True)["names"])
wv0 = np.array([w.get(n, 0.0) for n in names0], np.float64)
n_tot, s1, s2 = 0, 0.0, 0.0
for c in CHROMS:                                   # pass 1: genome-wide sd of the track
    v = chrom_track(c, wv0)
    n_tot += len(v); s1 += v.sum(); s2 += (v * v).sum()
sd = np.sqrt(s2 / n_tot - (s1 / n_tot) ** 2)
print(f"genome-wide sd before rescaling: {sd:.4f}", flush=True)
bw = None
tot = 0
for c in CHROMS:
    z = np.load(files[c], allow_pickle=True)
    names = list(z["names"])
    L = int(z["chrom_len"])
    if bw is None:
        wv = np.array([w.get(n, 0.0) for n in names], np.float64)
        print(f"{a.cell} {a.weighting}: {int((wv != 0).sum())} of {len(names)} motifs weighted; "
              f"top: {', '.join(n for n in np.array(names)[np.argsort(-wv)[:5]])}", flush=True)
        bw = pyBigWig.open(out, "w")
        sizes = {}
        for cc in CHROMS:
            sizes[cc] = int(np.load(files[cc], allow_pickle=True)["chrom_len"])
        bw.addHeader([(cc, sizes[cc]) for cc in CHROMS])
    val = (chrom_track(c, wv) / sd).astype(np.float32)
    # run-length encode, dropping zero runs
    change = np.nonzero(np.diff(val) != 0)[0] + 1
    starts = np.concatenate([[0], change])
    ends = np.concatenate([change, [L]])
    vals = val[starts]
    nz = vals != 0
    bw.addEntries([c] * int(nz.sum()), starts[nz].tolist(), ends=ends[nz].tolist(),
                  values=vals[nz].astype(float).tolist())
    tot += int((ends[nz] - starts[nz]).sum())
    print(f"  {c}: {int(nz.sum()):,} runs", flush=True)
bw.close()
print(f"wrote {out}; {tot / 1e6:.0f} Mb covered")
