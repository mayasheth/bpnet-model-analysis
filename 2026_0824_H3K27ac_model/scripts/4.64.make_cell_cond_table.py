#!/usr/bin/env python3
"""Per-cell conditioning table for train_multimodal_bpnet.py --cell-cond-tsv, from bulk chromVAR.

Input: the wide motif x cell chromVAR z matrix from 4.63 (ATAC only, no p300). Steps:
  1. drop motifs with a non-finite z in any cell (short motifs with no FIMO hits, 4.63 note);
  2. z-score each cell's vector across motifs, so library depth (about 25x across the panel)
     does not set the scale;
  3. PCA fitted on the TRAINING cells only (every cell except --holdout), then every cell,
     the holdout included, is projected onto those axes. The holdout therefore never shapes the
     axes; its row is computed from its own ATAC alone and used only at inference;
  4. keep k = min(--max-pcs, n_train - 1) components, all divided by ONE scale (the training-
     cell sd of pc1). Not per component: with k = n_train - 1, unit-variance components put
     every training cell at the same distance from every other (a whitened simplex), so the
     vector would carry identity and no similarity.

Output columns: cell, pc1..pck. Also writes <out>.loadings.tsv (top motifs per component) and
<out>.json (explained variance, cells, motifs kept).

Usage: 4.64.make_cell_cond_table.py --z-wide results/2026_1020_bulk_chromvar_panel_z_wide.tsv \
           --holdout K562 --out config/cell_cond_chromvar_pc.tsv
"""
import argparse
import json

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--z-wide", required=True, help="4.63 wide z matrix (motif rows, cell columns)")
ap.add_argument("--holdout", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--max-pcs", type=int, default=8)
ap.add_argument("--cells", default=None, help="comma list to restrict to (default: all columns)")
ap.add_argument("--panel-json", default=None,
                help="Write the panel JSON's cell names (e.g. MCF-7) rather than the file-safe column "
                     "names 4.63 uses (MCF_7); every panel cell must match exactly one column.")
a = ap.parse_args()

z = pd.read_csv(a.z_wide, sep="\t")
id_cols = [c for c in z.columns if c in ("motif", "motif_name")]
z = z.set_index(id_cols[0]).drop(columns=id_cols[1:])
if a.cells:
    z = z[a.cells.split(",")]
if a.panel_json:
    keys = list(json.load(open(a.panel_json)))
    safe = {k.replace("-", "_").replace(" ", "_"): k for k in keys}
    missing = [k for s_, k in safe.items() if s_ not in z.columns]
    if missing:
        raise SystemExit(f"error: panel cells with no chromVAR column: {missing}")
    z = z[list(safe)].rename(columns=safe)
if a.holdout not in z.columns:
    raise SystemExit(f"error: holdout {a.holdout} not in {list(z.columns)}")
ok = np.isfinite(z.to_numpy(float)).all(axis=1)
z = z[ok]
print(f"{ok.sum()} motifs kept, {(~ok).sum()} dropped as non-finite; cells {list(z.columns)}")

V = (z - z.mean(axis=0)) / z.std(axis=0, ddof=0)            # per cell, across motifs
train = [c for c in V.columns if c != a.holdout]
Xt = V[train].T.to_numpy(float)                               # cells x motifs
mu = Xt.mean(axis=0)
U, S, Wt = np.linalg.svd(Xt - mu, full_matrices=False)
k = min(a.max_pcs, len(train) - 1)
if k < 1:
    raise SystemExit("error: need at least 2 training cells")
W = Wt[:k].T                                                  # motifs x k
proj = (V.T.to_numpy(float) - mu) @ W                         # all cells x k
sd = proj[[V.columns.get_loc(c) for c in train], 0].std(ddof=0)
proj = proj / sd
ev = (S ** 2) / (S ** 2).sum()

out = pd.DataFrame(proj, columns=[f"pc{i + 1}" for i in range(k)])
out.insert(0, "cell", list(V.columns))
out.to_csv(a.out, sep="\t", index=False, float_format="%.6f")

load = []
for i in range(k):
    w = pd.Series(W[:, i], index=V.index)
    for sign, s in (("+", w.nlargest(8)), ("-", w.nsmallest(8))):
        load += [{"pc": f"pc{i + 1}", "sign": sign, "motif": m, "loading": float(v)} for m, v in s.items()]
pd.DataFrame(load).to_csv(a.out.replace(".tsv", ".loadings.tsv"), sep="\t", index=False)
json.dump({"z_wide": a.z_wide, "holdout": a.holdout, "training_cells": train, "k": k,
           "explained_variance_training": ev[:k].round(4).tolist(), "n_motifs": int(ok.sum())},
          open(a.out.replace(".tsv", ".json"), "w"), indent=2)

print(f"k = {k}; training-cell variance explained {', '.join(f'{v:.2f}' for v in ev[:k])}")
d = pd.DataFrame(proj, index=V.columns)
for c in V.columns:
    dist = np.sqrt(((d.loc[train] - d.loc[c]) ** 2).sum(axis=1)).drop(c, errors="ignore")
    print(f"  {c:<18} nearest training cell {dist.idxmin()} (distance {dist.min():.2f})")
