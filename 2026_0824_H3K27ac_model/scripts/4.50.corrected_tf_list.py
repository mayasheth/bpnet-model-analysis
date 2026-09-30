#!/usr/bin/env python3
"""Correct the ATAC-derived TF list for motifs that are open but do not recruit p300.

WHY. The deployable TF-activity channel (F-036) uses each cell's within-cell ATAC association
with CTCF removed by hand. In the training cells, 29 motifs have ATAC association more than 1 sd
above their true p300 importance in at least 4 of 5 cells: CTCF/CTCFL and GC-rich,
CpG-island-like motifs (KLF, SP, E2F, ZIC, ZBED4, ZNF417, PATZ1). Classic repressors (REST,
SNAI, ZEB1, YY1) are not among them: their ATAC association is already negative. This scores
corrections of the raw association, each using the TRAINING cells' p300 or the target's own
sequence only, leave one cell out, as 4.48 did.

    raw          within-cell ATAC association (4.47's assoc, standardised within cell)
    noctcf       CTCF-family motifs set to the minimum (the F-036 weighting)
    flagged      motifs whose association exceeds truth by > --flag-sd in >= all-but-one of the
                 OTHER four cells set to the minimum: CTCF removal generalised from data
    gap          raw - mean over the other cells of (association - truth), per motif.
                 NOTE: algebraically this is mean training truth + (target association - mean
                 training association), i.e. a reference-based specificity estimate plus a
                 prior; a target like one training cell loses a quarter of that cell's signal
    gc           raw with motif GC, CpG and GC^2 regressed out, within the target cell only
    gc+flagged   both

Scores against truth (4.43): Spearman, top-20 overlap, the departure from the training-cell mean
truth (the cell-specific part, F-034), CTCF's best rank, GATA in K562, K562's top 10.

Usage: 4.50.corrected_tf_list.py [--flag-sd 1.0]
       [--features results/2026_0928_learned_motif_importance.features.tsv.gz]
       [--out-prefix results/2026_0929_corrected_tf_list]
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
JASPAR = f"{P}/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"

ap = argparse.ArgumentParser()
ap.add_argument("--features", default="results/2026_0928_learned_motif_importance.features.tsv.gz")
ap.add_argument("--flag-sd", type=float, default=1.0)
ap.add_argument("--remove-to", choices=("zero", "min"), default="zero")
ap.add_argument("--out-prefix", default="results/2026_0929_corrected_tf_list")
a = ap.parse_args()

Fe = pd.read_csv(a.features, sep="\t")
CELLS = list(dict.fromkeys(Fe["cell"]))
mnames = list(Fe.loc[Fe["cell"] == CELLS[0], "motif"])
W = {c: Fe.loc[Fe["cell"] == c, "assoc"].values for c in CELLS}
T = {c: Fe.loc[Fe["cell"] == c, "truth"].values for c in CELLS}
gc = Fe.loc[Fe["cell"] == CELLS[0], "motif_gc"].values

from memelite.io import read_meme  # noqa: E402

pwm = read_meme(JASPAR)
cpg = np.array([float((pwm[n][1, :-1] * pwm[n][2, 1:]).sum() / max(pwm[n].shape[1] - 1, 1))
                for n in mnames])
ctcf = np.array([n.split(" ", 1)[1].upper().startswith("CTCF") for n in mnames])
X = np.column_stack([np.ones(len(gc)), gc, gc ** 2, cpg])
z = lambda v: (v - v.mean()) / v.std()


def floor(v, mask):
    """Remove motifs from the list: 0 (no information, as 4.45 does for CTCF in the track) or,
    with --remove-to min, below every other motif."""
    v = v.copy()
    v[mask] = 0.0 if a.remove_to == "zero" else v.min() - 1
    return v


def gc_resid(v):
    b = np.linalg.lstsq(X, v, rcond=None)[0]
    return z(v - X @ b)


EST = {k: {} for k in ("raw", "noctcf", "flagged", "gap", "gc", "gc+flagged")}
FLAG = {}
for c in CELLS:
    oth = [d for d in CELLS if d != c]
    G = np.array([W[d] - T[d] for d in oth])
    flag = (G > a.flag_sd).sum(0) >= len(oth) - 1
    FLAG[c] = flag
    EST["raw"][c] = W[c]
    EST["noctcf"][c] = floor(W[c], ctcf)
    EST["flagged"][c] = floor(W[c], flag | ctcf)
    EST["gap"][c] = W[c] - G.mean(0)
    EST["gc"][c] = gc_resid(W[c])
    EST["gc+flagged"][c] = floor(gc_resid(W[c]), flag | ctcf)
    print(f"{c}: {int(flag.sum())} motifs flagged from the other cells: "
          + ", ".join(sorted({mnames[i].split(' ', 1)[1] for i in np.nonzero(flag)[0]}))[:300],
          flush=True)
print(f"motif GC vs raw association (mean over cells): Spearman "
      f"{np.mean([spearmanr(gc, W[c])[0] for c in CELLS]):+.3f}; CpG {np.mean([spearmanr(cpg, W[c])[0] for c in CELLS]):+.3f}; "
      f"vs truth GC {np.mean([spearmanr(gc, T[c])[0] for c in CELLS]):+.3f}")

top = lambda v, k=20: set(np.argsort(-v)[:k])
rows = []
for name, E in EST.items():
    for c in CELLS:
        prior = np.mean([T[d] for d in CELLS if d != c], 0)
        e, t = E[c], T[c]
        rows.append(dict(estimate=name, cell=c, spearman=spearmanr(e, t)[0],
                         top20=len(top(e) & top(t)),
                         dep_spearman=spearmanr(e - prior, t - prior)[0],
                         dep_top20=len(top(e - prior) & top(t - prior)),
                         ctcf_best_rank=int(rankdata(-e)[ctcf].min())))
R = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== mean over the five held-out cells ===")
print(R.groupby("estimate", sort=False)[["spearman", "top20", "dep_spearman", "dep_top20"]].mean()
      .to_string(float_format=lambda v: f"{v:.3f}"))
print("\n=== per cell: Spearman vs truth / departure Spearman ===")
print(R.pivot(index="estimate", columns="cell", values="spearman").reindex(EST).round(3)
      .astype(str).add(" / ").add(R.pivot(index="estimate", columns="cell", values="dep_spearman")
                                  .reindex(EST).round(3).astype(str)).to_string())

gi = [i for i, n in enumerate(mnames) if "GATA" in n.upper() or "TRPS1" in n.upper()]
panel = [c for c in CELLS if c != "K562"]
dT = T["K562"] - np.mean([T[c] for c in panel], 0)
print("\n=== K562 ===")
print(f"  {'truth':<12} GATA {z(T['K562'])[gi].mean():+.2f}  top10: "
      + ", ".join(mnames[i].split(' ', 1)[1] for i in np.argsort(-T['K562'])[:10]))
for name, E in EST.items():
    e = E["K562"]
    dE = e - np.mean([E[c] for c in panel], 0)
    print(f"  {name:<12} GATA {z(e)[gi].mean():+.2f}  contrast vs true {spearmanr(dE, dT)[0]:+.3f}  "
          "top10: " + ", ".join(mnames[i].split(' ', 1)[1] for i in np.argsort(-e)[:10]))

wide = pd.DataFrame({"motif": mnames})
for name, E in EST.items():
    for c in CELLS:
        wide[f"{name}_{c}"] = E[c]
wide.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
R.to_csv(f"{a.out_prefix}.scores.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.motifs.tsv, .scores.tsv")
