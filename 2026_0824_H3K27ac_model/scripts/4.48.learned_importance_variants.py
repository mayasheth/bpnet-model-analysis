#!/usr/bin/env python3
"""Two ways to put the target cell's own signal back into 4.47's learned motif importance.

WHY. 4.47's leave-one-cell-out estimator beats raw ATAC association on overall Spearman in every
cell (0.48-0.73 vs 0.12-0.57), but mostly by learning which motifs matter for p300 beyond ATAC in
ANY cell: its K562-minus-panel contrast correlates +0.10 with the true contrast (raw association
+0.45) and it scores GATA in K562 at +0.84 (truth +3.87, raw association +3.22). Test 1 needs the
cell-specific part. Both variants use only the target cell's own ATAC and sequence at deployment,
plus quantities learned from the training cells' p300; neither compares the target cell's ATAC
to the training cells', so a new cell similar to a training cell does not cancel its own TFs.

    stack      truth ~ b1 * learned + b2 * assoc, fitted by least squares on INNER
               leave-one-out predictions (for held-out cell c, the learned score of each
               training cell d comes from a model fitted on the other three).
    departure  prior + predicted departure. prior(motif) = mean truth over the training cells.
               Departure model: gradient boosting on the same own-cell features, target =
               truth_d - mean truth over the OTHER training cells, so the cell-invariant part
               cannot explain it. Prediction for c = prior + departure model(c's features).

Also reported: departure-only (no prior), and how well each estimate's departure from the prior
tracks the true departure, the quantity 4.47's learned score gets wrong.

Reads 4.47's features table (features and truth already standardised within cell); no rescan.
Usage: 4.48.learned_importance_variants.py
       [--features results/2026_0928_learned_motif_importance.features.tsv.gz]
       [--out-prefix results/2026_0928_learned_importance_variants]
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

A = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
ANAT = f"{A}/results/2026_0925_p300_error_anatomy.features.tsv"
FEATS = ["assoc", "tail", "prom_frac", "prevalence", "ctcf_cooc", "fp_core", "fp_flank",
         "motif_gc"]

ap = argparse.ArgumentParser()
ap.add_argument("--features", default="results/2026_0928_learned_motif_importance.features.tsv.gz")
ap.add_argument("--out-prefix", default="results/2026_0928_learned_importance_variants")
a = ap.parse_args()

Fe = pd.read_csv(a.features, sep="\t")
CELLS = list(dict.fromkeys(Fe["cell"]))
mnames = list(Fe.loc[Fe["cell"] == CELLS[0], "motif"])
nm = len(mnames)
for c in CELLS:
    assert list(Fe.loc[Fe["cell"] == c, "motif"]) == mnames, f"motif order differs in {c}"
X = {c: Fe.loc[Fe["cell"] == c, FEATS].values for c in CELLS}
T = {c: Fe.loc[Fe["cell"] == c, "truth"].values for c in CELLS}
ASSOC = {c: Fe.loc[Fe["cell"] == c, "assoc"].values for c in CELLS}
LEARNED = {c: Fe.loc[Fe["cell"] == c, "pred"].values for c in CELLS}   # 4.47 leave-one-out


def gb():
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
                                         min_samples_leaf=20, random_state=0)


def fit_learned(train):
    m = gb()
    m.fit(np.vstack([X[d] for d in train]), np.concatenate([T[d] for d in train]))
    return m


def fit_departure(train):
    ys = [T[d] - np.mean([T[e] for e in train if e != d], 0) for d in train]
    m = gb()
    m.fit(np.vstack([X[d] for d in train]), np.concatenate(ys))
    return m


z = lambda v: (v - v.mean()) / v.std()
EST = {k: {} for k in ("stack", "departure", "departure only")}
COEF = {}
for c in CELLS:
    train = [d for d in CELLS if d != c]
    inner_l, inner_a, inner_t = [], [], []
    for d in train:
        m = fit_learned([e for e in train if e != d])
        inner_l.append(z(m.predict(X[d])))
        inner_a.append(ASSOC[d])
        inner_t.append(T[d])
    Z = np.column_stack([np.concatenate(inner_l), np.concatenate(inner_a)])
    b = np.linalg.lstsq(Z, np.concatenate(inner_t), rcond=None)[0]
    COEF[c] = b
    EST["stack"][c] = b[0] * z(LEARNED[c]) + b[1] * ASSOC[c]
    prior = np.mean([T[d] for d in train], 0)
    dep = fit_departure(train).predict(X[c])
    EST["departure"][c] = prior + dep
    EST["departure only"][c] = dep
    print(f"{c}: stack coefficients learned {b[0]:+.3f}, assoc {b[1]:+.3f}", flush=True)

ctcf = np.array([n.split(" ", 1)[1].upper().startswith("CTCF") for n in mnames])
top = lambda v, k=20: set(np.argsort(-v)[:k])
cands = {"raw ATAC association": ASSOC, "learned (4.47)": LEARNED,
         "stack: learned + assoc": EST["stack"], "departure: prior + model": EST["departure"],
         "departure model only": EST["departure only"]}
rows = []
for c in CELLS:
    prior = np.mean([T[d] for d in CELLS if d != c], 0)
    tdep = T[c] - prior
    for name, E in cands.items():
        e = E[c]
        rows.append(dict(cell=c, estimate=name, spearman=spearmanr(e, T[c])[0],
                         top20_overlap=len(top(e) & top(T[c])),
                         dep_spearman=spearmanr(e - prior, tdep)[0],
                         dep_top20=len(top(e - prior) & top(tdep)),
                         ctcf_best_rank=int(rankdata(-e)[ctcf].min())))
Rt = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== estimate vs truth, per held-out cell ===")
print("dep_* score (estimate - prior) against (truth - prior); prior = mean truth over the other "
      "four cells")
print(Rt.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\nmean over cells:")
print(Rt.groupby("estimate", sort=False)[["spearman", "top20_overlap", "dep_spearman",
                                          "dep_top20"]].mean()
      .to_string(float_format=lambda v: f"{v:.3f}"))

panel = [c for c in CELLS if c != "K562"]
dT = T["K562"] - np.mean([T[c] for c in panel], 0)
an = pd.read_csv(ANAT, sep="\t")
an = an[(an.cls == "nonpromoter") & an.feature.str.startswith("motif:")].copy()
rt = an.assign(motif=an.feature.str[6:]).set_index("motif")["coef_r_transfer"].reindex(mnames)
ok = rt.notna().values
print("\n=== K562 minus panel ===")
for name, E in cands.items():
    dE = E["K562"] - np.mean([E[c] for c in panel], 0)
    print(f"  {name:<26} vs true contrast {spearmanr(dE, dT)[0]:+.3f} (top-20 shared "
          f"{len(top(dE) & top(dT))});  vs 4.41 transfer residual "
          f"{spearmanr(dE[ok], rt.values[ok])[0]:+.3f}")
gi = [i for i, n in enumerate(mnames) if "GATA" in n.upper() or "TRPS1" in n.upper()]
print("\nGATA family, mean standardised score (each estimate standardised within cell):")
print(f"  {'truth':<26} " + ", ".join(f"{c} {z(T[c])[gi].mean():+.2f}" for c in CELLS))
for name, E in cands.items():
    print(f"  {name:<26} " + ", ".join(f"{c} {z(E[c])[gi].mean():+.2f}" for c in CELLS))
print("\nK562 top 10 by each estimate:")
for name, E in [("truth", T)] + list(cands.items()):
    print(f"  {name:<26} " + ", ".join(mnames[i].split(' ', 1)[1] for i in np.argsort(-E['K562'])[:10]))

wide = pd.DataFrame({"motif": mnames})
for key, E in (("stack", EST["stack"]), ("departure", EST["departure"])):
    for c in CELLS:
        wide[f"{key}_{c}"] = E[c]
wide.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
Rt.to_csv(f"{a.out_prefix}.scores.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.motifs.tsv (stack_<cell>, departure_<cell>), .scores.tsv")
