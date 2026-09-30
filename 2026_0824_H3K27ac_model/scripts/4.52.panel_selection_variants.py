#!/usr/bin/env python3
"""Greedy panel selection variants on 4.49's reliable family departures.

4.49's own greedy run picked two DNase samples first (Ishikawa, H1) on families that are
largely specific to the DNase samples as a group (FEZF2, OTX2, Prdm14), i.e. the assay. This
reruns the same facility-location greedy with the candidate pool restricted to ATAC samples,
with and without the GC-rich / zinc-finger families that 4.50 found to be open-but-not-p300
across the training cells. The departures themselves are 4.49's, unchanged.

Usage: 4.52.panel_selection_variants.py
       [--reliable results/2026_0929_panel_tf_diversity_loose.family_reliable.tsv]
       [--manifest config/panel_candidates_loose.tsv] [--k 6]
       [--out results/2026_0930_panel_selection_variants.tsv]
"""
import argparse

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--reliable", default="results/2026_0929_panel_tf_diversity_loose.family_reliable.tsv")
ap.add_argument("--manifest", default="config/panel_candidates_loose.tsv")
ap.add_argument("--base", default="GM12878,HepG2,A549,MCF-7")
ap.add_argument("--holdout", default="K562")
ap.add_argument("--k", type=int, default=6)
ap.add_argument("--out", default="results/2026_0930_panel_selection_variants.tsv")
a = ap.parse_args()

R = pd.read_csv(a.reliable, sep="\t", index_col=0)
man = pd.read_csv(a.manifest, sep="\t").set_index("cell")
base = a.base.split(",")
atac = [c for c in R.columns if man.loc[c, "assay"] == "ATAC-seq"]
dnase = [c for c in R.columns if c not in atac]
d = R[dnase].mean(1) - R[atac].mean(1)
print("families most specific to the DNase samples as a group: "
      + ", ".join(d.sort_values(ascending=False).head(8).index))
gcrich = R.index.str.upper().str.match(
    r"^(ZNF|ZFP|ZBTB|ZSCAN|ZKSCAN|ZBED|ZIC|PRDM|E2F|KLF|SP[0-9]|PLAG|HINFP|NRF1|NFY)")
rows = []
for label, fams in (("ATAC candidates, all families", R.index),
                    ("ATAC candidates, GC-rich/zinc-finger families removed", R.index[~gcrich])):
    P = R.loc[fams]
    pool = [c for c in atac if c not in base + [a.holdout]]
    cov = lambda S: float(P[S].max(1).sum())
    k = P[a.holdout].sort_values(ascending=False).head(10)
    hcov = lambda S: float(np.mean([P.loc[i, S].max() >= 0.5 * k[i] for i in k.index]))
    S = list(base)
    print(f"\n=== {label}: extending {', '.join(base)} (K562 top-10 covered "
          f"{100 * hcov(S):.0f}%; K562 top: {', '.join(k.index[:6])}) ===")
    rows.append(dict(run=label, step=0, cell="current panel", gain=np.nan, coverage=cov(S),
                     holdout_top10_covered=hcov(S), adds=""))
    for step in range(1, a.k + 1):
        g = {c: cov(S + [c]) - cov(S) for c in pool if c not in S}
        if not g:
            break
        b = max(g, key=g.get)
        add = (P[S + [b]].max(1) - P[S].max(1)).sort_values(ascending=False).head(4)
        S.append(b)
        why = ", ".join(add.index)
        print(f"  {step}. {b:<28} +{g[b]:5.1f}  K562 top-10 covered {100 * hcov(S):3.0f}%  adds {why}")
        rows.append(dict(run=label, step=step, cell=b, gain=g[b], coverage=cov(S),
                         holdout_top10_covered=hcov(S), adds=why))
pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
print(f"\nwrote {a.out}")
