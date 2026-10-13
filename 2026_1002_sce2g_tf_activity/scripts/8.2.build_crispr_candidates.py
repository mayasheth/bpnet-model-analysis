#!/usr/bin/env python3
"""8.2 Phase 1 of plans/8_crispr_feature_discovery.md: label-free candidate features on the CRISPR pairs.

SOURCE  the 7.1 library over all 155,134 Xu elements ($SCRATCH/feature_screen/) and the 7.2 p300 composite.
GROUPING, before any label is read:
  motifs       879 JASPAR count columns -> clusters of motifs whose element counts correlate (Spearman
               over all Xu elements) >= --motif-rho, average linkage; cluster value = max count of its
               members (similar motifs hit the same sites, so a sum would double count). Clusters with a
               site in < --min-cluster-frac of Xu elements are dropped: at ~470 positives a rarer motif
               cannot be resolved (rho 0.7 alone left 614 clusters, 2026-10-06)
  dinucleotide 16 frequencies -> top 3 principal components
ENCODING  every candidate = 1 + its percentile over all Xu elements. scE2G takes log(|x| + 0.01) of every
          feature, so the value must be positive and monotone; percentiles also remove outliers and put
          all families on one scale. Label-free and deployable (percentiles from the cell's own elements).
MAPPING   each CRISPR element (training-table chr/start/end, hg38) takes the features of the overlapping Xu
          element with the highest Xu ATAC, so one pair's features come from one element.

WRITES $SCRATCH/crispr_discovery/candidates.tsv.gz (pair keys + candidates), families.tsv, motif_clusters.tsv

Run in the tobias pixi env. Usage: 8.2.build_crispr_candidates.py [--motif-rho 0.7]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

SC = os.environ["SCRATCH"]
LIB = f"{SC}/feature_screen"
TAB = (f"{SC}/sce2g_runs/2026_1005_k562_xu_screen/K562_Xu_pl/multiome_n6_gc/"
       "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz")
KEY = ["chr", "start", "end"]
PAIR = KEY + ["TargetGene", "name"]
ap = argparse.ArgumentParser()
ap.add_argument("--motif-rho", type=float, default=0.7)
ap.add_argument("--min-present", type=float, default=0.005)
ap.add_argument("--min-cluster-frac", type=float, default=0.05,
                help="keep motif clusters with a site in at least this share of Xu elements (label-free power filter)")
ap.add_argument("--out-dir", default=f"{SC}/crispr_discovery")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

E = pd.read_csv(f"{LIB}/library.elements.tsv.gz", sep="\t")
M = np.load(f"{LIB}/library.motifs.npz", allow_pickle=True)
assert (M["chr"] == E["chr"].values).all() and (M["start"] == E["start"].values).all()
C = pd.read_csv(f"{LIB}/screen_composite.tsv.gz", sep="\t")
E = E.merge(C, on=KEY, how="left", validate="1:1")

FAMILIES = {
    "accessibility": ["atac_xu"],
    "composition": ["gc", "cpg_oe", "cpg_n", "at_run_max", "polyAT5_n"],
    "conservation": ["phylop_mean", "phylop_frac_gt2"],
    "repeats": [c for c in E if c.startswith("repeat_")],
    "cpg_island": ["cgi_frac", "cgi_log_dist"],
    "shape": ["width", "mappability"],
    "fragments": ["frag_short_frac", "frag_mono_frac", "frag_mean_log2_len", "ins_centre_frac"],
    "cells": ["cells_frac_open", "cells_mean_count"],
    "neighbourhood": ["neigh_n50k", "neigh_atac50k"],
    "footprints": ["ref_fp_any_tobias", "ref_ery_fp_tobias", "ref_ery_fp_print", "ref_ctcf_fp_print"],
    "p300_composite": ["screen_composite"],
}
F = E[sum(FAMILIES.values(), [])].copy()

# dinucleotide principal components (label-free)
D = E[[c for c in E if c.startswith("dinuc_")]].values
D = (D - D.mean(0)) / D.std(0)
u, sv, vt = np.linalg.svd(D, full_matrices=False)
for k in range(3):
    F[f"dinuc_pc{k + 1}"] = u[:, k] * sv[k]
FAMILIES["composition"] += [f"dinuc_pc{k + 1}" for k in range(3)]
print(f"dinucleotide PCs explain {np.round((sv[:3] ** 2) / (sv ** 2).sum(), 3).tolist()} of variance", flush=True)

# motif clusters (label-free)
names = np.array([str(x) for x in M["names"]])
N = M["n"].astype(np.float32)
keep = (N > 0).mean(0) >= a.min_present
N, names = N[:, keep], names[keep]
R = rankdata(N, axis=0).astype(np.float32)
R = (R - R.mean(0)) / R.std(0)
rho = (R.T @ R) / len(R)
np.fill_diagonal(rho, 1)
Z = linkage(squareform(np.clip(1 - rho, 0, 2), checks=False), method="average")
lab = fcluster(Z, t=1 - a.motif_rho, criterion="distance")
clusters, mcols = [], {}
for k in np.unique(lab):
    idx = np.where(lab == k)[0]
    lead = idx[np.argmax((N[:, idx] > 0).sum(0))]
    cname = f"motif:{names[lead]}" + (f"+{len(idx) - 1}" if len(idx) > 1 else "")
    v = N[:, idx].max(1)
    frac = float((v > 0).mean())
    clusters += [dict(cluster=cname, motif=names[i], n_members=len(idx), present_frac=frac,
                      kept=frac >= a.min_cluster_frac) for i in idx]
    if frac >= a.min_cluster_frac:
        mcols[cname] = v
F = pd.concat([F, pd.DataFrame(mcols, index=F.index)], axis=1)
FAMILIES["motifs"] = sorted(mcols)
pd.DataFrame(clusters).to_csv(f"{a.out_dir}/motif_clusters.tsv", sep="\t", index=False)
n_cl = len({c["cluster"] for c in clusters})
print(f"{keep.sum()} motifs -> {n_cl} clusters at rho >= {a.motif_rho}; {len(mcols)} with a site in >= "
      f"{a.min_cluster_frac:.0%} of Xu elements kept", flush=True)

# encoding: 1 + percentile over all Xu elements (NaN -> lowest)
cand = [c for fam in FAMILIES.values() for c in fam]
P = pd.DataFrame({c: 1 + rankdata(np.nan_to_num(F[c].values.astype(float), nan=-np.inf)) / len(F) for c in cand})
P[KEY] = E[KEY].values

# map CRISPR elements to the overlapping Xu element with the highest Xu ATAC
T = pd.read_csv(TAB, sep="\t", usecols=PAIR + ["Regulated"])
el = T[KEY].drop_duplicates().reset_index(drop=True)
best = np.full(len(el), -1)
for c, g in el.groupby("chr"):
    u_ = E.index[E["chr"] == c].values
    if len(u_) == 0:
        continue
    o = np.argsort(E.loc[u_, "start"].values)
    u_ = u_[o]
    us, ue, at = E.loc[u_, "start"].values, E.loc[u_, "end"].values, E.loc[u_, "atac_xu"].values
    assert (us[1:] >= ue[:-1]).all(), f"{c}: overlapping Xu elements"
    lo = np.searchsorted(ue, g["start"].values, side="right")
    hi = np.searchsorted(us, g["end"].values, side="left")
    for i, l, h in zip(g.index.values, lo, hi):
        if h > l:
            best[i] = u_[l + np.argmax(at[l:h])]
cov = (best >= 0).mean()
print(f"{len(el):,} CRISPR elements, {cov:.4f} overlap a Xu element", flush=True)
assert cov > 0.95, "CRISPR elements do not match the Xu elements"
vals = np.where(best[:, None] >= 0, P[cand].values[np.maximum(best, 0)], 1.0)
el = pd.concat([el, pd.DataFrame(vals, columns=cand, index=el.index)], axis=1)
out = T.merge(el, on=KEY, how="left", validate="m:1")
assert len(out) == len(T)
out.to_csv(f"{a.out_dir}/candidates.tsv.gz", sep="\t", index=False, float_format="%.6g")
fam = pd.DataFrame([dict(feature=c, family=f) for f, cs in FAMILIES.items() for c in cs])
fam.to_csv(f"{a.out_dir}/families.tsv", sep="\t", index=False)
print(fam["family"].value_counts().to_string())
print(f"wrote {a.out_dir}/candidates.tsv.gz: {len(out):,} pairs x {len(cand)} candidates in {len(FAMILIES)} families")
