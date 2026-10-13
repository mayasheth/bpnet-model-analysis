#!/usr/bin/env python3
"""10.1 Within-gene eQTL candidates for one scE2G cluster, scored by every plan 9 arm (plan 10).

Pre-registered in decisions [2026-10-13] (plan 10). Plan 9's enrichment readouts barely see the
activity term (F-065); here the RNA terms are held constant by keeping only expressed genes, and the
readout asks which of one gene's candidate elements carries its fine-mapped eQTL.

Per matched eQTL cell type of the biosample:
  candidates  scE2G pairs with gene in the universe, pseudobulk TPM >= 1, class != promoter,
              distance <= 1 Mb
  positive    the element contains a plan 9 variant (PIP >= 0.5, distal noncoding) for that gene
              in that cell type
  kept        only genes with >= 1 positive candidate, so the set is identical for every arm
Scores are 9.1's: each arm's logistic score on scE2G's quantile scale (rank over ALL pairs of the
cluster, zeros at 0, TPM < 1 -> 0), plus v3_published (the file's own E2G.Score, same rule).

CHECKS: the quantile + TPM rule reproduces E2G.Score.qnorm (as 9.1); elements do not overlap.

OUTPUT ($OUT/withingene/<biosample>.tsv.gz): tissue, eid, TargetGene, distance, positive, one
quantile column per arm; <biosample>.stats.json

Usage: 10.1.eqtl_withingene_cluster.py --biosample PBMC_11_CD14.Mono.1 [--out-dir $SCRATCH/eqtl_transfer]
"""
import argparse
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
import pysam
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("prep", os.path.join(HERE, "9.0.eqtl_transfer_prep.py"))
prep = importlib.util.module_from_spec(spec)
_argv, sys.argv = sys.argv, [sys.argv[0]]
spec.loader.exec_module(prep)
sys.argv = _argv

SC = os.environ["SCRATCH"]
TPM_MIN = 1.0
MAX_DIST = 1_000_000
ap = argparse.ArgumentParser()
ap.add_argument("--biosample", required=True)
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
a = ap.parse_args()
od = f"{a.out_dir}/withingene"
os.makedirs(od, exist_ok=True)

mt = pd.read_csv(f"{a.out_dir}/matches.tsv", sep="\t")
mt = mt[mt["biosample"] == a.biosample]
assert len(mt), f"{a.biosample} not in matches.tsv"
models = json.load(open(f"{a.out_dir}/models.json"))
COMPUTED = {"cgi_p1", "gata_p1", "ikzf2_p1"}
COLS = ["chr", "start", "end", "TargetGene", "class", "distance", "RNA_pseudobulkTPM", "E2G.Score", "E2G.Score.qnorm"]
NUM0 = len(COLS)
COLS += sorted({f for m in models.values() for f in m["features"]} - COMPUTED - set(COLS))
P = pd.read_csv(mt["pred_file"].iloc[0], sep="\t", usecols=COLS)
P["ubiqExpressed"] = P["ubiqExpressed"].astype(str).str.lower().map({"true": 1.0, "false": 0.0, "1": 1.0, "0": 0.0, "1.0": 1.0, "0.0": 0.0})
assert P["ubiqExpressed"].notna().all()
P = P.replace([np.inf, -np.inf], np.nan).fillna({c: 0 for c in COLS[NUM0:]})
print(f"{a.biosample}: {len(P):,} pairs", flush=True)


def quantile_tpm(score, tpm):
    """As 9.1: run_e2g_cv.calculate_quantiles, then 0 where TPM < TPM_MIN."""
    s = np.asarray(score, float)
    r = rankdata(np.where(s == 0, -np.inf, s), method="average")
    q = (r - 1) / (len(s) - 1)
    q[s == 0] = 0
    q[np.asarray(tpm) < TPM_MIN] = 0
    return q


q3 = quantile_tpm(P["E2G.Score"], P["RNA_pseudobulkTPM"])
rho = spearmanr(q3, P["E2G.Score.qnorm"])[0]
same0 = np.mean((q3 == 0) == (P["E2G.Score.qnorm"].values == 0))
assert rho > 0.999 and same0 > 0.999, f"quantile + TPM rule does not reproduce qnorm (rho {rho:.5f}, zeros {same0:.5f})"

el = P[["chr", "start", "end"]].drop_duplicates().sort_values(["chr", "start"]).reset_index(drop=True)
for c, g in el.groupby("chr"):
    assert (g["start"].values[1:] >= g["end"].values[:-1]).all(), f"{c}: overlapping elements"
el["mid"] = (el["start"] + el["end"]) // 2
G = pd.read_csv(prep.CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])
H = pd.read_csv(f"{a.out_dir}/motif_sites_gata_ikzf2.tsv.gz", sep="\t")
el = pd.concat([el, prep.element_features(el["chr"].values, el["mid"].values, G, H, pysam.FastaFile(prep.GENOME))], axis=1)
el["eid"] = np.arange(len(el))
P = P.merge(el[["chr", "start", "end", "eid", "cgi_p1", "gata_p1", "ikzf2_p1"]], on=["chr", "start", "end"], how="left", validate="m:1")

# quantiles over ALL pairs of the cluster (the deployed scale), before any restriction
Q = {"v3_published": q3}
for arm, m in models.items():
    X = np.log(np.abs(P[m["features"]].values.astype(float)) + m["eps"])
    Q[arm] = quantile_tpm(1 / (1 + np.exp(-(X @ np.array(m["coef"]) + m["intercept"]))), P["RNA_pseudobulkTPM"])
names = list(Q)

uni = set(pd.read_csv(prep.UNIVERSE, sep="\t", header=None, usecols=[3])[3])
cand = (P["TargetGene"].isin(uni) & (P["RNA_pseudobulkTPM"] >= TPM_MIN) & (P["class"] != "promoter")
        & (P["distance"] <= MAX_DIST)).values
C = P.loc[cand, ["eid", "TargetGene", "distance"]].assign(**{k: Q[k][cand] for k in names})
# a gene's pair with one element appears once (no overlapping elements); assert it
assert not C.duplicated(["eid", "TargetGene"]).any(), "duplicate (element, gene) candidates"

V = pd.read_csv(f"{a.out_dir}/variants.tsv.gz", sep="\t")
V = V[V["tissue"].isin(set(mt["tissue"]))].reset_index(drop=True)
V["eid"] = -1
for c, g in V.groupby("chr"):
    e = el[el["chr"] == c]
    if e.empty:
        continue
    k = np.searchsorted(e["start"].values, g["start"].values, side="right") - 1
    ok = (k >= 0) & (e["end"].values[np.maximum(k, 0)] > g["start"].values)
    V.loc[g.index.values[ok], "eid"] = e["eid"].values[k[ok]]
pos = V.loc[V["eid"] >= 0, ["tissue", "eid", "gene_hgnc"]].drop_duplicates().rename(columns={"gene_hgnc": "TargetGene"})

out, st_t = [], {}
for t in sorted(set(mt["tissue"])):
    pt = pos[pos["tissue"] == t].assign(positive=1).drop(columns="tissue")
    ct = C.merge(pt, on=["eid", "TargetGene"], how="left")
    ct["positive"] = ct["positive"].fillna(0).astype(np.int8)
    eg = set(ct.loc[ct["positive"] == 1, "TargetGene"])
    ct = ct[ct["TargetGene"].isin(eg)]
    n_eg_all = V.loc[V["tissue"] == t, "gene_hgnc"].nunique()
    st_t[t] = dict(n_egenes_variant_table=int(n_eg_all), n_egenes_kept=len(eg), n_candidates=int(len(ct)),
                   n_positive=int(ct["positive"].sum()), median_candidates_per_gene=float(ct.groupby("TargetGene").size().median()) if len(ct) else 0.0)
    out.append(ct.assign(tissue=t))
    print(f"  {t}: {st_t[t]}", flush=True)
O = pd.concat(out)[["tissue", "eid", "TargetGene", "distance", "positive"] + names]
O.to_csv(f"{od}/{a.biosample}.tsv.gz", sep="\t", index=False, float_format="%.6g")
st = dict(biosample=a.biosample, n_pairs=int(len(P)), n_candidates_all=int(cand.sum()), qnorm_check_rho=float(rho),
          qnorm_check_zeros=float(same0), per_tissue=st_t)
json.dump(st, open(f"{od}/{a.biosample}.stats.json", "w"), indent=1)
print("done", flush=True)
