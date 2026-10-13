#!/usr/bin/env python3
"""9.1 Score one scE2G cluster with each K562-trained arm and reduce it against its matched eQTLs.

Reads the cluster's unthresholded scE2G predictions (2025_0318_all_clusters, multiome_powerlaw_v3),
which carry every n6 feature. Adds the three features computed on the cluster's own elements (9.0
functions; percentiles over this cluster's elements only, so a cell type's scores never depend on
other cell types). Each arm's score is turned into scE2G's quantile scale exactly as
run_e2g_cv.py does before qnorm: quantile = (rank - 1) / (N - 1) over all pairs, zeros at 0, then 0
for genes with pseudobulk TPM < 1. qnorm only maps that quantile onto a fixed reference curve, so a
threshold on the quantile is a threshold on qnorm. Pairs are then filtered to the gene universe, as
eQTLEnrichment does.

CHECKS: the cluster's elements do not overlap (one element per variant); the quantile + TPM rule
applied to the published E2G.Score reproduces E2G.Score.qnorm (Spearman > 0.999, same zeros).

OUTPUT ($OUT/clusters/<biosample>.*):
  variants.tsv.gz  one row per matched (variant, eGene, cell type): per arm, the element quantile
                   (max over the containing element's genes) and the link quantile (that element's
                   pair with the eGene), -1 when absent
  bg.tsv           per arm, background SNPs inside elements whose quantile >= each grid threshold
  stats.json       counts and checks

Usage: 9.1.eqtl_transfer_cluster.py --biosample PBMC_11_CD14.Mono.1 [--out-dir $SCRATCH/eqtl_transfer]
"""
import argparse
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("prep", os.path.join(HERE, "9.0.eqtl_transfer_prep.py"))
prep = importlib.util.module_from_spec(spec)
import sys  # noqa: E402
_argv, sys.argv = sys.argv, [sys.argv[0]]
spec.loader.exec_module(prep)
sys.argv = _argv

SC = os.environ["SCRATCH"]
TPM_MIN = 1.0
GRID = np.unique(np.r_[np.linspace(0, 0.9, 91), 1 - np.logspace(-1, -5, 400)])
ap = argparse.ArgumentParser()
ap.add_argument("--biosample", required=True)
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
a = ap.parse_args()
od = f"{a.out_dir}/clusters"
os.makedirs(od, exist_ok=True)

mt = pd.read_csv(f"{a.out_dir}/matches.tsv", sep="\t")
mt = mt[mt["biosample"] == a.biosample]
assert len(mt), f"{a.biosample} not in matches.tsv"
models = json.load(open(f"{a.out_dir}/models.json"))
ARMS = list(models)
COMPUTED = {"cgi_p1", "gata_p1", "ikzf2_p1"}
COLS = ["chr", "start", "end", "TargetGene", "RNA_pseudobulkTPM", "E2G.Score", "E2G.Score.qnorm"]
COLS += sorted({f for m in models.values() for f in m["features"]} - COMPUTED - set(COLS))
P = pd.read_csv(mt["pred_file"].iloc[0], sep="\t", usecols=COLS)
P["ubiqExpressed"] = P["ubiqExpressed"].astype(str).str.lower().map({"true": 1.0, "false": 0.0, "1": 1.0, "0": 0.0, "1.0": 1.0, "0.0": 0.0})
assert P["ubiqExpressed"].notna().all()
P = P.replace([np.inf, -np.inf], np.nan).fillna({c: 0 for c in COLS[7:]})
print(f"{a.biosample}: {len(P):,} pairs", flush=True)


def quantile_tpm(score, tpm):
    """run_e2g_cv.calculate_quantiles, then 0 where TPM < TPM_MIN (filter_by_tpm)."""
    s = np.asarray(score, float)
    r = rankdata(np.where(s == 0, -np.inf, s), method="average")
    q = (r - 1) / (len(s) - 1)
    q[s == 0] = 0
    q[np.asarray(tpm) < TPM_MIN] = 0
    return q


# check the published score's conversion
q3 = quantile_tpm(P["E2G.Score"], P["RNA_pseudobulkTPM"])
rho = spearmanr(q3, P["E2G.Score.qnorm"])[0]
same0 = np.mean((q3 == 0) == (P["E2G.Score.qnorm"].values == 0))
assert rho > 0.999 and same0 > 0.999, f"quantile + TPM rule does not reproduce qnorm (rho {rho:.5f}, zeros {same0:.5f})"

# elements and their features
el = P[["chr", "start", "end"]].drop_duplicates().sort_values(["chr", "start"]).reset_index(drop=True)
for c, g in el.groupby("chr"):
    assert (g["start"].values[1:] >= g["end"].values[:-1]).all(), f"{c}: overlapping elements"
el["mid"] = (el["start"] + el["end"]) // 2
G = pd.read_csv(prep.CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])
H = pd.read_csv(f"{a.out_dir}/motif_sites_gata_ikzf2.tsv.gz", sep="\t")
import pysam  # noqa: E402
el = pd.concat([el, prep.element_features(el["chr"].values, el["mid"].values, G, H, pysam.FastaFile(prep.GENOME))], axis=1)
el["eid"] = np.arange(len(el))
P = P.merge(el[["chr", "start", "end", "eid", "cgi_p1", "gata_p1", "ikzf2_p1"]], on=["chr", "start", "end"], how="left", validate="m:1")

# arm scores on scE2G's quantile scale, plus the published score
Q = {"v3_published": q3}
n6_raw = None
for arm, m in models.items():
    X = np.log(np.abs(P[m["features"]].values.astype(float)) + m["eps"])
    s = 1 / (1 + np.exp(-(X @ np.array(m["coef"]) + m["intercept"])))
    if arm == "n6":
        n6_raw = s
    Q[arm] = quantile_tpm(s, P["RNA_pseudobulkTPM"])
rho_n6_v3 = spearmanr(n6_raw, P["E2G.Score"])[0]
uni = set(pd.read_csv(prep.UNIVERSE, sep="\t", header=None, usecols=[3])[3])
keep = P["TargetGene"].isin(uni).values
P = P[keep].reset_index(drop=True)
Q = {k: v[keep] for k, v in Q.items()}
names = list(Q)

# element quantile = max over its (universe) genes
emax = np.full((len(el), len(names)), 0.0)
for j, k in enumerate(names):
    np.maximum.at(emax[:, j], P["eid"].values, Q[k])

# background SNPs per element -> SNPs above each threshold, per arm
bg = np.load(f"{a.out_dir}/bg_snps.npz")
nbg = np.zeros(len(el), np.int64)
for c, g in el.groupby("chr"):
    if c not in bg.files:
        continue
    p = bg[c]
    nbg[g.index.values] = np.searchsorted(p, g["end"].values, side="left") - np.searchsorted(p, g["start"].values, side="left")
bgc = {k: [int(nbg[emax[:, j] >= t].sum()) for t in GRID] for j, k in enumerate(names)}
pd.DataFrame({"threshold": GRID, **bgc}).to_csv(f"{od}/{a.biosample}.bg.tsv", sep="\t", index=False)

# matched variants -> containing element -> element and link quantiles
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
inel = V["eid"].values >= 0
for j, k in enumerate(names):
    V[f"elem_{k}"] = np.where(inel, emax[np.maximum(V["eid"].values, 0), j], -1.0)
link = P[["eid", "TargetGene"]].assign(**{f"link_{k}": Q[k] for k in names})
link = link.groupby(["eid", "TargetGene"], as_index=False).max()
V = V.merge(link, left_on=["eid", "gene_hgnc"], right_on=["eid", "TargetGene"], how="left").drop(columns="TargetGene")
for k in names:
    V[f"link_{k}"] = V[f"link_{k}"].fillna(-1.0)
V["biosample"] = a.biosample
V.drop(columns=["eid"]).to_csv(f"{od}/{a.biosample}.variants.tsv.gz", sep="\t", index=False, float_format="%.6g")

stats = dict(biosample=a.biosample, n_pairs=int(len(keep)), n_pairs_universe=int(keep.sum()), n_elements=int(len(el)),
             qnorm_check_rho=float(rho), qnorm_check_zeros=float(same0), spearman_n6_refit_vs_v3=float(rho_n6_v3),
             n_variant_rows=int(len(V)), frac_variants_in_element=float(inel.mean()), n_bg_in_elements=int(nbg.sum()),
             feature_means={f: float(el[f].mean()) for f in ["cgi_log_dist", "gata_n", "ikzf2_n"]})
json.dump(stats, open(f"{od}/{a.biosample}.stats.json", "w"), indent=1)
print(json.dumps(stats), flush=True)
