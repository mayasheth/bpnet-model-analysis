#!/usr/bin/env python3
"""9.3 Add a published scE2G prediction file as an extra arm of the plan 9 eQTL benchmark.

Positive control for F-063 (Maya, 2026-10-13): does the benchmark see a known CRISPR-scale change?
scE2G scATAC (scATAC_powerlaw_v3) vs multiome predictions differ by about -0.05 AUPRC on K562
CRISPR. This reads one biosample's file for the given method column of config_predictions_hd_fine
(e.g. scE2G_ATAC), puts E2G.Score on scE2G's quantile scale (run_e2g_cv.calculate_quantiles, with
the TPM >= 1 rule only if that reproduces the file's own E2G.Score.qnorm; the rule used is
recorded), and writes the same reductions as 9.1 for that arm on its own elements. 9.2 merges them.

OUTPUT ($OUT/clusters/<biosample>.<arm>.variants.tsv.gz, .<arm>.bg.tsv, .<arm>.stats.json)

Usage: 9.3.eqtl_transfer_published_arm.py --biosample PBMC_25_NK --method scE2G_ATAC --arm scatac_published
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

SC = os.environ["SCRATCH"]
EQ = "/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment"
UNIVERSE = f"{EQ}/resources/genome_annotation/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6"
MATCH = f"{EQ}/config/mine/config_predictions_hd_fine.tsv"
GRID = np.unique(np.r_[np.linspace(0, 0.9, 91), 1 - np.logspace(-1, -5, 400)])  # as 9.1
ap = argparse.ArgumentParser()
ap.add_argument("--biosample", required=True)
ap.add_argument("--method", default="scE2G_ATAC")
ap.add_argument("--arm", default="scatac_published")
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
a = ap.parse_args()
od = f"{a.out_dir}/clusters"

cfg = pd.read_csv(MATCH, sep="\t").set_index("biosample")
path = cfg.loc[a.biosample, a.method]
mt = pd.read_csv(f"{a.out_dir}/matches.tsv", sep="\t")
mt = mt[mt["biosample"] == a.biosample]
P = pd.read_csv(path, sep="\t", usecols=["chr", "start", "end", "TargetGene", "RNA_pseudobulkTPM", "E2G.Score", "E2G.Score.qnorm"])
print(f"{a.biosample} {a.method}: {len(P):,} pairs from {path}", flush=True)


def quantiles(s):
    s = np.asarray(s, float)
    r = rankdata(np.where(s == 0, -np.inf, s), method="average")
    q = (r - 1) / (len(s) - 1)
    q[s == 0] = 0
    return q


q = quantiles(P["E2G.Score"])
checks = {}
for rule in ["tpm", "none"]:
    qq = np.where(P["RNA_pseudobulkTPM"].values < 1.0, 0, q) if rule == "tpm" else q
    checks[rule] = (float(spearmanr(qq, P["E2G.Score.qnorm"])[0]), float(np.mean((qq == 0) == (P["E2G.Score.qnorm"].values == 0))))
ok = [r for r, (rho, z) in checks.items() if rho > 0.999 and z > 0.999]
assert ok, f"no quantile rule reproduces qnorm: {checks}"
rule = ok[0]
if rule == "tpm":
    q = np.where(P["RNA_pseudobulkTPM"].values < 1.0, 0, q)

uni = set(pd.read_csv(UNIVERSE, sep="\t", header=None, usecols=[3])[3])
keep = P["TargetGene"].isin(uni).values
P, q = P[keep].reset_index(drop=True), q[keep]
el = P[["chr", "start", "end"]].drop_duplicates().sort_values(["chr", "start"]).reset_index(drop=True)
for c, g in el.groupby("chr"):
    assert (g["start"].values[1:] >= g["end"].values[:-1]).all(), f"{c}: overlapping elements"
el["eid"] = np.arange(len(el))
P = P.merge(el, on=["chr", "start", "end"], how="left", validate="m:1")
emax = np.zeros(len(el))
np.maximum.at(emax, P["eid"].values, q)

bg = np.load(f"{a.out_dir}/bg_snps.npz")
nbg = np.zeros(len(el), np.int64)
for c, g in el.groupby("chr"):
    if c in bg.files:
        p = bg[c]
        nbg[g.index.values] = np.searchsorted(p, g["end"].values, side="left") - np.searchsorted(p, g["start"].values, side="left")
pd.DataFrame({"threshold": GRID, a.arm: [int(nbg[emax >= t].sum()) for t in GRID]}).to_csv(f"{od}/{a.biosample}.{a.arm}.bg.tsv", sep="\t", index=False)

V = pd.read_csv(f"{a.out_dir}/variants.tsv.gz", sep="\t")
V = V[V["tissue"].isin(set(mt["tissue"]))].reset_index(drop=True)
V["eid"] = -1
for c, g in V.groupby("chr"):
    e = el[el["chr"] == c]
    if e.empty:
        continue
    k = np.searchsorted(e["start"].values, g["start"].values, side="right") - 1
    hit = (k >= 0) & (e["end"].values[np.maximum(k, 0)] > g["start"].values)
    V.loc[g.index.values[hit], "eid"] = e["eid"].values[k[hit]]
V[f"elem_{a.arm}"] = np.where(V["eid"].values >= 0, emax[np.maximum(V["eid"].values, 0)], -1.0)
link = P[["eid", "TargetGene"]].assign(lq=q).groupby(["eid", "TargetGene"], as_index=False)["lq"].max()
V = V.merge(link, left_on=["eid", "gene_hgnc"], right_on=["eid", "TargetGene"], how="left")
V[f"link_{a.arm}"] = V["lq"].fillna(-1.0)
V["biosample"] = a.biosample
V[["chr", "start", "end", "varID_hg38", "gene_hgnc", "tissue", "pip", "biosample", f"elem_{a.arm}", f"link_{a.arm}"]].to_csv(
    f"{od}/{a.biosample}.{a.arm}.variants.tsv.gz", sep="\t", index=False, float_format="%.6g")
st = dict(biosample=a.biosample, method=a.method, path=path, rule=rule, checks=checks, n_pairs=int(len(keep)),
          n_elements=int(len(el)), frac_variants_in_element=float((V["eid"] >= 0).mean()))
json.dump(st, open(f"{od}/{a.biosample}.{a.arm}.stats.json", "w"), indent=1)
print(json.dumps(st), flush=True)
