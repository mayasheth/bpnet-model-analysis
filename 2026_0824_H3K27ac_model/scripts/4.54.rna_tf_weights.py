#!/usr/bin/env python3
"""TF lists from RNA: weight each motif by how much of its TF is expressed in the cell.

WHY. The TF-activity track's gain is small and the ATAC-derived TF list only partly recovers
the cell-specific signal (F-034, F-038). DeepSCENIC (Partel et al. 2026, bioRxiv
2026.09.18.752607) models a TF's activity at an element as (sequence binding potential) x
(the TF's expression in the cell), i.e. the same sum our track computes, with RNA supplying the
cell-specific part. In the intended application (single-cell multiome) RNA is available in the
target cell, so an RNA-based TF list is deployable there.

EXPRESSION. ENCODE poly(A) RNA-seq, Gingeras lab, one untreated whole-cell experiment per cell
type (both replicates, RSEM gene quantifications, GENCODE V29), mean TPM, log2(TPM + 1). Gene IDs
mapped to symbols through the local GENCODE v43 GTF (versions stripped).

MOTIF -> TF. JASPAR names, upper-cased; a heterodimer (A::B) takes the smaller of its partners'
expression, since both are needed. JASPAR has one matrix per TF, but paralogs bind the same
site (GATA5's matrix is bound by GATA1 in K562, where GATA5 is not expressed), so the FAMILY
version takes, for each motif, the highest expression among the TFs of all motifs in its
co-occurrence family (motifs whose FIMO hits co-occur at Pearson r >= 0.5 across the union
candidate set, average linkage, as in 4.49). Unmapped motifs get NaN (weight 0 in a track).

WEIGHTINGS written (per cell, standardised within cell across scored motifs):
    rna_<cell>         z(log2 expression), motif's own TF
    rnafam_<cell>      z(log2 expression), family maximum
    rnaatac_<cell>     z(z(ATAC association, CTCF removed) + z(rnafam))
    rnagate_<cell>     ATAC association (CTCF removed) where the family is expressed
                       (family max TPM >= 1), else 0
Scored as 4.50 does against 4.43's truth (Spearman, cell-specific part, K562 GATA and top 10).
These motif-level scores are proxies and do not predict the benchmark reliably (F-030, F-036).

Usage: 4.54.rna_tf_weights.py [--rna-dir $SCRATCH/rnaseq] [--out-prefix results/2026_0930_rna_tf_weights]
"""
import argparse
import glob
import gzip
import os
import re

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
GTF = f"{D}/scE2G_temp/scE2G/resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
REC = f"{A}/results/2026_0926_tf_activity_recovery.motifs.tsv"
CACHE = f"{A}/data/tf_activity_recovery"
RNA = {"K562": ("ENCSR000AEM", ["ENCFF742CVV", "ENCFF222UVT"]),
       "GM12878": ("ENCSR000AED", ["ENCFF791MED", "ENCFF473KMX"]),
       "HepG2": ("ENCSR000CPE", ["ENCFF831QQF", "ENCFF168QKW"]),
       "A549": ("ENCSR000CON", ["ENCFF855AKQ", "ENCFF244DNJ"]),
       "MCF-7": ("ENCSR000CPT", ["ENCFF967AOT", "ENCFF921PJP"])}
CELLS = list(RNA)

ap = argparse.ArgumentParser()
ap.add_argument("--rna-dir", default=os.path.join(os.environ.get("SCRATCH", "/tmp"), "rnaseq"))
ap.add_argument("--family-r", type=float, default=0.5)
ap.add_argument("--out-prefix", default="results/2026_0930_rna_tf_weights")
a = ap.parse_args()

# ------------------------------------------------------------------ gene id -> symbol
sym = {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t", 9)
        if p[2] != "gene":
            continue
        gid = re.search(r'gene_id "([^"]+)"', p[8]).group(1).split(".")[0]
        sym[gid] = re.search(r'gene_name "([^"]+)"', p[8]).group(1).upper()
print(f"{len(sym):,} GENCODE v43 genes")

EXPR = {}
for c, (exp, accs) in RNA.items():
    tpm = []
    for acc in accs:
        t = pd.read_csv(f"{a.rna_dir}/{acc}.tsv", sep="\t", usecols=["gene_id", "TPM"])
        t["sym"] = t["gene_id"].str.split(".").str[0].map(sym)
        tpm.append(t.dropna(subset=["sym"]).groupby("sym")["TPM"].sum())
    e = pd.concat(tpm, axis=1).mean(1)
    EXPR[c] = e
    print(f"{c}: {exp}, {len(accs)} replicates, {int((e >= 1).sum()):,} genes TPM >= 1; "
          f"GATA1 {e.get('GATA1', np.nan):.1f}, TAL1 {e.get('TAL1', np.nan):.1f}, "
          f"HNF4A {e.get('HNF4A', np.nan):.1f}, ESR1 {e.get('ESR1', np.nan):.1f}, "
          f"CEBPB {e.get('CEBPB', np.nan):.1f} TPM", flush=True)

Rm = pd.read_csv(REC, sep="\t")
mnames = list(Rm["motif"])
tfs = [m.split(" ", 1)[1].upper() for m in mnames]
parts = [re.split(r"::", t) for t in tfs]


def motif_expr(e):
    out = np.full(len(mnames), np.nan)
    for i, ps in enumerate(parts):
        v = [e.get(p_) for p_ in ps]
        if all(x is not None for x in v):
            out[i] = min(v)
    return out


ME = {c: motif_expr(EXPR[c]) for c in CELLS}
unm = [tfs[i] for i in range(len(mnames)) if np.isnan(ME["K562"][i])]
print(f"motifs mapped to expressed-gene symbols: {len(mnames) - len(unm)} of {len(mnames)}; "
      f"unmapped e.g. {', '.join(unm[:12])}")

# ------------------------------------------------------------------ co-occurrence families
z = np.load(f"{CACHE}/union_fimo_presence.npz", allow_pickle=True)
M_all, names_all = z["presence"], list(z["names"])
M = M_all[:, [names_all.index(n) for n in mnames]]
rng = np.random.default_rng(0)
samp = M[rng.choice(len(M), min(300000, len(M)), replace=False)].astype(np.float32)
use = samp.sum(0) >= 50
r = np.nan_to_num(np.corrcoef(samp[:, use].T))
dist = np.clip(1 - r, 0, 2)
np.fill_diagonal(dist, 0)
lab = np.full(len(mnames), -1)
lab[use] = fcluster(linkage(squareform(dist, checks=False), "average"), 1 - a.family_r, "distance")
lab[~use] = np.arange((~use).sum()) + lab.max() + 1
print(f"{len(np.unique(lab))} co-occurrence families")


def family_max(v):
    s = pd.Series(v).groupby(lab).transform("max")
    return s.values


zs = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
ctcf = np.array([t.startswith("CTCF") for t in tfs])
W = {}
for c in CELLS:
    le = np.log2(ME[c] + 1)
    fam = family_max(le)
    atac = Rm[f"within_{c}"].values.astype(float).copy()
    atac[ctcf] = np.nan
    W[f"rna_{c}"] = zs(le)
    W[f"rnafam_{c}"] = zs(fam)
    W[f"rnaatac_{c}"] = zs(zs(atac) + zs(fam))
    g = zs(atac)
    W[f"rnagate_{c}"] = np.where(family_max(ME[c]) >= 1, g, 0.0)
    W[f"atac_{c}"] = zs(atac)
    for k in ("rna", "rnafam", "rnaatac", "rnagate"):
        W[f"{k}_{c}"][ctcf] = np.nan
T = {c: zs(Rm[f"truth_{c}"].values) for c in CELLS}

# ------------------------------------------------------------------ proxy scores
top = lambda v, k=20: set(np.argsort(-np.nan_to_num(v, nan=-9))[:k])
rows = []
for k in ("atac", "rna", "rnafam", "rnaatac", "rnagate"):
    for c in CELLS:
        e, t = np.nan_to_num(W[f"{k}_{c}"], nan=0.0), T[c]
        prior = np.mean([T[d] for d in CELLS if d != c], 0)
        rows.append(dict(estimate=k, cell=c, spearman=spearmanr(e, t)[0],
                         dep_spearman=spearmanr(e - prior, t - prior)[0], top20=len(top(e) & top(t))))
S = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== mean over the five cells (missing weights scored as 0) ===")
print(S.groupby("estimate", sort=False)[["spearman", "dep_spearman", "top20"]].mean().round(3).to_string())
print(S.pivot(index="estimate", columns="cell", values="dep_spearman")[CELLS].round(3).to_string())
gi = [i for i, t in enumerate(tfs) if "GATA" in t or "TRPS1" in t]
print(f"\nK562 truth       GATA {T['K562'][gi].mean():+.2f}  top10: "
      + ", ".join(tfs[i] for i in np.argsort(-T["K562"])[:10]))
for k in ("atac", "rna", "rnafam", "rnaatac", "rnagate"):
    e = np.nan_to_num(W[f"{k}_K562"], nan=0.0)
    print(f"K562 {k:<11} GATA {e[gi].mean():+.2f}  top10: "
          + ", ".join(tfs[i] for i in np.argsort(-e)[:10]))

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
out = pd.DataFrame({"motif": mnames, "tf": tfs, "family": lab})
for c in CELLS:
    out[f"log2tpm_{c}"] = np.log2(ME[c] + 1)
for k, v in W.items():
    out[k] = v
out.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
S.to_csv(f"{a.out_prefix}.scores.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.motifs.tsv, .scores.tsv")
