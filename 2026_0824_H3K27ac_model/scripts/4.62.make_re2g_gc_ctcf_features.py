#!/usr/bin/env python3
"""Element GC and CTCF motif count as extra ENCODE-rE2G features: the pre-registered test of F-055.

WHY. F-055 found that element GC (+0.0139 AUPRC) and CTCF motif count (+0.0083) improve scE2G on
the K562 CRISPR benchmark, both with negative coefficients, but both were picked from F-054's residuals on the same pairs.
This asks whether the same two features, defined the same way, help a different model family
(ENCODE-rE2G, ENCODE bulk ATAC, the July candidate regions) that played no part in choosing them.
Decision 2026-10-05 in .living/decisions.md fixes the comparisons and the pass rule before the run.

DEFINITIONS, as 2026_1002_sce2g_tf_activity/scripts/2.4 and 4.7:
    element_gc     GC fraction of the element midpoint +/- 250 bp (ACGT bases only)
    ctcf_motif_n   JASPAR MA0139.2 sites in the element, TOBIAS TFBScan at p < 1e-4, overlaps
                   resolved to the best site (BINDetect's scan and defaults); a site is in the
                   element if its midpoint is

The CTCF site BED comes from 4.62.submit_re2g_gc_ctcf_features.sh.

WRITES
    data/re2g_external/k562_gc_ctcf_enh.tsv.gz   chr start end TargetGene element_gc ctcf_motif_n
                                                 (join_by overlap: per element-gene pair)
    $E2G/config/mine/external_features_K562_GCCTCF.tsv
    $E2G/config/mine/feature_tables/ATAC_{GC,CTCF,GCCTCF}_pl_all.tsv     (from ATAC_pl_all)
    $E2G/config/mine/feature_tables/ATAC_H3K27ac_GCCTCF_pl_all.tsv       (from ATAC_H3K27ac_pl_all)

Usage: 4.62.make_re2g_gc_ctcf_features.py --ctcf-sites ctcf_sites.bed
"""
import argparse
import os

import numpy as np
import pandas as pd
import pysam

OAK = "/oak/stanford/groups/engreitz/Users/sheth"
E2G = f"{OAK}/ENCODE_rE2G"
A = f"{OAK}/EP300_BPNet/2026_0824_H3K27ac_model"
SRC = f"{E2G}/results/2026_0923_p300_features/ATAC_powerlaw/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
CAND = (f"{OAK}/ABC_working/ABC-Enhancer-Gene-Prediction/results/2026_0721_h3k27ac_counting_comparison/"
        "K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
GENOME = f"{OAK}/hg38_resources/hg38.fa"
OUT = f"{A}/data/re2g_external"
ENH = f"{OUT}/k562_gc_ctcf_enh.tsv.gz"
EFC = f"{E2G}/config/mine/external_features_K562_GCCTCF.tsv"
FT_ATAC = f"{E2G}/resources/feature_tables/ATAC_pl_all.tsv"
FT_ACT = f"{E2G}/resources/feature_tables/ATAC_H3K27ac_pl_all.tsv"
FT_DIR = f"{E2G}/config/mine/feature_tables"
HALF = 250
KEY = ["chr", "start", "end"]

ap = argparse.ArgumentParser()
ap.add_argument("--ctcf-sites", required=True)
a = ap.parse_args()
os.makedirs(OUT, exist_ok=True)
os.makedirs(FT_DIR, exist_ok=True)

P = pd.read_csv(SRC, sep="\t", usecols=KEY + ["TargetGene"])
assert not P.duplicated(KEY + ["TargetGene"]).any()
U = P[KEY].drop_duplicates().sort_values(KEY).reset_index(drop=True)
print(f"{len(P):,} element-gene pairs, {len(U):,} elements from {SRC}", flush=True)

C = pd.read_csv(CAND, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
assert len(U.merge(C, on=KEY)) == len(U), "prediction elements missing from the July candidate set"
# rE2G's overlap join uses closed intervals, so book-ended elements would share each other's values
for c, g in U.groupby("chr"):
    assert (g["start"].values[1:] > g["end"].values[:-1]).all(), f"touching or overlapping elements on {c}"

fa = pysam.FastaFile(GENOME)
mid = ((U["start"] + U["end"]) // 2).values
gc = np.empty(len(U))
for i, (ch, m) in enumerate(zip(U["chr"].values, mid)):
    s = fa.fetch(ch, max(0, m - HALF), m + HALF).upper()
    gc[i] = (s.count("C") + s.count("G")) / max(sum(s.count(b) for b in "ACGT"), 1)
U["element_gc"] = gc

S = pd.read_csv(a.ctcf_sites, sep="\t", header=None, usecols=[0, 1, 2, 3], names=KEY + ["motif"])
assert S["motif"].nunique() == 1 and "MA0139.2" in S["motif"].iloc[0], S["motif"].unique()
S["m"] = (S["start"] + S["end"]) // 2
n = np.zeros(len(U), int)
for c, g in U.groupby("chr"):
    sm = np.sort(S.loc[S["chr"] == c, "m"].values)
    n[g.index.values] = (np.searchsorted(sm, g["end"].values, side="left")
                         - np.searchsorted(sm, g["start"].values, side="left"))
U["ctcf_motif_n"] = n
in_el = n.sum()
print(f"{len(S):,} CTCF sites, {in_el:,} with their midpoint in a paired element", flush=True)
# candidates with no gene in range are absent from the predictions, so a few sites go unused
assert in_el > 0.95 * len(S), "CTCF sites outside the elements: wrong region set or coordinates"
print(U[["element_gc", "ctcf_motif_n"]].describe().round(3).to_string())
print(f"CTCF site in {(U['ctcf_motif_n'] > 0).mean():.3f} of elements", flush=True)

M = P.merge(U, on=KEY, how="left")
assert len(M) == len(P) and M[["element_gc", "ctcf_motif_n"]].notna().all().all()
M.to_csv(ENH, sep="\t", index=False, float_format="%.6g")
print(f"wrote {ENH} ({len(M):,} rows)")

# merge_external_features.R filters this table on input_col and renames source_col to it
pd.DataFrame([
    ("element_gc", "element_gc", "mean", "overlap", ENH),
    ("ctcf_motif_n", "ctcf_motif_n", "max", "overlap", ENH),
], columns=["input_col", "source_col", "aggregate_function", "join_by", "source_file"]).to_csv(
    EFC, sep="\t", index=False)
print(f"wrote {EFC}")

ROWS = {"GC": ("elementGC", "element_gc", "NA", "mean", "0", "GC fraction of E"),
        "CTCF": ("ctcfMotifN", "ctcf_motif_n", "NA", "max", "0", "CTCF motif sites in E")}
for name, src, keys in (("ATAC_GC", FT_ATAC, ["GC"]), ("ATAC_CTCF", FT_ATAC, ["CTCF"]),
                        ("ATAC_GCCTCF", FT_ATAC, ["GC", "CTCF"]),
                        ("ATAC_H3K27ac_GCCTCF", FT_ACT, ["GC", "CTCF"])):
    ft = pd.read_csv(src, sep="\t", keep_default_na=False)
    extra = pd.DataFrame([ROWS[k] for k in keys], columns=ft.columns)
    assert not set(extra["feature"]) & set(ft["feature"])
    # keep ABC.Score last, as in the shipped tables
    ft = pd.concat([ft[ft["feature"] != "ABC.Score"], extra, ft[ft["feature"] == "ABC.Score"]])
    path = f"{FT_DIR}/{name}_pl_all.tsv"
    ft.to_csv(path, sep="\t", index=False)
    print(f"wrote {path} ({len(ft)} features)")
