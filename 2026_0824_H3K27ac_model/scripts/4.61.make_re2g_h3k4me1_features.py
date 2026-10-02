#!/usr/bin/env python3
"""Observed K562 H3K4me1 at the enhancer and at the promoter as two extra rE2G features.

WHY. The rE2G arm ATAC_H3K4me1_powerlaw (4.32) puts H3K4me1 in the H3K27ac slot, so it asks
whether H3K4me1 can REPLACE H3K27ac. This builds the inputs for ATAC_H3K4me1sep_powerlaw: the
ten ATAC-only features (ATAC_pl_all) plus H3K4me1 at E and H3K4me1 at P as separate features,
added through rE2G's external_features_config. No combined-activity or ABC term uses H3K4me1.

VALUES. Taken from the swap arm's own ABC output (normalized_h3k27ac_enh / _prom computed from
the H3K4me1 BAMs), so the counting and normalisation are exactly ABC's, over the same staged
July candidate regions that every 4.32 dataset uses.

WRITES
    data/re2g_external/k562_h3k4me1_enh.tsv.gz    chr start end TargetGene normalized_h3k4me1_enh
                                                  (join_by overlap: per element-gene pair)
    data/re2g_external/k562_h3k4me1_prom.tsv      TargetGene normalized_h3k4me1_prom
                                                  (join_by TargetGene: one row per gene)
    $E2G/config/mine/external_features_K562_H3K4me1.tsv
    $E2G/config/mine/feature_tables/ATAC_H3K4me1sep_pl_all.tsv

Usage: 4.61.make_re2g_h3k4me1_features.py
"""
import os

import numpy as np
import pandas as pd

OAK = "/oak/stanford/groups/engreitz/Users/sheth"
E2G = f"{OAK}/ENCODE_rE2G"
A = f"{OAK}/EP300_BPNet/2026_0824_H3K27ac_model"
SRC = f"{E2G}/results/2026_0923_p300_features/ATAC_H3K4me1_powerlaw/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
OUT = f"{A}/data/re2g_external"
ENH = f"{OUT}/k562_h3k4me1_enh.tsv.gz"
PROM = f"{OUT}/k562_h3k4me1_prom.tsv"
EFC = f"{E2G}/config/mine/external_features_K562_H3K4me1.tsv"
FT_IN = f"{E2G}/resources/feature_tables/ATAC_pl_all.tsv"
FT_OUT = f"{E2G}/config/mine/feature_tables/ATAC_H3K4me1sep_pl_all.tsv"

os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.dirname(FT_OUT), exist_ok=True)

P = pd.read_csv(SRC, sep="\t", usecols=["chr", "start", "end", "TargetGene",
                                        "normalized_h3k27ac_enh", "normalized_h3k27ac_prom"])
print(f"{len(P):,} element-gene pairs, {P['TargetGene'].nunique():,} genes from {SRC}", flush=True)
assert P[["normalized_h3k27ac_enh", "normalized_h3k27ac_prom"]].notna().all().all()

enh = P[["chr", "start", "end", "TargetGene", "normalized_h3k27ac_enh"]].rename(
    columns={"normalized_h3k27ac_enh": "normalized_h3k4me1_enh"})
assert not enh.duplicated(["chr", "start", "end", "TargetGene"]).any()
# the enhancer value belongs to the element, so it must not vary across its genes
spread = enh.groupby(["chr", "start", "end"])["normalized_h3k4me1_enh"].agg(np.ptp)
assert (spread < 1e-9).all(), f"enhancer value varies by gene at {(spread >= 1e-9).sum()} elements"
enh.to_csv(ENH, sep="\t", index=False)

g = P.groupby("TargetGene")["normalized_h3k27ac_prom"]
spread = g.agg(np.ptp)
assert (spread < 1e-9).all(), f"promoter value varies within {(spread >= 1e-9).sum()} genes"
prom = g.first().rename("normalized_h3k4me1_prom").reset_index()
prom.to_csv(PROM, sep="\t", index=False)
print(f"wrote {ENH} ({len(enh):,} rows), {PROM} ({len(prom):,} genes)")
print("H3K4me1 at E quantiles:", np.round(enh["normalized_h3k4me1_enh"].quantile([.1, .5, .9]).values, 3))
print("H3K4me1 at P quantiles:", np.round(prom["normalized_h3k4me1_prom"].quantile([.1, .5, .9]).values, 3))

# merge_external_features.R filters this table on input_col and renames source_col to it
pd.DataFrame([
    ("normalized_h3k4me1_enh", "normalized_h3k4me1_enh", "mean", "overlap", ENH),
    ("normalized_h3k4me1_prom", "normalized_h3k4me1_prom", "mean", "TargetGene", PROM),
], columns=["input_col", "source_col", "aggregate_function", "join_by", "source_file"]).to_csv(
    EFC, sep="\t", index=False)

ft = pd.read_csv(FT_IN, sep="\t", keep_default_na=False)
extra = pd.DataFrame([
    ("normalizedH3K4me1_enh", "normalized_h3k4me1_enh", "NA", "mean", "0", "H3K4me1 signal at E"),
    ("normalizedH3K4me1_prom", "normalized_h3k4me1_prom", "NA", "mean", "0", "H3K4me1 signal at P"),
], columns=ft.columns)
# keep ABC.Score last, as in the shipped tables
ft = pd.concat([ft[ft["feature"] != "ABC.Score"], extra, ft[ft["feature"] == "ABC.Score"]])
ft.to_csv(FT_OUT, sep="\t", index=False)
print(f"wrote {EFC}\nwrote {FT_OUT} ({len(ft)} features)")
