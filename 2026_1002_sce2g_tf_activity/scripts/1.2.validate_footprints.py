#!/usr/bin/env python3
"""Do the 1.1 footprint scores detect TF binding beyond accessibility?

TEST 1, against ChIP. For CTCF motif hits (K562 CTCF ChIP, ENCSR000AKO optimal IDR) and GATA1 /
GATA1::TAL1 hits (K562 GATA1 ChIP peaks used by K562_GATA1_BPNet): label a hit bound if it
overlaps a ChIP peak, then compute the AUROC of FS for bound vs unbound WITHIN deciles of the
hit's insertion count, so accessibility cannot carry it. Reported beside the AUROC of insertion
count alone. An FS AUROC near 0.5 within strata means the score is not a footprint.

TEST 2, motif ranking. Mean FS per motif over well-covered hits (n_ins >= 50). K562 regulators
(CTCF, GATA, NFE2/MAF, TAL1) should rank above liver and gut TFs that K562 does not express
(HNF4A, HNF1A, FOXA).

Usage: 1.2.validate_footprints.py [--out-prefix results/2026_1002_footprint_validation]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
FP = os.path.join(os.environ["SCRATCH"], "footprints_k562_xu")
CHIP = {"CTCF": f"{D}/Data/ENCODE/K562/ENCFF519CXF.bed.gz",
        "GATA1": f"{D}/EP300_BPNet/K562_GATA1_BPNet/data/peaks_inliers.bed.gz"}
TESTS = [("MA0139.2 CTCF", "CTCF"), ("MA0035.5 GATA1", "GATA1"), ("MA0140.3 GATA1::TAL1", "GATA1")]
CONTROLS = ["MA0139.2 CTCF", "MA0035.5 GATA1", "MA0140.3 GATA1::TAL1", "MA0841.2 NFE2",
            "MA0501.2 MAF::NFE2", "MA0091.2 TAL1::TCF3", "MA0080.7 Spi1", "MA0114.5 HNF4A",
            "MA0046.3 HNF1A", "MA0148.5 FOXA1", "MA0047.4 FOXA2"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_1002_footprint_validation")
a = ap.parse_args()

names = pd.read_csv(f"{FP}/motif_names.tsv", sep="\t", header=None, index_col=0)[1].tolist()
idx = {n: i for i, n in enumerate(names)}
parts = []
for f in sorted(glob.glob(f"{FP}/fs_chr*.npz")):
    z = np.load(f)
    parts.append(pd.DataFrame({"chr": os.path.basename(f)[3:-4], **{k: z[k] for k in z.files}}))
H = pd.concat(parts, ignore_index=True)
H["name"] = np.array(names, dtype=object)[H["motif"].values]
print(f"{len(H):,} accessible motif hits on {H['chr'].nunique()} chromosomes", flush=True)


def overlaps(reg, iv):
    hit = np.zeros(len(reg), bool)
    for c, g in reg.groupby("chr"):
        e = iv[iv["chr"] == c].sort_values("start")
        if e.empty:
            continue
        es, ee = e["start"].values, np.maximum.accumulate(e["end"].values)
        k = np.searchsorted(es, g["end"].values, side="left") - 1
        ok = k >= 0
        hit[g.index.values[ok]] = ee[k[ok]] > g["start"].values[ok]
    return hit


def auroc(score, y):
    if y.sum() == 0 or (~y).sum() == 0:
        return np.nan
    r = rankdata(score)
    return (r[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum())


rows = []
for motif, tf in TESTS:
    g = H[H["name"] == motif].reset_index(drop=True)
    pk = pd.read_csv(CHIP[tf], sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    y = overlaps(g, pk)
    dec = pd.qcut(rankdata(g["n_ins"].values, method="ordinal"), 10, labels=False)
    w, auc_fs = [], []
    for d in range(10):
        m = dec == d
        auc_fs.append(auroc(g["fs"].values[m], y[m]))
        w.append(y[m].sum())
    auc_fs, w = np.array(auc_fs), np.array(w, float)
    ok = np.isfinite(auc_fs)
    rows.append(dict(motif=motif, chip=tf, n_hits=len(g), n_bound=int(y.sum()),
                     auroc_ins=auroc(g["n_ins"].values, y),
                     auroc_fs_within_ins_decile=np.average(auc_fs[ok], weights=w[ok]),
                     auroc_fs_top_decile=auc_fs[9]))
T1 = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== test 1: FS vs ChIP binding, within accessibility deciles ===")
print(T1.round(3).to_string(index=False))

cov = H[H["n_ins"] >= 50]
T2 = cov.groupby("name")["fs"].agg(["mean", "size"]).rename(columns={"size": "n"})
T2 = T2[T2["n"] >= 200].sort_values("mean", ascending=False)
T2["rank"] = np.arange(1, len(T2) + 1)
print(f"\n=== test 2: mean FS per motif, hits with n_ins >= 50 ({len(T2)} motifs with >= 200 hits) ===")
print("top 20:")
print(T2.head(20).round(3).to_string())
print("controls:")
print(T2.reindex([c for c in CONTROLS if c in T2.index]).round(3).to_string())

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
T1.to_csv(f"{a.out_prefix}.chip.tsv", sep="\t", index=False)
T2.to_csv(f"{a.out_prefix}.motif_rank.tsv", sep="\t")
print(f"\nwrote {a.out_prefix}.chip.tsv, .motif_rank.tsv")
