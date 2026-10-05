#!/usr/bin/env python3
"""Do TOBIAS footprint scores (2.2) detect TF binding beyond accessibility? Same tests as 1.2 (F-047).

TEST 1, against ChIP. For BINDetect sites of CTCF (K562 CTCF ChIP, ENCSR000AKO optimal IDR) and
GATA1 / GATA1::TAL1 (K562 GATA1 ChIP peaks used by K562_GATA1_BPNet): label a site bound if it
overlaps a ChIP peak, then compute the AUROC of the TOBIAS footprint score (`K562_Xu_score`) for
bound vs unbound WITHIN deciles of the site's accessibility (uncorrected insertions in the site
centre +/- 100 bp), beside the AUROC of accessibility alone. Near 0.5 within strata = no footprint.

TEST 2, motif ranking. BINDetect's mean footprint score per motif: K562 regulators (CTCF, GATA,
NFE2/MAF, TAL1) should rank above liver and gut TFs K562 does not express (HNF4A, HNF1A, FOXA).

Run in the tobias pixi env. Usage: 2.3.validate_tobias.py [--out-prefix results/2026_1002_tobias_validation]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
TB = os.path.join(os.environ["SCRATCH"], "tobias_k562_xu")
CHIP = {"CTCF": f"{D}/Data/ENCODE/K562/ENCFF519CXF.bed.gz",
        "GATA1": f"{D}/EP300_BPNet/K562_GATA1_BPNet/data/peaks_inliers.bed.gz"}
TESTS = [("MA0139.2", "CTCF"), ("MA0035.5", "GATA1"), ("MA0140.3", "GATA1")]
CONTROLS = ["MA0139.2", "MA0035.5", "MA0140.3", "MA0841.2", "MA0501.2", "MA0091.2", "MA0080.7",
            "MA0114.5", "MA0046.3", "MA0148.5", "MA0047.4"]
FLANK = 100

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_1002_tobias_validation")
ap.add_argument("--score", choices=["tobias", "print"], default="tobias",
                help="print: PRINT (3.1) scores at the same BINDetect sites; test 2 is TOBIAS-only")
a = ap.parse_args()
if a.score == "print":
    a.out_prefix = a.out_prefix.replace("tobias_validation", "print_validation")
from footprint_scores import site_scores  # noqa: E402
MID2NAME = {"MA0139.2": "CTCF", "MA0035.5": "GATA1", "MA0140.3": "GATA1::TAL1"}
assert os.path.exists(f"{TB}/tobias_complete.txt"), "TOBIAS 2.2 not complete"


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


bw = pyBigWig.open(f"{TB}/atacorrect/k562_xu_uncorrected.bw")
rows = []
for mid, tf in TESTS:
    (ov,) = glob.glob(f"{TB}/bindetect/*_{mid}/*_{mid}_overview.txt")
    g = site_scores(pd.read_csv(ov, sep="\t"), a.score, MID2NAME[mid])
    g = g.reset_index(drop=True).rename(columns={"TFBS_chr": "chr", "TFBS_start": "start", "TFBS_end": "end"})
    mid_pt = (g["start"].values + g["end"].values) // 2
    g["n_ins"] = [np.nansum(bw.values(c, max(0, m - FLANK), m + FLANK, numpy=True))
                  for c, m in zip(g["chr"].values, mid_pt)]
    pk = pd.read_csv(CHIP[tf], sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    y = overlaps(g, pk)
    s = g["K562_Xu_score"].values
    dec = pd.qcut(rankdata(g["n_ins"].values, method="ordinal"), 10, labels=False)
    auc, w = [], []
    for d in range(10):
        m = dec == d
        auc.append(auroc(s[m], y[m]))
        w.append(y[m].sum())
    auc, w = np.array(auc), np.array(w, float)
    ok = np.isfinite(auc)
    # 95% interval: resample sites within each decile, 200 reps
    rng = np.random.default_rng(0)
    boot = []
    for _ in range(200):
        bs, bw_ = [], []
        for d in range(10):
            i = rng.choice(np.flatnonzero(dec == d), (dec == d).sum())
            bs.append(auroc(s[i], y[i]))
            bw_.append(y[i].sum())
        bs, bw_ = np.array(bs), np.array(bw_, float)
        k = np.isfinite(bs)
        boot.append(np.average(bs[k], weights=bw_[k]))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    rows.append(dict(motif=os.path.basename(ov)[:-13], chip=tf, n_sites=len(g), n_bound=int(y.sum()),
                     auroc_ins=auroc(g["n_ins"].values, y), auroc_score=auroc(s, y),
                     auroc_score_within_ins_decile=np.average(auc[ok], weights=w[ok]),
                     ci_lo=lo, ci_hi=hi,
                     auroc_score_top_decile=auc[9],
                     frac_chip_in_tobias_bound=g.loc[y, "K562_Xu_bound"].mean(),
                     frac_unbound_in_tobias_bound=g.loc[~y, "K562_Xu_bound"].mean()))
T1 = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print("=== test 1: TOBIAS footprint score vs ChIP binding, within accessibility deciles ===")
print(T1.round(3).to_string(index=False))
f047 = "results/2026_1002_footprint_validation.chip.tsv"
if os.path.exists(f047):
    print("\nfor comparison, homemade footprints (1.2, F-047):")
    print(pd.read_csv(f047, sep="\t").round(3).to_string(index=False))

if a.score == "print":
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    T1.to_csv(f"{a.out_prefix}.chip.tsv", sep="\t", index=False)
    print(f"\nwrote {a.out_prefix}.chip.tsv (test 2 skipped: BINDetect ranking is TOBIAS-only)")
    raise SystemExit(0)
R = pd.read_csv(f"{TB}/bindetect/bindetect_results.txt", sep="\t")
R = R.sort_values("K562_Xu_mean_score", ascending=False).reset_index(drop=True)
R["rank"] = np.arange(1, len(R) + 1)
cols = ["name", "motif_id", "total_tfbs", "K562_Xu_mean_score", "K562_Xu_bound", "rank"]
print(f"\n=== test 2: BINDetect mean footprint score per motif ({len(R)} motifs) ===")
print("top 20:")
print(R[cols].head(20).round(3).to_string(index=False))
print("controls:")
print(R[R["motif_id"].isin(CONTROLS)][cols].round(3).to_string(index=False))

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
T1.to_csv(f"{a.out_prefix}.chip.tsv", sep="\t", index=False)
R[cols].to_csv(f"{a.out_prefix}.motif_rank.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.chip.tsv, .motif_rank.tsv")
