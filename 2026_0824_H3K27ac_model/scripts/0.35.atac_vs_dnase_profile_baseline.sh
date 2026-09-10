#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 3:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/atac_dnase_shape.%j.txt
#SBATCH -e log/atac_dnase_shape.%j.txt
#SBATCH --job-name=atac_dnase_shape
#
# THE COMPARATOR THE CONVERTER ACTUALLY HAS TO BEAT, WHICH THE PILOT DID NOT MEASURE.
#
# 2.31 scores the converter against the DNase inter-replicate ceiling, and the fold-0
# numbers look strong: 96.2% of ceiling in K562 at 1 bp on the top quintile, 71.8% on
# transfer to GM12878. But "fraction of the DNase ceiling" is the wrong denominator for
# the decision being made. The converter exists to REPLACE the ATAC track as the
# downstream model's accessibility input. So the question is not how close it gets to
# DNase, it is whether it is MORE DNase-like than the ATAC track already sitting there.
#
# ATAC and DNase are both cut-site assays over the same open chromatin. If the raw
# observed ATAC profile already correlates with the observed DNase profile at close to
# what the converter achieves, the converter is an expensive no-op and the downstream
# model would see nothing new.
#
# NOTE THE ATAC-ONLY MODEL ARM IS NOT THIS CONTROL. That arm is a trained ATAC->DNase
# mapping and has already learned a reweighting. This is the untouched input track.
#
# Model-free, so it runs on CPU and takes minutes.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

$PY - <<'PY'
import json
import numpy as np
import pandas as pd
import pyBigWig
import sys

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
sys.path.insert(0, f"{D}/scripts")
from train_multimodal_bpnet import load_peaks

FOLDS = f"{D}/reference/hg38_five_folds.json"
HW = 500
OUT_W = 2 * HW
BINS = [1, 5, 10, 25, 50, 100, 250]

CELLS = {
    "k562": dict(
        atac=f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
        dnase_p=f"{P}/data/k562_dnase_5p_plus.bw",
        dnase_m=f"{P}/data/k562_dnase_5p_minus.bw",
        elements=f"{D}/reference/K562_ATAC_candidate_elements.narrowPeak"),
    "gm12878": dict(
        atac=f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw",
        dnase_p=f"{P}/data/gm12878_dnase_5p_plus.bw",
        dnase_m=f"{P}/data/gm12878_dnase_5p_minus.bw",
        elements=f"{D}/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"),
}


def shape_corr(x, y, b):
    """Identical to 0.25 and 2.31, so the number lands on the same scale as both."""
    n, L = x.shape
    xb = x.reshape(n, L // b, b).sum(2)
    yb = y.reshape(n, L // b, b).sum(2)
    if xb.shape[1] < 3:
        return np.nan, 0
    xc = xb - xb.mean(1, keepdims=True)
    yc = yb - yb.mean(1, keepdims=True)
    num = (xc * yc).sum(1)
    den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1))
    good = den > 0
    if good.sum() == 0:
        return np.nan, 0
    return float((num[good] / den[good]).mean()), int(good.sum())


folds = json.load(open(FOLDS))
rows = []
for cell, C in CELLS.items():
    # Fold 0's validation chromosomes, matching what the converter was scored on, so the
    # baseline and the model number describe the same regions.
    els = load_peaks(C["elements"], folds["0"]["val"])
    ba = pyBigWig.open(C["atac"])
    bp = pyBigWig.open(C["dnase_p"]); bm = pyBigWig.open(C["dnase_m"])
    sizes = ba.chroms()
    A, Dn = [], []
    for _, r in els.iterrows():
        c = r["chr"]
        ctr = int(r["start"]) + int(r["summit"])
        s, e = ctr - HW, ctr + HW
        if c not in sizes or s < 0 or e > sizes[c]:
            continue
        a = ba.values(c, s, e, numpy=True)
        d1 = bp.values(c, s, e, numpy=True); d2 = bm.values(c, s, e, numpy=True)
        if a is None or d1 is None or d2 is None:
            continue
        if len(a) != OUT_W or len(d1) != OUT_W or len(d2) != OUT_W:
            continue
        A.append(np.nan_to_num(a, nan=0.0))
        Dn.append(np.nan_to_num(d1, nan=0.0) + np.nan_to_num(d2, nan=0.0))
    for b in (ba, bp, bm):
        b.close()
    A = np.stack(A); Dn = np.stack(Dn)
    # Stratify on the DNase total, the same quantity 2.31 strata on, so "topq" means the
    # same set of elements in both tables.
    tot = Dn.sum(1)
    q = np.quantile(tot, [0.2, 0.4, 0.6, 0.8])
    quint = np.digitize(tot, q)
    print(f"{cell}: {len(A):,} fold-0 held-out windows", flush=True)
    for b in BINS:
        for strat, mask in (("all", np.ones(len(tot), bool)), ("topq", quint == 4)):
            r_, n_ = shape_corr(Dn[mask], A[mask], b)
            rows.append({"cell": cell, "bin_bp": b, "stratum": strat,
                         "n_elements": n_, "shape_r_observed_atac_vs_dnase": r_})

df = pd.DataFrame(rows)
out = f"{P}/results/atac_vs_dnase_profile_baseline.tsv"
df.round(4).to_csv(out, sep="\t", index=False)
print("\nWrote", out)
print("\nObserved ATAC profile vs observed DNase profile, fold-0 held-out windows.")
print("This is what the converter's painted track must beat to be worth deploying.")
print(f"\n{'cell':>9}{'bin':>6}{'stratum':>9}{'r(ATAC, DNase)':>17}")
for _, r in df.iterrows():
    print(f"{r['cell']:>9}{int(r['bin_bp']):>6}{r['stratum']:>9}"
          f"{r['shape_r_observed_atac_vs_dnase']:>17.4f}")
PY
echo done
