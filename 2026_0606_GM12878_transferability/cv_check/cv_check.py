"""Held-out (CV) versus mean-of-folds Pearson r for every bar in Fig. C.2 and Fig. 4.4c,d.

CV prediction for an element = prediction of the fold model whose test chromosomes contain it
(reference/hg38_five_folds.json). Mean = average of the five fold models.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
GM = f"{P}/2026_0606_GM12878_transferability"
MM = f"{P}/2026_0529_multimodal_p300_model"
K = f"{P}/2025_0517_official_EP300_K562_model"
OUT = f"{GM}/cv_check"
os.makedirs(OUT, exist_ok=True)

folds = json.load(open(f"{P}/reference/hg38_five_folds.json"))
test_fold = {c: int(f) for f, d in folds.items() for c in d["test"]}

# (figure, eval cell type, label, kind, path)
# kind "wide": mean_predictions.tsv.gz with pred_fold0..4 (bpnet-refactor models)
# kind "long": cv_predictions.tsv.gz (fold column) + mean_predictions.tsv.gz (multimodal family)
MODELS = [
    ("C.2", "K562", "K562 BPNet (sequence only, in-cell-type)", "wide", f"{K}/predictions_mean/all_folds/mean_predictions.tsv.gz"),
    ("C.2", "K562", "K562 BPNet (ATAC only, in-cell-type)", "long", f"{MM}/predictions/atac_only"),
    ("C.2", "K562", "K562 BPNet (sequence + ATAC, in-cell-type)", "long", f"{MM}/predictions/atac"),
    ("C.2", "K562", "GM12878 BPNet (sequence only, cross-cell-type)", "wide", f"{GM}/predictions/gm12878_bpnet_on_k562/mean/all_folds/mean_predictions.tsv.gz"),
    ("C.2", "K562", "GM12878 BPNet (ATAC only, cross-cell-type)", "long", f"{GM}/predictions/gm12878_atac_only_on_k562"),
    ("C.2", "K562", "GM12878 BPNet (sequence + ATAC, cross-cell-type)", "long", f"{GM}/predictions/gm12878_multimodal_on_k562"),
    ("4.4c,d", "GM12878", "GM12878 BPNet (sequence only)", "wide", f"{GM}/predictions/gm12878_bpnet/mean/all_folds/mean_predictions.tsv.gz"),
    ("4.4c,d", "GM12878", "GM12878 BPNet (ATAC only)", "long", f"{GM}/predictions/gm12878_atac_only"),
    ("4.4c,d", "GM12878", "GM12878 BPNet (sequence + ATAC)", "long", f"{GM}/predictions/gm12878_multimodal_atac"),
    ("4.4c,d", "GM12878", "K562 BPNet (sequence only)", "wide", f"{GM}/predictions/k562_bpnet_v1/mean/all_folds/mean_predictions.tsv.gz"),
    ("4.4c,d", "GM12878", "K562 BPNet (ATAC only)", "long", f"{GM}/predictions/k562_atac_only"),
    ("4.4c,d", "GM12878", "K562 BPNet (sequence + ATAC)", "long", f"{GM}/predictions/k562_multimodal_atac"),
]

OC = "EP300_peak_overlap"


def r(x, y):
    return float(pearsonr(x, y)[0])


rows = []
for fig, cell, label, kind, path in MODELS:
    if kind == "wide":
        df = pd.read_csv(path, sep="\t")
        df = df[df["chrom"].isin(test_fold)].copy()
        f = df["chrom"].map(test_fold).values
        preds = df[[f"pred_fold{i}" for i in range(5)]].values
        df["cv_pred"] = preds[np.arange(len(df)), f]
        cv, mean = df.rename(columns={"cv_pred": "pred"}), df.rename(columns={"mean_pred_logcounts": "pred"})
        fold_check = "derived from chromosome"
    else:
        cv = pd.read_csv(f"{path}/cv_predictions.tsv.gz", sep="\t").rename(columns={"pred_logcounts": "pred"})
        mean = pd.read_csv(f"{path}/mean_predictions.tsv.gz", sep="\t")
        pcol = "mean_pred_logcounts" if "mean_pred_logcounts" in mean.columns else "pred_logcounts"
        mean = mean.rename(columns={pcol: "pred"})
        exp = cv["chrom"].map(test_fold)
        mism = int((exp.notna() & (exp != cv["fold"])).sum())
        fold_check = f"fold column vs JSON mismatches: {mism}"
    for subset, m in [("all", None), ("p300+", 1)]:
        c = cv if m is None else cv[cv[OC] == m]
        mm = mean if m is None else mean[mean[OC] == m]
        rows.append(dict(figure=fig, eval_cell=cell, model=label, subset=subset,
                         n_cv=len(c), r_cv=round(r(c["true_logcounts"], c["pred"]), 4),
                         n_mean=len(mm), r_mean=round(r(mm["true_logcounts"], mm["pred"]), 4),
                         check=fold_check))
    print(rows[-2], rows[-1], sep="\n", flush=True)

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/cv_vs_mean_accuracy.tsv", sep="\t", index=False)

# Guard (pass direction): our CV and mean values must reproduce the stored accuracy tables.
expected = {
    ("K562 BPNet (sequence only, in-cell-type)", "all"): (0.651, 0.6878),
    ("K562 BPNet (sequence only, in-cell-type)", "p300+"): (0.5207, 0.6492),
    ("K562 BPNet (sequence + ATAC, in-cell-type)", "all"): (0.785, None),
    ("GM12878 BPNet (sequence + ATAC)", "all"): (None, 0.821),
    ("K562 BPNet (sequence + ATAC)", "all"): (None, 0.7925),
}
bad = []
for (model, subset), (ecv, emean) in expected.items():
    row = res[(res.model == model) & (res.subset == subset)].iloc[0]
    if ecv is not None and abs(row.r_cv - ecv) > 0.002:
        bad.append((model, subset, "cv", row.r_cv, ecv))
    if emean is not None and abs(row.r_mean - emean) > 0.002:
        bad.append((model, subset, "mean", row.r_mean, emean))
print("GUARD:", "OK" if not bad else bad, flush=True)
sys.exit(1 if bad else 0)
