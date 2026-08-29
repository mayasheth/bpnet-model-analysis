#!/usr/bin/env python
"""Per-fold residual statistics with CIs, plus the artifact check the pooled table cannot do.

Two problems with the pooled residual_evaluation table:

1. It is POOLED. Between-fold sd in this project is 0.041-0.046, ~2-3x the run-to-run
   variance, so a pooled correlation reports these differences with no error bar. Project
   standard is mean +/- t-based 95% CI over the 5 folds.

2. residual_pearson has a free artifact channel. true_resid = obs - atac_pred, so ANY
   prediction that anti-correlates with atac_pred scores positive residual_pearson without
   containing information about what accessibility actually missed. A sequence model can
   partially predict accessibility, so this is a live risk, not a hypothetical -- and it
   would inflate exactly the number we would want to headline.

Checks that separate signal from artifact:
  r_out_vs_atac      r(model raw output, atac_pred). Large negative = artifact channel open.
  partial_r          r(true_resid, model_resid) with atac_pred partialled out of BOTH.
                     This is residual_pearson stripped of the shared -atac_pred term.
  incremental_r2     R2(obs, full_pred) - R2(obs, atac_pred). Out-of-sample gain in
                     predicting the OBSERVED signal -- immune to the artifact, since it
                     scores against obs rather than against a quantity containing atac_pred.
"""
import json, os, sys
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, t as tdist

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks, normalize_accessibility

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/data/atac.bw"
EL = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
TRIM, HW = 557, 500
OUT_W = 2 * HW
IN_W = OUT_W + 2 * TRIM
dev = "cuda" if torch.cuda.is_available() else "cpu"

MODELS = [
    ("atac5p",      "atac",       f"{P}/models/atac5p_hw500_clw10",            False),
    ("sequence5p",  "sequence",   f"{P}/models/sequence5p_hw500_clw10",        False),
    ("multimodal5p","multimodal", f"{P}/models/multimodal5p_hw500_clw10",      False),
    ("residual5p",  "sequence",   f"{P}/models/residual5pFIXED_hw500_clw10",   True),
]


def predict(model_dir, mode, fold, seqs, accs_raw):
    a = accs_raw
    if mode in ("multimodal", "atac"):
        st = json.load(open(f"{model_dir}/fold{fold}/acc_normalization_stats.json"))
        a = normalize_accessibility(accs_raw, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, a], axis=1) if mode == "multimodal"
         else seqs if mode == "sequence" else a).astype(np.float32)
    m = torch.load(f"{model_dir}/fold{fold}/multimodal_bpnet.torch",
                   map_location="cpu", weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = mode
    m = m.to(dev).eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), 256):
            _, lc = m(torch.from_numpy(X[i:i + 256]).to(dev))
            out.append(lc.squeeze(-1).cpu().numpy())
    m.to("cpu")
    del X
    return np.concatenate(out)


def partial_r(a, b, c):
    """r(a,b) with c linearly removed from both."""
    ra = a - np.polyval(np.polyfit(c, a, 1), c)
    rb = b - np.polyval(np.polyfit(c, b, 1), c)
    if ra.std() < 1e-9 or rb.std() < 1e-9:
        return np.nan
    return float(pearsonr(ra, rb)[0])


rows = []
folds_json = json.load(open(FOLDS))
for fold in range(5):
    val = folds_json[str(fold)]["val"]
    els = load_peaks(EL, val)
    seqs, sigs, accs, _ = extract_windows(
        els, GEN, f"{P}/data/h3k27ac_5p_plus.bw", f"{P}/data/h3k27ac_5p_minus.bw",
        ATAC, IN_W, OUT_W, 0, is_peak=True)
    obs = np.log1p(sigs.sum(axis=(1, 2)))
    del sigs

    raw = {}
    for label, mode, mdir, _ in MODELS:
        raw[label] = predict(mdir, mode, fold, seqs, accs)
    del seqs, accs

    atac_pred = raw["atac5p"]
    true_resid = obs - atac_pred
    # top signal quintile: the project's standard stratum, distinct from |resid| quintiles
    top = obs >= np.quantile(obs, 0.8)
    r2_atac = pearsonr(obs, atac_pred)[0] ** 2
    r2_atac_top = pearsonr(obs[top], atac_pred[top])[0] ** 2

    for label, mode, mdir, is_res in MODELS:
        if label == "atac5p":
            continue
        out = raw[label]
        full = out + atac_pred if is_res else out
        mres = out if is_res else out - atac_pred
        rows.append({
            "fold": fold, "config": label, "n": len(obs),
            "residual_pearson": pearsonr(true_resid, mres)[0],
            "partial_r": partial_r(true_resid, mres, atac_pred),
            "r_out_vs_atac": pearsonr(out, atac_pred)[0],
            "overall_pearson": pearsonr(obs, full)[0],
            "incremental_r2": pearsonr(obs, full)[0] ** 2 - r2_atac,
            "overall_pearson_topq": pearsonr(obs[top], full[top])[0],
            "incremental_r2_topq": pearsonr(obs[top], full[top])[0] ** 2 - r2_atac_top,
        })
    print(f"fold{fold}: n={len(obs):,} atac r={pearsonr(obs, atac_pred)[0]:.4f} done", flush=True)

df = pd.DataFrame(rows)
p1 = f"{P}/results/residual5p_per_fold.tsv"
df.round(4).to_csv(p1, sep="\t", index=False)
print("\nWrote", p1)
print(df.round(4).to_string(index=False))

TCRIT = tdist.ppf(0.975, df=4)
metrics = ["residual_pearson", "partial_r", "r_out_vs_atac", "overall_pearson",
           "incremental_r2", "overall_pearson_topq", "incremental_r2_topq"]
srows = []
for label, g in df.groupby("config", sort=False):
    r = {"config": label}
    for m in metrics:
        v = g[m].to_numpy()
        half = TCRIT * v.std(ddof=1) / np.sqrt(len(v))
        r[m] = f"{v.mean():.3f} [{v.mean()-half:.3f}, {v.mean()+half:.3f}]"
    srows.append(r)
summ = pd.DataFrame(srows)
p2 = f"{P}/results/residual5p_fold_summary.tsv"
summ.to_csv(p2, sep="\t", index=False)
print("\nWrote", p2)
print(summ.to_string(index=False))
