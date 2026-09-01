#!/usr/bin/env python
"""Per-fold metrics with CIs and paired tests, driven by a 2.4-style config JSON.

Generalises 2.8 / 2.11, which were copies with hardcoded cell-type paths. Any comparison
expressible as a 2.4 config -- in-cell grids, cross-cell-type transfer, input-definition
swaps -- can now be scored per fold without another copy.

Emits, for every compare entry against the config's baseline:
  overall_pearson       r(observed, full prediction). THE metric when the question is
                        "how well would this predict H3K27ac in a cell type where I have
                        none" -- it scores the quantity you actually want to produce.
  overall_pearson_topq  the same on the top signal quintile (project reporting standard).
  residual_pearson      r(observed - atac_pred, model_pred - atac_pred), the mechanistic
                        readout of what the model adds beyond the baseline.
  incremental_r2        R2(model) - R2(baseline) against the observed signal.

Residual-objective models are handled via `"residual": true` in the config entry: their
forward() emits the residual, so the baseline prediction is added back before scoring.

Usage: 2.15.perfold_from_config.py CONFIG_JSON OUT_PREFIX ELEMENTS [--pair A B]...
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, ttest_rel, t as tdist

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks, normalize_accessibility

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
TRIM, HW = 557, 500
OUT_W = 2 * HW
IN_W = OUT_W + 2 * TRIM
TC = tdist.ppf(0.975, df=4)
dev = "cuda" if torch.cuda.is_available() else "cpu"

ap = argparse.ArgumentParser()
ap.add_argument("config"); ap.add_argument("out_prefix"); ap.add_argument("elements")
ap.add_argument("--pair", nargs=2, action="append", default=[])
a = ap.parse_args()
spec = json.load(open(a.config))
entries = [spec["baseline"]] + spec["compare"]


def predict(cfg, fold, seqs, accs_raw):
    md, mode = cfg["model_dir"], cfg["mode"]
    marker = f"{md}/fold{fold}/training_complete.json"
    if not os.path.exists(marker):
        raise SystemExit(f"error: {marker} missing; refusing to score an unfinished fold.")
    x = accs_raw
    if mode in ("multimodal", "atac"):
        st = json.load(open(f"{md}/fold{fold}/acc_normalization_stats.json"))
        x = normalize_accessibility(accs_raw, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, x], axis=1) if mode == "multimodal"
         else seqs if mode == "sequence" else x).astype(np.float32)
    m = torch.load(f"{md}/fold{fold}/multimodal_bpnet.torch", map_location="cpu",
                   weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = mode
    m = m.to(dev).eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), 256):
            _, lc = m(torch.from_numpy(X[i:i + 256]).to(dev))
            out.append(lc.squeeze(-1).cpu().numpy())
    m.to("cpu"); del X
    return np.concatenate(out)


rows = []
folds_json = json.load(open(FOLDS))
for fold in range(5):
    els = load_peaks(a.elements, folds_json[str(fold)]["val"])
    b = spec["baseline"]
    seqs, sigs, accs, _ = extract_windows(
        els, GEN, b["signal_plus_bw"], b.get("signal_minus_bw"),
        b["accessibility_bw"], IN_W, OUT_W, 0, is_peak=True)
    obs = np.log1p(sigs.sum(axis=(1, 2))); del sigs
    base = predict(b, fold, seqs, accs)
    true_resid = obs - base
    top = obs >= np.quantile(obs, 0.8)
    r2b, r2bt = pearsonr(obs, base)[0] ** 2, pearsonr(obs[top], base[top])[0] ** 2
    for cfg in spec["compare"]:
        raw = predict(cfg, fold, seqs, accs)
        full = raw + base if cfg.get("residual") else raw
        mres = raw if cfg.get("residual") else raw - base
        rows.append({"fold": fold, "config": cfg["label"], "n": len(obs),
                     "overall_pearson": pearsonr(obs, full)[0],
                     "overall_pearson_topq": pearsonr(obs[top], full[top])[0],
                     "residual_pearson": pearsonr(true_resid, mres)[0],
                     "incremental_r2": pearsonr(obs, full)[0] ** 2 - r2b,
                     "incremental_r2_topq": pearsonr(obs[top], full[top])[0] ** 2 - r2bt})
    rows.append({"fold": fold, "config": b["label"], "n": len(obs),
                 "overall_pearson": pearsonr(obs, base)[0],
                 "overall_pearson_topq": pearsonr(obs[top], base[top])[0],
                 "residual_pearson": np.nan, "incremental_r2": 0.0,
                 "incremental_r2_topq": 0.0})
    del seqs, accs
    print(f"fold{fold}: n={len(obs):,}", flush=True)

df = pd.DataFrame(rows)
p1 = f"{P}/results/{a.out_prefix}per_fold.tsv"
df.round(4).to_csv(p1, sep="\t", index=False)
print("\nWrote", p1)

METRICS = ["overall_pearson", "overall_pearson_topq", "residual_pearson",
           "incremental_r2", "incremental_r2_topq"]
srows = []
for label, g in df.groupby("config", sort=False):
    r = {"config": label}
    for m in METRICS:
        v = g[m].to_numpy(float)
        if np.isnan(v).all():
            r[m] = "-"
            continue
        mu = v.mean(); half = TC * v.std(ddof=1) / np.sqrt(len(v))
        r[m] = f"{mu:.3f} [{mu-half:.3f}, {mu+half:.3f}]"
    srows.append(r)
summ = pd.DataFrame(srows)
p2 = f"{P}/results/{a.out_prefix}fold_summary.tsv"
summ.to_csv(p2, sep="\t", index=False)
print("Wrote", p2)
print(summ.to_string(index=False))

if a.pair:
    print("\n--- paired differences (A - B), within fold ---")
    piv = {m: df.pivot(index="fold", columns="config", values=m) for m in METRICS}
    for A, B in a.pair:
        for m in ["overall_pearson", "overall_pearson_topq"]:
            d = (piv[m][A] - piv[m][B]).to_numpy()
            mu = d.mean(); half = TC * d.std(ddof=1) / np.sqrt(len(d))
            p = ttest_rel(piv[m][A], piv[m][B]).pvalue
            print(f"  {m:<22} {A} - {B}: {mu:+.4f} [{mu-half:+.4f}, {mu+half:+.4f}]  p={p:.4f}")
