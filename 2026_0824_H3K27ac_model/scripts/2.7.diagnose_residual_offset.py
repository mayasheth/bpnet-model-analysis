#!/usr/bin/env python
"""Empirically confirm what an offset-trained model's forward() actually returns.

Reading the code says MultiModalBPNet.forward() never adds the count offset -- it is
applied only inside fit() when scoring the loss -- so a model trained with
--count-offset-model should emit the RESIDUAL (observed - atac_pred), not logcounts.
That reading decides whether 2.4.evaluate_residual.py must subtract atac_pred or not,
and getting it backwards produces a plausible wrong number rather than an error. So
check it against the data instead of trusting the read.

Discriminating prediction, on fold 0 held-out elements:
  if forward() returns logcounts -> mean(pred) ~ mean(observed logcounts), several units
  if forward() returns residual  -> mean(pred) ~ 0, and pred+atac_pred ~ observed
"""
import json, os, sys
import numpy as np
import torch
from scipy.stats import pearsonr

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks, normalize_accessibility

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/data/atac.bw"
EL = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
TRIM, HW, FOLD = 557, 500, 0
OUT_W = 2 * HW
IN_W = OUT_W + 2 * TRIM
dev = "cuda" if torch.cuda.is_available() else "cpu"

val = json.load(open(FOLDS))[str(FOLD)]["val"]
els = load_peaks(EL, val)
print(f"fold{FOLD}: {len(els):,} held-out elements on {val}")

seqs, sigs, accs, _ = extract_windows(
    els, GEN, f"{P}/data/h3k27ac_5p_plus.bw", f"{P}/data/h3k27ac_5p_minus.bw",
    ATAC, IN_W, OUT_W, 0, is_peak=True)
obs = np.log1p(sigs.sum(axis=(1, 2)))
print(f"observed log1p counts: mean {obs.mean():.3f}  sd {obs.std():.3f}")


def predict(model_dir, mode, acc_raw):
    mp = f"{model_dir}/fold{FOLD}/multimodal_bpnet.torch"
    a = accs
    if mode in ("multimodal", "atac"):
        st = json.load(open(f"{model_dir}/fold{FOLD}/acc_normalization_stats.json"))
        a = normalize_accessibility(acc_raw, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, a], axis=1) if mode == "multimodal"
         else seqs if mode == "sequence" else a).astype(np.float32)
    m = torch.load(mp, map_location="cpu", weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = mode
    m = m.to(dev).eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), 256):
            _, lc = m(torch.from_numpy(X[i:i + 256]).to(dev))
            out.append(lc.squeeze(-1).cpu().numpy())
    m.to("cpu")
    return np.concatenate(out)


atac_pred = predict(f"{P}/models/atac5p_hw500_clw10", "atac", accs)
seq_pred = predict(f"{P}/models/sequence5p_hw500_clw10", "sequence", accs)
res_pred = predict(f"{P}/models/residual5pFIXED_hw500_clw10", "sequence", accs)

print()
print(f"{'model':<24} {'mean':>9} {'sd':>9} {'r(obs,raw)':>12} {'r(obs,raw+atac)':>17}")
for name, p in [("atac5p", atac_pred), ("sequence5p", seq_pred),
                ("residual5pFIXED", res_pred)]:
    r_raw = pearsonr(obs, p)[0]
    r_off = pearsonr(obs, p + atac_pred)[0]
    print(f"{name:<24} {p.mean():9.3f} {p.std():9.3f} {r_raw:12.4f} {r_off:17.4f}")

true_resid = obs - atac_pred
print()
print(f"true residual (obs - atac_pred): mean {true_resid.mean():.3f} sd {true_resid.std():.3f}")
print(f"r(true_resid, residual_raw_output)      = {pearsonr(true_resid, res_pred)[0]:.4f}")
print(f"r(true_resid, residual_minus_atac_pred) = {pearsonr(true_resid, res_pred - atac_pred)[0]:.4f}")

# verdict
d_logcounts = abs(res_pred.mean() - obs.mean())
d_residual = abs(res_pred.mean() - 0.0)
print()
if d_residual < d_logcounts:
    print(f"VERDICT: residual5pFIXED emits a RESIDUAL "
          f"(|mean-0|={d_residual:.3f} < |mean-obs_mean|={d_logcounts:.3f}). "
          f"2.4 must NOT subtract atac_pred for it.")
else:
    print(f"VERDICT: residual5pFIXED emits LOGCOUNTS "
          f"(|mean-obs_mean|={d_logcounts:.3f} < |mean-0|={d_residual:.3f}). "
          f"Revert the 'residual' flag in 2.4 -- the code reading was wrong.")
# sanity: the plain sequence model must look like logcounts, proving the test discriminates
print(f"CONTROL sequence5p: |mean-obs_mean|={abs(seq_pred.mean()-obs.mean()):.3f} "
      f"vs |mean-0|={abs(seq_pred.mean()):.3f} "
      f"-> {'logcounts (test discriminates)' if abs(seq_pred.mean()-obs.mean()) < abs(seq_pred.mean()) else 'UNEXPECTED - test is not discriminating'}")
