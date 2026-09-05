#!/usr/bin/env python3
"""Do p300 models predict the elements the H3K27ac model misses?

The H3K27ac model's under-predicted tail is 6.6x EP300-peak-enriched, 3.1x H3K4me1, zero
H3K27me3, CTCF-depleted -- canonical active enhancers treated as background. That makes p300
look like the missing signal. But observed p300 being present there says nothing about
whether a p300 MODEL can predict it from ATAC + sequence, and only the latter is deployable.

This is the prerequisite for both p300 ideas (multi-head, and stacking predicted p300):
  if predicted p300 IS elevated at those elements, p300 is learnable there and an auxiliary
  p300 task should transfer that signal into the H3K27ac model;
  if predicted p300 is NOT elevated, the p300 models fail on the same elements for the same
  reason, and neither idea helps -- drop both.

Everything is reported as fold-elevation against each track's OWN genome-wide median, because
painted-prediction bigwigs and observed read-count tracks are on different scales and raw
cross-track comparison is meaningless (a mistake made once already in 4.12).
"""
import os
import numpy as np
import pandas as pd
import pyBigWig

D = "/oak/stanford/groups/engreitz/Users/sheth"
E = f"{D}/EP300_BPNet"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS_K27 = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt"
PRED_K27 = f"{ABC}/2026_0903_predicted_activity/pred_k562_multimodal/Neighborhoods/EnhancerList.txt"
PRED_DIR = f"{E}/2026_0824_H3K27ac_model/data/abc_predicted"
P300_OBS = [f"{E}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig",
            f"{E}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"]
P300_CTRL = [f"{E}/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_plus.bigWig",
             f"{E}/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_minus.bigWig"]
COLS = ["chr", "start", "end", "H3K27ac.RPM", "ATAC.RPM", "class"]


def load(path, tag):
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df.index = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    return df.rename(columns={"H3K27ac.RPM": tag})


d = load(OBS_K27, "k27_obs").join(load(PRED_K27, "k27_pred")[["k27_pred"]], how="inner")
print(f"{len(d):,} ABC regions")


def bw_sum(paths, df, label):
    """Sum of bigwig values over each region, added across the listed tracks."""
    tot = np.zeros(len(df))
    for path in paths:
        if not os.path.exists(path):
            print(f"   MISSING {label}: {path}")
            return None
        bw = pyBigWig.open(path)
        chroms = bw.chroms()
        v = np.zeros(len(df))
        for i, (c, s, e) in enumerate(zip(df["chr"], df["start"], df["end"])):
            if c in chroms and e <= chroms[c]:
                v[i] = bw.stats(c, int(s), int(e), type="sum", exact=True)[0] or 0.0
        bw.close()
        tot += v
    return tot


d["p300_obs"] = bw_sum(P300_OBS, d, "observed p300")
d["p300_ctrl"] = bw_sum(P300_CTRL, d, "p300 input control")
for arm, tag in (("multimodal", "p300_pred_mm"), ("atac", "p300_pred_atac")):
    f = f"{PRED_DIR}/predp300_{arm}.bw"
    d[tag] = bw_sum([f], d, f"predicted p300 ({arm})")

tracks = [t for t in ("k27_obs", "k27_pred", "p300_obs", "p300_ctrl",
                      "p300_pred_mm", "p300_pred_atac") if d.get(t) is not None]

# Strata from the H3K27ac model's accessibility-residualised error, as in 4.10
d["latac"] = np.log1p(d["ATAC.RPM"])
A = np.vstack([np.ones(len(d)), d["latac"].to_numpy()]).T
for c in ("k27_obs", "k27_pred"):
    beta, *_ = np.linalg.lstsq(A, np.log1p(d[c]).to_numpy(), rcond=None)
    d["r_" + c] = np.log1p(d[c]).to_numpy() - A @ beta
d["err"] = d["r_k27_pred"] - d["r_k27_obs"]
q = d["err"].quantile([0.01, 0.05, 0.25, 0.75, 0.95, 0.99])
STRATA = [("H3K27ac under-pred 1%", d["err"] <= q[0.01]),
          ("H3K27ac under-pred 5%", d["err"] <= q[0.05]),
          ("typical (middle 50%)", (d["err"] > q[0.25]) & (d["err"] < q[0.75])),
          ("H3K27ac over-pred 5%", d["err"] >= q[0.95]),
          ("H3K27ac over-pred 1%", d["err"] >= q[0.99])]

med = {t: max(d[t].median(), 1e-9) for t in tracks}
print("\nFold elevation vs each track's OWN genome-wide median")
print(f"{'stratum':<24}{'n':>7}" + "".join(f"{t:>15}" for t in tracks))
for lab, mask in STRATA:
    s = d[mask]
    print(f"{lab:<24}{len(s):>7,}" +
          "".join(f"{s[t].median()/med[t]:>15.2f}" for t in tracks))

print("\nVERDICT on whether p300 is learnable where H3K27ac fails")
u = d[STRATA[0][1]]
obs_fold = u["p300_obs"].median() / med["p300_obs"]
if "p300_pred_mm" in tracks:
    pred_fold = u["p300_pred_mm"].median() / med["p300_pred_mm"]
    ctrl_fold = u["p300_ctrl"].median() / med["p300_ctrl"] if "p300_ctrl" in tracks else float("nan")
    print(f"   observed p300 elevation {obs_fold:.2f}x, input control {ctrl_fold:.2f}x, "
          f"predicted p300 {pred_fold:.2f}x")
    if obs_fold / max(ctrl_fold, 1e-9) < 1.5:
        print("   observed p300 barely exceeds its own input control here: the peak-based")
        print("   enrichment may be an accessibility artefact, so treat it cautiously.")
    elif pred_fold >= 0.5 * obs_fold:
        print("   p300 IS predicted at these elements. An auxiliary p300 task can carry")
        print("   signal the H3K27ac objective misses -- multi-head is worth building.")
    else:
        print("   p300 is NOT predicted here even though it is present. The p300 models fail")
        print("   on the same elements, so neither multi-head nor stacking helps. Drop both.")
