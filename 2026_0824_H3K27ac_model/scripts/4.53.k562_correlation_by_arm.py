#!/usr/bin/env python3
"""Correlation of every ABC arm's K562 activity track with observed K562 p300.

WHY. The headline question of report 6 is whether the changes that raise correlation also raise
the CRISPR benchmark. This scores every arm on the quantity ABC consumes: the activity value of
each of the 153,545 K562 candidate elements (EnhancerList `H3K27ac.RPM`, which holds the
predicted p300 in these runs), against observed K562 p300 counted the same way (the
`p300obs_k562` arm).

    spearman   Spearman rho over all candidates
    partial    non-promoter candidates: Spearman of predicted with observed after removing ATAC
               (rank and rank^2) from both, i.e. what the track adds beyond accessibility

Usage: 4.53.k562_correlation_by_arm.py [--out results/2026_0930_k562_correlation_by_arm.tsv]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ABC = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
DIRS = ["2026_0905_p300_activity", "2026_0922_p300_multicell", "2026_0924_p300_rank",
        "2026_0925_p300_ranksweep", "2026_0924_p300_nbr", "2026_0924_p300_nbrx",
        "2026_0926_p300_gataoracle", "2026_0927_p300_motifact", "2026_0928_p300_motifact2",
        "2026_0928_p300_teach", "2026_0929_p300_seeds"]
KEY = ["chr", "start", "end"]
ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/2026_0930_k562_correlation_by_arm.tsv")
a = ap.parse_args()


def load(d):
    return pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t",
                       usecols=KEY + ["class", "H3K27ac.RPM", "ATAC.RPM"])


obs = load(f"{ABC}/2026_0905_p300_activity/p300obs_k562").rename(
    columns={"H3K27ac.RPM": "obs", "ATAC.RPM": "atac"})
obs["np"] = ~obs["class"].isin(["tss", "promoter"])
ro = rankdata(obs["obs"]) / len(obs)
npm = obs["np"].values
atac_np = rankdata(obs.loc[npm, "atac"]) / npm.sum()
Z = np.column_stack([np.ones(npm.sum()), atac_np, atac_np ** 2])
Q, _ = np.linalg.qr(Z)
res = lambda v: v - Q @ (Q.T @ v)
yo = res(rankdata(obs.loc[npm, "obs"]) / npm.sum())
rows = []
for d in DIRS:
    for arm in sorted(glob.glob(f"{ABC}/{d}/p300pred_*")):
        p = load(arm)
        m = obs[KEY].merge(p[KEY + ["H3K27ac.RPM"]], on=KEY, how="left", validate="one_to_one")
        v = m["H3K27ac.RPM"].values
        assert np.isfinite(v).all(), arm
        rho = np.corrcoef(rankdata(v) / len(v), ro)[0, 1]
        pr = np.corrcoef(res(rankdata(v[npm]) / npm.sum()), yo)[0, 1]
        rows.append(dict(arm=os.path.basename(arm), results_dir=d, spearman=rho, partial=pr))
        print(f"{os.path.basename(arm):<60} rho {rho:.4f}  partial {pr:.4f}", flush=True)
pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
print("wrote", a.out)
