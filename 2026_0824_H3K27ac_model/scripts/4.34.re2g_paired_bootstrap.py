#!/usr/bin/env python3
"""Paired bootstrap of AUPRC differences between ENCODE-rE2G models trained by 4.32/4.33.

rE2G's own performance table gives each model an independent bootstrap interval, which is
too wide to compare models that are scored on the same CRISPR pairs. This joins every
model's leave-one-chromosome-out predictions (the E2G.Score column of
model_results/<model>/model/training_predictions.tsv; E2G.Score_full is the model trained on
all chromosomes and would leak) on the element-gene pair, and resamples pairs once per draw
for all models.

EVALUATED EXACTLY AS rE2G's compare_all_models.py DOES, which a first version of this script
did not and read ~0.045 too high for every model: (1) CRISPR pairs whose element falls
outside every candidate region are ADDED BACK with score 0 (about 11.6% of pairs, 31 of them
regulated), since a model that cannot score a regulated pair has missed it; (2) AUPRC is
rE2G's own statistic_aupr (precision_recall_curve_modified + auc_mod), imported from the
checkout. Pairs are outer-joined across models with a missing score set to 0, which is the
same rule applied per model. A self-check reproduces each model's AUPRC in rE2G's table.

Usage: 4.34.re2g_paired_bootstrap.py RESULTS_DIR --pair A B [--pair ...] [--out-tsv T]
"""
import argparse
import os

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/workflow/scripts/"
                   "model_training")
from training_functions import statistic_aupr  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("results_dir", help="rE2G results dir, e.g. .../2026_0923_p300_features")
ap.add_argument("--pair", nargs=2, action="append", required=True, metavar=("A", "B"))
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--out-tsv", default=None)
ap.add_argument("--group", action="append", default=[], metavar="NAME=A,B,...",
                help="a replicate group: its AUPRC is the mean of its members' AUPRCs, computed "
                     "within each bootstrap draw; NAME can then be used in --pair")
a = ap.parse_args()
GROUPS = {g_.split("=", 1)[0]: g_.split("=", 1)[1].split(",") for g_ in a.group}

KEY = ["chr", "start", "end", "TargetGene"]
MR = os.path.join(a.results_dir, "model_results")


def aupr(y, s):
    return statistic_aupr(np.asarray(y).astype(np.int64), np.asarray(s, dtype=float))


models = sorted({m_ for x in {m for pr in a.pair for m in pr} for m_ in GROUPS.get(x, [x])})
reported = pd.read_csv(os.path.join(MR, "performance_across_models.tsv"), sep="\t")
reported = dict(zip(reported["model"], reported["AUPRC"]))

frames = []
print("self-check against rE2G's table (each model on its own pairs):")
for m in models:
    d = pd.read_csv(os.path.join(MR, m, "model", "training_predictions.tsv"), sep="\t",
                    usecols=KEY + ["Regulated", "E2G.Score"])
    miss = pd.read_csv(os.path.join(MR, m,
                                    "merged_CRISPR_dataset.missing_from_features.NAfilled.tsv.gz"),
                       sep="\t", usecols=KEY + ["Regulated"]).assign(**{"E2G.Score": 0.0})
    d = pd.concat([d, miss], ignore_index=True)
    own = aupr(d["Regulated"].astype(bool), d["E2G.Score"])
    print(f"  {m:<36} recomputed {own:.4f}   rE2G {reported.get(m, float('nan')):.4f}   "
          f"n={len(d):,}")
    frames.append(d.rename(columns={"E2G.Score": m, "Regulated": f"reg_{m}"}))

# Duplicate keys: a CRISPR pair can appear twice when the element overlaps two candidate
# regions. rE2G keeps every row, so do the same by numbering duplicates within each key.
for f in frames:
    f["_dup"] = f.groupby(KEY).cumcount()
w = frames[0]
for f in frames[1:]:
    w = w.merge(f, on=KEY + ["_dup"], how="outer")
regcols = [c for c in w.columns if c.startswith("reg_")]
lab = w[regcols].bfill(axis=1).iloc[:, 0]
agree = w[regcols].apply(lambda r: r.dropna().nunique() <= 1, axis=1)
assert agree.all(), "Regulated labels disagree between models"
for m in models:
    w[m] = w[m].fillna(0.0)
y = lab.astype(bool).values
print(f"\n{len(w):,} element-gene pairs in the union; {int(y.sum())} regulated "
      f"({100 * y.mean():.2f}%). A pair a model did not score counts as score 0, as in rE2G.")

obs = {m: aupr(y, w[m].values) for m in models}
for m in sorted(models, key=lambda k: -obs[k]):
    print(f"  {m:<36} {obs[m]:.4f}")

rng = np.random.default_rng(a.seed)
n = len(w)
boot = {m: np.empty(a.n_boot) for m in models}
for b in range(a.n_boot):
    idx = rng.integers(0, n, n)
    yb = y[idx]
    for m in models:
        boot[m][b] = aupr(yb, w[m].values[idx])

# replicate groups: mean AUPRC over members, draw by draw
for g_, mem in GROUPS.items():
    obs[g_] = float(np.mean([obs[m_] for m_ in mem]))
    boot[g_] = np.mean([boot[m_] for m_ in mem], axis=0)
    print(f"group {g_}: " + ", ".join(f"{obs[m_]:.4f}" for m_ in mem)
          + f"  -> mean {obs[g_]:.4f}, sd {np.std([obs[m_] for m_ in mem], ddof=1):.4f}")

rows = []
print(f"\npaired bootstrap, {a.n_boot:,} resamples (same resample for every model)")
for A, B in a.pair:
    dd = boot[A] - boot[B]
    delta = obs[A] - obs[B]
    lo, hi = np.percentile(dd, [2.5, 97.5])
    kept = float(np.mean(np.sign(dd) == np.sign(delta)))
    print(f"  {A} - {B}\n      {delta:+.4f}  [{lo:+.4f}, {hi:+.4f}]  sign kept {100 * kept:.1f}%")
    rows.append(dict(model_a=A, model_b=B, auprc_a=obs[A], auprc_b=obs[B], delta=delta,
                     ci_lo=lo, ci_hi=hi, sign_kept_pct=100 * kept, n_pairs=n,
                     n_regulated=int(y.sum()), n_boot=a.n_boot))
if a.out_tsv:
    pd.DataFrame(rows).to_csv(a.out_tsv, sep="\t", index=False)
    print("wrote", a.out_tsv)
