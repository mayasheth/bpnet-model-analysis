#!/usr/bin/env python3
"""Paired difference between two 2.15 runs over the same configs and folds.

For an A/B where the CHANGE IS IN THE EVALUATOR rather than in the models -- test-time
reverse-complement averaging being the case in hand -- the two arms live in separate output
tables, so `--pair` inside a single run cannot express the comparison. This pairs the tables
on (fold, config) and reports mean difference, t-based 95% CI and a paired t-test per metric
per config.

Pairing on (fold, config) is what makes it sensitive: between-fold sd on this project is
0.041-0.046, so an unpaired comparison of a ~0.005 effect sees nothing.

Usage: 2.20.paired_delta_between_runs.py A.tsv B.tsv --label-a rc --label-b plain
Reports A - B, so a positive number means A is better (except profile_jsd, lower is better).
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, t as tdist

ap = argparse.ArgumentParser()
ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--label-a", default="A"); ap.add_argument("--label-b", default="B")
args = ap.parse_args()

A = pd.read_csv(args.a, sep="\t").sort_values(["config", "fold"]).reset_index(drop=True)
B = pd.read_csv(args.b, sep="\t").sort_values(["config", "fold"]).reset_index(drop=True)

for c in ("fold", "config"):
    assert A[c].equals(B[c]), f"tables are not over the same {c} values"
assert A["n"].equals(B["n"]), "region counts differ; the two runs are not comparable"

metrics = [c for c in A.columns if c not in ("fold", "config", "n")]
nf = A["fold"].nunique()
tc = tdist.ppf(0.975, df=nf - 1)
print(f"{args.label_a} - {args.label_b}   ({nf} folds, paired within fold)")
print(f"{'config':<22}{'metric':<24}{'delta':>10}{'95% CI':>24}{'p':>9}")
for cfg, ga in A.groupby("config", sort=False):
    gb = B[B["config"] == cfg]
    for m in metrics:
        x = ga[m].to_numpy(float); y = gb[m].to_numpy(float)
        if np.isnan(x).all() or np.isnan(y).all():
            continue
        d = x - y
        mu = d.mean()
        half = tc * d.std(ddof=1) / np.sqrt(len(d)) if d.std(ddof=1) > 0 else 0.0
        p = ttest_rel(x, y).pvalue if d.std(ddof=1) > 0 else 1.0
        star = " *" if p < 0.05 else ""
        print(f"{cfg:<22}{m:<24}{mu:>+10.4f}"
              f"{f'[{mu-half:+.4f}, {mu+half:+.4f}]':>24}{p:>9.4f}{star}")
