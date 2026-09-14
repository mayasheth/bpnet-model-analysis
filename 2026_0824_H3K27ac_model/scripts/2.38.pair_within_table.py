#!/usr/bin/env python3
"""Paired within-fold differences between configs that are ALREADY in one 2.15 per-fold table.

WHY THIS EXISTS. `2.15 --pair` computes the same thing, but only while it is running
inference, so recovering a paired comparison from a finished run means re-scoring every fold
on a GPU. `2.20` pairs ACROSS two tables and asserts identical config sets, so it cannot
express "GATE - ungated" when both labels sit in the same table. This reads the table.

Pairing on fold is what makes it sensitive: between-fold sd on this project is 0.041-0.046,
so an unpaired comparison of a ~0.005 effect sees nothing, which is exactly the regime the
gate and asymmetric-loss arms live in.

Usage: 2.38.pair_within_table.py PERFOLD.tsv --pair A B [--pair C B] [--metrics m1 m2]
Reports A - B, so positive means A is better, except profile_jsd where lower is better.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, t as tdist

ap = argparse.ArgumentParser()
ap.add_argument("table")
ap.add_argument("--pair", nargs=2, action="append", required=True)
ap.add_argument("--metrics", nargs="+", default=None,
                help="default: every numeric column except fold/config/n")
a = ap.parse_args()

T = pd.read_csv(a.table, sep="\t")
metrics = a.metrics or [c for c in T.columns if c not in ("fold", "config", "n")]
missing = [m for m in metrics if m not in T.columns]
if missing:
    raise SystemExit(f"no such column(s) in {a.table}: {missing}")

W = max(len(f"{A} - {B}") for A, B in a.pair) + 2
print(f"{a.table}\n{'comparison':<{W}}{'metric':<24}{'delta':>10}{'95% CI':>24}{'p':>9}")
for A, B in a.pair:
    for label in (A, B):
        if label not in set(T["config"]):
            raise SystemExit(f"config {label!r} is not in the table; have "
                             f"{sorted(set(T['config']))}")
    ga = T[T["config"] == A].sort_values("fold")
    gb = T[T["config"] == B].sort_values("fold")
    if not ga["fold"].to_list() == gb["fold"].to_list():
        raise SystemExit(f"{A} and {B} were not scored on the same folds")
    # Same fold means the same held-out regions, so a differing n means the two arms did not
    # see the same windows and the pairing would be comparing different denominators.
    if "n" in T.columns and not np.array_equal(ga["n"].to_numpy(), gb["n"].to_numpy()):
        raise SystemExit(f"{A} and {B} disagree on n per fold; not comparable")
    nf = len(ga)
    tc = tdist.ppf(0.975, df=nf - 1)
    for m in metrics:
        x, y = ga[m].to_numpy(float), gb[m].to_numpy(float)
        if np.isnan(x).all() or np.isnan(y).all():
            continue
        d = x - y
        mu, sd = d.mean(), d.std(ddof=1)
        if nf < 2 or sd == 0:
            print(f"{A + ' - ' + B:<{W}}{m:<24}{mu:>+10.4f}{'(no CI)':>24}{'-':>9}")
            continue
        half = tc * sd / np.sqrt(nf)
        p = ttest_rel(x, y).pvalue
        star = " *" if p < 0.05 else ""
        print(f"{A + ' - ' + B:<{W}}{m:<24}{mu:>+10.4f}"
              f"{f'[{mu - half:+.4f}, {mu + half:+.4f}]':>24}{p:>9.4f}{star}")
