#!/usr/bin/env python3
"""Paired bootstrap on CRISPR-benchmark AUPRC differences between arms.

WHY THIS IS NEEDED. `performance_summary.txt` reports a bootstrap CI per predictor
independently. Those intervals are ~+/-0.05 wide while the differences under discussion are
0.03-0.10, so every pair of arms has overlapping intervals and no comparison can be made from
them -- exactly the situation where an unpaired interval understates resolution. The arms are
scored on the IDENTICAL element-gene pair set, so resampling pairs once and recomputing every
arm on that same resample removes the shared pair-sampling variance, which is most of it.

This is the same argument as the within-fold pairing used for the model metrics (Report 2),
applied to the downstream benchmark.

Reports, for each requested comparison: the observed delta, a percentile CI, and the fraction
of resamples in which the sign is preserved. That last number is the one to quote -- it is a
direct statement about how often the ordering would reproduce.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("merged", help="expt_pred_merged_annot.txt.gz")
ap.add_argument("--pair", nargs=2, action="append", required=True,
                metavar=("A", "B"), help="report AUPRC(A) - AUPRC(B). Repeatable.")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out-tsv", default=None,
                help="Write the deltas, CIs and sign-retention to a TSV so the report's "
                     "headline numbers are registerable in the numbers manifest instead of "
                     "being read off a log and hand-typed.")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

d = pd.read_csv(a.merged, sep="\t", low_memory=False,
                usecols=["name", "measuredGeneSymbol", "Regulated", "pred_uid", "pred_value"])
d["pair"] = d["name"].astype(str) + "|" + d["measuredGeneSymbol"].astype(str)
d["y"] = d["Regulated"].astype(str).str.upper().eq("TRUE").astype(int)

# Wide: one row per element-gene pair, one column per arm. inner join semantics via dropna
# below, so every arm is scored on exactly the same pairs.
score = d.pivot_table(index="pair", columns="pred_uid", values="pred_value", aggfunc="max")
label = d.groupby("pair")["y"].max()
score = score.loc[label.index]

arms = sorted({x for p in a.pair for x in p})
missing = [x for x in arms if x not in score.columns]
if missing:
    raise SystemExit(f"arm(s) not in the merged table: {missing}\n"
                     f"available: {list(score.columns)}")

sub = score[arms].dropna()
y = label.loc[sub.index].to_numpy()
print(f"{len(sub):,} element-gene pairs scored by every requested arm; "
      f"{int(y.sum()):,} regulated ({100 * y.mean():.2f}%)\n")
if y.sum() < 20:
    raise SystemExit("too few positives for a meaningful AUPRC bootstrap")


def auprc(y_true, s):
    """Average precision, computed directly so no sklearn dependency is introduced."""
    o = np.argsort(-s, kind="mergesort")
    yt = y_true[o]
    tp = np.cumsum(yt)
    prec = tp / np.arange(1, len(yt) + 1)
    n_pos = yt.sum()
    return float((prec * yt).sum() / n_pos) if n_pos else np.nan


X = {c: sub[c].to_numpy(dtype=float) for c in arms}
obs = {c: auprc(y, X[c]) for c in arms}
print("observed AUPRC on the shared pair set")
for c in sorted(arms, key=lambda k: -obs[k]):
    print(f"   {c:<44}{obs[c]:.4f}")

rng = np.random.RandomState(a.seed)
n = len(y)
boot = {c: np.empty(a.n_boot) for c in arms}
for b in range(a.n_boot):
    idx = rng.randint(0, n, n)
    yb = y[idx]
    if yb.sum() < 5:                      # degenerate resample, redraw
        while yb.sum() < 5:
            idx = rng.randint(0, n, n)
            yb = y[idx]
    for c in arms:
        boot[c][b] = auprc(yb, X[c][idx])

print(f"\npaired bootstrap, {a.n_boot:,} resamples of the pair set (same resample for every arm)")
print(f"{'comparison':<62}{'delta':>9}{'95% CI':>22}{'sign kept':>11}")
for A, B in a.pair:
    dl = boot[A] - boot[B]
    lo, hi = np.percentile(dl, [2.5, 97.5])
    frac = float((dl > 0).mean()) if obs[A] >= obs[B] else float((dl < 0).mean())
    print(f"{A + ' - ' + B:<62}{obs[A] - obs[B]:>+9.4f}"
          f"{f'[{lo:+.4f}, {hi:+.4f}]':>22}{100 * frac:>10.1f}%")

print("\n'sign kept' is the fraction of resamples preserving the observed ordering. Read it as")
print("the reproducibility of the ranking, not as a p-value; the CI is the effect size.")

if a.out_tsv:
    rows = []
    for A, B in a.pair:
        dl = boot[A] - boot[B]
        lo, hi = np.percentile(dl, [2.5, 97.5])
        frac = float((dl > 0).mean()) if obs[A] >= obs[B] else float((dl < 0).mean())
        rows.append({"arm_a": A, "arm_b": B,
                     "auprc_a": obs[A], "auprc_b": obs[B],
                     "delta": obs[A] - obs[B], "ci_lo": lo, "ci_hi": hi,
                     "sign_kept_frac": frac, "sign_kept_pct": 100.0 * frac,
                     "n_pairs": len(sub), "n_regulated": int(y.sum()),
                     "n_boot": a.n_boot})
    pd.DataFrame(rows).round(6).to_csv(a.out_tsv, sep="\t", index=False)
    print(f"\nwrote {a.out_tsv}")
