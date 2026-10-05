#!/usr/bin/env python3
"""Paired bootstrap of CRISPR AUPRC differences between scE2G models trained in the same run.

Each model's `training_predictions.tsv` holds held-out-chromosome predictions (`E2G.Score`) for
the same CRISPR E-G pairs (`Regulated` = label). Pairs are resampled jointly, so a difference's
interval reflects the models disagreeing, not the benchmark's sampling noise. AUPRC = area under
the precision-recall curve (sklearn precision_recall_curve + auc), checked against scE2G's own
`performance_across_models.tsv` point estimate. CRISPR pairs whose element is not a candidate
(`missing.*.tsv.gz`, identical across models of a run) are added as score 0, as scE2G does.

Usage: 4.2.paired_bootstrap.py --run DIR --models a,b,c --pairs b:a,c:a,c:b [--extra name=path] [--out f.tsv]
"""
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True, help="results_dir/<cluster>")
ap.add_argument("--models", required=True)
ap.add_argument("--pairs", required=True, help="comma list of new:reference")
ap.add_argument("--extra", default="", help="name=path/to/training_predictions.tsv from another run")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out", required=True)
a = ap.parse_args()

KEY = ["chr", "start", "end", "TargetGene", "name"]
paths = {m: f"{a.run}/{m}/model/training_predictions.tsv" for m in a.models.split(",")}
if a.extra:
    for kv in a.extra.split(","):
        k, v = kv.split("=")
        paths[k] = v
D = None
for m, p in paths.items():
    t = pd.read_csv(p, sep="\t", usecols=KEY + ["Regulated", "E2G.Score"]).rename(columns={"E2G.Score": m})
    assert not t.duplicated(KEY).any(), m
    D = t if D is None else D.merge(t, on=KEY + ["Regulated"], how="inner", validate="1:1")
MISS = "missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"
miss = None
for m in a.models.split(","):
    t = pd.read_csv(f"{a.run}/{m}/{MISS}", sep="\t", usecols=["name", "Regulated"]).sort_values("name")
    if miss is None:
        miss = t
    else:
        assert (t.values == miss.values).all(), f"{m}: missing pairs differ"
y = np.concatenate([D["Regulated"].astype(bool).values, miss["Regulated"].astype(bool).values])
print(f"{len(D):,} pairs, {y.sum():,} positives; models: {list(paths)}", flush=True)
for m in paths:
    assert len(D) == len(pd.read_csv(paths[m], sep="\t", usecols=["Regulated"])), f"{m}: pairs lost in join"


def auprc(s, yy):
    p, r, _ = precision_recall_curve(yy, s)
    return auc(r, p)


S = {m: np.concatenate([np.nan_to_num(D[m].values, nan=0.0), np.zeros(len(miss))]) for m in paths}
print(f"+ {len(miss):,} missing pairs ({int(miss['Regulated'].astype(bool).sum())} positives) at score 0", flush=True)
point = {m: auprc(S[m], y) for m in paths}
rng = np.random.default_rng(0)
idx = [rng.integers(0, len(y), len(y)) for _ in range(a.n_boot)]
boot = {m: np.array([auprc(S[m][i], y[i]) for i in idx]) for m in paths}
rows = []
for m in paths:
    lo, hi = np.percentile(boot[m], [2.5, 97.5])
    rows.append(dict(comparison=m, estimate=point[m], ci_lo=lo, ci_hi=hi, p_le_0=np.nan))
for pr in a.pairs.split(","):
    new, ref = pr.split(":")
    d = boot[new] - boot[ref]
    lo, hi = np.percentile(d, [2.5, 97.5])
    rows.append(dict(comparison=f"{new} - {ref}", estimate=point[new] - point[ref], ci_lo=lo, ci_hi=hi,
                     p_le_0=float((d <= 0).mean())))
T = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(T.round(4).to_string(index=False))
perf = os.path.join(os.path.dirname(a.run.rstrip("/")), "performance_across_models.tsv")
if os.path.exists(perf):
    P = pd.read_csv(perf, sep="\t").set_index("model")["AUPRC"]
    print("\nscE2G's own AUPRC:", P.round(4).to_dict())
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
T.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
