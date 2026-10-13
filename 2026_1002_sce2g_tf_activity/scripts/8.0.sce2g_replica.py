#!/usr/bin/env python3
"""8.0 Offline replica of scE2G's model training, for fast feature discovery (plan:
plans/8_crispr_feature_discovery.md, Phase 0).

scE2G (ENCODE_rE2G/workflow/scripts/model_training/train_model.py) fits
LogisticRegression(**training_params.pkl) on log(|x| + 0.01) features over the CRISPR pairs that have a
candidate element, with one model per held-out chromosome (E2G.Score). AUPRC as 4.2: sklearn
precision_recall_curve + auc, with the CRISPR pairs that miss every candidate added at score 0.

This refits every model of a run from that run's own training table and checks it against the
pipeline: per-pair LOCO scores and AUPRC. Checked 2026-10-06 on the 2026_1005 screen run: AUPRC within
0.0005 for all 4 models; per-pair scores differ by up to 0.04 (unpenalised lbfgs at tol 1e-4 stops at
slightly different points under this env's sklearn 1.7.2), so compare models within the replica. The library functions are imported by 8.1.

Usage: 8.0.sce2g_replica.py [--run DIR] [--tol 0.002]
"""
import argparse
import os
import pickle

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, precision_recall_curve

SC = os.environ.get("SCRATCH", "")
TAB = "for_training.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"
MISS = "missing.EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"
KEY = ["chr", "start", "end", "TargetGene", "name"]
EPS = 0.01


def load_run(run, model):
    d = f"{run}/{model}"
    df = pd.read_csv(f"{d}/{TAB}", sep="\t")
    miss = pd.read_csv(f"{d}/{MISS}", sep="\t", usecols=["Regulated"])
    with open(f"{d}/model/training_params.pkl", "rb") as fh:
        params = {**pickle.load(fh), "verbose": 0}
    return df, miss["Regulated"].astype(bool).values, params


def loco(df, features, params):
    """Leave-one-chromosome-out probabilities, exactly as train_model.py."""
    X = np.log(np.abs(df[features].values.astype(float)) + EPS)
    y = df["Regulated"].values.astype(np.int64)
    out = np.full(len(df), np.nan)
    for c in np.unique(df["chr"]):
        te = (df["chr"] == c).values
        out[te] = LogisticRegression(**params).fit(X[~te], y[~te]).predict_proba(X[te])[:, 1]
    return out


def auprc(y, s, n_miss_pos=0, n_miss=0):
    """AUPRC over scored pairs plus missing pairs at score 0 (4.2 / scE2G convention)."""
    y = np.r_[np.asarray(y, bool), np.ones(n_miss_pos, bool), np.zeros(n_miss - n_miss_pos, bool)]
    s = np.r_[np.asarray(s, float), np.zeros(n_miss)]
    p, r, _ = precision_recall_curve(y, s)
    return auc(r, p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1005_k562_xu_screen/K562_Xu_pl")
    ap.add_argument("--tol", type=float, default=0.002)
    a = ap.parse_args()
    rows = []
    for m in sorted(x for x in os.listdir(a.run) if os.path.exists(f"{a.run}/{x}/model/training_predictions.tsv")):
        df, miss, params = load_run(a.run, m)
        feats = pd.read_csv(f"{a.run}/{m}/model/model_coefficients.tsv", sep="\t").query("test_chr == 'none'")["feature"].tolist()
        pipe = pd.read_csv(f"{a.run}/{m}/model/training_predictions.tsv", sep="\t", usecols=KEY + ["E2G.Score"])
        mine = df[KEY].assign(rep=loco(df, feats, params)).merge(pipe, on=KEY, how="inner", validate="1:1")
        assert len(mine) == len(df), f"{m}: pairs lost joining to the pipeline's predictions"
        y = df.set_index(KEY).loc[pd.MultiIndex.from_frame(mine[KEY]), "Regulated"].astype(bool).values
        a_rep = auprc(y, mine["rep"], miss.sum(), len(miss))
        a_pipe = auprc(y, mine["E2G.Score"], miss.sum(), len(miss))
        # negative control: the same check on a wrong model (ARC.E2G.Score dropped) must fail
        assert "ARC.E2G.Score" in feats
        neg = auprc(df["Regulated"].astype(bool).values, loco(df, [f for f in feats if f != "ARC.E2G.Score"], params),
                    miss.sum(), len(miss))
        rows.append(dict(model=m, n_features=len(feats), params=str(params),
                         max_abs_score_diff=float(np.abs(mine["rep"] - mine["E2G.Score"]).max()),
                         auprc_replica=a_rep, auprc_pipeline=a_pipe, auprc_wrong_model=neg))
    R = pd.DataFrame(rows)
    print(R.drop(columns="params").round(5).to_string(index=False))
    print("params:", R["params"].iloc[0])
    bad = R[(R["auprc_replica"] - R["auprc_pipeline"]).abs() > a.tol]
    assert bad.empty, f"replica differs from the pipeline by > {a.tol}:\n{bad}"
    assert ((R["auprc_wrong_model"] - R["auprc_pipeline"]).abs() > a.tol).all(), "check cannot tell a wrong model apart"
    print(f"replica matches the pipeline within {a.tol} AUPRC for all {len(R)} models")
