#!/usr/bin/env python
"""Add the missing `n_topq` column to a 2.15 per-fold table.

2.15 records `n = len(obs)` and nothing else, but half its columns (`*_topq`) are computed
on `obs >= quantile(obs, 0.8)`. Figures that annotate a top-quintile panel with the table's
`n` therefore overstate the sample size roughly fivefold. This script recovers the real
count so the annotation can be correct without re-scoring anything.

NOTHING IS PREDICTED HERE. `top` depends only on the observed signal, so the whole
computation is the element loading and window extraction from 2.15 with the model forward
passes deleted. Checkpoints are still opened, because `model.trimming` sets IN_W_MAX and
IN_W_MAX sets which regions near chromosome ends survive -- but only on CPU and only for
that attribute.

THE `n` COLUMN IS THE CHECK. If this script reconstructs a different element set from the
one that was scored -- wrong elements file, wrong config, a fold definition that has moved
-- the recomputed `n` will not match the stored `n` and the script exits without writing.
That is the whole reason it recomputes a value it already knows.

DO NOT DERIVE n_topq AS n/5. `np.quantile(obs, 0.8)` is a VALUE threshold, so ties at the
threshold all fall inside the top quintile. H3K27ac windows are frequently zero, so the
tie mass is real and the count is not n/5. (2.2 is different: it ranks with `pd.qcut`, so
its `top_quintile` stratum genuinely is a fifth. The two must not be conflated.)

Usage: 2.30.backfill_n_topq.py CONFIG_JSON OUT_PREFIX ELEMENTS [--dry-run]
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import torch

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
HW = 500
OUT_W = 2 * HW

ap = argparse.ArgumentParser()
ap.add_argument("config"); ap.add_argument("out_prefix"); ap.add_argument("elements")
ap.add_argument("--dry-run", action="store_true",
                help="Report the counts and the n check, write nothing.")
a = ap.parse_args()

spec = json.load(open(a.config))
entries = [spec["baseline"]] + spec["compare"]

table = f"{P}/results/{a.out_prefix}per_fold.tsv"
if not os.path.exists(table):
    raise SystemExit(f"error: {table} not found")
df = pd.read_csv(table, sep="\t")
if "n_topq" in df.columns:
    raise SystemExit(f"{table} already has n_topq; nothing to do")

# Same window resolution as 2.15: read the receptive field off each entry's own fold0
# checkpoint and extract at the largest, so the valid-region mask matches what was scored.
for cfg in entries:
    m0 = torch.load(f'{cfg["model_dir"]}/fold0/multimodal_bpnet.torch',
                    map_location="cpu", weights_only=False)
    cfg["_in_window"] = OUT_W + 2 * m0.trimming
    if "accessibility_bw" not in cfg:
        cfg["accessibility_bw"] = spec["baseline"].get("accessibility_bw")
    del m0
IN_W_MAX = max(c["_in_window"] for c in entries)
print(f"extracting at in_window={IN_W_MAX}")

b = spec["baseline"]
folds_json = json.load(open(FOLDS))
counts = {}
for fold in range(5):
    els = load_peaks(a.elements, folds_json[str(fold)]["val"])
    # One accessibility spec is enough: 2.15 asserts every spec yields the same valid mask,
    # so the baseline's reproduces the element set the metrics were computed on.
    _, sigs, _, _ = extract_windows(
        els, GEN, b["signal_plus_bw"], b.get("signal_minus_bw"),
        b["accessibility_bw"], IN_W_MAX, OUT_W, 0, is_peak=True)
    obs = np.log1p(sigs.sum(axis=(1, 2)))
    top = obs >= np.quantile(obs, 0.8)
    counts[fold] = (len(obs), int(top.sum()))
    print(f"  fold{fold}: n={len(obs):,}  n_topq={int(top.sum()):,} "
          f"({100 * top.sum() / len(obs):.2f}%, n/5 would be {len(obs) / 5:,.0f})",
          flush=True)
    del sigs, obs, top

stored = df.groupby("fold")["n"].agg(["min", "max"])
bad = [f for f in counts if not (stored.loc[f, "min"] == stored.loc[f, "max"] == counts[f][0])]
if bad:
    for f in bad:
        print(f"fold{f}: stored n={stored.loc[f, 'min']}..{stored.loc[f, 'max']}, "
              f"recomputed {counts[f][0]}", file=sys.stderr)
    raise SystemExit("error: recomputed element counts disagree with the stored `n`; "
                     "the config/elements pair does not reproduce the scored set. "
                     "Nothing written.")

df.insert(df.columns.get_loc("n") + 1, "n_topq", df["fold"].map(lambda f: counts[f][1]))
if a.dry_run:
    print("\n--dry-run: not writing")
else:
    df.to_csv(table, sep="\t", index=False)
    print(f"\nwrote {table} (+n_topq)")
