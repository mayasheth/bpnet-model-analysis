#!/usr/bin/env python3
"""Compare two 2.15 per-fold tables as a regression check on the evaluator.

Byte-identity is the WRONG invariant here. cuDNN convolution is not bit-reproducible across
GPU models, so re-scoring the same models on a different node moves metrics in the 4th
decimal place even when the code is untouched. A `diff` gate therefore fails for reasons
that have nothing to do with the change under test.

Split the invariant by what is actually deterministic:

  EXACT   the set of folds, the set of config labels, and n per fold. These come from window
          extraction and bounds checking, which is pure numpy on file contents. This is
          precisely what a geometry or region-set change would break, so it must match
          exactly.
  TOLERANT  every metric, to --tol (default 1e-3). Well below this project's between-fold sd
          of 0.041-0.046 and below any effect size under discussion, so a real behavioural
          change -- scoring the wrong input, cropping off-centre, misaligning rows -- moves
          numbers by far more than the float noise this absorbs.

Usage: 2.18.compare_perfold_tables.py NEW.tsv REFERENCE.tsv [--tol 1e-3]
Exit 0 on pass, 1 on fail.
"""
import argparse, sys
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("new"); ap.add_argument("ref")
ap.add_argument("--tol", type=float, default=1e-3)
a = ap.parse_args()

new = pd.read_csv(a.new, sep="\t")
ref = pd.read_csv(a.ref, sep="\t")
fail = []

if list(new.columns) != list(ref.columns):
    fail.append(f"columns differ: {list(new.columns)} vs {list(ref.columns)}")
    print("FAIL:", fail[0]); sys.exit(1)

key = ["fold", "config"]
new = new.sort_values(key).reset_index(drop=True)
ref = ref.sort_values(key).reset_index(drop=True)

for c in key:
    if not new[c].equals(ref[c]):
        fail.append(f"{c} column differs between tables")
if new.shape != ref.shape:
    fail.append(f"shape differs: {new.shape} vs {ref.shape}")

if not fail and not new["n"].equals(ref["n"]):
    bad = new.loc[new["n"] != ref["n"], key + ["n"]]
    fail.append("n (region count) differs -- the scored region set changed:\n"
                + bad.to_string(index=False))

metrics = [c for c in new.columns if c not in key + ["n"]]
worst = 0.0
if not fail:
    for m in metrics:
        d = (new[m].to_numpy(float) - ref[m].to_numpy(float))
        d = np.abs(np.nan_to_num(d, nan=0.0))
        i = int(d.argmax())
        if d[i] > worst:
            worst = d[i]
        if d[i] > a.tol:
            fail.append(f"{m}: max |diff| {d[i]:.6f} > tol {a.tol:g} at "
                        f"fold{new.loc[i,'fold']}/{new.loc[i,'config']} "
                        f"({new.loc[i,m]} vs {ref.loc[i,m]})")

print(f"compared {len(new)} rows x {len(metrics)} metrics")
print(f"region counts per fold: {'IDENTICAL' if not any('n (' in f for f in fail) else 'DIFFER'}")
print(f"largest metric difference: {worst:.6f}  (tolerance {a.tol:g})")
if fail:
    print("\nREGRESSION FAIL")
    for f in fail:
        print("  -", f)
    sys.exit(1)
print("\nREGRESSION PASS: region sets identical, metrics within float/GPU noise")
