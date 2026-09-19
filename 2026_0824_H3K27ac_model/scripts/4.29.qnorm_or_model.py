#!/usr/bin/env python3
"""Is F-021's accessibility-slope error a property of the models, or of ABC's qnorm?

F-021 measured the ratio of transferred to in-cell predicted activity by ATAC decile and found
it ran 0.64 to 2.16, calling it a monotone accessibility-slope error in the transferred model.
That was measured on `normalized_h3k27ac_enh` in the ABC putative predictions, restricted to
CRISPR-tested element-gene pairs. `4.28` then measured the same ratio on the RAW predicted
bigwigs over all candidate regions and got 0.475 to 0.516, a spread of 1.20x rather than 3.4x.

Both cannot describe the same thing, so something between the bigwig and the putative table
creates the slope. Two candidates, and this separates them:
  QNORM         ABC rank-normalizes each arm's activity onto a reference distribution, so two
                tracks with differently shaped marginals are remapped differently.
  CRISPR SUBSET the tested element-gene pairs are not a random sample of candidate regions;
                they are biased toward accessible, gene-proximal elements.

Reading the post-qnorm values over ALL elements distinguishes them. If the 3.4x swing survives
over all elements, qnorm produces it. If it disappears, the CRISPR restriction produces it.
Either way F-021's mechanism sentence needs rewriting, because it attributes the slope to what
the model learned.
"""
import sys
import numpy as np
import pandas as pd

A = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
PF = "Predictions/EnhancerPredictionsAllPutative.tsv.gz"
ARMS = {
    "in_cell":    f"{A}/2026_0905_p300_activity/p300pred_k562_multimodal/{PF}",
    "transferred": f"{A}/2026_0906_p300_transfer/p300pred_gm12878_multimodal/{PF}",
}
KEEP = ["chr", "start", "end", "normalized_atac_enh", "normalized_h3k27ac_enh"]

el = {}
for label, path in ARMS.items():
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=KEEP, chunksize=1_000_000,
                          low_memory=False):
        parts.append(ch)
    d = pd.concat(parts)
    d["ekey"] = (d["chr"].astype(str) + ":" + d["start"].astype(str) + "-"
                 + d["end"].astype(str))
    # normalized_*_enh are element-level, repeated across that element's gene pairs.
    d = d.drop_duplicates("ekey").set_index("ekey")
    el[label] = d
    print(f"{label:<12} {len(d):,} distinct elements")

keys = sorted(set(el["in_cell"].index) & set(el["transferred"].index))
t = pd.DataFrame({
    "atac": el["in_cell"].loc[keys, "normalized_atac_enh"].to_numpy(),
    "in_cell": el["in_cell"].loc[keys, "normalized_h3k27ac_enh"].to_numpy(),
    "transferred": el["transferred"].loc[keys, "normalized_h3k27ac_enh"].to_numpy(),
}, index=keys).dropna()
print(f"\n{len(t):,} elements shared by both arms\n")

t["dec"] = pd.qcut(t["atac"].rank(method="first"), 10, labels=False) + 1
rows = []
for dec, g in t.groupby("dec"):
    rows.append({"dec": int(dec), "n": len(g), "med_atac": g["atac"].median(),
                 "med_in_cell": g["in_cell"].median(),
                 "med_transferred": g["transferred"].median(),
                 "tx_over_in": g["transferred"].median() / max(g["in_cell"].median(), 1e-12)})
out = pd.DataFrame(rows)
print("POST-QNORM ACTIVITY, ALL ELEMENTS (not just CRISPR-tested pairs)")
print(out.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))
v = out["tx_over_in"].to_numpy()
print(f"\ntx_over_in: decile 1 {v[0]:.3f} -> decile 10 {v[-1]:.3f}, "
      f"spread max/min x{v.max() / max(v.min(), 1e-12):.2f}")
print()
print("For reference, the same ratio measured two other ways:")
print("  raw bigwigs, all candidate regions (4.28)        0.475 -> 0.516, spread x1.20")
print("  post-qnorm, CRISPR-tested pairs only (F-021)     0.64  -> 2.16,  spread x3.40")
print()
# Decompose rather than pick a winner. The first version of this script printed a binary
# verdict on a >2x threshold and mislabelled a result that is roughly half and half.
raw, allq, crispr = 1.20, v.max() / max(v.min(), 1e-12), 3.40
import math
lr, la, lc = math.log(raw), math.log(allq), math.log(crispr)
print("DECOMPOSITION of the 3.4x swing, on a log scale, since these compose multiplicatively:")
print(f"  the models themselves       x{raw:.2f}   {100 * lr / lc:4.0f}% of it")
print(f"  ABC's qnorm adds            x{allq / raw:.2f}   {100 * (la - lr) / lc:4.0f}%")
print(f"  the CRISPR subset adds      x{crispr / allq:.2f}   {100 * (lc - la) / lc:4.0f}%")
print()
print("So the accessibility slope is mostly MANUFACTURED DOWNSTREAM of the model, by rank")
print("normalization against a reference distribution and then by scoring on a biased,")
print("accessible, gene-proximal subset of elements. F-021's numbers are correct for the")
print("quantity ABC uses where the benchmark looks, but its mechanism sentence attributes")
print("them to what the transferred model learned, and that is not supported.")
print()
print("CONSEQUENCE FOR THE PROPOSED FIX: a monotone recalibration fitted on raw predictions")
print("in the source cell type cannot repair a defect the raw predictions barely have.")
