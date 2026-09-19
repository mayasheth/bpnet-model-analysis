#!/usr/bin/env python3
"""Does changing the predicted track's scale survive ABC's qnorm? Structurally it should not.

`run_qnorm(qnorm_method="rank")` in ABC's neighborhoods.py maps each element's WITHIN-ARM rank
through an interpolation onto a reference distribution, separately for promoters and
nonpromoters. So the post-qnorm activity is a function of the predicted track's per-element
RANKING and nothing else: any transformation that preserves the per-element ordering is
discarded exactly.

The 2026-09-17 accessibility renormalisation changed the p300 tracks by x1.5 (in-cell) and x1.8
(transferred) in scale while moving the per-element ranking by only 1.6% and 3.7%
(Spearman 0.9969 and 0.9853). The prediction is therefore that post-qnorm activity is nearly
unchanged and the benchmark barely moves, for reasons of ABC's architecture rather than biology.

This measures it at two levels, because they can disagree:
  post-qnorm activity   should be near-identical, since it is rank-determined
  ABC.Score             can differ more, since the geomean arms multiply the predicted track
                        against real ATAC before qnorm sees it, and contact scaling follows

A null result here is informative, not a failure: it says the input-normalisation asymmetry,
real as it is (0.73 sd against 0.38 sd), cannot be what costs the transferred arm 0.0458 AUPRC.
"""
import numpy as np
import pandas as pd

A = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
PF = "Predictions/EnhancerPredictionsAllPutative.tsv.gz"
PAIRS = [
    ("in-cell",
     "2026_0905_p300_activity/p300pred_k562_multimodal",
     "2026_0917_accnorm_activity/p300pred_k562_multimodal_accnorm"),
    ("transferred",
     "2026_0906_p300_transfer/p300pred_gm12878_multimodal",
     "2026_0917_accnorm_activity/p300pred_gm12878_multimodal_accnorm"),
]
KEEP = ["chr", "start", "end", "normalized_h3k27ac_enh", "activity_base_enh",
        "ABC.Score", "TargetGene"]


def load(path):
    parts = [c for c in pd.read_csv(f"{A}/{path}/{PF}", sep="\t", usecols=KEEP,
                                    chunksize=1_000_000, low_memory=False)]
    x = pd.concat(parts)
    x["ek"] = (x["chr"].astype(str) + ":" + x["start"].astype(str) + "-"
               + x["end"].astype(str))
    return x


for label, p_orig, p_accn in PAIRS:
    o, n = load(p_orig), load(p_accn)

    # Element level: normalized_h3k27ac_enh is element-level, repeated across gene pairs.
    e1 = o.drop_duplicates("ek").set_index("ek")["normalized_h3k27ac_enh"]
    e2 = n.drop_duplicates("ek").set_index("ek")["normalized_h3k27ac_enh"]
    k = sorted(set(e1.index) & set(e2.index))
    v1, v2 = e1.loc[k].to_numpy(), e2.loc[k].to_numpy()
    sp_act = np.corrcoef(pd.Series(v1).rank(), pd.Series(v2).rank())[0, 1]
    print(f"{label}")
    print(f"  post-qnorm activity   n={len(k):,}  identical at "
          f"{100 * np.isclose(v1, v2, rtol=1e-9, atol=1e-12).mean():.2f}% of elements  "
          f"spearman={sp_act:.6f}  median {np.median(v1):.4f} -> {np.median(v2):.4f}")

    # Pair level.
    o["pk"] = o["ek"] + "|" + o["TargetGene"].astype(str)
    n["pk"] = n["ek"] + "|" + n["TargetGene"].astype(str)
    s1 = o.drop_duplicates("pk").set_index("pk")["ABC.Score"]
    s2 = n.drop_duplicates("pk").set_index("pk")["ABC.Score"]
    kk = sorted(set(s1.index) & set(s2.index))
    a1, a2 = s1.loc[kk].to_numpy(), s2.loc[kk].to_numpy()
    print(f"  ABC.Score             n={len(kk):,}  identical at "
          f"{100 * np.isclose(a1, a2, rtol=1e-9).mean():.2f}% of pairs  "
          f"spearman={np.corrcoef(pd.Series(a1).rank(), pd.Series(a2).rank())[0, 1]:.6f}")
