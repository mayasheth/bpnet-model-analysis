#!/usr/bin/env python3
"""Config for scoring the K562 hold-out of the multi-cell p300 panel.

THE QUESTION: does training on four cell types beat training on one, when the target cell
type was seen by neither? Everything is scored on K562 p300 over K562's own EP300 peaks,
with K562 ATAC as the accessibility input, so the arms differ only in what they were
trained on.

BASELINE IS K562's OWN ACCESSIBILITY-ONLY p300 MODEL, which is the assay-matched floor
F-016 requires: the same assay as the models it bounds, in the target cell type. Using an
H3K27ac-derived or cross-assay floor is the error F-016 caught and F-010 was withdrawn for.

The two single-source GM12878 arms are both included on purpose. The corrupted one is what
F-021 benchmarked; the corrected one is the F-025 retrain. Keeping both means the
multi-cell comparison is readable against either.
"""
import json, os, sys

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
Q = f"{D}/2026_0529_multimodal_p300_model"
SIG_P = f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig"
SIG_M = f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"
ACC = f"{Q}/data/atac.bw"

# K562 ATAC over THESE regions (the EP300 peaks), from 4.26. Not the ABC-candidate-region
# value of 3.7232/1.1428: EP300 peaks are the more accessible subset, and the statistic has
# to describe the windows actually scored.
TGT_MEAN, TGT_STD = 4.2275, 1.3051


def e(label, model_dir, mode="multimodal", target_norm=False):
    d = {"label": label, "mode": mode, "half_window": 500, "model_dir": model_dir,
         "accessibility_bw": ACC, "signal_plus_bw": SIG_P, "signal_minus_bw": SIG_M}
    if target_norm:
        d["acc_mean"], d["acc_std"] = TGT_MEAN, TGT_STD
    return d

cfg = {
    "baseline": e("K562_atac_only_floor", f"{Q}/models/atac_only", mode="atac"),
    "compare": [
        e("K562_incell",            f"{Q}/models/atac"),
        e("GM_single_corrupted",    f"{D}/2026_0606_GM12878_transferability/"
                                    f"GM12878_multimodal_BPNet/models/atac"),
        e("GM_single_corrected",    f"{P}/models/gm12878_p300_corrected"),
        # POLICY-MATCHED COMPARATOR. The multi-cell model has no single training statistic,
        # so it MUST be normalized on the target. The single-source models have one and
        # conventionally use it, which is how F-021 benchmarked them. Comparing the two
        # directly would therefore mix "trained on four cell types" with "normalized
        # differently". This arm is the corrected single-source model scored under the
        # multi-cell model's normalization policy, so the pair against it isolates the
        # training set.
        e("GM_single_corrected_tgtnorm", f"{P}/models/gm12878_p300_corrected",
          target_norm=True),
        e("multicell_holdK562",     f"{P}/models/p300_multicell_holdK562_hw500_clw10",
          target_norm=True),
    ],
}
missing = []
for entry in [cfg["baseline"]] + cfg["compare"]:
    for k in ("model_dir", "accessibility_bw", "signal_plus_bw", "signal_minus_bw"):
        if not os.path.exists(entry[k]):
            missing.append(f"{entry['label']}: {k} -> {entry[k]}")
if missing:
    print("MISSING:\n  " + "\n  ".join(missing), file=sys.stderr)
    raise SystemExit(1)

out = f"{P}/config/p300_multicell_k562_configs.json"
with open(out, "w") as f:
    json.dump(cfg, f, indent=2)
print(f"wrote {out}")
print(f"  baseline: {cfg['baseline']['label']}")
for c in cfg["compare"]:
    print(f"  compare : {c['label']}")
