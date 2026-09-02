#!/usr/bin/env python3
"""Emit 2.15 configs for the wide-receptive-field and fragment-channel comparisons.

The baseline in every config is the NARROW ATAC-only model. It is the same model in both
arms of each paired test, so it cancels out of the comparison and only sets the reference
for residual_pearson and incremental_r2.

`sequence_narrow` points at sequence5p_hw500_clw10, not an _accs5p directory: sequence mode
never reads the accessibility input, so no _accs5p sequence model exists or is needed.
"""
import json, os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"
M = f"{K27}/models"
CFG = f"{K27}/config"

CELLS = {
    "k562": dict(pfx="", acc=f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
                 sp=f"{K27}/data/h3k27ac_5p_plus.bw",
                 sm=f"{K27}/data/h3k27ac_5p_minus.bw"),
    "gm12878": dict(pfx="gm12878_", acc=f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw",
                    sp=f"{K27}/data/gm12878_h3k27ac_5p_plus.bw",
                    sm=f"{K27}/data/gm12878_h3k27ac_5p_minus.bw"),
}

FRAG = ",".join([f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw"] +
                [f"{K27}/data/atac_{k}5p.bw" for k in ("sub", "mono", "di", "poly")])


def entry(label, mode, model_dir, c, acc=None):
    return {"label": label, "mode": mode, "half_window": 500,
            "model_dir": model_dir,
            "accessibility_bw": acc or c["acc"],
            "signal_plus_bw": c["sp"], "signal_minus_bw": c["sm"]}


for cell, c in CELLS.items():
    p = c["pfx"]
    spec = {
        "baseline": entry("atac_narrow", "atac", f"{M}/{p}atac5p_accs5p_hw500_clw10", c),
        "compare": [
            entry("sequence_narrow", "sequence", f"{M}/{p}sequence5p_hw500_clw10", c),
            entry("sequence_WIDE", "sequence",
                  f"{M}/{p}sequence5p_accs5p_wide_hw500_clw10", c),
            entry("multimodal_narrow", "multimodal",
                  f"{M}/{p}multimodal5p_accs5p_hw500_clw10", c),
            entry("multimodal_WIDE", "multimodal",
                  f"{M}/{p}multimodal5p_accs5p_wide_hw500_clw10", c),
        ],
    }
    out = f"{CFG}/wide_{cell}_configs.json"
    json.dump(spec, open(out, "w"), indent=1)
    print("wrote", out)

# Fragment channels: K562 only, since fragment length needs the PE BAMs.
c = CELLS["k562"]
spec = {
    "baseline": entry("atac_flat", "atac", f"{M}/atac5p_accs5p_hw500_clw10", c),
    "compare": [
        entry("multimodal_flat", "multimodal", f"{M}/multimodal5p_accs5p_hw500_clw10", c),
        entry("multimodal_FRAGCHAN", "multimodal",
              f"{M}/multimodal5p_fragchan_hw500_clw10", c, acc=FRAG),
    ],
}
out = f"{CFG}/fragchan_k562_configs.json"
json.dump(spec, open(out, "w"), indent=1)
print("wrote", out)

for f in ("wide_k562_configs.json", "wide_gm12878_configs.json",
          "fragchan_k562_configs.json"):
    s = json.load(open(f"{CFG}/{f}"))
    for e in [s["baseline"]] + s["compare"]:
        for k in ("accessibility_bw", "signal_plus_bw", "signal_minus_bw"):
            for path in e[k].split(","):
                assert os.path.exists(path), f"{f}: {e['label']}: missing {path}"
    print(f"  {f}: all input files present")
