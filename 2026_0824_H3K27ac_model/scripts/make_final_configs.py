#!/usr/bin/env python3
"""Configs for the four comparisons whose models all finished on 2026-09-05.

  1. atac_WIDE          does the receptive-field gain need sequence at all?
  2. fragchan x wide     do the two accessibility-side gains combine, or read the same thing?
  3. ATAC elements       does training on the element set ABC scores us on help?
  4. GM12878 fragchan    does the fragment-channel gain replicate in a second cell type?

All four are scored with --no-rc-average (see 2.23) so they are directly comparable to the
existing single-pass figures. RC averaging is the 2.15 default now, and 2.22 converts every
table together; mixing the two inside one report is the failure mode to avoid.
"""
import json, os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"
CF, M = f"{K27}/config", f"{K27}/models"

K5 = dict(acc=f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
          sp=f"{K27}/data/h3k27ac_5p_plus.bw", sm=f"{K27}/data/h3k27ac_5p_minus.bw")
GM = dict(acc=f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw",
          sp=f"{K27}/data/gm12878_h3k27ac_5p_plus.bw",
          sm=f"{K27}/data/gm12878_h3k27ac_5p_minus.bw")
FRAG_K5 = ",".join([K5["acc"]] + [f"{K27}/data/atac_{k}5p.bw"
                                  for k in ("sub", "mono", "di", "poly")])
FRAG_GM = ",".join([GM["acc"]] + [f"{K27}/data/gm12878_atac_{k}5p.bw"
                                  for k in ("sub", "mono", "di", "poly")])


def e(label, mode, model_dir, c, acc=None):
    return {"label": label, "mode": mode, "half_window": 500, "model_dir": model_dir,
            "accessibility_bw": acc or c["acc"],
            "signal_plus_bw": c["sp"], "signal_minus_bw": c["sm"]}


specs = {}

# 1. add the ATAC-only wide arm to both wide configs
for cell, c, p in (("k562", K5, ""), ("gm12878", GM, "gm12878_")):
    specs[f"wide_{cell}_configs.json"] = {
        "baseline": e("atac_narrow", "atac", f"{M}/{p}atac5p_accs5p_hw500_clw10", c),
        "compare": [
            e("atac_WIDE", "atac", f"{M}/{p}atac5p_accs5p_wide_hw500_clw10", c),
            e("sequence_narrow", "sequence", f"{M}/{p}sequence5p_hw500_clw10", c),
            e("sequence_WIDE", "sequence", f"{M}/{p}sequence5p_accs5p_wide_hw500_clw10", c),
            e("multimodal_narrow", "multimodal", f"{M}/{p}multimodal5p_accs5p_hw500_clw10", c),
            e("multimodal_WIDE", "multimodal", f"{M}/{p}multimodal5p_accs5p_wide_hw500_clw10", c),
        ]}

# 2. the 2x2: receptive field x fragment channels
specs["fragwide_k562_configs.json"] = {
    "baseline": e("atac_flat", "atac", f"{M}/atac5p_accs5p_hw500_clw10", K5),
    "compare": [
        e("mm_narrow_flat", "multimodal", f"{M}/multimodal5p_accs5p_hw500_clw10", K5),
        e("mm_narrow_FRAG", "multimodal", f"{M}/multimodal5p_fragchan_hw500_clw10", K5, FRAG_K5),
        e("mm_WIDE_flat", "multimodal", f"{M}/multimodal5p_accs5p_wide_hw500_clw10", K5),
        e("mm_WIDE_FRAG", "multimodal", f"{M}/multimodal5p_fragchan_wide_hw500_clw10", K5, FRAG_K5),
    ]}

# 3. element derivation, scored on the ATAC-derived set both models will face in ABC
specs["atacel_k562_configs.json"] = {
    "baseline": e("atac_flat", "atac", f"{M}/atac5p_accs5p_hw500_clw10", K5),
    "compare": [
        e("mm_trained_DNase_elements", "multimodal", f"{M}/multimodal5p_accs5p_hw500_clw10", K5),
        e("mm_trained_ATAC_elements", "multimodal", f"{M}/multimodal5p_atacel_hw500_clw10", K5),
    ]}

# 4. fragment channels in GM12878
specs["fragchan_gm12878_configs.json"] = {
    "baseline": e("atac_flat", "atac", f"{M}/gm12878_atac5p_accs5p_hw500_clw10", GM),
    "compare": [
        e("multimodal_flat", "multimodal", f"{M}/gm12878_multimodal5p_accs5p_hw500_clw10", GM),
        e("multimodal_FRAGCHAN", "multimodal", f"{M}/gm12878_multimodal5p_fragchan_hw500_clw10", GM, FRAG_GM),
    ]}

bad = []
for name, spec in specs.items():
    json.dump(spec, open(f"{CF}/{name}", "w"), indent=1)
    print(f"\n{name}")
    for x in [spec["baseline"]] + spec["compare"]:
        n = sum(os.path.exists(f"{x['model_dir']}/fold{i}/training_complete.json")
                for i in range(5))
        nch = len(x["accessibility_bw"].split(","))
        print(f"   {x['label']:<26} {n}/5 folds  {nch}ch  {os.path.basename(x['model_dir'])}")
        if n < 5:
            bad.append(f"{name}:{x['label']} {n}/5")
        for path in x["accessibility_bw"].split(",") + [x["signal_plus_bw"], x["signal_minus_bw"]]:
            if not os.path.exists(path):
                bad.append(f"{name}:{x['label']}: missing {path}")
print("\n" + ("ALL READY" if not bad else "PROBLEMS: " + "; ".join(bad)))
