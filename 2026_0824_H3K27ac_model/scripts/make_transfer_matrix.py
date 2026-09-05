#!/usr/bin/env python3
"""Transfer matrix: rank the candidate architectures by CROSS-CELL-TYPE performance.

The in-cell-type ranking is not the deployment ranking. Wide gained +0.028 in K562 and
+0.016 in GM12878 in-cell-type, but only +0.012 and +0.005 on transfer, neither significant.
Fragment channels have never been tested on transfer at all, and they are the arm most at
risk: fragment-size distributions are library properties (depth, size selection), so those
five channels may carry batch signal that does not travel.

Each table holds the transferred arms AND the target cell type's own models, so the transfer
drop reads off one table without cross-referencing. Baseline is the target's own ATAC-only
model in both directions, identical across arms, so it cancels from the paired tests.

fragchan_wide exists for K562 only, so the GM12878->K562 direction has three arms not four.
"""
import json, os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"
CF, M = f"{K27}/config", f"{K27}/models"

K5 = dict(acc=f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
          sp=f"{K27}/data/h3k27ac_5p_plus.bw", sm=f"{K27}/data/h3k27ac_5p_minus.bw",
          frag=",".join([f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw"] +
                        [f"{K27}/data/atac_{k}5p.bw" for k in ("sub","mono","di","poly")]))
GM = dict(acc=f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw",
          sp=f"{K27}/data/gm12878_h3k27ac_5p_plus.bw",
          sm=f"{K27}/data/gm12878_h3k27ac_5p_minus.bw",
          frag=",".join([f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw"] +
                        [f"{K27}/data/gm12878_atac_{k}5p.bw" for k in ("sub","mono","di","poly")]))

# model dirs by (training cell, architecture); None where the model does not exist
ARCH = {
    ("k562", "narrow_flat"): f"{M}/multimodal5p_accs5p_hw500_clw10",
    ("k562", "WIDE_flat"):   f"{M}/multimodal5p_accs5p_wide_hw500_clw10",
    ("k562", "narrow_FRAG"): f"{M}/multimodal5p_fragchan_hw500_clw10",
    ("k562", "WIDE_FRAG"):   f"{M}/multimodal5p_fragchan_wide_hw500_clw10",
    ("gm12878", "narrow_flat"): f"{M}/gm12878_multimodal5p_accs5p_hw500_clw10",
    ("gm12878", "WIDE_flat"):   f"{M}/gm12878_multimodal5p_accs5p_wide_hw500_clw10",
    ("gm12878", "narrow_FRAG"): f"{M}/gm12878_multimodal5p_fragchan_hw500_clw10",
}
LOCAL_ATAC = {"k562": f"{M}/atac5p_accs5p_hw500_clw10",
              "gm12878": f"{M}/gm12878_atac5p_accs5p_hw500_clw10"}
CELL = {"k562": K5, "gm12878": GM}


def e(label, mode, md, tgt, frag=False):
    return {"label": label, "mode": mode, "half_window": 500, "model_dir": md,
            "accessibility_bw": tgt["frag"] if frag else tgt["acc"],
            "signal_plus_bw": tgt["sp"], "signal_minus_bw": tgt["sm"]}


bad = []
for src, dst in (("k562", "gm12878"), ("gm12878", "k562")):
    tgt = CELL[dst]
    compare = []
    for arch in ("narrow_flat", "WIDE_flat", "narrow_FRAG", "WIDE_FRAG"):
        frag = arch.endswith("FRAG")
        md_src, md_dst = ARCH.get((src, arch)), ARCH.get((dst, arch))
        if md_src:
            compare.append(e(f"{arch}_transferred", "multimodal", md_src, tgt, frag))
        if md_dst:
            compare.append(e(f"{arch}_local", "multimodal", md_dst, tgt, frag))
    spec = {"baseline": e("atac_LOCAL", "atac", LOCAL_ATAC[dst], tgt), "compare": compare}
    out = f"{CF}/txmatrix_{src}_to_{dst}_configs.json"
    json.dump(spec, open(out, "w"), indent=1)
    print(f"\n{os.path.basename(out)}")
    for x in [spec["baseline"]] + compare:
        n = sum(os.path.exists(f"{x['model_dir']}/fold{i}/training_complete.json")
                for i in range(5))
        print(f"   {x['label']:<24} {n}/5  {len(x['accessibility_bw'].split(','))}ch  "
              f"{os.path.basename(x['model_dir'])}")
        if n < 5:
            bad.append(f"{x['label']} {n}/5")
        for path in x["accessibility_bw"].split(","):
            if not os.path.exists(path):
                bad.append(f"{x['label']}: missing {path}")
print("\n" + ("ALL READY" if not bad else "PROBLEMS: " + "; ".join(bad)))
