#!/usr/bin/env python3
"""Generate the combined wide-receptive-field + fragment-channel script.

Both changes gained on the top quintile separately -- +0.027 for n_layers 10 and +0.0135 for
the fragment channels -- and both act on the accessibility side, so they may be reading the
same neighbourhood structure and may not sum. This is the run that decides what the deployed
model should be.

Derived from 1.15 (fragment channels) with EXACTLY THE THREE WIDE EDITS used for 1.13/1.14:
n_layers, the derived trimming sum, and the output directory. The accessibility input stays
the 5-channel fragment set, so the two comparison arms are
  models/multimodal5p_fragchan_hw500_clw10            (narrow + fragment channels)
  models/multimodal5p_accs5p_wide_hw500_clw10          (wide  + flat channel)
and the new model differs from each by exactly one of the two changes.
"""
import os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"

HDR = """#
# ===== AUTO-DERIVED: fragment channels PLUS the ~4.2 kb receptive field =====
# Generated from {src} by scripts/make_fragwide_scripts.py.
# EXACTLY THREE changes vs the original: n_layers 8 -> 10, the trimming sum extended to
# match, and the output directory. The 5-channel fragment accessibility input is unchanged.
#
# Tests whether the two accessibility-side gains combine. Separately, on the K562 top
# quintile: n_layers 10 gave +0.027 and the fragment channels gave +0.0135, both against
# multimodal5p_accs5p_hw500_clw10. If they read the same neighbourhood structure they will
# not sum.
#
"""

EDITS = [
    ("N_LAYERS=8", "N_LAYERS=10"),
    ("TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))",
     "TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 + 512 + 1024 ))"),
    ('OUT_DIR="$PROJ/models/${MODE}5p_fragchan_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
     'OUT_DIR="$PROJ/models/${MODE}5p_fragchan_wide_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"'),
]

src = f"{K27}/scripts/1.15.submit_training_fragchan.sh"
dst = f"{K27}/scripts/1.16.submit_training_fragchan_wide.sh"
with open(src) as f:
    text = f.read()
for old, new in EDITS:
    n = text.count(old)
    assert n == 1, f"expected 1 occurrence of {old!r}, found {n}"
    text = text.replace(old, new)
lines = text.split("\n")
assert lines[0].startswith("#!")
text = lines[0] + "\n" + HDR.format(src=os.path.basename(src)) + "\n".join(lines[1:])
with open(dst, "w") as f:
    f.write(text)
os.chmod(dst, 0o755)
print("wrote", dst)
