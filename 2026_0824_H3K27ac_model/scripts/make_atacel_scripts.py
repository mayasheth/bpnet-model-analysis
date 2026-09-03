#!/usr/bin/env python3
"""Generate the ATAC-derived-element training script: EXACTLY TWO changes from 1.11.

Element derivation is the one variable. Everything else -- 5' target, 5' ATAC accessibility
input, count_loss_weight 10, n_layers 8, window geometry, negatives, folds -- is byte
identical, so a difference is attributable to which assay called the regions.

Why it matters beyond tidiness: the ABC experiment scores our models on the ATAC-derived
regions while they were trained on the DNase-derived ones. Scoring BOTH models on the
ATAC-derived set answers whether closing that mismatch is worth anything.
"""
import os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"

HDR = """#
# ===== AUTO-DERIVED: elements switched from DNase-called to ATAC-called =====
# Generated from {src} by scripts/make_atacel_scripts.py.
# EXACTLY TWO changes vs the original: the element file, and the output directory.
#
# Old: reference/K562_DNase_candidate_elements.narrowPeak       (150,528 DNase-called)
# New: reference/K562_ATAC_candidate_elements.narrowPeak        (153,545 ATAC-called)
# The new file is the rE2G ATAC_H3K27ac_powerlaw candidate regions, byte-identical to the
# candidate regions ABC scores, so this model is trained on the regions it will be evaluated
# on. Both files centre windows on the region midpoint, so centring is unchanged.
#
"""

EDITS = [
    ('ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"',
     'ELEMENTS="$PROJECT_DIR/reference/K562_ATAC_candidate_elements.narrowPeak"'),
    ('OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
     'OUT_DIR="$PROJ/models/${MODE}5p_atacel_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"'),
]

src = f"{K27}/scripts/1.11.submit_training_5prime_accs5p.sh"
dst = f"{K27}/scripts/1.17.submit_training_atac_elements.sh"
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
