#!/usr/bin/env python3
"""Training scripts for the gate x loss factorial.

Four arms against the existing multimodal5p_accs5p baseline:

  gate         sequence gates the accessibility branch, symmetric loss
  asym         symmetric architecture, over-prediction weighted 3x
  gate_asym    both -- the intended experiment
  (indicator)  built separately by 1.22, since it needs new input tracks

`gate` alone is expected to do little: the gate has no gradient pressure to close while the
loss barely charges for over-prediction. That is a claim, so it gets its own arm rather than
an assumption. `asym` alone isolates how much of any gain is the loss rather than the gate.

Derived from 1.11 with exactly two edits per arm (the training flags, and the output
directory), so the baseline is the same 5 folds already trained.
"""
import os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"
SRC = f"{K27}/scripts/1.11.submit_training_5prime_accs5p.sh"

HDR = """#
# ===== AUTO-DERIVED: {desc} =====
# Generated from 1.11.submit_training_5prime_accs5p.sh by scripts/make_gate_scripts.py.
# EXACTLY TWO changes vs the original: the training flags below, and the output directory.
#
# Motivation: every model that sees ATAC over-predicts H3K27ac at accessible-but-
# unacetylated elements (CpG-island promoters, CTCF sites) by 7-8x relative to its own
# genome-wide median, while the sequence-only model elevates them only 1.6x and observed
# H3K27ac not at all. The signal is in sequence; the additive trunk gives it no way to veto
# accessibility. And log1pMSE is ~400x more sensitive to missing a strong enhancer than to
# inventing a weak one, so nothing pushes the model to fix it.
#
"""

ARMS = [
    ("gate", "sequence gates accessibility, symmetric loss",
     "--gate-accessibility"),
    ("asym", "over-prediction weighted 3x, ungated",
     "--overprediction-weight 3.0"),
    ("gate_asym", "sequence gates accessibility AND over-prediction weighted 3x",
     "--gate-accessibility --overprediction-weight 3.0"),
]

with open(SRC) as f:
    base = f.read()

ANCHOR = "    --negative-ratio 0.1"
OUT_OLD = ('OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}'
           '_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"')

for tag, desc, flags in ARMS:
    text = base
    assert text.count(ANCHOR) == 1
    text = text.replace(ANCHOR, ANCHOR + " \\\n    " + flags)
    out_new = ('OUT_DIR="$PROJ/models/${MODE}5p_%s_hw${HALF_WINDOW}'
               '_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"' % tag)
    assert text.count(OUT_OLD) == 1
    text = text.replace(OUT_OLD, out_new)
    lines = text.split("\n")
    assert lines[0].startswith("#!")
    text = lines[0] + "\n" + HDR.format(desc=desc) + "\n".join(lines[1:])
    dst = f"{K27}/scripts/1.19.submit_training_{tag}.sh"
    with open(dst, "w") as f:
        f.write(text)
    os.chmod(dst, 0o755)
    print(f"wrote {os.path.basename(dst)}  ->  models/multimodal5p_{tag}_hw500_clw10")
