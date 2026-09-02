#!/usr/bin/env python3
"""Generate wide-receptive-field training scripts by copying the accs5p submit scripts and
changing EXACTLY THREE things: n_layers, the derived trimming sum, and the output directory.

WHY. H3K27ac signal is still at 28% of maximum at +/-2 kb, but the 8-layer model sees only
47 + sum(2^i, i=1..8) = 557 bp either side of each output position (~1.1 kb receptive
field). Ten layers give 47 + sum(2^i, i=1..10) = 2093 bp either side (~4.2 kb), which
covers the range over which the mark is actually spread.

Geometry follows from n_layers alone:
  n_layers=8   trimming=557    in_window=1000+2*557 =2114
  n_layers=10  trimming=2093   in_window=1000+2*2093=5186
train_multimodal_bpnet.py asserts in_window - 2*model.trimming == out_window, so a
mismatch fails immediately rather than training something subtly wrong.

Everything else -- count_loss_weight 10, n_filters 64, n_acc_filters 8, max_jitter 50,
batch size 64, out_window 1000, elements, negatives, folds, and the 5' ATAC accessibility
input -- is byte-identical to the accs5p runs, so those are the matched baseline.

CAVEAT for evaluation: a 5186 bp input window drops regions within 2593 bp of a chromosome
end, where 2114 bp kept them. Compare the wide and narrow models on the INTERSECTION of
their valid regions, not on each model's own region set.
"""
import os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"

HDR = """#
# ===== AUTO-DERIVED: receptive field widened from ~1.1 kb to ~4.2 kb =====
# Generated from {src} by scripts/make_wide_scripts.py.
# EXACTLY THREE changes vs the original: n_layers 8 -> 10, the trimming sum extended to
# match, and the output directory. Every other hyperparameter is byte-identical, so a
# difference in results is attributable to the receptive field and nothing else.
#
# n_layers=8   trimming=557   in_window=2114  (~1.1 kb receptive field)
# n_layers=10  trimming=2093  in_window=5186  (~4.2 kb receptive field)
#
# Accessibility input stays the ChromBPNet 5' insertion track, so the matched baseline is
# the *_accs5p_hw500_clw10 model of the same mode, cell type and fold.
#
"""

EDITS = [
    ("N_LAYERS=8", "N_LAYERS=10"),
    ("TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))",
     "TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 + 512 + 1024 ))"),
]

JOBS = [
    (f"{K27}/scripts/1.11.submit_training_5prime_accs5p.sh",
     f"{K27}/scripts/1.13.submit_training_5prime_wide.sh",
     EDITS + [('OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
               'OUT_DIR="$PROJ/models/${MODE}5p_accs5p_wide_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"')]),

    (f"{K27}/scripts/1.12.submit_training_gm12878_accs5p.sh",
     f"{K27}/scripts/1.14.submit_training_gm12878_wide.sh",
     EDITS + [('OUT_DIR="$PROJ/models/gm12878_${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
               'OUT_DIR="$PROJ/models/gm12878_${MODE}5p_accs5p_wide_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"')]),
]

for src, dst, edits in JOBS:
    with open(src) as f:
        text = f.read()
    for old, new in edits:
        n = text.count(old)
        assert n == 1, f"{src}: expected 1 occurrence of {old!r}, found {n}"
        text = text.replace(old, new)
    # insert the provenance header immediately after the shebang
    lines = text.split("\n")
    assert lines[0].startswith("#!"), src
    text = lines[0] + "\n" + HDR.format(src=os.path.basename(src)) + "\n".join(lines[1:])
    with open(dst, "w") as f:
        f.write(text)
    os.chmod(dst, 0o755)
    print("wrote", dst)
