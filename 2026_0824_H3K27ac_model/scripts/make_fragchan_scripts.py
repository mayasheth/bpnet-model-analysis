#!/usr/bin/env python3
"""Generate the fragment-size-stratified training script from the accs5p script, changing
EXACTLY THREE things: the accessibility input, the output directory, and the allowed modes.

Accessibility input becomes five channels instead of one:
    [all, sub, mono, di, poly]
`all` is atac_5p.bw itself and sub/mono/di/poly exhaustively partition fragments by length
(0.23.make_atac_fragment_5p_channels.sh), all in the same single-base 5' insertion
convention. The first convolution is a linear map across channels, so the model can
reproduce the accs5p input exactly by zeroing the four stratified channels -- this is a
strict superset of the baseline input and cannot regress except through optimisation noise.

n_acc_filters stays 8, so the merged representation is still 64+8=72 wide and the only
difference from the accs5p baseline is what the accessibility branch is allowed to see.

`sequence` mode is rejected: it ignores the accessibility input entirely, so it would
silently duplicate the existing sequence5p_accs5p run under a different name.
"""
import os

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
K27 = f"{D}/2026_0824_H3K27ac_model"

HDR = """#
# ===== AUTO-DERIVED: flat ATAC replaced by 5 fragment-size-stratified 5' channels =====
# Generated from {src} by scripts/make_fragchan_scripts.py.
# EXACTLY THREE changes vs the original: the accessibility input (1 channel -> 5), the
# output directory, and the rejection of `sequence` mode. Every other hyperparameter is
# byte-identical, so the matched baseline is the *_accs5p_hw500_clw10 model of the same
# mode and fold.
#
# Channels: [all, sub(<=139), mono(140-329), di(330-620), poly(>=621)] -- an exhaustive
# partition plus the flat track, every one a single-base Tn5 insertion count. Flat ATAC
# says where the chromatin is open; fragment length says where the nucleosomes are, and
# H3K27ac can only exist on a nucleosome.
#
"""

CHANS = ('ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw,'
         '$PROJ/data/atac_sub5p.bw,$PROJ/data/atac_mono5p.bw,'
         '$PROJ/data/atac_di5p.bw,$PROJ/data/atac_poly5p.bw"')

JOBS = [
    (f"{K27}/scripts/1.11.submit_training_5prime_accs5p.sh",
     f"{K27}/scripts/1.15.submit_training_fragchan.sh",
     [('ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"', CHANS),
      ('OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
       'OUT_DIR="$PROJ/models/${MODE}5p_fragchan_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"'),
      ('case "$MODE" in sequence|multimodal|atac) ;; *) echo "bad MODE \'$MODE\'" >&2; exit 1 ;; esac',
       'case "$MODE" in multimodal|atac) ;; *) echo "MODE must be multimodal or atac '
       '(sequence ignores the accessibility input)" >&2; exit 1 ;; esac')]),
]

for src, dst, edits in JOBS:
    with open(src) as f:
        text = f.read()
    for old, new in edits:
        n = text.count(old)
        assert n == 1, f"{src}: expected 1 occurrence of {old!r}, found {n}"
        text = text.replace(old, new)
    lines = text.split("\n")
    assert lines[0].startswith("#!"), src
    text = lines[0] + "\n" + HDR.format(src=os.path.basename(src)) + "\n".join(lines[1:])
    with open(dst, "w") as f:
        f.write(text)
    os.chmod(dst, 0o755)
    print("wrote", dst)
