#!/usr/bin/env python3
"""Generate ATAC-input-5' retraining scripts by copying each original submit script and
changing EXACTLY TWO things: the accessibility bigwig, and the output directory.

Everything else -- count_loss_weight, n_filters, n_acc_filters, max_jitter, batch size,
window geometry, elements, negatives, folds -- must stay byte-identical, or the comparison
measures the hyperparameters as well as the input. p300 in particular uses different
settings from the H3K27ac models (clw 1.0, n_filters 64, max_jitter 50), so copying rather
than re-authoring is the only safe route.

Naming: `_accs5p` marks the ACCESSIBILITY input as 5'-end insertion counts. Note the
existing `5p` in H3K27ac model names refers to the TARGET, which is unchanged here.
"""
import os, sys

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P300 = f"{D}/2026_0529_multimodal_p300_model"
K27 = f"{D}/2026_0824_H3K27ac_model"

HDR = """#
# ===== AUTO-DERIVED: ATAC input switched to ChromBPNet-convention 5' insertion counts =====
# Generated from {src} by scripts/make_accs5p_scripts.py.
# EXACTLY TWO changes vs the original: the accessibility bigwig, and the output directory.
# Every other hyperparameter is byte-identical, so a difference in results is attributable
# to the input definition and nothing else.
#
# Old input: genomecov -bg over the FULL tagAlign interval -> coverage scales with read
#            length (94-95 bp K562/GM12878, 35-36 bp TeloHAEC).
# New input: genomecov -bg -5 -> single-base 5' insertion counts, read-length independent,
#            matching chrombpnet/helpers/preprocessing/reads_to_bigwig.py.
# Verified: sum(old)/sum(new) equals the mean read length in all six tracks.
#
"""

JOBS = [
    # (src, dst, [(old, new), ...])
    (f"{P300}/scripts/1.1.submit_training_atac.sh",
     f"{P300}/scripts/1.1b.submit_training_atac_accs5p.sh",
     [('ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"',
       'ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"'),
      ('OUTPUT_DIR="$SCRIPT_DIR/../models/atac/fold${FOLD}"',
       'OUTPUT_DIR="$SCRIPT_DIR/../models/atac_accs5p/fold${FOLD}"')]),

    (f"{P300}/scripts/1.3.submit_training_atac_only.sh",
     f"{P300}/scripts/1.3b.submit_training_atac_only_accs5p.sh",
     [('OUTPUT_DIR=$THIS_DIR/models/atac_only/fold${FOLD}',
       'OUTPUT_DIR=$THIS_DIR/models/atac_only_accs5p/fold${FOLD}')]),

    (f"{K27}/scripts/1.4.submit_training_5prime.sh",
     f"{K27}/scripts/1.11.submit_training_5prime_accs5p.sh",
     [('ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"',
       'ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"'),
      ('OUT_DIR="$PROJ/models/${MODE}5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
       'OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"')]),

    (f"{K27}/scripts/1.8.submit_training_gm12878.sh",
     f"{K27}/scripts/1.12.submit_training_gm12878_accs5p.sh",
     [('ATAC_BW="$TRANS/data/atac.bw"', 'ATAC_BW="$TRANS/data/atac_5p.bw"'),
      ('OUT_DIR="$PROJ/models/gm12878_${MODE}5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"',
       'OUT_DIR="$PROJ/models/gm12878_${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"')]),
]

for src, dst, subs in JOBS:
    s = open(src).read()
    orig = s
    for old, new in subs:
        assert old.strip(), "empty target in " + src
        n = s.count(old)
        assert n == 1, "in %s: expected 1 occurrence of %r, found %d" % (src, old[:60], n)
        s = s.replace(old, new)
    # 1.3 references the ATAC bigwig via a differently-named variable; catch it generically
    if "atac_5p.bw" not in s:
        cand = [l for l in s.splitlines() if "data/atac.bw" in l]
        assert len(cand) == 1, "in %s: cannot locate the atac.bw reference, found %r" % (src, cand)
        s = s.replace("data/atac.bw", "data/atac_5p.bw")
        print("  (%s: patched atac.bw via generic rule)" % os.path.basename(src))
    assert "data/atac_5p.bw" in s, "no 5' ATAC path in " + dst
    assert "/data/atac.bw" not in s, "stale atac.bw remains in " + dst
    assert s != orig
    # insert the provenance header after the shebang
    lines = s.split("\n")
    lines.insert(1, HDR.format(src=os.path.basename(src)))
    open(dst, "w").write("\n".join(lines))
    os.chmod(dst, 0o755)
    print("wrote", dst)
