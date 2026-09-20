#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mceval.%j.txt
#SBATCH -e log/mceval.%j.txt
#SBATCH --job-name=mc_k562_eval
#
# Score the K562 hold-out of the multi-cell p300 panel: does training on four cell types
# beat training on one, when neither has seen the target?
#
# Every arm is scored on K562 p300 over K562's own EP300 peaks with K562 ATAC as the
# accessibility input, so the arms differ only in what they were trained on. The baseline
# is K562's own accessibility-only p300 model, the assay-matched floor F-016 requires.
#
# THE DECISIVE PAIR IS THE FIRST ONE. multicell_holdK562 against
# GM_single_corrected_tgtnorm compares four source cell types against one under the SAME
# normalisation policy. The multi-cell model has no single training statistic so it must be
# normalised on the target; comparing it against a single-source model using its own
# training statistics would mix the training set with the normalisation. The second pair,
# against GM_single_corrected, is that unmatched comparison, kept so the size of the
# policy effect is visible rather than assumed.
#
# out_prefix takes no directory and needs its own trailing underscore: 2.15 writes
# f"{P}/results/{out_prefix}per_fold.tsv".
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
cd "$P"

"$D/.pixi/envs/multimodal/bin/python" scripts/2.15.perfold_from_config.py \
    config/p300_multicell_k562_configs.json \
    p300_multicell_k562_ \
    "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --pair multicell_holdK562 GM_single_corrected_tgtnorm \
    --pair multicell_holdK562 GM_single_corrected \
    --pair multicell_holdK562 GM_single_corrupted \
    --pair multicell_holdK562 K562_atac_only_floor \
    --pair GM_single_corrected K562_atac_only_floor \
    --pair K562_incell multicell_holdK562
