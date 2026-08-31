#!/bin/bash
#SBATCH -p gpu
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/residual_grid_eval.%j.txt
#SBATCH -e log/residual_grid_eval.%j.txt
#SBATCH --job-name=k27_resgrid
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Evaluate the complete residual grid: 3 standard models + 3 residual-objective models.
#
# The question this answers, which one residual model could not: is the gain from
# residual TRAINING, or from the residual being sequence-predictable at all?
#   residual_multimodal  if it beats jointly-trained multimodal (0.551 residual r), the
#                        objective matters even when ATAC is already an input.
#   residual_atac        negative control, expected ~0: it predicts an ATAC model's own
#                        errors from the same ATAC input. Clearly non-zero would mean the
#                        residual metric measures something other than we think, which
#                        would put the sequence result back in question.
#
# On the CURRENT atac.bw, deliberately -- these must be comparable to the existing 0.459 /
# 0.551 numbers. The ChromBPNet-style 5' ATAC rebuild is separate and invalidates the whole
# ATAC-input set at once.
#
# gpu partition only, not owners: 2 of the 10 training folds were preempted there.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
ATAC=$D/2026_0529_multimodal_p300_model/data/atac.bw
FOLDS=$D/reference/hg38_five_folds.json
EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

echo "### 1/2 pooled residual evaluation, 5 models ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/residual_grid_eval_configs.json \
    --out-prefix residual_grid_ \
    --elements "$EL" --genome "$GEN" \
    --signal-bw "$P/data/h3k27ac_5p_plus.bw" \
    --accessibility-bw "$ATAC" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo
echo "### 2/2 per-fold CIs + artifact controls ###"
$PY scripts/2.8.residual_perfold_and_artifact.py residual_grid_

echo "ALL_DONE"
