#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -o log/residual5p_eval.%j.txt
#SBATCH -e log/residual5p_eval.%j.txt
#SBATCH --job-name=k27_res5p_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Evaluate the 5'-target residual models (models/residual5pFIXED_hw500_clw10), which were
# trained but never scored. Open question behind F-002: an independently-trained sequence
# model captures almost none of the residual left by ATAC (r=0.100). Does a model trained
# EXPLICITLY on (observed - atac_pred) recover a sequence complement that the jointly
# trained model missed, or is the residual genuinely not sequence-predictable?
#
# Step 1 confirms empirically that the offset-trained model emits a residual rather than
# logcounts, since 2.4's arithmetic depends on that and a wrong guess fails silently.

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

echo "### 1/2 confirm offset-model output semantics ###"
$PY scripts/2.7.diagnose_residual_offset.py

echo
echo "### 2/2 residual evaluation ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/residual5p_eval_configs.json \
    --out-prefix residual5p_ \
    --elements "$EL" --genome "$GEN" \
    --signal-bw "$P/data/h3k27ac_5p_plus.bw" \
    --accessibility-bw "$ATAC" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "ALL_DONE"
