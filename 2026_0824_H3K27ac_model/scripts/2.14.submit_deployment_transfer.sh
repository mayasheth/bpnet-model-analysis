#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/deploy_transfer.%j.txt
#SBATCH -e log/deploy_transfer.%j.txt
#SBATCH --job-name=k27_deploy
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# THE DEPLOYMENT QUESTION: you have ATAC in a new cell type but NO H3K27ac. What is the best
# way to predict H3K27ac there?
#
# Everything transfers from the source cell type; only ATAC is available locally:
#
#     ATAC-only        source_atac_model(target ATAC)
#     multimodal       source_multimodal_model(target sequence + target ATAC)
#     sequence-only    source_sequence_model(target sequence)
#     residual         source_atac_model(target ATAC) + source_residual_model(target seq)
#
# This differs from 2.13, which used the TARGET cell type's own ATAC model as the offset.
# That model is trained on target H3K27ac, which does not exist in this scenario -- so 2.13
# is an optimistic upper bound, not a deployable configuration. Running both is deliberate:
# the gap between them is what you lose by being unable to fit the accessibility baseline
# locally.
#
# METRIC. `overall_pearson` is the one that matters here -- correlation with OBSERVED
# H3K27ac in the target cell type, directly comparable across all four approaches, since
# that is the quantity you are trying to produce. `residual_pearson` is reported alongside
# as the mechanistic readout of whether the sequence component itself survived transfer, but
# it is measured against a baseline that is now itself a transferred model, so read it as
# "beyond what the transferred ATAC model predicts", not as an absolute.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLDS=$D/reference/hg38_five_folds.json
cd "$P"

echo "### 1/2 DEPLOY: K562-trained -> GM12878 (only GM12878 ATAC available) ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/deploy_k562_to_gm_configs.json \
    --out-prefix deploy_k562_to_gm_ \
    --elements "$TR/reference/GM12878_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/gm12878_h3k27ac_5p_plus.bw" \
    --accessibility-bw "$TR/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo
echo "### 2/2 DEPLOY: GM12878-trained -> K562 (only K562 ATAC available) ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/deploy_gm_to_k562_configs.json \
    --out-prefix deploy_gm_to_k562_ \
    --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/h3k27ac_5p_plus.bw" \
    --accessibility-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "ALL_DONE"
