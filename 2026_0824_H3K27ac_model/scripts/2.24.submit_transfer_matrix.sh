#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/txmatrix.%j.txt
#SBATCH -e log/txmatrix.%j.txt
#SBATCH --job-name=k27_txmatrix
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Which architecture transfers best? The in-cell-type ranking is not the deployment ranking.
# Each pairing below asks a distinct question:
#   *_transferred vs atac_LOCAL   is the transferred model better than the target's own ATAC?
#   WIDE vs narrow (transferred)  does the receptive-field gain travel?
#   FRAG vs flat (transferred)    do fragment channels travel, or carry library batch signal?
#   local vs transferred          the transfer drop, per architecture
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

echo "########## K562 -> GM12878 ##########"
$PY scripts/2.15.perfold_from_config.py config/txmatrix_k562_to_gm12878_configs.json \
    txmatrix_k562_to_gm12878_ "$GMEL" \
    --pair WIDE_flat_transferred narrow_flat_transferred \
    --pair narrow_FRAG_transferred narrow_flat_transferred \
    --pair WIDE_FRAG_transferred narrow_flat_transferred \
    --pair narrow_flat_local narrow_flat_transferred \
    --pair WIDE_flat_local WIDE_flat_transferred \
    --pair narrow_FRAG_local narrow_FRAG_transferred
echo
echo "########## GM12878 -> K562 ##########"
$PY scripts/2.15.perfold_from_config.py config/txmatrix_gm12878_to_k562_configs.json \
    txmatrix_gm12878_to_k562_ "$K5EL" \
    --pair WIDE_flat_transferred narrow_flat_transferred \
    --pair narrow_FRAG_transferred narrow_flat_transferred \
    --pair narrow_flat_local narrow_flat_transferred \
    --pair WIDE_flat_local WIDE_flat_transferred \
    --pair narrow_FRAG_local narrow_FRAG_transferred
