#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/gcmatch.%j.txt
#SBATCH -e log/gcmatch.%j.txt
#SBATCH --job-name=k27_gcmatch
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does GC-matching the training negatives help? Paired within fold against the matched-settings
# baseline, which differs only in the negatives file (1.11 vs 1.20).
#
# Scored in-cell K562 AND transferred to GM12878. The transfer arm is the one the hypothesis
# predicts: if the unmatched pool taught the sequence branch to be a GC detector, the cost
# should show up most where accessibility does not generalise.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "########## in-cell K562 ##########"
$PY scripts/2.15.perfold_from_config.py config/gcmatch_configs.json \
    gcmatch_k562_ "$K5EL" --pair negatives_GCMATCHED negatives_UNMATCHED
echo
echo "########## K562 -> GM12878 (deployment) ##########"
$PY scripts/2.15.perfold_from_config.py config/gcmatch_transfer_configs.json \
    gcmatch_k562_to_gm_ "$GMEL" \
    --pair negatives_GCMATCHED_transferred negatives_UNMATCHED_transferred \
    --pair negatives_GCMATCHED_transferred atac_LOCAL
