#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/p300depthtest.%j.txt
#SBATCH -e log/p300depthtest.%j.txt
#SBATCH --job-name=p300depthtest
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does training-signal volume explain the p300 transfer asymmetry?
#
# The K562-trained model transferred to GM12878 keeping 83% of its local advantage; the
# GM12878-trained one transferred to K562 keeping none, while both fit their own cell type
# equally well. K562 carries 2.30x the reads in peaks. This scores a K562 model retrained on
# GM12878's budget (30.0M reads, 21,068 peaks) alongside the full-depth original.
#
# The decisive pairing is the transferred one, on GM12878:
#   K562_DEPTHMATCHED_transferred - GMatac_local   does the thinner K562 model still clear
#                                                  the target cell type's own accessibility?
# If it does not, volume explains the asymmetry. If it does, something else about K562 as a
# training cell type does, and a third cell type is required to say what.
#
# The in-cell arm is the control that keeps the transfer arm interpretable: a depth-matched
# model that is simply worse everywhere says nothing about portability specifically.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "########## evaluated on K562 (in-cell control) ##########"
$PY scripts/2.15.perfold_from_config.py config/p300depth_k562_configs.json \
    p300depth_k562_ "$K5EL" \
    --pair K562p300_local_DEPTHMATCHED K562p300_local_FULL \
    --pair K562p300_local_DEPTHMATCHED K562atac_local
echo
echo "########## evaluated on GM12878 (the transfer test) ##########"
$PY scripts/2.15.perfold_from_config.py config/p300depth_gm_configs.json \
    p300depth_gm_ "$GMEL" \
    --pair K562p300_transferred_DEPTHMATCHED GMatac_local \
    --pair K562p300_transferred_FULL GMatac_local \
    --pair K562p300_transferred_DEPTHMATCHED K562p300_transferred_FULL
