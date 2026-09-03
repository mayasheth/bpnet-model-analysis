#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/widetx.%j.txt
#SBATCH -e log/widetx.%j.txt
#SBATCH --job-name=k27_widetx
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does the wider receptive field still pay off across cell types?
#
# The in-cell-type gain is +0.027 (K562) and +0.014 (GM12878) on the top quintile, and the
# deployment conclusion -- transfer the multimodal model -- was measured on the narrow one.
# Each table carries four arms so the transferred and local versions of both geometries sit
# side by side, making the transfer drop readable without cross-referencing another file.
#
# RC averaging is ON (the 2.15 default since 2026-09-03) for every arm here, so these
# numbers are not comparable to the single-pass transfer tables from 2026-09-01.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

echo "########## K562 -> GM12878 ##########"
$PY scripts/2.15.perfold_from_config.py config/widetx_k562_to_gm12878_configs.json \
    widetx_k562_to_gm12878_ "$GMEL" \
    --pair k562_mm_WIDE_to_gm12878 k562_mm_narrow_to_gm12878 \
    --pair gm12878_mm_WIDE_local gm12878_mm_narrow_local

echo
echo "########## GM12878 -> K562 ##########"
$PY scripts/2.15.perfold_from_config.py config/widetx_gm12878_to_k562_configs.json \
    widetx_gm12878_to_k562_ "$K5EL" \
    --pair gm12878_mm_WIDE_to_k562 gm12878_mm_narrow_to_k562 \
    --pair k562_mm_WIDE_local k562_mm_narrow_local
