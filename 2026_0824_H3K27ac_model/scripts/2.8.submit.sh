#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -o log/residual5p_perfold.%j.txt
#SBATCH -e log/residual5p_perfold.%j.txt
#SBATCH --job-name=k27_res5p_pf
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $D/2026_0824_H3K27ac_model
$D/.pixi/envs/multimodal/bin/python scripts/2.8.residual_perfold_and_artifact.py
echo "ALL_DONE"
