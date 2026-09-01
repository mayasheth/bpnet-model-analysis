#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/gm12878_residual_grid.%j.txt
#SBATCH -e log/gm12878_residual_grid.%j.txt
#SBATCH --job-name=gm_resgrid
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# GM12878 residual grid, the same 2x2-plus-control run as K562.
#
# In K562: residual training gained +0.310 for sequence-only, COST -0.037 for multimodal
# (all five folds), and the ATAC control scored -0.003. The question is whether that is a
# property of the objective or of K562's unusually tight ATAC-H3K27ac coupling. GM12878
# has lower coupling (0.409 vs 0.510) so more residual is available, but a higher ceiling
# (0.832 vs 0.760) so less of it is noise -- the two pull opposite ways.
#
# Uses the same evaluator with only cell-type paths substituted, so any difference is data,
# not method.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLDS=$D/reference/hg38_five_folds.json
cd "$P"

echo "### 1/2 pooled residual evaluation, GM12878 ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/gm12878_residual_grid_configs.json \
    --out-prefix gm12878_residual_grid_ \
    --elements "$TR/reference/GM12878_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/gm12878_h3k27ac_5p_plus.bw" \
    --accessibility-bw "$TR/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo
echo "### 2/2 per-fold CIs + artifact controls, GM12878 ###"
$PY scripts/2.11.gm12878_residual_perfold.py gm12878_residual_grid_

echo "ALL_DONE"
