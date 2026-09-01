#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/transfer_perfold.%j.txt
#SBATCH -e log/transfer_perfold.%j.txt
#SBATCH --job-name=k27_tx_perfold
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Per-fold CIs and paired tests for the four transfer runs. The pooled tables show
# multimodal-vs-residual gaps of 0.001-0.016, which is below this project's between-fold sd
# (0.041-0.046) -- exactly the regime where only a within-fold paired test can say whether
# the difference is real.
#
# Two variants per direction:
#   transfer_*  offset from the TARGET cell type's own ATAC model. Optimistic: that model is
#               trained on target H3K27ac, which the deployment scenario lacks.
#   deploy_*    offset from the SOURCE ATAC model, transferred. Only target ATAC is used,
#               which is the real application.
#
# Headline metric is overall_pearson against OBSERVED H3K27ac -- the quantity you are trying
# to produce in a cell type where you have none.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

run () {  # cfg out_prefix elements  mmLabel resMMLabel resSeqLabel
  echo "########## $2 ##########"
  $PY scripts/2.15.perfold_from_config.py "config/$1" "$2" "$3" \
      --pair "$5" "$4" --pair "$6" "$4"
  echo
}

run transfer_k562_to_gm_configs.json  transfer_k562_to_gm_  "$GMEL" K562mm_to_GM   K562resMM_to_GM   K562resSeq_to_GM
run deploy_k562_to_gm_configs.json    deploy_k562_to_gm_    "$GMEL" K562mm_to_GM   K562resMM_to_GM   K562resSeq_to_GM
run transfer_gm_to_k562_configs.json  transfer_gm_to_k562_  "$K5EL" GMmm_to_K562   GMresMM_to_K562   GMresSeq_to_K562
run deploy_gm_to_k562_configs.json    deploy_gm_to_k562_    "$K5EL" GMmm_to_K562   GMresMM_to_K562   GMresSeq_to_K562

echo "ALL_DONE"
