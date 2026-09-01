#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/residual_transfer.%j.txt
#SBATCH -e log/residual_transfer.%j.txt
#SBATCH --job-name=k27_res_transfer
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Cross-cell-type transfer of RESIDUAL-objective models vs TOTAL-target models.
#
# THE DESIGN. A residual model predicts `observed - atac_pred`, so transferring it requires
# choosing whose atac_pred. We use the TARGET cell type's own ATAC-only model, reading the
# TARGET cell type's observed ATAC. Nothing about the source cell type's accessibility
# crosses over; only the sequence-derived complement transfers. Concretely, for a
# K562-trained residual model scored on GM12878:
#
#     prediction = K562_residual_model(GM12878 sequence)
#                + gm12878_atac5p_model(GM12878 observed ATAC)
#
# The same atac_pred defines the metric baseline (`true_resid = observed - atac_pred`) and
# the model's offset, so residual_pearson stays internally consistent -- we are evaluating
# against the true residual in the target cell type, which is the metric that matters.
#
# Baseline is the ATAC-only MODEL's prediction, not the raw ATAC track. That was decided on
# 2026-08-25 (.living/decisions.md): with the raw track the residual absorbs everything a
# linear read of the track misses, which flatters any model that merely learns a better
# ATAC transform.
#
# WHY THIS MIGHT FAVOUR THE RESIDUAL DESIGN. A multimodal model carries an internal
# ATAC->H3K27ac mapping that is cell-type-specific (model-free coupling is 0.510 in K562 vs
# 0.409 in GM12878), and that mapping travels with it. The residual design never transfers
# that mapping -- accessibility is re-measured locally and only sequence moves. If sequence
# grammar is more portable than the accessibility-to-acetylation relationship, the residual
# models should retain more. If they do not, the sequence component is the cell-type-specific
# part, which is the more interesting outcome.
#
# In-cell references already exist (residual_grid_* and gm12878_residual_grid_*), giving the
# four evaluations F-003 requires for any cross-cell-type claim.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLDS=$D/reference/hg38_five_folds.json
cd "$P"

echo "### 1/2 K562-trained -> GM12878 ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/transfer_k562_to_gm_configs.json \
    --out-prefix transfer_k562_to_gm_ \
    --elements "$TR/reference/GM12878_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/gm12878_h3k27ac_5p_plus.bw" \
    --accessibility-bw "$TR/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo
echo "### 2/2 GM12878-trained -> K562 ###"
$PY scripts/2.4.evaluate_residual.py \
    --config-json config/transfer_gm_to_k562_configs.json \
    --out-prefix transfer_gm_to_k562_ \
    --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/h3k27ac_5p_plus.bw" \
    --accessibility-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "ALL_DONE"
