#!/bin/bash
#SBATCH -p gpu
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -o log/accs5p_compare.%j.txt
#SBATCH -e log/accs5p_compare.%j.txt
#SBATCH --job-name=k27_accs5p_cmp
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does the ChromBPNet-convention 5' accessibility input change model performance?
#
# Each config pairs OLD (genomecov -bg over the full tagAlign interval, so coverage scales
# with read length) against 5PRIME (genomecov -bg -5, single-base insertion counts) for the
# same mode, same target, same elements, same folds. The two arms differ ONLY in the
# accessibility bigwig -- the training scripts were byte-identical copies apart from the
# ATAC path and output dir -- so any difference is attributable to the input definition.
#
# 2.2 reads `accessibility_bw` per config entry, so both arms are scored with their own
# input and their own saved normalisation statistics in a single paired run.
#
# Expected: little or no difference WITHIN a cell type, because read length is constant
# there and z-normalisation absorbs a pure scale change. The value of the switch is
# ACROSS cell types (TeloHAEC at 36 bp vs K562/GM12878 at 94-95 bp), which this job does
# not test. A large within-cell-type difference would be the surprise and would mean the
# smear was doing something beyond scaling.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
P3=$D/2026_0529_multimodal_p300_model
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLDS=$D/reference/hg38_five_folds.json
cd "$P"

echo "########## 1/3 K562 H3K27ac ##########"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/accs5p_k562_configs.json --out-prefix accs5p_k562_ \
    --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/h3k27ac_5p_plus.bw" \
    --accessibility-bw "$P3/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "########## 2/3 GM12878 H3K27ac ##########"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/accs5p_gm12878_configs.json --out-prefix accs5p_gm12878_ \
    --elements "$TR/reference/GM12878_candidate_elements.narrowPeak" --genome "$GEN" \
    --signal-bw "$P/data/gm12878_h3k27ac_5p_plus.bw" \
    --accessibility-bw "$TR/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "########## 3/3 p300 ##########"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/accs5p_p300_configs.json --out-prefix accs5p_p300_ \
    --elements "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" --genome "$GEN" \
    --signal-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --accessibility-bw "$P3/data/atac.bw" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "ALL_DONE"
