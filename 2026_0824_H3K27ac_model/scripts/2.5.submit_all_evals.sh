#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/evals_perfold.%j.txt
#SBATCH -e log/evals_perfold.%j.txt
#SBATCH --job-name=k27_evals_pf
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Re-run all three stratified evaluations so each emits PER-FOLD statistics
# (*_stratified_per_fold.tsv) and a mean +/- 95% CI summary
# (*_stratified_fold_summary.tsv) alongside the pooled table.
#
# Why per-fold: each chromosome-holdout fold is an independent replicate of the
# train-and-evaluate procedure, so the spread across folds is the only uncertainty
# estimate available. Measured run-to-run variance is ~0.018 on a single fold, and
# several comparisons in this analysis rest on differences of that size — a pooled
# correlation reports them with no error bar at all.
#
# Usage: sbatch 2.5.submit_all_evals.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
ATAC=$D/2026_0529_multimodal_p300_model/data/atac.bw
FOLDS=$D/reference/hg38_five_folds.json
EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak

cd "$P"

echo "### 1/3 five-prime three-mode comparison ###"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/fiveprime_eval_configs.json --out-prefix fiveprime_ \
    --elements "$EL" --genome "$GEN" --signal-bw data/h3k27ac_5p_plus.bw \
    --accessibility-bw "$ATAC" --fold-json "$FOLDS" \
    --ceiling results/fiveprime_replicate_ceiling_by_window.tsv \
    --outdir results --figdir figures

echo "### 2/3 p300 vs H3K27ac ###"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/compare_configs.json --out-prefix p300_vs_h3k27ac_ \
    --elements "$EL" --genome "$GEN" \
    --signal-bw /oak/stanford/groups/engreitz/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw \
    --accessibility-bw "$ATAC" --fold-json "$FOLDS" \
    --outdir results --figdir figures

echo "### 3/3 GM12878 transfer ###"
$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/gm_transfer_configs.json --out-prefix gm12878_transfer_ \
    --elements "$GMEL" --genome "$GEN" \
    --signal-bw data/gm12878_h3k27ac_frag250.bw \
    --accessibility-bw "$D/2026_0606_GM12878_transferability/data/atac.bw" \
    --fold-json "$FOLDS" --outdir results --figdir figures

echo "ALL_DONE"
