#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 6:00:00
#SBATCH --mem=96G
#SBATCH -o log/conv_eval.%j.txt
#SBATCH -e log/conv_eval.%j.txt
#SBATCH --job-name=conv_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Score the ATAC->DNase converter pilot. TWO SCORING PASSES, BECAUSE THEY ANSWER DIFFERENT
# QUESTIONS AND ONLY ONE OF THEM IS THE POINT.
#
#   STEP 1, 2.15: counts and smoothed-profile metrics, with the ATAC-ONLY arm as the
#   baseline. This is the control that decides whether there is a converter here at all.
#   If `incremental_r2` for the sequence+ATAC arms is ~0, the model is a learned rescaling
#   of its own input and nothing downstream will see anything ATAC normalisation could not
#   already provide.
#
#   STEP 2, 2.31: RAW within-element shape correlation across bin sizes, against the
#   inter-replicate ceiling computed on the same held-out windows. This is the deliverable.
#   2.15's profile_pearson is Gaussian-smoothed (kernel_sigma=7, width=81) and cannot be
#   read against 0.25's raw 1 bp ceiling of 0.848 -- doing so would report the converter as
#   closer to the ceiling than it is.
#
# Memory is 96G to match 2.19: profiles for every held-out element are held alive.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL_ATAC=$D/reference/K562_ATAC_candidate_elements.narrowPeak
cd "$P"

echo "########## STEP 1: counts + paired comparison vs the ATAC-only control ##########"
$PY scripts/2.15.perfold_from_config.py config/converter_k562_configs.json \
    converter_k562_ "$K5EL_ATAC" --no-rc-average --folds 0 \
    --pair seq_atac_clw10 atac_only_clw10 \
    --pair seq_atac_clw1 atac_only_clw10 \
    --pair seq_atac_clw1 seq_atac_clw10

echo
echo "########## STEP 2: raw profile shape vs the inter-replicate ceiling ##########"
run_shape () {  # label, model_dir, mode
    echo "===== $1 ====="
    $PY scripts/2.31.converter_profile_eval.py \
        --model-dir "$P/models/$2" --fold 0 --label "$1" --cell k562 --mode "$3"
    echo
}
run_shape seq_atac_clw10 multimodal5p_conv_atac2dnase_atacel_hw500_clw10 multimodal
run_shape seq_atac_clw1  multimodal5p_conv_atac2dnase_atacel_hw500_clw1  multimodal
run_shape atac_only_clw10 atac5p_conv_atac2dnase_atacel_hw500_clw10      atac

echo "done"
