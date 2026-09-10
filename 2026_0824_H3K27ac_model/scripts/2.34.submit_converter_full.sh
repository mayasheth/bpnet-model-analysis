#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -o log/conv_full.%j.txt
#SBATCH -e log/conv_full.%j.txt
#SBATCH --job-name=conv_full
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# The five-fold converter result: everything the fold-0 pilot showed, with folds to put a
# CI on it and a paired test behind it.
#
#   1  counts + smoothed profile, paired against the ATAC-only control (2.15, all folds)
#   2  raw shape vs the DNase ceiling, in-cell and transferred, both arms (2.31)
#   3  the deployment comparator: raw observed ATAC vs observed DNase, per fold (0.35)
#
# Step 3 is the one that makes step 2 interpretable. Fraction-of-DNase-ceiling says how
# close the converter gets to DNase; only the comparison against the untouched ATAC track
# says whether swapping it in changes what the downstream model reads.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL_ATAC=$D/reference/K562_ATAC_candidate_elements.narrowPeak
cd "$P"

echo "########## STEP 1: counts, paired over 5 folds ##########"
$PY scripts/2.15.perfold_from_config.py config/converter_k562_5fold_configs.json \
    converter5f_k562_ "$K5EL_ATAC" --no-rc-average \
    --pair seq_atac_clw10 atac_only_clw10

echo
echo "########## STEP 2: raw shape vs ceiling, 5 folds, in-cell and transferred ##########"
shape () {  # label, model_dir, mode, cell
    echo "===== $1 on $4 ====="
    $PY scripts/2.31.converter_profile_eval.py \
        --model-dir "$P/models/$2" --fold 0 --fold 1 --fold 2 --fold 3 --fold 4 \
        --label "$1" --cell "$4" --mode "$3"
    echo
}
shape converter5f  multimodal5p_conv_atac2dnase_atacel_hw500_clw10 multimodal k562
shape ataconly5f   atac5p_conv_atac2dnase_atacel_hw500_clw10       atac       k562
shape converter5f  multimodal5p_conv_atac2dnase_atacel_hw500_clw10 multimodal gm12878
shape ataconly5f   atac5p_conv_atac2dnase_atacel_hw500_clw10       atac       gm12878
echo done
