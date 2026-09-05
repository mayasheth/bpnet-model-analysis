#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/finalcmp.%j.txt
#SBATCH -e log/finalcmp.%j.txt
#SBATCH --job-name=k27_finalcmp
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# The four comparisons whose models finished on 2026-09-05.
#
# --no-rc-average throughout, deliberately: RC averaging became the 2.15 default on
# 2026-09-03, but every existing figure is single-pass, and 2.22 converts the whole set
# together. Mixing single-pass and RC-averaged tables inside one report is undetectable by a
# reader, so these stay single-pass until that conversion happens.
#
# The ATAC-element comparison is scored on the ATAC-DERIVED element set, which is the set ABC
# uses. Models trained on different element sets cannot be compared on their own sets, so
# both arms are scored on the one that matters downstream.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
K5EL_ATAC=$D/reference/K562_ATAC_candidate_elements.narrowPeak
cd "$P"

echo "########## 1. receptive field, all three input modes ##########"
for CELL in k562 gm12878; do
    [[ "$CELL" == k562 ]] && EL="$K5EL" || EL="$GMEL"
    echo "===== $CELL ====="
    $PY scripts/2.15.perfold_from_config.py "config/wide_${CELL}_configs.json" \
        "wide_${CELL}_" "$EL" --no-rc-average \
        --pair atac_WIDE atac_narrow \
        --pair sequence_WIDE sequence_narrow \
        --pair multimodal_WIDE multimodal_narrow
    echo
done

echo "########## 2. receptive field x fragment channels (2x2) ##########"
$PY scripts/2.15.perfold_from_config.py config/fragwide_k562_configs.json \
    fragwide_k562_ "$K5EL" --no-rc-average \
    --pair mm_narrow_FRAG mm_narrow_flat \
    --pair mm_WIDE_flat mm_narrow_flat \
    --pair mm_WIDE_FRAG mm_narrow_flat \
    --pair mm_WIDE_FRAG mm_WIDE_flat \
    --pair mm_WIDE_FRAG mm_narrow_FRAG
echo

echo "########## 3. element derivation, scored on ATAC-derived elements ##########"
$PY scripts/2.15.perfold_from_config.py config/atacel_k562_configs.json \
    atacel_k562_ "$K5EL_ATAC" --no-rc-average \
    --pair mm_trained_ATAC_elements mm_trained_DNase_elements
echo

echo "########## 4. fragment channels in GM12878 ##########"
$PY scripts/2.15.perfold_from_config.py config/fragchan_gm12878_configs.json \
    fragchan_gm12878_ "$GMEL" --no-rc-average \
    --pair multimodal_FRAGCHAN multimodal_flat
