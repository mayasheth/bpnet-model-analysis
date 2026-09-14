#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 16:00:00
#SBATCH --mem=96G
#SBATCH -o log/eval2x2.%j.txt
#SBATCH -e log/eval2x2.%j.txt
#SBATCH --job-name=eval_2x2
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# The 2x2: four accessibility inputs x two training cell types, each scored in both cell
# types. Four runs, each with the ATAC-input arm of that same direction as its baseline, so
# every delta reads as "is this input better than plain ATAC, here, in this direction".
#
# A transferred model is fed the TARGET cell type's version of its own input variant. Feeding
# it the source cell type's track would measure nothing anyone would deploy.
#
# WHY BOTH DIRECTIONS. F-009 found p300 transfer strongly asymmetric: K562 -> GM12878 retained
# most of its advantage and GM12878 -> K562 retained none. Nothing says H3K27ac under a
# DNase-family input behaves the same way in both directions, and one direction cannot tell
# the difference between "this input transfers" and "K562 transfers".
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

run () {  # train_cell target_cell elements
    echo "########## trained on $1, scored on $2 ##########"
    $PY scripts/2.15.perfold_from_config.py \
        "config/h27_2x2_$1_on_$2_configs.json" "h27_2x2_$1_on_$2_" "$3" \
        --no-rc-average \
        --pair real_dnase atac_input \
        --pair dnase_smooth250 atac_input \
        --pair converted_pow atac_input \
        --pair real_dnase converted_pow \
        --pair real_dnase dnase_smooth250
    echo
}

run k562    k562    "$K5EL"
run k562    gm12878 "$GMEL"
run gm12878 gm12878 "$GMEL"
run gm12878 k562    "$K5EL"
echo done
