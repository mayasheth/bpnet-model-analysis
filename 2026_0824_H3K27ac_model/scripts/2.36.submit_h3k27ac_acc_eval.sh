#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/h27acceval.%j.txt
#SBATCH -e log/h27acceval.%j.txt
#SBATCH --job-name=h27_acc_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does a better accessibility input improve H3K27ac prediction, and is it shape or magnitude?
#
# Four arms on one element set, baseline = ATAC input:
#   real_dnase       F-010, the target to match
#   converted_dnase  the converter's first test where its profile can matter
#   dnase_smooth250  real DNase, structure erased, magnitude kept: THE CONTROL
#
# READ THE CONTROL FIRST. If smooth250 keeps real DNase's advantage, shape never mattered and
# the converter cannot work as an input no matter how good its profile gets. If smooth250
# loses it, shape is the thing and a converted arm below real DNase is a converter-quality
# problem rather than a dead idea.
#
# Elements are the DNase-derived set, matching 1.11 and 1.22 so every arm is comparable to
# its own anchors.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

$PY scripts/2.15.perfold_from_config.py config/h3k27ac_accvariants_configs.json \
    h27acc_k562_ "$EL" --no-rc-average \
    --pair real_dnase atac_input \
    --pair converted_dnase atac_input \
    --pair dnase_smooth250 atac_input \
    --pair real_dnase converted_dnase \
    --pair real_dnase dnase_smooth250
echo done
