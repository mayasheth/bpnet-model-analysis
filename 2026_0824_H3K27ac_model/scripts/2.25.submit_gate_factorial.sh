#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/gatefac.%j.txt
#SBATCH -e log/gatefac.%j.txt
#SBATCH --job-name=k27_gatefac
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# The 2x2 factorial for the two fixes aimed at accessibility over-reliance:
#   GATE       sequence-gated accessibility (X_acc *= sigmoid(conv(X_seq)))
#   ASYM       overprediction_weight 3.0 in the count loss
#   GATE_ASYM  both
# Baseline is the matched ungated symmetric multimodal model: same accessibility bigwig
# (atac_5p.bw), same n_layers 8, same clw 10. The only differences are the two flags.
#
# Scored twice, because the two questions have different answers:
#   in-cell K562   does the fix cost anything where accessibility is trustworthy?
#   K562->GM12878  does the fix help where over-reliance actually bites? This is the
#                  deployment question -- the gate exists to stop the model leaning on
#                  accessibility that does not generalise.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "########## in-cell K562 ##########"
$PY scripts/2.15.perfold_from_config.py config/gate_factorial_configs.json \
    gatefac_k562_ "$K5EL" \
    --pair GATE multimodal_ungated \
    --pair ASYM multimodal_ungated \
    --pair GATE_ASYM multimodal_ungated \
    --pair GATE_ASYM GATE
echo
echo "########## K562 -> GM12878 (deployment) ##########"
$PY scripts/2.15.perfold_from_config.py config/gate_transfer_configs.json \
    gatefac_k562_to_gm_ "$GMEL" \
    --pair GATE_transferred ungated_transferred \
    --pair ASYM_transferred ungated_transferred \
    --pair GATE_ASYM_transferred ungated_transferred \
    --pair GATE_transferred atac_LOCAL \
    --pair ungated_transferred atac_LOCAL
