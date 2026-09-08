#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/dnaseswap.%j.txt
#SBATCH -e log/dnaseswap.%j.txt
#SBATCH --job-name=dnaseswap
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does DNase beat ATAC as the model's accessibility input, in-cell type and on transfer?
#
# Each entry carries its OWN accessibility bigwig, so a DNase-trained model is always fed
# DNase and an ATAC-trained model always ATAC, including when transferred. The comparison is
# therefore assay-versus-assay end to end rather than a model being fed the wrong input.
#
# The transfer arm is the one that decides whether the ATAC->DNase converter is worth
# building: a converter only pays if DNase is the better representation to be in, and it only
# pays for DEPLOYMENT if the DNase advantage survives the move between cell types.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "########## in-cell K562 ##########"
$PY scripts/2.15.perfold_from_config.py config/dnaseswap_k562_configs.json \
    dnaseswap_k562_ "$K5EL" --pair DNase_input ATAC_input
echo
echo "########## GM12878: in-cell and K562-transferred, both assays ##########"
$PY scripts/2.15.perfold_from_config.py config/dnaseswap_gm_configs.json \
    dnaseswap_gm_ "$GMEL" \
    --pair DNase_input_local ATAC_input_local \
    --pair DNase_input_K562transferred ATAC_input_K562transferred \
    --pair DNase_input_K562transferred atac_LOCAL_floor \
    --pair ATAC_input_K562transferred atac_LOCAL_floor
