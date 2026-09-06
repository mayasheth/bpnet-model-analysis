#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/p300tx.%j.txt
#SBATCH -e log/p300tx.%j.txt
#SBATCH --job-name=p300tx
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# The control the p300 benchmark transfer result is missing.
#
# The CRISPR benchmark showed a GM12878-trained p300 model applied to K562 sitting at the
# ATAC-only floor (+0.009 [-0.004, +0.022]) while the K562-trained model clears it by +0.055.
# That drop conflates TWO things: what is lost by changing cell type, and the possibility that
# the GM12878 p300 model is simply weaker -- it comes from a different experiment
# (ENCSR000DZG against ENCSR000EGE), a different project, and has never been scored in its own
# cell type.
#
# The benchmark cannot separate them, because CRISPR data exists only for K562. Correlation
# metrics can, by scoring both models in BOTH cell types -- the same local-versus-transferred
# 2x2 the H3K27ac transfer matrix (2.24) used throughout.
#
#   GM local vs K562 local     is the GM12878 p300 model comparably good at home?
#   transferred vs local       the transfer drop, measured per direction on identical data
#   transferred vs atac_LOCAL  does transferring beat the target's own accessibility?
#
# If the GM12878 model is comparably good in GM12878, the benchmark drop is about transfer.
# If it is much worse at home, the benchmark drop is partly model quality and the transfer
# claim has to be weakened.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GMEL=$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "########## evaluated on K562 (target: K562 EP300) ##########"
$PY scripts/2.15.perfold_from_config.py config/p300tx_k562_configs.json \
    p300tx_k562_ "$K5EL" \
    --pair GMp300_transferred K562p300_local \
    --pair GMp300_transferred K562atac_local \
    --pair K562p300_local K562atac_local
echo
echo "########## evaluated on GM12878 (target: GM12878 EP300) ##########"
$PY scripts/2.15.perfold_from_config.py config/p300tx_gm_configs.json \
    p300tx_gm_ "$GMEL" \
    --pair K562p300_transferred GMp300_local \
    --pair K562p300_transferred GMatac_local \
    --pair GMp300_local GMatac_local
