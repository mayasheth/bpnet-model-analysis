#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -o log/reciprocal.%j.txt
#SBATCH -e log/reciprocal.%j.txt
#SBATCH --job-name=k27_reciprocal
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Both transfer directions on the 5'-end target, plus the in-cell-type reference each
# needs. Retention can only be read against the model's own cell type, so all four
# evaluations must exist on the same target processing:
#
#   1. GM-trained  -> GM12878   in-cell-type reference for the reciprocal direction
#   2. GM-trained  -> K562      RECIPROCAL transfer
#   3. K562-trained-> GM12878   forward transfer, re-derived on 5' (was fragment)
#   (K562-trained -> K562 already exists as fiveprime_*)
#
# Usage: sbatch 2.6.submit_reciprocal.sh

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
T=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLDS=$D/reference/hg38_five_folds.json
KEL=$D/reference/K562_DNase_candidate_elements.narrowPeak
GEL=$T/reference/GM12878_candidate_elements.narrowPeak
cd "$P"

echo "### 1/3 GM-trained on GM12878 (in-cell-type reference) ###"
$PY scripts/2.2.evaluate_stratified.py --config-json config/gm_incell_configs.json \
    --out-prefix gm_incell_ --elements "$GEL" --genome "$GEN" \
    --signal-bw data/gm12878_h3k27ac_5p_plus.bw --accessibility-bw "$T/data/atac.bw" \
    --fold-json "$FOLDS" --ceiling results/gm12878_fiveprime_replicate_ceiling_by_window.tsv \
    --outdir results --figdir figures

echo "### 2/3 GM-trained on K562 (RECIPROCAL transfer) ###"
$PY scripts/2.2.evaluate_stratified.py --config-json config/gm_to_k562_configs.json \
    --out-prefix gm_to_k562_ --elements "$KEL" --genome "$GEN" \
    --signal-bw data/h3k27ac_5p_plus.bw \
    --accessibility-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw" \
    --fold-json "$FOLDS" --ceiling results/fiveprime_replicate_ceiling_by_window.tsv \
    --outdir results --figdir figures

echo "### 3/3 K562-trained on GM12878, re-derived on 5-prime ###"
$PY scripts/2.2.evaluate_stratified.py --config-json config/k562_to_gm_5p_configs.json \
    --out-prefix k562_to_gm_5p_ --elements "$GEL" --genome "$GEN" \
    --signal-bw data/gm12878_h3k27ac_5p_plus.bw --accessibility-bw "$T/data/atac.bw" \
    --fold-json "$FOLDS" --ceiling results/gm12878_fiveprime_replicate_ceiling_by_window.tsv \
    --outdir results --figdir figures

echo ALL_DONE
