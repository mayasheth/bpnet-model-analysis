#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o log/wide_eval.%j.txt
#SBATCH -e log/wide_eval.%j.txt
#SBATCH --job-name=k27_wide_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# 2.15 was generalised to handle models with different receptive fields and different
# accessibility inputs. Two things must happen in this order.
#
# STEP 1, REGRESSION. Re-score an existing all-8-layer comparison and require the per-fold
# table to be byte-identical to the stored one. For 8-layer models the new code path
# resolves in_window to 2114 and the centre-crop is a no-op, so any difference at all means
# the generalisation changed behaviour it should not have touched.
#
# STEP 2, the wide comparison itself. Baseline is the narrow ATAC-only model in both arms,
# so it contributes identically and the paired tests isolate the receptive field.
# Note the sequence baseline is sequence5p_hw500_clw10: `sequence` mode never reads the
# accessibility input, so there is no _accs5p sequence model to pair against.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

echo "########## STEP 1: backward-compatibility regression ##########"
$PY scripts/2.15.perfold_from_config.py config/transfer_k562_to_gm_configs.json \
    REGRESSION_transfer_k562_to_gm_ "$GMEL" \
    --pair K562resMM_to_GM K562mm_to_GM --pair K562resSeq_to_GM K562mm_to_GM
if diff -q results/REGRESSION_transfer_k562_to_gm_per_fold.tsv \
           results/transfer_k562_to_gm_per_fold.tsv; then
    echo "REGRESSION PASS: per-fold table unchanged"
    rm -f results/REGRESSION_transfer_k562_to_gm_{per_fold,fold_summary}.tsv
else
    echo "REGRESSION FAIL: 2.15 changed behaviour on all-8-layer configs" >&2
    diff results/REGRESSION_transfer_k562_to_gm_per_fold.tsv \
         results/transfer_k562_to_gm_per_fold.tsv | head -40 >&2
    exit 1
fi

echo
echo "########## STEP 2: wide vs narrow receptive field ##########"
for CELL in k562 gm12878; do
    [[ "$CELL" == k562 ]] && EL="$K5EL" || EL="$GMEL"
    echo "===== $CELL ====="
    $PY scripts/2.15.perfold_from_config.py "config/wide_${CELL}_configs.json" \
        "wide_${CELL}_" "$EL" \
        --pair sequence_WIDE sequence_narrow --pair multimodal_WIDE multimodal_narrow
    echo
done
