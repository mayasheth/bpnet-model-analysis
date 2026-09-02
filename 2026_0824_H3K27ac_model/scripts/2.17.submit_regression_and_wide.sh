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
# accessibility inputs. Two things happen here, in this order.
#
# STEP 1, REGRESSION. Re-score an existing all-8-layer comparison and compare against the
# stored table via 2.18: region counts per fold must match EXACTLY (that is what a geometry
# or region-set change would break) and every metric must match within 1e-3 (cuDNN is not
# bit-reproducible across GPU models, so re-scoring on a different node moves the 4th
# decimal place with the code untouched). The wide numbers are not produced unless this
# passes.
#
# STEP 2, the comparison. Baseline is the narrow ATAC-only model in both arms, so it
# contributes identically and the paired tests isolate the variable under test.
# Note the sequence baseline is sequence5p_hw500_clw10: `sequence` mode never reads the
# accessibility input, so there is no _accs5p sequence model to pair against.
#
# Usage: sbatch 2.17.submit_regression_and_wide.sh [SET]
#   SET = wide_seqonly   sequence arm only (runnable as soon as those 10 folds finish)
#         wide           sequence + multimodal (needs all 20 wide folds)
#         fragchan       fragment-size channels, K562 only
#   Default: wide

set -euo pipefail
export PYTHONUNBUFFERED=1
SET=${1:-wide}
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
if $PY scripts/2.18.compare_perfold_tables.py \
        results/REGRESSION_transfer_k562_to_gm_per_fold.tsv \
        results/transfer_k562_to_gm_per_fold.tsv --tol 1e-3; then
    rm -f results/REGRESSION_transfer_k562_to_gm_{per_fold,fold_summary}.tsv
else
    echo "2.15 changed behaviour on all-8-layer configs; not producing wide numbers" >&2
    exit 1
fi

echo
echo "########## STEP 2: $SET ##########"
case "$SET" in
  wide_seqonly) CELLS="k562 gm12878"; PAIRS=(--pair sequence_WIDE sequence_narrow) ;;
  wide)         CELLS="k562 gm12878"; PAIRS=(--pair sequence_WIDE sequence_narrow
                                             --pair multimodal_WIDE multimodal_narrow) ;;
  fragchan)     CELLS="k562";         PAIRS=(--pair multimodal_FRAGCHAN multimodal_flat) ;;
  *) echo "unknown SET '$SET'" >&2; exit 1 ;;
esac

for CELL in $CELLS; do
    [[ "$CELL" == k562 ]] && EL="$K5EL" || EL="$GMEL"
    echo "===== $CELL ====="
    $PY scripts/2.15.perfold_from_config.py "config/${SET}_${CELL}_configs.json" \
        "${SET}_${CELL}_" "$EL" "${PAIRS[@]}"
    echo
done
