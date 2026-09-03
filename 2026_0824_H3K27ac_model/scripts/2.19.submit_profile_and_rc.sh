#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/prof_rc.%j.txt
#SBATCH -e log/prof_rc.%j.txt
#SBATCH --job-name=k27_prof_rc
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Two questions, one job.
#
# 1. IS THE PROFILE HEAD LEARNING ANYTHING? Every metric this project has ever reported
#    comes from the counts head; 2.15 called the model as `_, lc = m(X)` and threw the
#    profile away. profile_pearson and profile_jsd are now computed with the same bpnetlite
#    call the training loop uses, so they are comparable to the Validation Profile Pearson
#    column in the logs. Run over the full K562 residual grid so all three input modes and
#    both objectives are covered at once.
#
# 2. DOES TEST-TIME REVERSE-COMPLEMENT AVERAGING HELP? Training-time RC augmentation is
#    already on; averaging predictions over a sequence and its reverse complement at
#    inference is free and untested here. Because the change is in the evaluator and not the
#    models, the two arms are separate tables, so 2.20 pairs them on (fold, config).
#
# The regression gate runs first, as in 2.17. 2.18 now permits ADDED columns (an older
# reference table predates the profile metrics) while still failing on missing ones.
#
# Memory is 96G rather than 64G: keeping the observed profiles and a predicted profile array
# alive costs ~2 x (n, 2, 1000) float32 on top of the windows.

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
    REGRESSION_transfer_k562_to_gm_ "$GMEL"
if $PY scripts/2.18.compare_perfold_tables.py \
        results/REGRESSION_transfer_k562_to_gm_per_fold.tsv \
        results/transfer_k562_to_gm_per_fold.tsv --tol 1e-3; then
    rm -f results/REGRESSION_transfer_k562_to_gm_{per_fold,fold_summary}.tsv
else
    echo "2.15 changed behaviour on the counts metrics; stopping" >&2
    exit 1
fi

echo
echo "########## STEP 2: profile metrics, K562 residual grid ##########"
$PY scripts/2.15.perfold_from_config.py config/prof_residual_grid_configs.json \
    prof_residual_grid_ "$K5EL" \
    --pair multimodal5p sequence5p --pair residual_multimodal multimodal5p

echo
echo "########## STEP 3: same grid with test-time RC averaging ##########"
$PY scripts/2.15.perfold_from_config.py config/prof_residual_grid_configs.json \
    prof_rc_residual_grid_ "$K5EL" --rc-average

echo
echo "########## STEP 4: what RC averaging is worth (RC - plain) ##########"
$PY scripts/2.20.paired_delta_between_runs.py \
    results/prof_rc_residual_grid_per_fold.tsv \
    results/prof_residual_grid_per_fold.tsv --label-a RC --label-b plain

echo
echo "########## STEP 5: profile metrics for the wide sequence models ##########"
for CELL in k562 gm12878; do
    [[ "$CELL" == k562 ]] && EL="$K5EL" || EL="$GMEL"
    echo "===== $CELL ====="
    $PY scripts/2.15.perfold_from_config.py "config/wide_seqonly_${CELL}_configs.json" \
        "prof_wide_seqonly_${CELL}_" "$EL" --pair sequence_WIDE sequence_narrow
done
