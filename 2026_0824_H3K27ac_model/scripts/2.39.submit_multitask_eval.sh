#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/mteval.%j.txt
#SBATCH -e log/mteval.%j.txt
#SBATCH --job-name=multitask_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Does giving the profile head a LEARNABLE target make the counts head better, and is the
# answer just the stopping rule?
#
# Three arms, paired within fold against the status-quo baseline:
#   dnase_profile_head      profile head on K562 DNase, counts head on H3K27ac (1.28)
#   epoch_matched_baseline  the baseline given the multi-task arm's 100-epoch budget (1.29)
#
# WHY THE CONTROL. Early stopping watches the TOTAL validation loss, and for the multi-task
# arm that includes the reweighted DNase profile term, so the arms stop on different
# objectives: 53-99 epochs against 32-55, every fold longer. Without the control, "a
# learnable auxiliary task improves the trunk" and "this arm trained longer" fit the result
# equally well.
#
# HOW TO READ IT. The two paired deltas are only meaningful together:
#   control ~ multi-task  -> the win was the stopping rule, and the auxiliary task is moot
#   control ~ baseline    -> the auxiliary task is doing real work
# The additional pair `dnase_profile_head - epoch_matched_baseline` states it directly, which
# is the comparison the whole control exists to make.
#
# READ ONLY THE COUNTS COLUMNS. overall_pearson, overall_pearson_topq, residual_pearson and
# incremental_r2 score the counts head against observed H3K27ac. The `profile_*` columns
# score every arm's profile head against H3K27ac, so for the multi-task arm they measure how
# well a DNase-trained head happens to predict H3K27ac shape. A drop there is by
# construction and carries no information.
#
# --no-rc-average keeps these comparable to every other table in the project.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

# Fail before burning a GPU hour: a preempted fold leaves a checkpoint that looks finished.
for m in multimodal5p_accs5p_hw500_clw10 \
         multimodal5p_accs5p_dnaseprof00561_hw500_clw10 \
         multimodal5p_accs5p_ep100_hw500_clw10; do
    for f in 0 1 2 3 4; do
        [[ -s "$P/models/$m/fold$f/training_complete.json" ]] || {
            echo "ERROR: $m/fold$f has no training_complete.json, so it is still running "\
"or was preempted. Scoring it would mix a best-so-far checkpoint into the comparison." >&2
            exit 1
        }
    done
done

$PY scripts/2.15.perfold_from_config.py config/multitask_k562_configs.json \
    multitask_k562_ "$EL" --no-rc-average \
    --pair dnase_profile_head atac_input \
    --pair epoch_matched_baseline atac_input \
    --pair dnase_profile_head epoch_matched_baseline
echo MULTITASK_EVAL_DONE
