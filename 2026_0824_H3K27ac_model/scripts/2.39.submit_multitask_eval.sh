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
# Does giving the profile head a LEARNABLE target make the counts head better?
#
# Two arms, paired within fold: multimodal5p_accs5p (profile head on H3K27ac, the status
# quo) against multimodal5p_accs5p_dnaseprof00561 (profile head on K562 DNase). Everything
# else is identical by construction, 1.28 being 1.11 plus two flags.
#
# READ ONLY THE COUNTS COLUMNS. overall_pearson, overall_pearson_topq, residual_pearson and
# incremental_r2 all score the counts head against observed H3K27ac, which is the quantity
# of interest and the thing the hypothesis is about. The `profile_*` columns score the
# profile head against H3K27ac for BOTH arms, so for the multi-task arm they measure how
# well a head trained on DNase happens to predict H3K27ac shape. That is not what it was
# trained for and a drop there is expected rather than informative.
#
# Elements are the DNase-derived set, matching 1.11 and therefore both arms.
# --no-rc-average keeps these numbers comparable to every other table in the project; RC
# averaging is still opt-in project-wide (report Fig. 13).
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

# Fail before burning a GPU hour: a preempted fold leaves a checkpoint that looks finished.
for m in multimodal5p_accs5p_hw500_clw10 multimodal5p_accs5p_dnaseprof00561_hw500_clw10; do
    for f in 0 1 2 3 4; do
        [[ -s "$P/models/$m/fold$f/training_complete.json" ]] || {
            echo "ERROR: $m/fold$f has no training_complete.json, so it either did not "\
"finish or was preempted. Scoring it would silently mix a best-so-far checkpoint into "\
"the comparison." >&2
            exit 1
        }
    done
done

$PY scripts/2.15.perfold_from_config.py config/multitask_k562_configs.json \
    multitask_k562_ "$EL" --no-rc-average \
    --pair dnase_profile_head atac_input
echo MULTITASK_EVAL_DONE
