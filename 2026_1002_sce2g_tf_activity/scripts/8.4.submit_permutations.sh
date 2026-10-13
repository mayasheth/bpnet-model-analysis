#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 1:30:00
#SBATCH --mem=24G
#SBATCH -c 6
#SBATCH --array=1-100%25
#SBATCH -o log/nested_fperm.%A_%a.txt
#SBATCH --job-name=nested_fperm
# 8.4 Phase 3 null (plans/8_crispr_feature_discovery.md): rerun 8.3's selection procedures with the
# candidate vectors shuffled among pairs within TSS-distance bins (labels and base untouched), one seed per
# array task. Replaces the within-gene label shuffle (invalid, 2026-10-06). Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/tobias/bin/python scripts/8.3.nested_selection.py --procedures forward,l1,family \
    --feature-perm-seed $SLURM_ARRAY_TASK_ID --n-boot 0 --jobs 5
