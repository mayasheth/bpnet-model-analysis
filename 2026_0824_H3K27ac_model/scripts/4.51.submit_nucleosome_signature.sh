#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 16:00:00
#SBATCH --mem=24G
#SBATCH -c 2
#SBATCH -o log/nucsig.%A_%a.txt
#SBATCH --job-name=nucsig
#SBATCH --array=0-4
#
# 4.51 features, one cell type per task (fragment scan of that cell's ATAC BAMs). The analyze
# stage is submitted separately with --dependency=afterok on this array.
set -euo pipefail
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
CELLS=(K562 GM12878 HepG2 A549 MCF-7)
cd "$P"
$PY scripts/4.51.nucleosome_signature.py features --cell "${CELLS[${SLURM_ARRAY_TASK_ID:-0}]}"
