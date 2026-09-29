#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH --array=0-4
#SBATCH -o log/acttrack.%A_%a.txt
#SBATCH --job-name=acttrack
set -euo pipefail
CELLS=(K562 GM12878 HepG2 A549 MCF-7)
/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python \
    scripts/4.45.build_activity_track.py --cell ${CELLS[$SLURM_ARRAY_TASK_ID]} --weighting ${WEIGHTING:-oracle}
