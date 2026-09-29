#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 3:00:00
#SBATCH --mem=48G
#SBATCH -c 16
#SBATCH --array=0-23
#SBATCH -o log/motifscan.%A_%a.txt
#SBATCH --job-name=motifscan
set -euo pipefail
CH=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python \
    scripts/4.44.scan_genome_motifs.py --chrom ${CH[$SLURM_ARRAY_TASK_ID]}
