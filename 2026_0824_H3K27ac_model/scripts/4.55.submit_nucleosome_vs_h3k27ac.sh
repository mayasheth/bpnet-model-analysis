#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 3:00:00
#SBATCH --mem=24G
#SBATCH -c 2
#SBATCH -o log/nuck27.%j.txt
#SBATCH --job-name=nuck27
#
# 4.55: nucleosome features (4.51) against H3K27ac vs p300 in K562 and GM12878.
set -euo pipefail
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
cd "$P"
$PY scripts/4.55.nucleosome_vs_h3k27ac.py
