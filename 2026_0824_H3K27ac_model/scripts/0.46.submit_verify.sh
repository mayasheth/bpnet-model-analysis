#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 1:00:00
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH -o log/verify046.%j.txt
#SBATCH -e log/verify046.%j.txt
#SBATCH --job-name=verify_tracks
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $D/2026_0824_H3K27ac_model
$D/.pixi/envs/multimodal/bin/python scripts/0.46.verify_panel_tracks.py   --panel-dir data/p300_panel --rebuilt-prefix GM12878
