#!/bin/bash
#SBATCH -p normal
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH -o logs/composite_figure_%j.out
#SBATCH -e logs/composite_figure_%j.err
#SBATCH -J composite_figure

set -eo pipefail

cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet

mkdir -p logs

module load devel pixi/0.53.0

pixi run python scripts/plot_finemo_composite_figure.py
