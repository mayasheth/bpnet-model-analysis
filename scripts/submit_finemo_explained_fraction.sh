#!/bin/bash
#SBATCH -p normal
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH -o logs/finemo_explained_fraction_%j.out
#SBATCH -e logs/finemo_explained_fraction_%j.err
#SBATCH -J finemo_explained_fraction

set -eo pipefail

cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet

mkdir -p logs

source ~/.bashrc
conda activate analysis

python3 scripts/plot_finemo_explained_fraction.py
