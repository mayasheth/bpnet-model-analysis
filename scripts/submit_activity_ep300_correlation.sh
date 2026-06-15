#!/bin/bash
#SBATCH -p normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH -o logs/activity_ep300_corr_%j.out
#SBATCH -e logs/activity_ep300_corr_%j.err
#SBATCH -J activity_ep300_corr

set -eo pipefail

cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet

mkdir -p logs

source ~/.bashrc
conda activate analysis

python3 scripts/plot_activity_ep300_correlation.py
