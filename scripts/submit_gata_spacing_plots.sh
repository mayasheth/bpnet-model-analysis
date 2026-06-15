#!/bin/bash
#SBATCH -J gata_spacing_plots
#SBATCH -p normal
#SBATCH --time=00:30:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/motif_spacing/slurm_spacing_plots_%j.out
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/motif_spacing/slurm_spacing_plots_%j.err

set -euo pipefail

PROJ=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet

module load devel pixi/0.53.0

cd "$PROJ"

echo "Generating GATA n-copy spacing plots (n=2,3,4) ..."
pixi run -e ism python scripts/plot_gata_ncopy_spacing.py --n 2 3 4

echo "All done."
