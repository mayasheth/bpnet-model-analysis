#!/bin/bash
#SBATCH -J finemo_tf_corr
#SBATCH -p normal
#SBATCH --time=02:00:00
#SBATCH --mem=48GB
#SBATCH --cpus-per-task=4
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/tf_correlation/slurm_%j.out
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/tf_correlation/slurm_%j.err

set -euo pipefail

PROJ=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
FINEMO_DIR="$PROJ/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2"

module load devel pixi/0.53.0

mkdir -p "$FINEMO_DIR/tf_correlation"

cd "$PROJ"
pixi run python scripts/finemo_tf_correlation.py \
    --finemo-dir "$FINEMO_DIR"
