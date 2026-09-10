#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 3:00:00
#SBATCH --mem=48G
#SBATCH -o log/backfill_n_topq.%j.txt
#SBATCH -e log/backfill_n_topq.%j.txt
#SBATCH --job-name=n_topq
#SBATCH -n 1
#
# Recover `n_topq` for the four 2.15 tables that feed a figure with a top-quintile panel.
# CPU only -- no model is evaluated, only `model.trimming` is read (see 2.30's docstring).
#
# The figures these correct:
#   prof_residual_grid_  fig13, both panels are top-quintile metrics
#   fragchan_k562_       fig12 panel b
#   wide_k562_           fig10, the two K562 panels
#   wide_gm12878_        fig10, the two GM12878 panels
#
# Every table is written in place only after the recomputed `n` matches the stored `n`.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

run () {  # config, out_prefix, elements
    echo "########## $2 ##########"
    $PY scripts/2.30.backfill_n_topq.py "config/$1" "$2" "$3" ${DRY:-}
    echo
}

run prof_residual_grid_configs.json prof_residual_grid_ "$K5EL"
run fragchan_k562_configs.json      fragchan_k562_      "$K5EL"
run wide_k562_configs.json          wide_k562_          "$K5EL"
run wide_gm12878_configs.json       wide_gm12878_       "$GMEL"

echo "done"
