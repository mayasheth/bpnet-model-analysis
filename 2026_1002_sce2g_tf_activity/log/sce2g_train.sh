#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity/log/sce2g_train.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity/log/sce2g_train.%j.txt
#SBATCH --job-name=sce2g_train
# scE2G main (7cb2af7) training-mode run: chromosome-CV predictions on the K562 CRISPR benchmark.
# CFG overrides the config (default: the baseline).
set -euo pipefail
# --rerun-triggers mtime: frag_to_tagAlign's params hold workflow.source_path(), a random
# /tmp source-cache path, so the default params trigger re-runs it (and everything after it)
# in every invocation, including inside each SLURM job.
TRIG="--rerun-triggers mtime"
export PYTHONUNBUFFERED=1
set +u; source ~/.bashrc
conda activate run_snakemake; set -u
CFG=${CFG:-/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity/config/config_training.yaml}
cd /oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G_main
snakemake -s workflow/Snakefile_training --configfile $CFG --use-conda --profile $HOME/.config/snakemake/slurm --unlock || true
snakemake -s workflow/Snakefile_training --configfile $CFG --use-conda --profile $HOME/.config/snakemake/slurm $TRIG -n -q
snakemake -s workflow/Snakefile_training --configfile $CFG --use-conda --profile $HOME/.config/snakemake/slurm $TRIG \
    --rerun-incomplete --keep-going --max-jobs-per-second 0.5 --max-status-checks-per-second 5
echo "scE2G training complete"
