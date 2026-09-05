#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH -c 16
#SBATCH -o log/crispr.%j.txt
#SBATCH -e log/crispr.%j.txt
#SBATCH --job-name=crispr_cmp
#
# Benchmark the eleven predicted-activity ABC arms against CRISPR data.
#
# RUN FROM workflow/, NOT THE REPO ROOT. Snakemake resolves .snakemake/conda relative to the
# working directory, and this repo's built R env lives at
#   workflow/.snakemake/conda/a5c95d47c3f0a63f7a24cdab52a16531_   (R 4.1.1)
# Running from the repo root makes .snakemake/conda look empty, which is why an earlier
# attempt concluded the env had to be rebuilt. It also explains the relative
# `../resources/cell_type_mapping/...` paths in this repo's older configs: relative to
# workflow/, `..` is the repo root.
#
# ENVIRONMENT.
#   snakemake -> run_snakemake (7.32.4), called by ABSOLUTE PATH so its bin never lands on
#                PATH and cannot shadow the rule env's R or python
#   R         -> supplied by --use-conda from the prebuilt env above; nothing R-related goes
#                on PATH by hand
#   mamba     -> exposed through a shim dir APPENDED to PATH. --use-conda checks for a conda
#                frontend even when the env already exists, and appending means it cannot
#                shadow anything. The shim lives on Oak, not local /tmp, so rule jobs on any
#                node can see it.
#
# PREFLIGHT first: assert the prebuilt env is the one Snakemake will pick and that it can
# load every R package the pipeline needs.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
C=/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison
SM_ENV=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
SM=$SM_ENV/bin/snakemake
CONDA_ENV=$C/workflow/.snakemake/conda/a5c95d47c3f0a63f7a24cdab52a16531_
CFG=$C/config/config_predicted_activity.yml
mkdir -p "$D/2026_0824_H3K27ac_model/log"

SHIM=$D/2026_0824_H3K27ac_model/.mamba_shim
mkdir -p "$SHIM"
ln -sfn "$SM_ENV/bin/mamba" "$SHIM/mamba"
ln -sfn "$SM_ENV/bin/conda" "$SHIM/conda"
export PATH="$PATH:$SHIM"
command -v mamba >/dev/null || { echo "ERROR: mamba not on PATH via $SHIM" >&2; exit 1; }

cd "$C/workflow"

echo "=== preflight ==="
echo "node:      $(hostname)"
echo "snakemake  $($SM --version)"
echo "workdir    $(pwd)"
[[ -d "$CONDA_ENV" ]] || { echo "ERROR: prebuilt env missing: $CONDA_ENV" >&2; exit 1; }
echo "rule env   $CONDA_ENV"
"$CONDA_ENV/bin/Rscript" -e 'libs <- c("tidyverse","data.table","ROCR","caTools","boot","DT","rmarkdown","bookdown","optparse","cowplot","ggpubr");
            miss <- libs[!sapply(libs, requireNamespace, quietly=TRUE)];
            if (length(miss)) { cat("MISSING R PACKAGES:", paste(miss, collapse=", "), "\n"); quit(status=1) };
            cat("all", length(libs), "R packages load OK in the prebuilt env\n")'

echo
echo "=== dry run (same flags as the run) ==="
$SM --configfile "$CFG" --use-conda -n -q

echo
echo "=== run ==="
$SM --configfile "$CFG" --use-conda -j 4 --rerun-incomplete --nolock
echo "CRISPR comparison complete"
