#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 00:10:00
#SBATCH --mem=4G
#SBATCH -c 1
#SBATCH -o log/shimtest.%j.txt
#SBATCH -e log/shimtest.%j.txt
#SBATCH --job-name=shimtest
#
# Reproduce, on a COMPUTE NODE, the two things that killed create_neighborhoods:
#   1. mamba must be resolvable from the shared shim (driver-local /tmp was invisible here)
#   2. the rule conda env's python must win, and must import pyranges
# Costs seconds; the alternative is another 37-job round trip per hypothesis.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
SM_ENV=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
CONDA_ENV=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/.snakemake/conda/a0febb704d007e3c12f577f28bdeccb8_
SHIM=$D/2026_0824_H3K27ac_model/.mamba_shim
export PATH="$PATH:$SHIM"

echo "node: $(hostname)"
echo "1) mamba on PATH: $(command -v mamba || echo MISSING)"
echo "   conda on PATH: $(command -v conda || echo MISSING)"
echo "2) bare python:   $(command -v python || echo none)"
echo "3) env python:    $CONDA_ENV/bin/python"
"$CONDA_ENV/bin/python" -c 'import pyranges, pyBigWig, pandas; print("   imports OK, pyranges", pyranges.__version__)'
echo "4) after conda activate, which python wins:"
set +u
source "$SM_ENV/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "$CONDA_ENV" 2>/dev/null || echo "   (conda activate unavailable; snakemake uses its own activation)"
echo "   python -> $(command -v python)"
python -c 'import pyranges; print("   pyranges importable via activated env")' 2>&1 | tail -1
set -u
echo "ALL SHIM CHECKS DONE"
