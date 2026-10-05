#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/scprinter_fetch.%j.txt
#SBATCH --job-name=scprinter_fetch
# 3.0 Fetch scPrinter's hg38 (v2) Tn5 bias model (hg38_bias_v2.h5, 12.8 GB) and GENCODE v41 fasta
# into $SCRATCH/scprinter_cache through scPrinter's own pooch registry (md5-checked). Maya OK'd 2026-10-03.
set -euo pipefail
export PYTHONUNBUFFERED=1
export PIXI_CACHE_DIR=$SCRATCH/.cache/rattler UV_CACHE_DIR=$SCRATCH/.cache/uv
export SCPRINTER_DATA=$SCRATCH/scprinter_cache
mkdir -p $SCPRINTER_DATA
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/scprinter/bin/python - <<'EOF'
import scprinter as scp
g = scp.genome.hg38
print("bias:", g.fetch_bias())
print("fasta:", g.fetch_fa())
EOF
ls -la $SCPRINTER_DATA
du -sh $SCPRINTER_DATA
