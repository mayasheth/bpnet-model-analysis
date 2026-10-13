#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH -o log/feature_screen.%j.txt
#SBATCH --job-name=feature_screen
# 7.1 library then 7.2 screen (plans/7_feature_screen.md). Run from the A2 dir.
# SMOKE=chr21+chr22 ("+": sbatch --export splits on commas) builds and screens those chromosomes only, into $SCRATCH/feature_screen/smoke_*.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
PY=$A2/.pixi/envs/tobias/bin/python
export PATH=$A2/.pixi/envs/tobias/bin:$PATH
cd $A2
if [[ -n "${SMOKE:-}" ]]; then
    SMOKE=${SMOKE//+/,}
    LIB=$SCRATCH/feature_screen/smoke_${SMOKE//,/_}
    $PY scripts/7.1.build_feature_library.py --chroms $SMOKE
    $PY scripts/7.2.screen_features.py --lib-dir $LIB --out-prefix $LIB/screen --n-stab 5 --n-boot 20
else
    [[ -e $SCRATCH/feature_screen/library.elements.tsv.gz ]] && { echo "library exists; not overwriting" >&2; exit 1; }
    $PY scripts/7.1.build_feature_library.py
    $PY scripts/7.2.screen_features.py
fi
