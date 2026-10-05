#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/insbam.%A_%a.txt
#SBATCH --array=1-23
#SBATCH --job-name=insbam
# 2.1 per chromosome into $SCRATCH; merge with: MERGE=1 sbatch --array=1 --dependency=afterok:<id> 2.1...
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
FRAG=/oak/stanford/groups/engreitz/Projects/scE2G/data/K562_260122/K562_Xu/atac_fragments.tsv.gz
OUT=$SCRATCH/tobias_k562_xu; mkdir -p $OUT
export PATH=$A2/.pixi/envs/tobias/bin:$PATH
cd $A2
if [[ "${MERGE:-0}" == 1 ]]; then
    $A2/.pixi/envs/tobias/bin/samtools merge -f -@ 2 $OUT/k562_xu_insertions.bam $(ls $OUT/chr*.bam | sort -V)
    $A2/.pixi/envs/tobias/bin/samtools index $OUT/k562_xu_insertions.bam
    $A2/.pixi/envs/tobias/bin/samtools idxstats $OUT/k562_xu_insertions.bam | head -n 25
    exit 0
fi
CHR=($(seq -f "chr%g" 1 22) chrX)
C=${CHR[$(( SLURM_ARRAY_TASK_ID - 1 ))]}
$A2/.pixi/envs/tobias/bin/python scripts/2.1.fragments_to_insertion_bam.py --frag $FRAG --chroms $C --out $OUT/$C.bam
