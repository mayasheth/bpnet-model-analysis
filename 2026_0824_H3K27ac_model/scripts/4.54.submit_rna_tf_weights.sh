#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH -o log/rna_tf.%j.txt
#SBATCH --job-name=rna_tf
#
# Download the ENCODE poly(A) RNA-seq gene quantifications (10 files, ~110 MB; approved by Maya
# 2026-09-30) into $SCRATCH/rnaseq, then build and score RNA-based TF lists (4.54).
set -euo pipefail
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
OUT=$SCRATCH/rnaseq; mkdir -p "$OUT"; cd "$P"
for acc in ENCFF742CVV ENCFF222UVT ENCFF791MED ENCFF473KMX ENCFF831QQF ENCFF168QKW ENCFF855AKQ ENCFF244DNJ ENCFF967AOT ENCFF921PJP; do
    [[ -s "$OUT/$acc.tsv" ]] || curl -sSfL -o "$OUT/$acc.tsv" "https://www.encodeproject.org/files/$acc/@@download/$acc.tsv"
    echo "$acc $(wc -l < "$OUT/$acc.tsv") lines"
done
$PY scripts/4.54.rna_tf_weights.py --rna-dir "$OUT"
