#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 6:00:00
#SBATCH --mem=128G
#SBATCH -c 4
#SBATCH -o log/panel_div.%j.txt
#SBATCH --job-name=panel_div
#
# Download the candidates' ENCODE peak files (config/panel_candidates*.tsv; 47 + 77 MB, approved by
# Maya 2026-09-29) into $SCRATCH, then rank cell types for the p300 training panel (4.49).
set -euo pipefail
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
MAN=${MANIFEST:-config/panel_candidates_loose.tsv}
PK=$SCRATCH/panel_candidates
cd "$P"; mkdir -p "$PK" log
tail -n +2 "$MAN" | cut -f4 | while read acc; do
    f=$PK/$acc.bed.gz
    [[ -s "$f" ]] || curl -sSfL -o "$f" "https://www.encodeproject.org/files/$acc/@@download/$acc.bed.gz"
    echo "$acc $(stat -c %s "$f") bytes"
done
$PY scripts/4.49.panel_tf_diversity.py --manifest "$MAN" --peak-dir "$PK" --out-prefix ${OUT_PREFIX:-results/2026_0929_panel_tf_diversity_loose}
