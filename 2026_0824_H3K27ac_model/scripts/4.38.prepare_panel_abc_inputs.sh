#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --array=0-3
#SBATCH -o log/panel_abc_inputs.%A_%a.txt
#SBATCH --job-name=panel_abc_in
#
# Filtered read files for running ABC in the four non-K562 cell types of the p300 panel.
#
# WHY. The deployable neighbour model (4.36 fitted OUTSIDE K562) needs, in each training cell
# type, ABC's own candidate regions with ATAC and observed p300 counted over them, so the
# features it is fitted on are defined exactly as they are in K562. ENCODE ships these as
# unfiltered alignments; the same filters as 0.45 are applied (EP300: -F 1804 -q 30; ATAC:
# additionally -f 2), so ABC counts the reads the training tracks were built from.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
ST=$D/.pixi/envs/multimodal/bin/samtools
CTS=(GM12878 HepG2 A549 MCF-7)
CT=${CTS[$SLURM_ARRAY_TASK_ID]}
OUT=$D/2026_0824_H3K27ac_model/data/panel_abc/$CT
mkdir -p "$OUT"
for kind in ATAC EP300; do
    for b in $(ls "$E/$CT/$kind"/*.bam | grep -v filtered.sorted); do
        o="$OUT/${kind}_$(basename "${b%.bam}").filt.bam"
        [[ -s "$o.bai" ]] && { echo "exists: $o"; continue; }
        flags="-F 1804 -q 30"; [[ $kind == ATAC ]] && flags="-f 2 -F 1804 -q 30"
        $ST view -@ 8 -b $flags "$b" | $ST sort -@ 8 -m 2G -T "$SCRATCH/st_${CT}_$$" -o "$o" -
        $ST index -@ 8 "$o"
        echo "$CT $kind $(basename "$b") -> $($ST view -c -@ 8 "$o") reads"
    done
done
echo "done $CT"
