#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/errannot.%j.txt
#SBATCH -e log/errannot.%j.txt
#SBATCH --job-name=errannot
#
# Test the CTCF hypothesis directly against the prediction-error strata.
#
# The over-predicted tail is 10x more accessible than typical, barely more acetylated
# (ATAC/K27ac 13.25 vs 2.20), GC 0.59 and CpG o/e 0.55, and 2.4x promoter-enriched. That is
# a CpG-island/CTCF phenotype by composition. These are the actual ENCODE peak calls, so it
# is a direct test rather than a PWM proxy:
#   CTCF      ENCSR000AKO  ENCFF519CXF  optimal IDR-thresholded
#   EP300     ENCSR000EGE  ENCFF702XPO  IDR-thresholded
#   H3K4me1   ENCSR000AKS  ENCFF135ZLM  replicated
#   H3K27me3  ENCSR000AKQ  ENCFF323WOT  pseudoreplicated
#   H3K27ac   ENCSR000AKP  ENCFF544LXB  pseudoreplicated
#
# Prediction if the CTCF explanation holds: over-predicted enriched for CTCF and depleted for
# EP300; under-predicted the reverse, since EP300 marks the coactivator-bound enhancers whose
# acetylation is not readable from accessibility.

set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
W=${SCRATCH}/errannot_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

"$D/.pixi/envs/multimodal/bin/python" scripts/4.10.characterize_prediction_errors.py --emit-beds "$W" > /dev/null
ls -la "$W"/*.bed

printf '%-24s %8s' stratum n
for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac; do printf ' %9s' "$m"; done; printf '\n'

declare -A PK=( [CTCF]=ENCFF519CXF.bed.gz [EP300]=ENCFF702XPO.bed.gz
                [H3K4me1]=ENCFF135ZLM.bed.gz [H3K27me3]=ENCFF323WOT.bed.gz
                [H3K27ac]=ENCFF544LXB.bed.gz )
for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac; do
    zcat "$E/${PK[$m]}" | cut -f1-3 | sort -k1,1 -k2,2n -T "$W" > "$W/$m.bed"
done

for f in "$W"/stratum_*.bed; do
    lab=$(basename "$f" .bed | sed 's/^stratum_//')
    n=$(wc -l < "$f")
    printf '%-24s %8d' "$lab" "$n"
    for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac; do
        ov=$(bedtools intersect -u -a "$f" -b "$W/$m.bed" | wc -l)
        printf ' %8.1f%%' "$(echo "scale=4; 100*$ov/$n" | bc)"
    done
    printf '\n'
done
