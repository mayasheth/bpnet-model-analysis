#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -o log/ep300cmp.%j.txt
#SBATCH -e log/ep300cmp.%j.txt
#SBATCH --job-name=ep300cmp
#
# Why does a K562-trained p300 model transfer to GM12878 (+0.207 over the target's floor)
# while a GM12878-trained one transfers to K562 not at all (+0.006, p=0.72), when both fit
# their own cell type equally well (+0.313 and +0.309 over matched floors)?
#
# The H3K27ac explanation for its own transfer asymmetry -- that GM12878 is the harder cell
# type -- cannot apply, because here both cell types are equally predictable locally. The
# leading remaining candidate is that the K562 EP300 experiment supports a more portable
# sequence model than GM12878's despite equal local fit, which would show up as a difference
# in the training signal's quality rather than its predictability.
#
# Measures, for ENCSR000EGE (K562) and ENCSR000DZG (GM12878): peak count, peak width, mapped
# depth, and FRiP. FRiP is the informative one -- it is the fraction of reads inside peaks, so
# it bounds how much of the training signal is signal rather than background, and a model
# trained on a low-FRiP target can fit its own cell type by memorising accessibility while
# learning little portable sequence.
#
# This is diagnostic, not decisive: two experiments cannot establish a general relationship
# between FRiP and portability. It either supports the hypothesis enough to justify a third
# cell type, or removes it and sends the question elsewhere.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
W=${SCRATCH}/ep300cmp_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$D/2026_0824_H3K27ac_model"

report () {
    local label=$1 peaks=$2; shift 2
    local bams=("$@")
    # normalise the peak file to plain sorted BED regardless of gz/narrowPeak
    if [[ "$peaks" == *.gz ]]; then zcat "$peaks"; else cat "$peaks"; fi \
        | cut -f1-3 | sort -k1,1 -k2,2n > "$W/$label.bed"
    local n width
    n=$(wc -l < "$W/$label.bed")
    width=$(awk '{s += ($3-$2)} END {printf "%.0f", s/NR}' "$W/$label.bed")
    echo "=========== $label"
    printf '  peaks            %d\n  mean peak width  %s bp\n' "$n" "$width"
    printf '  total peak bp    %s\n' "$(awk '{s += ($3-$2)} END {print s}' "$W/$label.bed")"
    local tot=0 inpk=0
    for b in "${bams[@]}"; do
        [[ -s "$b" ]] || { echo "  MISSING BAM: $b" >&2; continue; }
        local t i
        t=$(samtools view -c -F 0x400 -@ 8 "$b")
        i=$(samtools view -c -F 0x400 -@ 8 -L "$W/$label.bed" "$b")
        printf '  %-28s mapped %12d  in peaks %12d  FRiP %.4f\n' \
            "$(basename "$b" | cut -c1-28)" "$t" "$i" "$(echo "scale=6; $i/$t" | bc)"
        tot=$((tot + t)); inpk=$((inpk + i))
    done
    printf '  POOLED           mapped %12d  in peaks %12d  FRiP %.4f\n' \
        "$tot" "$inpk" "$(echo "scale=6; $inpk/$tot" | bc)"
}

report K562_ENCSR000EGE \
    "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    "$E/K562/ENCFF466WKF.filtered.sorted.bam" \
    "$E/K562/ENCFF163FSR.filtered.sorted.bam"

report GM12878_ENCSR000DZG \
    "$E/GM12878/EP300/ENCFF926AKK.bed.gz" \
    "$E/GM12878/EP300/ENCFF515HYM.filtered.sorted.bam" \
    "$E/GM12878/EP300/ENCFF215GSQ.filtered.sorted.bam"

echo
echo "A higher FRiP means more of the training signal sits in peaks. If K562's is materially"
echo "higher, that supports the portability hypothesis and justifies a third cell type; if the"
echo "two are comparable, the asymmetry is not explained by target quality and the next"
echo "candidate is the peak set used to define training positives."
