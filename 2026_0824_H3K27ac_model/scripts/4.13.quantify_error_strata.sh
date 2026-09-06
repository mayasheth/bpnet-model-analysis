#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -o log/errquant.%j.txt
#SBATCH -e log/errquant.%j.txt
#SBATCH --job-name=errquant
#
# QUANTITATIVE signal in the prediction-error strata, to sit alongside 4.11's peak overlaps.
#
# Peak calls are thresholded and the two measures demonstrably diverge here: the
# over-predicted stratum has H3K27ac peaks called at ~2x background while its H3K27ac RPM is
# AT background. There is substantial signal without a call and calls without much signal, so
# neither measure alone is trustworthy and both get reported.
#
# Counts reads from the local filtered BAMs rather than downloading portal bigwigs -- same
# data, no transfer. Normalised as RPKM (per million mapped, per kb) so strata of different
# size and elements of different width are comparable.

set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
W=${SCRATCH}/errquant_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

"$D/.pixi/envs/multimodal/bin/python" scripts/4.10.characterize_prediction_errors.py \
    --emit-beds "$W" > /dev/null

# one representative replicate per mark; adding replicates would average out, not change rank
declare -A BAM=(
  [CTCF]=ENCFF216XRV.filtered.sorted.bam
  [EP300]=ENCFF466WKF.filtered.sorted.bam
  [H3K4me1]=ENCFF204MWI.filtered.sorted.bam
  [H3K27me3]=ENCFF351YGP.filtered.sorted.bam
  [H3K27ac]=ENCFF790GFL.se.filtered.sorted.bam
  [IgG]=ENCFF396DTD.sorted.bam
)

echo "mapped-read totals (RPKM denominators)"
declare -A TOT
for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac IgG; do
    b="$E/${BAM[$m]}"
    [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }
    TOT[$m]=$(samtools view -c -F 0x400 -@ 4 "$b")
    printf '  %-10s %12d\n' "$m" "${TOT[$m]}"
done

echo
OUT_TSV="$P/results/error_strata_rpkm.tsv"
{ printf 'stratum\tn'; for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac IgG; do printf '\t%s' "$m"; done; printf '\n'; } > "$OUT_TSV"
{ printf 'mark\tmapped_reads\n'; for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac IgG; do printf '%s\t%s\n' "$m" "${TOT[$m]}"; done; } > "$P/results/error_strata_rpkm_denominators.tsv"
printf '%-26s %8s' stratum n
for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac IgG; do printf ' %9s' "$m"; done; printf '\n'
for f in "$W"/stratum_*.bed; do
    lab=$(basename "$f" .bed | sed 's/^stratum_//')
    n=$(wc -l < "$f")
    kb=$(awk '{s += ($3-$2)} END {printf "%.6f", s/1000}' "$f")
    printf '%-26s %8d' "$lab" "$n"
    printf '%s\t%s' "$lab" "$n" >> "$OUT_TSV"
    for m in CTCF EP300 H3K4me1 H3K27me3 H3K27ac IgG; do
        c=$(samtools view -c -F 0x400 -@ 4 -L "$f" "$E/${BAM[$m]}")
        v=$(echo "scale=6; 1000000*$c/(${TOT[$m]}*$kb)" | bc)
        printf ' %9.2f' "$v"
        printf '\t%s' "$v" >> "$OUT_TSV"
    done
    printf '\n'
    printf '\n' >> "$OUT_TSV"
done
echo "wrote $OUT_TSV"
echo
echo "RPKM per stratum. IgG is the background control: a mark whose pattern tracks IgG is"
echo "reporting chromatin accessibility to the antibody rather than the mark itself."
