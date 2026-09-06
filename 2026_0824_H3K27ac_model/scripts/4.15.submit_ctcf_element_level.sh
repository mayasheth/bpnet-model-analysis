#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/ctcfelem.%j.txt
#SBATCH -e log/ctcfelem.%j.txt
#SBATCH --job-name=ctcfelem
#
# Per-element CTCF/IgG counts over all ABC candidate regions, then the element-level tests
# that the stratum-level mean enrichment in 4.13 cannot support. Same BAMs and the same
# mapped-read totals as 4.13, so the two are directly comparable.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562
PY="$D/.pixi/envs/multimodal/bin/python"
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
W=${SCRATCH}/ctcfelem_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

EL=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt
awk 'NR>1 {print $1"\t"$2"\t"$3}' "$EL" | sort -k1,1 -k2,2n > "$W/regions.bed"
echo "regions: $(wc -l < "$W/regions.bed")"

declare -A BAM=(
  [CTCF]=ENCFF216XRV.filtered.sorted.bam
  [IgG]=ENCFF396DTD.sorted.bam
)
ORDER=(CTCF IgG)
BAMS=(); TOT=()
for m in "${ORDER[@]}"; do
  BAMS+=("$E/${BAM[$m]}")
  n=$(samtools view -c -F 260 "$E/${BAM[$m]}")
  echo "$m mapped=$n"
  TOT+=("$m=$n")
done

bedtools multicov -bams "${BAMS[@]}" -bed "$W/regions.bed" > "$W/counts.bed"
echo "multicov done: $(wc -l < "$W/counts.bed") rows"

IFS=,; MARKS="${ORDER[*]}"; TOTALS="${TOT[*]}"; unset IFS
$PY scripts/4.15.ctcf_element_level.py "$W/counts.bed" "$MARKS" "$TOTALS"
