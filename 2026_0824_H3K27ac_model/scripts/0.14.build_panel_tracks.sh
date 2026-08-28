#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 12:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/panel_tracks.%j.txt
#SBATCH -e log/panel_tracks.%j.txt
#SBATCH --job-name=k27_panel_tracks
#
# Build the tracks needed to extend the analysis to additional cell types.
#
# Two different requirements, because target and input are not held to the same standard:
#   H3K27ac is the TARGET  -> needs >=2 replicates of one processing so an inter-replicate
#                             ceiling can be computed. Built as stranded 5'-end counts
#                             (genomecov -5 -dz), merged plus per-replicate.
#   DNase is an INPUT      -> only needs a usable track; a single replicate is fine, since
#                             no ceiling is computed on an input. Built as read coverage
#                             (genomecov -bg), matching the accessibility convention used
#                             for ATAC, with -pc for paired-end libraries.
#
# Usage: sbatch 0.14.build_panel_tracks.sh CELLTYPE ASSAY RUNTYPE BAM[,BAM...]
#   e.g. sbatch 0.14.build_panel_tracks.sh H9 H3K27ac se /path/a.bam,/path/b.bam

set -euo pipefail
export PYTHONUNBUFFERED=1

CT=${1:?usage: CELLTYPE ASSAY RUNTYPE BAMS}
ASSAY=${2:?}
RUNTYPE=${3:?}
BAMS_CSV=${4:?}

D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT="$P/data/panel"; mkdir -p "$OUT" "$P/log"
WORK="${SCRATCH}/panel_${CT}_${ASSAY}_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT

IFS=',' read -ra BAMS <<< "$BAMS_CSV"
for b in "${BAMS[@]}"; do [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }; done
SAFE=$(echo "$CT" | tr -d ' ' | tr '-' '_')

# --- H3K27ac target: stranded 5' ends -------------------------------------
make_5p () {
    local bam="$1" prefix="$2"
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_5p_${tag}.bw"
        [[ -s "$out" ]] && { echo "  exists: $(basename $out)"; continue; }
        bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
          | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
          | sort -k1,1 -k2,2n -T "$WORK" > "$WORK/x.bg"
        bedGraphToBigWig "$WORK/x.bg" "$CHR" "$out"; rm -f "$WORK/x.bg"
        echo "  wrote $(basename $out)"
    done
}

# --- DNase input: read coverage -------------------------------------------
make_cov () {
    local bam="$1" out="$2"
    [[ -s "$out" ]] && { echo "  exists: $(basename $out)"; return; }
    local pc=()
    [[ "$RUNTYPE" == *pe* && "$RUNTYPE" != *se* ]] && pc=(-pc)   # fragment coverage for PE
    bedtools genomecov -ibam "$bam" -bg "${pc[@]+${pc[@]}}" \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
      | sort -k1,1 -k2,2n -T "$WORK" > "$WORK/x.bg"
    bedGraphToBigWig "$WORK/x.bg" "$CHR" "$out"; rm -f "$WORK/x.bg"
    echo "  wrote $(basename $out)"
}

echo "=== $CT $ASSAY (run_type=$RUNTYPE, ${#BAMS[@]} replicate file(s)) ==="
if [[ ${#BAMS[@]} -gt 1 ]]; then
    samtools merge -@ 4 -f "$WORK/merged.bam" "${BAMS[@]}"
    samtools index "$WORK/merged.bam"
    MERGED="$WORK/merged.bam"
else
    MERGED="${BAMS[0]}"
fi

if [[ "$ASSAY" == "H3K27ac" ]]; then
    make_5p "$MERGED" "${SAFE}_h3k27ac"
    i=1
    for b in "${BAMS[@]}"; do make_5p "$b" "${SAFE}_h3k27ac_rep${i}"; i=$((i+1)); done
else
    make_cov "$MERGED" "$OUT/${SAFE}_dnase.bw"
fi
echo "done: $CT $ASSAY"
