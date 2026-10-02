#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/k4me1_tracks.%A_%a.txt
#SBATCH -e log/k4me1_tracks.%A_%a.txt
#SBATCH --job-name=k4me1_tracks
#SBATCH --array=1-5
#
# Build H3K4me1 stranded 5'-end tracks for the five-cell panel (K562, GM12878, HepG2, A549,
# MCF-7), merged and per replicate, by 0.14's make_5p code path so they are comparable with the
# H3K27ac targets.
#
# FILTERING MATCHES K562 AND GM12878. Their H3K4me1 BAMs were filtered by
# Data/scripts/process_single_ended_bam.sh (-F 780 -q 30). The 2026-10-01 downloads are
# ENCODE4 unfiltered single-end alignments with no duplicate flags (checked on the first 200k
# records of each), so the same -F 780 -q 30 is applied here; -F 1804 would be identical on
# these files. Duplicates are not removed in any cell.
#
# Filtered BAMs for the new cells stay on $SCRATCH beside the raw downloads (Oak is 97% full);
# the raw files are re-downloadable by accession. Tracks go to data/h3k4me1_panel/.
#
# Usage: sbatch 0.47.build_h3k4me1_tracks.sh
#        CELLS=K562 sbatch --array=1 0.47.build_h3k4me1_tracks.sh   # rebuild one
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
DATA=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
RAW=$SCRATCH/h3k4me1
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT="$P/data/h3k4me1_panel"; mkdir -p "$OUT" "$P/log"
cd "$P"

IFS=',' read -ra CELL_LIST <<< "${CELLS:-K562,GM12878,HepG2,A549,MCF-7}"
CT=${CELL_LIST[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
SAFE=$(echo "$CT" | tr '-' '_' | tr '[:upper:]' '[:lower:]')
WORK="$SCRATCH/k4me1_${SAFE}_${SLURM_JOB_ID:-$$}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
echo "=== $CT -> $OUT (work $WORK) ==="
[[ -s "$CHR" ]] || { echo "ERROR: missing chrom sizes $CHR" >&2; exit 1; }

case "$CT" in
    K562)    ACCS=(ENCFF524BOJ ENCFF204MWI ENCFF665JSC) ;;
    GM12878) ACCS=(ENCFF757IRH ENCFF579EEP) ;;
    HepG2)   ACCS=(ENCFF748FEE ENCFF699BSB) ;;
    A549)    ACCS=(ENCFF563VWO ENCFF050COQ ENCFF162FGL) ;;
    MCF-7)   ACCS=(ENCFF987GWM ENCFF379DBF) ;;
    *) echo "ERROR: unknown cell $CT" >&2; exit 1 ;;
esac

BAMS=()
for acc in "${ACCS[@]}"; do
    pre="$DATA/$CT/$acc.filtered.sorted.bam"
    if [[ -s "$pre" ]]; then
        BAMS+=("$pre"); echo "  using pre-filtered $(basename "$pre")"; continue
    fi
    raw="$RAW/$CT/$acc.bam"
    f="$RAW/$CT/$acc.filtered.sorted.bam"
    [[ -s "$raw" ]] || { echo "ERROR: missing $raw" >&2; exit 1; }
    if [[ ! -s "$f.bai" ]]; then
        samtools view -@ 4 -u -F 780 -q 30 "$raw" | samtools sort -@ 4 -m 2G -T "$WORK/st" -o "$f.tmp" -
        mv "$f.tmp" "$f"; samtools index -@ 4 "$f"
    fi
    BAMS+=("$f")
    echo "  filtered $acc: $(samtools view -c -@ 4 "$f") of $(samtools view -c -@ 4 "$raw") reads"
done

make_5p () {  # make_5p <bam> <prefix>
    local bam="$1" prefix="$2"
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_5p_${tag}.bw"
        [[ -s "$out" ]] && { echo "  exists: $(basename "$out")"; continue; }
        bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
          | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
          | sort -k1,1 -k2,2n -S 4G -T "$WORK" > "$WORK/x.bg"
        bedGraphToBigWig "$WORK/x.bg" "$CHR" "$out.tmp"; mv "$out.tmp" "$out"; rm -f "$WORK/x.bg"
        echo "  wrote $(basename "$out")"
    done
}

if [[ ${#BAMS[@]} -gt 1 ]]; then
    samtools merge -@ 4 -f "$WORK/merged.bam" "${BAMS[@]}"
    samtools index -@ 4 "$WORK/merged.bam"
    MERGED="$WORK/merged.bam"
else
    MERGED="${BAMS[0]}"
fi
echo "  merged: $(samtools view -c -@ 4 "$MERGED") reads"

make_5p "$MERGED" "${SAFE}_h3k4me1"
i=1
for b in "${BAMS[@]}"; do make_5p "$b" "${SAFE}_h3k4me1_rep${i}"; i=$((i+1)); done

echo "=== done $CT ==="
ls -l "$OUT" | grep "${SAFE}_" || true
