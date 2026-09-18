#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 24:00:00
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH -o log/p300panel.%A_%a.txt
#SBATCH -e log/p300panel.%A_%a.txt
#SBATCH --job-name=p300_tracks
#SBATCH --array=1-3
#
# Build EP300 target tracks and ATAC input tracks for the new panel cell types.
#
# CONSTRUCTION MUST MATCH THE EXISTING CELL TYPES OR THE PANEL IS CONFOUNDED. K562 and
# GM12878 are not rebuilt here, so every choice below is made to reproduce how theirs were
# built rather than to be independently sensible. 0.31 made exactly this argument when it
# built a full-depth control track by the same code path as its subsampled one: without it,
# a comparison confounds the thing under test with how the tracks were made.
#
# EP300, the TARGET. Stranded 5'-end counts, `genomecov -5 -dz -strand`, matching 0.14's
# treatment of H3K27ac. ENCODE ships these experiments as UNFILTERED alignments, and
# reference/DATA_INVENTORY.md's convention is to filter locally, so -F 1804 -q 30 is applied
# here; the existing K562 and GM12878 EP300 BAMs were filtered the same way before use.
#
# ATAC, the INPUT. The existing p300 models read `atac.bw`, built by
# 2026_0529.../0.1.make_accessibility_bigwig.sh from *.tn5.sorted.tagAlign.gz. Those
# tagAligns are Tn5-shifted, and 0.22 calibrated the shift that reproduces them from raw PE
# BAMs as plus +4 / minus -5, at Pearson r = 1.0000 exactly, independently in K562 and
# GM12878. The downloads here are BAMs, so the shift is applied on conversion. Getting it
# wrong would put every new cell type a few bp out of phase with the two existing ones.
#
# Both a coverage track and a 5'-insertion track are written. The p300 models use coverage
# (`atac.bw`); the 5' track is built because it costs one extra genomecov and the H3K27ac
# line uses that form.
#
# Usage: sbatch 0.45.build_p300_panel_tracks.sh          # array over the three new cells
#        CELLS=K562 sbatch --array=1 0.45...             # rebuild one, for the control
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
DATA=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT="$P/data/p300_panel"; mkdir -p "$OUT" "$P/log"
cd "$P"

IFS=',' read -ra CELL_LIST <<< "${CELLS:-A549,HepG2,MCF-7}"
CT=${CELL_LIST[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
SAFE=$(echo "$CT" | tr '-' '_')
WORK="${SCRATCH:-/tmp}/p300panel_${SAFE}_${SLURM_JOB_ID:-$$}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
echo "=== $CT -> $OUT (work $WORK) ==="

[[ -s "$CHR" ]] || { echo "ERROR: missing chrom sizes $CHR" >&2; exit 1; }

bg_to_bw () {  # bg_to_bw <bedgraph> <out.bw>
    sort -k1,1 -k2,2n -S 4G -T "$WORK" "$1" > "$1.sorted"
    bedGraphToBigWig "$1.sorted" "$CHR" "$2"
    rm -f "$1" "$1.sorted"
    echo "  wrote $(basename "$2")"
}

# ---------------------------------------------------------------- EP300 target
EP_BAMS=($(ls "$DATA/$CT/EP300"/*.bam 2>/dev/null | grep -v 'filtered.sorted' || true))
[[ ${#EP_BAMS[@]} -gt 0 ]] || { echo "ERROR: no EP300 bams for $CT" >&2; exit 1; }
echo "EP300: ${#EP_BAMS[@]} replicate file(s)"

FILT=()
for b in "${EP_BAMS[@]}"; do
    o="$WORK/$(basename "${b%.bam}").filt.bam"
    # -F 1804: unmapped, mate-unmapped, secondary, QC-fail, duplicate. -q 30: unique.
    samtools view -@ 4 -b -F 1804 -q 30 "$b" | samtools sort -@ 4 -T "$WORK/st" -o "$o" -
    FILT+=("$o")
    echo "  filtered $(basename "$b") -> $(samtools view -c -@ 4 "$o") reads"
done
if [[ ${#FILT[@]} -gt 1 ]]; then
    samtools merge -@ 4 -f "$WORK/ep300.bam" "${FILT[@]}"
else
    cp "${FILT[0]}" "$WORK/ep300.bam"
fi
samtools index -@ 4 "$WORK/ep300.bam"
echo "  merged EP300: $(samtools view -c -@ 4 "$WORK/ep300.bam") reads"

for strand in + -; do
    [[ "$strand" == "+" ]] && tag=plus || tag=minus
    out="$OUT/${SAFE}_ep300_5p_${tag}.bw"
    if [[ -s "$out" ]]; then echo "  exists: $(basename "$out")"; continue; fi
    bedtools genomecov -ibam "$WORK/ep300.bam" -5 -dz -strand "$strand" \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
      > "$WORK/ep.bg"
    bg_to_bw "$WORK/ep.bg" "$out"
done

# ---------------------------------------------------------------- ATAC input
AT_BAMS=($(ls "$DATA/$CT/ATAC"/*.bam 2>/dev/null || true))
if [[ ${#AT_BAMS[@]} -eq 0 ]]; then
    echo "no ATAC bams for $CT; skipping the input track"
else
    echo "ATAC: ${#AT_BAMS[@]} replicate file(s)"
    TA="$WORK/${SAFE}.tn5.tagAlign"
    : > "$TA"
    for b in "${AT_BAMS[@]}"; do
        # -f 2 keeps properly paired reads; ENCODE's ATAC `alignments` output is already
        # filtered, so -F 1804 -q 30 is belt-and-braces rather than a second opinion.
        # The Tn5 shift is +4 on the plus strand and -5 on the minus, which 0.22 measured
        # against the existing tagAligns at r = 1.0000.
        samtools view -@ 4 -b -f 2 -F 1804 -q 30 "$b" \
          | bedtools bamtobed -i - \
          | awk 'BEGIN{OFS="\t"}
                 {if($6=="+"){$2=$2+4}else{$3=$3-5}
                  if($2<0)$2=0; if($3<=$2)$3=$2+1;
                  print $1,$2,$3,"N",1000,$6}' >> "$TA"
        echo "  converted $(basename "$b")"
    done
    echo "  tagAlign reads: $(wc -l < "$TA")"

    out_cov="$OUT/${SAFE}_atac.bw"
    if [[ -s "$out_cov" ]]; then echo "  exists: $(basename "$out_cov")"; else
        bedtools genomecov -i <(sort -k1,1 -k2,2n -S 4G -T "$WORK" "$TA") -g "$CHR" -bg \
          > "$WORK/cov.bg"
        bg_to_bw "$WORK/cov.bg" "$out_cov"
    fi
    out_5p="$OUT/${SAFE}_atac_5p.bw"
    if [[ -s "$out_5p" ]]; then echo "  exists: $(basename "$out_5p")"; else
        bedtools genomecov -i <(sort -k1,1 -k2,2n -S 4G -T "$WORK" "$TA") -g "$CHR" -bg -5 \
          > "$WORK/5p.bg"
        bg_to_bw "$WORK/5p.bg" "$out_5p"
    fi
fi

# ---------------------------------------------------------------- peaks
PK_SRC=$(ls "$DATA/$CT/EP300"/*.bed.gz 2>/dev/null | head -1 || true)
if [[ -n "$PK_SRC" ]]; then
    zcat "$PK_SRC" | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
      | sort -k1,1 -k2,2n > "$OUT/${SAFE}_ep300_peaks.narrowPeak"
    echo "  peaks: $(wc -l < "$OUT/${SAFE}_ep300_peaks.narrowPeak") from $(basename "$PK_SRC")"
fi

echo "=== done $CT ==="
ls -lh "$OUT" | grep -i "$SAFE" || true
