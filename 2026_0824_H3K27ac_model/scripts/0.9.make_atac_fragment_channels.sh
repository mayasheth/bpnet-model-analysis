#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/atac_frag.%j.txt
#SBATCH -e log/atac_frag.%j.txt
#SBATCH --job-name=atac_frag
#
# Fragment-size-stratified ATAC BigWigs — supplies nucleosome POSITIONING, which flat
# ATAC coverage discards. H3K27ac can only exist where a nucleosome exists, so this is
# the mechanistically missing input.
#
# Two channels, the standard ATAC split:
#   sub      < 100 bp   nucleosome-free / TF-bound
#   mono   180-247 bp   mono-nucleosomal
#
# Built from the paired-end BAMs on $SCRATCH (fragment length = TLEN, only available in
# the PE BAMs; the tagAlign files are per-read and lose it). Oak is at 98%, so the BAMs
# stay on SCRATCH and only the small BigWigs land on Oak.
#
# MultiModalBPNet already accepts arbitrary accessibility channels via n_acc_filters, so
# no architecture change is needed — only more input channels.

set -euo pipefail
export PYTHONUNBUFFERED=1
PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
export PATH="$PROJECT_DIR/.pixi/envs/multimodal/bin:$PATH"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
BAM_DIR="${SCRATCH}/atac_pe"
mkdir -p "$PROJ/data" "$PROJ/log"
WORK="${SCRATCH}/atac_frag_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT

BAMS=("$BAM_DIR/ENCFF077FBI.pe.bam" "$BAM_DIR/ENCFF128WZG.pe.bam" "$BAM_DIR/ENCFF534DCE.pe.bam")
for b in "${BAMS[@]}"; do [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }; done

# fragment-length histogram, for the record — the split points should be checked against it
echo "=== fragment length distribution (chr1 sample, rep1) ==="
samtools view -f 0x2 "${BAMS[0]}" chr1 2>/dev/null | head -500000 \
  | awk '{if($9>0 && $9<1000) print int($9/10)*10}' | sort -n | uniq -c \
  | awk '{printf "%5d bp: %d\n", $2, $1}' | head -40 \
  > "$PROJ/results/atac_fragment_length_hist.txt" || true
head -12 "$PROJ/results/atac_fragment_length_hist.txt" || true

make_channel () {
    local name="$1" lo="$2" hi="$3"
    local out="$PROJ/data/atac_${name}.bw"
    if [[ -s "$out" ]]; then echo "  $out exists"; return; fi
    echo "=== channel $name (${lo}-${hi} bp) ==="
    # Properly-paired reads, TLEN in range; emit the FRAGMENT interval once per pair
    # (positive TLEN only, so each pair is counted a single time).
    for b in "${BAMS[@]}"; do
        samtools view -@ 2 -f 0x2 "$b" \
          | awk -v lo="$lo" -v hi="$hi" 'BEGIN{OFS="\t"}
              $9 > 0 && $9 >= lo && $9 <= hi { print $3, $4-1, $4-1+$9 }'
    done \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR_SIZES" - \
      | sort -k1,1 -k2,2n -T "$WORK" \
      | bedtools genomecov -i stdin -g "$CHR_SIZES" -bg \
      | sort -k1,1 -k2,2n -T "$WORK" > "$WORK/${name}.bedGraph"
    bedGraphToBigWig "$WORK/${name}.bedGraph" "$CHR_SIZES" "$out"
    rm -f "$WORK/${name}.bedGraph"
    echo "  wrote $out"
}

make_channel sub  1   99
make_channel mono 180 247

ls -la "$PROJ/data"/atac_sub.bw "$PROJ/data"/atac_mono.bw
