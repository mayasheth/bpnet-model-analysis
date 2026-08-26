#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 6:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/gm_5p_bw.%j.txt
#SBATCH -e log/gm_5p_bw.%j.txt
#SBATCH --job-name=gm_5p_bw
#
# GM12878 H3K27ac 5'-end BigWigs, for cross-cell-type transfer.
# Same processing as the K562 target (0.5.make_5prime_bigwigs.sh) so transfer numbers
# are not confounded by a processing difference between train and test cell type.
#
# GM12878 H3K27ac: ENCFF645BAL, ENCFF865OOP (single-end), per
# $OAK/Users/sheth/Data/ENCODE/GM12878/log.sh.
#
# Usage: sbatch 0.7.make_gm12878_5prime.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
BIN="$PROJECT_DIR/.pixi/envs/multimodal/bin"
export PATH="$BIN:$PATH"

GM_DIR="/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
REP1="$GM_DIR/ENCFF645BAL.filtered.sorted.bam"
REP2="$GM_DIR/ENCFF865OOP.filtered.sorted.bam"

OUT_DIR="$PROJ/data"
mkdir -p "$OUT_DIR" "$PROJ/log"
WORK="${SCRATCH:-/tmp}/gm_5p_$$"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT

for f in "$REP1" "$REP2"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

make_5p () {
    local bam="$1" prefix="$2"
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT_DIR/${prefix}_5p_${tag}.bw"
        if [[ -s "$out" ]]; then echo "  $out exists, skipping"; continue; fi
        local bg="$WORK/${prefix}_${tag}.bedGraph"
        bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
            | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' \
                  OFS='\t' "$CHR_SIZES" - \
            | sort -k1,1 -k2,2n -T "$WORK" > "$bg"
        bedGraphToBigWig "$bg" "$CHR_SIZES" "$out"
        rm -f "$bg"
        echo "  wrote $out"
    done
}

echo "=== merged (evaluation target) ==="
MERGED="$WORK/gm_merged.bam"
samtools merge -@ 4 -f "$MERGED" "$REP1" "$REP2"
samtools index "$MERGED"
make_5p "$MERGED" "gm12878_h3k27ac"
rm -f "$MERGED" "$MERGED.bai"

echo "=== per-replicate (GM12878 ceiling) ==="
make_5p "$REP1" "gm12878_h3k27ac_rep1"
make_5p "$REP2" "gm12878_h3k27ac_rep2"

ls -la "$OUT_DIR"/gm12878_*_5p_*.bw
