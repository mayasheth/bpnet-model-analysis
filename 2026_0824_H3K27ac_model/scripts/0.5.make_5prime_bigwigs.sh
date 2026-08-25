#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 8:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/5prime_bw.%j.txt
#SBATCH -e log/5prime_bw.%j.txt
#SBATCH --job-name=k27_5p_bw
#
# Build 5'-END H3K27ac BigWigs, matching the p300 convention.
#
# Why: the existing target (Data/share/IGV/ENCSR000AKP_coverage.bw) extends 36 bp reads
# to a fixed 250 bp fragment. That was written for IGV display and was never tuned for
# modelling, and it has two costs:
#   1. It smears signal across the nucleosome-free center, flattening the very
#      bimodality that motivates a wide counting window (observed dip is only ~9% below
#      the shoulders, which is far shallower than the biology should give).
#   2. It makes the profile task ill-posed. bpnetlite's profile loss is MNLL, which
#      expects multinomial READ COUNTS. Fragment-extended coverage inflates per-window
#      totals ~250x and makes adjacent positions near-identical, which is why profile
#      MNLL sat at ~2800 here versus ~500 for the p300 runs.
# 5'-end counting fixes both: one count per read at a single base, which is what MNLL
# is built for, and no bleed across the NFR or into neighbouring elements.
#
# Matches scripts/0.3.make_training_bw.sh (the p300 path): genomecov -5 -dz, stranded.
# Stranded costs nothing here — bpnetlite's count target sums all channels, so counts are
# identical either way, but the profile head gets a more meaningful target and the setup
# mirrors p300 exactly for comparison.
#
# Builds merged (training target) and per-replicate (inter-replicate ceiling) tracks.
#
# Usage: sbatch 0.5.make_5prime_bigwigs.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
BIN="$PROJECT_DIR/.pixi/envs/multimodal/bin"
export PATH="$BIN:$PATH"

ENCODE_DIR="/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
REP1="$ENCODE_DIR/ENCFF790GFL.se.filtered.sorted.bam"
REP2="$ENCODE_DIR/ENCFF817HMW.se.filtered.sorted.bam"

OUT_DIR="$PROJ/data"
mkdir -p "$OUT_DIR" "$PROJ/log"

WORK="${SCRATCH:-/tmp}/k27_5p_$$"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT

# Emit stranded 5'-end BigWigs for one BAM (or a merged BAM).
make_5p () {
    local bam="$1" prefix="$2"
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT_DIR/${prefix}_5p_${tag}.bw"
        if [[ -s "$out" ]]; then echo "  $out exists, skipping"; continue; fi
        local bg="$WORK/${prefix}_${tag}.bedGraph"
        # -5: count only the 5' end of each read.  -dz: per-base, 0-based, non-zero only.
        bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
            | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' \
                  OFS='\t' "$CHR_SIZES" - \
            | sort -k1,1 -k2,2n -T "$WORK" \
            > "$bg"
        # NOTE: 0.3.make_training_bw.sh does NOT sort its experiment bedGraphs.
        # bedGraphToBigWig requires chrom+start order, and genomecov follows BAM header
        # order, which need not match. Sorting explicitly here.
        bedGraphToBigWig "$bg" "$CHR_SIZES" "$out"
        rm -f "$bg"
        echo "  wrote $out"
    done
}

echo "=== merged (training target) ==="
MERGED="$WORK/merged.bam"
samtools merge -@ 4 -f "$MERGED" "$REP1" "$REP2"
samtools index "$MERGED"
make_5p "$MERGED" "h3k27ac"
rm -f "$MERGED" "$MERGED.bai"

echo "=== per-replicate (for the ceiling) ==="
make_5p "$REP1" "h3k27ac_rep1"
make_5p "$REP2" "h3k27ac_rep2"

echo "=== done ==="
ls -la "$OUT_DIR"/*_5p_*.bw
