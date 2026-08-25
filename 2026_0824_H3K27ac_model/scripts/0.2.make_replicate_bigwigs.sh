#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 8:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/rep_bigwigs.%j.txt
#SBATCH -e log/rep_bigwigs.%j.txt
#SBATCH --job-name=k27_rep_bw
#
# Build PER-REPLICATE H3K27ac coverage BigWigs.
#
# The merged track (Data/share/IGV/ENCSR000AKP_coverage.bw) is the training target.
# These per-replicate tracks exist only to measure the inter-replicate ceiling — the
# upper bound on any model's counts correlation — as a function of counting window.
#
# Settings match the merged track exactly so the ceiling is comparable to the target:
# raw counts (no normalization), single-end with fragment extension to 250 bp, via
# bedtools genomecov -bg -fs 250. Same as $OAK/Users/sheth/Data/scripts/bam_to_bigWig.sh
# with -r SINGLE; reimplemented here per replicate to keep the fragment size explicit.
#
# Usage: sbatch 0.2.make_replicate_bigwigs.sh

set -euo pipefail

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
BIN="$PROJECT_DIR/.pixi/envs/multimodal/bin"
export PATH="$BIN:$PATH"

ENCODE_DIR="/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
FRAG_SIZE=250   # must match the merged training target

declare -A REPS=(
  [rep1]="$ENCODE_DIR/ENCFF790GFL.se.filtered.sorted.bam"
  [rep2]="$ENCODE_DIR/ENCFF817HMW.se.filtered.sorted.bam"
)

OUT_DIR="$PROJ/data"
mkdir -p "$OUT_DIR" "$PROJ/log"

# Scratch, not Oak: these intermediates are large and transient.
WORK="${SCRATCH:-/tmp}/k27_rep_bw_$$"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT

for rep in rep1 rep2; do
  BAM="${REPS[$rep]}"
  OUT_BW="$OUT_DIR/h3k27ac_${rep}.bw"

  if [[ -s "$OUT_BW" ]]; then
    echo "$OUT_BW already exists — skipping"
    continue
  fi

  echo "=== $rep: $(basename "$BAM") ==="
  BG="$WORK/${rep}.bedGraph"

  # Restrict to canonical chroms, extend to fragment size, per-base coverage
  bedtools genomecov -ibam "$BAM" -bg -fs "$FRAG_SIZE" \
    | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR_SIZES" - \
    | sort -k1,1 -k2,2n -T "$WORK" \
    > "$BG"

  bedGraphToBigWig "$BG" "$CHR_SIZES" "$OUT_BW"
  rm -f "$BG"
  echo "Wrote $OUT_BW"
done

echo "All done."
ls -la "$OUT_DIR"/h3k27ac_rep*.bw
