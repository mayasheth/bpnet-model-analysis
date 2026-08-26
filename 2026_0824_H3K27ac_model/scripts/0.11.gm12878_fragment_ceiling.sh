#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 6:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/gm_frag_ceiling.%j.txt
#SBATCH -e log/gm_frag_ceiling.%j.txt
#SBATCH --job-name=gm_frag_ceil
#
# GM12878 per-replicate FRAGMENT-EXTENDED tracks + the ceiling on them.
#
# Why this exists: the cross-cell-type transfer evaluation (job 40874283) scored K562
# models against the fragment-extended GM12878 target, because that is what those models
# were trained on. But the GM12878 ceiling computed so far (job 40876333) is on the 5'
# target. In K562 the two ceilings were identical to three decimals (0.7601 vs 0.7605
# top-quintile at +/-500), so substituting one for the other was defensible — but it was
# an assumption, not a measurement, in GM12878. This measures it.
#
# Processing matches the merged GM12878 fragment target (0.8) and the K562 target
# exactly: raw counts, genomecov -bg -fs 250.
#
# Usage: sbatch 0.11.gm12878_fragment_ceiling.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
BIN="$PROJECT_DIR/.pixi/envs/multimodal/bin"
export PATH="$BIN:$PATH"

GM_DIR="/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
ELEMENTS="$PROJECT_DIR/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"
FRAG_SIZE=250

declare -A REPS=( [rep1]="$GM_DIR/ENCFF645BAL.filtered.sorted.bam"
                  [rep2]="$GM_DIR/ENCFF865OOP.filtered.sorted.bam" )

mkdir -p "$PROJ/data" "$PROJ/log" "$PROJ/results" "$PROJ/figures"
WORK="${SCRATCH:-/tmp}/gm_frag_ceil_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT

for rep in rep1 rep2; do
  OUT="$PROJ/data/gm12878_h3k27ac_${rep}_frag250.bw"
  if [[ -s "$OUT" ]]; then echo "$OUT exists, skipping"; continue; fi
  echo "=== $rep ==="
  BG="$WORK/${rep}.bedGraph"
  bedtools genomecov -ibam "${REPS[$rep]}" -bg -fs "$FRAG_SIZE" \
    | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR_SIZES" - \
    | sort -k1,1 -k2,2n -T "$WORK" > "$BG"
  bedGraphToBigWig "$BG" "$CHR_SIZES" "$OUT"
  rm -f "$BG"
  echo "wrote $OUT"
done

echo "=== ceiling on the fragment-extended GM12878 target ==="
"$BIN/python" "$PROJ/scripts/0.3.replicate_ceiling_by_window.py" \
  --elements "$ELEMENTS" \
  --rep1-bw "$PROJ/data/gm12878_h3k27ac_rep1_frag250.bw" \
  --rep2-bw "$PROJ/data/gm12878_h3k27ac_rep2_frag250.bw" \
  --outdir "$PROJ/results" --figdir "$PROJ/figures" \
  --label gm12878_frag250_
