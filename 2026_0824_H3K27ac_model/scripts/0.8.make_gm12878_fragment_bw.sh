#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/gm_frag_bw.%j.txt
#SBATCH -e log/gm_frag_bw.%j.txt
#SBATCH --job-name=gm_frag_bw
#
# GM12878 H3K27ac 250bp-fragment-extended BigWig.
#
# Needed because the current K562 models were TRAINED on the fragment-extended target.
# Transfer must be evaluated against a GM12878 target processed the same way, or the
# transfer number confounds cell-type difference with target-definition difference.
# Once the K562 models are retrained on 5' ends, the 5' GM12878 tracks (already built by
# 0.7) become the matching pair and this track is only needed for the current comparison.
#
# Matches $OAK/Users/sheth/Data/scripts/bam_to_bigWig.sh -r SINGLE: genomecov -bg -fs 250.

set -euo pipefail
export PYTHONUNBUFFERED=1
PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
export PATH="$PROJECT_DIR/.pixi/envs/multimodal/bin:$PATH"
GM_DIR="/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878"
CHR_SIZES="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes"
OUT="$PROJ/data/gm12878_h3k27ac_frag250.bw"
mkdir -p "$PROJ/data" "$PROJ/log"
WORK="${SCRATCH:-/tmp}/gm_frag_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT

if [[ -s "$OUT" ]]; then echo "$OUT exists"; exit 0; fi
samtools merge -@ 4 -f "$WORK/m.bam" "$GM_DIR/ENCFF645BAL.filtered.sorted.bam" \
                                     "$GM_DIR/ENCFF865OOP.filtered.sorted.bam"
samtools index "$WORK/m.bam"
bedtools genomecov -ibam "$WORK/m.bam" -bg -fs 250 \
  | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR_SIZES" - \
  | sort -k1,1 -k2,2n -T "$WORK" > "$WORK/m.bedGraph"
bedGraphToBigWig "$WORK/m.bedGraph" "$CHR_SIZES" "$OUT"
echo "wrote $OUT"; ls -la "$OUT"
