#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 12:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/telohaec_tracks.%j.txt
#SBATCH -e log/telohaec_tracks.%j.txt
#SBATCH --job-name=k27_telo_tracks
#
# Build TeloHAEC ATAC + H3K27ac tracks for the four conditions.
#
# EXPLICIT FILE LISTS, NEVER A GLOB. Three separate traps in this directory:
#   1. TeloHAEC_ctrl/ATAC holds 9 bams/tagAligns of which 3 are cell_line = Eahy926
#      (SRR20809434/435/436) -- a different endothelial line under the same sample_name.
#   2. There are duplicate `initial_*.tn5.sorted.tagAlign.gz` files alongside the real ones.
#   3. SRR20809416.filtered.sorted.dedup.tagAlign.gz is 20 bytes -- a truncated artifact.
#
# ATAC is built from Tn5-shift-corrected tagAlign via 0.1.make_accessibility_bigwig.sh,
# matching exactly how K562's data/atac.bw was built. Using BAMs instead would make
# TeloHAEC's input differ from K562's in construction as well as cell type.
#
# H3K27ac 5' ends: TeloHAEC's ChIP is PAIRED-end, while K562 and GM12878 are SINGLE-end.
# A naive `genomecov -5` on PE data counts TWO 5' ends per fragment, one at each end, which
# is not the same quantity the SE tracks measure. Rather than assume, build BOTH and compare
# (Maya, 2026-08-29) -- the same way the 5'-vs-fragment-extension question was settled:
#   *_r1_5p_*    R1 only (-f 64): one 5' end per fragment, matches the SE convention.
#   *_both_5p_*  both mates: all reads, but signal at both fragment ends.
# Per-replicate tracks are built for each variant so the inter-replicate ceiling can be
# compared too -- a variant that raises the ceiling is measuring a more reproducible
# quantity, which is the criterion that settled the 5' switch.
#
# Usage: sbatch 0.15.build_telohaec_tracks.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
T=/oak/stanford/groups/engreitz/Projects/E2G/endothelial_cells/TeloHAEC_GSE210489_GSE210491
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT="$P/data/panel"; mkdir -p "$OUT" "$P/log"
WORK="${SCRATCH}/telo_tracks_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT

# --- TeloHAEC-only accessions, Eahy926 deliberately excluded -----------------
ATAC_ctrl="SRR20809416 SRR20809419 SRR20809420 SRR20809430 SRR20809431 SRR20809433"
ATAC_IL1b="SRR20809411 SRR20809414 SRR20809415"
ATAC_TNFa="SRR20809409 SRR20809410 SRR20809412"
ATAC_VEGF="SRR20809407 SRR20809408 SRR20809413"

K27_ctrl="SRR20810532 SRR20810533 SRR20810544 SRR20810545"
K27_IL1b="SRR20810530 SRR20810531"
K27_TNFa="SRR20810528 SRR20810529"
K27_VEGF="SRR20810526 SRR20810527"

declare -A DIR=( [ctrl]=TeloHAEC_ctrl [IL1b]=TeloHAEC_IL1b
                 [TNFa]=TeloHAEC_TNFa [VEGF]=TeloHAEC_VEGF_ctrl )

build_atac () {
    local cond="$1"; shift
    local srrs="$1"
    local out="$OUT/TeloHAEC_${cond}_atac.bw"
    [[ -s "$out" ]] && { echo "  exists: $(basename $out)"; return; }
    local files=()
    for s in $srrs; do
        f="$T/${DIR[$cond]}/ATAC/${s}.tn5.sorted.tagAlign.gz"
        [[ -s "$f" ]] || { echo "ERROR missing/empty $f" >&2; exit 1; }
        # guard against the 20-byte truncated artifact pattern
        [[ $(stat -c %s "$f") -gt 1000000 ]] || { echo "ERROR suspiciously small $f" >&2; exit 1; }
        files+=("$f")
    done
    echo "  ATAC $cond: ${#files[@]} tagAlign files"
    bash $D/2026_0529_multimodal_p300_model/scripts/0.1.make_accessibility_bigwig.sh \
        --input "${files[@]}" --output "$out" --chrom-sizes "$CHR" --type atac
}

# stranded 5' ends. variant=r1 keeps read 1 only; variant=both keeps both mates.
make_5p () {
    local bam="$1" prefix="$2" variant="$3"
    local flags
    if [[ "$variant" == "r1" ]]; then flags=(-f 64 -F 3852); else flags=(-F 3852); fi
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_${variant}_5p_${tag}.bw"
        [[ -s "$out" ]] && { echo "  exists: $(basename $out)"; continue; }
        samtools view -b "${flags[@]}" "$bam" \
          | bedtools genomecov -ibam stdin -5 -dz -strand "$strand" \
          | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
          | sort -k1,1 -k2,2n -T "$WORK" > "$WORK/x.bg"
        bedGraphToBigWig "$WORK/x.bg" "$CHR" "$out"; rm -f "$WORK/x.bg"
        echo "  wrote $(basename $out)"
    done
}

build_k27 () {
    local cond="$1"; shift
    local srrs="$1"
    local bams=()
    for s in $srrs; do
        f="$T/${DIR[$cond]}/H3K27ac/${s}.filtered.sorted.dedup.bam"
        [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }
        bams+=("$f")
    done
    echo "  H3K27ac $cond: ${#bams[@]} bams"
    local merged
    if [[ ${#bams[@]} -gt 1 ]]; then
        merged="$WORK/k27_${cond}.bam"
        samtools merge -@ 4 -f "$merged" "${bams[@]}"
        samtools index "$merged"
    else
        merged="${bams[0]}"
    fi
    for variant in r1 both; do
        make_5p "$merged" "TeloHAEC_${cond}_h3k27ac" "$variant"
        local i=1
        for b in "${bams[@]}"; do
            make_5p "$b" "TeloHAEC_${cond}_h3k27ac_rep${i}" "$variant"; i=$((i+1))
        done
    done
    rm -f "$merged" "$merged.bai"
}

echo "=== TeloHAEC ATAC ==="
build_atac ctrl "$ATAC_ctrl"
build_atac IL1b "$ATAC_IL1b"
build_atac TNFa "$ATAC_TNFa"
build_atac VEGF "$ATAC_VEGF"

echo "=== TeloHAEC H3K27ac (5' ends: r1 and both variants) ==="
build_k27 ctrl "$K27_ctrl"
build_k27 IL1b "$K27_IL1b"
build_k27 TNFa "$K27_TNFa"
build_k27 VEGF "$K27_VEGF"

echo "ALL_DONE"
