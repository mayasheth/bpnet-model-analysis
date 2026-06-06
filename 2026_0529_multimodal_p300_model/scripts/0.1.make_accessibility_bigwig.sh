#!/bin/bash
# Convert ATAC tagAlign or DNase BAM files to 1-bp resolution BigWig.
#
# ATAC mode:  accepts one or more .tagAlign.gz files (already Tn5-shift
#             corrected); merges them and runs bedtools genomecov.
# DNase mode: accepts one or more indexed BAM files; merges with samtools,
#             then runs bedtools genomecov (raw counts; training normalizes).
#
# Usage:
#   bash 0.1.make_accessibility_bigwig.sh \
#       --input FILE1 [FILE2 ...] \
#       --output OUTPUT.bw \
#       --chrom-sizes CHROM_SIZES \
#       --type {atac|dnase}
#
# Requirements: bedtools, samtools, bedGraphToBigWig (pixi multimodal env)

set -euo pipefail

INPUTS=()
OUTPUT=""
CHROM_SIZES=""
TYPE="atac"

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                INPUTS+=("$1")
                shift
            done
            ;;
        --output)       OUTPUT=$2;       shift 2 ;;
        --chrom-sizes)  CHROM_SIZES=$2;  shift 2 ;;
        --type)         TYPE=$2;         shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ ${#INPUTS[@]} -eq 0 || -z "$OUTPUT" || -z "$CHROM_SIZES" ]] && {
    echo "Error: --input, --output, and --chrom-sizes are required"; exit 1
}

mkdir -p "$(dirname "$OUTPUT")"
WORK_TMP="/tmp/bw_$$"
mkdir -p "$WORK_TMP"
trap "rm -rf $WORK_TMP" EXIT

if [[ "$TYPE" == "atac" ]]; then
    echo "Generating ATAC BigWig from tagAlign files (merged, 1bp resolution)..."
    echo "  Inputs: ${INPUTS[*]}"

    MERGED_BG="$WORK_TMP/merged.bedGraph"

    # Merge tagAlign files, keep only canonical chroms, compute per-base coverage
    # Use local /tmp for sort temp files to avoid slow network filesystem I/O
    zcat "${INPUTS[@]}" \
        | awk 'NR==FNR{chroms[$1]=1; next} ($1 in chroms)' \
            "$CHROM_SIZES" - \
        | sort -k1,1 -k2,2n -T "$WORK_TMP" \
        | bedtools genomecov -i stdin -g "$CHROM_SIZES" -bg \
        | sort -k1,1 -k2,2n -T "$WORK_TMP" \
        > "$MERGED_BG"

    bedGraphToBigWig "$MERGED_BG" "$CHROM_SIZES" "$OUTPUT"

elif [[ "$TYPE" == "dnase" ]]; then
    echo "Generating DNase BigWig from BAM files (raw counts, 1bp resolution)..."
    echo "  Inputs: ${INPUTS[*]}"

    MERGED_BG="$WORK_TMP/merged.bedGraph"

    if [[ ${#INPUTS[@]} -gt 1 ]]; then
        MERGED_BAM="$WORK_TMP/merged.bam"
        samtools merge -f "$MERGED_BAM" "${INPUTS[@]}"
        samtools index "$MERGED_BAM"
        INPUT_BAM="$MERGED_BAM"
    else
        INPUT_BAM="${INPUTS[0]}"
    fi

    bedtools genomecov -ibam "$INPUT_BAM" -bg \
        | awk 'NR==FNR{chroms[$1]=1; next} ($1 in chroms)' "$CHROM_SIZES" - \
        | sort -k1,1 -k2,2n -T "$WORK_TMP" \
        > "$MERGED_BG"

    bedGraphToBigWig "$MERGED_BG" "$CHROM_SIZES" "$OUTPUT"

else
    echo "Error: --type must be atac or dnase"; exit 1
fi

echo "Done: $OUTPUT"
