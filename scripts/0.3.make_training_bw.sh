#!/bin/bash

# Note: input bam files are assumed to be sorted and indexed!
# The script will create the the following files in the output directory: control_minus.bw, control_plus.bw, minus.bw, plus.bw 

# Initialize variables
OUT_DIR=""
CHR_SIZES=""
BAM_FILES=()
CONTROL_BAM_FILES=()

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -i|--input)
      shift
      while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
        BAM_FILES+=("$1")
        shift
      done
      ;;
    -c|--control)
      shift
      while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
        CONTROL_BAM_FILES+=("$1")
        shift
      done
      ;;
    -o|--outdir)
      OUT_DIR="$2"
      shift 2
      ;;
    -g|--chromsizes)
      CHR_SIZES="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# Create output directories
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/interm"

## Merge and re-index BAM files if necessary
# Experiment
if [ ${#BAM_FILES[@]} -gt 1 ]; then
  MERGED_INPUT_UNFILT="$OUT_DIR/interm/merged.unfiltered.bam"
  MERGED_INPUT="$OUT_DIR/interm/merged.bam"
  samtools merge -f "$MERGED_INPUT_UNFILT" "${BAM_FILES[@]}"
  samtools index "$MERGED_INPUT_UNFILT"
  samtools view -b "$MERGED_INPUT_UNFILT" $(cut -f 1 $CHR_SIZES) -o "$MERGED_INPUT"
  samtools index "$MERGED_INPUT"
  rm "$MERGED_INPUT_UNFILT" "$MERGED_INPUT_UNFILT".bai
else
  MERGED_INPUT="$OUT_DIR/interm/merged.bam"
  samtools view -b "${BAM_FILES[0]}" $(cut -f 1 $CHR_SIZES) -o "$MERGED_INPUT"
  samtools index "$MERGED_INPUT"
fi
echo "Merged and indexed experiment BAM files."

# Control
if [ ${#CONTROL_BAM_FILES[@]} -gt 1 ]; then
  MERGED_CONTROL="$OUT_DIR/interm/control_merged.bam"
  samtools merge -f "$MERGED_CONTROL".unfiltered "${CONTROL_BAM_FILES[@]}"
  samtools view -b "$MERGED_CONTROL".unfiltered $(cut -f 1 $CHR_SIZES) -o "$MERGED_CONTROL"
  samtools index "$MERGED_CONTROL"
  rm "$MERGED_CONTROL".unfiltered
else
  MERGED_CONTROL="$OUT_DIR/interm/control_merged.bam"
  samtools view -b "${CONTROL_BAM_FILES[0]}" $(cut -f 1 $CHR_SIZES) -o "$MERGED_CONTROL"
  samtools index "$MERGED_CONTROL"
fi
echo "Merged and indexed control BAM files."

## Create bedgraph files
# NOTE: bedGraphToBigWig requires chrom+start order. genomecov follows BAM
# header order, which is not guaranteed to match chrom.sizes, so sort
# explicitly. The control bedGraphs below already did; the experiment ones
# did not (fixed 2026-08-25).
# Experiment
BG_PLUS="$OUT_DIR/interm/plus.bedGraph"
BG_MINUS="$OUT_DIR/interm/minus.bedGraph"
bedtools genomecov -ibam "$MERGED_INPUT" -5 -dz -strand + | awk '{print $1, $2, $2+1, $3}' OFS='\t' | sort -k1,1 -k2,2n > "$BG_PLUS"
bedtools genomecov -ibam "$MERGED_INPUT" -5 -dz -strand - | awk '{print $1, $2, $2+1, $3}' OFS='\t' | sort -k1,1 -k2,2n > "$BG_MINUS"
echo "Made experiment bedgraph files."

# Control
CONTROL_BG_PLUS="$OUT_DIR/interm/control_plus.bedGraph"
CONTROL_BG_MINUS="$OUT_DIR/interm/control_minus.bedGraph"
bedtools genomecov -ibam "$MERGED_CONTROL" -5 -dz -strand + | awk '{print $1, $2, $2+1, $3}' OFS='\t' | sort -k1,1 -k2,2n > "$CONTROL_BG_PLUS"
bedtools genomecov -ibam "$MERGED_CONTROL" -5 -dz -strand - | awk '{print $1, $2, $2+1, $3}' OFS='\t' | sort -k1,1 -k2,2n > "$CONTROL_BG_MINUS"
echo "Made control bedgraph files."

## Bedgraph to BigWig
# Experiment
BW_PLUS="$OUT_DIR/plus.bw"
BW_MINUS="$OUT_DIR/minus.bw"
bedGraphToBigWig "$BG_PLUS" "$CHR_SIZES" "$BW_PLUS"
bedGraphToBigWig "$BG_MINUS" "$CHR_SIZES" "$BW_MINUS"
echo "Made experiment bw files!"

# Control
CONTROL_BW_PLUS="$OUT_DIR/control_plus.bw"
CONTROL_BW_MINUS="$OUT_DIR/control_minus.bw"
bedGraphToBigWig "$CONTROL_BG_PLUS" "$CHR_SIZES" "$CONTROL_BW_PLUS"
bedGraphToBigWig "$CONTROL_BG_MINUS" "$CHR_SIZES" "$CONTROL_BW_MINUS"
echo "Made control bw files!"

# Remove temporary files
rm -r "$OUT_DIR/interm"
echo "All done."