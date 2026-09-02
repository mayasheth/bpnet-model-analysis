#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/tn5cal.%j.txt
#SBATCH -e log/tn5cal.%j.txt
#SBATCH --job-name=tn5cal
#
# Calibrate the Tn5 shift needed to reproduce atac_5p.bw from the paired-end BAMs.
#
# WHY. atac_5p.bw was built from *.tn5.sorted.tagAlign.gz, which is already shifted, so
# 0.20 applied `genomecov -bg -5` with no further offset. The fragment-size channels must
# be built from the PE BAMs instead (TLEN exists only there), and those are UNSHIFTED
# alignments. If the offset is guessed wrong, the size-stratified channels sit a few bp
# off the flat channel and the model sees a phase artefact rather than nucleosome
# positioning.
#
# METHOD. On chr20, dump unshifted 5' ends separately by strand, then score every
# (plus_delta, minus_delta) pair in -8..+8 by Pearson r against atac_5p.bw. The argmax is
# the shift the tagAlign pipeline applied. Reported alongside the r at (0,0) so it is
# obvious whether the peak is real.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
PY="$D/.pixi/envs/multimodal/bin/python"
BAM_DIR="${SCRATCH}/atac_pe"
WORK="${SCRATCH}/tn5cal_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT
mkdir -p "$P/log" "$P/results"

CHROM=chr20
for b in ENCFF077FBI ENCFF128WZG ENCFF534DCE; do
    samtools view -b -f 0x2 -F 0x400 "$BAM_DIR/$b.pe.bam" "$CHROM" \
      | bedtools bamtobed -i stdin
done | awk -v OFS='\t' '{if($6=="+") print "P", $2; else print "M", $3-1}' > "$WORK/ends.txt"
wc -l "$WORK/ends.txt"

$PY "$P/scripts/0.22.calibrate_tn5_shift.py" "$WORK/ends.txt" "$CHROM" \
  | tee "$P/results/tn5_shift_calibration.txt"
