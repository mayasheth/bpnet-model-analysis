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
# Cell type selects the BAMs and the reference bigwig. Default k562 so prior runs are
# reproduced byte for byte.
CELL=${1:-k562}
case "$CELL" in
  k562)
    BAMS=("${SCRATCH}/atac_pe/ENCFF077FBI.pe.bam" "${SCRATCH}/atac_pe/ENCFF128WZG.pe.bam"
          "${SCRATCH}/atac_pe/ENCFF534DCE.pe.bam")
    REF_BW="$D/2026_0529_multimodal_p300_model/data/atac_5p.bw" ;;
  gm12878)
    G=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878/ATAC
    BAMS=("$G/ENCFF440GRZ.bam" "$G/ENCFF962FMH.bam" "$G/ENCFF981FXV.bam")
    REF_BW="$D/2026_0606_GM12878_transferability/data/atac_5p.bw" ;;
  *) echo "unknown CELL '$CELL'" >&2; exit 1 ;;
esac
for b in "${BAMS[@]}"; do [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }; done
[[ -s "$REF_BW" ]] || { echo "ERROR missing $REF_BW" >&2; exit 1; }

# This script reads ONE chromosome, which needs an index. The K562 BAMs on $SCRATCH were
# indexed when they were downloaded; the GM12878 ones on Oak were not, and under
# `set -euo pipefail` the missing index kills the run five seconds in with samtools' error
# swallowed by the 2>/dev/null on the region query. Build the index if it is absent.
for b in "${BAMS[@]}"; do
    if [[ ! -s "$b.bai" && ! -s "${b%.bam}.bai" ]]; then
        echo "indexing $(basename "$b") (no .bai found)"
        samtools index -@ 4 "$b"
    fi
done

WORK="${SCRATCH}/tn5cal_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT
mkdir -p "$P/log" "$P/results"

CHROM=chr20
echo "cell=$CELL ref=$REF_BW"
for b in "${BAMS[@]}"; do
    samtools view -b -f 0x2 -F 0x400 "$b" "$CHROM" \
      | bedtools bamtobed -i stdin
done | awk -v OFS='\t' '{if($6=="+") print "P", $2; else print "M", $3-1}' > "$WORK/ends.txt"
wc -l "$WORK/ends.txt"

$PY "$P/scripts/0.22.calibrate_tn5_shift.py" "$WORK/ends.txt" "$CHROM" "$REF_BW" \
  | tee "$P/results/tn5_shift_calibration_${CELL}.txt"
