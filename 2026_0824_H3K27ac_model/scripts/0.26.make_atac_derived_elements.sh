#!/bin/bash
#
# Convert the ATAC-derived candidate regions to the narrowPeak layout load_peaks expects.
#
# Source is the rE2G ATAC_H3K27ac_powerlaw run, which is byte-identical (md5
# 7d5995ce17fbaad18958f19d5f0b6e1b) to the candidate regions of the July 2026 ABC run, so a
# model trained on this file is trained on exactly the regions ABC scores.
#
# load_peaks reads 10 columns and centres each window on start + summit. The existing
# DNase element file already uses the region MIDPOINT as its summit (e.g. 500 bp region,
# summit 250), so setting summit to the midpoint here makes the two files differ only in
# which assay called the regions -- no change in centring convention.
#
# Usage: bash 0.26.make_atac_derived_elements.sh
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
SRC=/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2025_0226_ATAC_powerlaw_models/ATAC_H3K27ac_powerlaw/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
OUT=$D/reference/K562_ATAC_candidate_elements.narrowPeak

[[ -s "$SRC" ]] || { echo "ERROR missing $SRC" >&2; exit 1; }
echo "source md5: $(md5sum "$SRC" | cut -d' ' -f1)  (expect 7d5995ce17fbaad18958f19d5f0b6e1b)"

awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
        mid = int(($3 - $2) / 2)
        print $1, $2, $3, $1":"$2"-"$3, 0, ".", 0, -1, -1, mid
     }' "$SRC" > "$OUT"

echo "wrote $OUT"
wc -l "$OUT" "$D/reference/K562_DNase_candidate_elements.narrowPeak"
echo "--- format check against the DNase file ---"
head -2 "$OUT"
awk '{if (NF != 10) {print "ERROR: line "NR" has "NF" fields" > "/dev/stderr"; exit 1}}' "$OUT"
echo "all rows have 10 fields"
