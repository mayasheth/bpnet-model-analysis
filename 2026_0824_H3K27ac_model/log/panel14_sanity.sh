#!/bin/bash
# Sanity check for the panel_14 biosamples, with the three 2026-09-29 tissues as reference rows:
# p300 peak count, fraction of p300 peaks overlapping the biosample's own ATAC IDR peaks
# (bedtools intersect -u, any overlap), total ATAC alignments from results/panel_depth_audit.tsv
# (0.44's unit), and that every track 1.33 points at opens and has data. Flags < 3,000 p300
# peaks or < 20M ATAC reads. Ends with TRACKS_VERIFIED only when every track passes.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
cd "$P"
PD=data/p300_panel
AUD=results/panel_depth_audit.tsv
NEW=SK-N-SH_RA,esophagus_muscularis_mucosa,gastroesophageal_sphincter,sigmoid_colon,tibial_nerve,upper_lobe_of_left_lung
OLD=stomach,transverse_colon,ovary
bad=0
printf "%-28s %-4s %8s %9s %14s %s\n" cell set p300_pk frac_atac atac_reads flags
for spec in "new:$NEW:$SCRATCH/p300_panel14" "ref:$OLD:$SCRATCH/p300_tissues"; do
    IFS=: read -r tag cells root <<< "$spec"
    IFS=, read -ra CL <<< "$cells"
    for c in "${CL[@]}"; do
        s=${c//-/_}
        pk="$PD/${s}_ep300_peaks.narrowPeak"
        atpk=$(ls "$root/$c/ATAC"/*.bed.gz)
        n=$(wc -l < "$pk")
        ov=$(cut -f1-3 "$pk" | bedtools intersect -u -a - -b <(zcat "$atpk" | cut -f1-3) | wc -l)
        frac=$(awk -v a="$ov" -v b="$n" 'BEGIN{printf "%.3f", (b>0)?a/b:0}')
        reads=$(awk -F'\t' -v c="$c" '$1==c && $2=="ATAC"{s+=$5} END{print s+0}' "$AUD")
        fl=""
        [[ "$n" -lt 3000 ]] && fl="$fl LOW_PEAKS"
        if [[ "$reads" -eq 0 ]]; then fl="$fl NO_AUDIT_ROW"
        elif [[ "$reads" -lt 20000000 ]]; then fl="$fl LOW_ATAC"; fi
        if [[ "$tag" == new ]]; then
            for t in ep300_5p_plus ep300_5p_minus atac atac_5p; do
                bw="$PD/${s}_${t}.bw"
                if ! python3 -c 'import sys, pyBigWig
b = pyBigWig.open(sys.argv[1]); h = b.header()
assert h["nBasesCovered"] > 0 and h["sumData"] > 0, h
assert "chr1" in b.chroms()' "$bw" 2>/dev/null; then
                    fl="$fl BAD_TRACK:$t"; bad=$((bad + 1))
                fi
            done
        fi
        printf "%-28s %-4s %8d %9s %14d %s\n" "$c" "$tag" "$n" "$frac" "$reads" "${fl:- ok}"
    done
done
[[ $bad -eq 0 ]] && echo TRACKS_VERIFIED || { echo "TRACKS_FAILED: $bad"; exit 1; }
