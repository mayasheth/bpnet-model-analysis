#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=24G
#SBATCH -c 8
#SBATCH -o log/re2g_gcctcf_feats.%j.txt
#SBATCH --job-name=re2g_gcctcf_feats
# CTCF (JASPAR MA0139.2) sites over the July rE2G candidate regions with TOBIAS TFBScan, the scan
# BINDetect ran for the scE2G features (2.2, defaults p < 1e-4), then 4.62 builds GC + CTCF
# features. Run from 2026_0824_H3K27ac_model: sbatch scripts/4.62.submit_re2g_gc_ctcf_features.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
OAK=/oak/stanford/groups/engreitz/Users/sheth
T=$OAK/EP300_BPNet/2026_1002_sce2g_tf_activity/.pixi/envs/tobias/bin
CAND=$OAK/ABC_working/ABC-Enhancer-Gene-Prediction/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
GENOME=$OAK/hg38_resources/hg38.fa
MOTIFS=$OAK/EP300_BPNet/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
OUT=data/re2g_external
mkdir -p $OUT

# MEME header plus the MA0139.2 block only
awk '/^MOTIF/ {seen = 1; keep = ($2 == "MA0139.2")} !seen || keep' $MOTIFS > $OUT/ctcf_MA0139.2.meme
[[ $(grep -c '^MOTIF' $OUT/ctcf_MA0139.2.meme) == 1 ]] || { echo "motif extraction failed" >&2; exit 1; }
cut -f1-3 $CAND | sort -k1,1 -k2,2n > $OUT/re2g_july_candidates.bed
echo "candidates: $(wc -l < $OUT/re2g_july_candidates.bed)"

$T/TOBIAS TFBScan --motifs $OUT/ctcf_MA0139.2.meme --fasta $GENOME --regions $OUT/re2g_july_candidates.bed \
    --outfile $OUT/re2g_july_candidates.ctcf_MA0139.2.bed --pvalue 1e-4 --cores $SLURM_CPUS_PER_TASK
$T/python scripts/4.62.make_re2g_gc_ctcf_features.py --ctcf-sites $OUT/re2g_july_candidates.ctcf_MA0139.2.bed
