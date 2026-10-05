#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 16:00:00
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -o log/tobias.%j.txt
#SBATCH --job-name=tobias
# 2.2 TOBIAS (Bentsen et al. 2020) footprints on the Xu K562 pseudobulk: ATACorrect -> ScoreBigwig
# -> BINDetect over all JASPAR 2024 non-redundant motifs (the 4.44 set), restricted to scE2G's
# candidate elements so scores join the feature table directly. Insertion BAM from 2.1 is already
# Tn5-shifted, hence --read_shift 0 0. Run from the A2 dir: sbatch scripts/2.2.submit_tobias.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
T=$A2/.pixi/envs/tobias/bin
PEAKS=${PEAKS:-$SCRATCH/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed}
BAM=$SCRATCH/tobias_k562_xu/k562_xu_insertions.bam
GENOME=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
BLACKLIST=/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G_main/ENCODE_rE2G/ABC/reference/hg38/GRCh38_unified_blacklist.bed
MOTIFS=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
OUT=$SCRATCH/tobias_k562_xu
for f in "$PEAKS" "$BAM" "$BAM.bai" "$GENOME" "$BLACKLIST" "$MOTIFS"; do
    [[ -s "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done
echo "peaks: $PEAKS ($(wc -l < "$PEAKS") regions)"
# keep only the chromosomes present in the BAM; TOBIAS errors on contigs missing from the header
$T/samtools idxstats $BAM | awk '$3 > 0 {print $1}' > $OUT/bam_chroms.txt
awk 'NR==FNR {k[$1]; next} ($1 in k)' $OUT/bam_chroms.txt "$PEAKS" | cut -f1-3 | sort -k1,1 -k2,2n > $OUT/peaks.bed
echo "peaks on BAM chroms: $(wc -l < $OUT/peaks.bed)"

$T/TOBIAS ATACorrect --bam $BAM --genome $GENOME --peaks $OUT/peaks.bed --blacklist $BLACKLIST \
    --read_shift 0 0 --prefix k562_xu --outdir $OUT/atacorrect --cores $SLURM_CPUS_PER_TASK
$T/TOBIAS ScoreBigwig --signal $OUT/atacorrect/k562_xu_corrected.bw --regions $OUT/peaks.bed \
    --output $OUT/k562_xu_footprints.bw --cores $SLURM_CPUS_PER_TASK
$T/TOBIAS BINDetect --motifs $MOTIFS --signals $OUT/k562_xu_footprints.bw --genome $GENOME \
    --peaks $OUT/peaks.bed --cond-names K562_Xu --skip-excel --outdir $OUT/bindetect \
    --cores $SLURM_CPUS_PER_TASK
date > $OUT/tobias_complete.txt
echo "tobias complete"
