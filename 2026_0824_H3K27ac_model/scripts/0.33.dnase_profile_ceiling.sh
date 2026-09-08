#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/dnaseceil.%j.txt
#SBATCH -e log/dnaseceil.%j.txt
#SBATCH --job-name=dnaseceil
#
# Does DNase have a reproducible base-resolution profile? This gates the ATAC->DNase converter.
#
# The converter has to emit something the downstream model can eat, and that model's
# accessibility input is a base-resolution 5' insertion track read over a 2,114 bp window. A
# converter predicting only per-element DNase COUNTS and painting them flat would strip the
# base-resolution structure the accessibility branch uses, so it has to predict the PROFILE
# and the painted track has to be profile x counts at base resolution.
#
# That inverts which head matters. Every result in this project so far comes from the counts
# head, and the profile head has been treated as near-useless because H3K27ac has no
# reproducible base-resolution shape (inter-replicate ceiling 0.21 at 1 bp, rising to 0.72 at
# 50 bp). DNase should differ: the cut sites ARE the signal, rather than a smear over the
# nucleosomes flanking the element. If DNase's 1 bp ceiling is also ~0.2, the converter cannot
# produce a usable input track and the idea dies here for the cost of this job.
#
# RUN-TYPE CAVEAT. K562's ENCSR000EOT is listed as mixed pe/se, and a recorded constraint says
# HCT116's H3K27ac was excluded from a ceiling for exactly that reason. The objection is
# weaker under 5' counting, which takes one cut site per read regardless of length and uses
# read 1 only for paired data, so a PE and an SE library are comparable at the cut-site level.
# GM12878's ENCSR000EMT is two clean SE replicates and is the control: if the two cell types
# give the same shape of answer, the mixed run type is not distorting K562's.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
PY=$D/.pixi/envs/multimodal/bin/python
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT=$P/data
mkdir -p "$OUT" "$P/log" "$P/results"
W=${SCRATCH}/dnaseceil_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

# Stranded per-replicate 5' tracks, read 1 only where the BAM is paired.
make_rep () {
    local bam=$1 prefix=$2
    local n_paired filt
    n_paired=$(samtools view -c -f 0x1 -@ 8 "$bam")
    if [[ "$n_paired" -gt 0 ]]; then filt="-F 0x900 -f 0x40"; else filt="-F 0x900"; fi
    printf '  %-30s paired_reads=%d  filter="%s"\n' "$(basename "$bam")" "$n_paired" "$filt"
    samtools view -b $filt -@ 8 "$bam" > "$W/f.bam"
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_${tag}.bw"
        [[ -s "$out" ]] && { echo "    exists $(basename "$out")"; continue; }
        bedtools genomecov -ibam "$W/f.bam" -5 -dz -strand "$strand" \
            | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
            | LC_COLLATE=C sort -k1,1 -k2,2n -T "$W" > "$W/x.bg"
        bedGraphToBigWig "$W/x.bg" "$CHR" "$out"; rm -f "$W/x.bg"
        echo "    wrote $(basename "$out")"
    done
    rm -f "$W/f.bam"
}

echo "=== K562 DNase replicates (ENCSR000EOT, listed pe+se) ==="
make_rep "$E/K562/ENCFF205FNC.filtered.sorted.bam"  k562_dnase_rep1_5p
make_rep "$E/K562/ENCFF860XAE.filtered.sorted.bam"  k562_dnase_rep2_5p
echo "=== GM12878 DNase replicates (ENCSR000EMT, se) ==="
make_rep "$E/GM12878/ENCFF467CXY_sorted.bam" gm12878_dnase_rep1_5p
make_rep "$E/GM12878/ENCFF940NSD_sorted.bam" gm12878_dnase_rep2_5p

echo
echo "########## K562 DNase profile ceiling ##########"
$PY scripts/0.25.profile_ceiling_by_binsize.py --label k562_dnase \
    --rep1-plus data/k562_dnase_rep1_5p_plus.bw --rep1-minus data/k562_dnase_rep1_5p_minus.bw \
    --rep2-plus data/k562_dnase_rep2_5p_plus.bw --rep2-minus data/k562_dnase_rep2_5p_minus.bw \
    --elements $D/reference/K562_DNase_candidate_elements.narrowPeak
echo
echo "########## GM12878 DNase profile ceiling (clean SE control) ##########"
$PY scripts/0.25.profile_ceiling_by_binsize.py --label gm12878_dnase \
    --rep1-plus data/gm12878_dnase_rep1_5p_plus.bw \
    --rep1-minus data/gm12878_dnase_rep1_5p_minus.bw \
    --rep2-plus data/gm12878_dnase_rep2_5p_plus.bw \
    --rep2-minus data/gm12878_dnase_rep2_5p_minus.bw \
    --elements $D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
echo
echo "Compare against H3K27ac, whose 1 bp top-quintile ceiling is 0.21 (K562) and 0.18"
echo "(GM12878), rising to 0.72 and 0.70 at 50 bp. A DNase 1 bp ceiling near those numbers"
echo "means the converter cannot emit a usable base-resolution track."
