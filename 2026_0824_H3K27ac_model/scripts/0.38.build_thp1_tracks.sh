#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/thp1tracks.%j.txt
#SBATCH -e log/thp1tracks.%j.txt
#SBATCH --job-name=thp1_tracks
#
# THP-1 tracks: the third cell type for the DNase-input panel.
#
# WHAT THP-1 CAN AND CANNOT DO. DNase is ONE replicate, so no DNase inter-replicate ceiling
# exists and THP-1 cannot be a converter TARGET (0.25's shape ceiling needs two). H3K27ac is
# two paired-end replicates of the same run type, so the H3K27ac ceiling IS computable, which
# is the one the panel is read against. No ATAC at all, so THP-1 cannot train a converter
# either. It is a DNase-input panel member and nothing more, and that is still worth having:
# it is the first cell type beyond K562 and GM12878 that can carry a DNase-input model.
#
# READ 1 ONLY for the H3K27ac targets. Both replicates are paired-end, and both mates would
# contribute a 5-prime end and double-count each fragment. This is the project's standing rule
# and it has bitten before.
#
# Pooling is by MERGING BAMs then calling genomecov once, never by combining bedGraphs.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
# macs2 is NOT in the pixi env. The only installation on Oak is inside ABC's built snakemake
# conda env, which is where the project's other MACS2 calls came from. Appended, not
# prepended, so it cannot shadow the pixi python or bedtools.
MACS2_ENV=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/.snakemake/conda/a0febb704d007e3c12f577f28bdeccb8_
export PATH="$PATH:$MACS2_ENV/bin"
MACS2="$MACS2_ENV/bin/python $MACS2_ENV/bin/macs2"  # its own python: PATH alone leaves MACS2 unimportable
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
E2G=/oak/stanford/groups/engreitz/Projects/E2G/THP1/THP1_PRJNA830917
DNASE=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/THP1/DNase/AG81591.filtered.bam
K27_R1=$E2G/THP1_macrophages_0000min_R1/H3K27ac/SRR18899252.pe.filtered.sorted.dedup.bam
K27_R2=$E2G/THP1_macrophages_0000min_R2/H3K27ac/SRR18899253.pe.filtered.sorted.dedup.bam
OUT=$P/data
mkdir -p "$OUT" "$P/log"
W=${SCRATCH}/thp1_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

for f in "$DNASE" "$K27_R1" "$K27_R2"; do
    [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }
done

# Filter to primary reads, and to read 1 where the BAM is paired. Pairedness is detected per
# BAM rather than assumed, the same way 0.32 and 0.34 do it.
filt () {  # in.bam out.bam
    local n_paired
    n_paired=$(samtools view -c -f 0x1 -@ 8 "$1")
    if [[ "$n_paired" -gt 0 ]]; then
        echo "    $(basename "$1"): paired -> read 1 only"
        samtools view -b -F 0x900 -f 0x40 -@ 8 "$1" > "$2"
    else
        echo "    $(basename "$1"): single-end"
        samtools view -b -F 0x900 -@ 8 "$1" > "$2"
    fi
    samtools index -@ 8 "$2"
}

main_chrom_reads () {  # bam -> count on the contigs $CHR keeps
    samtools idxstats "$1" \
      | awk 'NR==FNR {k[$1]=1; next} ($1 in k) {s += $3} END {print s+0}' "$CHR" -
}

# Unstranded 5-prime accessibility input, matching 0.32's convention exactly.
echo "=== THP-1 DNase 5-prime accessibility input ==="
if [[ -s "$OUT/thp1_dnase_5p.bw" ]]; then
    echo "  exists"
else
    filt "$DNASE" "$W/dnase.bam"
    exp=$(main_chrom_reads "$W/dnase.bam")
    bedtools genomecov -ibam "$W/dnase.bam" -bg -5 \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
      | LC_COLLATE=C sort -k1,1 -k2,2n -T "$W" > "$W/d.bg"
    got=$(awk '{s += ($3-$2)*$4} END {printf "%.0f", s}' "$W/d.bg")
    printf '  main-chrom reads %s, track carries %s\n' "$exp" "$got"
    awk -v g="$got" -v e="$exp" 'BEGIN{r=g/e; if (r<0.995||r>1.005){
        printf "ERROR: track sum / read count = %.4f, expected 1.000\n", r > "/dev/stderr"; exit 1}}'
    bedGraphToBigWig "$W/d.bg" "$CHR" "$OUT/thp1_dnase_5p.bw"
    echo "  wrote thp1_dnase_5p.bw"
fi

# Stranded H3K27ac target: pooled for training, per-replicate for the ceiling.
echo "=== THP-1 H3K27ac 5-prime stranded targets ==="
filt "$K27_R1" "$W/k1.bam"
filt "$K27_R2" "$W/k2.bam"
samtools merge -@ 8 -f "$W/k_pool.bam" "$W/k1.bam" "$W/k2.bam"
samtools index -@ 8 "$W/k_pool.bam"

make_stranded () {  # bam prefix
    local tot=0
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${2}_${tag}.bw"
        [[ -s "$out" ]] && { echo "    exists $(basename "$out")"; continue; }
        bedtools genomecov -ibam "$1" -5 -dz -strand "$strand" \
          | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
          | LC_COLLATE=C sort -k1,1 -k2,2n -T "$W" > "$W/x.bg"
        local s; s=$(awk '{s += ($3-$2)*$4} END {printf "%.0f", s}' "$W/x.bg")
        tot=$(( tot + s ))
        bedGraphToBigWig "$W/x.bg" "$CHR" "$out"; rm -f "$W/x.bg"
        echo "    wrote $(basename "$out") (sum $s)"
    done
    local exp; exp=$(main_chrom_reads "$1")
    printf '  %s: main-chrom reads %s, plus+minus %s\n' "$2" "$exp" "$tot"
}

make_stranded "$W/k_pool.bam" thp1_h3k27ac_5p
make_stranded "$W/k1.bam"     thp1_h3k27ac_rep1_5p
make_stranded "$W/k2.bam"     thp1_h3k27ac_rep2_5p

echo
echo "=== candidate elements from the DNase, matching the K562 DNase-derived convention ==="
EL=$D/reference/THP1_DNase_candidate_elements.narrowPeak
if [[ -s "$EL" ]]; then
    echo "  exists: $(wc -l < "$EL") elements"
else
    [[ -s "$W/dnase.bam" ]] || filt "$DNASE" "$W/dnase.bam"
    $MACS2 callpeak -t "$W/dnase.bam" -f BAM -g hs -n thp1_dnase \
        --outdir "$W/macs2" --nomodel --shift -75 --extsize 150 -q 0.01 2>&1 | tail -3
    # load_peaks wants 10 columns and centres on start + summit; the K562 file uses the region
    # MIDPOINT as its summit (0.26), so the same convention is applied here rather than the
    # MACS2 summit, or the two element sets would differ in centring as well as in caller.
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
            mid = int(($3 - $2) / 2)
            print $1, $2, $3, $1":"$2"-"$3, 0, ".", 0, -1, -1, mid
         }' "$W/macs2/thp1_dnase_peaks.narrowPeak" \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
      | LC_COLLATE=C sort -k1,1 -k2,2n > "$EL"
    echo "  wrote $EL ($(wc -l < "$EL") elements)"
fi

echo
echo "=== H3K27ac inter-replicate ceiling (the one this panel is read against) ==="
"$D/.pixi/envs/multimodal/bin/python" scripts/0.25.profile_ceiling_by_binsize.py \
    --label thp1_h3k27ac \
    --rep1-plus "$OUT/thp1_h3k27ac_rep1_5p_plus.bw" \
    --rep1-minus "$OUT/thp1_h3k27ac_rep1_5p_minus.bw" \
    --rep2-plus "$OUT/thp1_h3k27ac_rep2_5p_plus.bw" \
    --rep2-minus "$OUT/thp1_h3k27ac_rep2_5p_minus.bw" \
    --elements "$EL"
echo done
