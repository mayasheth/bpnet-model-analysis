#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH -o log/dnase5p.%j.txt
#SBATCH -e log/dnase5p.%j.txt
#SBATCH --job-name=dnase5p
#
# DNase 5'-insertion accessibility inputs for K562 and GM12878, to test swapping DNase for
# ATAC as the model's accessibility channel.
#
# WHY: DNase couples to H3K27ac and to enhancer activity better than ATAC does, and observed
# DNase already beats observed ATAC inside ABC. If a DNase-input model beats its ATAC
# counterpart in-cell or on transfer, that also decides whether the ATAC->DNase converter is
# worth building, since a converter is only worth having if DNase is the better representation
# to be in.
#
# CONVENTION: `genomecov -bg -5`, unstranded, matching how atac_5p.bw was built by 0.20 so the
# swap changes the assay and nothing else. The 5' convention is what makes this safe across
# experiments with different read lengths and run types: one count per read regardless of
# length, which is why K562's mixed pe/se DNase can be pooled with no read-length confound.
#
# PAIRED-END HANDLING: for a paired BAM, both mates would contribute a 5' end and double-count
# each fragment, so paired BAMs are filtered to read 1. Pairedness is detected per BAM rather
# than taken from the metadata table, because that table records the EXPERIMENT's run type and
# ENCSR000EOT is listed as "pe, se".
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT=$P/data
mkdir -p "$OUT" "$P/log"
W=${SCRATCH}/dnase5p_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

# Replicates are pooled by MERGING THE BAMs and calling genomecov once, not by summing
# per-replicate bedGraphs.
#
# The bedGraph route is wrong and fails silently. `genomecov -bg` emits runs of constant
# coverage, so adjacent intervals are bookended (one's end is the next's start). `bedtools
# merge` collapses bookended intervals by default, so summing that way fuses long runs into
# single intervals and adds their values: the first attempt produced a K562 track with sum
# 10,290,635,309 from 301,152,633 reads, 34x too high, over 165M covered bp instead of ~288M.
# A 5' insertion track must have sum == read count, which is the assertion at the end here.
build () {
    local out=$1; shift
    local bams=("$@")
    [[ -s "$out" ]] && { echo "  exists: $(basename "$out")"; return; }
    local filtered=() expected=0
    for b in "${bams[@]}"; do
        [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }
        local n_paired total tag kept f
        n_paired=$(samtools view -c -f 0x1 -@ 8 "$b")
        total=$(samtools view -c -F 0x900 -@ 8 "$b")
        f="$W/$(basename "$b" .bam).f.bam"
        if [[ "$n_paired" -gt 0 ]]; then
            tag="paired -> read 1 only"
            samtools view -b -F 0x900 -f 0x40 -@ 8 "$b" > "$f"
        else
            tag="single-end"
            samtools view -b -F 0x900 -@ 8 "$b" > "$f"
        fi
        # The denominator must count exactly the contigs the track keeps, which is whatever
        # is in $CHR. Two earlier versions got this wrong in opposite directions: a whole-BAM
        # count failed on GM12878, whose DNase is 9.5% chrM against K562's 0.075%, and a
        # hardcoded chr1-22,X,Y regex then failed the other way because $CHR includes chrM.
        # Reading the contig list from $CHR is the only version that cannot drift from it.
        samtools index -@ 8 "$f"
        kept=$(samtools idxstats "$f" \
            | awk 'NR==FNR {keep[$1]=1; next} ($1 in keep) {s += $3} END {print s+0}' \
                  "$CHR" -)
        local allp; allp=$(samtools view -c -@ 8 "$f")
        printf '  %-34s %-22s total %11d  primary %11d  main-chrom %11d\n' \
            "$(basename "$b")" "$tag" "$total" "$allp" "$kept"
        filtered+=("$f"); expected=$(( expected + kept ))
    done
    if [[ ${#filtered[@]} -gt 1 ]]; then
        samtools merge -@ 8 -f "$W/pool.bam" "${filtered[@]}"
    else
        cp "${filtered[0]}" "$W/pool.bam"
    fi
    bedtools genomecov -ibam "$W/pool.bam" -bg -5 \
        | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
        | LC_COLLATE=C sort -k1,1 -k2,2n -T "$W" > "$W/pool.bg"
    bedGraphToBigWig "$W/pool.bg" "$CHR" "$out"
    # sum == reads on the main chromosomes, so a silent merge bug cannot pass again
    local got
    got=$(awk '{s += ($3-$2)*$4} END {printf "%.0f", s}' "$W/pool.bg")
    printf '  %-34s main-chrom reads %d, track carries %s\n' "sanity" "$expected" "$got"
    awk -v g="$got" -v e="$expected" 'BEGIN {r = g/e; if (r < 0.995 || r > 1.005) {
        printf "ERROR: track sum / main-chrom read count = %.4f, expected 1.000. Both sides now\n", r;
        print  "  count the same reads, so any gap is a pooling or interval bug." > "/dev/stderr"; exit 1}}'
    rm -f "$W/pool.bam" "$W/pool.bg" "${filtered[@]}"
    echo "  wrote $(basename "$out")"
}

echo "=== K562 DNase (ENCSR000EOT pe+se, ENCSR000EKS se) ==="
build "$OUT/k562_dnase_5p.bw" \
    "$E/K562/ENCFF205FNC.filtered.sorted.bam" \
    "$E/K562/ENCFF860XAE.filtered.sorted.bam" \
    "$E/K562/ENCFF987IUK.filtered.sorted.bam"

echo "=== GM12878 DNase (ENCSR000EMT se) ==="
build "$OUT/gm12878_dnase_5p.bw" \
    "$E/GM12878/ENCFF467CXY_sorted.bam" \
    "$E/GM12878/ENCFF940NSD_sorted.bam"

echo
echo "=== totals against the ATAC inputs they replace ==="
"$D/.pixi/envs/multimodal/bin/python" - <<'PY'
import pyBigWig
D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
for lab, fp in [
    ("K562 ATAC 5'",    f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw"),
    ("K562 DNase 5'",   f"{P}/data/k562_dnase_5p.bw"),
    ("GM12878 ATAC 5'", f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw"),
    ("GM12878 DNase 5'", f"{P}/data/gm12878_dnase_5p.bw"),
]:
    try:
        b = pyBigWig.open(fp)
        h = b.header()
        print(f"  {lab:<18} sum {h['sumData']:>14,.0f}   covered bp {h['nBasesCovered']:>13,}")
        b.close()
    except Exception as e:
        print(f"  {lab:<18} UNREADABLE: {e}")
PY
