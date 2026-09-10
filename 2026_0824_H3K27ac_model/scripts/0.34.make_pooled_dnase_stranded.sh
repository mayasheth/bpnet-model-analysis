#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/dnase_stranded.%j.txt
#SBATCH -e log/dnase_stranded.%j.txt
#SBATCH --job-name=dnase_str
#
# Pooled STRANDED DNase 5' tracks: the training target for the ATAC->DNase converter.
#
# What exists already and why neither file will do:
#   0.32  pooled but UNSTRANDED (k562_dnase_5p.bw). That is the model's accessibility
#         INPUT. A target has to be stranded to match how every other target in this
#         project is built (0.5, 0.3), because bpnetlite's profile head is trained per
#         strand; the counts head sums the channels, so counts are unaffected either way.
#   0.33  stranded but PER REPLICATE, and only the two ENCSR000EOT runs for K562. Training
#         on one replicate throws away half the depth and, worse, would put the converter's
#         target on a different read set from the DNase input the swap experiment used.
#
# So: same BAM sets as 0.32, same read-1 filtering, pooled the same way, split by strand.
#
# POOL BY MERGING BAMs, NEVER BY COMBINING bedGraphs. `genomecov` emits bookended intervals
# and `bedtools merge` fuses and sums them, which produced a K562 track 34x too hot with
# exit 0 and a valid-looking bigwig. Ten GPU folds were trained against it before the totals
# gave it away. The guard is the assertion at the end: a 5' track carries exactly one count
# per read, so plus + minus must equal the main-chromosome read count.
#
# THE DENOMINATOR COMES FROM $CHR, NOT A REGEX. GRCh38.main.chrom.sizes INCLUDES chrM, and
# GM12878's DNase is 9.5% chrM against K562's 0.075%. A hardcoded chr1-22,X,Y denominator
# fails on GM12878; a whole-BAM count fails the other way. Reading the contig list from the
# same file the track is built against is the only version that cannot drift.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT=$P/data
mkdir -p "$OUT" "$P/log"
W=${SCRATCH}/dnase_stranded_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

build () {
    local prefix=$1; shift
    local bams=("$@")
    local plus="$OUT/${prefix}_5p_plus.bw" minus="$OUT/${prefix}_5p_minus.bw"
    if [[ -s "$plus" && -s "$minus" ]]; then
        echo "  exist: $(basename "$plus"), $(basename "$minus")"; return
    fi
    local filtered=() expected=0
    for b in "${bams[@]}"; do
        [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }
        local n_paired tag kept f
        n_paired=$(samtools view -c -f 0x1 -@ 8 "$b")
        f="$W/$(basename "$b" .bam).f.bam"
        # Pairedness is detected per BAM, not read off the metadata table: that table records
        # the EXPERIMENT's run type and ENCSR000EOT is listed as "pe, se". Both mates of a
        # pair would each contribute a 5' end and double-count the fragment.
        if [[ "$n_paired" -gt 0 ]]; then
            tag="paired -> read 1 only"
            samtools view -b -F 0x900 -f 0x40 -@ 8 "$b" > "$f"
        else
            tag="single-end"
            samtools view -b -F 0x900 -@ 8 "$b" > "$f"
        fi
        samtools index -@ 8 "$f"
        kept=$(samtools idxstats "$f" \
            | awk 'NR==FNR {keep[$1]=1; next} ($1 in keep) {s += $3} END {print s+0}' \
                  "$CHR" -)
        printf '  %-34s %-22s main-chrom %11d\n' "$(basename "$b")" "$tag" "$kept"
        filtered+=("$f"); expected=$(( expected + kept ))
    done
    if [[ ${#filtered[@]} -gt 1 ]]; then
        samtools merge -@ 8 -f "$W/pool.bam" "${filtered[@]}"
    else
        cp "${filtered[0]}" "$W/pool.bam"
    fi

    local got_total=0
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_5p_${tag}.bw"
        # -5 -dz -strand, matching 0.5 and 0.3 so the converter's target is built exactly
        # the way every other target in this project is.
        bedtools genomecov -ibam "$W/pool.bam" -5 -dz -strand "$strand" \
            | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' OFS='\t' "$CHR" - \
            | LC_COLLATE=C sort -k1,1 -k2,2n -T "$W" > "$W/x.bg"
        local s; s=$(awk '{s += ($3-$2)*$4} END {printf "%.0f", s}' "$W/x.bg")
        got_total=$(( got_total + s ))
        printf '    %-28s carries %11d\n' "$tag" "$s"
        bedGraphToBigWig "$W/x.bg" "$CHR" "$out"; rm -f "$W/x.bg"
    done

    printf '  %-34s main-chrom reads %d, plus+minus %d\n' "sanity" "$expected" "$got_total"
    awk -v g="$got_total" -v e="$expected" 'BEGIN {r = g/e; if (r < 0.995 || r > 1.005) {
        printf "ERROR: (plus+minus) / main-chrom read count = %.4f, expected 1.000.\n", r > "/dev/stderr";
        print  "  Both sides count the same reads, so any gap is a pooling or strand bug." > "/dev/stderr";
        exit 1}}'
    rm -f "$W/pool.bam" "${filtered[@]}" "${filtered[@]/%/.bai}"
    echo "  wrote ${prefix}_5p_{plus,minus}.bw"
}

echo "=== K562 DNase (ENCSR000EOT pe+se, ENCSR000EKS se) -- same BAMs as 0.32 ==="
build k562_dnase \
    "$E/K562/ENCFF205FNC.filtered.sorted.bam" \
    "$E/K562/ENCFF860XAE.filtered.sorted.bam" \
    "$E/K562/ENCFF987IUK.filtered.sorted.bam"

echo "=== GM12878 DNase (ENCSR000EMT se) -- same BAMs as 0.32 ==="
build gm12878_dnase \
    "$E/GM12878/ENCFF467CXY_sorted.bam" \
    "$E/GM12878/ENCFF940NSD_sorted.bam"

echo
echo "=== stranded targets against the unstranded input they were pooled from ==="
"$D/.pixi/envs/multimodal/bin/python" - <<'PY'
import pyBigWig
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
for cell in ("k562", "gm12878"):
    tot = 0
    for tag in ("plus", "minus"):
        b = pyBigWig.open(f"{P}/data/{cell}_dnase_5p_{tag}.bw")
        s = b.header()["sumData"]; tot += s; b.close()
        print(f"  {cell} {tag:<6} sum {s:>14,.0f}")
    b = pyBigWig.open(f"{P}/data/{cell}_dnase_5p.bw")
    u = b.header()["sumData"]; b.close()
    # The unstranded 0.32 track was pooled from the same reads, so this is a cross-script
    # check that the two agree and not merely an internal one.
    flag = "OK" if abs(tot - u) / u < 0.005 else "MISMATCH"
    print(f"  {cell} plus+minus {tot:>14,.0f}  vs 0.32 unstranded {u:>14,.0f}  [{flag}]")
PY
