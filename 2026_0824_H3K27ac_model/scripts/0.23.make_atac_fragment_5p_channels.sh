#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/atacfrag5p.%j.txt
#SBATCH -e log/atacfrag5p.%j.txt
#SBATCH --job-name=atacfrag5p
#
# Fragment-size-stratified ATAC channels in the ChromBPNet 5' convention (K562).
#
# WHY THIS REPLACES 0.9. The channels built by 0.9.make_atac_fragment_channels.sh are
# `genomecov -bg` over the full fragment interval, the same read-length-dependent smear
# that 0.20 removed from the flat ATAC track. Feeding those alongside atac_5p.bw would mix
# the two incompatible accessibility conventions inside a single input tensor. These
# channels are single-base Tn5 insertion counts, so every accessibility channel the model
# sees is on the same footing. atac_sub.bw / atac_mono.bw are left in place but are
# superseded; do not use them with a 5'-input model.
#
# WHY IT CANNOT REGRESS AGAINST accs5p. The model input is
#   [all, sub, mono, di, poly]
# where `all` is atac_5p.bw itself and sub/mono/di/poly exhaustively partition fragments by
# length. The first convolution is a linear map across channels, so it can reproduce the
# accs5p input exactly by zeroing the four stratified channels. Any difference in the
# fitted model is added information, not a changed input.
#
# BIN EDGES from results/atac_fragment_length_hist.txt (K562 rep1, chr1): NFR mode ~40 bp,
# trough ~140, mono-nucleosome mode ~200, trough ~330, di-nucleosome mode ~390. The 0.9
# edges (<100 and 180-247) sat inside the modes and covered only half the fragments; these
# sit in the troughs and cover all of them.
#   sub   <= 139     nucleosome-free / TF-bound
#   mono  140-329    mono-nucleosomal
#   di    330-620    di-nucleosomal
#   poly  >= 621     higher-order
#
# SHIFT. The PE BAMs are unshifted alignments; atac_5p.bw came from already-shifted
# tagAligns. PLUS_DELTA/MINUS_DELTA must therefore be set to the offsets measured by
# 0.22.calibrate_tn5_shift.sh -- there is no default, because guessing puts the stratified
# channels a few bp out of phase with the flat channel.
#
# Usage: PLUS_DELTA=4 MINUS_DELTA=-5 sbatch 0.23.make_atac_fragment_5p_channels.sh

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
BAM_DIR="${SCRATCH}/atac_pe"
PLUS_DELTA=${PLUS_DELTA:?set PLUS_DELTA from results/tn5_shift_calibration.txt}
MINUS_DELTA=${MINUS_DELTA:?set MINUS_DELTA from results/tn5_shift_calibration.txt}

WORK="${SCRATCH}/atacfrag5p_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT
mkdir -p "$P/data" "$P/log" "$P/results"

BAMS=("$BAM_DIR/ENCFF077FBI.pe.bam" "$BAM_DIR/ENCFF128WZG.pe.bam" "$BAM_DIR/ENCFF534DCE.pe.bam")
for b in "${BAMS[@]}"; do [[ -s "$b" ]] || { echo "ERROR missing $b" >&2; exit 1; }; done
echo "plus_delta=$PLUS_DELTA minus_delta=$MINUS_DELTA"

# One pass over the BAMs: emit `class<TAB>chrom<TAB>pos` for both insertion sites of every
# properly-paired, non-duplicate fragment. TLEN>0 selects the leftmost mate, so each
# fragment is written exactly once.
echo "=== dumping insertion sites ==="
for b in "${BAMS[@]}"; do
    samtools view -@ 2 -f 0x2 -F 0x400 "$b"
done | awk -v OFS='\t' -v dp="$PLUS_DELTA" -v dm="$MINUS_DELTA" '
    $9 > 0 {
        L = $9; s = $4 - 1; e = s + L
        if      (L <= 139) k = "sub"
        else if (L <= 329) k = "mono"
        else if (L <= 620) k = "di"
        else               k = "poly"
        n[k]++
        left = s + dp; right = e - 1 + dm
        if (left  >= 0) print k, $3, left
        if (right >= 0) print k, $3, right
    }
    END { for (k in n) printf "FRAGCOUNT\t%s\t%d\n", k, n[k] > "/dev/stderr" }
' 2> "$WORK/counts.txt" > "$WORK/all_sites.txt"
grep FRAGCOUNT "$WORK/counts.txt" | tee "$P/results/atac_fragment_class_counts.txt"

for k in sub mono di poly; do
    out="$P/data/atac_${k}5p.bw"
    echo "=== channel $k -> $(basename $out) ==="
    awk -v k="$k" -v OFS='\t' '$1 == k { print $2, $3, $3 + 1 }' "$WORK/all_sites.txt" \
      | awk 'NR==FNR{c[$1]=$2; next} ($1 in c) && $3 <= c[$1]' "$CHR" - \
      | sort -k1,1 -k2,2n -T "$WORK" \
      | bedtools genomecov -i stdin -g "$CHR" -bg \
      | LC_COLLATE=C sort -k1,1 -k2,2n -T "$WORK" > "$WORK/$k.bg"
    bedGraphToBigWig "$WORK/$k.bg" "$CHR" "$out"
    rm -f "$WORK/$k.bg"
    echo "  wrote $out"
done

ls -la "$P/data"/atac_{sub,mono,di,poly}5p.bw
