#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH -o log/atac_5p.%j.txt
#SBATCH -e log/atac_5p.%j.txt
#SBATCH --job-name=k27_atac5p
#
# Build ChromBPNet-convention accessibility tracks: single-base 5' insertion counts.
#
# WHY. Our existing atac.bw files use `bedtools genomecov -bg` over the FULL tagAlign
# interval, so coverage is proportional to read length. K562 and GM12878 tagAligns are
# 94-95 bp; TeloHAEC's are 35-36 bp. K562/GM12878 match each other, which is why their
# transfer work was unaffected -- but TeloHAEC differs by ~2.6x from read length alone,
# on top of ~2x depth. Same defect one level up as the 250 bp fragment extension rejected
# for the H3K27ac target: a read-length-dependent smear that breaks cross-sample
# comparability.
#
# ChromBPNet (chrombpnet/helpers/preprocessing/reads_to_bigwig.py), verbatim:
#   awk (strand shift) | bedtools genomecov -bg -5 -i stdin -g <genome>
#   ATAC:  plus_shift_delta, minus_shift_delta = 4-plus_shift, -4-minus_shift
#   DNase: plus_shift_delta, minus_shift_delta = -plus_shift, 1-minus_shift
# Single-base 5' ends, unstranded, `-5` applied AFTER the shift. Their pipeline auto-detects
# an existing shift; our inputs are already Tn5-shifted (*.tn5.sorted.tagAlign.gz), so the
# deltas are zero and the correct build reduces to adding `-5`. No further shift is applied
# here -- applying +4/-4 again would double-shift.
#
# NON-DESTRUCTIVE: writes *_5p.bw. Never overwrite atac.bw --
# 2026_0529_multimodal_p300_model/data/atac.bw is shared with the p300 models, and every
# current ATAC-input result is on it.
#
# Usage: sbatch 0.20.make_atac_5prime_bigwigs.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
DATA=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
T=/oak/stanford/groups/engreitz/Projects/E2G/endothelial_cells/TeloHAEC_GSE210489_GSE210491
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
WORK="${SCRATCH}/atac5p_$$"; mkdir -p "$WORK"; trap "rm -rf $WORK" EXIT
mkdir -p "$P/log" "$P/data/panel"

build () {
    local out="$1"; shift
    if [[ -s "$out" ]]; then echo "  exists: $(basename $out)"; return; fi
    local files=("$@")
    for f in "${files[@]}"; do
        [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }
    done
    echo "  building $(basename $out) from ${#files[@]} tagAlign file(s)"
    zcat "${files[@]}" \
      | awk 'NR==FNR{c[$1]=1; next} ($1 in c)' "$CHR" - \
      | sort -k1,1 -k2,2n -T "$WORK" \
      | bedtools genomecov -bg -5 -i stdin -g "$CHR" \
      | LC_COLLATE=C sort -k1,1 -k2,2n -T "$WORK" > "$WORK/x.bg"
    bedGraphToBigWig "$WORK/x.bg" "$CHR" "$out"
    rm -f "$WORK/x.bg"
    echo "  wrote $(basename $out)"
}

echo "=== K562 (94-95 bp reads) ==="
build "$D/2026_0529_multimodal_p300_model/data/atac_5p.bw" \
  "$DATA/K562/ENCFF077FBI.tn5.sorted.tagAlign.gz" \
  "$DATA/K562/ENCFF128WZG.tn5.sorted.tagAlign.gz" \
  "$DATA/K562/ENCFF534DCE.tn5.sorted.tagAlign.gz"

echo "=== GM12878 (94-95 bp reads) ==="
build "$D/2026_0606_GM12878_transferability/data/atac_5p.bw" \
  "$DATA/GM12878/ATAC/ENCFF440GRZ.tn5.sorted.tagAlign.gz" \
  "$DATA/GM12878/ATAC/ENCFF962FMH.tn5.sorted.tagAlign.gz" \
  "$DATA/GM12878/ATAC/ENCFF981FXV.tn5.sorted.tagAlign.gz"

# TeloHAEC-only accessions; the 3 Eahy926 files in TeloHAEC_ctrl/ATAC are excluded.
echo "=== TeloHAEC ctrl (35-36 bp reads; Eahy926 excluded) ==="
build "$P/data/panel/TeloHAEC_ctrl_atac_5p.bw" \
  $T/TeloHAEC_ctrl/ATAC/SRR20809416.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_ctrl/ATAC/SRR20809419.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_ctrl/ATAC/SRR20809420.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_ctrl/ATAC/SRR20809430.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_ctrl/ATAC/SRR20809431.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_ctrl/ATAC/SRR20809433.tn5.sorted.tagAlign.gz

echo "=== TeloHAEC IL1b ==="
build "$P/data/panel/TeloHAEC_IL1b_atac_5p.bw" \
  $T/TeloHAEC_IL1b/ATAC/SRR20809411.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_IL1b/ATAC/SRR20809414.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_IL1b/ATAC/SRR20809415.tn5.sorted.tagAlign.gz

echo "=== TeloHAEC TNFa ==="
build "$P/data/panel/TeloHAEC_TNFa_atac_5p.bw" \
  $T/TeloHAEC_TNFa/ATAC/SRR20809409.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_TNFa/ATAC/SRR20809410.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_TNFa/ATAC/SRR20809412.tn5.sorted.tagAlign.gz

echo "=== TeloHAEC no-VEGF ==="
build "$P/data/panel/TeloHAEC_VEGF_atac_5p.bw" \
  $T/TeloHAEC_VEGF_ctrl/ATAC/SRR20809407.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_VEGF_ctrl/ATAC/SRR20809408.tn5.sorted.tagAlign.gz \
  $T/TeloHAEC_VEGF_ctrl/ATAC/SRR20809413.tn5.sorted.tagAlign.gz

echo "ALL_DONE"
