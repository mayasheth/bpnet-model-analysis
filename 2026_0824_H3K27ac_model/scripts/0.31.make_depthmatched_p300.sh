#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH -o log/depthmatch.%j.txt
#SBATCH -e log/depthmatch.%j.txt
#SBATCH --job-name=depthmatch
#
# Build a K562 EP300 target matched to GM12878's data budget, to test whether training-signal
# volume explains the transfer asymmetry.
#
# THE ASYMMETRY. K562-trained p300 transferred to GM12878 keeps +0.207 over that cell type's
# own accessibility model (83% of its local advantage); GM12878-trained transferred to K562
# keeps +0.006 (p=0.72). Both models are equally strong at home, so it is not model quality.
# The clearest difference between the two experiments is volume: K562 has 2.30x the reads in
# peaks (3,590,540 against 1,562,462), from 1.70x the depth over 1.35x the peaks. FRiP does
# NOT separate them once replicate spread is accounted for, so the hypothesis under test is
# signal quantity, not purity.
#
# THE TEST. Give K562 GM12878's budget and retrain. If the subsampled K562 model stops
# transferring, volume explains the asymmetry. If it still transfers, something else about
# K562 as a training cell type does.
#
# TRACK-CONSTRUCTION CONTROL. The existing K562 target (2025_0703/data/ENCSR000EGE_*.bigWig)
# was built by an earlier script. If the subsampled track were built differently, the
# comparison would confound depth with construction. So this builds a FULL-DEPTH track by the
# same code path as the subsampled one and reports its agreement with the existing file; only
# if they agree can the existing model serve as the full-depth arm.
#
# PEAK SUBSAMPLING is random with a fixed seed rather than top-N by score. Random ablates the
# number of positives while leaving the score distribution alone. Re-calling peaks on the
# subsampled BAM would be more faithful to how GM12878 ended up with 21,068, and is the
# follow-up if this test comes back positive.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
BIN=$D/.pixi/envs/multimodal/bin
export PATH="$BIN:$PATH"
E=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562
CHR_SIZES=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
OUT=$P/data/p300_depthmatched
mkdir -p "$OUT" "$P/log"
W=${SCRATCH}/depthmatch_$$; mkdir -p "$W"; trap "rm -rf $W" EXIT
cd "$P"

REP1=$E/ENCFF466WKF.filtered.sorted.bam
REP2=$E/ENCFF163FSR.filtered.sorted.bam
# Targets measured by 0.30: K562 51,127,010 mapped over 28,532 peaks; GM12878 30,001,681
# over 21,068.
TARGET_READS=30001681
TARGET_PEAKS=21068
SEED=42

k562_total=$(( $(samtools view -c -F 0x400 -@ 8 "$REP1") + $(samtools view -c -F 0x400 -@ 8 "$REP2") ))
FRAC=$(echo "scale=6; $TARGET_READS / $k562_total" | bc)
echo "K562 mapped $k562_total, target $TARGET_READS, keeping fraction $FRAC"

# --- merged full-depth and subsampled BAMs -----------------------------------
samtools merge -@ 8 -f "$W/full.bam" "$REP1" "$REP2"
samtools index -@ 8 "$W/full.bam"
samtools view -@ 8 -b -s "${SEED}${FRAC#0}" "$W/full.bam" > "$W/sub.bam"
samtools index -@ 8 "$W/sub.bam"
echo "full  $(samtools view -c -F 0x400 -@ 8 "$W/full.bam")"
echo "sub   $(samtools view -c -F 0x400 -@ 8 "$W/sub.bam")"

# --- stranded 5'-end bigwigs, same convention as 0.5 -------------------------
make_5p () {
    local bam=$1 prefix=$2
    for strand in + -; do
        local tag; [[ "$strand" == "+" ]] && tag=plus || tag=minus
        local out="$OUT/${prefix}_${tag}.bw"
        [[ -s "$out" ]] && { echo "  $out exists"; continue; }
        bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
            | awk 'NR==FNR{c[$1]=1; next} ($1 in c){print $1, $2, $2+1, $3}' \
                  OFS='\t' "$CHR_SIZES" - \
            | sort -k1,1 -k2,2n -T "$W" > "$W/${prefix}_${tag}.bedGraph"
        bedGraphToBigWig "$W/${prefix}_${tag}.bedGraph" "$CHR_SIZES" "$out"
        echo "  wrote $out"
    done
}
make_5p "$W/full.bam" ep300_fulldepth
make_5p "$W/sub.bam"  ep300_depthmatched

# --- peaks subsampled to GM12878's count -------------------------------------
PK=$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak
shuf --random-source=/dev/zero -n "$TARGET_PEAKS" "$PK" | sort -k1,1 -k2,2n \
    > "$OUT/ep300_peaks_depthmatched.narrowPeak"
echo "peaks $(wc -l < "$PK") -> $(wc -l < "$OUT/ep300_peaks_depthmatched.narrowPeak")"

# --- construction control ----------------------------------------------------
echo
echo "=== does the full-depth rebuild match the existing target? ==="
"$BIN/python" - <<'PY'
import numpy as np, pyBigWig
D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
OUT = f"{D}/2026_0824_H3K27ac_model/data/p300_depthmatched"
pairs = [("plus", f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig"),
         ("minus", f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig")]
ok = True
for tag, old_fp in pairs:
    new_fp = f"{OUT}/ep300_fulldepth_{tag}.bw"
    o, n = pyBigWig.open(old_fp), pyBigWig.open(new_fp)
    so, sn = o.header()["sumData"], n.header()["sumData"]
    v_o = np.nan_to_num(np.array(o.values("chr1", 1_000_000, 4_000_000), dtype=float))
    v_n = np.nan_to_num(np.array(n.values("chr1", 1_000_000, 4_000_000), dtype=float))
    r = np.corrcoef(v_o, v_n)[0, 1]
    print(f"  {tag:<6} existing sum {so:,.0f}  rebuilt sum {sn:,.0f}  "
          f"ratio {sn/so:.4f}  r(chr1 3Mb) {r:.4f}")
    if abs(sn/so - 1) > 0.02 or r < 0.99:
        ok = False
    o.close(); n.close()
print()
print("MATCH: the existing K562 p300 model can serve as the full-depth arm; train only the"
      "\n  subsampled arm." if ok else
      "MISMATCH: construction differs, so BOTH arms must be retrained from these tracks.")
PY
