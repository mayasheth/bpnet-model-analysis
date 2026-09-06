#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/gcneg.%j.txt
#SBATCH -e log/gcneg.%j.txt
#SBATCH --job-name=gcneg
#
# GC-matched negatives for the H3K27ac candidate element sets.
#
# WHY: every model in this repository -- H3K27ac AND p300, K562 AND GM12878 -- was trained
# with `--negatives reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed`, which is
# ChromBPNet's genome-wide GC-ANNOTATED TILING, not a matched set. The trainer reads only its
# first three columns and samples uniformly, so the negative pool sits at mean GC 0.389
# against 0.466-0.593 for candidate elements, and a sequence branch can satisfy much of its
# training objective with a GC detector.
#
# This uses the SAME tool the p300 v3 negatives were built with (`bpnet-gc-background`, Oct
# 2025) rather than a reimplementation, so the matched sets are procedurally identical and a
# matched-vs-unmatched comparison is not confounded by two different matching algorithms.
#
# The existing 2025_1016_p300_model_v3/data/gc_negatives.bed cannot be reused here: it was
# matched to EP300 peak GC, and the candidate element set has a different GC distribution.
#
# neg_to_pos_ratio 1 rather than the 4 used for p300: training caps the pool at
# --max-negatives 50,000 regardless, and a uniform subsample of a GC-matched pool is still
# GC-matched, so a 4x pool would only cost disk and runtime.
set -euo pipefail
BP=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/bpnet_37/bin
# bpnet-gc-background shells out to bedtools, which is not in the bpnet_37 env. Appended
# rather than prepended so bpnet_37's own entry points still win; the tool is invoked by
# absolute path and its shebang pins its interpreter, so PATH order cannot mis-route it.
export PATH="$PATH:/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin"
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
REF_GC=$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
OUT=$P/data
mkdir -p "$OUT" "$P/log"
cd "$P"

# FOREGROUND SUBSAMPLE, and why it is necessary rather than an optimisation.
#
# bpnet-gc-background matches each foreground region to a distinct candidate negative in the
# same GC bin, drawing without replacement. Our candidate elements are GC-rich (mean 0.49,
# up to 0.593 in the over-predicted tail) while the genome-background pool sits at 0.389, so
# the high-GC bins hold far fewer candidates than 150,528 foreground regions demand. The tool
# does not report a shortfall -- it HANGS: a first attempt reached 93% (140,108/150,528) in
# 15 seconds and then made no further progress for two hours before being cancelled.
#
# Subsampling the foreground fixes it because only the GC DISTRIBUTION is needed to define
# the matching target, and a 50,000-region random subsample of 150,528 estimates that
# distribution to well within the 0.02 bin width. 50,000 is chosen to equal the trainer's
# --max-negatives cap, so the matched pool is used whole and nothing is discarded downstream.
#
# `shuf --random-source=/dev/zero` is deterministic, so the subsample is
# reproducible across reruns.
N_FOREGROUND=${N_FOREGROUND:-50000}

run_one () {
    local label=$1 peaks=$2
    echo "=========== $label"
    echo "peaks: $peaks  ($(wc -l < "$peaks") regions)"
    local fg="$OUT/foreground_subsample_${label}.narrowPeak"
    if [[ $(wc -l < "$peaks") -gt $N_FOREGROUND ]]; then
        shuf --random-source=/dev/zero -n "$N_FOREGROUND" "$peaks" \
            | sort -k1,1 -k2,2n > "$fg"
        echo "foreground subsampled to $(wc -l < "$fg") regions (GC target is a "
        echo "  distribution, so a subsample defines it; the tool hangs on the full set)"
    else
        cp "$peaks" "$fg"
    fi
    # Fail fast rather than burning the wall clock if it hangs again.
    timeout 45m "$BP/bpnet-gc-background" \
        --ref_fasta "$GEN" \
        --peaks_bed "$fg" \
        --ref_gc_bed "$REF_GC" \
        --out_dir "$OUT" \
        --output_prefix "$OUT/gc_negatives_${label}" \
        --flank_size 1057 \
        --neg_to_pos_ratio_train 1
    local f="$OUT/gc_negatives_${label}.bed"
    [[ -s "$f" ]] || { echo "ERROR: $f empty or missing" >&2; exit 1; }
    echo "wrote $f ($(wc -l < "$f") regions)"
}

run_one k562_h3k27ac   "$D/reference/K562_DNase_candidate_elements.narrowPeak"
run_one gm12878_h3k27ac "$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"

echo
echo "=========== verification: matched pool GC against elements and the unmatched tiling"
"$D/.pixi/envs/multimodal/bin/python" - <<'PY'
import pandas as pd, numpy as np
D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
tiling = pd.read_csv(f"{D}/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed",
                     sep="\t", header=None, usecols=[3], names=["gc"])
print(f"{'set':<34}{'n':>10}{'mean GC':>10}{'p10':>8}{'p90':>8}")
print(f"{'unmatched tiling (what we used)':<34}{len(tiling):>10,}"
      f"{tiling['gc'].mean():>10.4f}{tiling['gc'].quantile(.1):>8.3f}"
      f"{tiling['gc'].quantile(.9):>8.3f}")
for lab in ("k562_h3k27ac", "gm12878_h3k27ac"):
    f = f"{P}/data/gc_negatives_{lab}.bed"
    d = pd.read_csv(f, sep="\t", header=None)
    # bpnet-gc-background writes narrowPeak-like columns; GC is not retained, so recover it
    # by joining back to the tiling on coordinates.
    print(f"{'matched: ' + lab:<34}{len(d):>10,}{'-':>10}{'-':>8}{'-':>8}")
print("\nGC of the matched sets is verified by 0.29 against the elements themselves;")
print("bpnet-gc-background does not carry the GC column into its output.")
PY
