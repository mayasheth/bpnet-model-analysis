#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -o log/testaccnorm.%j.txt
#SBATCH -e log/testaccnorm.%j.txt
#SBATCH --job-name=test_accnorm
#
# Regression gate for the 2026-09-17 additions to 4.1: --acc-mean/--acc-std and
# --count-offset-model.
#
# WHAT IT PROVES, and why each check is needed.
#
# 1. PASSING THE MODEL'S OWN STORED STATISTICS EXPLICITLY MUST REPRODUCE THE DEFAULT EXACTLY.
#    This is the whole safety argument for the eleven arms that already exist. The new flag
#    reroutes the normalization call, and the same patch also moved the np.concatenate of the
#    predictions above the `del accs` that frees the offset model's input. If either change
#    perturbed anything, two runs that should be identical will not be. Bit-identity is the
#    right bar: same weights, same input, same order of operations.
#    NOTE the stats are PER FOLD (4.155 to 4.232 across the five p300 folds), so the test
#    reads the stored pair for the fold that actually owns the test chromosome rather than
#    assuming fold0.
#
# 2. A DIFFERENT MEAN MUST ACTUALLY CHANGE THE OUTPUT. A flag that is silently ignored would
#    pass check 1 perfectly, so the negative control is what makes check 1 meaningful.
#
# 3. THE OFFSET MUST BE ADDITIVE IN LOG SPACE. 4.1 paints expm1(logcounts)/width, so adding an
#    offset of L to the logcounts multiplies the painted value by exp(L). The check recomputes
#    the offset model's own prediction independently and confirms the residual arm's track
#    equals expm1(residual_logits + offset_logits) to within float tolerance, which is the one
#    piece of arithmetic that would silently produce a plausible-looking but wrong track.
#
# chr22 only: 3,520 regions, enough to compare tracks and small enough to run four times.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
PY=$D/.pixi/envs/multimodal/bin/python
REG=$A/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
SIZES=$A/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
M=$D/2026_0529_multimodal_p300_model/models/atac
ACC=$D/2026_0529_multimodal_p300_model/data/atac.bw
T=${SCRATCH:-/tmp}/accnorm_test.$$
mkdir -p "$T" "$P/log"
cd "$P"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "=== 1. default path ==="
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$M" \
    --mode multimodal --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --out-bw "$T/default.bw" | tee "$T/default.log"

FOLD=$(grep -o 'chr22: fold[0-9]' "$T/default.log" | grep -o '[0-9]$')
echo "chr22 is held out in fold${FOLD}"
read -r MEAN STD < <($PY -c "
import json; st=json.load(open('$M/fold${FOLD}/acc_normalization_stats.json'))
print(st['acc_mean'], st['acc_std'])")
echo "stored stats for fold${FOLD}: mean=$MEAN std=$STD"

echo "=== 2. same statistics passed explicitly (must be bit-identical) ==="
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$M" \
    --mode multimodal --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --acc-mean "$MEAN" --acc-std "$STD" --out-bw "$T/explicit.bw" >/dev/null

echo "=== 3. negative control: a different mean must move the output ==="
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$M" \
    --mode multimodal --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --acc-mean 3.7232 --acc-std 1.1428 --out-bw "$T/regions.bw" >/dev/null

echo "=== 4. offset arithmetic, on a residual model ==="
RM=$P/models/residual5p_multimodal_hw500_clw10
OM=$P/models/atac5p_hw500_clw10
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$RM" \
    --mode multimodal --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --count-offset-model "$OM" --out-bw "$T/resid_with.bw" >/dev/null
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$RM" \
    --mode multimodal --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --out-bw "$T/resid_without.bw" >/dev/null
$PY scripts/4.1.predict_h3k27ac_for_abc.py --regions "$REG" --model-dir "$OM" \
    --mode atac --accessibility-bw "$ACC" --chrom-sizes "$SIZES" \
    --chroms chr22 --no-rc-average --out-bw "$T/offset_only.bw" >/dev/null

$PY - "$T" <<'PYIN'
import sys, numpy as np, pyBigWig
T = sys.argv[1]

def vals(name):
    b = pyBigWig.open(f"{T}/{name}.bw")
    iv = b.intervals("chr22")
    b.close()
    iv = sorted(iv)
    return (np.array([x[0] for x in iv]), np.array([x[1] for x in iv]),
            np.array([x[2] for x in iv], dtype=np.float64))

s0, e0, d = vals("default")
s1, e1, x = vals("explicit")
_, _, r = vals("regions")
assert (s0 == s1).all() and (e0 == e1).all(), "interval sets differ"

fail = []
if not np.array_equal(d, x):
    n = int((d != x).sum())
    fail.append(f"CHECK 1 FAILED: explicit stored stats changed {n}/{len(d)} values, "
                f"max abs diff {np.abs(d - x).max():.3e}")
else:
    print(f"CHECK 1 PASS: {len(d):,} intervals bit-identical to the default path")

if np.array_equal(d, r):
    fail.append("CHECK 2 FAILED: region statistics produced an identical track, so "
                "--acc-mean is being ignored")
else:
    frac = float(np.mean(r / np.maximum(d, 1e-12)))
    print(f"CHECK 2 PASS: region statistics moved the track; mean ratio "
          f"region/default = {frac:.4f}, corr = {np.corrcoef(d, r)[0, 1]:.6f}")

sw, ew, w = vals("resid_with")
_, _, wo = vals("resid_without")
so, eo, o = vals("offset_only")
if not ((sw == so).all() and (ew == eo).all()):
    fail.append("CHECK 3 SKIPPED: interval sets differ between residual and offset runs")
else:
    width = (ew - sw).astype(float)
    # 4.1 paints expm1(lc).clip(min=0)/width, so logcounts are recoverable ONLY where the
    # clip did not fire. A residual model predicts observed MINUS accessibility, so its raw
    # logcounts are negative wherever the accessibility model already over-explains the
    # signal, and the no-offset track is floored at 0 there. Those intervals carry no
    # information about lc_res and cannot be used to verify the sum; excluding them is not
    # leniency, it is the only well-posed comparison. The production path is unaffected
    # because there the clip applies to residual+offset, which is the quantity that must be
    # non-negative.
    usable = (wo > 0) & (o > 0)
    lc_res = np.log1p(wo[usable] * width[usable])
    lc_off = np.log1p(o[usable] * width[usable])
    expect = np.expm1(lc_res + lc_off) / width[usable]
    got = w[usable]
    bad = ~np.isclose(got, expect, rtol=2e-4, atol=1e-8)
    print(f"  {int(usable.sum()):,} of {len(w):,} intervals have an unclipped residual "
          f"track; {int((~usable).sum()):,} were floored at 0 by the clip")
    if int(usable.sum()) < 100:
        fail.append("CHECK 3 INCONCLUSIVE: too few unclipped intervals to verify")
    elif bad.any():
        i = int(np.argmax(np.abs(got - expect)))
        fail.append(f"CHECK 3 FAILED: {int(bad.sum())}/{len(got)} unclipped intervals off; "
                    f"worst got {got[i]:.6g} expected {expect[i]:.6g}")
    else:
        print(f"CHECK 3 PASS: offset is additive in log space across "
              f"{int(usable.sum()):,} unclipped intervals; median painted value "
              f"{np.median(wo[usable]):.4g} -> {np.median(got):.4g} "
              f"(x{np.median(got) / max(np.median(wo[usable]), 1e-12):.2f} from the offset)")
    # Independent of the clip: the offset must RAISE the track everywhere, since the offset
    # model's logcounts are positive. This check needs no reconstruction at all.
    if (w < wo - 1e-9).any():
        fail.append(f"CHECK 4 FAILED: the offset lowered the track at "
                    f"{int((w < wo - 1e-9).sum())} intervals")
    else:
        # Report the ratio on the unclipped subset only. Over all intervals the median of
        # `wo` is 0, because most residuals are negative and get floored, so an overall
        # ratio divides by the 1e-12 guard and prints a meaningless 1e10.
        print(f"CHECK 4 PASS: the offset raised or held every interval; median ratio "
              f"x{np.median(got) / max(np.median(wo[usable]), 1e-12):.2f} on the "
              f"{int(usable.sum()):,} unclipped intervals")

print()
if fail:
    print("\n".join(fail))
    raise SystemExit("ACCNORM_OFFSET_TESTS_FAILED")
print("ACCNORM_OFFSET_TESTS_PASSED")
PYIN
rm -rf "$T"
