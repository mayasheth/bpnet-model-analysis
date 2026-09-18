#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH -o log/testmulticell.%j.txt
#SBATCH -e log/testmulticell.%j.txt
#SBATCH --job-name=test_multicell
#
# Regression gate for the 2026-09-18 multi-cell-type additions to
# train_multimodal_bpnet.py: --cell-types-json and --holdout-cell.
#
# THE CHANGE TOUCHED THE SINGLE-CELL PATH, which is what makes a gate necessary rather
# than nice to have. Extraction, the negative cap, the profile-target checks and
# accessibility normalization all moved inside a per-cell-type loop, and the provenance
# dump moved after panel resolution. Every existing model in this project was trained by
# the code being edited.
#
# CHECKS
#   1. OLD vs NEW on the same single-cell-type arguments. Runs the pre-change script from
#      git alongside the current one and compares window counts and the accessibility
#      normalization statistics. Those are deterministic given the inputs, unlike the
#      weights, which cuDNN does not reproduce bit-for-bit across nodes (decision
#      2026-09-03). If the data path moved, these numbers move.
#   2. A ONE-CELL PANEL equals the single-cell path. Same comparison, but driving the new
#      code through --cell-types-json. This is what proves the panel machinery is a
#      generalisation and not a second implementation.
#   3. A TWO-CELL PANEL pools. Window counts must be the sum of the parts, and each cell
#      type's normalization statistics must be its OWN, not a shared pooled value; a
#      single mean and std across cell types is the failure this design exists to avoid.
#   4. --holdout-cell drops exactly the named cell type and errors on an unknown name.
#
# One epoch, one fold. Nothing here depends on the model learning anything.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
T=${SCRATCH:-/tmp}/multicell_test.$$
mkdir -p "$T" "$P/log"
cd "$D"

git show HEAD:scripts/train_multimodal_bpnet.py > "$T/train_OLD.py"
# The old copy lives outside scripts/, and python only puts the SCRIPT's own directory on
# sys.path, so it cannot find its sibling module multimodal_bpnet. Both runs get the real
# scripts/ on PYTHONPATH so they import the identical model definition; that module is not
# what is being tested here.
export PYTHONPATH="$D/scripts${PYTHONPATH:+:$PYTHONPATH}"
echo "old script: $(wc -l < "$T/train_OLD.py") lines; new: $(wc -l < scripts/train_multimodal_bpnet.py) lines"

K_PEAKS=$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak
K_SIG_P=$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
K_SIG_M=$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig
K_ACC=$D/2026_0529_multimodal_p300_model/data/atac.bw
G_PEAKS=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz
G_SIG_P=$D/2026_0606_GM12878_transferability/data/ENCFF960OFK_plus.bw
G_SIG_M=$D/2026_0606_GM12878_transferability/data/ENCFF941MGK_minus.bw
G_ACC=$D/2026_0606_GM12878_transferability/data/atac.bw
NEG=$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
FOLD=$D/reference/hg38_five_folds.json
echo "fold file: $FOLD"

COMMON=(--negatives "$NEG" --genome "$GEN" --fold "$FOLD" --fold-key 0
        --mode multimodal --max-epochs 1 --early-stopping 1
        --max-negatives 2000 --batch-size 32 --device cpu --n-workers 2)

cat > "$T/panel_one.json" <<JSON
{"K562": {"peaks": "$K_PEAKS", "signal_plus_bw": "$K_SIG_P",
          "signal_minus_bw": "$K_SIG_M", "accessibility_bw": "$K_ACC"}}
JSON
cat > "$T/panel_two.json" <<JSON
{"K562": {"peaks": "$K_PEAKS", "signal_plus_bw": "$K_SIG_P",
          "signal_minus_bw": "$K_SIG_M", "accessibility_bw": "$K_ACC"},
 "GM12878": {"peaks": "$G_PEAKS", "signal_plus_bw": "$G_SIG_P",
             "signal_minus_bw": "$G_SIG_M", "accessibility_bw": "$G_ACC"}}
JSON

run () {  # run <label> <script> [extra args...]
    local label=$1 script=$2; shift 2
    echo "=== running $label ==="
    $PY "$script" --output-dir "$T/$label" \
        --peaks "$K_PEAKS" --signal-plus-bw "$K_SIG_P" --signal-minus-bw "$K_SIG_M" \
        --accessibility-bw "$K_ACC" "${COMMON[@]}" "$@" > "$T/$label.log" 2>&1 \
        || { echo "FAILED: $label"; tail -25 "$T/$label.log"; return 1; }
    grep -E 'Extracted|normalization:' "$T/$label.log" | sed 's/^/    /'
}

run old "$T/train_OLD.py"
run new_single scripts/train_multimodal_bpnet.py
run new_panel1 scripts/train_multimodal_bpnet.py --cell-types-json "$T/panel_one.json"
run new_panel2 scripts/train_multimodal_bpnet.py --cell-types-json "$T/panel_two.json"
run new_hold   scripts/train_multimodal_bpnet.py --cell-types-json "$T/panel_two.json" \
                                                 --holdout-cell GM12878

echo "=== unknown --holdout-cell must fail ==="
if $PY scripts/train_multimodal_bpnet.py --output-dir "$T/bad" \
      --peaks "$K_PEAKS" --signal-plus-bw "$K_SIG_P" --accessibility-bw "$K_ACC" \
      "${COMMON[@]}" --cell-types-json "$T/panel_two.json" --holdout-cell NOPE \
      > "$T/bad.log" 2>&1; then
    echo "  CHECK 4a FAILED: an unknown holdout name was accepted"; BAD4A=1
else
    grep -q 'is not in' "$T/bad.log" && echo "  CHECK 4a PASS: rejected with a useful message" \
        || { echo "  CHECK 4a FAILED: rejected for the wrong reason"; tail -3 "$T/bad.log"; BAD4A=1; }
fi

$PY - "$T" <<'PYIN'
import json, re, sys, os
T = sys.argv[1]
def counts(lbl):
    """Windows actually used, as (train peaks, negatives, val peaks).

    A panel logs one 'Extracted N' per cell type per region set, so a two-cell panel emits
    six numbers and not three. Sum them per position rather than reading them positionally,
    which is what the first version of this test got wrong: it compared a six-element list
    against a three-element one and called a correctly pooled panel a failure.
    """
    txt = open(f"{T}/{lbl}.log").read()
    n = [int(x) for x in re.findall(r"Extracted (\d+)", txt)]
    assert len(n) % 3 == 0, f"{lbl}: {len(n)} extraction lines is not a multiple of 3"
    return [sum(n[i::3]) for i in range(3)]
def stats(lbl):
    p = f"{T}/{lbl}/acc_normalization_stats.json"
    return json.load(open(p)) if os.path.exists(p) else None

old, new_s = counts("old"), counts("new_single")
p1, p2, hold = counts("new_panel1"), counts("new_panel2"), counts("new_hold")
so, sn, s1, s2 = stats("old"), stats("new_single"), stats("new_panel1"), stats("new_panel2")
fail = []

print(f"\nwindow counts (train peaks, negatives, val peaks)")
for lbl, c in (("old", old), ("new_single", new_s), ("new_panel1", p1),
               ("new_panel2", p2), ("new_hold", hold)):
    print(f"  {lbl:<12} {c}")

if old != new_s:
    fail.append(f"CHECK 1 FAILED: single-cell window counts moved, {old} -> {new_s}")
elif so != sn:
    fail.append(f"CHECK 1 FAILED: normalization stats moved, {so} -> {sn}")
else:
    print(f"\nCHECK 1 PASS: old and new agree on the single-cell path, counts {old}, "
          f"stats {sn}")

# A one-cell panel must match the single-cell path in BOTH window counts and the shape
# and values of the statistics file, since its models are consumed by the same 4.1.
if p1 != new_s or s1 != sn:
    fail.append(f"CHECK 2 FAILED: a one-cell panel differs from the single-cell path; "
                f"counts {p1} vs {new_s}, stats {s1} vs {sn}")
else:
    print("CHECK 2 PASS: a one-cell panel reproduces the single-cell path exactly")

# The pooled set must be strictly larger than either cell type alone on peaks and
# validation. Negatives are capped PER cell type, so they scale with the panel too.
if all(p2[i] > new_s[i] for i in range(3)):
    print(f"CHECK 3a PASS: the two-cell panel pools, {new_s} -> {p2}")
else:
    fail.append(f"CHECK 3a FAILED: two-cell panel did not pool: {new_s} -> {p2}")
if isinstance(s2, dict) and "per_cell_type" in s2:
    pc = s2["per_cell_type"]
    if len(pc) == 2 and pc.get("K562") != pc.get("GM12878"):
        print(f"CHECK 3b PASS: per-cell-type statistics kept separate: "
              + ", ".join(f"{k} mean={v['acc_mean']:.4f}" for k, v in pc.items()))
    else:
        fail.append(f"CHECK 3b FAILED: expected two distinct per-cell stats, got {pc}")
else:
    fail.append(f"CHECK 3b FAILED: multi-cell stats file has the wrong shape: {s2}")

if hold == new_s:
    print("CHECK 4 PASS: holding out GM12878 leaves exactly the K562 single-cell set")
else:
    fail.append(f"CHECK 4 FAILED: holdout gave {hold}, expected the K562 set {new_s}")

tr = json.load(open(f"{T}/new_panel2/training_target.json"))
if sorted(tr.get("trained_on", [])) == ["GM12878", "K562"] and "panel" in tr:
    print("CHECK 5 PASS: training_target.json names the cell types and their tracks")
else:
    fail.append(f"CHECK 5 FAILED: provenance missing; trained_on={tr.get('trained_on')}")

print()
if fail:
    print("\n".join(fail)); raise SystemExit("MULTICELL_TESTS_FAILED")
print("MULTICELL_TESTS_PASSED")
PYIN
echo "artifacts in $T"
