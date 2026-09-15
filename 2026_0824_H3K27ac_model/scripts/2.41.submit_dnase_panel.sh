#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -o log/dnasepanel.%j.txt
#SBATCH -e log/dnasepanel.%j.txt
#SBATCH --job-name=dnase_panel
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Three-way DNase-input transfer panel: K562, GM12878, THP-1.
#
# THE QUESTION TWO FINDINGS ARE STUCK ON. F-009 (p300) and F-014 (DNase inputs) both report
# the same asymmetry: a K562-trained model keeps its advantage in GM12878, a GM12878-trained
# model loses it in K562. With two cell types that is unattributable, because "K562 is a good
# TRAINING cell type" and "K562->GM12878 is a good PAIR" predict identical numbers. A third
# target separates them:
#   * if K562-trained models travel well to BOTH GM12878 and THP-1, K562 is the good trainer
#   * if K562->GM12878 works but K562->THP-1 does not, it was the pair
#   * if nothing transfers into K562 from either source, the asymmetry is about K562 being a
#     hard TARGET rather than about any source being weak
#
# ONE PASS PER TARGET CELL TYPE, because the element set and the observed mark differ per
# target and 2.15 scores one element set at a time. Each pass has that target's own
# DNase-only model as its baseline, so residual_pearson and incremental_r2 measure what
# sequence adds beyond locally available accessibility. Configs come from 2.40, which
# applies the source/target rule mechanically: model from the SOURCE, accessibility and
# observed mark and elements from the TARGET.
#
# THP-1 CARRIES A CAVEAT THE OTHER TWO DO NOT: one DNase replicate, so no DNase shape
# ceiling, and its H3K27ac is DEEPER than K562's (mean 208.2 per 1 kb window against 51.1)
# while its DNase is SHALLOWER (281.4 against 910.2). Depth differences cut both ways, so
# neither cell type is straightforwardly the easier target and a transfer drop into THP-1 is
# not automatically evidence about the source model.
#
# --no-rc-average keeps these comparable to every other table in the project.
# Env: DEPTH_MATCHED=1 scores the depth-matched variant instead, where all three DNase
#      libraries are thinned to GM12878's 50.3M reads (0.41/0.42). Element sets are
#      unchanged, so the two panels are directly comparable region for region.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

# Regenerate with --check so a preempted fold stops this before it burns GPU time. 2.40
# refuses to write unless all five folds of all twelve arms left a completion marker.
DM_ARG=(); SUFFIX=""; PREFIX="dnasepanel_"
if [[ -n "${DEPTH_MATCHED:-}" ]]; then
    DM_ARG=(--depth-matched); SUFFIX="_depthmatched"; PREFIX="dnasepaneldm_"
    echo "DEPTH-MATCHED VARIANT: all three DNase inputs thinned to 50.3M reads"
fi
$PY scripts/2.40.make_dnase_panel_configs.py "${DM_ARG[@]}" --check

for target in k562 gm12878 thp1; do
    CFG=config/dnase_panel${SUFFIX}_on_${target}_configs.json
    EL=$($PY -c "import json;print(json.load(open('$CFG'))['_elements'])")
    PAIRS=()
    for lbl in $($PY -c "import json;print(' '.join(c['label'] for c in json.load(open('$CFG'))['compare']))"); do
        PAIRS+=(--pair "$lbl" dnase_only_LOCAL)
    done
    # Also state the transfer question directly: each transferred arm against the local one.
    LOCAL=$($PY -c "import json;print(json.load(open('$CFG'))['compare'][0]['label'])")
    for lbl in $($PY -c "import json;print(' '.join(c['label'] for c in json.load(open('$CFG'))['compare'][1:]))"); do
        PAIRS+=(--pair "$lbl" "$LOCAL")
    done
    echo "########## DNase panel scored on ${target} ##########"
    echo "elements: $EL"
    $PY scripts/2.15.perfold_from_config.py "$CFG" "${PREFIX}${target}_" "$EL" \
        --no-rc-average "${PAIRS[@]}"
done
echo DNASE_PANEL_DONE
