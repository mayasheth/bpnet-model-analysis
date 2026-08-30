#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH -o log/telohaec_coupling.%j.txt
#SBATCH -e log/telohaec_coupling.%j.txt
#SBATCH --job-name=k27_telo_coup
#
# Model-free ATAC vs H3K27ac coupling for the four TeloHAEC conditions.
#
# Run BEFORE any TeloHAEC model. Cross-cell-type transfer numbers are uninterpretable
# without it: the K562->GM12878 ATAC-only "transfer drop" turned out to be entirely
# explained by weaker ATAC-H3K27ac coupling in GM12878 (0.409 vs 0.510), with no model
# degradation at all. Without the model-free reference we would have read that as a
# generalisation failure.
#
# Uses the r1 5' target per the 2026-08-29 decision, and the merged (not per-replicate)
# tracks. Appends to results/atac_vs_h3k27ac_by_celltype.tsv alongside K562 and GM12878.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

declare -A EL=( [ctrl]=TeloHAEC_ctrl [IL1b]=TeloHAEC_IL1b
                [TNFa]=TeloHAEC_TNFa [VEGF]=TeloHAEC_VEGF_ctrl )

for cond in ctrl IL1b TNFa VEGF; do
  echo "########## TeloHAEC_$cond ##########"
  $PY scripts/0.12.atac_vs_h3k27ac.py \
      --label "TeloHAEC_${cond}" \
      --elements "$P/reference/celltype_elements/${EL[$cond]}_ATAC_candidate_elements.narrowPeak" \
      --atac-bw "$P/data/panel/TeloHAEC_${cond}_atac.bw" \
      --h3k27ac-bw "$P/data/panel/TeloHAEC_${cond}_h3k27ac_r1_5p_plus.bw,$P/data/panel/TeloHAEC_${cond}_h3k27ac_r1_5p_minus.bw" \
      --half-window 500 --outdir results --figdir figures
done
echo "ALL_DONE"
