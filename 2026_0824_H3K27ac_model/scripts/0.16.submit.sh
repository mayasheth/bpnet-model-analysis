#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 8:00:00
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH -o log/telohaec_pe_compare.%j.txt
#SBATCH -e log/telohaec_pe_compare.%j.txt
#SBATCH --job-name=k27_telo_cmp
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

# condition -> element file stem (VEGF's element set is named VEGF_ctrl)
declare -A EL=( [ctrl]=TeloHAEC_ctrl [IL1b]=TeloHAEC_IL1b
                [TNFa]=TeloHAEC_TNFa [VEGF]=TeloHAEC_VEGF_ctrl )

for cond in ctrl IL1b TNFa VEGF; do
  echo "########## $cond ##########"
  $PY scripts/0.16.compare_pe_5prime_variants.py \
      --panel-dir "$P/data/panel" \
      --elements "$P/reference/celltype_elements/${EL[$cond]}_ATAC_candidate_elements.narrowPeak" \
      --condition "$cond" --outdir results --figdir figures
done
echo "ALL_DONE"
