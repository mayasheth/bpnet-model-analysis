#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/slopeprof.%j.txt
#SBATCH -e log/slopeprof.%j.txt
#SBATCH --job-name=slope_prof
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
cd $D/2026_0824_H3K27ac_model
P=data/abc_predicted
$D/.pixi/envs/multimodal/bin/python scripts/4.28.slope_profile.py \
  --regions $A/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
  --accessibility-bw $D/2026_0529_multimodal_p300_model/data/atac.bw \
  --track in_cell=$P/predp300_multimodal.bw \
  --track transferred=$P/predp300gm_gm12878_multimodal.bw \
  --track in_cell_an=$P/predp300an_k562_multimodal.bw \
  --track transferred_an=$P/predp300an_gm12878_multimodal.bw \
  --track h3k27ac_in_cell=$P/predk27ac_k562_multimodal.bw \
  --track h3k27ac_transferred=$P/predk27ac_gm12878_multimodal.bw \
  --track h3k27ac_resid_in_cell=$P/predk27ac_k562_residual.bw \
  --track h3k27ac_resid_transferred=$P/predk27ac_gm12878_residual.bw \
  --ratio transferred/in_cell --ratio transferred_an/in_cell_an \
  --ratio h3k27ac_transferred/h3k27ac_in_cell \
  --ratio h3k27ac_resid_transferred/h3k27ac_resid_in_cell
