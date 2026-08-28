#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH -o log/regen_figs.%j.txt
#SBATCH -e log/regen_figs.%j.txt
#SBATCH --job-name=k27_regen_figs
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

echo "### profile comparison (Fig 2) ###"
$PY scripts/0.6.compare_profiles.py \
  --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" \
  --track "ATAC:$D/2026_0529_multimodal_p300_model/data/atac.bw" \
  --track "H3K27ac:$P/data/h3k27ac_5p_plus.bw,$P/data/h3k27ac_5p_minus.bw" \
  --stratify-by "H3K27ac" --outdir results --figdir figures

echo "### K562 ceiling by window (Fig 3) ###"
$PY scripts/0.3.replicate_ceiling_by_window.py \
  --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" \
  --rep1-bw "$P/data/h3k27ac_rep1_5p_plus.bw,$P/data/h3k27ac_rep1_5p_minus.bw" \
  --rep2-bw "$P/data/h3k27ac_rep2_5p_plus.bw,$P/data/h3k27ac_rep2_5p_minus.bw" \
  --outdir results --figdir figures --label fiveprime_

echo "### ATAC fragment lengths (Fig 7) ###"
$PY scripts/0.10.plot_fragment_lengths.py \
  --bams "$SCRATCH/atac_pe/ENCFF077FBI.pe.bam" "$SCRATCH/atac_pe/ENCFF128WZG.pe.bam" \
         "$SCRATCH/atac_pe/ENCFF534DCE.pe.bam" \
  --outdir results --figdir figures

echo "### ATAC vs H3K27ac coupling, both cell types (new Fig 8) ###"
$PY scripts/0.12.atac_vs_h3k27ac.py --label K562 \
  --elements "$D/reference/K562_DNase_candidate_elements.narrowPeak" \
  --atac-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw" \
  --h3k27ac-bw "$P/data/h3k27ac_5p_plus.bw,$P/data/h3k27ac_5p_minus.bw" \
  --outdir results --figdir figures
$PY scripts/0.12.atac_vs_h3k27ac.py --label GM12878 \
  --elements "$D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak" \
  --atac-bw "$D/2026_0606_GM12878_transferability/data/atac.bw" \
  --h3k27ac-bw "$P/data/gm12878_h3k27ac_5p_plus.bw,$P/data/gm12878_h3k27ac_5p_minus.bw" \
  --outdir results --figdir figures
echo ALL_DONE
