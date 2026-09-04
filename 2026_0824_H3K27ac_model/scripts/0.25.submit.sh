#!/bin/bash
#SBATCH -p engreitz,normal,owners
#
# PARTITION: `engreitz` first -- the lab-owned partition, 9 nodes at 24+ cores and 192 GB+
# with a 7-day limit and no GPUs, so it is the right home for CPU work and does not compete
# with the general GPU queues. `normal` and `owners` follow as fallbacks; SLURM starts the
# job wherever a slot frees first.
#SBATCH -t 2:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/profceil.%j.txt
#SBATCH -e log/profceil.%j.txt
#SBATCH --job-name=profceil
#
# Profile-shape ceiling vs bin size, K562 and GM12878. Data only, no GPU, no training.
# Gates the binned-profile-target architecture change: if the 1 bp ceiling is near zero and
# rises with bin size, the profile head is currently fitting noise and the saturation point
# is the resolution worth predicting.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

echo "########## K562 ##########"
$PY scripts/0.25.profile_ceiling_by_binsize.py --label k562 \
    --rep1-plus data/h3k27ac_rep1_5p_plus.bw --rep1-minus data/h3k27ac_rep1_5p_minus.bw \
    --rep2-plus data/h3k27ac_rep2_5p_plus.bw --rep2-minus data/h3k27ac_rep2_5p_minus.bw \
    --elements $D/reference/K562_DNase_candidate_elements.narrowPeak

echo
echo "########## GM12878 ##########"
$PY scripts/0.25.profile_ceiling_by_binsize.py --label gm12878 \
    --rep1-plus data/gm12878_h3k27ac_rep1_5p_plus.bw \
    --rep1-minus data/gm12878_h3k27ac_rep1_5p_minus.bw \
    --rep2-plus data/gm12878_h3k27ac_rep2_5p_plus.bw \
    --rep2-minus data/gm12878_h3k27ac_rep2_5p_minus.bw \
    --elements $D/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
