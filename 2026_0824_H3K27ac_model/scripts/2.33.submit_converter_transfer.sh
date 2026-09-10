#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 4:00:00
#SBATCH --mem=96G
#SBATCH -o log/conv_tx.%j.txt
#SBATCH -e log/conv_tx.%j.txt
#SBATCH --job-name=conv_tx
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# THE DECISIVE TEST: the K562-trained converter applied to GM12878.
#
# The fold-0 pilot reached 96.2% of the K562 inter-replicate ceiling at 1 bp on the top
# quintile. That is in-cell, and in-cell is exactly where this project has been fooled
# before -- fragment channels were "real in-cell, exactly null transferred both
# directions", and the whole reason the DNase input result mattered was that it was the
# one change that survived transfer. A converter that only works in the cell type it was
# trained on is worth nothing: the point is to run where DNase does not exist.
#
# The ATAC-only arm travels too, as the control. If the sequence margin over ATAC collapses
# on transfer while ATAC-only holds up, the converter is reading K562-specific sequence
# associations rather than a general ATAC->DNase mapping.
#
# GM12878's own DNase ceiling is lower (1 bp top-quintile 0.686 against K562's 0.848), so
# read the fraction of ceiling, not the raw correlation, across cell types.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

run () {  # label, model_dir, mode
    echo "===== $1 -> GM12878 ====="
    $PY scripts/2.31.converter_profile_eval.py \
        --model-dir "$P/models/$2" --fold 0 --label "$1" --cell gm12878 --mode "$3"
    echo
}
run seq_atac_clw10_tx  multimodal5p_conv_atac2dnase_atacel_hw500_clw10 multimodal
run atac_only_clw10_tx atac5p_conv_atac2dnase_atacel_hw500_clw10       atac
echo done
