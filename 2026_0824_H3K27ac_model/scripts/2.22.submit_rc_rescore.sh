#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 16:00:00
#SBATCH --mem=96G
#SBATCH -o log/rc_rescore.%j.txt
#SBATCH -e log/rc_rescore.%j.txt
#SBATCH --job-name=k27_rc_rescore
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Re-score every model-based report table with RC averaging on, IN ONE BATCH.
#
# RC averaging became the 2.15 default on 2026-09-03. It shifts every number by +0.002 to
# +0.016, so a report containing both single-pass and RC-averaged tables would be internally
# inconsistent in a way no reader could detect. This job produces the full RC-averaged set so
# the switch happens once.
#
# Outputs are written under an `rc_` prefix rather than overwriting, so the single-pass
# tables survive for comparison and the switch stays reversible.
#
# NOT re-scored here, because no model is involved and RC is meaningless for them:
# the coupling panel, the fragment-length distribution, both ceilings.

set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
TR=$D/2026_0606_GM12878_transferability
PY=$D/.pixi/envs/multimodal/bin/python
GMEL=$TR/reference/GM12878_candidate_elements.narrowPeak
K5EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
P300EL=$D/reference/K562_DNase_candidate_elements.narrowPeak
cd "$P"

run () {  # config out_prefix elements [pair args...]
  local cfg="$1" pfx="$2" el="$3"; shift 3
  echo "########## $pfx ##########"
  $PY scripts/2.15.perfold_from_config.py "config/$cfg" "$pfx" "$el" "$@"
  echo
}

# Fig 1 headline, and Fig 5 residual grids
run rc_fiveprime_configs.json           rc_fiveprime_            "$K5EL" \
    --pair multimodal5p sequence5p
run rc_residual_grid_k562_configs.json  rc_residual_grid_k562_   "$K5EL" \
    --pair residual_sequence sequence5p --pair residual_multimodal multimodal5p
run rc_residual_grid_gm12878_configs.json rc_residual_grid_gm12878_ "$GMEL" \
    --pair residual_sequence sequence5p --pair residual_multimodal multimodal5p

# Fig 6 accessibility input definition
run rc_accs5p_k562_configs.json         rc_accs5p_k562_          "$K5EL" \
    --pair atac_5PRIME atac_OLD --pair multimodal_5PRIME multimodal_OLD
run rc_accs5p_gm12878_configs.json      rc_accs5p_gm12878_       "$GMEL" \
    --pair atac_5PRIME atac_OLD --pair multimodal_5PRIME multimodal_OLD
run rc_accs5p_p300_configs.json         rc_accs5p_p300_          "$P300EL" \
    --pair atac_5PRIME atac_OLD --pair multimodal_5PRIME multimodal_OLD

# Fig 7 transfer and deployment
run transfer_k562_to_gm_configs.json    rc_transfer_k562_to_gm_  "$GMEL" \
    --pair K562resMM_to_GM K562mm_to_GM --pair K562resSeq_to_GM K562mm_to_GM
run deploy_k562_to_gm_configs.json      rc_deploy_k562_to_gm_    "$GMEL" \
    --pair K562resMM_to_GM K562mm_to_GM --pair K562resSeq_to_GM K562mm_to_GM
run transfer_gm_to_k562_configs.json    rc_transfer_gm_to_k562_  "$K5EL" \
    --pair GMresMM_to_K562 GMmm_to_K562 --pair GMresSeq_to_K562 GMmm_to_K562
run deploy_gm_to_k562_configs.json      rc_deploy_gm_to_k562_    "$K5EL" \
    --pair GMresMM_to_K562 GMmm_to_K562 --pair GMresSeq_to_K562 GMmm_to_K562

# Figs 10 and 12 architecture
run wide_k562_configs.json              rc_wide_k562_            "$K5EL" \
    --pair sequence_WIDE sequence_narrow --pair multimodal_WIDE multimodal_narrow
run wide_gm12878_configs.json           rc_wide_gm12878_         "$GMEL" \
    --pair sequence_WIDE sequence_narrow --pair multimodal_WIDE multimodal_narrow
run fragchan_k562_configs.json          rc_fragchan_k562_        "$K5EL" \
    --pair multimodal_FRAGCHAN multimodal_flat

echo "All tables re-scored with RC averaging. Re-run 3.7 and every 3.x figure script "
echo "against the rc_ prefixes before touching the report."
