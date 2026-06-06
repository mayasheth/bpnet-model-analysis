#!/bin/bash
# Log of commands for the multimodal p300 BPNet project
# Architecture: 5-channel input (DNA seq + accessibility), middle fusion
# See scripts/multimodal_bpnet.py for model definition
# See scripts/train_multimodal_bpnet.py for training script
#
# Environment: pixi run -e multimodal (see pixi.toml)
#   module load devel pixi/0.53.0

# ============================================================
# Stage 0: Preprocess accessibility BigWigs
# ============================================================
# Signal BigWigs (stranded, n_outputs=2) already exist from retrained p300 model:
#   plus:  2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
#   minus: 2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig

CHROM_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes

# ATAC: tagAlign files (Tn5-shift corrected) -> merged BigWig
# Files: ENCFF077FBI, ENCFF128WZG, ENCFF534DCE (K562 ATAC-seq, ENCODE)
module load devel pixi/0.53.0
pixi run -e multimodal bash scripts/0.1.make_accessibility_bigwig.sh \
    --input \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF077FBI.tn5.sorted.tagAlign.gz \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF128WZG.tn5.sorted.tagAlign.gz \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF534DCE.tn5.sorted.tagAlign.gz \
    --output 2026_0529_multimodal_p300_model/data/atac.bw \
    --chrom-sizes $CHROM_SIZES \
    --type atac

# DNase: BAM files -> raw-count BigWig (training normalizes accessibility)
# Files: ENCFF205FNC, ENCFF860XAE (K562 DNase-seq, ENCODE)
pixi run -e multimodal bash scripts/0.1.make_accessibility_bigwig.sh \
    --input \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF205FNC.filtered.sorted.bam \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF860XAE.filtered.sorted.bam \
    --output 2026_0529_multimodal_p300_model/data/dnase.bw \
    --chrom-sizes $CHROM_SIZES \
    --type dnase

# ============================================================
# Stage 1: Train model (ATAC and DNase variants, 5 folds each)
# ============================================================
cd $OAK/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model

for FOLD in 0 1 2 3 4; do
    sbatch scripts/1.1.submit_training_atac.sh $FOLD
    sbatch scripts/1.2.submit_training_dnase.sh $FOLD
done

# ============================================================
# Stage 2: SHAP attributions
# ============================================================
# sbatch scripts/2.1.submit_shap.sh atac fold0
