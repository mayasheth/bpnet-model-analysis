#!/bin/bash
# GM12878 cross-cell-type transferability analysis
# Goal: evaluate whether K562-trained p300 models generalize to GM12878
# Models tested:
#   (1) GM12878 EP300 BPNet (ENCODE/TFAtlas) — in-cell-type ceiling
#   (2) K562 p300 BPNet v1 — sequence-only (K562-trained)
#   (3) K562 multimodal ATAC BPNet — DNA + ATAC accessibility (K562-trained)
# See HANDOVER.md for full status and context.
# Date: 2026-06-06

# --- FILE PATHS --- #

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
SCRIPTS_DIR=$PROJECT_DIR/scripts
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
REF_DIR=$THIS_DIR/reference
DATA_DIR=$THIS_DIR/data
LOG_DIR=$THIS_DIR/log
PRED_DIR=$THIS_DIR/predictions

# GM12878 candidate elements (154,224 regions; derived from H3K27ac megamap)
GM12878_ELEMENTS=$REF_DIR/GM12878_candidate_elements.narrowPeak  # DONE

# GM12878 p300 peaks (21,068 peaks)
GM12878_PEAKS=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz

# GM12878 EP300 signal BigWigs (5' end, merged replicates) — output of stage 0.3
EP300_PLUS=$DATA_DIR/EP300_plus.bw
EP300_MINUS=$DATA_DIR/EP300_minus.bw

# ENCODE processed signal BigWigs (fold-change over control) — for sanity checking only
ENCODE_PLUS=$DATA_DIR/ENCFF960OFK_plus.bw
ENCODE_MINUS=$DATA_DIR/ENCFF941MGK_minus.bw

# GM12878 ATAC BigWig (for multimodal model) — generated in stage 0.2
ATAC_BW=$DATA_DIR/atac.bw

# Reference genome
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa

# K562-trained models
V1_MODEL_PATTERN="$PROJECT_DIR/2025_0517_official_EP300_K562_model/models/release_run_1/fold{fold}/ENCSR000EGE/ENCSR000EGE_split000"
MULTIMODAL_MODEL_DIR=$PROJECT_DIR/2026_0529_multimodal_p300_model

# GM12878 BPNet model (in-cell-type ceiling)
# Experiment accession: ENCSR000DZG; extracted from ENCFF778VIB.tar.gz
# TF SavedModel format; note fold_ (underscore) and accession in path vs K562 fold{N}/model_split000
GM12878_MODEL_DIR=$THIS_DIR/GM12878_EP300_BPNet
GM12878_MODEL_PATTERN="$GM12878_MODEL_DIR/models/fold_{fold}/ENCSR000DZG_split000"


# ============================================================
# STAGE 0: PREPARE INPUT DATA
# ============================================================

## [0.1] GM12878 candidate elements narrowPeak — DONE 2026-06-06
# 154,224 regions; source: .../GM12878_H3K27ac_megamap/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
# bash $SCRIPTS_DIR/0.2.bed_to_narrowPeak.sh \
#   .../macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
#   $GM12878_ELEMENTS 1057


## [0.2] Generate GM12878 ATAC BigWig from tagAlign files
# (once tagAlign files are available at $OAK/Users/sheth/Data/ENCODE/GM12878/ATAC)
# sbatch --partition=owners --time=8:00:00 --mem=64G --cpus-per-task=4 \
#   --job-name=gm_atac_bw --output=$LOG_DIR/atac_bw.%j.txt --error=$LOG_DIR/atac_bw.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal bash \
#     $PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/0.1.make_accessibility_bigwig.sh \
#     --input <tagAlign1> <tagAlign2> ... \
#     --output $ATAC_BW --chrom-sizes $CHR_SIZES --type atac"


## [0.3] Generate GM12878 EP300 BigWigs (merged reps, 5' end stranded) — job 27925774
# Input BAMs: ENCFF515HYM.filtered.sorted.bam (14.4M reads), ENCFF215GSQ.filtered.sorted.bam (15.6M reads)
# Location: $OAK/Users/sheth/Data/ENCODE/GM12878/EP300/
# Outputs: $EP300_PLUS, $EP300_MINUS
# Note: job 27925622 failed (missing samtools index before view); fixed in 27925774


## [0.4] Download ENCODE BPNet signal BigWigs — job 27925868
# Plus:  ENCFF960OFK → $ENCODE_PLUS
# Minus: ENCFF941MGK → $ENCODE_MINUS
# These are the actual signal BigWigs from the ENCODE BPNet annotation object —
# the same files used to train the GM12878 BPNet model (not a standard fold-change track).
# Purpose: sanity-check against our BAM-derived BigWigs; could also be used directly as model input.


## [0.5] Download and extract GM12878 EP300 BPNet model — DONE (job 27925831)
# Source: https://www.encodeproject.org/files/ENCFF778VIB/@@download/ENCFF778VIB.tar.gz
# Extracted to $GM12878_MODEL_DIR/models/fold_{0-4}/ENCSR000DZG_split000/ (TF SavedModel)
# Tars extracted in-place: for fold in 0 1 2 3 4; do
#   tar -xf $GM12878_MODEL_DIR/models/fold_${fold}/model.fold_${fold}.ENCSR000DZG.tar \
#     -C $GM12878_MODEL_DIR/models/fold_${fold}/; done
# Model path pattern: $GM12878_MODEL_PATTERN (uses fold_ underscore; no script changes needed)


# ============================================================
# STAGE 1: PREDICT WITH GM12878 BPNET (IN-CELL-TYPE CEILING)
# ============================================================
# Model: $GM12878_MODEL_DIR (5-fold, TF SavedModel format — same arch as K562 v1)
# Loci: GM12878 candidate elements
# Requires: GM12878 EP300 signal BigWigs + config JSON

GM12878_V1_PRED_DIR=$PRED_DIR/gm12878_bpnet
GM12878_V1_CONFIG=$THIS_DIR/config/input_data_predict.gm12878_bpnet.json

# SUBMITTED — job 27974981 (array 0-4), log: $LOG_DIR/gm12878_pred_f*.27974981.txt
# Note: uses K562 controls as bias proxy (GM12878 controls not downloaded; minor effect on counts Pearson)
# sbatch \
#   --array=0-4 \
#   --output=$LOG_DIR/gm12878_pred_f%a.%A.txt \
#   --error=$LOG_DIR/gm12878_pred_f%a.%A.txt \
#   --job-name=gm_bpnet_pred \
#   $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
#     "$GM12878_MODEL_PATTERN" \
#     $CHR_SIZES $GENOME \
#     $GM12878_V1_PRED_DIR/mean \
#     $GM12878_V1_CONFIG


# ============================================================
# STAGE 2: PREDICT WITH K562 p300 BPNET V1 (SEQUENCE-ONLY)
# ============================================================

K562_V1_PRED_DIR=$PRED_DIR/k562_bpnet_v1
K562_V1_CONFIG=$THIS_DIR/config/input_data_predict.k562_v1.json

## [2.1] Config created: $K562_V1_CONFIG
# Signal: $EP300_PLUS / $EP300_MINUS (GM12878 BAM-derived 5' end BigWigs)
# Loci: $GM12878_ELEMENTS
# Bias: K562 controls (ENCSR000EGE_control_plus/minus.bigWig) — correct for this K562-trained model

## [2.2] Submit mean predictions — SUBMITTED job 28041042 (array 0-4)
# log: $LOG_DIR/k562_v1_pred_f*.28041042.txt
# Note: first attempt (27974979) failed — model path was missing /ENCSR000EGE_split000
# sbatch \
#   --array=0-4 \
#   --output=$LOG_DIR/k562_v1_pred_f%a.%A.txt \
#   --error=$LOG_DIR/k562_v1_pred_f%a.%A.txt \
#   --job-name=gm_v1_pred \
#   $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
#     "$V1_MODEL_PATTERN" \
#     $CHR_SIZES $GENOME \
#     $K562_V1_PRED_DIR/mean \
#     $K562_V1_CONFIG


# ============================================================
# STAGE 3: PREDICT WITH K562 MULTIMODAL ATAC BPNET
# ============================================================

MULTIMODAL_PRED_DIR=$PRED_DIR/k562_multimodal_atac

## [3.1] Submit multimodal predictions — SUBMITTED job 28041421
# Note: 2.1.submit_predict_multimodal.sh has K562 paths hardcoded; submitted directly to predict_multimodal.py
# log: $LOG_DIR/multimodal_pred.28041421.txt
# sbatch --partition=owners,engreitz,normal \
#   --time=6:00:00 --mem=64G --cpus-per-task=4 \
#   --job-name=gm_mm_pred \
#   --output=$LOG_DIR/multimodal_pred.%j.txt \
#   --error=$LOG_DIR/multimodal_pred.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $MULTIMODAL_MODEL_DIR/scripts/2.1.predict_multimodal.py \
#     --elements $GM12878_ELEMENTS \
#     --genome $GENOME \
#     --signal-plus-bw $EP300_PLUS \
#     --signal-minus-bw $EP300_MINUS \
#     --accessibility-bw $ATAC_BW \
#     --model-dir $MULTIMODAL_MODEL_DIR/models/atac \
#     --fold-json $PROJECT_DIR/reference/hg38_five_folds.json \
#     --peaks $GM12878_PEAKS \
#     --output-dir $MULTIMODAL_PRED_DIR \
#     --batch-size 512 --device cpu"


# ============================================================
# STAGE 4: EVALUATE TRANSFERABILITY
# ============================================================
# Compare predicted vs. observed GM12878 p300 log counts (Pearson r)
# Use same script as K562: scripts/2.3.compute_prediction_performance.py

## [4.1] GM12878 BPNet (in-cell-type ceiling) — SUBMITTED job 28079524
# log: $LOG_DIR/eval_gm12878_bpnet.28079524.txt
# sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
#   --job-name=gm_eval_gm12878 \
#   --output=$LOG_DIR/eval_gm12878_bpnet.%j.txt \
#   --error=$LOG_DIR/eval_gm12878_bpnet.%j.txt \
#   --wrap="source ~/.bashrc && conda activate tfmodisco && \
#     python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
#       --mean-pred-dir $GM12878_V1_PRED_DIR/mean \
#       --peaks $GM12878_PEAKS \
#       --overlap-col EP300_peak_overlap \
#       --h5-name ENCSR000DZG_split000_predictions.h5 \
#       --output-dir $GM12878_V1_PRED_DIR"

## [4.2] K562 v1 BPNet (sequence-only cross-cell-type) — SUBMITTED job 28079525
# log: $LOG_DIR/eval_k562_v1.28079525.txt
# sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
#   --job-name=gm_eval_k562v1 \
#   --output=$LOG_DIR/eval_k562_v1.%j.txt \
#   --error=$LOG_DIR/eval_k562_v1.%j.txt \
#   --wrap="source ~/.bashrc && conda activate tfmodisco && \
#     python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
#       --mean-pred-dir $K562_V1_PRED_DIR/mean \
#       --peaks $GM12878_PEAKS \
#       --overlap-col EP300_peak_overlap \
#       --h5-name ENCSR000EGE_split000_predictions.h5 \
#       --output-dir $K562_V1_PRED_DIR"

## [4.3] K562 multimodal ATAC BPNet (multimodal cross-cell-type) — SUBMITTED job 28079529
# log: $LOG_DIR/eval_multimodal.28079529.txt
# sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
#   --job-name=gm_eval_mm \
#   --output=$LOG_DIR/eval_multimodal.%j.txt \
#   --error=$LOG_DIR/eval_multimodal.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $MULTIMODAL_MODEL_DIR/scripts/2.2.plot_prediction_accuracy.py \
#     --predictions-dir $MULTIMODAL_PRED_DIR \
#     --peaks $GM12878_PEAKS \
#     --max-counts 10"


# ============================================================
# NOTE: eval jobs 28079524/525/529 (2026-06-06) all failed.
# Fixed and resubmitted 2026-06-07:
#   [4.1] job 28116257 — gm_eval_gm12878 — COMPLETED
#   [4.2] job 28116346 — gm_eval_k562v1  — COMPLETED
#   [4.3] job 28116354 — gm_eval_mm      — COMPLETED
# Fixes: pixi run -e ism instead of conda; --cv-pred-dir made optional
# in 2.3.compute_prediction_performance.py; seaborn added to ism pixi env;
# multimodal eval corrected to use 2.2.plot_prediction_accuracy.py without --peaks.


# ============================================================
# STAGE 5: SHAP ON GM12878 EP300 BPNET (peaks only, counts head)
# ============================================================
# Compute per-fold SHAP scores on GM12878 p300 peaks (21,068 regions),
# then merge per-chromosome, then average across folds.
# Uses: conda bpnet_37 + cuda/11.1.1 cudnn/8.1.1.33 (via 3.1.submit_mean_shap_one_fold.sh)
# Signal BigWigs: ENCODE processed (ENCFF960OFK_plus.bw, ENCFF941MGK_minus.bw) — match training
# Bias: K562 control BigWigs as proxy (GM12878 controls not available)

SHAP_DIR=$THIS_DIR/shap_peaks
CONFIG=$THIS_DIR/config/input_data_gm12878_shap.json
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

## [5.1] Submit SHAP per-fold — COMPLETE (all 5 folds, 24/24 chroms each) 2026-06-08
# NOTE: per-chromosome array approach (jobs 28131714–28131722) was cancelled — overkill for 21k peaks.
# Replaced with one SLURM job per fold looping all chromosomes internally (~2h/fold vs 25 array tasks).
# Script: $THIS_DIR/scripts/3.1.submit_gm12878_shap.sh
# chrY has 0 GM12878 p300 peaks — DONE.txt created manually; chrY removed from submit script.

# First attempt: job array 28132251 (2026-06-07, 4h limit) — timed out at 6–19/24 chroms/fold
# Second attempt: job array 28324270 (2026-06-08, 8h limit) — all folds complete
# sbatch --array=0-4 $THIS_DIR/scripts/3.1.submit_gm12878_shap.sh

## [5.2+5.3] Merge per-chromosome h5s and compute mean across folds — SUBMITTED job 28394331 2026-06-08
# Script: $THIS_DIR/scripts/5.2.submit_merge_mean_shap.sh
# Outputs: $SHAP_DIR/fold{0-4}/shap_counts_merged.h5, $SHAP_DIR/all_folds/counts_mean_shap_scores.h5
# sbatch $THIS_DIR/scripts/5.2.submit_merge_mean_shap.sh


# ============================================================
# STAGE 6: TFMODISCO ON GM12878 MEAN SHAP
# ============================================================
# Parameters mirror K562 shap_peaks MoDISCo run.

MODISCO_DIR=$THIS_DIR/modisco
MOTIF_REF=$PROJECT_DIR/reference/MotifCompendium-Database-Human.meme.txt
SHAP_MEAN_H5=$SHAP_DIR/all_folds/counts_mean_shap_scores.h5
MODISCO_WIDTH=400

## [6.1] Run TF-MoDISco — submit after step 5.3 completes
# sbatch $SCRIPTS_DIR/4.1.submit_counts_modisco.sh \
#   $MODISCO_DIR/max_seqlets_250k_30_10_0 \
#   $SHAP_MEAN_H5 \
#   250000 \
#   30 \
#   10 \
#   0 \
#   $MOTIF_REF \
#   $MODISCO_WIDTH


# ============================================================
# STAGE 7: PLOT TRANSFERABILITY RESULTS
# ============================================================

## [7.1] Fig 2d + S2a/b — transferability bar chart and scatter plots — job 28131732 COMPLETED (8m27s)
# Outputs: $THIS_DIR/figures/transferability_{bar,scatter_all,scatter_peaks}.pdf
# Note: initial outputs went to $PROJECT_DIR/figures/ (wrong); fixed --output-dir default in script.
# sbatch --partition=owners,engreitz,normal --time=0:30:00 --mem=32G \
#   --job-name=plot_transfer \
#   --output=$LOG_DIR/plot_transferability.%j.txt \
#   --error=$LOG_DIR/plot_transferability.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
#     $SCRIPTS_DIR/plot_transferability.py"


# ============================================================
# STAGE 8: GM12878 MULTIMODAL BPNET (in-cell-type ceiling)
# ============================================================
# Train same architecture as K562 ATAC multimodal, but using GM12878 EP300 signal + GM12878 ATAC.
# Purpose: in-cell-type multimodal ceiling to compare against cross-cell-type K562 transfer.
# Config: $THIS_DIR/config/input_data_gm12878_multimodal.json
# Script: $THIS_DIR/scripts/1.1.submit_training_gm12878_multimodal.sh
# Output: $THIS_DIR/GM12878_multimodal_BPNet/models/atac/fold{0-4}/

GM12878_MM_DIR=$THIS_DIR/GM12878_multimodal_BPNet
GM12878_MM_CONFIG=$THIS_DIR/config/input_data_gm12878_multimodal.json

## [8.1] Submit GM12878 multimodal training — SUBMITTED 2026-06-08 as job array 28359131
# Folds 0/1/2 running, 3/4 pending; 24h time limit
# Logs: $THIS_DIR/log/gm12878_mm_f{N}.28359131.txt
# sbatch --array=0-4 $THIS_DIR/scripts/1.1.submit_training_gm12878_multimodal.sh

## [8.2] Predict on 154k GM12878 elements (after training completes)
# Use same predict_multimodal.py as K562 but with GM12878 signal/ATAC BigWigs
# sbatch --partition=owners,engreitz,normal --time=8:00:00 --mem=64G \
#   --job-name=gm12878_mm_pred \
#   --output=$LOG_DIR/gm12878_mm_pred.%j.txt \
#   --error=$LOG_DIR/gm12878_mm_pred.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/2.1.predict_multimodal.py \
#     --elements $GM12878_ELEMENTS \
#     --genome $GENOME \
#     --signal-plus-bw $EP300_PLUS \
#     --signal-minus-bw $EP300_MINUS \
#     --accessibility-bw $ATAC_BW \
#     --model-dir $GM12878_MM_DIR/models/atac \
#     --fold-json $PROJECT_DIR/reference/hg38_five_folds.json \
#     --peaks $GM12878_PEAKS \
#     --output-dir $THIS_DIR/predictions/gm12878_multimodal_atac \
#     --batch-size 512 --device cpu"

## [8.2] GM12878 multimodal predictions — COMPLETE job 28387044 2026-06-08
# Output: $THIS_DIR/predictions/gm12878_multimodal_atac/{cv,mean}_predictions.tsv.gz
# sbatch $THIS_DIR/scripts/2.1.submit_predict_gm12878_multimodal.sh

## [8.3] Evaluate GM12878 multimodal predictions — SUBMITTED job 28438419 2026-06-08
# Note: predict script outputs TSVs (not h5); use --from-tsv + --mean-output-dir
# Output: $THIS_DIR/predictions/gm12878_multimodal_atac/prediction_accuracy.tsv
# sbatch --partition=owners,engreitz,normal --time=0:30:00 --mem=32G \
#   --job-name=gm12878_mm_eval \
#   --output=$THIS_DIR/log/gm12878_mm_eval.%j.txt \
#   --error=$THIS_DIR/log/gm12878_mm_eval.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
#     $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
#     --mean-output-dir $THIS_DIR/predictions/gm12878_multimodal_atac \
#     --overlap-col EP300_peak_overlap \
#     --from-tsv"

# ============================================================
# Stage 9: ATAC-only model (GM12878 in-cell-type + K562 transfer)
# ============================================================

## [9.1] GM12878 ATAC-only training — COMPLETE 2026-06-08 (job array 28443287); all 5 folds done in 2-9 min
# Script: scripts/1.2.submit_training_gm12878_atac_only.sh
# Logs: log/train_gm12878_atac_only.28443287_{0-4}.txt
# Output: GM12878_ATAC_only_BPNet/models/atac_only/fold{0-4}/multimodal_bpnet.torch
# sbatch --array=0-4 $THIS_DIR/scripts/1.2.submit_training_gm12878_atac_only.sh

## [9.2] K562 ATAC-only predictions on GM12878 elements (cross-cell-type transfer) — SUBMITTED job 28443298 2026-06-08
# Uses K562-trained ATAC-only model; predicts on 154k GM12878 elements
# Output: predictions/k562_atac_only/
# sbatch --partition=owners,engreitz,normal --time=8:00:00 --mem=64G \
#   --job-name=k562_atac_gm12878 \
#   --output=$THIS_DIR/log/predict_k562_atac_only_on_gm12878.%j.txt \
#   --error=$THIS_DIR/log/predict_k562_atac_only_on_gm12878.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/2.1.predict_multimodal.py \
#     --mode atac \
#     --elements $THIS_DIR/reference/GM12878_candidate_elements.narrowPeak \
#     --signal-plus-bw $THIS_DIR/data/EP300_plus.bw \
#     --signal-minus-bw $THIS_DIR/data/EP300_minus.bw \
#     --accessibility-bw $THIS_DIR/data/atac.bw \
#     --model-dir $PROJECT_DIR/2026_0529_multimodal_p300_model/models/atac_only \
#     --fold-json $PROJECT_DIR/reference/hg38_five_folds.json \
#     --peaks /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz \
#     --output-dir $THIS_DIR/predictions/k562_atac_only \
#     --batch-size 512 --device cpu"

## [9.3] Evaluate K562 ATAC-only transfer predictions — SUBMITTED job 28489944 2026-06-08
# sbatch --partition=owners,engreitz,normal --time=0:30:00 --mem=32G \
#   --job-name=k562_atac_gm12878_eval \
#   --output=$THIS_DIR/log/eval_k562_atac_only.%j.txt \
#   --error=$THIS_DIR/log/eval_k562_atac_only.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
#     $PROJECT_DIR/scripts/2.3.compute_prediction_performance.py \
#     --mean-pred-dir $THIS_DIR/predictions/k562_atac_only \
#     --mean-output-dir $THIS_DIR/predictions/k562_atac_only \
#     --overlap-col EP300_peak_overlap \
#     --from-tsv"

## [9.4] GM12878 ATAC-only predictions on GM12878 elements (in-cell-type) — COMPLETE job 28454690 2026-06-08
## [9.5] Evaluate GM12878 ATAC-only in-cell-type predictions — COMPLETE job 28499108 2026-06-08
# K562 ATAC-only → GM12878 eval results (job 28489944): Pearson=0.717 all, 0.467 p300+
# GM12878 ATAC-only in-cell-type eval results (job 28499108): Pearson=0.683 all, 0.579 p300+
# sbatch --partition=owners,engreitz,normal --time=8:00:00 --mem=64G \
#   --job-name=gm12878_atac_pred \
#   --output=$THIS_DIR/log/predict_gm12878_atac_only.%j.txt \
#   --error=$THIS_DIR/log/predict_gm12878_atac_only.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/2.1.predict_multimodal.py \
#     --mode atac \
#     --elements $THIS_DIR/reference/GM12878_candidate_elements.narrowPeak \
#     --signal-plus-bw $THIS_DIR/data/ENCFF960OFK_plus.bw \
#     --signal-minus-bw $THIS_DIR/data/ENCFF941MGK_minus.bw \
#     --accessibility-bw $THIS_DIR/data/atac.bw \
#     --model-dir $THIS_DIR/GM12878_ATAC_only_BPNet/models/atac_only \
#     --fold-json $PROJECT_DIR/reference/hg38_five_folds.json \
#     --peaks /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz \
#     --output-dir $THIS_DIR/predictions/gm12878_atac_only \
#     --batch-size 512 --device cpu"

# ============================================================
# [10] Reverse direction: evaluate GM12878-trained models on K562 elements
# ============================================================

## [10.1] GM12878 ATAC-only BPNet predictions on K562 elements — COMPLETE 2026-07-09 (job 33320529, 2:58:20)
# Script: scripts/6.1.submit_predict_gm_atac_only_on_k562.sh
# Output: predictions/gm12878_atac_only_on_k562/{cv,mean}_predictions.tsv.gz
# sbatch scripts/6.1.submit_predict_gm_atac_only_on_k562.sh

## [10.2] GM12878 multimodal BPNet predictions on K562 elements — job 33320530 FAILED (8:41, exit 1); resubmitted as 33328385, COMPLETE 2026-07-09 (2:39:10)
# Script: scripts/6.2.submit_predict_gm_multimodal_on_k562.sh
# Output: predictions/gm12878_multimodal_on_k562/{cv,mean}_predictions.tsv.gz
# sbatch scripts/6.2.submit_predict_gm_multimodal_on_k562.sh

## [10.3] GM12878 sequence-only BPNet predictions on K562 elements — COMPLETE 2026-07-09 (job array 33320531_{0-4}, ~1h/fold)
# Script: scripts/6.3.submit_predict_gm_bpnet_on_k562.sh
# Output: predictions/gm12878_bpnet_on_k562/mean/fold{0-4}/ENCSR000DZG_split000_predictions.h5
# bash scripts/6.3.submit_predict_gm_bpnet_on_k562.sh

## [10.4] Evaluate all 3 GM->K562 predictions (aggregate CV Pearson/Spearman across folds) — SUBMITTED 2026-07-09 (jobs 33357106, 33357107, 33357108)
# Script: scripts/6.4.submit_eval_gm_on_k562.sh
# [1] 33357106: gm_bpnet_on_k562 eval -> FAILED (unrecognized arg --output-dir; 2.3.compute_prediction_performance.py
#     only accepts --cv-output-dir/--mean-output-dir). Fixed script (dropped --output-dir) and resubmitted as 33357837.
# [2] 33357107: gm_atac_only_on_k562 eval -> predictions/gm12878_atac_only_on_k562/ (RUNNING)
# [3] 33357108: gm_multimodal_on_k562 eval -> predictions/gm12878_multimodal_on_k562/ (RUNNING)
# bash scripts/6.4.submit_eval_gm_on_k562.sh

## [10.5] Plot reverse-direction transferability bar chart + scatter panels — PENDING (waiting on 10.4 to finish)
# Script: scripts/plot_transferability_on_k562.py (already written, mirrors scripts/plot_transferability.py)
# 7 models: K562 seq-only (in-cell ceiling), GM12878 seq-only/ATAC-only/multimodal (cross-cell-type),
#           K562 ATAC-only/multimodal (in-cell-type), K562 inter-replicate ceiling
# Outputs: figures/transferability_bar_on_k562.pdf (+ _all_elements/_p300plus split),
#          figures/transferability_scatter_on_k562_{all,peaks}.pdf
# pixi run -e ism python scripts/plot_transferability_on_k562.py --output-dir 2026_0606_GM12878_transferability/figures/

## [10.4-result] Eval results (2026-07-09):
#   GM12878 seq-only BPNet -> K562:   Pearson = 0.5352 all, 0.3374 p300+ (job 33357837)
#   GM12878 ATAC-only BPNet -> K562:  Pearson = 0.5969 all, 0.3783 p300+ (job 33357107)
#   GM12878 multimodal BPNet -> K562: Pearson = 0.6838 all, 0.4507 p300+ (job 33357108)

## [10.5] Plot reverse-direction transferability bar chart + scatter panels — SUBMITTED job 33358993 2026-07-09
# Script: scripts/plot_transferability_on_k562.py (already existed, mirrors scripts/plot_transferability.py)
# Submitted as sbatch job (not run on login node, per Sherlock policy) — COMPLETE 2026-07-09 (16:59, 98.5% CPU eff);
# sbatch --partition=owners,engreitz,normal --time=0:30:00 --mem=32G \
#   --job-name=plot_transfer_k562 \
#   --output=$LOG_DIR/plot_transferability_on_k562.%j.txt \
#   --error=$LOG_DIR/plot_transferability_on_k562.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
#     $SCRIPTS_DIR/plot_transferability_on_k562.py --output-dir $THIS_DIR/figures/"
