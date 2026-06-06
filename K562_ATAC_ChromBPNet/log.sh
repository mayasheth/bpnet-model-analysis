
### INPUT FILES & PARAMS ###
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/K562_ATAC_ChromBPNet
DATA_DIR=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/ATAC/ENCSR467RSV

## reference files
ABC_peaks=$OAK/Users/sheth/ENCODE_rE2G/results/2025_0226_ATAC_powerlaw_models/ATAC_H3K27ac_powerlaw/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
chr_list=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.txt
genome=$OAK/Users/sheth/hg38_resources/hg38.fa

## derived paths
pred_dir=$THIS_DIR/predictions
peaks_use=$THIS_DIR/data/K562_candidate_elements.narrowPeak
shap_dir=$THIS_DIR/shap
model_path_pattern="$DATA_DIR/models/fold_{fold}/model.chrombpnet_nobias.fold_{fold}.ENCSR868FGK.h5"


### STEP 1: DOWNLOAD MODELS [DONE] ###
# Models downloaded from ENCODE (ENCSR868FGK), 5 folds
# Output: $DATA_DIR/models/fold_{0-4}/model.chrombpnet_nobias.fold_{fold}.ENCSR868FGK.h5

cd $DATA_DIR
wget https://www.encodeproject.org/files/ENCFF231MHU/@@download/ENCFF231MHU.tar.gz
mkdir models & cd models
tar -xvf ../ENCFF574YLK.tar.gz


### STEP 2: PREPARE REGIONS [DONE] ###
# Reformat ABC candidate elements to narrowPeak format
# Output: $peaks_use (~154k regions)

mkdir -p $THIS_DIR/data
sh $PROJECT_DIR/scripts/0.2.bed_to_narrowPeak.sh $ABC_peaks $peaks_use 1057


### STEP 3: PREDICT ON CANDIDATE ELEMENTS [DONE - needs rerun, see bottom] ###
# Per-fold predictions of log counts AND per-base profile for all candidate elements.
# Uses pixi ism env (TF-GPU 2.17) for GPU inference.
# Output: $pred_dir/fold{0-4}_predictions.{tsv,h5}
#   h5 keys: coords/{chrom,center,start,end}, predictions/{logcounts,profile (N,1000)}
#
# Note: original runs (2026-04-23) did not include profile output. Re-run required.

sbatch $PROJECT_DIR/scripts/2.4.submit_mean_predict_chrombpnet.sh \
    "$model_path_pattern" \
    $chr_sizes \
    $genome \
    $pred_dir \
    $peaks_use


### STEP 4: MEAN PREDICTIONS ACROSS FOLDS [DONE - needs rerun, see bottom] ###
# Average log-counts and profile across the 5 folds.
# Reads per-fold h5 files (not TSV).
# Output: $pred_dir/mean_predictions.tsv  (153,545 regions, log-counts)
#         $pred_dir/mean_predictions.h5   (pred_logcounts, pred_prof (N,1000))
#
# Note: original run (2026-04-23) used TSV and produced no h5/profile. Re-run required.

pixi run python $PROJECT_DIR/scripts/2.5.mean_predictions_chrombpnet.py \
    --pred-dir $pred_dir \
    --output $pred_dir/mean_predictions.tsv \
    --output-h5 $pred_dir/mean_predictions.h5


### STEP 5: COMPUTE SHAP SCORES (counts head) [DONE] ###
# DeepSHAP contribution scores per fold, split by chromosome for parallelism.
# Each fold submits 24 chromosome jobs (GPU), then a merge job runs on dependency.
# Output: $shap_dir/fold{0-4}/shap_counts_merged.h5
#
# Debugging history:
#   Initial runs (2026-04-10/11) failed because the pixi ism environment has
#   standalone Keras 3 (3.10.0) + TF 2.17, which broke shap.TFDeepExplainer:
#     1. tf.reduce_sum() cannot accept KerasTensor (fixed with Lambda layer)
#     2. tf.compat.v1.keras.backend.get_session() removed in Keras 3
#     3. KerasTensor has no .op attribute (used internally by TFDeepExplainer)
#   Fix: switched 3.2.submit_shap_chrombpnet.sh to use conda bpnet_37
#   (TF 2.4.1 + Keras 2), which is compatible with shap.TFDeepExplainer.
#   The merge job was also switched from pixi to bpnet_37.
#
# Status (2026-04-13):
#   Folds 0-2: running with bpnet_37 (jobs 21440181/21440204/21440214)
#              merge jobs 21440184/21440204/21440214 pending dependency
#              merge jobs for folds 1+2 hit QOS limit; resubmit after chr jobs finish
#   Folds 3-4: not yet submitted; resubmit with updated script when queue opens

for fold in 0 1 2; do
    bash $PROJECT_DIR/scripts/3.2.submit_shap_chrombpnet.sh \
        $fold \
        "$model_path_pattern" \
        $peaks_use \
        $shap_dir \
        $genome
done

# Folds 3 and 4 submitted 2026-04-13 (jobs 21452258/21452261):
for fold in 3 4; do
    bash $PROJECT_DIR/scripts/3.2.submit_shap_chrombpnet.sh \
        $fold \
        "$model_path_pattern" \
        $peaks_use \
        $shap_dir \
        $genome
done

# All 5 folds resubmitted 2026-04-15 after fixing GPU issue:
#   Root cause: cuda/11.1.1 lacks libcusolver.so.10 (renamed to .so.11 in CUDA 11.1+),
#   so TF 2.4 silently fell back to CPU (~100 regions/hr vs. thousands/hr on GPU).
#   Fix: switched 3.2.submit_shap_chrombpnet.sh to module load cuda/11.2.0 (has libcusolver.so.10).
#   GPU confirmed working on A100 (~400 regions/min).
#   Prior runs only completed chrY for folds 0, 3, 4; nothing for folds 1, 2.
#
# Fold 0: resubmitted 2026-04-15, 23 chr jobs (chrY already done), merge job 21608037
# Fold 1: resubmitted 2026-04-15, 24 chr jobs, merge job 21609863
# Fold 2: resubmitted 2026-04-15, 24 chr jobs, merge job 21609870
# Folds 3-4: resubmit after folds 0-2 clear the queue

for fold in 0 1 2; do
    bash $PROJECT_DIR/scripts/3.2.submit_shap_chrombpnet.sh \
        $fold \
        "$model_path_pattern" \
        $peaks_use \
        $shap_dir \
        $genome
done

# Merge jobs for folds 0-2 failed initially because 'conda activate' does not work in
# SLURM --wrap context. Fix: switched to 'conda run -n bpnet_37' in 3.2.submit_shap_chrombpnet.sh.
# Merge jobs resubmitted manually 2026-04-15 (jobs 21646699/21646700/21646701).
# All three completed successfully: 153,545 regions, shape (153545, 2114, 4).

# Folds 3-4 submitted 2026-04-15 (23 chr jobs each, chrY already done):
#   Fold 3 chr jobs: 21648240, merge: 21648241
#   Fold 4 chr jobs: 21648243, merge: 21648244
# All 5 folds completed 2026-04-15: 153,545 regions, shape (153545, 2114, 4) per fold.

for fold in 3 4; do
    bash $PROJECT_DIR/scripts/3.2.submit_shap_chrombpnet.sh \
        $fold \
        "$model_path_pattern" \
        $peaks_use \
        $shap_dir \
        $genome
done


### STEP 6: MEAN SHAP ACROSS FOLDS [DONE] ###
# Average hyp_scores across all 5 per-fold merged h5 files.
# Output: $shap_dir/all_folds/counts_mean_shap_scores.h5  (hyp_scores (153545, 2114, 4))
#         $shap_dir/all_folds/counts_peaks_valid_scores.bed
# Completed 2026-04-23.

bash $PROJECT_DIR/scripts/3.5.submit_mean_shap_chrombpnet.sh $shap_dir


### STEP 7: BIGWIGS [TODO - see bottom] ###
# Generate BigWig files from mean SHAP scores and mean profile predictions.
# Output: $shap_dir/all_folds/counts_scores.bw
#         $shap_dir/all_folds/counts_scores.stats.txt
#         $pred_dir/mean_predictions_profile.bw
#
# Requires mean_predictions.h5 from step 4 (re-run) before profile BigWig can be made.

bash $PROJECT_DIR/scripts/4.0.submit_make_bigwigs_chrombpnet.sh \
    $shap_dir \
    $pred_dir \
    $chr_sizes