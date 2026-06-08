# ---  BPNET ENV --- # 
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

# --- FILE PATHS --- #
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
RESULTS_DIR=$PROJECT_DIR/K562_GATA1_BPNet
SCRIPTS_DIR=$PROJECT_DIR/scripts

DATA_CONFIG=$RESULTS_DIR/config/input_data_predict.json
CANDIDATE_ELEMENTS=$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak
MODEL_PATH_PATTERN="$RESULTS_DIR/models/fold{fold}/model_split000"
MEAN_PRED_DIR=$RESULTS_DIR/predictions_mean; #mkdir -p $MEAN_PRED_DIR
CV_PRED_DIR=$RESULTS_DIR/predictions_cv; #mkdir -p $CV_PRED_DIR

# --- REFERENCE FILES --- #
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
CHR_LIST=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.txt
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
BLACKLIST=$OAK/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_unified_blacklist.bed

### --- [1] DOWNLOAD --- ### 

# input bws and peaks, trained models
sh $SCRIPTS_DIR/0.1.download_tfatlas_bpnet.sh ENCSR000EWM $RESULTS_DIR

### --- [2] PREDICT --- ###

## --- [2A] CV PREDICTIONS --- ##

ARRAY="0-4"
sbatch \
  --job-name=cvp_f%a \
  --output=log/cvp_f%a.%A.txt \
  --error=log/cvp_f%a.%A.txt \
  --array=$ARRAY \
  $SCRIPTS_DIR/2.1.submit_cv_predict.sh \
    $MODEL_PATH_PATTERN \
    $CHR_SIZES \
    $GENOME \
    $CV_PRED_DIR \
    $DATA_CONFIG

## --- [2B] MEAN PREDICTIONS --- ##

# --- [2B.i] GENERATE PREDICTIONS FROM EACH FOLD MODEL --- #

ARRAY="0-4"
sbatch \
  --job-name=mp_f%a \
  --output=log/mp_f%a.%A.txt \
  --error=log/mp_f%a.%A.txt \
  --array=$ARRAY \
  $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
    $MODEL_PATH_PATTERN \
    $CHR_SIZES \
    $GENOME \
    $MEAN_PRED_DIR \
    $DATA_CONFIG

# --- [2B.ii] CALCULATE MEAN PREDICTIONS --- #

conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

CHROM_ANNOT_FILE=/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv
COLS_KEEP="chrom,start,end,EP300.RPM,EP300_peak_overlap"

mkdir -p $MEAN_PRED_DIR/all_folds
PRED_H5=model_split000_predictions.h5
PRED_H5_LIST="$MEAN_PRED_DIR/fold0/$PRED_H5,$MEAN_PRED_DIR/fold1/$PRED_H5,$MEAN_PRED_DIR/fold2/$PRED_H5,$MEAN_PRED_DIR/fold3/$PRED_H5,$MEAN_PRED_DIR/fold4/$PRED_H5"

python $BPNET_DIR/utils/mean_predictions.py \
  --prediction_h5s $PRED_H5_LIST \
  --chrom_sizes $CHR_SIZES \
  --output_dir $MEAN_PRED_FILE/all_folds \
  --chrom_annot $CHROM_ANNOT_FILE \
  --chrom_cols $COLS_KEEP

## --- [2C] PREDICTION METRICS --- ##

conda activate tfmodisco

python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
  --mean_pred_dir $MEAN_PRED_DIR \
  --cv_pred_dir $CV_PRED_DIR \
  --overlap_col "EP300_peak_overlap"

### --- [2D] PREDICTION RESUBMISSION (2026-06-07) --- ###
#
# Root cause: jobs were hanging, not slow. --threads 2 spawns 2 parallel generator
# workers. When one worker crashes (BigWig read error), the stealer thread blocks
# forever waiting for batches that never arrive.
#
# Mean predictions: chrM:15952-16568 is the last candidate element (row 150527/150528).
# Its 2114bp window (end=17317) exceeds the GATA1 BigWig chrM size (16569 bp),
# causing pyBigWig RuntimeError and deadlocking the generator.
#
# CV fold 2: cause unknown (chr4/11/12/15/Y have no out-of-bounds elements);
# running --threads 1 to surface exact element with a clean error.
#
# Fix: GPU constraint, chrM excluded from mean predictions, --threads 1 for CV fold 2.
# Scripts: resubmit_cv_fold2.sh, resubmit_mean_predict.sh

cd $RESULTS_DIR
sbatch resubmit_cv_fold2.sh           # job 28137653 — CV fold 2
sbatch --array=0-4 resubmit_mean_predict.sh  # job 28137654 — mean predictions folds 0-4

# After mean h5s complete, calculate mean across folds then run performance metrics:
# (same as [2B.ii] / [2C] above, using pixi run -e ism instead of conda)
module load devel pixi/0.53.0

PRED_H5=model_split000_predictions.h5
PRED_H5_LIST="$MEAN_PRED_DIR/fold0/$PRED_H5,$MEAN_PRED_DIR/fold1/$PRED_H5,$MEAN_PRED_DIR/fold2/$PRED_H5,$MEAN_PRED_DIR/fold3/$PRED_H5,$MEAN_PRED_DIR/fold4/$PRED_H5"

conda activate bpnet_37
python $BPNET_DIR/utils/mean_predictions.py \
  --prediction_h5s $PRED_H5_LIST \
  --chrom_sizes $CHR_SIZES \
  --output_dir $MEAN_PRED_DIR/all_folds

pixi run -e ism python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
  --mean-pred-dir $MEAN_PRED_DIR \
  --cv-pred-dir $CV_PRED_DIR \
  --peaks $PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak \
  --overlap-col "GATA1_peak_overlap" \
  --h5-name model_split000_predictions.h5

### --- [3] SHAP (ALL ELEMENTS) --- ###

## --- [3A] COMPUTE SHAP PER FOLD --- ##
# folds complete: 0, 1, 2
# folds submitted: 3, 4

FOLD=4
sh $SCRIPTS_DIR/3.1.submit_mean_shap_one_fold.sh \
   $FOLD \
   $MODEL_PATH_PATTERN \
   $DATA_CONFIG \
   $CANDIDATE_ELEMENTS \
   $RESULTS_DIR/shap \
   $GENOME \
   $CHR_SIZES \
   $BPNET_DIR


### --- [6] MOTIF SPACING --- ###

## --- [6A] 2,3,4 MOTIFS --- ## 
# took just over 9h
qsjob -j sp_g234 --cpu -m 160G -t 18:00:00 --env bpnet_37 "ml cuda/11.1.1 cudnn/8.1.1.33; python $OAK/Users/sheth/EP300_BPNet/scripts/6.0.motif_spacing.one_motif.py --analysis-id GATA_50bp_n234 --motif-name GATA --motif-seq GATAA --motif-counts 2 3 4 --narrow-peak-type all --model-type GATA1 --out-dir $OAK/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/motif_spacing"

python $OAK/Users/sheth/EP300_BPNet/scripts/plot_motif_spacing.py $OAK/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv --out-dir $OAK/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234 --motif-name GATA

############################################################
#### SHAP (COUNTS ONLY) ####

# 1) run each fold on all chr, parallelize jobs across chr (chr1 requires 12h, the rest only need 5h)(ARRAY="0-23")
# folds done: 0, 1, 2, 3, 4
# run on just p300 peaks: 0, 1, 2, 3, 4

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
project_dir=$OAK/Users/sheth/EP300_BPNet
pred_data_config=$project_dir/2025_0703_retrain_p300_model/config/input_data_predict.json
peaks_use=$project_dir/2025_0325_K562_BPNet/data/input_peaks.narrowPeak
shap_dir=$project_dir/2025_0703_retrain_p300_model/shap
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
model_path_pattern="$project_dir/2025_0703_retrain_p300_model/models/fold{fold}/model_split000"

peak_shap_dir=$project_dir/2025_0703_retrain_p300_model/shap_peaks
train_data_config=$project_dir/2025_0703_retrain_p300_model/config/input_data.json
p300_peaks_use=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0325_K562_BPNet/data/peak_inliers.bed

FOLD=4
sh $SCRIPTS_DIR/4.1.submit_mean_shap_one_fold.sh \
   $FOLD \
   $model_path_pattern \
   $train_data_config \
   $p300_peaks_use \
   $peak_shap_dir \
   $genome \
   $chr_sizes \
   $BPNET_DIR

## 2) combine h5 across chromosomes across folds
# done: 0, ..., 4
# peaks only: 0, ..., 4

conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

SHAP_ID=shap_peaks
FOLD=0
python $BPNET_DIR/utils/merge_shap_across_chrom.py \
  --input-dir $OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/${SHAP_ID}/fold${FOLD} \
  --h5-filename counts_scores.h5 \
  --output-file $OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/${SHAP_ID}/fold${FOLD}/shap_counts_merged.h5

# 3) compute mean shap across folds 
# outputs: $out_dir/counts_mean_shap_scores.h5 and $out_dir/counts_peaks_valid_scores.bed
project_dir=$OAK/Users/sheth/EP300_BPNet
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
SHAP_ID=shap_peaks

shap_dir=$project_dir/2025_0703_retrain_p300_model/${SHAP_ID}
out_dir=$shap_dir/all_folds; mkdir -p $out_dir
shap_h5=shap_counts_merged.h5
shap_h5_list="$shap_dir/fold0/$shap_h5,$shap_dir/fold1/$shap_h5,$shap_dir/fold2/$shap_h5,$shap_dir/fold3/$shap_h5,$shap_dir/fold4/$shap_h5"

python $BPNET_DIR/utils/mean_shap_plus_peaks.py \
  --counts_shaps $shap_h5_list \
  --output_dir $out_dir 


# 4) make merged shap bw
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
project_dir=$OAK/Users/sheth/EP300_BPNet
out_dir=$project_dir/2025_0703_retrain_p300_model/shap/all_folds

merged_shap_h5=$out_dir/counts_mean_shap_scores.h5
merged_peaks=$out_dir/counts_peaks_valid_scores.bed

python $BPNET_DIR/utils/run_importance_hdf5_to_bigwig.py \
  --hdf5_path $merged_shap_h5 \
  --regions_path $merged_peaks \
  --chrom_sizes $chr_sizes \
  --outfile $out_dir/counts_scores.bw \
  --outstats $out_dir/counts_scores.stats.txt

#### MODISCO (COUNTS ONLY) ####
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0703_retrain_p300_model
motif_reference=$base_dir/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
max_seqlets=250000
modisco_width=400

SHAP_ID=shap_peaks
MODISCO_ID=modisco_peaks
shap_counts_h5=$this_dir/${SHAP_ID}/all_folds/counts_mean_shap_scores.h5

## longer motifs
modisco_counts_dir=$this_dir/${MODISCO_ID}/max_seqlets_250k_20_5_10
trim_size=20
initial_flank_to_add=5
final_flank_to_add=10

## shorter motifs
modisco_counts_dir=$this_dir/${MODISCO_ID}/max_seqlets_250k_30_10_0
trim_size=30
initial_flank_to_add=10
final_flank_to_add=0

## submit
sbatch $SCRIPTS_DIR/5.1.submit_counts_modisco.sh \
    $modisco_counts_dir \
    $shap_counts_h5 \
    $max_seqlets \
    $trim_size \
    $initial_flank_to_add \
    $final_flank_to_add \
    $motif_reference \
    $modisco_width


#### FINEMO ####

## 1) call hits
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model
modisco_width=400

SHAP_ID=shap_peaks
MODISCO_ID=modisco_peaks
FINEMO_ID=finemo_peaks
shap_counts_h5=$this_dir/${SHAP_ID}/all_folds/counts_mean_shap_scores.h5
peaks_use=$this_dir/${SHAP_ID}/all_folds/counts_peaks_valid_scores.bed

run_id=max_seqlets_250k_30_10_0
run_id=max_seqlets_250k_20_5_10

peak_width=500
out_dir=$this_dir/${FINEMO_ID}/pkw_${peak_width}_$run_id
modisco_counts_h5=$this_dir/${MODISCO_ID}/$run_id/counts_scores.h5

sbatch $SCRIPTS_DIR/6.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $peak_width \
  $modisco_counts_h5 \
  $modisco_width

## 2) annotate motifs to create config file...

## 3) call hits with curated motifs
conda activate analysis

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model
modisco_width=400
SHAP_ID=shap
FINEMO_ID=finemo

peaks_use=$this_dir/${SHAP_ID}/all_folds/counts_peaks_valid_scores.bed
peak_width=500
out_dir=$this_dir/${FINEMO_ID}/pkw_${peak_width}_curated_motifs
modisco_path=$OAK/Users/sheth/EP300_BPNet/reference/curated_motif_data_for_finemo

sbatch $SCRIPTS_DIR/6.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $peak_width \
  $modisco_path \
  $modisco_width



## 4) format finemo hits and plot
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model
pred_file=$this_dir/predictions_mean/all_folds/mean_predictions_counts.tsv.gz
motif_file=$OAK/Users/sheth/EP300_BPNet/reference/motif_annotations.tsv
target_name="p300 v2"

# 30 10 0
run_id="250k_30_10_0"
peak_width=500
finemo_dir=$this_dir/finemo/pkw_${peak_width}_max_seqlets_$run_id
out_dir=$finemo_dir/annotated_motifs

# 20 5 10
run_id="250k_20_5_10"
peak_width=500
finemo_dir=$this_dir/finemo/max_seqlets_$run_id
out_dir=$finemo_dir/annotated_motifs

python $SCRIPTS_DIR/6.2.format_finemo_hits.py \
  --finemo_dir $finemo_dir \
  --pred $pred_file \
  --motifs $motif_file \
  --target_name "$target_name" \
  --overlap_col "EP300_peak_overlap" \
  --modisco_id $run_id \
  --out_dir $out_dir

# curated motifs
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model
pred_file=$this_dir/predictions_mean/all_folds/mean_predictions_counts.tsv.gz

peak_width=400
finemo_dir=$this_dir/finemo/pkw_500_curated_motifs
out_dir=$finemo_dir/annotated_motifs

python $SCRIPTS_DIR/6.2.format_finemo_hits.py \
--finemo_dir $finemo_dir \
--pred $pred_file \
--target_name "$target_name" \
--overlap_col "EP300_peak_overlap" \
--out_dir $out_dir \
--stratify_hits

## 5) filter finemo hits?
# conda activate finemo_gpu

# all_hits=$finemo_dir/annotated_motifs/hits_renamed.bed.gz
# threshold=0.01
# $SCRIPTS_DIR/6.3.filter_finemo_hits.sh $all_hits $threshold

# output: $finemo_dir/annotated_motifs/hits_renamed.filtered_{threshold}.bed.gz