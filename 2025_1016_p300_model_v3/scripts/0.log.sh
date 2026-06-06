#### GOAL ####
# retrain p300 BPNet model with negatives = accessible but not p300-binding elements AND genomic negatives

#### BPNET ENV #### 
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

#### FILE PATHS ####
project_dir=$OAK/Users/sheth/EP300_BPNet
results_dir=$OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3

## REFERENCE FILES
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
chr_list=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.txt
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
blacklist=$OAK/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_unified_blacklist.bed

## INPUT DATA (peaks = IDR-thresholded from ENCODE)
# EP300: ENCSR000EGE, control: ENCSR000EHI
ep300_peaks=$OAK/Users/sheth/Data/ENCODE/K562/ENCFF702XPO.narrowPeak
ABC_peaks=$OAK/Users/sheth/ENCODE_rE2G_main/ENCODE_rE2G/results/2025_0313_ENCODE_K562_DNase/Stam_ENCFF860XAE_ENCFF205FNC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed

#### PREPROCESSING ####

## 1) MAKE STRANDED BW FILES FOR EXPERIMENT AND CONTROL ##

# Using bigWig files from Vivek, because the ones I made with this script miss entire chromosomes!?
# Downloaded to /data from here: https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/ENCSR000EGE/

sh $project_dir/scripts/make_training_bw.sh -i $ep300_bam1 $ep300_bam2 -c $control_bam -o $results_dir/data -g $chr_sizes

## 2) PREPARE INPUT POSITIVE PEAKS ## 

outliers_config=$OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/config/input_outliers.json
ep300_peaks_inliers=$OAK/Users/sheth/EP300_BPNet/reference/ENCSR000EGE_peaks_inliers.narrowPeak

# remove outliers: went from 28710 to 28516 peaks
bpnet-outliers \
        --input-data $outliers_config  \
        --quantile 0.99 \
        --quantile-value-scale-factor 1.2 \
        --task 0 \
        --chrom-sizes $chr_sizes \
        --chroms $(paste -s -d ' ' $chr_list) \
        --sequence-len 1057 \
        --blacklist $blacklist \
        --global-sample-weight 1.0 \
        --output-bed $ep300_peaks_inliers

## 3) PREPARE NEGATIVES FROM ACCESSIBLE ELEMENTS + GC-MATCHED BACKGROUND ##
# generate gc reference
bpnet-gc-reference \
        --ref_fasta $genome \
        --chrom_sizes $chr_sizes \
        --output_prefix $OAK/Users/sheth/EP300_BPNet/reference/genomewide_gc_stride_1000_flank_size_1057.gc \
        --inputlen 2114 \
        --stride 1000
    
# get gc-matched background regions 
bpnet-gc-background \
        --ref_fasta $genome \
        --peaks_bed $ep300_peaks_inliers \
        --out_dir $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data \
        --ref_gc_bed $OAK/Users/sheth/EP300_BPNet/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed \
        --output_prefix $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/gc_negatives \
        --flank_size 1057 \
        --neg_to_pos_ratio_train 4

# combine with all accessible regions ($ABC_peaks or below file)
gc_negatives=$OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/gc_negatives.bed # 114128
accessible_elements=$OAK/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/analysis_v2/select_p300_negatives/K562_candidate_elements.no_EP300_peak_overlap.narrowPeak  # 124713

cat $gc_negatives $accessible_elements | sort -k1,1 -k2,2n | bedtools merge -i stdin| wc -l # 201883
# remove gc negatives that overlap accessible regions...
bedtools intersect -a $gc_negatives -b $accessible_elements -v > $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/gc_negatives.filtered_accessible_regions.bed # 93361
cat $accessible_elements $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/gc_negatives.filtered_accessible_regions.bed | sort -k1,1 -k2,2n > $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/combined_negatives.bed
final_negatives=$OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/data/combined_negatives.bed

## 4) GET OPTIMAL COUNTS LOSS WEIGHT ##

bpnet-counts-loss-weight --input-data $results_dir/config/input_data.json
# 94.01

#### TRAINING ####
# 22h for fold0, 18h for the rest
# submit job for each fold
# done: 0, 1, 2, 2, 4

project_dir=$OAK/Users/sheth/EP300_BPNet
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
results_dir=$OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3
model_out_dir=$results_dir/models
input_config=$results_dir/config/input_data.json
params_config=$results_dir/config/bpnet_params.json
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes

FOLD=1
splits=$project_dir/reference/fold${FOLD}.json 

sbatch \
  --job-name=train_f${FOLD} \
  --output=slurm_logs/train_f${FOLD}.%j.txt \
  --error=slurm_logs/train_f${FOLD}.%j.txt \
  --time=22:00:00 \
  $SCRIPTS_DIR/1.1.submit_training_one_fold.sh \
  $FOLD \
  $model_out_dir \
  $splits \
  $input_config \
  $params_config \
  $genome \
  $chr_sizes


#### CV PREDICTIONS ####
# generate predictions on held-out chromosomes from each fold
# submitted

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
project_dir=$OAK/Users/sheth/EP300_BPNet
pred_data_config=$project_dir/2025_1016_p300_model_v3/config/input_data_predict.json
cv_pred_dir=$project_dir/2025_1016_p300_model_v3/predictions_cv
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
model_path_pattern="$project_dir/2025_1016_p300_model_v3/models/fold{fold}/model_split000"

ARRAY="0-4"
sbatch \
  --job-name=cvp_f%a \
  --output=log/cvp_f%a.%A.txt \
  --error=log/cvp_f%a.%A.txt \
  --array=$ARRAY \
  $SCRIPTS_DIR/2.1.submit_cv_predict.sh \
    $model_path_pattern \
    $chr_sizes \
    $genome \
    $cv_pred_dir \
    $pred_data_config


#### MEAN PREDICTIONS ####
## 1) generate predictions with models from each fold on all elements
# done

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
project_dir=$OAK/Users/sheth/EP300_BPNet
pred_data_config=$project_dir/2025_1016_p300_model_v3/config/input_data_predict.json
mean_pred_dir=$project_dir/2025_1016_p300_model_v3/predictions_mean
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
model_path_pattern="$project_dir/2025_1016_p300_model_v3/models/fold{fold}/model_split000"

ARRAY="0-4"
sbatch \
  --job-name=mp_f%a \
  --output=log/mp_f%a.%A.txt \
  --error=log/mp_f%a.%A.txt \
  --array=$ARRAY \
  $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
    $model_path_pattern \
    $chr_sizes \
    $genome \
    $mean_pred_dir \
    $pred_data_config


## 2) compute mean predictions
# done
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
project_dir=$OAK/Users/sheth/EP300_BPNet
mean_pred_dir=$project_dir/2025_1016_p300_model_v3/predictions_mean
chrom_annot_file=/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv
cols_keep="chrom,start,end,EP300.RPM,EP300_peak_overlap"

out_dir=$mean_pred_dir/all_folds; mkdir -p $out_dir
pred_h5=model_split000_predictions.h5
pred_h5_list="$mean_pred_dir/fold0/$pred_h5,$mean_pred_dir/fold1/$pred_h5,$mean_pred_dir/fold2/$pred_h5,$mean_pred_dir/fold3/$pred_h5,$mean_pred_dir/fold4/$pred_h5"

python $BPNET_DIR/utils/mean_predictions.py \
  --prediction_h5s $pred_h5_list \
  --chrom_sizes $chr_sizes \
  --output_dir $out_dir \
  --chrom_annot $chrom_annot_file \
  --chrom_cols $cols_keep \
  --generate_bigwigs \
  --negate_minus_strand

#### PREDICTION METRICS ####

conda activate tfmodisco

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
project_dir=$OAK/Users/sheth/EP300_BPNet
mean_pred_dir=$project_dir/2025_1016_p300_model_v3/predictions_mean
cv_pred_dir=$project_dir/2025_1016_p300_model_v3/predictions_cv
overlap_col="EP300_peak_overlap"

python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
  --mean_pred_dir $mean_pred_dir \
  --cv_pred_dir $cv_pred_dir \
  --overlap_col $overlap_col

#### SHAP (COUNTS ONLY) ####

# 1) run each fold on all chr, parallelize jobs across chr (chr1 requires 12h, the rest only need 5h)(ARRAY="0-23")
# folds complete: 0, 1, 2, 3, 4

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
project_dir=$OAK/Users/sheth/EP300_BPNet
pred_data_config=$project_dir/2025_1016_p300_model_v3/config/input_data_predict.json
#peaks_use=$project_dir/2025_0325_K562_BPNet/data/input_peaks.narrowPeak
peaks_use=$project_dir/reference/K562_DNase_candidate_elements.narrowPeak
shap_dir=$project_dir/2025_1016_p300_model_v3/shap
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
model_path_pattern="$project_dir/2025_1016_p300_model_v3/models/fold{fold}/model_split000"

FOLD=3
sh $SCRIPTS_DIR/3.1.submit_mean_shap_one_fold.sh \
   $FOLD \
   $model_path_pattern \
   $pred_data_config \
   $peaks_use \
   $shap_dir \
   $genome \
   $chr_sizes \
   $BPNET_DIR


## 2) combine h5 across chromosomes across folds
# done: 0, 1, 2, 3, 4

conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

SHAP_ID=shap
FOLD=3
python $BPNET_DIR/utils/merge_shap_across_chrom.py \
  --input-dir $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/${SHAP_ID}/fold${FOLD} \
  --h5-filename counts_scores.h5 \
  --output-file $OAK/Users/sheth/EP300_BPNet/2025_1016_p300_model_v3/${SHAP_ID}/fold${FOLD}/shap_counts_merged.h5

# 3) compute mean shap across folds 
# outputs: $out_dir/counts_mean_shap_scores.h5 and $out_dir/counts_peaks_valid_scores.bed
project_dir=$OAK/Users/sheth/EP300_BPNet
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
SHAP_ID=shap

shap_dir=$project_dir/2025_1016_p300_model_v3/${SHAP_ID}
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
out_dir=$project_dir/2025_1016_p300_model_v3/shap/all_folds

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
this_dir=$base_dir/2025_1016_p300_model_v3
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