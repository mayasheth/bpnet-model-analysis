#### GOAL ####
# intepret official ENCODE p300 BPNet model 

#### BPNET ENV #### 
conda activate bpnet_37
module load cuda/11.1.1
module load cudnn/8.1.1.33

#### FILE PATHS ####

base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model

ABC_peaks=$OAK/Users/sheth/ENCODE_rE2G_main/ENCODE_rE2G/results/2025_0313_ENCODE_K562_DNase/Stam_ENCFF860XAE_ENCFF205FNC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
chr_list=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.txt
genome=$OAK/Users/sheth/hg38_resources/hg38.fa

peaks_use=$base_dir/2025_0325_K562_BPNet/data/input_peaks.narrowPeak
pred_data_config=$this_dir/config/input_data_predict.json

#### GET MODELS ####

# 1) download models for each fold from mitra
wget -r -np -nH --cut-dirs=4 -R "index.html*" https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold0/ENCSR000EGE/
wget -r -np -nH --cut-dirs=4 -R "index.html*" https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold1/ENCSR000EGE/
wget -r -np -nH --cut-dirs=4 -R "index.html*" https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold2/ENCSR000EGE/
wget -r -np -nH --cut-dirs=4 -R "index.html*" https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold3/ENCSR000EGE/
wget -r -np -nH --cut-dirs=4 -R "index.html*" https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold4/ENCSR000EGE/

# 2) extract tar files to $this_dir/models/release_run_1/fold${n_fold}/ENCSR000EGE/ENCSR000EGE_split000
this_shap_dir=$shap_dir/fold${n_fold}; mkdir -p $this_shap_dir
model_tar_path=$this_dir/models/release_run_1/fold${n_fold}/ENCSR000EGE/ENCSR000EGE_split000.tar
model_path_dir=$this_dir/models/release_run_1/fold${n_fold}/ENCSR000EGE
(
  cd "$model_path_dir"  && \
  tar -xvf "$model_tar_path"
)


#### PREPROCESSING ####

# reformat candidate elements to narrowPeak format
sh $base_dir/scripts/reformat_input_peaks.sh $ABC_peaks $peaks_use 1057


#### CV PREDICTIONS ####

# predict on candidate elements (submits one job per fold)

base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
pred_dir=$this_dir/predictions_cv; mkdir -p $pred_dir

sbatch scripts/1.1.submit_predict_cv.sh



#### MEAN PREDICTIONS #### 

# 1) predict with all folds on all chr
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
pred_dir=$this_dir/predictions_mean; mkdir -p $pred_dir

sbatch scripts/1.2.submit_predict_all_chr.sh

# 2) average as follows:
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
project_dir=$OAK/Users/sheth/EP300_BPNet
mean_pred_dir=$project_dir/2025_0517_official_EP300_K562_model/predictions_mean
chrom_annot_file=/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv
cols_keep="chrom,start,end,EP300.RPM,EP300_peak_overlap"

out_dir=$mean_pred_dir/all_folds; mkdir -p $out_dir
pred_h5=ENCSR000EGE_split000_predictions.h5
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
mean_pred_dir=$project_dir/2025_0517_official_EP300_K562_model/predictions_mean
cv_pred_dir=$project_dir/2025_0517_official_EP300_K562_model/predictions_cv
peaks=$project_dir/reference/ENCSR000EGE_peaks_inliers.narrowPeak
overlap_col="EP300_peak_overlap"

# Note: corrected metric — log1p(expm1(+strand) + expm1(-strand)) instead of
# log1p(+strand) + log1p(-strand). Requires h5 files with coords/ group.
# Peak overlap computed directly via interval overlap against p300 peaks (not chromatin annot join).
python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
  --cv-pred-dir $cv_pred_dir \
  --mean-pred-dir $mean_pred_dir \
  --peaks $peaks \
  --overlap-col $overlap_col

# Inter-replicate correlation stratified by p300 peak overlap
chromatin_annot=$OAK/Users/sheth/TF_analysis/2026_0603_EP300_reps/EnhancerList.chromatin_annotations.tsv
python $SCRIPTS_DIR/compute_replicate_correlations.py \
  --chromatin-annot $chromatin_annot \
  --rep1-col EP300_R1.RPM \
  --rep2-col EP300_R2.RPM \
  --overlap-col $overlap_col \
  --output $project_dir/2025_0517_official_EP300_K562_model/replicate_correlations.tsv


#### SHAP ####
# run shap for all folds on all chrom to get mean shap scores (counts only)

mkdir -p $this_dir/shap_peaks

## 1) run each fold on all chr, parallelize jobs across chr
# folds done: 0, 1, 2, 3, 4

ARRAY="1-23"
FOLD=4
sbatch \
  --job-name=shap_f${FOLD}_chr \
  --output=log/shap_f${FOLD}_chr.%a.%A.txt \
  --error=log/shap_f${FOLD}_chr.%a.%A.txt \
  --array=$ARRAY \
  scripts/2.2.submit_mean_shap.sh $FOLD

ARRAY="0-1"
FOLD=4
sbatch \
  --job-name=shap_f${FOLD}_chr \
  --output=log/shap_f${FOLD}_chr.%a.%A.txt \
  --error=log/shap_f${FOLD}_chr.%a.%A.txt \
  --array=$ARRAY \
  --time=10:00:00 \
  scripts/2.2.submit_mean_shap.sh $FOLD

## 2) combine h5 across chromosomes across folds
# folds completed: 0, 1, 2, 3, 4

conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

FOLD=4
python $BPNET_DIR/utils/merge_shap_across_chrom.py \
  --input-dir $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/shap/fold${FOLD} \
  --h5-filename counts_scores.h5 \
  --output-file $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/shap/fold${FOLD}/shap_counts_merged.h5

## 3) compute mean shap across folds 
# outputs: $out_dir/counts_mean_shap_scores.h5 and $out_dir/counts_peaks_valid_scores.bed
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

shap_dir=$this_dir/shap
out_dir=$shap_dir/all_folds; mkdir -p $out_dir
shap_h5=shap_counts_merged.h5
shap_h5_list="$shap_dir/fold0/$shap_h5,$shap_dir/fold1/$shap_h5,$shap_dir/fold2/$shap_h5,$shap_dir/fold3/$shap_h5,$shap_dir/fold4/$shap_h5"

python $BPNET_DIR/utils/mean_shap_plus_peaks.py \
  --counts_shaps $shap_h5_list \
  --output_dir $out_dir 

## 4) make merged shap bw
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes

base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
out_dir=$this_dir/shap/all_folds

merged_shap_h5=$out_dir/counts_mean_shap_scores.h5
merged_peaks=$out_dir/counts_peaks_valid_scores.bed

python $BPNET_DIR/utils/run_importance_hdf5_to_bigwig.py \
  --hdf5_path $merged_shap_h5 \
  --regions_path $merged_peaks \
  --chrom_sizes $chr_sizes \
  --outfile $out_dir/counts_scores.bw \
  --outstats $out_dir/counts_scores.stats.txt

#### SHAP FOR PEAKS ONLY
# 1) run each fold on all chr, parallelize jobs across chr (chr1 requires 12h, the rest only need 5h)(ARRAY="0-23")
# folds complete: 0, 1, 2, 3, 4

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
project_dir=$OAK/Users/sheth/EP300_BPNet
pred_data_config=$project_dir/2025_0517_official_EP300_K562_model/config/input_data_predict.json
peaks_use=$project_dir/reference/ENCSR000EGE_peaks_inliers.narrowPeak
shap_dir=$project_dir/2025_0517_official_EP300_K562_model/shap_peaks
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
model_path_pattern="$project_dir/2025_0517_official_EP300_K562_model/models/release_run_1/fold{fold}/ENCSR000EGE/ENCSR000EGE_split000"

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

SHAP_ID=shap_peaks
FOLD=4
python $BPNET_DIR/utils/merge_shap_across_chrom.py \
  --input-dir $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/${SHAP_ID}/fold${FOLD} \
  --h5-filename counts_scores.h5 \
  --output-file $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/${SHAP_ID}/fold${FOLD}/shap_counts_merged.h5

## 3) compute mean shap across folds 
# outputs: $out_dir/counts_mean_shap_scores.h5 and $out_dir/counts_peaks_valid_scores.bed
project_dir=$OAK/Users/sheth/EP300_BPNet
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
SHAP_ID=shap_peaks

shap_dir=$project_dir/2025_0517_official_EP300_K562_model/${SHAP_ID}
out_dir=$shap_dir/all_folds; mkdir -p $out_dir
shap_h5=shap_counts_merged.h5
shap_h5_list="$shap_dir/fold0/$shap_h5,$shap_dir/fold1/$shap_h5,$shap_dir/fold2/$shap_h5,$shap_dir/fold3/$shap_h5,$shap_dir/fold4/$shap_h5"

python $BPNET_DIR/utils/mean_shap_plus_peaks.py \
  --counts_shaps $shap_h5_list \
  --output_dir $out_dir 


#### MODISCO (COUNTS ONLY) ####
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
shap_counts_h5=$this_dir/shap/all_folds/counts_mean_shap_scores.h5
motif_reference=$base_dir/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
max_seqlets=250000
modisco_width=400

## longer motifs
modisco_counts_dir=$this_dir/modisco/max_seqlets_250k_20_5_10
trim_size=20
initial_flank_to_add=5
final_flank_to_add=10

## shorter motifs
modisco_counts_dir=$this_dir/modisco/max_seqlets_250k_30_10_0
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


#### MODISCO (COUNTS ONLY) (IN PEAKS)  ####

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
motif_reference=$base_dir/reference/MotifCompendium-Database-Human.meme.txt
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
sbatch $SCRIPTS_DIR/4.1.submit_counts_modisco.sh \
    $modisco_counts_dir \
    $shap_counts_h5 \
    $max_seqlets \
    $trim_size \
    $initial_flank_to_add \
    $final_flank_to_add \
    $motif_reference \
    $modisco_width


#### FINEMO ####

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
shap_counts_h5=$this_dir/shap/all_folds/counts_mean_shap_scores.h5
peaks_use=$this_dir/shap/all_folds/counts_peaks_valid_scores.bed

## official model version (width=500?) (done)
modisco_width=500
modisco_id=encode_run
out_dir=$this_dir/finemo/$modisco_id
modisco_counts_h5=$OAK/Users/sheth/EP300_BPNet/2025_0407_MoDISCo/ENCSR000EGE_trim20_flank5_0/modisco_results_counts.h5

## shorter motifs (to run)
modisco_width=400
modisco_id=max_seqlets_250k_30_10_0
out_dir=$this_dir/finemo/$modisco_id
modisco_counts_h5=$this_dir/modisco/$modisco_id/counts_scores.h5

## longer motifs (to run)
modisco_width=400
modisco_id=max_seqlets_250k_20_5_10
out_dir=$this_dir/finemo/$modisco_id
modisco_counts_h5=$this_dir/modisco/$modisco_id/counts_scores.h5

## 1) call hits
sbatch $SCRIPTS_DIR/6.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $modisco_counts_h5 \
  $modisco_width

## 2) annotate motifs to create config file...

## 3) call hits with curated motifs
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
modisco_width=400
SHAP_ID=shap
FINEMO_ID=finemo

shap_counts_h5=$this_dir/${SHAP_ID}/all_folds/counts_mean_shap_scores.h5
peaks_use=$this_dir/${SHAP_ID}/all_folds/counts_peaks_valid_scores.bed
peak_width=500
out_dir=$this_dir/${FINEMO_ID}/pkw_${peak_width}_curated_motifs_v2
modisco_path=$OAK/Users/sheth/EP300_BPNet/reference/curated_motif_data_for_finemo_v2

sbatch $SCRIPTS_DIR/5.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $peak_width \
  $modisco_path \
  $modisco_width


## 4) format motifs and such
conda activate analysis

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
pred_file=$this_dir/predictions_mean/all_folds/mean_predictions_counts.tsv.gz
motif_file=$OAK/Users/sheth/EP300_BPNet/reference/motif_annotations.tsv
target_name="p300 v1"

# 30 10 0
run_id="250k_30_10_0"
peak_width=400
finemo_dir=$this_dir/finemo/max_seqlets_$run_id
out_dir=$finemo_dir/annotated_motifs

# 20 5 10
run_id="250k_20_5_10"
peak_width=400
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
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
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


#### FIMO FOR MOTIFS ####

## 1) make fasta from candidate elements and p300 inliers (named chr:start-end)
candidate_elements=$OAK/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak
p300_peaks=$OAK/Users/sheth/EP300_BPNet/reference/ENCSR000EGE_peaks_inliers.narrowPeak
genome=$OAK/Users/sheth/hg38_resources/hg38.fa
ml biology bedtools

bedtools getfasta -fi $genome -bed $candidate_elements > $OAK/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.fa
bedtools getfasta -fi $genome -bed $p300_peaks > $OAK/Users/sheth/EP300_BPNet/reference/ENCSR000EGE_peaks_inliers.fa

## 2) prepare meme file of relevant motifs

## 3) run FIMO
singularity pull docker://memesuite/memesuite:latest # singularity-izes the docker image, DONE

singularity shell docker://memesuite/memesuite:latest
meme=$OAK/Users/sheth/EP300_BPNet/reference/MotifCompendium-test.meme.txt
fasta=$OAK/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.fa
out_dir=$OAK/Users/sheth/EP300_BPNet/FIMO/2025_1112_test
fimo --oc $out_dir --qv-thresh --thresh 1e-3 --max-stored-scores 100000 $meme $fasta 


#### FINEMO TF CORRELATION ANALYSIS (2026-06-10) ####

this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
finemo_dir=$this_dir/finemo/pkw_500_curated_motifs
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts

## Compute Spearman correlations between each motif's per-peak hit activity
## and every TF ChIP-seq signal RPM in the K562 enhancer activity feature file.
## conda activate analysis
# sbatch $SCRIPTS_DIR/submit_finemo_tf_correlation.sh
# Output: $finemo_dir/tf_correlation/correlations.tsv.gz
# Output: $finemo_dir/tf_correlation/top10_per_motif.pdf
# Output: $finemo_dir/tf_correlation/correlation_density.pdf

## Compute fraction of total |SHAP| signal explained by FiNeMo motif hits
## conda activate analysis
# sbatch $SCRIPTS_DIR/submit_finemo_explained_fraction.sh
# Results: overall 8.4%, p300+ 11.4%, p300− 7.6%
# Output: $finemo_dir/annotated_motifs/motif_explained_fraction.pdf
# Output: $finemo_dir/annotated_motifs/motif_explained_fraction.motifs_only.pdf

## Correlate all 522 TF ChIP RPM signals with EP300 RPM (no motifs involved)
## conda activate analysis
# sbatch $SCRIPTS_DIR/submit_activity_ep300_correlation.sh
# Top: TBL1XR1 r=0.784, JUND r=0.774, RCOR1 r=0.737, SMARCA4 r=0.695
# Output: $finemo_dir/activity_ep300_correlation/top_corr_ep300.pdf
# Output: $finemo_dir/activity_ep300_correlation/correlations.tsv.gz

## 3-motif composite figure (CTCF / AP-1 / GATA) with p300 interaction markers
## conda activate analysis
# sbatch $SCRIPTS_DIR/submit_composite_figure.sh
# Output: $finemo_dir/tf_correlation/composite_motif_tf_figure.pdf
