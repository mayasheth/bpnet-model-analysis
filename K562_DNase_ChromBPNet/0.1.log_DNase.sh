
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/K562_DNase_ChromBPNet


### DOWNLOAD MODELS ###
wget https://www.encodeproject.org/files/ENCFF574YLK/@@download/ENCFF574YLK.tar.gz
tar -xvzf ENCFF574YLK.tar.gz # extracted to K562_DNase_ChromBPNet/models


# extract tar files for each fold model to $this_dir/models/fold${n_fold}/
n_fold=4
model_tar_path=$this_dir/models/fold_${n_fold}/model.chrombpnet_nobias.fold_${n_fold}.ENCSR000EOT.tar
target_model_dir=$this_dir/models/fold${n_fold}
model_path_dir=$this_dir/models/fold${n_fold}/chrombpnet_wo_bias
(
  mkdir -p "$target_model_dir"  && cd "$target_model_dir" && tar -xvf "$model_tar_path"
)

# actually we want the h5 file, which we have examples of how to load?
# move some stuff around...
n_fold=4
model_path_orig=$this_dir/models/fold_${n_fold}/model.chrombpnet_nobias.fold_${n_fold}.ENCSR000EOT.h5
model_path_target=$this_dir/models/fold${n_fold}/model.chrombpnet_nobias.h5

mv $model_path_orig $model_path_target


### CALL HITS FROM DNASE-SEQ CHROMBPNET MODEL MODISCO RESULTS ###
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/K562_DNase_ChromBPNet
data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase
out_dir=$this_dir/finemo; mkdir -p $out_dir

# 1) some preprocessing of downloaded files
# get regions for modisco input
conda activate ucsc_tools

base_dir=$OAK/Users/sheth/EP300_BPNet
data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase

# wrong regions (170k)
bigBedToBed  $data_dir/ENCFF229FNM.bigBed $data_dir/ENCFF229FNM.bed
sh $base_dir/scripts/reformat_input_peaks.sh $data_dir/ENCFF229FNM.bed $data_dir/ENCFF229FNM.narrowPeak 0

# right regions (190k)
peaks_all_base=$data_dir/ENCFF739NDX/logs.seq_contrib.counts.ENCSR000EOT/logs.seq_contrib.counts.input_regions.modisco.ENCSR000EOT
zcat ${peaks_all_base}.bed.gz > ${peaks_all_base}.bed
sh $base_dir/scripts/reformat_input_peaks.sh ${peaks_all_base}.bed ${peaks_all_base}.narrowPeak 0

# 2) run tf-modisco to get html report and look at motifs
conda activate tfmodisco

data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase
base_dir=$OAK/Users/sheth/EP300_BPNet
motif_ref=$base_dir/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

modisco report -i $data_dir/ENCFF404TKA/counts/tfmodisco.raw_output.counts.ENCSR000EOT.hd5 \
    -o $data_dir/ENCFF404TKA/counts/counts_report/ \
    -s $data_dir/ENCFF404TKA/counts/counts_report/ \
    -m $motif_ref

# 3) run finemo

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
out_dir=$this_dir/finemo_DNase
shap_counts_h5=$data_dir/ENCFF739NDX/seq_contrib.counts.fold_mean.modisco_input.ENCSR000EOT.h5
peaks_use=$data_dir/ENCFF739NDX/logs.seq_contrib.counts.ENCSR000EOT/logs.seq_contrib.counts.input_regions.modisco.ENCSR000EOT.narrowPeak
modisco_counts_h5=$data_dir/ENCFF404TKA/counts/tfmodisco.raw_output.counts.ENCSR000EOT.hd5
modisco_width=500

## 1) call hits
sbatch $SCRIPTS_DIR/6.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $modisco_counts_h5 \
  $modisco_width


#sbatch scripts/3.3.1.finemo_call_hits_dnase.sh

# 4) rerun reporting with relabeled motifs
conda activate finemo_gpu

base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
out_dir=$this_dir/finemo_DNase
regions_npz=$out_dir/contribution_regions.npz
motif_labels=$out_dir/motif_label_map.txt

finemo report \
        --regions $regions_npz \
        --hits $out_dir \
        --out-dir $out_dir \
        --modisco-region-width 400 \
        --motif-labels $motif_labels



### RERUN MODISCO ###
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
base_dir=$OAK/Users/sheth/EP300_BPNet
this_dir=$base_dir/2025_0517_official_EP300_K562_model
data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase
shap_counts_h5=$data_dir/ENCFF739NDX/seq_contrib.counts.fold_mean.modisco_input.ENCSR000EOT.h5
motif_reference=$base_dir/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
max_seqlets=250000
modisco_width=400

## longer motifs (running)
modisco_counts_dir=$this_dir/modisco_DNase/max_seqlets_250k_20_5_10
trim_size=20
initial_flank_to_add=5
final_flank_to_add=10

## shorter motifs (running?)
modisco_counts_dir=$this_dir/modisco_DNase/max_seqlets_250k_30_10_0
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


### RERUN FINEMO ON NEW MODISCO RUNS ###

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
data_dir=$OAK/Users/sheth/Data/ENCODE/K562/BPNet/DNase
shap_counts_h5=$data_dir/ENCFF739NDX/seq_contrib.counts.fold_mean.modisco_input.ENCSR000EOT.h5
peaks_use=$data_dir/ENCFF739NDX/logs.seq_contrib.counts.ENCSR000EOT/logs.seq_contrib.counts.input_regions.modisco.ENCSR000EOT.narrowPeak

## official model version (width=500?) (done)
modisco_width=500
modisco_id=encode_run
out_dir=$this_dir/finemo_DNase/$modisco_id
modisco_counts_h5=$OAK/Users/sheth/EP300_BPNet/2025_0407_MoDISCo/ENCSR000EGE_trim20_flank5_0/modisco_results_counts.h5

## shorter motifs (to run - took to long/didn't work..)
modisco_width=400
modisco_id=max_seqlets_250k_30_10_0
out_dir=$this_dir/finemo_DNase/$modisco_id
modisco_counts_h5=$this_dir/modisco_DNase/$modisco_id/counts_scores.h5

## longer motifs (to run took to long/didn't work)
modisco_width=400
modisco_id=max_seqlets_250k_20_5_10
out_dir=$this_dir/finemo_DNase/$modisco_id
modisco_counts_h5=$this_dir/modisco_DNase/$modisco_id/counts_scores.h5

## 1) call hits
sbatch $SCRIPTS_DIR/6.1.submit_finemo.sh \
  $out_dir \
  $shap_counts_h5 \
  $peaks_use \
  $modisco_counts_h5 \
  $modisco_width


### REFORMAT FINEMO HITS ###

SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
this_dir=$OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model
motif_file=$OAK/Users/sheth/EP300_BPNet/reference/motif_annotations.tsv

target_name="dnase"
modisco_id="encode_run"

finemo_dir=$this_dir/finemo_DNase/${modisco_id}
out_dir=$finemo_dir/annotated_motifs

python $SCRIPTS_DIR/6.2.format_finemo_hits.py \
  --finemo_dir $finemo_dir \
  --pred "NONE" \
  --motifs $motif_file \
  --target_name "$target_name" \
  --modisco_id $modisco_id \
  --out_dir $out_dir