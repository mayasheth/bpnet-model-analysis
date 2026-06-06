import numpy as np
import os
import tensorflow as tf
from argparse import Namespace
from bpnet.cli.predict import predict


## general file paths
MY_OAK = "/oak/stanford/groups/engreitz/Users/sheth"
base_dir = os.path.join(MY_OAK, "EP300_BPNet")
this_dir = os.path.join(base_dir, "2025_0517_official_EP300_K562_model")
pred_dir = os.path.join(this_dir, "predictions_and_metrics"); os.makedirs(pred_dir, exist_ok = True)

peaks_use = os.path.join(base_dir, "2025_0325_K562_BPNet/data/input_peaks.narrowPeak")
chr_sizes = os.path.join(MY_OAK, "hg38_resources/GRCh38.main.chrom.sizes")
genome = os.path.join(MY_OAK, "hg38_resources/hg38.fa")
pred_data_config= os.path.join(this_dir, "config/input_data_predict.json")

n_fold = 0
model_path = f"/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models/release_run_1/fold{n_fold}/ENCSR000EGE/ENCSR000EGE_model/ENCSR000EGE_split000"
test_chr = "chr1 chr3 chr6"
this_pred_dir = os.path.join(pred_dir, f"fold{n_fold}"); os.makedirs(this_pred_dir, exist_ok = True)
this_shap_dir = os.path.join(shap_dir, f"fold{n_fold}"); os.makedirs(this_shap_dir, exist_ok = True)


# Define dimensions
num_sequences = 2 # or any number of test sequences
input_seq_length = 2114
output_length = 1000
num_strands = 2
num_bases = 4 # A, C, G, T

# Create dummy zero arrays
sequence = np.zeros((num_sequences, input_seq_length, num_bases), dtype=np.float32)

profile_bias_input = np.zeros((num_sequences, output_length, num_strands), dtype=np.float32)
counts_bias_input = np.zeros((num_sequences, num_strands), dtype=np.float32)

# Load model from the untarred folder
n_fold = 0
model_path = f"/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models/release_run_1/fold{n_fold}/ENCSR000EGE/ENCSR000EGE_model/ENCSR000EGE_split000"
model = tf.saved_model.load(model_path)

# Run inference
predictions = model.signatures['serving_default'](**{
  'profile_bias_input_0': profile_bias_input,
  'counts_bias_input_0': counts_bias_input,
  'sequence': sequence
})