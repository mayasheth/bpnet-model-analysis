#!/bin/bash

target=GATA1
experiment_id=ENCSR000EWM

out_dir=$OAK/Users/sheth/EP300_BPNet/K562_${target}_BPNet
mkdir -p $out_dir

bw_plus=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/${experiment_id}/${experiment_id}_plus.bigWig
bw_minus=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/${experiment_id}/${experiment_id}_minus.bigWig
peaks=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/peaks_inliers.bed.gz
model_dirs=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold{0,1,2,3,4}/${experiment_id}/${experiment_id}_model/${experiment_id}_split000

# Create directory structure
mkdir -p $out_dir/data
mkdir -p $out_dir/models/fold{0,1,2,3,4}

# Download bigWig files
echo "Downloading bigWig files..."
wget -O $out_dir/data/plus.bigWig $bw_plus
wget -O $out_dir/data/minus.bigWig $bw_minus

# Download peaks file
echo "Downloading peaks file..."
wget -O $out_dir/data/peaks.bed.gz $peaks

# Download model directories for each fold
echo "Downloading model directories..."
for fold in {0,1,2,3,4}; do
    echo "Downloading fold $fold..."
    model_url=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold${fold}/${experiment_id}/${experiment_id}_model/${experiment_id}_split000
    wget -r -np -nH --cut-dirs=6 -P $out_dir/models/fold${fold}/ $model_url
    # Move the downloaded content to the expected location
    if [ -d "$out_dir/models/fold${fold}/${experiment_id}_split000" ]; then
        mv $out_dir/models/fold${fold}/${experiment_id}_split000 $out_dir/models/fold${fold}/model_split000
    fi
done

echo "Download complete!"
echo "Data saved to: $out_dir/data/"
echo "Models saved to: $out_dir/models/fold{0,1,2,3,4}/model_split000"