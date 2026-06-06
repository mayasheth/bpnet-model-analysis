#!/bin/bash

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <experiment_id> <out_dir>"
    echo "Example: $0 ENCSR000EWM /path/to/output"
    exit 1
fi

# Parse command line arguments
experiment_id=$1
out_dir=$2

# Save a copy of this script to the output directory
script_name="download.sh"
script_path="$(readlink -f "$0")"
mkdir -p "$out_dir"
cp "$script_path" "$out_dir/$script_name"
echo "Script saved to: $out_dir/$script_name"

# Set up URLs and paths
bw_plus=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/${experiment_id}/${experiment_id}_plus.bigWig
bw_minus=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/${experiment_id}/${experiment_id}_minus.bigWig
peaks=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/processed_data/peaks_inliers.bed.gz

# Create directory structure
echo "Creating directory structure..."
mkdir -p $out_dir/data
mkdir -p $out_dir/models/fold{0,1,2,3,4}
mkdir -p $out_dir/config

# Download bigWig files
echo "Downloading bigWig files..."
wget -O $out_dir/data/plus.bigWig $bw_plus
wget -O $out_dir/data/minus.bigWig $bw_minus

# Download peaks file
echo "Downloading peaks file..."
wget -O $out_dir/data/peaks_inliers.bed.gz $peaks

# Convert peaks to narrowPeak
zcat $out_dir/data/peaks_inliers.bed.gz > $out_dir/data/peaks_inliers.bed
sh $OAK/Users/sheth/EP300_BPNet/scripts/0.2.bed_to_narrowPeak.sh $out_dir/data/peaks_inliers.bed $out_dir/data/peaks_inliers.narrowPeak 1057
rm $out_dir/data/peaks_inliers.bed

# Download model directories for each fold
echo "Downloading model directories..."
for fold in {0,1,2,3,4}; do
    echo "Downloading fold $fold..."
    model_url=https://mitra.stanford.edu/kundaje/oak/vir/tfatlas/models/release_run_1/fold${fold}/${experiment_id}/${experiment_id}_model/${experiment_id}_split000/
    wget -r --reject "index.html*" -np -nH --cut-dirs=9 -P $out_dir/models/fold${fold}/ $model_url
    # Move the downloaded content to the expected location
    if [ -d "$out_dir/models/fold${fold}/${experiment_id}_split000" ]; then
        mv $out_dir/models/fold${fold}/${experiment_id}_split000 $out_dir/models/fold${fold}/model_split000
    fi
done

# Create input_data.json config file
echo "Creating input_data.json config file..."
cat > $out_dir/config/input_data.json << EOF
{
    "0": {
        "signal": {
            "source": ["$out_dir/data/plus.bigWig",
                       "$out_dir/data/minus.bigWig"]
        },
        "loci": {
            "source": ["$out_dir/data/peaks_inliers.narrowPeak"]
        },
        "bias": {
            "source": ["/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_plus.bigWig",
                       "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_minus.bigWig"],
            "smoothing": [null, null]
        }
    }
}
EOF

# Create input_data_predict.json config file
    # K562 candidate elements previously created from $OAK/Users/sheth/ENCODE_rE2G_main/ENCODE_rE2G/results/2025_0313_ENCODE_K562_DNase/Stam_ENCFF860XAE_ENCFF205FNC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
    # 
echo "Creating input_data.json config file..."
cat > $out_dir/config/input_data.json << EOF
{
    "0": {
        "signal": {
            "source": ["$out_dir/data/plus.bigWig",
                       "$out_dir/data/minus.bigWig"]
        },
        "loci": {
            "source": ["$OAK/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak"]
        },
        "bias": {
            "source": ["/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_plus.bigWig",
                       "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_control_minus.bigWig"],
            "smoothing": [null, null]
        }
    }
}
EOF

echo "Download complete!"
echo "Data saved to: $out_dir/data/"
echo "Models saved to: $out_dir/models/fold{0,1,2,3,4}/model_split000"