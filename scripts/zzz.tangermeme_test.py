### ENV ###
# conda activate tfmodisco
# ml biology
# ml htslib/1.16

### IMPORTS ###
import os, subprocess
import h5py
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn; seaborn.set_style('whitegrid')

from bpnetlite.bpnet import BasePairNet
from bpnetlite.bpnet import ControlWrapper
from bpnetlite.bpnet import CountWrapper


from tangermeme.utils import random_one_hot
from tangermeme.predict import predict

import torch
print("PyTorch version:", torch.__version__)
print("CUDA compiled:", torch.version.cuda)
print("CUDA available at runtime:", torch.cuda.is_available())

### HELPER FUNCTIONS ###
def get_counts_model(model_path):
    model = BasePairNet.from_bpnet(model_path)
    model = CountWrapper(ControlWrapper(model))
    return model

def df_to_dict(df, key_col, value_col):
    df_uniq = df[[key_col, value_col]].drop_duplicates()
    df_map = dict(zip(df_uniq[key_col], df_uniq[value_col]))

    return df_map
    
# returns: reformatted hits_df
def read_finemo_hits(hits_path, motif_data, rename = False):
    hits = pd.read_table(hits_path)
    hits["peak_name"] = hits["chr"] + ":" + hits["start"].astype(str) + "-" + hits["end"].astype(str)

    if rename:
        # rename motif_name in hits from motif_name_orig --> motif_name
        motif_name_map = df_to_dict(motif_data, "motif_name_orig", "motif_name")
        hits["motif_name"] = [motif_name_map[m] for m in hits["motif_name"]]
        uniq_names = list(set(motif_data["motif_name"].to_list()))
        motif_idx_map = dict(zip(uniq_names, range(len(uniq_names))))
        motif_names = dict(zip(range(len(uniq_names)), uniq_names))
    else:
        # motif index, same for + and - strand
        motif_idx_map = df_to_dict(motif_data, "motif_name_orig", "motif_idx")
        motif_names = df_to_dict(motif_data, "motif_idx", "motif_name_orig")
        
    hits["motif_idx"] = [motif_idx_map[m] for m in hits["motif_name"]]

    # add other columns for compatibility
    hits["sequence_name"] = hits["peak_id"]
    hits["score"] = hits["hit_coefficient_global"]
    hits["p-value"] = 0
    hits["attribution"] = hits["hit_importance"]

    return hits, motif_names

def add_label_to_hits(hits, motif_data):
    # add index based on original names 
    motif_idx_map = df_to_dict(motif_data, "motif_name_orig", "motif_idx")
    motif_label_map = df_to_dict(motif_data, "motif_name_orig", "motif_name")
    hits["motif_idx"] = [motif_idx_map[m] for m in hits["motif_name"]]

    # label: (orig_idx) motif_renamed
    hits["motif_label"] = [f"({motif_idx_map[n]}) {motif_label_map[n]}" for n in hits["motif_name"]] # new name

    return hits

def save_hits_to_bed(hits_path, motif_data, out_prefix): 
    tsv_path = out_prefix + ".tsv"
    bed_path = out_prefix + ".bed"

    hits = pd.read_table(hits_path)
    hits = add_label_to_hits(hits, motif_data)
    hits.to_csv(tsv_path, sep="\t", index=False, header=True)
    
    hits_bed = hits[["chr", "start", "end", "motif_label", "hit_importance", "strand"]]
    bed_path = out_prefix + ".bed"
    hits_bed.to_csv(bed_path, sep="\t", index=False, header=False)
    subprocess.run(["bgzip", "-f", bed_path], check=True)
    subprocess.run(["tabix", "-p", "bed", f"{bed_path}.gz"], check=True)

    print(f"Saved and indexed hits at {bed_path}.gz!")

def compile_predictions_and_hits(pred_path, finemo_dir, motif_data, out_file):
    ## finemo peaks = pred peaks +/- 250 bp

    ## summarize hits per peak
    # read in all peaks for coordinate info
    peaks_path = os.path.join(finemo_dir, "peaks_qc.tsv") # peak_id, chr, peak_region_start
    peaks = pd.read_table(peaks_path)
    peaks["chr_start"] = peaks["chr"] + ":" + (peaks["peak_region_start"] - 250).astype(str) 
    peaks = peaks[['peak_id', 'chr_start']]

    # read in hits and summarize
    hits_path = os.path.join(finemo_dir, "hits.tsv") # peak_id
    hits = pd.read_table(hits_path)
    hits = add_label_to_hits(hits, motif_data)
    hits_smry = (
        hits.groupby('peak_id').agg(
            n_hits=('motif_name', 'count'),
            motif_list=('motif_label', lambda x: ' | '.join(x)),
            motif_strand_list=('strand', lambda x: ' | '.join(x))).reset_index()
    )
    print(hits_smry.head())

    # merge hits with peaks and fill missing values
    hits_smry = peaks.merge(hits_smry, on = "peak_id", how = "left")
    print(hits_smry.head())

    hits_smry['n_hits'] = hits_smry['n_hits'].fillna(0).astype(int)
    hits_smry[['motif_list', 'motif_strand_list']] = hits_smry[['motif_list', 'motif_strand_list']].fillna("")
    print(hits_smry.head())

    ### merge hits summary with predictions
    # read in predictions for all peaks
    pred = pd.read_table(pred_path, compression = "gzip")
    pred["chr_start"] = pred["chrom"] + ":" + pred["start"].astype(str) # for matching
    pred["peak_name"] = pred["chrom"] + ":" + pred["start"].astype(str) + "-" + pred["end"].astype(str)
    #pred = pred.drop(columns = ['chrom', 'start', 'end'])
    print(pred.head())

    # merge hits with pred
    pred = pred.merge(hits_smry, on='chr_start', how='left')
    pred = pred.drop(columns = "chr_start")
    print(pred.head())

    print(f"Predictions without peaks from finemo: {pred['n_hits'].isnull().sum()}")

    # save
    pred.to_csv(out_file, sep="\t", index=False)


def torch_cosine_normalize(counts):
    # equivalent to cooccurrence from finemo
    # normalized values = cosine similarity of binary motif presence; 1 meaning they always cooccur and 0 meaning they never do

    diag = torch.diag(counts)
    inv_sqrt = 1 / torch.sqrt(diag)
    norm = inv_sqrt[:, None] * inv_sqrt[None, :]
    result = counts * norm

    return result

def plot_cooccurence(hits, motif_names, out_file):
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import pdist, squareform
    from tangermeme.annotate import count_annotations
    from tangermeme.annotate import pairwise_annotations

    cooccur_counts = pairwise_annotations((hits['sequence_name'].astype(int), hits['motif_idx']), symmetric = True)
    cooccur_norm = torch_cosine_normalize(cooccur_counts)

    motif_labels = list(motif_names.values())

    # Step 3: Plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    # Raw counts
    m0 = ax[0].imshow(cooccur_counts, cmap='Blues', interpolation = "nearest", vmin=0, vmax=len(hits) / 15)
    ax[0].set_xticks(range(len(motif_labels)))
    ax[0].set_xticklabels(motif_labels, rotation=90, fontsize=6)
    ax[0].set_yticks(range(len(motif_labels)))
    ax[0].set_yticklabels(motif_labels, fontsize=6)
    ax[0].set_title("Co-occurrence peak counts")
    ax[0].grid(False)
    cbar0 = fig.colorbar(m0, ax=ax[0], extend='max', shrink = 0.6)

    # Normalized
    m1 = ax[1].imshow(cooccur_norm, cmap='Purples', interpolation = "nearest", vmin=0, vmax=1)
    ax[1].set_xticks(range(len(motif_labels)))
    ax[1].set_xticklabels(motif_labels, rotation=90, fontsize=6)
    ax[1].set_yticks(range(len(motif_labels)))
    ax[1].set_yticklabels(motif_labels, fontsize=6)
    ax[1].set_title("Normalized co-occurrence")
    ax[1].grid(False)
    cbar1 = fig.colorbar(m1, ax=ax[1], shrink = 0.6)

    plt.tight_layout()
    plt.savefig(out_file)
    plt.show()


### FILE PATHS AND PARAMS ###
folds = ["fold0", "fold1", "fold2", "fold3", "fold4"]

out_dir = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/tangermeme_all_chr"
os.makedirs(out_dir, exist_ok=True)

model_dir = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models_h5"
model_paths = [os.path.join(model_dir, f, "ENCSR000EGE_model.h5") for f in folds]

finemo_dir = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/finemo_all_chr"
finemo_hits_path = os.path.join(finemo_dir, "hits.tsv")
motif_data_path = os.path.join(finemo_dir, "motif_data_renamed.tsv")

pred_dir = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/predictions_and_metrics_all_chr"
pred_counts_path = os.path.join(pred_dir, "all_folds/mean_predictions_counts.tsv.gz")

motif_data = pd.read_table(motif_data_path)
motif_data["motif_idx"] = [x // 2 for x in motif_data["motif_id"]] # regardless of strand

compile_predictions_and_hits(pred_counts_path, finemo_dir, motif_data, os.path.join(finemo_dir, "hits_per_peak.with_predictions.tsv.gz"))

save_hits_to_bed(finemo_hits_path, motif_data, os.path.join(finemo_dir, "hits_renamed"))

# hits, motif_names = read_finemo_hits(finemo_hits_path, motif_data, rename = True)
# motif_id = "grouped_motifs"

# out_file = os.path.join(out_dir, f"cooccurence.no_filter.{motif_id}.pdf")
# plot_cooccurence(hits, motif_names, out_file)

# hits = hits.loc[hits["attribution"] > 0.25]
# out_file = os.path.join(out_dir, f"cooccurence.attrib_over_25.{motif_id}.pdf")
# plot_cooccurence(hits, motif_names, out_file)

# hits = hits.loc[hits["attribution"] > 0.5]
# out_file = os.path.join(out_dir, f"cooccurence.attrib_over_50.{motif_id}.pdf")
# plot_cooccurence(hits, motif_names, out_file)


if False:
    hits = hits.sort_values("attribution")

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].plot(hits['attribution'].values)
    ax[0].set_xlabel("Hit (Sorted by Contribution)")
    ax[0].set_ylabel("Contribution Sum")
    ax[0].set_title("Contribution Score per Hit")

    hits = hits.sort_values("score")
    ax[1].plot(hits['score'].values)
    ax[1].set_xlabel("Hit (Sorted by Score)")
    ax[1].set_ylabel("Score Sum")
    ax[1].set_title("Score per Hit")

    plt.show()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "attributions.pdf"))





# X = random_one_hot((10, 4, 2114), random_state = 0)
# X = X.float()

# model = get_counts_model(model_paths[0])

# # example without gpu
# predict(model, X, device = "cpu", verbose = True)


## imports
import os
import numpy as np
import pandas as pd
import bpnetlite.bpnet as bpnet
import torch
import tangermeme
from tangermeme.io import extract_loci

## file paths
my_oak = "/oak/stanford/groups/engreitz/Users/sheth"
this_dir = os.path.join(my_oak, "EP300_BPNet","2025_0421_tangermeme"); os.makedirs(this_dir, exist_ok = True)

model_path = os.path.join("")
gc_neg_path = os.path.join(my_oak, "EP300_BPNet",  "2025_0325_K562_BPNet/data/gc_negatives.bed")

hg38_fasta = os.path.join(my_oak, "hg38_resources", "hg38.fa")


test = bpnet.from_bpnet(model_path)

print(test)

gc_neg = pd.read_csv(gc_neg_path, sep="\t", usecols=(0, 1, 2), names=['chrom', 'start', 'end'])
gc_neg_seq = extract_loci(gc_neg, hg38_fasta, verbose = True, in_window = 2114).float()
print(gc_neg_seq.shape)

X = gc_neg_seq
X = X[X.sum(dim=(1, 2)) == X.shape[-1]]
X.shape