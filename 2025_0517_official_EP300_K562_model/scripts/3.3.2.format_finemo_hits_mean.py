import os
import numpy as np
import pandas as pd
import subprocess

# file paths and set up
BASE_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model"
OUT_DIR = os.path.join(BASE_DIR, "finemo_DNase")


input_peaks_file = os.path.join("/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/BPNet/DNase",
    "ENCFF739NDX/logs.seq_contrib.counts.ENCSR000EOT/logs.seq_contrib.counts.input_regions.modisco.ENCSR000EOT.narrowPeak")
finemo_hits_file = os.path.join(OUT_DIR, "hits.tsv")

NARROWPEAK_SCHEMA = ["chr", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand", "peak_signal", "peak_pval", "peak_qval", "peak_summit"]

def save_renamed_hits(hits_path, mapping_path)

# goal: 
assembled_hits = []
formatted_peaks = []

for this_fold in FOLDS: 
    print(f"Working on {this_fold}...")

    # load peaks
    peaks = pd.read_table(os.path.join(BASE_DIR, "finemo", this_fold, "input_peaks.narrowPeak"),
                        names = NARROWPEAK_SCHEMA, header = None, comment = "#")
    peaks["peak_id"] = range(len(peaks))  # index starting with 0
    peaks["peak_name"] = peaks["chr"] + ":" + peaks["peak_start"].astype(str) + "-" + peaks["peak_end"].astype(str)
    peaks["fold"] = this_fold
    
    # add true and predicted counts
    for col in ["true_logcounts_0", "true_logcounts_1", "pred_logcounts_0", "pred_logcounts_1"]:
        peaks[col] = np.load(os.path.join(BASE_DIR, "predictions_and_metrics", this_fold, f"{col}.npy"))
    peaks["true_logcounts_total"] = peaks["true_logcounts_0"] + peaks["true_logcounts_1"]
    peaks["pred_logcounts_total"] = peaks["pred_logcounts_0"] + peaks["pred_logcounts_1"]

    # map motif names for finemo hits
    finemo_hits = pd.read_table(os.path.join(BASE_DIR, "finemo", this_fold, "hits.tsv"))
    finemo_hits['motif_name'] = finemo_hits['motif_name'].map(MOTIF_DICT)
    assembled_hits.append(finemo_hits)

    # summarize hits per peak
    motif_summary = (
        finemo_hits.groupby('peak_id').agg(
            n_hits=('motif_name', 'count'),
            motif_names=('motif_name', lambda x: ' | '.join(sorted(set(x))))).reset_index()
    )

    # merge into the full peaks dataframe (preserve all rows in `peaks`)
    peaks = peaks.merge(motif_summary, on='peak_id', how='left')

    # fill missing values for peaks with no motif hits
    peaks['n_hits'] = peaks['n_hits'].fillna(0).astype(int)
    peaks['motif_names'] = peaks['motif_names'].fillna("")

    formatted_peaks.append(peaks)

## combine formatted peaks, add other annotations
all_peaks = pd.concat(formatted_peaks, axis = 0, ignore_index=True)
all_peaks.to_csv(os.path.join(OUT_DIR, "all_peaks_with_model_annot.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")

# annot = pd.read_table(categorized_elements,
#     usecols = ["cell_type", "chr", "start", "end", "element_category"])
# annot_K562 = annot[annot["cell_type"] == "K562"].drop(columns = 'cell_type').rename(columns={'start': 'peak_start', 'end': 'peak_end'})
# all_peaks = pd.merge(all_peaks, annot_K562, on=["chr", "peak_start", "peak_end"], how = "left")

# tf_annot = pd.read_table(annotated_elements).drop(columns = 'elementName').rename(columns={'start': 'peak_start', 'end': 'peak_end'})
# all_peaks = pd.merge(all_peaks, tf_annot, on=["chr", "peak_start", "peak_end"], how = "left")

# all_peaks.to_csv(os.path.join(OUT_DIR, "all_peaks_with_chromatin_annotations.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")

## combine all hits, sort, save table
all_hits = pd.concat(assembled_hits, axis = 0, ignore_index = True).sort_values(
    by = ["chr", "start"])
all_hits.to_csv(os.path.join(OUT_DIR, "hits_renamed.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")

# save as bed file, bgzip and index
all_hits_bed = all_hits[["chr", "start", "end", "motif_name", "hit_importance", "strand"]]
bed_path = os.path.join(OUT_DIR, "hits_renamed.bed")
all_hits_bed.to_csv(bed_path, sep="\t", index=False, header=False)
subprocess.run(["bgzip", "-f", bed_path], check=True)
subprocess.run(["tabix", "-p", "bed", f"{bed_path}.gz"], check=True)

## plot predicted versus observed for each fold, positive and negative strands
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
for i in range(len(FOLDS)):
    this_fold = FOLDS[i]
    # pos strand
    plot_scatter(all_peaks[all_peaks["fold"] == this_fold]["true_logcounts_0"],
            all_peaks[all_peaks["fold"] == this_fold]["pred_logcounts_0"],
            ax = axes[0, i],
            xlabel="Observed log counts",
            ylabel=f"Predicted log counts",
            title=f"EP300 counts ({this_fold}, + strand)")
    # neg strand
    plot_scatter(all_peaks[all_peaks["fold"] == this_fold]["true_logcounts_1"],
        all_peaks[all_peaks["fold"] == this_fold]["pred_logcounts_1"],
        ax = axes[1, i],
        xlabel="Observed log counts",
        ylabel=f"Predicted log counts",
        title=f"EP300 counts ({this_fold}, - strand)")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pred_vs_observed_by_fold.png"))

## plot NONZERO predicted versus observed for each fold, positive and negative strands
nonzero_peaks_0 = all_peaks[all_peaks["true_logcounts_0"] > 0]
nonzero_peaks_1 = all_peaks[all_peaks["true_logcounts_1"] > 0]

fig, axes = plt.subplots(2, 5, figsize=(20, 8))
for i in range(len(FOLDS)):
    this_fold = FOLDS[i]
    # pos strand
    plot_scatter(nonzero_peaks_0[nonzero_peaks_0["fold"] == this_fold]["true_logcounts_0"],
            nonzero_peaks_0[nonzero_peaks_0["fold"] == this_fold]["pred_logcounts_0"],
            ax = axes[0, i],
            xlabel="Observed log counts",
            ylabel=f"Predicted log counts",
            title=f"EP300 counts ({this_fold}, + strand, nonzero)")
    # neg strand
    plot_scatter(nonzero_peaks_1[nonzero_peaks_1["fold"] == this_fold]["true_logcounts_1"],
        nonzero_peaks_1[nonzero_peaks_1["fold"] == this_fold]["pred_logcounts_1"],
        ax = axes[1, i],
        xlabel="Observed log counts",
        ylabel=f"Predicted log counts",
        title=f"EP300 counts ({this_fold}, - strand, nonzero)")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pred_vs_observed_by_fold.nonzero.png"))

## plot summary across folds
fig, axes = plt.subplots(2, 2, figsize=(9, 9))

# pos strand
plot_scatter(all_peaks["true_logcounts_0"],
        all_peaks["pred_logcounts_0"],
        ax = axes[0, 0],
        xlabel="Observed log counts",
        ylabel=f"Predicted log counts",
        title=f"EP300 counts (All folds, + strand)")

# neg strand
plot_scatter(all_peaks["true_logcounts_1"],
    all_peaks["pred_logcounts_1"],
    ax = axes[1, 0],
    xlabel="Observed log counts",
    ylabel=f"Predicted log counts",
    title=f"EP300 counts (All folds, - strand)")

# pos strand, nonzero
plot_scatter(nonzero_peaks_0["true_logcounts_0"],
        nonzero_peaks_0["pred_logcounts_0"],
        ax = axes[0, 1],
        xlabel="Observed log counts",
        ylabel=f"Predicted log counts",
        title=f"EP300 counts (All folds, + strand, nonzero)")
# neg strand, nonzero
plot_scatter(nonzero_peaks_1["true_logcounts_1"],
    nonzero_peaks_1["pred_logcounts_1"],
    ax = axes[1, 1],
    xlabel="Observed log counts",
    ylabel=f"Predicted log counts",
    title=f"EP300 counts (All folds, - strand, nonzero)")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pred_vs_observed.all_folds.png"))

## plots hits per peak (all folds)
fig, axes = plt.subplots(1, 2, figsize=(8, 4))

sns.set_style("white")      # no grid
# build integer bin edges
mn = 0
mx = 15
bins = np.arange(mn, mx+2)

ax0 = axes[0]
sns.histplot(
    all_peaks['n_hits'],
    kde=False,
    bins = bins,
    color="#b778b3",
    shrink=0.8,      # controls bar width relative to bin size
    discrete = True,
    ax=ax0
)

# remove any remaining grid
ax0.grid(False)
ax0.set_xlabel('Hits per peak')
ax0.set_ylabel('Occurences')
ax0.set_title("Hits per peak (all)")

# nonzero only
ax1 = axes[1]
sns.histplot(
    nonzero_peaks_0['n_hits'],
    kde=False,
    bins = bins,
    color="#b778b3",
    shrink=0.8,      # controls bar width relative to bin size
    discrete = True,
    ax=ax1
)

# remove any remaining grid
ax1.grid(False)
ax1.set_xlabel('Hits per peak')
ax1.set_ylabel('Occurences')
ax1.set_title("Hits per peak (nonzero)")


sns.despine(trim = True)         # remove top/right spines
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "motif_hits_per_peak.all_folds.pdf"))

## plot observed or predicted counts by number of hits 
fig, axes = plt.subplots(1, 2, figsize=(8, 4))

# observed
plot_simple_boxplot(nonzero_peaks_0["n_hits"],
        nonzero_peaks_0["true_logcounts_0"],
        ax = axes[0],
        xlabel="# motif hits per peak",
        ylabel=f"Observed log counts",
        title=f"Observed EP300 by motif hits\nAll folds, + strand, nonzero")

# predicted
plot_simple_boxplot(nonzero_peaks_0["n_hits"],
        nonzero_peaks_0["pred_logcounts_0"],
        ax = axes[1],
        xlabel="# motif hits per peak",
        ylabel=f"Predicted log counts",
        title=f"Predicted EP300 counts by motif hits\nAll folds, + strand, nonzero")

sns.despine(trim=True)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "motif_hits_by_counts.all_folds.nonzero.pdf"))