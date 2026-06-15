#!/bin/bash
#SBATCH -J finemo_v2_plots
#SBATCH -p normal
#SBATCH --time=01:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=2
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/slurm_summary_plots_%j.out
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/slurm_summary_plots_%j.err

set -euo pipefail

PROJ=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
FINEMO_V2="$PROJ/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2"
ANNOT_DIR="$FINEMO_V2/annotated_motifs"
PLOT_DIR="$FINEMO_V2/plot_annotated_motifs"
CHROMATIN_ANNOT="/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"

mkdir -p "$ANNOT_DIR" "$PLOT_DIR"

module load devel pixi/0.53.0

cd "$PROJ"

# ── 1. prepare R input files ──────────────────────────────────────────────────
# hits_renamed.tsv: v2 hits.tsv filtered to exclude FOS_JUN
# hits_per_peak.with_predictions.tsv.gz: peaks + true/predicted logcounts + n_hits

echo "Preparing R input files ..."
pixi run -e ism python - <<'PYEOF'
import os, sys
import pandas as pd
import numpy as np

PROJ = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
FINEMO_V2 = f"{PROJ}/2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2"
ANNOT_DIR = f"{FINEMO_V2}/annotated_motifs"
CHROMATIN_ANNOT = (
    "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/"
    "2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"
)
EXCLUDE = {"FOS_JUN", "GATA_TAL1_8BP"}

# hits_renamed.tsv
hits = pd.read_csv(f"{FINEMO_V2}/hits.tsv", sep="\t")
hits = hits[~hits["motif_name"].isin(EXCLUDE)]
hits.to_csv(f"{ANNOT_DIR}/hits_renamed.tsv", sep="\t", index=False)
print(f"Wrote hits_renamed.tsv  ({len(hits):,} rows)", flush=True)

# hits_per_peak.with_predictions.tsv.gz
peaks = pd.read_csv(f"{FINEMO_V2}/peaks_qc.tsv", sep="\t")[["peak_id", "chr", "peak_region_start"]]
peaks["start"] = peaks["peak_region_start"]
peaks["end"]   = peaks["peak_region_start"] + 500

annot = pd.read_csv(CHROMATIN_ANNOT, sep="\t",
                    usecols=["peak_id", "true_logcounts", "mean_pred_logcounts"])

n_hits = hits.groupby("peak_id").size().reset_index(name="n_hits")

df = (peaks
      .merge(annot, on="peak_id", how="left")
      .merge(n_hits, on="peak_id", how="left"))
df["n_hits"] = df["n_hits"].fillna(0).astype(int)

out = f"{ANNOT_DIR}/hits_per_peak.with_predictions.tsv.gz"
df.to_csv(out, sep="\t", index=False, compression="gzip")
print(f"Wrote hits_per_peak.with_predictions.tsv.gz  ({len(df):,} rows)", flush=True)
PYEOF

# ── 2. fraction explained plot ────────────────────────────────────────────────
echo "Running plot_finemo_explained_fraction.py ..."
pixi run -e ism python scripts/plot_finemo_explained_fraction.py \
    --finemo-dir "$FINEMO_V2"

# ── 3. violin + frequency/importance bar chart ────────────────────────────────
echo "Running plot_finemo_hit_importance_summary.py ..."
pixi run -e ism python scripts/plot_finemo_hit_importance_summary.py \
    --finemo-dir "$FINEMO_V2"

# ── 4. R bubble plot and related plots ────────────────────────────────────────
echo "Running 5.3.plot_finemo_hits.R ..."
conda run -n analysis Rscript scripts/5.3.plot_finemo_hits.R \
    "$ANNOT_DIR/hits_renamed.tsv" \
    "$ANNOT_DIR/hits_per_peak.with_predictions.tsv.gz" \
    "$PLOT_DIR"

echo "All done."
