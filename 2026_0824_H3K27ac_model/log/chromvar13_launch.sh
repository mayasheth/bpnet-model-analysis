#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -J chromvar13
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/log/chromvar13.%j.txt
set -euo pipefail
A=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
cd $A
PC=$SCRATCH/panel_candidates
for f in $SCRATCH/p300_panel14/*/ATAC/ENCFF*.bed.gz $SCRATCH/p300_tissues/*/ATAC/ENCFF*.bed.gz; do
    [[ -s $PC/$(basename $f) ]] || cp -v "$f" $PC/
done
for acc in ENCFF429FOV ENCFF915FZC; do
    [[ -s $PC/$acc.bed.gz ]] || curl -sfL -o $PC/$acc.bed.gz "https://www.encodeproject.org/files/$acc/@@download/$acc.bed.gz"
done
for acc in $(cut -f4 config/bulk_chromvar_peaks_panel13.tsv | tail -n +2); do
    [[ -s $PC/$acc.bed.gz ]] || { echo "missing $acc"; exit 1; }
done
echo "== default path regression (8 cells, cached) =="
$PY scripts/4.63.bulk_chromvar_panel.py --procs 8 --seeds 0 --out-prefix $SCRATCH/bulk_chromvar/regress8 > $SCRATCH/bulk_chromvar/regress8.log
$PY - <<PY
import pandas as pd, numpy as np
a = pd.read_csv("results/2026_1020_bulk_chromvar_panel_z_wide.tsv", sep="\t")
b = pd.read_csv("$SCRATCH/bulk_chromvar/regress8_z_wide.tsv", sep="\t")
assert list(a.columns) == list(b.columns) and a["motif"].equals(b["motif"])
d = np.nanmax(np.abs(a.iloc[:, 2:].to_numpy(float) - b.iloc[:, 2:].to_numpy(float)))
print(f"max |z difference| vs committed 8-cell output: {d:.3g}")
assert d < 1e-3, "default path changed"
PY
echo "== panel 13 =="
$PY scripts/4.63.bulk_chromvar_panel.py --procs 8 --panel-json config/p300_panel_13.json \
    --manifest config/bulk_chromvar_peaks_panel13.tsv \
    --cache $SCRATCH/bulk_chromvar/panel13_counts_gc_annot.npz \
    --elements data/bulk_chromvar/elements_atac_idr_union_500bp_panel13.bed.gz \
    --out-prefix results/2026_1021_bulk_chromvar_panel13
$PY scripts/4.64.make_cell_cond_table.py --z-wide results/2026_1021_bulk_chromvar_panel13_z_wide.tsv \
    --holdout K562 --panel-json config/p300_panel_13.json --out config/cell_cond_chromvar13pc.tsv
