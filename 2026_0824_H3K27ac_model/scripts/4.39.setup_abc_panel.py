#!/usr/bin/env python3
"""ABC biosample table for the four non-K562 p300-panel cell types, observed p300 as activity.

Each arm calls its OWN candidate regions from its own ATAC (these are different cell types,
so there is no shared region set to preserve, unlike every K562 comparison). The filtered
reads come from 4.38. The predicted-p300 arms for these cell types are painted over the
regions this run produces and reuse its Peaks, the same scheme 4.3 uses in K562.

Usage: 4.39.setup_abc_panel.py [--results-dir results/2026_0924_p300_panel_obs]
"""
import argparse
import glob
import os

ABC = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction"
IN = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/data/panel_abc"
COLS = ["biosample", "DHS", "ATAC", "H3K27ac", "default_accessibility_feature",
        "HiC_file", "HiC_type", "HiC_resolution", "alt_TSS", "alt_genes"]
ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results/2026_0924_p300_panel_obs")
a = ap.parse_args()

rows = []
for ct in ("GM12878", "HepG2", "A549", "MCF-7"):
    atac = sorted(glob.glob(f"{IN}/{ct}/ATAC_*.filt.bam"))
    ep = sorted(glob.glob(f"{IN}/{ct}/EP300_*.filt.bam"))
    for f in atac + ep:
        assert os.path.exists(f + ".bai"), f"not indexed (4.38 unfinished?): {f}"
    assert atac and ep, f"{ct}: missing filtered reads (ATAC {len(atac)}, EP300 {len(ep)})"
    r = dict.fromkeys(COLS, "")
    r.update(biosample=f"p300obs_{ct.replace('-', '_')}", ATAC=",".join(atac),
             H3K27ac=",".join(ep), default_accessibility_feature="ATAC")
    rows.append(r)
    print(f"{ct}: {len(atac)} ATAC, {len(ep)} EP300 files")

tsv = f"{ABC}/config/mine/config_biosamples_p300_panel_obs.tsv"
with open(tsv, "w") as f:
    f.write("\t".join(COLS) + "\n")
    for r in rows:
        f.write("\t".join(r[c] for c in COLS) + "\n")
base = open(f"{ABC}/config/config.yaml").read()
base = base.replace('biosamplesTable: "config/config_biosamples_chr22.tsv"',
                    'biosamplesTable: "config/mine/config_biosamples_p300_panel_obs.tsv"')
base = base.replace('results_dir: "results/"', f'results_dir: "{a.results_dir}"')
cfg = f"{ABC}/config/mine/config_p300_panel_obs.yaml"
open(cfg, "w").write(base)
print("wrote", tsv, "and", cfg)
