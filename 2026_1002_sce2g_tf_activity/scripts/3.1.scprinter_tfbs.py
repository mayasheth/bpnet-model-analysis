#!/usr/bin/env python3
"""PRINT TF binding scores (scPrinter; Hu et al. 2025) at the TOBIAS motif sites, Xu K562 pseudobulk.

INPUT. scE2G's filtered fragments (the cells scE2G scores, `fragments_filtered.tsv.gz`), imported
once into a printer object; all cells form one pseudobulk group. Tn5 bias: hg38 v2
(`hg38_bias_v2.h5`); pretrained TFBS model `TFBS_model_py.pt`; default dispersion model.

SITES. Every BINDetect site (2.2) of the 2.4 motifs: GATA1, GATA1::TAL1, NFE2, MAF::NFE2,
TAL1::TCF3, CTCF, HNF4A. scPrinter assumes equal-width regions and scores only positions >=
contextRadius (100 bp) from the edge, and 62,796 candidate regions are > 1 kb, so each unique
site centre gets its own 300 bp window and is scored in 10 bp tiles; the site's score is the tile
holding the centre, passed as a manual site. Output: one row per site with TOBIAS's `K562_Xu_score` beside PRINT's, so 2.3
and 2.4 can be rerun on PRINT scores at the same sites.

Usage: 3.1.scprinter_tfbs.py [--max-sites N] [--n-jobs 16]
"""
import argparse
import glob
import os
import time

import numpy as np
import pandas as pd
import scprinter as scp

SCR = os.environ["SCRATCH"]
TB = f"{SCR}/tobias_k562_xu"
OUT = f"{SCR}/scprinter_k562_xu"
FRAG = f"{SCR}/sce2g_runs/2026_1002_k562_xu/K562_Xu_pl/fragments_filtered.tsv.gz"
MOTIFS = {"GATA1": "MA0035.5", "GATA1::TAL1": "MA0140.3", "NFE2": "MA0841.2",
          "MAF::NFE2": "MA0501.2", "TAL1::TCF3": "MA0091.2", "CTCF": "MA0139.2", "HNF4A": "MA0114.5"}
WIDTH, RADIUS, TILE = 300, 100, 10

ap = argparse.ArgumentParser()
ap.add_argument("--max-sites", type=int, default=0, help="0 = all; >0 = smoke test on the first N site centres")
ap.add_argument("--n-shards", type=int, default=1)
ap.add_argument("--shard", type=int, default=None)
ap.add_argument("--sub-shards", type=int, default=1, help="split one shard further (rerun of a failed shard)")
ap.add_argument("--sub-shard", type=int, default=None)
ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 16)))
a = ap.parse_args()
# the SLURM array index fills --sub-shard when splitting a shard, else --shard
tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
if a.sub_shards > 1:
    assert a.shard is not None, "give --shard explicitly with --sub-shards"
    a.sub_shard = tid if a.sub_shard is None else a.sub_shard
else:
    a.shard = tid if a.shard is None else a.shard
    a.sub_shard = 0
assert os.environ.get("SCPRINTER_DATA"), "set SCPRINTER_DATA so models load from the $SCRATCH cache"
os.makedirs(f"{OUT}/tmp", exist_ok=True)
tag = f"n{a.max_sites}" if a.max_sites else "all"
if a.n_shards > 1:
    assert not a.max_sites and 0 <= a.shard < a.n_shards
    tag = f"shard{a.shard}of{a.n_shards}"
    if a.sub_shards > 1:
        assert 0 <= a.sub_shard < a.sub_shards
        tag += f"_sub{a.sub_shard}of{a.sub_shards}"
genome = scp.genome.hg38
print("bias:", genome.fetch_bias(), flush=True)

# ---------------------------------------------------------------- printer (imported once)
pfile = f"{OUT}/k562_xu_printer.h5ad"
t0 = time.time()
if a.n_shards > 1:
    # binding-score runs write into the printer file, so each shard copies the raw files first
    # and only ever opens its own copy
    import shutil
    assert os.path.exists(f"{OUT}/import_complete.txt"), "import the printer first (unsharded run)"
    sd = f"{OUT}/{tag}"
    if not os.path.exists(f"{sd}/copy_complete.txt"):
        os.makedirs(sd, exist_ok=True)
        shutil.copy(pfile, sd)
        shutil.copytree(pfile.replace(".h5ad", "_supp"), f"{sd}/k562_xu_printer_supp", dirs_exist_ok=True)
        open(f"{sd}/copy_complete.txt", "w").write(time.ctime() + "\n")
    printer = scp.load_printer(f"{sd}/k562_xu_printer.h5ad", genome)
elif os.path.exists(f"{OUT}/import_complete.txt"):
    printer = scp.load_printer(pfile, genome)
else:
    printer = scp.pp.import_fragments(path_to_frags=FRAG, barcodes=None, savename=pfile, genome=genome,
                                      min_num_fragments=1, min_tsse=0, sorted_by_barcode=False,
                                      tempdir=f"{OUT}/tmp")
    open(f"{OUT}/import_complete.txt", "w").write(time.ctime() + "\n")
cells = list(printer.insertion_file.obs_names)
print(f"printer: {len(cells):,} cells ({time.time() - t0:.0f} s)", flush=True)
printer.load_disp_model()
printer.load_bindingscore_model("TF", scp.datasets.TFBS_model())

# ---------------------------------------------------------------- sites -> 300 bp windows
parts = []
for name, mid in MOTIFS.items():
    (ov,) = glob.glob(f"{TB}/bindetect/*_{mid}/*_{mid}_overview.txt")
    s = pd.read_csv(ov, sep="\t", usecols=["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand",
                                           "peak_chr", "peak_start", "peak_end", "K562_Xu_score"])
    s["motif"] = name
    parts.append(s)
S = pd.concat(parts, ignore_index=True)
S["centre"] = (S["TFBS_start"] + S["TFBS_end"]) // 2
W = (S[["TFBS_chr", "centre", "TFBS_strand"]].drop_duplicates()
     .sort_values(["TFBS_chr", "centre", "TFBS_strand"]).reset_index(drop=True))
if a.max_sites:
    W = W.iloc[:a.max_sites]
if a.n_shards > 1:
    W = np.array_split(W, a.n_shards)[a.shard].reset_index(drop=True)
    if a.sub_shards > 1:
        W = np.array_split(W, a.sub_shards)[a.sub_shard].reset_index(drop=True)
W["Start"] = W["centre"] - WIDTH // 2
W = W[W["Start"] >= 0].reset_index(drop=True)
W["region_name"] = (W["TFBS_chr"] + ":" + W["Start"].astype(str) + "-" + (W["Start"] + WIDTH).astype(str)
                    + ":" + W["TFBS_strand"])
regions = pd.DataFrame({"Chromosome": W["TFBS_chr"], "Start": W["Start"], "End": W["Start"] + WIDTH,
                        "region_name": W["region_name"]})
# one site per window: the 10 bp at the centre, on the motif's strand; Score = 1 as in tile mode, so
# the PRINT score carries footprint only, no motif-match strength. Passing sites also avoids
# pyranges.window(), which fails in this env (sorted_nearest 0.0.41: NameError 'long').
sites = [pd.DataFrame({"Chromosome": [c], "Start": [m - TILE // 2], "End": [m + TILE // 2],
                       "Strand": [st], "TF": [""], "Score": [1.0]})
         for c, m, st in zip(W["TFBS_chr"], W["centre"], W["TFBS_strand"])]
print(f"{len(S):,} sites, {len(regions):,} unique site windows", flush=True)

# ---------------------------------------------------------------- PRINT binding scores
key = f"TFBS_sites_{tag}"
t0 = time.time()
# backed=True: scPrinter writes each finished chunk to <printer>_supp/<key>.h5ad, so a timeout keeps
# partial results; reading the scores back with h5py is seconds, where per-key AnnData obsm access
# took ~6 h for 75k windows
scp.tl.get_binding_score(printer, cell_grouping=[cells], group_names=["K562_Xu"],
                         regions=regions, region_width=None, model_key="TF",
                         n_jobs=a.n_jobs, contextRadius=RADIUS, save_key=key,
                         backed=True, open_anndata=False, overwrite=True, manual_sites=sites,
                         verbose=False)
print(f"binding scores: {time.time() - t0:.0f} s for {len(regions):,} windows", flush=True)
import h5py
t0 = time.time()
with h5py.File(os.path.join(printer.file_path, f"{key}.h5ad"), "r") as h:
    obsm = h["obsm"]
    first = obsm[W["region_name"].iloc[0]][()]
    print("per-window score array shape:", first.shape, flush=True)
    assert first.reshape(-1).shape[0] == 1, f"expected 1 score per window, got {first.shape}"
    missing = [r for r in W["region_name"] if r not in obsm]
    assert not missing, f"{len(missing)} windows without a score, e.g. {missing[:3]}"
    W["print_score"] = [float(obsm[r][()].reshape(-1)[0]) for r in W["region_name"]]
print(f"read scores: {time.time() - t0:.0f} s", flush=True)
print("PRINT score summary:", W["print_score"].describe().round(3).to_dict(), flush=True)

S = S.merge(W[["TFBS_chr", "centre", "TFBS_strand", "print_score"]], on=["TFBS_chr", "centre", "TFBS_strand"], how="inner")
f = f"{OUT}/print_tfbs_sites.{tag}.tsv.gz"
S.to_csv(f, sep="\t", index=False)
print(S.groupby("motif")[["K562_Xu_score", "print_score"]].corr(method="spearman").round(3), flush=True)
print(f"wrote {f} ({len(S):,} sites)")
