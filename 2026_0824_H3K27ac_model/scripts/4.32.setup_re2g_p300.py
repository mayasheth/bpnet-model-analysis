#!/usr/bin/env python3
"""Configure ENCODE-rE2G training runs that put p300 tracks in the activity slot.

WHAT IS BEING TESTED. ABC uses the activity track through a fixed formula, and on this
project's benchmark the whole distance from the ATAC-only floor to observed p300 is 0.0993
AUPRC. rE2G instead learns a logistic regression over features, several of which are built
from the activity track (enhancer and promoter activity, combined activity, ABC score). This
asks whether a supervised model can get more out of a PREDICTED p300 track than ABC can.

FIVE DATASETS, all trained on the K562 CRISPR data with rE2G's own chromosome CV:

    ATAC_powerlaw               ATAC only                       10 features, ATAC_pl_all
    ATAC_H3K27ac_powerlaw       + observed H3K27ac (reference)  14 features, ATAC_H3K27ac_pl_all
    ATAC_p300obs_powerlaw       + observed p300 (ceiling)       14
    ATAC_p300incell_powerlaw    + predicted p300, K562 in-cell  14
    ATAC_p300multicell_powerlaw + predicted p300, K562 held out 14   <- deployment-legal

The p300 tracks go in the H3K27ac column because that is the column rE2G's feature tables
read activity from; nothing about the H3K27ac features is H3K27ac-specific except the name.
Every model uses the ALL-features table rather than the finalised subset, because the subset
was selected for observed H3K27ac and a predicted track may be useful through different
features.

ONE CANDIDATE REGION SET FOR ALL FIVE. The predicted bigwigs were painted over the July ABC
candidate regions, so each dataset's Peaks/ is seeded from that run and stamped newer than
every input in rule-chain order. Otherwise rE2G's ABC would call peaks per dataset, and the
2025 rE2G ATAC run already differs from the July set by 2 regions, which is enough to leave
a predicted region unpainted and read as zero. The ATAC-only and observed-H3K27ac datasets
are re-run on the shared set too, rather than reusing the 2025 numbers, so every comparison
is paired.

Usage: 4.32.setup_re2g_p300.py [--results-dir results/2026_0923_p300_features] [--stage-peaks]
"""
import argparse
import os
import shutil
import time

OAK = "/oak/stanford/groups/engreitz/Users/sheth"
E2G = f"{OAK}/ENCODE_rE2G"
ABC = f"{OAK}/ABC_working/ABC-Enhancer-Gene-Prediction"
JULY_PEAKS = f"{ABC}/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks"
PRED = f"{OAK}/EP300_BPNet/2026_0824_H3K27ac_model/data/abc_predicted"
DATA = f"{OAK}/Data/ENCODE/K562"

ATAC = ",".join(f"{DATA}/{x}.tn5.sorted.tagAlign.gz"
                for x in ("ENCFF077FBI", "ENCFF128WZG", "ENCFF534DCE"))
H3K27AC_OBS = ",".join(f"{DATA}/{x}.se.filtered.sorted.bam"
                       for x in ("ENCFF790GFL", "ENCFF817HMW"))
P300_OBS = ",".join(f"{DATA}/{x}.filtered.sorted.bam" for x in ("ENCFF466WKF", "ENCFF163FSR"))

DATASETS = [
    ("ATAC_powerlaw",               ""),
    ("ATAC_H3K27ac_powerlaw",       H3K27AC_OBS),
    ("ATAC_p300obs_powerlaw",       P300_OBS),
    ("ATAC_p300incell_powerlaw",    f"{PRED}/predp300_multimodal.bw"),
    ("ATAC_p300multicell_powerlaw", f"{PRED}/predp300mc_multicell_holdK562.bw"),
    # ranking-loss multi-cell models (1.36), added 2026-09-24
    ("ATAC_p300mcrank2_powerlaw",   f"{PRED}/predp300mc_multicell_holdK562_rank2.bw"),
    ("ATAC_p300mcrank10_powerlaw",  f"{PRED}/predp300mc_multicell_holdK562_rank10.bw"),
    # ranking-loss sweep, added 2026-09-25
    ("ATAC_p300mcrank20_powerlaw",  f"{PRED}/predp300mc_multicell_holdK562_rank20.bw"),
    ("ATAC_p300mcrank50_powerlaw",  f"{PRED}/predp300mc_multicell_holdK562_rank50.bw"),
    ("ATAC_p300mcrank100_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_rank100.bw"),
    # Test 1, TF-activity channel with oracle weights (1.37), added 2026-09-27
    ("ATAC_p300mcmotifact_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle.bw"),
    # Test 1 follow-ups, added 2026-09-28: deployable channel; oracle channel + ranking loss
    ("ATAC_p300mcmotifctcf_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_assoc_noctcf.bw"),
    ("ATAC_p300mcmotifrank_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_rank10_motifact_oracle.bw"),
    # teaching the model to use the oracle channel (1.38), added 2026-09-28
    ("ATAC_p300mcteachhead_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_head.bw"),
    ("ATAC_p300mcteachpair_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_pair.bw"),
    ("ATAC_p300mcteachheadpair_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_headpair.bw"),
    # fragment-size ATAC channels (1.40), added 2026-09-30
    ("ATAC_p300mcflat5p_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_flat5p.bw"),
    ("ATAC_p300mcflat5pseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_flat5p_seed1.bw"),
    ("ATAC_p300mcflat5pseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_flat5p_seed2.bw"),
    ("ATAC_p300mcfrag5p_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_frag5p.bw"),
    ("ATAC_p300mcfrag5pseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_frag5p_seed1.bw"),
    ("ATAC_p300mcfrag5pseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_frag5p_seed2.bw"),
    # RNA-based TF lists (4.54), added 2026-09-30
    ("ATAC_p300mcrnafam_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnafam.bw"),
    ("ATAC_p300mcrnafamseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnafam_seed1.bw"),
    ("ATAC_p300mcrnafamseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnafam_seed2.bw"),
    ("ATAC_p300mcrnaatac_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnaatac.bw"),
    ("ATAC_p300mcrnaatacseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnaatac_seed1.bw"),
    ("ATAC_p300mcrnaatacseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_rnaatac_seed2.bw"),
    # frozen two-stage models (1.39), added 2026-09-29
    ("ATAC_p300mctwostageoracle_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_twostage.bw"),
    ("ATAC_p300mctwostagectcf_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_assoc_noctcf_twostage.bw"),
    # seed replicates (1.37 SEED=1,2), added 2026-09-29
    ("ATAC_p300mcseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_seed1.bw"),
    ("ATAC_p300mcseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_seed2.bw"),
    ("ATAC_p300mcmotifactseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_seed1.bw"),
    ("ATAC_p300mcmotifactseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_oracle_seed2.bw"),
    ("ATAC_p300mcmotifctcfseed1_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_assoc_noctcf_seed1.bw"),
    ("ATAC_p300mcmotifctcfseed2_powerlaw", f"{PRED}/predp300mc_multicell_holdK562_motifact_assoc_noctcf_seed2.bw"),
    # p300 + H3K27ac combined (4.58) and corrected transferred H3K27ac (4.57), added 2026-09-30.
    # Observed tracks repainted over the candidates so all go through the same bigwig path.
    ("ATAC_k27obspaint_powerlaw", f"{PRED}/predk27frag_obs_k27.bw"),
    ("ATAC_p300obspaint_powerlaw", f"{PRED}/predk27frag_obs_p300.bw"),
    ("ATAC_p300k27obs_powerlaw", f"{PRED}/predk27frag_obs_p300k27.bw"),
    ("ATAC_p300k27pred_powerlaw", f"{PRED}/predk27frag_pred_p300k27.bw"),
    ("ATAC_k27txrepaint_powerlaw", f"{PRED}/predk27frag_tx_repaint.bw"),
    ("ATAC_k27txbase_powerlaw", f"{PRED}/predk27frag_tx_base.bw"),
    ("ATAC_k27txbreadthfrag_powerlaw", f"{PRED}/predk27frag_tx_breadth_frag.bw"),
    # observed H3K4me1 in the activity slot, added 2026-10-01
    ("ATAC_H3K4me1_powerlaw", ",".join(f"{DATA}/{x}.filtered.sorted.bam"
                                       for x in ("ENCFF524BOJ", "ENCFF204MWI", "ENCFF665JSC"))),
    # ATAC features plus observed H3K4me1 at E and at P as two separate features, through
    # external_features_config (4.61), added 2026-10-01
    ("ATAC_H3K4me1sep_powerlaw", ""),
]
# Datasets with features beyond the pipeline's own: name -> (feature table, external features
# config). The tables and configs are written by 4.61.
EXTERNAL = {"ATAC_H3K4me1sep_powerlaw": (
    "config/mine/feature_tables/ATAC_H3K4me1sep_pl_all.tsv",
    f"{E2G}/config/mine/external_features_K562_H3K4me1.tsv")}
FT_ATAC = "resources/feature_tables/ATAC_pl_all.tsv"
FT_ACT = "resources/feature_tables/ATAC_H3K27ac_pl_all.tsv"

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results/2026_0923_p300_features")
ap.add_argument("--stage-peaks", action="store_true",
                help="Seed every dataset's Peaks/ from the July ABC run. Run once the "
                     "results dir is about to be used; re-running re-stamps.")
ap.add_argument("--stage-new", action="store_true",
                help="Stage Peaks only for datasets that have none, and leave an existing "
                     "chrom-sizes BED alone. Re-stamping finished datasets or the BED "
                     "would make their outputs stale and re-run them.")
a = ap.parse_args()

inputs = set(ATAC.split(","))
for _, act in DATASETS:
    inputs.update(x for x in act.split(",") if x)
inputs.update((f"{E2G}/{FT_ATAC}", f"{E2G}/{FT_ACT}", JULY_PEAKS))
for ft, efc in EXTERNAL.values():
    inputs.update((f"{E2G}/{ft}", efc))
    if not os.path.exists(efc):
        continue
    with open(efc) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        inputs.update(line.rstrip("\n").split("\t")[hdr.index("source_file")] for line in fh if line.strip())
missing = [x for x in inputs if not os.path.exists(x)]
if missing:
    raise SystemExit("missing inputs:\n  " + "\n  ".join(sorted(set(missing))))

# ---------------------------------------------------------------- dataset config
cols = ["biosample", "crispr_cell_type", "DHS", "ATAC", "H3K27ac",
        "default_accessibility_feature", "HiC_file",
        "HiC_type", "HiC_resolution", "HiC_gamma", "HiC_scale", "alt_TSS", "alt_genes",
        "external_features_config", "model_dir"]
ds_path = f"{E2G}/config/mine/config_datasets_K562_p300.tsv"
with open(ds_path, "w") as f:
    f.write("\t".join(cols) + "\n")
    for name, act in DATASETS:
        r = dict.fromkeys(cols, "")
        r.update(biosample=name, crispr_cell_type="K562", ATAC=ATAC, H3K27ac=act,
                 default_accessibility_feature="ATAC",
                 external_features_config=EXTERNAL.get(name, ("", ""))[1])
        f.write("\t".join(r[c] for c in cols) + "\n")
print("wrote", ds_path)

# ---------------------------------------------------------------- model config
md_path = f"{E2G}/config/mine/config_models_p300.tsv"
with open(md_path, "w") as f:
    # crispr_dataset names a key under crispr_dataset: in the training yaml. This checkout
    # (train_multiple_celltypes) requires the column; the 2025 configs predate it.
    f.write("\t".join(["model", "dataset", "crispr_dataset", "feature_table", "polynomial",
                       "override_params"]) + "\n")
    for name, act in DATASETS:
        ft = EXTERNAL[name][0] if name in EXTERNAL else (FT_ACT if act else FT_ATAC)
        f.write("\t".join([f"{name}_all", name, "training", ft,
                           "False", ""]) + "\n")
print("wrote", md_path)

# ---------------------------------------------------------------- training yaml
# Built from the checkout's DEFAULT training yaml, not the 2025 config/mine one: the current
# Snakefile_training needs keys (max_memory_allocation_mb, the dict form of crispr_dataset)
# that the older file lacks.
base = open(f"{E2G}/config/config_training.yaml").read()
for old, new in (('model_config: "config/config_models_test.tsv"',
                  'model_config: "config/mine/config_models_p300.tsv"'),
                 ('dataset_config: "config/config_biosamples_chr22.tsv"',
                  'dataset_config: "config/mine/config_datasets_K562_p300.tsv"'),
                 ('results_dir: "results"', f'results_dir: "{a.results_dir}"')):
    assert base.count(old) == 1, f"training yaml anchor not found: {old}"
    base = base.replace(old, new)
ty_path = f"{E2G}/config/mine/config_training_p300.yaml"
open(ty_path, "w").write(base)
print("wrote", ty_path)

# ---------------------------------------------------------------- peaks
# Rule-chain order, same scheme as 2026_0824_H3K27ac_model/scripts/4.3: each stage is stamped
# later than the one feeding it, and all of them later than every input file.
ORDER = [
    lambda f: f in ("macs2_peaks.xls", "macs2_summits.bed", "macs2_peaks.narrowPeak"),
    lambda f: f == "macs2_peaks.narrowPeak.sorted",
    lambda f: f.endswith("Counts.bed"),
    lambda f: f.endswith("candidateRegions.bed"),
    lambda f: f in ("candidateRegions.qc.txt", "params.txt"),
]
if a.stage_peaks or a.stage_new:
    newest = max(os.path.getmtime(p) for _, act in DATASETS
                 for p in (ATAC + ("," + act if act else "")).split(","))
    base_t = max(newest, time.time()) + 60
    # abc_generate_chrom_sizes_bed_file feeds abc_sort_narrowpeaks, so a missing or newer
    # chrom-sizes BED re-triggers region calling for every dataset. Write it exactly as the
    # rule does (ABC/workflow/rules/macs2.smk: chrom, 0, size) and stamp it BEFORE the Peaks.
    sizes = f"{E2G}/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv"
    tmpd = f"{E2G}/{a.results_dir}/tmp"
    os.makedirs(tmpd, exist_ok=True)
    bed = f"{tmpd}/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed"
    if a.stage_new and os.path.exists(bed):
        print(f"kept {bed}")
    else:
        with open(sizes) as fi, open(bed, "w") as fo:
            for line in fi:
                f_ = line.rstrip("\n").split("\t")
                fo.write(f"{f_[0]}\t0\t{f_[1]}\n" if line.strip() else line)
        t_bed = max(base_t - 30, os.path.getmtime(sizes) + 1)
        os.utime(bed, (t_bed, t_bed))
        print(f"wrote {bed}")
    for name, _ in DATASETS:
        dst = f"{E2G}/{a.results_dir}/{name}/Peaks"
        if a.stage_new and os.path.isdir(dst):
            print(f"kept {dst}")
            continue
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(JULY_PEAKS, dst)
        files = sorted(os.listdir(dst))
        for fn in files:
            stage = next((i for i, test in enumerate(ORDER) if test(fn)), len(ORDER))
            t = base_t + 10 * stage
            os.utime(os.path.join(dst, fn), (t, t))
        os.utime(dst, (base_t + 10 * (len(ORDER) + 1),) * 2)
        print(f"staged {dst} ({len(files)} files)")
else:
    print("Peaks NOT staged (pass --stage-peaks)")
