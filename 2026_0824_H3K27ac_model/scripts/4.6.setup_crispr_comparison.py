#!/usr/bin/env python3
"""Build the CRISPR_comparison config for the eleven predicted-activity arms.

Eleven arms: the nine predicted-activity runs plus the two July baselines, reused rather
than recomputed. Every arm scores the SAME candidate regions (md5 7d5995ce), so differences
are attributable to the activity term alone.

ARM FAMILIES, and what each answers:
  ATAC_only            floor: accessibility alone, no H3K27ac
  ATAC_H3K27ac         ceiling: observed H3K27ac, the standard ABC activity
  pred_k562_*          geomean(observed ATAC, predicted H3K27ac), K562-trained
  pred_gm12878_*       the same with GM12878-trained models, i.e. no target-cell H3K27ac
                       anywhere in the pipeline -- the actual deployment scenario
  k27only_k562_*       predicted H3K27ac alone as activity, no geometric mean

The floor and ceiling are what make the middle interpretable: a predicted arm landing near
ATAC_H3K27ac means the model substitutes for the assay, and one landing near ATAC_only means
it adds nothing over accessibility.

Colours: greys for baselines, one hue per model family, darkening with input richness.
"""
import argparse
import os

_ap = argparse.ArgumentParser()
_ap.add_argument("--arms",
                 choices=("h3k27ac", "p300", "p300all", "accessibility", "dnaseinput"),
                 default="h3k27ac",
                 help="Which activity target's arms to benchmark. The p300 set reuses the "
                      "same July floor and observed-H3K27ac ceiling, so the two "
                      "comparisons are read on one scale.")
_a = _ap.parse_args()

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
NEW = f"{ABC}/2026_0903_predicted_activity"
JULY = f"{ABC}/2026_0721_h3k27ac_counting_comparison"
P300 = f"{ABC}/2026_0905_p300_activity"
P300TX = f"{ABC}/2026_0906_p300_transfer"
ACC = f"{ABC}/2026_0910_accessibility_activity"
CC = f"{D}/CRISPR_comparison_v3/CRISPR_comparison"
PRED_FILE = "Predictions/EnhancerPredictionsAllPutative.tsv.gz"

# (arm id, results dir, long name for plots, colour)
ARMS = [
    ("K562_ATAC_only",          JULY, "ATAC only (floor)",                 "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (ceiling)", "#404040"),
    ("pred_k562_atac",          NEW, "ATAC x predicted H3K27ac (K562 ATAC model)",       "#c6dbef"),
    ("pred_k562_sequence",      NEW, "ATAC x predicted H3K27ac (K562 sequence model)",   "#6baed6"),
    ("pred_k562_multimodal",    NEW, "ATAC x predicted H3K27ac (K562 multimodal model)", "#2171b5"),
    ("pred_gm12878_atac",       NEW, "ATAC x predicted H3K27ac (GM12878 ATAC model)",       "#fcae91"),
    ("pred_gm12878_sequence",   NEW, "ATAC x predicted H3K27ac (GM12878 sequence model)",   "#fb6a4a"),
    ("pred_gm12878_multimodal", NEW, "ATAC x predicted H3K27ac (GM12878 multimodal model)", "#cb181d"),
    ("k27only_k562_atac",       NEW, "Predicted H3K27ac alone (K562 ATAC model)",       "#c7e9c0"),
    ("k27only_k562_sequence",   NEW, "Predicted H3K27ac alone (K562 sequence model)",   "#74c476"),
    ("k27only_k562_multimodal", NEW, "Predicted H3K27ac alone (K562 multimodal model)", "#238b45"),
]

# p300 as the activity term. Same two anchors as above, so an AUPRC here is directly
# comparable to the H3K27ac run, plus observed p300 as this concept's own ceiling -- without
# it a predicted-p300 number has nothing to be read against.
P300_ARMS = [
    ("K562_ATAC_only",           JULY, "ATAC only (floor)",                 "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac",           "#404040"),
    ("p300obs_k562",             P300, "ATAC x observed p300 (ceiling)",     "#54278f"),
    ("p300pred_k562_multimodal", P300, "ATAC x predicted p300 (multimodal)", "#807dba"),
    ("p300pred_k562_atac",       P300, "ATAC x predicted p300 (ATAC model)", "#bcbddc"),
    ("p300only_k562_multimodal", P300, "Predicted p300 alone (multimodal)",  "#d94801"),
    ("p300only_k562_atac",       P300, "Predicted p300 alone (ATAC model)",  "#fdae6b"),
]

# Everything p300 in one comparison, so local and transferred arms are scored on the
# IDENTICAL pair set and the paired bootstrap can compare them directly. Splitting them
# across two runs would leave transfer-vs-local uncomparable except through the anchors.
P300_ALL_ARMS = P300_ARMS + [
    ("p300pred_gm12878_multimodal", P300TX,
     "ATAC x predicted p300 (GM12878 model -> K562)", "#08519c"),
    ("p300only_gm12878_multimodal", P300TX,
     "Predicted p300 alone (GM12878 model -> K562)", "#6baed6"),
]

# The accessibility assay as the whole activity term. Asks what the ABC activity slot is
# worth before any H3K27ac model is involved: real ATAC against real DNase, and then whether
# DNase predicted from ATAC recovers the difference.
#
# NO JULY ANCHORS HERE, DELIBERATELY. The July ATAC-only (0.457) and DHS-only (0.576) arms
# count reads through count_bam/count_tagalign while these four count bigwig values, so
# mixing them into one comparison would put a counting-path difference inside the very
# contrast the run exists to measure. The two real-track arms below are this run's own
# floor and ceiling; the July numbers stay external reference points.
ACC_ARMS = [
    ("acc_real_atac",    ACC, "Real ATAC (floor)",                          "#bdbdbd"),
    ("acc_real_dnase",   ACC, "Real DNase (ceiling for this concept)",      "#404040"),
    ("acc_conv_k562",    ACC, "DNase predicted from ATAC (K562 converter)", "#762A83"),
    ("acc_conv_gm12878", ACC, "DNase predicted from ATAC (GM12878 converter -> K562)",
     "#9970ab"),
]

# The DNase-INPUT model as the activity term. F-010 gave the model +0.037 in-cell K562 by
# swapping DNase for ATAC as its accessibility input, and it has never been benchmarked.
# F-012 separately showed real DNase beats real ATAC in the activity slot by +0.0709, the
# largest downstream effect in the project, so the question is whether a MODEL that reads
# DNase adds anything on top of that assay swap.
#
# WHICH CONTRASTS IN THIS TABLE ARE LEGITIMATE, because the arms do not all share a counting
# path and a careless pairing would measure count_bigwig vs count_tagalign instead of the
# thing of interest:
#   pred_dnase_atacacc  vs  pred_k562_multimodal   VALID. Both geomean(ATAC tagAligns,
#                           predicted bigwig); only the predictor's input differs. THE
#                           headline contrast.
#   k27only_dnase       vs  k27only_k562_multimodal  VALID. Both a single predicted bigwig.
#   pred_dnase_dhsacc   vs  acc_real_dnase         VALID. Both count_bigwig throughout.
#                           Answers "does predicted H3K27ac add anything on top of real
#                           DNase in the activity slot".
#   anything            vs  the two July anchors   READ AS REFERENCE ONLY. Those count reads
#                           through count_bam/count_tagalign. They are included because the
#                           floor and ceiling are what make any middle number readable, and
#                           because the existing h3k27ac family already reports against
#                           them, but a small difference against them is not attributable.
DNASE_IN = f"{ABC}/2026_0915_dnase_input_activity"
DNASE_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor, reference)",          "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (ceiling, reference)",
     "#404040"),
    ("acc_real_dnase",            ACC,  "Real DNase alone (F-012 winner)",       "#0F6E56"),
    ("pred_k562_multimodal",      NEW,  "ATAC x predicted H3K27ac (ATAC-input model)",
     "#2171b5"),
    ("pred_dnase_atacacc",        DNASE_IN,
     "ATAC x predicted H3K27ac (DNase-input model)",                             "#08306b"),
    ("k27only_k562_multimodal",   NEW,  "Predicted H3K27ac alone (ATAC-input model)",
     "#74c476"),
    ("k27only_dnase",             DNASE_IN,
     "Predicted H3K27ac alone (DNase-input model)",                              "#238b45"),
    ("pred_dnase_dhsacc",         DNASE_IN,
     "Real DNase x predicted H3K27ac (DNase-input model)",                       "#762A83"),
]

if _a.arms == "dnaseinput":
    ARMS = DNASE_ARMS
elif _a.arms == "accessibility":
    ARMS = ACC_ARMS
elif _a.arms == "p300":
    ARMS = P300_ARMS
elif _a.arms == "p300all":
    ARMS = P300_ALL_ARMS
TAG = {"h3k27ac": "predicted_activity", "p300": "p300_activity",
       "p300all": "p300_all", "accessibility": "accessibility_activity",
       "dnaseinput": "dnase_input_activity"}[_a.arms]
RUN = {"h3k27ac": "2026_0904_predicted_activity", "p300": "2026_0905_p300_activity",
       "p300all": "2026_0906_p300_all",
       "accessibility": "2026_0910_accessibility_activity",
       "dnaseinput": "2026_0915_dnase_input_activity"}[_a.arms]

BASELINES = [
    ("distToTSS",      "FALSE", "mean", "Inf", "TRUE",  "Distance to TSS",      "#c5cad7"),
    ("nearestTSS",     "TRUE",  "max",  "0",   "FALSE", "Nearest TSS",          "#6e788d"),
    ("nearestGene",    "TRUE",  "max",  "0",   "FALSE", "Nearest gene",         "#435369"),
    ("within100kbTSS", "TRUE",  "max",  "0",   "FALSE", "Within 100kb of TSS",  "#1c2a43"),
]

missing = []
for arm, root, _n, _c in ARMS:
    p = f"{root}/{arm}/{PRED_FILE}"
    if not os.path.exists(p):
        missing.append(p)
if missing:
    raise SystemExit("missing prediction files:\n  " + "\n  ".join(missing))

# --- pred_config -------------------------------------------------------------
os.makedirs(f"{CC}/resources/pred_config", exist_ok=True)
pc = f"{CC}/resources/pred_config/pred_config_{TAG}.tsv"
with open(pc, "w") as f:
    f.write("pred_id\tpred_col\tboolean\talpha\taggregate_function\tfill_value\t"
            "inverse_predictor\tpred_name_long\tcolor\n")
    for col, boolean, agg, fill, inv, name, colour in BASELINES:
        f.write(f"baseline\t{col}\t{boolean}\t\t{agg}\t{fill}\t{inv}\t{name}\t{colour}\n")
    for arm, _root, name, colour in ARMS:
        f.write(f"{arm}\tABC.Score\tFALSE\t\tmax\t0\tFALSE\t{name}\t{colour}\n")
print("wrote", pc)

# --- cell type mappings ------------------------------------------------------
os.makedirs(f"{CC}/resources/cell_type_mapping", exist_ok=True)
for arm, _root, _n, _c in ARMS:
    m = f"{CC}/resources/cell_type_mapping/{arm}.txt"
    if not os.path.exists(m):
        with open(m, "w") as f:
            f.write(f"experiment\tpredictions\nK562\t{arm}\n")
print(f"cell type mappings present for all {len(ARMS)} arms")

# --- comparison config -------------------------------------------------------
GENE_U = f"{D}/scE2G_temp/scE2G/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed6"
TSS_U = f"{D}/scE2G_temp/scE2G/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6"
EXPT = f"{CC}/resources/crispr_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
for p in (GENE_U, TSS_U, EXPT):
    if not os.path.exists(p):
        raise SystemExit(f"missing reference: {p}")

cfg = f"{CC}/config/config_{TAG}.yml"
with open(cfg, "w") as f:
    f.write("# ABC activity term benchmarked against CRISPR data.\n")
    f.write(f"# All {len(ARMS)} arms score the identical candidate-region set (md5 7d5995ce).\n\n")
    f.write(f"comparisons:\n  {RUN}:\n    pred:\n")
    for arm, root, _n, _c in ARMS:
        f.write(f"      {arm}: {root}/{arm}/{PRED_FILE}\n")
    f.write(f"    expt: {EXPT}\n")
    f.write(f'    gene_universe: "{GENE_U}"\n')
    f.write(f'    tss_universe:  "{TSS_U}"\n')
    f.write(f'    pred_config: "{pc}"\n')
    f.write("    cell_type_mapping:\n")
    for arm, _root, _n, _c in ARMS:
        f.write(f"      {arm}: {CC}/resources/cell_type_mapping/{arm}.txt\n")
    f.write("    dist_bins_kb: [0, 10, 100, 2500]\n")
    f.write("    include_col: Null\n    gene_features: Null\n")
    f.write("    enh_features: Null\n    enh_assays: Null\n")
print("wrote", cfg)
print(f"\n{len(ARMS)} arms wired for {_a.arms}: floor=K562_ATAC_only, "
      f"{len(ARMS) - 2} non-anchor arms, run name {RUN}")
