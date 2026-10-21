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
import subprocess

_ap = argparse.ArgumentParser()
_ap.add_argument("--arms",
                 choices=("h3k27ac", "p300", "p300all", "accessibility", "dnaseinput",
                          "multitask", "slopefix", "p300fix", "multicell", "oracle",
                          "rank", "nbr", "nbrx", "ranksweep", "gata", "motifact", "motifact2", "teach", "seeds", "rna", "frag", "k27frag", "noatac", "tissue", "panelmean", "film13"),
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

# The multi-task model (F-015): DNase profile head, H3K27ac counts head, ATAC ALONE at
# inference. The only remaining candidate that meets the real deployment constraint, since
# F-018's winner needs DNase at prediction time and the ATAC-to-DNase converter is closed
# (F-014).
#
# THE DECISIVE CONTRAST IS pred_multitask_atacacc vs pred_ep100_atacacc, NOT vs
# pred_k562_multimodal. Changing the loss changes the stopping rule too: the multi-task arm
# trained 53-99 epochs against the early-stopped baseline's 32-55, so a comparison against
# the latter cannot separate "the auxiliary task helped" from "it trained longer". F-015 had
# to exclude exactly this confound upstream, where the epoch-matched baseline moved +0.0000
# (p=0.92); the same control belongs here. pred_k562_multimodal is kept as the published
# reference point.
#
# Secondary: pred_multitask_atacacc vs pred_k562_multimodal. Same architecture,
# same ATAC input, same counts target; the only difference is that one trained its profile
# head against DNase instead of against H3K27ac. Both are geomean(ATAC tagAligns, predicted
# bigwig), so they are matched on counting path.
#
# pred_dnase_atacacc is carried as the REFERENCE CEILING FOR THIS QUESTION: it is what the
# same architecture achieves when it gets to READ DNase rather than merely be supervised by
# it (+0.0403 over the floor, F-018). If the multi-task arm approaches it, DNase supervision
# recovers some of what DNase input provides without needing the assay. If it sits at the
# ATAC-input arm, the auxiliary task buys nothing downstream either.
MULTITASK = f"{ABC}/2026_0916_multitask_activity"
MULTITASK_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)",                      "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (reference)",     "#404040"),
    ("pred_k562_multimodal",      NEW,  "ATAC x predicted H3K27ac (plain ATAC-input model)",
     "#c6dbef"),
    ("pred_multitask_atacacc",    MULTITASK,
     "ATAC x predicted H3K27ac (multi-task, DNase-supervised profile head)",       "#2171b5"),
    ("k27only_k562_multimodal",   NEW,  "Predicted H3K27ac alone (plain ATAC-input model)",
     "#c7e9c0"),
    ("k27only_multitask",         MULTITASK,
     "Predicted H3K27ac alone (multi-task)",                                       "#238b45"),
    ("pred_ep100_atacacc",        MULTITASK,
     "ATAC x predicted H3K27ac (epoch-matched baseline, no DNase supervision)",    "#9ecae1"),
    ("k27only_ep100",             MULTITASK,
     "Predicted H3K27ac alone (epoch-matched baseline)",                           "#a1d99b"),
    ("pred_dnase_atacacc",        DNASE_IN,
     "ATAC x predicted H3K27ac (DNase-INPUT model, needs DNase at inference)",     "#762A83"),
]

if _a.arms == "multitask":
    ARMS = MULTITASK_ARMS
elif _a.arms == "dnaseinput":
    ARMS = DNASE_ARMS
# SLOPEFIX: everything needed to test F-021's two cheap strategies, in ONE comparison so the
# paired bootstrap can compare them against the arms they are meant to improve.
#
# WHY ONE RUN. F-021's numbers came from two separate comparisons (p300 in 2026_0906_p300_all,
# H3K27ac in 2026_0904_predicted_activity), which means they were only comparable through the
# shared July anchors and could not be paired against each other. The two strategies are each
# worth about 0.01-0.05, the same size as the between-comparison slop, so they have to be
# scored on one pair set or the result is unreadable. Listing an already-finished arm's
# directory here re-scores it; it does not re-run its ABC.
#
# TWO STRATEGIES, FOUR NEW ARMS:
#   *_accnorm      the same p300 models re-predicted with the accessibility input standardized
#                  on the PREDICTION regions instead of each model's training windows. Removes
#                  a 0.73 sd (transferred) against 0.38 sd (in-cell) centring asymmetry that
#                  had nothing to do with transfer. BOTH arms move, so the pair stays honest.
#   *_residual     the residual objective, where the model predicts H3K27ac counts minus an
#                  accessibility-only model's counts, so the accessibility slope is not a
#                  channel it can get wrong. F-005 judged this objective in-cell only.
#
# COUNTING PATHS. Every non-anchor arm here is geomean(real ATAC tagAligns, predicted bigwig)
# or a single predicted bigwig, so they are matched. The two July anchors count reads through
# count_bam/count_tagalign and are reference only, as everywhere else in this script.
ACCNORM_D = f"{ABC}/2026_0917_accnorm_activity"
RESID_D = f"{ABC}/2026_0917_residual_transfer"
SLOPEFIX_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor, reference)",        "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (ceiling, reference)",
     "#404040"),
    # p300, geomean form: the F-021 headline pair and its renormalised twin
    ("p300pred_k562_multimodal",    P300,   "ATAC x predicted p300 (K562)",            "#807dba"),
    ("p300pred_gm12878_multimodal", P300TX, "ATAC x predicted p300 (GM12878 -> K562)", "#08519c"),
    ("p300pred_k562_multimodal_accnorm", ACCNORM_D,
     "ATAC x predicted p300 (K562, region-normalised)", "#bcbddc"),
    ("p300pred_gm12878_multimodal_accnorm", ACCNORM_D,
     "ATAC x predicted p300 (GM12878 -> K562, region-normalised)", "#6baed6"),
    # p300, prediction-alone form. F-019 found this form more sensitive than the geomean,
    # which multiplies the prediction by real ATAC and compresses differences between models.
    ("p300only_k562_multimodal",    P300,   "Predicted p300 alone (K562)",            "#d94801"),
    ("p300only_gm12878_multimodal", P300TX, "Predicted p300 alone (GM12878 -> K562)", "#fdae6b"),
    ("p300only_k562_multimodal_accnorm", ACCNORM_D,
     "Predicted p300 alone (K562, region-normalised)", "#fdd0a2"),
    ("p300only_gm12878_multimodal_accnorm", ACCNORM_D,
     "Predicted p300 alone (GM12878 -> K562, region-normalised)", "#feedde"),
    # H3K27ac, plain against residual objective, in-cell and transferred
    ("pred_k562_multimodal",    NEW, "ATAC x predicted H3K27ac (K562)",            "#2171b5"),
    ("pred_gm12878_multimodal", NEW, "ATAC x predicted H3K27ac (GM12878 -> K562)", "#cb181d"),
    ("pred_k562_residual",    RESID_D,
     "ATAC x predicted H3K27ac (K562, residual objective)", "#238b45"),
    ("pred_gm12878_residual", RESID_D,
     "ATAC x predicted H3K27ac (GM12878 -> K562, residual objective)", "#74c476"),
]

# F-025 re-test: the corrected GM12878 p300 model against the corrupted one it replaces,
# the in-cell model, and the floor, all on one pair set so every delta is paired.
P300FIX_D = f"{ABC}/2026_0919_p300_corrected"
P300FIX_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor, reference)", "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (ceiling, reference)",
     "#404040"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "ATAC x predicted p300 (K562 in-cell)", "#807dba"),
    ("p300pred_gm12878_multimodal", P300TX,
     "ATAC x predicted p300 (GM12878 -> K562, CORRUPTED target)", "#cb181d"),
    ("p300pred_gm12878_corrected", P300FIX_D,
     "ATAC x predicted p300 (GM12878 -> K562, corrected target)", "#238b45"),
    ("p300only_gm12878_multimodal", P300TX,
     "Predicted p300 alone (GM12878 -> K562, CORRUPTED target)", "#fdae6b"),
    ("p300only_gm12878_corrected", P300FIX_D,
     "Predicted p300 alone (GM12878 -> K562, corrected target)", "#74c476"),
]

# F-026 downstream: the multi-cell-type p300 model on the K562 CRISPR benchmark.
#
# The model never saw K562. Upstream it clears K562's assay-matched accessibility-only floor
# by +0.0614 top-quintile Pearson, the project's first positive transfer result; the CRISPR
# benchmark has disagreed with Pearson five times, so the claim is not settled until it runs
# here.
#
# THE COMPARATOR IS POLICY-MATCHED. A multi-cell model has no single training statistic for
# accessibility, so it must be predicted with target-derived statistics. Scoring it against
# p300pred_gm12878_corrected, which used its own training statistics, would confound the
# training set with the normalisation policy. The _tgtnorm arm is the same corrected
# single-source model re-predicted under the multi-cell model's statistics, and IT is the
# comparator the headline rests on. The own-statistics arm is carried as a reference so the
# normalisation effect itself stays visible.
MULTICELL_D = f"{ABC}/2026_0922_p300_multicell"
MULTICELL_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor, reference)", "#bdbdbd"),
    ("K562_ATAC_H3K27ac_element", JULY, "ATAC x observed H3K27ac (ceiling, reference)",
     "#404040"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "ATAC x predicted p300 (K562 in-cell)", "#807dba"),
    ("p300pred_gm12878_corrected", f"{ABC}/2026_0919_p300_corrected",
     "ATAC x predicted p300 (GM12878 -> K562, own-stats norm)", "#cb181d"),
    ("p300pred_gm12878_corrected_tgtnorm", MULTICELL_D,
     "ATAC x predicted p300 (GM12878 -> K562, target-derived norm)", "#fd8d3c"),
    ("p300pred_multicell_holdK562", MULTICELL_D,
     "ATAC x predicted p300 (4 cell types, K562 held out)", "#238b45"),
    ("p300only_gm12878_corrected_tgtnorm", MULTICELL_D,
     "Predicted p300 alone (GM12878 -> K562, target-derived norm)", "#fdd0a2"),
    ("p300only_multicell_holdK562", MULTICELL_D,
     "Predicted p300 alone (4 cell types, K562 held out)", "#74c476"),
]

# Activity-loss decomposition (4.31). For each model, three arms on one pair set:
#   A  predicted everywhere          (existing)
#   B  observed ranks at tested elements, predicted elsewhere   (oracle)
#   C  observed everywhere           (existing ceiling)
# B - A is the AUPRC owed to ranking the tested elements; C - B is the AUPRC owed to every
# other candidate element.
ORACLE_D = f"{ABC}/2026_0923_p300_oracle"
ORACLE_ARMS = [
    ("K562_ATAC_only",              JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",                P300, "C: observed p300 everywhere", "#54278f"),
    ("p300pred_multicell_holdK562", f"{ABC}/2026_0922_p300_multicell",
     "A: multi-cell predicted everywhere", "#238b45"),
    ("p300pred_oracle_multicell_holdK562", ORACLE_D,
     "B: multi-cell, observed ranks at tested elements", "#74c476"),
    ("p300pred_k562_multimodal",    P300, "A': in-cell predicted everywhere", "#807dba"),
    ("p300pred_oracle_k562_multimodal", ORACLE_D,
     "B': in-cell, observed ranks at tested elements", "#bcbddc"),
]

# Ranking-loss models (1.36) against the MSE-only multi-cell model they were derived from.
# Same panel, folds, normalisation and stopping rule; stopping epochs 36.4 (w=2), 34.0 (w=10)
# against 34.6, so training length is not a confound. Headline contrasts:
#   p300pred_multicell_holdK562_rank{2,10}  vs  p300pred_multicell_holdK562
#   p300only_multicell_holdK562_rank{2,10}  vs  p300only_multicell_holdK562
RANK_D = f"{ABC}/2026_0924_p300_rank"
RANK_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "ATAC x predicted p300 (K562 in-cell)", "#807dba"),
    ("p300pred_multicell_holdK562", MULTICELL_D,
     "ATAC x predicted p300 (multi-cell, MSE only)", "#238b45"),
    ("p300pred_multicell_holdK562_rank2", RANK_D,
     "ATAC x predicted p300 (multi-cell, ranking loss w=2)", "#6baed6"),
    ("p300pred_multicell_holdK562_rank10", RANK_D,
     "ATAC x predicted p300 (multi-cell, ranking loss w=10)", "#08519c"),
    ("p300only_multicell_holdK562", MULTICELL_D,
     "Predicted p300 alone (multi-cell, MSE only)", "#74c476"),
    ("p300only_multicell_holdK562_rank2", RANK_D,
     "Predicted p300 alone (multi-cell, ranking loss w=2)", "#fdae6b"),
    ("p300only_multicell_holdK562_rank10", RANK_D,
     "Predicted p300 alone (multi-cell, ranking loss w=10)", "#d94801"),
]

# Neighbour upper bound (4.36/4.37). Headline contrast: boostnbr vs boostnonbr for each
# model, both refit to K562 labels the same way, so the difference is the neighbour terms.
NBR_D = f"{ABC}/2026_0924_p300_nbr"
NBR_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "ATAC x predicted p300 (K562 in-cell)", "#807dba"),
    ("p300pred_nbr_incell_boostnonbr", NBR_D,
     "in-cell, boosting refit without neighbours", "#bcbddc"),
    ("p300pred_nbr_incell_boostnbr", NBR_D,
     "in-cell, boosting refit with neighbours", "#3f007d"),
    ("p300pred_multicell_holdK562", MULTICELL_D,
     "ATAC x predicted p300 (multi-cell, MSE only)", "#238b45"),
    ("p300pred_nbr_multicell_boostnonbr", NBR_D,
     "multi-cell, boosting refit without neighbours", "#a1d99b"),
    ("p300pred_nbr_multicell_boostnbr", NBR_D,
     "multi-cell, boosting refit with neighbours", "#00441b"),
    ("p300pred_multicell_holdK562_rank10", RANK_D,
     "multi-cell, ranking loss w=10 (reference)", "#08519c"),
]

# Cross-class neighbour upper bound, scored together with the same-class arms and their
# no-neighbour comparators so every contrast is paired. The atac arm's comparator is the floor.
NBRX_D = f"{ABC}/2026_0924_p300_nbrx"
NBRX_ARMS = NBR_ARMS + [
    ("p300pred_nbrx_multicell_boostnbr", NBRX_D,
     "multi-cell, boosting with cross-class neighbours", "#006d2c"),
    ("p300pred_nbrx_incell_boostnbr", NBRX_D,
     "in-cell, boosting with cross-class neighbours", "#54278f"),
    ("p300pred_nbrx_atac_boostnbr", NBRX_D,
     "no p300 model: ATAC + neighbour ATAC, boosting", "#737373"),
]

# Ranking-loss dose-response: MSE-only, then weights 2, 10, 20, 50, 100, one pair set.
RANKSWEEP_D = f"{ABC}/2026_0925_p300_ranksweep"
RANKSWEEP_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_multicell_holdK562", MULTICELL_D, "multi-cell, MSE only", "#238b45"),
    ("p300pred_multicell_holdK562_rank2", RANK_D, "ranking loss w=2", "#c6dbef"),
    ("p300pred_multicell_holdK562_rank10", RANK_D, "ranking loss w=10", "#6baed6"),
    ("p300pred_multicell_holdK562_rank20", RANKSWEEP_D, "ranking loss w=20", "#2171b5"),
    ("p300pred_multicell_holdK562_rank50", RANKSWEEP_D, "ranking loss w=50", "#08519c"),
    ("p300pred_multicell_holdK562_rank100", RANKSWEEP_D, "ranking loss w=100", "#08306b"),
]

# GATA oracle. Headline: gata_incell - multicell (what K562 GATA biology is worth), and
# gata_incell - ctrl_incell (the part specific to GATA rather than to replacing N regions).
GATA_D = f"{ABC}/2026_0926_p300_gataoracle"
GATA_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "in-cell model", "#807dba"),
    ("p300pred_multicell_holdK562", MULTICELL_D, "multi-cell model", "#238b45"),
    ("p300pred_gataoracle_gata_incell", GATA_D, "multi-cell, GATA elements at in-cell ranks", "#d94801"),
    ("p300pred_gataoracle_gata_obs", GATA_D, "multi-cell, GATA elements at observed ranks", "#7f2704"),
    ("p300pred_gataoracle_ctrl_incell", GATA_D, "multi-cell, matched non-GATA at in-cell ranks", "#fdae6b"),
    ("p300pred_gataoracle_ctrl_obs", GATA_D, "multi-cell, matched non-GATA at observed ranks", "#fdd0a2"),
]

# Test 1. Headline: motifact_oracle vs the MSE-only multi-cell model (same training except
# the channel). The GATA oracle arms and rank w=10 are carried as reference points.
MOTIFACT_D = f"{ABC}/2026_0927_p300_motifact"
MOTIFACT_ARMS = [
    ("K562_ATAC_only",            JULY, "ATAC only (floor)", "#bdbdbd"),
    ("p300obs_k562",              P300, "ATAC x observed p300 (ceiling)", "#54278f"),
    ("p300pred_k562_multimodal",  P300, "in-cell model", "#807dba"),
    ("p300pred_multicell_holdK562", MULTICELL_D, "multi-cell, MSE only", "#238b45"),
    ("p300pred_multicell_holdK562_rank10", RANK_D, "multi-cell, ranking loss w=10", "#08519c"),
    ("p300pred_multicell_holdK562_motifact_oracle", MOTIFACT_D,
     "multi-cell + TF-activity channel (oracle weights)", "#d94801"),
    ("p300pred_gataoracle_gata_obs", GATA_D, "GATA elements at observed ranks (reference)", "#7f2704"),
]

# Test 1 follow-ups (2026-09-28): the deployable channel and oracle channel + ranking loss.
MOTIFACT2_D = f"{ABC}/2026_0928_p300_motifact2"
MOTIFACT2_ARMS = [a_ for a_ in MOTIFACT_ARMS if not a_[0].startswith("p300pred_gataoracle")] + [
    ("p300pred_multicell_holdK562_motifact_assoc_noctcf", MOTIFACT2_D,
     "multi-cell + TF-activity channel (ATAC association, CTCF removed)", "#fd8d3c"),
    ("p300pred_multicell_holdK562_rank10_motifact_oracle", MOTIFACT2_D,
     "multi-cell + oracle channel + ranking loss w=10", "#7f2704"),
]

# Teaching the model to use the oracle channel (1.38): comparator is the oracle channel model.
TEACH_D = f"{ABC}/2026_0928_p300_teach"
TEACH_ARMS = MOTIFACT2_ARMS + [
    ("p300pred_multicell_holdK562_motifact_oracle_head", TEACH_D,
     "oracle channel, additive motif head", "#6baed6"),
    ("p300pred_multicell_holdK562_motifact_oracle_pair", TEACH_D,
     "oracle channel, cross-cell paired loss", "#3182bd"),
    ("p300pred_multicell_holdK562_motifact_oracle_headpair", TEACH_D,
     "oracle channel, head + paired loss", "#08306b"),
]

# Seed replicates (1.37 SEED=1,2): the original run of each model plus two replicates.
SEEDS_D = f"{ABC}/2026_0929_p300_seeds"
SEEDS_ARMS = [a_ for a_ in MOTIFACT2_ARMS if "rank10_motifact" not in a_[0]] + [
    (f"p300pred_multicell_holdK562{m}_seed{n}", SEEDS_D, f"{lab}, seed {n}", col)
    for m, lab, col in (("", "multi-cell, MSE only", "#238b45"),
                        ("_motifact_oracle", "oracle channel", "#d94801"),
                        ("_motifact_assoc_noctcf", "assoc_noctcf channel", "#fd8d3c"))
    for n in (1, 2)] + [
    ("p300pred_multicell_holdK562_motifact_oracle_twostage", SEEDS_D,
     "two-stage: frozen DNA+ATAC, oracle channel correction", "#6a51a3"),
    ("p300pred_multicell_holdK562_motifact_assoc_noctcf_twostage", SEEDS_D,
     "two-stage: frozen DNA+ATAC, assoc_noctcf channel correction", "#9e9ac8")]

# RNA-based TF lists (4.54), three runs each; the seed arms ride along for seed-averaged pairs.
RNA_D = f"{ABC}/2026_0930_p300_rna"
RNA_ARMS = [x for x in SEEDS_ARMS if "twostage" not in x[0]] + [
    (f"p300pred_multicell_holdK562_motifact_{w}{s_}", RNA_D, f"{lab}{' ' + s_[1:] if s_ else ''}", col)
    for w, lab, col in (("rnafam", "TF list from RNA", "#3182bd"),
                        ("rnaatac", "TF list from RNA + ATAC", "#6baed6"))
    for s_ in ("", "_seed1", "_seed2")]

# Fragment-size ATAC channels (1.40), 3 runs each; the MSE-only seed arms ride along.
FRAG_D = f"{ABC}/2026_0930_p300_frag"
FRAG_ARMS = [x for x in SEEDS_ARMS if "motifact" not in x[0]] + [
    (f"p300pred_multicell_holdK562_{w}{s_}", FRAG_D, f"{lab}{' ' + s_[1:] if s_ else ''}", col)
    for w, lab, col in (("flat5p", "flat 5' ATAC input", "#74c476"),
                        ("frag5p", "5' ATAC + 4 fragment-size channels", "#006d2c"))
    for s_ in ("", "_seed1", "_seed2")]

# Transferred H3K27ac corrected with within-cell-ranked fragment features (4.57, F-039).
K27FRAG_D = f"{ABC}/2026_0930_k27frag"
K27FRAG_ARMS = [x for x in ARMS if x[0] in ("K562_ATAC_only", "K562_ATAC_H3K27ac_element",
                                          "pred_k562_multimodal", "pred_gm12878_multimodal")] + [
    ("k27frag_tx_repaint", K27FRAG_D, "GM12878 model -> K562, repainted (control)", "#fcbba1"),
    ("k27frag_tx_base", K27FRAG_D, "+ recalibrated in GM12878, no fragments (control)", "#fb6a4a"),
    ("k27frag_tx_frag", K27FRAG_D, "+ ranked fragment features (linear)", "#a50f15"),
    ("k27frag_tx_frag_gbm", K27FRAG_D, "+ ranked fragment features (boosted trees)", "#67000d"),
    ("k27frag_atac_frag", K27FRAG_D, "ATAC + fragment features, no sequence model", "#737373"),
    ("k27frag_tx_breadth", K27FRAG_D, "+ ATAC within 1 kb, no fragment lengths (control)", "#fd8d3c"),
    ("k27frag_tx_breadth_frag", K27FRAG_D, "+ ATAC within 1 kb + fragment features", "#7f2704"),
    # 4.58: p300 and H3K27ac combined, observed (ceilings) and predicted (deployable)
    ("k27frag_obs_k27", K27FRAG_D, "ATAC x observed H3K27ac (painted)", "#252525"),
    ("k27frag_obs_p300", K27FRAG_D, "ATAC x observed p300 (painted)", "#3f007d"),
    ("k27frag_obs_p300k27", K27FRAG_D, "ATAC x observed p300 x observed H3K27ac", "#000000"),
    ("k27frag_obsonly_p300k27", K27FRAG_D, "observed p300 x observed H3K27ac, no ATAC", "#636363"),
    ("k27frag_obsonly_p300", K27FRAG_D, "observed p300 alone, no ATAC", "#6a51a3"),
    ("k27frag_obsonly_k27", K27FRAG_D, "observed H3K27ac alone, no ATAC", "#969696"),
    ("k27frag_pred_p300k27", K27FRAG_D, "ATAC x predicted p300 x corrected predicted H3K27ac", "#54278f"),
    ("k27frag_pred_p300ch_k27", K27FRAG_D, "same, p300 with TF-activity channel", "#807dba"),
    ("k27frag_pred_p300_k27tx", K27FRAG_D, "ATAC x predicted p300 x uncorrected H3K27ac", "#bcbddc")]
# predicted p300 alone, so the predicted combinations can be paired against it
K27FRAG_ARMS += [x for x in SEEDS_ARMS if x[0] in ("p300pred_multicell_holdK562",
                 "p300pred_multicell_holdK562_motifact_assoc_noctcf")]
assert len(K27FRAG_ARMS) == 22, len(K27FRAG_ARMS)

# Activity without the ATAC term (F-043 follow-up) and observed H3K4me1; the geomean versions
# of the same tracks ride along from their own runs for pairing.
NOATAC_D = f"{ABC}/2026_1001_noatac"
_seen, NOATAC_ARMS = set(), []
for x in ARMS + P300_ARMS + SEEDS_ARMS:
    if x[0] in ("K562_ATAC_only", "K562_ATAC_H3K27ac_element", "p300obs_k562",
                "p300pred_multicell_holdK562", "p300pred_multicell_holdK562_seed1",
                "p300pred_multicell_holdK562_seed2", "p300pred_multicell_holdK562_motifact_assoc_noctcf",
                "p300pred_multicell_holdK562_motifact_assoc_noctcf_seed1",
                "p300pred_multicell_holdK562_motifact_assoc_noctcf_seed2") and x[0] not in _seen:
        _seen.add(x[0])
        NOATAC_ARMS.append(x)
NOATAC_ARMS += [(f"p300only_multicell_holdK562{m}{s_}", NOATAC_D, f"predicted p300 alone{lab}{s_}", col)
                for m, lab, col in (("", "", "#74c476"), ("_motifact_assoc_noctcf", ", ATAC TF list", "#fd8d3c"))
                for s_ in ("", "_seed1", "_seed2")] + [
    ("k4me1obs_k562", NOATAC_D, "ATAC x observed H3K4me1", "#8c6bb1"),
    ("k4me1only_k562", NOATAC_D, "observed H3K4me1 alone", "#bcbddc"),
    ("k27obs_k562", NOATAC_D, "ATAC x observed H3K27ac (read counts)", "#525252"),
    ("k27obsonly_k562", NOATAC_D, "observed H3K27ac alone (read counts)", "#969696"),
    ("p300obsonly_k562", NOATAC_D, "observed p300 alone (read counts)", "#807dba")]
assert len(NOATAC_ARMS) == 20, len(NOATAC_ARMS)

# Multi-cell MSE-only model trained with three ENCODE tissues added (1.41), 3 runs; the MSE-only
# and assoc_noctcf seed arms ride along from their own run for seed-averaged pairs.
TISSUE_D = f"{ABC}/2026_1002_p300_tissue"
TISSUE_ARMS = [x for x in SEEDS_ARMS if "twostage" not in x[0] and "oracle" not in x[0]] + [
    (f"p300pred_multicell_tissues_holdK562{s_}", TISSUE_D, f"multi-cell + 3 tissues{' ' + s_[1:] if s_ else ''}", "#08519c")
    for s_ in ("", "_seed1", "_seed2")]

# Panel-mean (Corgi idea 1) model (1.42), 3 runs, and the panel mean alone; the MSE-only
# and assoc_noctcf seed arms ride along for seed-averaged pairs (decisions [2026-10-13]).
PANELMEAN_D = f"{ABC}/2026_1013_p300_panelmean"
PANELMEAN_ARMS = [x for x in SEEDS_ARMS if "twostage" not in x[0] and "oracle" not in x[0]] + [
    (f"p300pred_multicell_panelmean_holdK562{s_}", PANELMEAN_D, f"panel-mean model{' ' + s_[1:] if s_ else ''}", "#7a0177")
    for s_ in ("", "_seed1", "_seed2")] + [
    ("p300pred_panelmeanonly_holdK562", PANELMEAN_D, "panel-mean p300 alone (no model)", "#f768a1")]

# Bolero idea 2 (1.43, decisions [2026-10-21]): FILM (K562 row), CLOSEST (GM12878 row) and the
# matched BASE on panel 13, 3 runs each; Corgi 5-cell (F-067) read from its own ABC run, so
# BASE - Corgi5 is paired on this pair set.
FILM13_D = f"{ABC}/2026_1021_p300_film13"
_F13 = "p300pred_multicell_panelmean_film_chromvar13pc_panel13_holdK562"
FILM13_ARMS = P300_ARMS[:2] + [
    x for x in PANELMEAN_ARMS if x[0].startswith("p300pred_multicell_panelmean_holdK562")] + [
    (f"{_F13}_row{r}{s_}", FILM13_D, f"{lab}{' ' + s_[1:] if s_ else ''}", col)
    for r, lab, col in (("K562", "FiLM panel 13, K562 row", "#006eae"),
                        ("GM12878", "FiLM panel 13, GM12878 row (closest)", "#9bcae9"))
    for s_ in ("", "_seed1", "_seed2")] + [
    (f"p300pred_multicell_panelmean_panel13_holdK562{s_}", FILM13_D,
     f"Corgi panel 13, no conditioning{' ' + s_[1:] if s_ else ''}", "#c5373d")
    for s_ in ("", "_seed1", "_seed2")]
assert len(FILM13_ARMS) == 14, len(FILM13_ARMS)

if _a.arms == "motifact":
    ARMS = MOTIFACT_ARMS
elif _a.arms == "motifact2":
    ARMS = MOTIFACT2_ARMS
elif _a.arms == "teach":
    ARMS = TEACH_ARMS
elif _a.arms == "seeds":
    ARMS = SEEDS_ARMS
elif _a.arms == "rna":
    ARMS = RNA_ARMS
elif _a.arms == "frag":
    ARMS = FRAG_ARMS
elif _a.arms == "k27frag":
    ARMS = K27FRAG_ARMS
elif _a.arms == "noatac":
    ARMS = NOATAC_ARMS
elif _a.arms == "tissue":
    ARMS = TISSUE_ARMS
elif _a.arms == "panelmean":
    ARMS = PANELMEAN_ARMS
elif _a.arms == "film13":
    ARMS = FILM13_ARMS
elif _a.arms == "gata":
    ARMS = GATA_ARMS
elif _a.arms == "ranksweep":
    ARMS = RANKSWEEP_ARMS
elif _a.arms == "nbrx":
    ARMS = NBRX_ARMS
elif _a.arms == "nbr":
    ARMS = NBR_ARMS
elif _a.arms == "rank":
    ARMS = RANK_ARMS
elif _a.arms == "oracle":
    ARMS = ORACLE_ARMS
elif _a.arms == "multicell":
    ARMS = MULTICELL_ARMS
elif _a.arms == "p300fix":
    ARMS = P300FIX_ARMS
elif _a.arms == "slopefix":
    ARMS = SLOPEFIX_ARMS
elif _a.arms == "accessibility":
    ARMS = ACC_ARMS
elif _a.arms == "p300":
    ARMS = P300_ARMS
elif _a.arms == "p300all":
    ARMS = P300_ALL_ARMS
TAG = {"h3k27ac": "predicted_activity", "p300": "p300_activity",
       "p300all": "p300_all", "accessibility": "accessibility_activity",
       "dnaseinput": "dnase_input_activity",
       "multitask": "multitask_activity",
       "slopefix": "slopefix",
       "p300fix": "p300_corrected",
       "multicell": "p300_multicell",
       "oracle": "p300_oracle",
       "rank": "p300_rank",
       "nbr": "p300_nbr",
       "nbrx": "p300_nbrx",
       "ranksweep": "p300_ranksweep",
       "gata": "p300_gataoracle",
       "motifact": "p300_motifact",
       "motifact2": "p300_motifact2",
       "teach": "p300_teach",
       "seeds": "p300_seeds",
       "rna": "p300_rna",
       "frag": "p300_frag",
       "k27frag": "k27frag_activity",
       "noatac": "noatac_activity",
       "tissue": "p300_tissue",
       "panelmean": "p300_panelmean",
       "film13": "p300_film13"}[_a.arms]
RUN = {"h3k27ac": "2026_0904_predicted_activity", "p300": "2026_0905_p300_activity",
       "p300all": "2026_0906_p300_all",
       "accessibility": "2026_0910_accessibility_activity",
       "dnaseinput": "2026_0915_dnase_input_activity",
       "multitask": "2026_0916_multitask_activity",
       "slopefix": "2026_0917_slopefix",
       "p300fix": "2026_0919_p300_corrected",
       "multicell": "2026_0922_p300_multicell",
       "oracle": "2026_0923_p300_oracle",
       "rank": "2026_0924_p300_rank",
       "nbr": "2026_0924_p300_nbr",
       "nbrx": "2026_0924_p300_nbrx",
       "ranksweep": "2026_0925_p300_ranksweep",
       "gata": "2026_0926_p300_gataoracle",
       "motifact": "2026_0927_p300_motifact",
       "motifact2": "2026_0928_p300_motifact2",
       "teach": "2026_0928_p300_teach",
       "seeds": "2026_0929_p300_seeds",
       "rna": "2026_0930_p300_rna",
       "frag": "2026_0930_p300_frag",
       "k27frag": "2026_0930_k27frag",
       "noatac": "2026_1001_noatac",
       "tissue": "2026_1002_p300_tissue",
       "panelmean": "2026_1013_p300_panelmean",
       "film13": "2026_1021_p300_film13"}[_a.arms]

BASELINES = [
    ("distToTSS",      "FALSE", "mean", "Inf", "TRUE",  "Distance to TSS",      "#c5cad7"),
    ("nearestTSS",     "TRUE",  "max",  "0",   "FALSE", "Nearest TSS",          "#6e788d"),
    ("nearestGene",    "TRUE",  "max",  "0",   "FALSE", "Nearest gene",         "#435369"),
    ("within100kbTSS", "TRUE",  "max",  "0",   "FALSE", "Within 100kb of TSS",  "#1c2a43"),
]

missing, truncated = [], []
for arm, root, _n, _c in ARMS:
    p = f"{root}/{arm}/{PRED_FILE}"
    if not os.path.exists(p):
        missing.append(p)
        continue
    # EXISTENCE IS NOT COMPLETENESS. On 2026-09-17 this check passed two arms whose
    # create_predictions job was still streaming into the file: 71 MB present, gzip
    # truncated. The comparison would then have scored partial data silently.
    if subprocess.call(["gzip", "-t", p], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL) != 0:
        truncated.append(p)
if missing:
    raise SystemExit("missing prediction files:\n  " + "\n  ".join(missing))
if truncated:
    raise SystemExit("prediction files present but TRUNCATED (still being written?):\n  "
                     + "\n  ".join(truncated))

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
