#!/bin/bash
# Command log — H3K27ac sequence vs sequence+ATAC models (K562)
# Append every command actually run, in order. This file is the provenance record.

D=$OAK/Users/sheth/EP300_BPNet
PROJ=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python   # pixi not on PATH on Sherlock; call env python directly

### 2026-08-24 — Step 0.1: pick the counting window empirically ###
# H3K27ac sits on flanking nucleosomes, so the p300 window (out=1000 / in=2114) may
# be too narrow. Profile signal vs distance from DNase candidate element centers and
# weigh signal gained against neighbour-element contamination.
#
# Data: ENCSR000AKP (K562 H3K27ac), 2 replicate BAMs merged into a raw-count,
# 250bp-fragment-extended coverage BigWig by $OAK/Users/sheth/Data/scripts/bam_to_bigWig.sh
# (see $OAK/Users/sheth/Data/ENCODE/log.sh). Not normalized — correct for training.
H3K27AC_BW=$OAK/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw
ELEMENTS=$D/reference/K562_DNase_candidate_elements.narrowPeak

cd $PROJ
$PY scripts/0.1.profile_signal_vs_distance.py \
  --elements  $ELEMENTS \
  --signal-bw $H3K27AC_BW \
  --outdir results --figdir figures
# -> results/signal_vs_distance.tsv, results/window_tradeoff.tsv
# -> figures/h3k27ac_signal_vs_distance.{pdf,png}
# Result: bimodal, shoulders at -275/+275 bp, central dip, distal plateau 17.1
#         (3.48x enrichment), signal ~background by +/-2000 bp.
#         Neighbour contamination: 0% at +/-500, 9.9% at +/-750, 19.9% at +/-1000.
# Peak memory 264 MB — safe on a login node. Do NOT raise --flank without keeping the
# binned extraction; a full-resolution +/-10kb matrix OOMs.

### 2026-08-24 — Trainer change: unstranded target support ###
# scripts/train_multimodal_bpnet.py: --signal-minus-bw and --accessibility-bw are now
# optional. n_outputs is DERIVED (2 when a minus track is given, 1 when not) rather than
# hardcoded to 2, and --max-negatives caps the negative pool.
# Verified the stranded p300 path is bit-identical to before (regression test compared
# extract_windows output element-by-element against the pre-patch file). RC augmentation
# needed no change: torch.flip(yi,[0,1]) on a 1-channel target is a no-op in dim 0.

### 2026-08-24 — Step 0.2/0.3: inter-replicate ceiling (runs alongside training) ###
sbatch scripts/0.2.make_replicate_bigwigs.sh          # job 40697417
# -> data/h3k27ac_rep1.bw, data/h3k27ac_rep2.bw  (raw counts, -fs 250, matching target)
# then:
# $PY scripts/0.3.replicate_ceiling_by_window.py \
#     --elements $D/reference/K562_DNase_candidate_elements.narrowPeak \
#     --rep1-bw data/h3k27ac_rep1.bw --rep2-bw data/h3k27ac_rep2.bw \
#     --outdir results --figdir figures

### 2026-08-24 — Step 1.1: training ###
# Smoke test first — 30 jobs failing identically is a waste of queue time.
sbatch --job-name=k27_smoke_seq_hw500 scripts/1.1.submit_training.sh sequence 0 500  # job 40697384
# Once that clears data extraction and starts epoch 1, submit the rest:
#   bash scripts/1.2.submit_all.sh --dry-run   # 30 jobs: 3 modes x 2 windows x 5 folds
#   bash scripts/1.2.submit_all.sh
#
# Windows under test: hw=500 (out 1000 / in 2114, zero neighbour contamination) and
# hw=1000 (out 2000 / in 3114, ~20% contamination but captures most of the enriched
# region). COUNT_LOSS_WEIGHT defaults to 10 to down-weight the profile head; CHECK the
# reported profile vs count losses in the log and retune — the profile loss roughly
# doubles going from out_window 1000 to 2000, so the same weight is not equivalent
# across the two windows.
