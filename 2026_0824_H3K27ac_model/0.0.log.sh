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
