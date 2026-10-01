# TODO

Tasks only. Results live in `.living/findings/`, decisions in `.living/decisions.md`. Notes from the
pre-standard list (2026-10-01) moved to `.living/decisions.md` under "[2026-10-01] Notes migrated
from todo/TODOLIST.md"; the old file is kept as `todo/TODOLIST.pre-standard.md`.

## Now
<!-- At most 5 items actually being worked on this week. -->
- [ ] P1 Read out the no-ATAC and H3K4me1 benchmark and write F-044 #benchmark (F-043)
  - ABC half done: `2026_0824_H3K27ac_model/results/2026_1001_noatac_auprc_pairs.tsv`; rE2G 46228392 -> `re2g_k4_boot` 46228394 pending
  - launcher `2026_0824_H3K27ac_model/log/noatac_launch.sh`
- [ ] P1 Process the H3K4me1 panel: filter the BAMs, build 5' bigwigs and the training panel #data
  - `0.47` array 46238895 -> `data/h3k4me1_panel/`; filtered new-cell BAMs stay in `$SCRATCH/h3k4me1/<cell>/`
- [ ] P1 Measure H3K4me1 beyond ATAC in the five cells and its cross-cell consistency, before training #h3k4me1
  - `4.60` job 46239127, after 46238895 -> `results/2026_1001_h3k4me1_beyond_atac.*`
- [ ] P1 Add the tissue p300/ATAC panel (stomach, transverse colon, ovary) to multi-cell training; score on K562, 3 runs #data (F-038)
  - stomach is the only candidate with GATA activity (4.49/4.52); data `$SCRATCH/p300_tissues`, `reference/p300_tissue_manifest.tsv`

## Next
<!-- Ordered: top item is the next thing to start. -->
- [ ] P2 Train H3K4me1 models (in-cell K562, GM12878 -> K562, multi-cell holding out K562) if the beyond-ATAC check shows signal #h3k4me1
  - bar is high: observed K562 H3K4me1 in ABC is no better than ATAC alone and 0.060 below observed H3K27ac (F-044)
- [ ] P2 Add H3K4me1 as a separate rE2G feature (`external_features_config`) beside ATAC and p300 #benchmark
- [ ] P3 Build a two-target model (p300 + H3K27ac heads) only if an H3K27ac head would transfer beyond ATAC #model (F-043)
- [ ] P3 Run the rE2G promoter-feature ablation (drop `normalizedH3K27ac_prom`, `combinedActivity_prom`) #benchmark (F-029)

## Later
- [ ] P3 Run motif syntax analysis (SHAP / TF-MoDISco / FiNeMo) on the residual-trained model, +/-500 bp window #motif
- [ ] P3 Run SHAP on the K562 ATAC multimodal model (separate from the residual-model run) #motif
- [ ] P3 Compare the 26 GM12878 MoDISCo motifs (`modisco/max_seqlets_250k_30_10_0/logos/`) to the K562 set #motif
- [ ] P3 Confirm the identity of `REPEAT_G`, `NF2L_NFE` and `ELF` in the v2 FiNeMo motif set #motif
- [ ] P3 Update `P300_INTERACTORS` in `plot_finemo_composite_figure.py` after a BioGRID review #motif
- [ ] P3 Redefine p300+ across figures (about 40% of `EP300_peak_overlap` elements have `true_logcounts = 0`) #report
  - candidate: top 20% by observed p300 counts; affected figures listed in the migrated notes
- [ ] P3 Make the two BioRender schematics, Fig 1a (BPNet) and Fig 2a (multimodal) #report
- [ ] P3 Compute ATAC ChromBPNet true CV performance for Fig 1d (placeholder r = 0.70, no p300+ number) #report
- [ ] P3 Finish section 3 figures: MoDISCo logos and FiNeMo hit panels #report
  - done already: `scripts/plot_individual_motif_insertions.py`, `scripts/plot_motif_pair_heatmaps.py`
- [ ] P3 Add an accessibility-stratified metric (Pearson and predicted/observed by ATAC decile) to the `2.x` eval #benchmark (F-021)
- [ ] P3 Measure counts accuracy on CRISPR-tested elements per arm, not only genome-wide #benchmark (F-021)
- [ ] P3 Use partial rho given ATAC (`4.35`) as the upstream screen, over all candidates, never tested elements #benchmark (F-029)
- [ ] P3 Test removing qnorm in ABC as its own arm (`use_qnorm: False`) #benchmark
- [ ] P3 Decide whether to make test-time RC averaging the default and re-score the tables as a set #benchmark
- [ ] P3 Re-measure the dynamic range of the DNase-input H3K27ac prediction (`4.9`/`4.10` on `predk27ac_k562_dnase.bw`) #h3k27ac (F-018)
- [ ] P3 Test whether the DNase-input ABC arm transfers (GM12878 DNase-input model onto K562, one `4.4` call) #h3k27ac (F-018)
- [ ] P3 Run the indicator-channel control: GC and CpG tracks as extra channels (`0.27` builds them) #h3k27ac
- [ ] P3 Characterise the 12 threshold-level false negatives (observed H3K27ac 39x the genome median) #h3k27ac (F-007)
- [ ] P3 Decide the profile head (drop / keep 1 bp / bin to ~50 bp), jointly with a powered clw sweep #h3k27ac
- [ ] P3 Re-run the `count_loss_weight` sweep with >= 3 folds and >= 2 seeds #h3k27ac
- [ ] P3 Re-check +/-500 vs +/-1000 (and the counting window) on the 5' target at the next retrain #h3k27ac
- [ ] P3 Train on the ATAC-derived K562 element set (md5 7d5995ce; path in the migrated notes) #h3k27ac
- [ ] P3 Switch the accessibility input to 5' counts as a set (~35 fold-jobs) #h3k27ac
- [ ] P3 Decide whether to re-run the transfer and deployment comparisons at `n_layers` 10 #h3k27ac
- [ ] P3 Reweight the accessible-but-unacetylated quadrant (up-weight, not relabel) #h3k27ac
- [ ] P3 Sweep `negative_ratio` (0.1) and `--max-negatives` (50,000) for H3K27ac #h3k27ac
- [ ] P3 Re-run the depth thinning at a second seed (`0.41 --seed 1`) if F-017 goes into a figure #h3k27ac (F-017)
- [ ] P3 Predict the composite geomean(DHS, H3K27ac) directly; design caution in the migrated notes #h3k27ac
- [ ] P3 Repeat the residual-objective comparison for p300 #model
- [ ] P3 Train a pretrained sequence trunk (Borzoi or Enformer embeddings, frozen then fine-tuned) #model
- [ ] P3 Try a multiplicative two-tower model (sequence potential x accessibility availability) #model
- [ ] P3 Add gradient reversal on a cell-type classifier over the sequence embedding #model
- [ ] P3 Train multi-task across marks available in six cell types (H3K4me1, H3K27me3, CTCF) #model
- [ ] P3 Train TeloHAEC and score transfer, then +/-IL1b / +/-TNFa / -VEGF inference #data
  - caveats (ATAC-derived elements, 36 bp reads, depth) in the migrated notes
- [ ] P3 Assess hashFrag chromosome splits for train/test sequence-similarity leakage #data
- [ ] P3 Write `training_target.json` for the June GM12878 p300 model (`2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac`) #infra
- [ ] P3 Make submit scripts refuse to overwrite a completed fold directory #infra
- [ ] P3 Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history is confirmed #infra
- [ ] P3 Remove the stray directory `Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` #infra
- [ ] P3 Sort experiment bedGraphs before merging in `scripts/0.3.make_training_bw.sh` #infra

## Waiting
<!-- Blocked on a job, a person or data. Say on what. -->

## Inbox
<!-- Untriaged: questions, ideas, links. One line each, dated. -->
- 2026-09-07 question: is a mostly-accessibility model useful for the intended application, or does the goal need the sequence component?

## Settled
<!-- Decisions not to revisit, one line each, with the finding that settled it. -->
- Target = 5' ends, +/-500 bp window, `count_loss_weight = 10`; PE H3K27ac targets use read 1 only
- Panel and deployment input are ATAC only; DNase-input models and a DNase panel are off-direction (F-016)
- The K562 CRISPR benchmark is the only benchmark; no search for a second cell type (Maya, 2026-10-01)
- AlphaGenome on hold: every released model trained on K562 EP300 and H3K27ac (Maya, 2026-09-23)
- TF-channel and TF-list variants plateau within run noise; no further estimators, teaching losses, seed replicates, finemo TF lists or leave-one-out Test 1 (F-034 to F-038, F-041; Maya, 2026-10-01)
- RNA TF lists do not beat the ATAC-derived list; per-TF factorised model and MotifCompendium remap dropped (F-041)
- Neighbour-graph re-ranking is null in ABC (F-031)
- Fragment-size input channels are null on transfer; do not combine, replicate or retry (F-005, F-042)
- Flank nucleosome fragments are not longer at active elements (F-040)
- Correcting transferred H3K27ac with fragment features does not reach the benchmark (F-043)
- Predicted p300 keeps the ATAC term in ABC: dropping it costs 0.021 [0.005, 0.038]; no no-ATAC re-scoring of deployable arms (F-044)
- Residual objective gives no transfer gain (F-023)
- Accessibility-normalisation fix is not the mechanism and is mildly harmful; do not adopt (F-022, F-023)
- Monotone recalibration of the predicted track dropped; reopen only as a study of distribution shape (F-021)
- Training longer (100 epochs) is not attributable; do not adopt (F-019)
- ATAC-to-DNase converter is finished as a model input; converter count accuracy not pursued (F-014)
- Sequence gate and asymmetric loss closed; no p300 gate (decision 2026-09-05)
- Predicted p300 as an input channel to an H3K27ac model is redundant; not pursued
- Wide receptive field and fragment channels are sub-additive; prefer wide + flat (2026-09-14)
- Depth-match inputs before any cross-cell-type comparison (F-017)
- A floor must use the same assay as the model it bounds (F-016)
- Any change to the loss needs an epoch-matched comparator (F-015)
- Use explicit accession lists, never a directory glob (TeloHAEC_ctrl/ATAC holds 3 EA.hy926 files)

## Done
<!-- Last two weeks only; older items live in git and findings. -->
- [x] P1 Benchmark observed and predicted p300 x H3K27ac, with no-ATAC controls (F-043)
- [x] P1 Benchmark fragment-size ATAC channels, 3 runs (F-042)
- [x] P1 Benchmark RNA TF lists, 3 runs (F-041)
- [x] P2 Test flank fragment length at H3K27ac-high elements (F-040)
- [x] P2 Fragment features vs H3K27ac and the profile figure (F-039)
- [x] P1 Seed replicates and frozen two-stage model (F-038)
- [x] P2 Teach the model to use the TF channel: head, pair, headpair (F-037)
- [x] P1 Deployable TF channel `assoc_noctcf` (F-036)
- [x] P1 TF-activity track with oracle weights, Test 1 (F-035)
- [x] P2 Learned TF-list estimators (F-034)
- [x] P2 TF-activity recovery from ATAC (F-033)
- [x] P2 Error anatomy and GATA oracle (F-032)
- [x] P2 Ranking loss, weights 2 to 100, through ABC and rE2G (F-030)
- [x] P2 Oracle decomposition of the transfer gap (F-028)
- [x] P2 rE2G with p300 in the activity slot (F-029)
- [x] P2 Multi-cell joint training holding out K562 (F-027)
- [x] P3 Quantify run-to-run variance (F-038)
