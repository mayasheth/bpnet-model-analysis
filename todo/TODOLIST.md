# TODO

Tasks only. Results live in `.living/findings/`, decisions in `.living/decisions.md`. Notes from the
pre-standard list (2026-10-01) moved to `.living/decisions.md` under "[2026-10-01] Notes migrated
from todo/TODOLIST.md"; the old file is kept as `todo/TODOLIST.pre-standard.md`.

## Now
<!-- At most 5 items actually being worked on this week. -->
- [x] P1 Baseline scE2G training run on Xu K562 multiome #sce2g
  - driver 46428971 (2026-10-03): CRISPR AUPRC multiome_powerlaw_n6 0.628 [0.575, 0.672], scATAC_powerlaw 0.579 [0.524, 0.623], distance 0.414 [0.363, 0.461]; 10.4% of CRISPR elements missing; `$SCRATCH/sce2g_runs/2026_1002_k562_xu/performance_across_models.tsv`
  - needed: #118 `source_path` fix (local branch `local-with-pr118`), `--rerun-triggers mtime`, and no empty `external_features_config` column (R NA bug)
- [x] P1 Feature B: published-method TF footprints (TOBIAS, PRINT) as scE2G features #sce2g: CLOSED, null on CRISPR
  - both pass ChIP (GATA1 0.674 TOBIAS / 0.605 PRINT within accessibility deciles) and the p300 gate (~0.1 beyond motif), neither moves AUPRC: TOBIAS -0.0036 [-0.0090, +0.0019], PRINT -0.0020 [-0.0061, +0.0026] vs baseline (F-048 to F-051)
  - reopen only with a footprint several times stronger against p300 (seq2PRINT trained on this sample); the ABC-activity version is the item below
  - homemade footprints (`scripts/1.1`, `1.2`) fail on GATA1 (F-047); BINDetect per-motif means are GC-biased
- [x] P1 Feature A: peak-TF Kendall (`scripts/5.0`, `5.1`), gate vs p300 #sce2g: considered and DEFERRED 2026-10-04, too weak for scE2G
  - strict gate (beyond ATAC, motif count, permuted-TF and any-TF controls) 0.015 [0.009, 0.019], + ENCODE ATAC 0.028 [0.022, 0.034]; a fifth of the footprint signal that was null on CRISPR (F-052)
  - reopen only in a heterogeneous multiome sample (K562 has little TF expression variation across cells)
- [x] P1 Footprints inside the ABC activity term (`scripts/4.9`, exact reweighting) #sce2g: CLOSED 2026-10-05, null on CRISPR
  - lam 1: TOBIAS +0.0011 [-0.0063, +0.0087], PRINT +0.0013 [-0.0061, +0.0088]; lam 4 -0.010 / -0.011 (F-053)
  - reopen only if the PRINT lam 1 point estimate replicates in a second sample or cell
- [ ] P1 Confirm the GC / CTCF gains (F-055) on pairs they were not chosen on: ENCODE-rE2G with GC + CTCF features (pre-registered), and another CRISPR set if one exists #sce2g
  - done: GC within non-promoter elements, AUPRC 0.686 -> 0.700 (F-055 caveats, 6.1 rerun)
  - rE2G test considered and DEFERRED 2026-10-05 (Maya: stay on scE2G evals); driver 46625263 cancelled mid-run. Built and ready: `4.62` features (153,360 July elements, 41% with a CTCF site), 4 arms in `4.32`, `log/gcctcf_launch.sh`. Reopen with `log/gcctcf_launch.sh` minus the features job; the driver needs `--rerun-incomplete` after the cancel
  - pass rule fixed in decision [2026-10-05]: ATAC_GCCTCF - ATAC interval above 0 and both coefficients negative; shares ~470 positives with F-055, so a pass shows transfer across model families, not freedom from selection
  - 2026-10-08: in scE2G, GC adds nothing once CpG-island distance and the motif features are in (F-059); an rE2G rerun should test CpG-island distance alongside GC
- [x] P1 Observed-p300 ceiling and exact-reweight rerun of footprints in ABC #sce2g: observed p300 +0.056 to +0.060 (F-056); footprints in ABC null, lam 1 +0.001 (F-053, exact)
- [ ] P1 Find deployable evidence for "accessible but inactive" elements, the part of the p300 ceiling (F-056) that matters: e.g. predicted p300 from the in-cell / multi-cell models as a scE2G feature, more sequence composition terms #sce2g
  - predicted p300 CLOSED 2026-10-05, null (F-057): multi-cell (deployable) +0.0012 [-0.0056, +0.0081], in-cell +0.0092 [+0.0008, +0.0180], ENCODE-ATAC control -0.0049; observed p300 is +0.058 above the deployable arm. Reopen only with a p300 model whose in-cell version clears ~+0.03 in scE2G
  - next candidate source: the feature screen below
- [x] P1 Feature screen on observed p300 beyond Xu ATAC, then one pre-registered scE2G run (`A2/plans/7_feature_screen.md`, `scripts/7.0`-`7.3`) #sce2g: CLOSED 2026-10-05, fails the rule (F-058)
  - p300: short-fragment share partial rho 0.19 [0.18, 0.20]; trees recover 0.25 of p300 beyond ATAC + composition
  - CRISPR vs n6 + GC: short-fragment share +0.0054 [-0.0001, +0.0109], composite -0.0011 [-0.0018, -0.0004]
  - reopen only with a second cell or CRISPR set (for the short-fragment share), or a selection target closer to CRISPR than p300
- [x] P1 Direct CRISPR feature discovery with selection inside nested chromosome CV + within-gene label-permutation null (`A2/plans/8_crispr_feature_discovery.md`) #sce2g: DONE 2026-10-08, pre-registered PASS (F-059)
  - Phase 0 done 2026-10-06: offline replica `scripts/8.0` matches scE2G AUPRC within 0.0005 for 4 models (per-pair scores within 0.04; compare models within the replica)
  - 8.1 (2026-10-06): ARC split into ABC + exp(Kendall) is within 0.004 of n6 [-0.015, +0.006]; ARC keeps Kendall's sign (log ARC = log ABC + w tau), exp(tau) matches; base for Phase 2
  - Phase 1 done: `scripts/8.2`, 564 candidates in 12 families (528 motif clusters), `$SCRATCH/crispr_discovery/`
  - Phase 2 done (`scripts/8.3`, job 46825675): nested gain vs base 0.627: forward +0.027 [+0.011, +0.042], family +0.023 [+0.011, +0.034], L1 +0.024 [+0.000, +0.049], trees +0.008; consistent picks: CpG-island log distance 5/5, GATA6 motif cluster 5/5, IKZF2 motif 4/5; families alone: CpG island +0.023, composition +0.017, rest within +/-0.006, all motifs -0.034
  - Phase 3: within-gene label shuffle INVALID (Maya, 2026-10-06: 218 of 429 positives in genes it cannot change; moves positives distal, median 28 -> 45 kb); replaced by shuffling candidate vectors within TSS-distance bins (`8.4`, array 46836062, `8.5` -> `results/2026_1006_nested_null.tsv`)
  - null result: forward +0.027 vs null max +0.002, family (CpG island) +0.023 vs max +0.004, composition family +0.017 vs max +0.003; all p = 1/101; base unchanged in every null run
  - distance control (46913125, `8.3 --base-logdist`): gains hold, forward +0.032 [+0.015, +0.048], CpG-island family +0.025; same picks: CpG-island log distance 5/5, GATA cluster (GATA6/Gata3/TRPS1) 5/5, IKZF2 4/5; all three raise the regulated rate (GATA site 0.043 -> 0.111)
  - promoter-side CpG (Maya's idea, `8.7`, replica LOCO, hypothesis fixed in the docstring): null; CGI-island promoter overlap +0.001 [-0.004, +0.006], 3 promoter features -0.004 [-0.010, +0.003], on top of the 3 element features +0.005 [-0.003, +0.013]; CGI promoters slightly MORE often regulated (0.049 vs 0.043), opposite of the hypothesis
  - Phase 4 PASS (decisions [2026-10-07]; `8.6`, driver 47016567 after the Oak-outage SIGBUS of 46915679, bootstrap 47016570 -> `results/2026_1007_discovery_paired_bootstrap.tsv`): n6 + the three vs n6 +0.036 [+0.022, +0.049] (0.632 -> 0.668); CpG-island distance alone +0.024 [+0.014, +0.033]; GC adds nothing on top of the three (+0.036 either way); in-sample, so quote the nested ~+0.021 over n6 as expected gain
- [x] P1 External test of the three F-059 features (CpG-island log distance, GATA cluster, IKZF2) as scE2G defaults #sce2g: DONE 2026-10-12, pre-registered FAIL (F-063)
  - eQTL transfer (plan 9, `9.0`-`9.2`, decisions [2026-10-12]): 62 matched immune pairs, the three R1 -0.34 [-1.15, +0.11], R2 +0.11 [-0.18, +0.43]; CGI alone WORSE (R1 -0.56 [-1.14, -0.13]); no gain on K562's own eQTLs either. Not adopted; reopen only with a CRISPR set in a second cell showing a gain
- [x] P1 Positive control for the eQTL transfer benchmark #sce2g: DONE 2026-10-13 (F-064): scATAC - multiome R1 -2.63 [-3.48, -1.79], R2 -1.72 [-2.27, -1.12] for -0.049 AUPRC on CRISPR; F-063's upper bounds exclude the gain F-059 predicts
  - [x] activity-only control (F-065, 2026-10-13): ARC-E2G -> ABC moves R2 only -0.35 [-0.47, -0.20] (R1 -0.40 [-0.78, +0.01]); the scATAC drop is mostly RNA; F-063 is inconclusive for activity effects of F-059's size, not a refutation
- [x] P1 Find a transfer benchmark that weighs element activity: a CRISPR set in a second cell, or an eQTL readout restricted to expressed genes / elements where the activity term decides the ranking (calibrate with the ABC arm, F-065) #sce2g: eQTL route CLOSED 2026-10-13, not calibrated (F-066); the second-cell CRISPR route is in Waiting
  - 2026-10-13 Maya chose the eQTL route: plan 10 pre-registered (decisions [2026-10-13]): expressed genes, non-promoter candidates <= 1 Mb, E1 pooled AUPRC + E2 within-gene AP; calibration on n6_abc with n6_three / n6_cgi blinded (`10.1`, `10.2`)
  - result (clusters 47687538, gather 47687541, diagnostics `10.3`): ABC swap E1 -0.0002 [-0.0023, +0.0021], E2 +0.0007 [-0.0027, +0.0040]; n6 and ABC rank a gene's candidates alike (within-gene Spearman 0.994); eQTL positives are the proximal top element (median 22 kb). n6_three / n6_cgi stay blinded
  - reopen only with a readout that a CRISPR-scale activity swap moves before any test arm is read
- [ ] P1 Corgi idea 1: multi-cell p300 model that predicts the change from the panel mean (`train_multimodal_bpnet.py --panel-mean`, `A/scripts/1.42`; pre-registered decisions [2026-10-13]) #model
  - leave-one-out panel mean in training (all four cells at K562); inputs ATAC z, ATAC z - panel mean z, panel-mean p300; count head through a log1p panel-mean offset. Judge: seed-averaged ABC vs MSE-only (0.4828), rE2G secondary (0.5579); panel mean alone as a descriptive arm
  - 2026-10-13: helper verified against direct recomputation + invalid-window guard; smoke 47696615 submitted; downstream family `panelmean` wired into 4.1/4.3/4.4/4.6/4.32; scoring chain `A/log/panelmean_score_launch.sh`
- [ ] P2 Plan 9 eQTL rescoring as a cheap screen for K562-CRISPR-selected features: catches harm (CGI alone lowers R1) and RNA-side effects, but cannot confirm activity gains below about 0.05 AUPRC (F-065) #sce2g
  - the GATA term is likely erythroid-specific; a non-erythroid cell tests whether CpG-island distance alone carries over
  - label-free proxy done 2026-10-09 (`8.8`, F-060, Maya's idea): CGI distance tracks p300 beyond ATAC in K562 +0.233 and 3 of 4 other cells (pre-set rule: shared; MCF-7 +0.059 below half); the shared part is the near end (nearest decile under-acetylated relative to ATAC in all five, and depleted for regulated K562 pairs, 0.026 vs 0.031-0.060)
  - [x] near-CGI encoding (`8.9`, F-061, 2026-10-09): nearest-decile flag keeps 23% of the gain, +0.005 [-0.001, +0.012] vs +0.023 [+0.014, +0.033]; keep the percentile encoding for the external test and report by cell, since the far-end trend that carries the gain is cell-dependent (F-060)
  - [x] gene-density control (`8.10`, F-062, 2026-10-09): CGI distance keeps 57% of its gain on top of gene density, +0.013 [+0.005, +0.022]; gene density adds nothing once CGI distance is in. No need to carry it into the external test
  - report 7 Part C (Figs 12-15) records F-057 to F-062 (2026-10-09)
- [ ] P2 Decide the next direction now that multiome TF evidence (A, B, B-in-ABC) is null: e.g. observed p300 as a scE2G ceiling, or back to the H3K4me1 item in Next #sce2g (Maya)

## Next
<!-- Ordered: top item is the next thing to start. -->
- [ ] P1 Corgi idea 2: RNA TF list weighted by expression z-score across the panel + K562 instead of log TPM (`4.54` -> `4.45`, 3 seeds); after idea 1 is submitted #model
  - Corgi's per-cell TF importance tracks expression z-score (39-68% of variance); the log-TPM list was null (F-041). Same weighting could replace the hand-picked erythroid motif set in scE2G
- [ ] P2 Train H3K4me1 models (in-cell K562, GM12878 -> K562, multi-cell holding out K562) if the beyond-ATAC check shows signal #h3k4me1
  - bar is high: observed K562 H3K4me1 in ABC is no better than ATAC alone and 0.060 below observed H3K27ac (F-044)
  - beyond-ATAC gate passed (F-045), but observed H3K4me1 as separate rE2G features adds only +0.020 [+0.006, +0.031] over ATAC (F-044): Maya to decide
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
- [ ] Review of [EngreitzLab/scE2G#118](https://github.com/EngreitzLab/scE2G/pull/118) (`source_path` fix), requested 2026-10-03; run on `local-with-pr118` until merged (on reviewer)
- [ ] Is there a second CRISPR set (other cell type or held-out K562 screens) to confirm GC / CTCF on? (Maya)
- [ ] Lift the Settled "K562 CRISPR is the only benchmark" rule? A second-cell CRISPR set is now the only test of activity features in unseen cells (F-065, F-066) (Maya)

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
- Panel composition is null: adding three ENCODE tissues (incl. GATA-active stomach) +0.0015 [-0.0077, +0.0104]; no more cell types (F-046)
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
- [x] P1 Tissue p300/ATAC panel in multi-cell training, 3 runs: null, +0.0015 (F-046)
- [x] P1 rE2G with ATAC + observed H3K4me1 at E and P as separate features (F-044)
- [x] P1 Measure H3K4me1 beyond ATAC in five cells and across cells (`4.60`, F-045)
- [x] P1 Benchmark predicted p300 without the ATAC term and observed H3K4me1 in ABC and rE2G (F-044)
- [x] P1 Build H3K4me1 5' tracks for the five panel cells (`0.47`, `data/h3k4me1_panel/`)
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
