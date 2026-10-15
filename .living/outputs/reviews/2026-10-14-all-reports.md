# Review — EP300_BPNet reports 1-7 and h3k27ac_model_report — 2026-10-14

**Scope**: the eight rendered analysis reports in `2026_0824_H3K27ac_model/` (`h3k27ac_model_report`, `report1`-`report7`, .qmd sources and HTML), the scripts that produce their figures and numbers (`A/scripts`, `A2/scripts`, `scripts/train_multimodal_bpnet.py`, `scripts/multimodal_bpnet.py`), and the `.living/` record they cite. Paths: `A/` = `2026_0824_H3K27ac_model/`, `A2/` = `2026_1002_sce2g_tf_activity/`.
**Files reviewed**: 8 reports, 66 figures, about 60 scripts
**Sub-agents run**: 9 (figure-visual split into two parts, reports 1-5 and 6-7). vault-fidelity self-skipped: no Vault pointer note exists (F63).

## Key decisions in this analysis

- **Verdict metric** — K562 CRISPR AUPRC (10,342 pairs, 466 regulated) through ABC and rE2G is the only verdict; correlation screens are not. See F2, F13, F34.
- **Downstream uncertainty** — paired percentile bootstrap over element-gene pairs, i.i.d.; "resolvable" = 95% interval excludes 0. See F7.
- **Seed replication** — verdicts need >= 3 runs against a matched comparator; most report 5/6 rows are still single runs against the lowest-scoring MSE run. See F3, F4.
- **Upstream split** — 5 chromosome folds; models early-stop on `val`, and every upstream metric is also scored on `val`; `test` is never read. See F13.
- **Genome-wide track assembly** — each chromosome painted by the lowest fold where it is val or test. See F15.
- **Held-out cell** — K562 excluded from multi-cell training, but its observed p300 guided the GATA target, CTCF removal and panel selection. See F14.
- **Accessibility input and normalisation** — log1p then per-cell z-score; full-interval coverage `atac.bw` for all p300 models; target stats over ABC candidates downstream, over K562 EP300 peaks upstream. See F16, F23.
- **AUPRC estimator** — step-wise AP (4.17, ABC), rE2G's `statistic_aupr` (4.34), sklearn trapezoid (A2 scE2G). See F42.
- **Promoter definition** — one collapsed TSS per gene (GENCODE v43, ±250 bp). See F20.
- **Mark comparison** — observed p300 vs observed H3K27ac as ABC activity without depth matching; this re-targeted the project. See F19.
- **Pre-registration** — Part C scE2G and eQTL tests carry pass rules written before scoring. See F2.
- **Multiple comparisons** — none applied across about 50 arms on one benchmark; informational, see F5 and F9.
- **Figure style** — `nature_style.py` for reports 1-4 figures; per-script styles for reports 5-7. See F78.

## Questions for the analyst

- Which reports are meant to be read by others now (lab, collaborators, a paper draft), and which are working records? That decides whether F1-F22 need fixing in place or a single "superseded" banner.
- Is the K562 CRISPR set the confirmation set as well as the selection set by design, or would you hold back a split (half the chromosomes) for adopted arms? (F5, F14)
- Should the borderline verdicts (lower bounds +0.0001 to +0.0012, including Corgi's +0.0010 from 2026-10-13) be re-read under a gene-level bootstrap before anything is adopted? (F7)
- Is the observed-p300-over-H3K27ac result (F-024) worth a depth-matched rerun, given that it set the project's target? (F19)
- The upstream tables were all scored on early-stopping chromosomes: rescore on `test`, or relabel them as validation numbers? (F13)

## Findings

### Statistics & causal inference

#### Major

##### F1. Report 5's "upstream does not predict downstream" rests on n = 5, a withdrawn row and difference-in-significance

`A/report5_multicell_p300_transfer.qmd:67, 100-114`
```markdown
upstream top-quintile Pearson and downstream AUPRC are uncorrelated (**r = -0.12, n = 5**). Four of the five disagree about whether the change helped.
| 100 epochs instead of early stopping | +0.0033 [-0.0064, +0.0130] | +0.0086 [+0.0005, +0.0175] | no |
```
**Why it matters here**: r = -0.12 has a 95% interval of about [-0.91, +0.85], and by sign alone 4 of 5 rows agree. The 100-epoch row is the only example of "discards changes that work", but F-019/F-020 withdrew the epoch attribution on 2026-09-17, six days before the report. Report 2 (L120) re-quotes the result as a scoring rule. Line 260 says run variance is unmeasured; F-038 has since measured it. (stats-causal, llm-failure-modes and doc-fidelity all flagged this.)
**Fix**: Restate as "no relationship shown (n = 5; 4/5 share a sign)", relabel the 100-epoch row "second training run, cause unresolved (F-019/F-020)", and carry the same wording into report 2.

##### F2. Report 7 replaces the pre-registered eQTL verdict with a post-hoc reading and leaves CGI-alone WORSE out of the summary

`A/report7_tf_evidence_multiome.qmd:2, 80, 589, 596`
```markdown
title: "... a gain eQTLs in unseen immune cells can neither confirm nor rule out"
# .living/decisions.md:1745  FAIL or WORSE -> the F-059 gain is K562-specific;
```
**Why it matters here**: The pre-registered verdict was FAIL (the three features) and WORSE (CGI distance alone, R1 -0.56 [-1.14, -0.13]). The reading "can neither confirm nor rule out" comes from a calibration arm added afterwards and extrapolated linearly. F-066 now supports the conclusion that eQTLs cannot calibrate activity effects. Still, on R1 the CGI-alone arm moved more than the calibration arm did, and in the harmful direction.
**Fix**: Lead with the pre-registered verdict and consequence, present F-065/F-066 as post-hoc calibration, and put CGI-alone WORSE on R1 in the summary.

##### F3. The seed-averaged ATAC TF-list gain includes the run that got it selected

`A/report7_tf_evidence_multiome.qmd:76, 205, 212`
```markdown
seed-averaged the ATAC TF list gains +0.0059 [+0.0001, +0.0132] in ABC and +0.0116 [+0.0038, +0.0202] in ENCODE-rE2G.
# results/2026_0930_seeds_auprc_pairs.tsv: ctcf runs 0.4910 / 0.4880 / 0.4872; mse runs 0.4810 / 0.4832 / 0.4842
```
**Why it matters here**: The arm's original run is its highest of three, and the comparator's original run is its lowest. On the two fresh seeds alone the gain is about +0.004 ABC and +0.006 rE2G. "A seventh of the gap" has a lower bound of +0.0001.
**Fix**: Report the seeds-1-2-only contrast beside the 3-run mean, and quote the fraction of the gap with its range or drop it from the title.

##### F4. Single-run arms are called "real gains" against the lowest MSE run

`A/report2_evaluation_methodology.qmd:225-231`
```markdown
**First prospective test: the screen missed a real gain.** ... Both gained.
In ABC, weight 10 beat the MSE-only model by +0.0087 [+0.0006, +0.0179]; in ENCODE-rE2G,
```
**Why it matters here**: Weight 10 is the best of 5 weights in 2 frameworks, from one run each. The comparator's original run sits 0.0027 (ABC) and 0.006 (rE2G) below the mean of its other two seeds. Adding the seed sd (F-038) to the bootstrap SE puts both lower bounds below 0. The same inflation applies to every single-run row in report 6 Fig. 1. Report 2 builds a methodological rule on this result. (stats-causal and data-pipeline both flagged this.)
**Fix**: Mark the rank-loss result provisional until it has 3 seeds against the 3-seed MSE mean, compare the remaining single-run rows against that mean, and drop "real gain".

##### F5. Unresolved intervals are written up as nulls and moved to "Closed, do not reopen"

`A/report3_design_decisions.qmd:98, 539, 1218`
```markdown
**Fragment channels do not travel**, despite a real in-cell gain (*p*=0.002).
- Fragment channels for deployment: null transferred in both directions.
| narrow + fragments | -0.0012 [-0.0185, +0.0162] *p*=0.86 | -0.0029 [-0.0267, +0.0208] *p*=0.75 |
```
**Why it matters here**: The transferred intervals contain the full in-cell gain (+0.0135). The same pattern recurs in several places:
- report 4's title "does not improve ABC" for +0.0208 [-0.0049, +0.0441];
- report 2:332 "clean null" for [-0.004, +0.022];
- report 6:79 "lose to the flat input" for -0.0046 [-0.0127, +0.0033];
- report 3:143 "refuted", while report 3:1090 says the effect is "not excluded".

With about 50 arms tried on one benchmark, items closed this way will not be revisited.
**Fix**: Write "not resolved; the interval includes <named effect>", and move such items from Closed to Deferred (underpowered), each with a reopening condition.

##### F6. Report 5 says no transferred model clears the floor, but its own single-source interval excludes zero

`A/report5_multicell_p300_transfer.qmd:66, 176, 253`
```markdown
- Both transferred models are close to the ATAC-only floor: multi-cell **+0.0130 [-0.0041, +0.0287]**, single-source **+0.0149 [+0.0012, +0.0278]**.
- **No transferred model has cleared an assay-matched floor downstream by a resolvable margin.**
```
**Why it matters here**: By the project's own rule, the target-normalised single-source arm clears the floor. The report also quotes two values for "single-source minus floor" upstream (+0.0106 in Fig. 1, -0.0022 in the table) without saying they use different normalisations. "The p300 gain is lost on transfer" (report 3:96) describes about 27% of the in-cell gain kept, with an interval spanning roughly 0-50%. (stats-causal and doc-fidelity both flagged this.)
**Fix**: Say transferred models keep about a quarter of the in-cell gain, with intervals at or just above zero, and label the normalisation on each single-source row.

##### F7. Every CRISPR interval resamples pairs independently, ignoring gene and locus clustering

`A/scripts/4.17.paired_auprc_bootstrap.py:82-83`
```python
for b in range(a.n_boot):
    idx = rng.randint(0, n, n)
# same in A/scripts/4.34.re2g_paired_bootstrap.py:96 and A2/scripts/4.2.paired_bootstrap.py:63
```
**Why it matters here**: The 429 regulated pairs sit in 297 genes, with tiled elements per locus. Within-gene correlation of paired differences makes these intervals too narrow. The borderline calls are the most exposed: the ATAC TF list (+0.0001), rank-loss w10 (+0.0006), the single-source floor (+0.0012), and the 2026-10-13 Corgi result (+0.0010 ABC, +0.0002 rE2G). The project already uses cluster bootstraps elsewhere (eGene, 1 Mb blocks).
**Fix**: Add a gene-level (or locus-block) resampling option to 4.17, 4.34 and A2/4.2, and re-read the borderline verdicts under it.

#### Minor

##### F8. Statistically identical estimates get opposite verdicts (RNA + ATAC vs ATAC TF list)

`A/report6_improving_transfer.qmd:269, 273`
```markdown
## RNA TF lists and fragment channels do not help; observed p300 and H3K27ac do best without ATAC in ABC
# rnaatac vs mse 0.00583 [-8.4e-05, 0.013402]; ctcf vs mse 0.005945 [0.00012, 0.013237]
```
**Why it matters here**: The verdicts flip on the third decimal of a lower bound, and the heading contradicts its own "ties the ATAC list".
**Fix**: Describe both lists as "about +0.006, at the edge of resolution", and keep "does not help" for the RNA-only list.

##### F9. A nominal p = 0.039 among six contrasts is headlined as "significantly worse"

`A/report3_design_decisions.qmd:1120-1130`
```markdown
**A:** Neither helps, and together they hurt on transfer: ... the combination is **significantly
| both | -0.0032 [-0.0157, +0.0092] *p*=0.51 | **-0.0128 [-0.0246, -0.0010]** *p*=0.039 |
```
**Why it matters here**: Bonferroni-adjusted across the 2x2 factorial in two settings, this is about p = 0.23, and it comes from one run.
**Fix**: Drop "significantly worse" and "hurt", or flag it as one nominal p among six.

##### F10. "Fixing GATA recovers half the gap" uses the observed-rank oracle, and the fig29 title disagrees with its bars

`A/report7_tf_evidence_multiome.qmd:107, 111, 118`
```markdown
## The transferred model under-ranks GATA elements, and fixing only those would recover half the gap
GATA-specific gain (GATA minus matched): +0.0103 [-0.0001, +0.0213] with in-cell ranks, +0.0217 [+0.0060, +0.0385] with observed ranks.
```
**Why it matters here**: The matched in-cell-rank oracle gives about a quarter of the gap, with an interval touching zero. In `A/figures/fig29_error_anatomy_gata.png` the plotted GATA bar reaches the full in-cell line; the "half" is GATA minus matched, which is not drawn. This oracle motivates all of Part A. (stats-causal and figure-visual both flagged this.)
**Fix**: Quote the in-cell-rank oracle (about 25%) as the share of the gap, show the observed-rank oracle as an upper bound, and plot or annotate the GATA-minus-matched difference.

##### F11. Mechanistic conclusions from correlational or single-run comparisons

`A/report7_tf_evidence_multiome.qmd:223, 545`
```markdown
Its correction cannot see DNA, and it gains no more, so that shortcut does not explain the small gain; the limit is the information the track itself carries.
Agreement across cells shows a shared mechanism, not a shared CRISPR gain;
```
**Why it matters here**: The two-stage comparison is single-run. Cross-cell agreement is consistent with a shared confounder (unannotated promoters, F20). Report 3:548 "better models lose more" was never tested as a difference between drops.
**Fix**: Use "consistent with" wording, and test the difference in drops before attributing it to anything.

##### F12. Ratios of point estimates quoted as precise ("a seventh", "57%", "83%")

`A/report7_tf_evidence_multiome.qmd:558`
```markdown
**A:** The far range: flagging only the nearest decile keeps 23% of the gain (+0.005 of +0.023), and with gene density in the model CGI distance still adds +0.013 [+0.005, +0.022], 57% of its gain.
```
**Why it matters here**: The numerator interval alone puts 57% anywhere from about 22% to 96%, yet the pre-set 50% bar was applied to the point estimate.
**Fix**: Bootstrap each ratio within the same paired draws, or print the numerator interval next to each percentage.

### Data pipeline & leakage

#### Major

##### F13. Every upstream number is scored on the chromosomes that chose the checkpoint; `test` is never read

`A/scripts/2.15.perfold_from_config.py:238`
```python
for fold in FOLDS_TO_SCORE:
    els = load_peaks(a.elements, folds_json[str(fold)]["val"])
# repo_scripts/multimodal_bpnet.py:542  if valid_loss < best_loss: torch.save(...)
# hg38_five_folds.json fold 0: "test": ["chr1","chr3","chr6"], "val": ["chr8","chr20"]
```
**Why it matters here**: All chromosome-holdout Pearsons, top-quintile deltas and paired p-values in reports 1-5 and the H3K27ac report are validation-set numbers. Arms whose loss changes the stopping rule are selected and compared on the same set. Only 10 of 24 chromosomes are ever evaluated. Verified: neither 2.15 nor 2.2 reads `"test"`.
**Fix**: Score on `folds_json[str(fold)]["test"]` in 2.15 and 2.2, relabel the existing tables as validation numbers, and check that the headline paired deltas survive.

##### F14. K562, the "unseen" cell, informed the design of the transfer fixes

`A/report6_improving_transfer.qmd:88`
```markdown
The "multi-cell model" is trained on four cell types (A549, GM12878, HepG2, MCF-7) and applied to K562, which it never saw.
# A/scripts/4.49.panel_tf_diversity.py:23  SPLIT-HALF RELIABILITY (added after the first run ... ranked noisy zinc-finger and GC-rich families above K562's GATA/TAL1).
```
**Why it matters here**: Four choices were made by looking at K562's own observed p300:
- the GATA error target;
- dropping CTCF from the deployable TF list;
- revising panel selection until it recovered GATA/TAL1;
- picking the tissue panel to cover the GATA error.

K562 benchmark gains from these choices are therefore not clean transfer estimates. The gains are small, so current conclusions hold, but the claim "never saw" overstates it. The leave-one-cell-out motif rule (report 7 Fig. 3) gives CTCF removal a legal justification.
**Fix**: State in reports 6 and 7 which choices used K562 information, and score design-derived changes on a second held-out cell using the other four leave-one-out models.

#### Minor

##### F15. ABC tracks take val chromosomes from the fold that early-stopped on them

`A/scripts/4.1.predict_h3k27ac_for_abc.py:107-111`
```python
# chromosome -> the lowest-numbered fold that did not train on it
for k in sorted(folds):
    for c in folds[k]["val"] + folds[k]["test"]:
        chrom_fold.setdefault(c, k)
```
**Why it matters here**: For in-cell K562 arms, checkpoint choice saw K562 labels on chr7, 8, 12, 17, 20, 21 and 22, slightly inflating the in-cell anchor and the 0.042 gap. Transferred arms are unaffected. 4.20 copies the rule.
**Fix**: Assign each chromosome from a fold where it is `test` first, falling back to `val` only for uncovered chromosomes (there are none).

##### F16. Upstream K562 scoring normalises ATAC over K562's own EP300 peaks

`A/scripts/2.42.make_multicell_k562_config.py:26-36`
```python
# K562 ATAC over THESE regions (the EP300 peaks), from 4.26. Not the ABC-candidate-region
TGT_MEAN, TGT_STD = 4.2275, 1.3051
```
**Why it matters here**: A new cell has no p300 peaks, so this statistic is not deployable. The policy is worth +0.0128 upstream, the same order as the report 5 effects. The downstream ABC path uses label-free statistics and is unaffected.
**Fix**: Compute the upstream statistics over label-free windows, or report both versions and say which is deployable.

##### F17. Report 7 summary bullets quote in-sample gains beside out-of-sample ones; fig40 does the same

`A/report7_tf_evidence_multiome.qmd:78-79`
```markdown
Element GC and CTCF sites, chosen from an error analysis, add +0.015 together (Fig. 11).
CpG-island distance alone gives +0.024 and absorbs the GC gain (Fig. 13).
```
**Why it matters here**: These sit next to the nested +0.021. In `A/figures/fig40_tf_evidence_headline.png` the "(pre-registered)" row at +0.036 reads as unbiased, and the selection optimism appears only in the caption. (data-pipeline and figure-visual both flagged this.)
**Fix**: Tag both bullets "(in-sample)", and add fig47c's nested-estimate marker to fig40b or label the row "(in-sample)".

##### F18. An open question proposes screening on the CRISPR-tested elements, which report 2 forbids

`A/report6_improving_transfer.qmd:319`
```markdown
Would an error screen restricted to the 3,406 CRISPR-tested elements agree better with the benchmark?
```
**Why it matters here**: Run as written, it puts benchmark labels into model choice, which report2:215-217 rules out.
**Fix**: Reword it as a diagnostic only, or restrict it to a held-out CRISPR split.

### Bioinformatics

#### Major

##### F19. "Observed p300 beats observed H3K27ac" compares libraries of very different depth

`A/report3_design_decisions.qmd:139, 1179`
```markdown
| **Observed p300 instead of observed H3K27ac as the activity term** | **+0.038 [+0.019, +0.058]** | use p300 where measured |
p300 is now the project's default target (F-024).
```
**Why it matters here**: The K562 H3K27ac libraries are 36 bp single-end, 2011-era, with about 25 reads per 500 bp element; the p300 library is several-fold deeper (findings file: 17.8x per 1 kb window). Poisson noise in the H3K27ac ranks could produce the gap by itself. This result re-targeted the project, and its caveats never mention depth.
**Fix**: Binomially thin the observed p300 BAMs to the pooled H3K27ac read count (0.41 already exists), rerun the observed-p300 ABC arm, and pair it against observed H3K27ac.

##### F20. "Non-promoter" means not within ±250 bp of one collapsed TSS per gene

`A2/scripts/8.8.cgi_distance_activity_panel.py:62-63, 123`
```python
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
U["prom"] = overlaps(U, tss)
```
**Why it matters here**: Alternative-transcript and lncRNA promoters count as distal elements. Those promoters would produce exactly the "near-CGI, accessible, low p300" signature that Part C reads as a CpG-island effect shared across five cells. The same issue undercuts the claim that GC is not a promoter proxy. The caveat appears only in a fold-out (L545); the summary states a shared mechanism.
**Fix**: Rerun Fig. 14 and the GC check excluding all GENCODE v43 transcript TSSs, lncRNA included, ±500 bp (optionally also K562 H3K4me3/CAGE peaks), and report how much of the near-decile depletion remains.

#### Minor

##### F21. The H3K27ac meta-profile sums strand 5' ends without a fragment shift

`A/scripts/0.5.make_5prime_tracks.sh:61`
```bash
bedtools genomecov -ibam "$bam" -5 -dz -strand "$strand" \
```
**Why it matters here**: Each nucleosome shows up as two strand peaks offset by about ±d/2, which fills the central dip and smears phasing. That undermines the report 1 statements "flanking nucleosome at ±250", "centre carries ~94%" and "broad envelope". The counting-window choice still stands.
**Fix**: Plot per strand, or shift 5' ends by ±d/2 before summing, then restate the positions.

##### F22. DNase 1 bp reproducibility and the converter's sequence gain are likely cleavage bias

`A/report4_atac_to_dnase_converter.qmd:90-91, 142`
```markdown
## Sequence adds +0.229 over an ATAC-only model, so the converter does more than rescale ATAC
```
**Why it matters here**: Base-level DNase I cleavage is dominated by shared sequence preference, and no bias model is used. Bias is the likely reason the converted track was a worse input (F-014).
**Fix**: Name enzyme bias as the leading explanation in reports 1 and 4, and use a bias-factorised target if the converter is ever reopened.

##### F23. The multi-cell panel uses the coverage ATAC input that the 2026-08-30 decision rejected

`A/scripts/1.33.make_panel_json.py:35, 48`
```python
"accessibility_bw": f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
"accessibility_bw": f"{PANEL_DIR}/{safe}_atac.bw",
```
**Why it matters here**: Coverage scales with read length, which per-cell z-scoring cannot remove, and the panel cells' read lengths are not recorded. F-042 bounds the cost (flat5p minus coverage +0.0041 [-0.0028, +0.0104]), but every report 5-7 number uses this input, and the reports never say so.
**Fix**: Record each panel cell's read length, and state the input convention and the F-042 bound in report 5 Methods.

##### F24. The "IKZF2 motif" is a GGAA-core 6-mer, and its definition changed between selection and the eQTL test

`A2/scripts/9.0.eqtl_transfer_prep.py:60-61`
```python
IKZF2 = "IKZF2_MA2326.1"
IKZF2_KMERS = ["AGGAAG", "CTTCCT"]
```
**Why it matters here**: A 6-bp GGAA core cannot be attributed to IKZF2 rather than ETS factors. Selection used BINDetect counts, while the transfer test used exact k-mers, and no correlation between the two versions is reported (GATA's is 0.957).
**Fix**: Call it a GGAA/ETS-like 6-mer, report the k-mer vs BINDetect correlation, and note the substitution in the Fig. 16 legend.

### LLM coding antipatterns

#### Minor

##### F25. 8.3 silences every warning, including ConvergenceWarnings from the selection fits (verify)

`A2/scripts/8.3.nested_selection.py:42, 171`
```python
warnings.filterwarnings("ignore")
fit = lambda C, i: lr(X[i], y[i], penalty="l1", solver="liblinear", C=C, max_iter=5000)  # noqa: E731
```
**Why it matters here**: Forward selection, L1 over about 1,700 columns and the 100-run null all depend on these fits. A non-converging fit would change the picks without showing it, and the out-of-sample +0.021 and the PASS rest on this script.
**Fix**: Scope suppression with `warnings.catch_warnings()` around known-benign calls, or count ConvergenceWarnings per fold and write the count out.

##### F26. Report 7's footprint range leaves out the two largest effects

`A/report7_tf_evidence_multiome.qmd:77`
```markdown
... yet change CRISPR AUPRC by -0.004 to +0.001 as features or inside the activity term (Figs 1b, 9).
```
**Why it matters here**: The weight-4 activity arms are -0.0098 and -0.0114, and fig40b plots them, so the bullet reads as a clean zero. (llm-failure-modes and figure-visual both flagged this.)
**Fix**: Write "-0.011 to +0.001".

##### F27. "A third" and "two thirds" of the GATA error use different baselines

`A/report7_tf_evidence_multiome.qmd:170`
```markdown
## Given the true TF list, the model removes a third of the GATA error; the cross-cell 'pair' loss removes two thirds
```
**Why it matters here**: +5.65 to +3.15 is a 44% cut, and 'pair' (to +2.14) is 62% measured from zero. Against the in-cell value of -2.6 the cuts are 31% and 43%. As written, 'pair' looks twice as good. F-035 repeats "a third".
**Fix**: Use one baseline, e.g. "cuts the GATA error from +5.7 to +3.1 and +2.1 points", and fix F-035 to match.

##### F28. Report 3's summary still says ATAC-input H3K27ac does not clear the floor

`A/report3_design_decisions.qmd:95`
```markdown
Predicted H3K27ac from the same architecture reading ATAC does not clear the floor (+0.025 [-0.003, +0.052]); reading DNase, it does, by +0.0403
```
**Why it matters here**: F-019 and the report's own L574 note show ATAC-input retrains clearing the floor (+0.0208 [+0.0039, +0.0373]).
**Fix**: Add "(an ATAC-input retrain clears it by +0.021, cause unresolved; F-019/F-020)".

##### F29. 9.1 replaces non-finite model features with 0 without counting them

`A2/scripts/9.1.eqtl_transfer_cluster.py:63`
```python
P = P.replace([np.inf, -np.inf], np.nan).fillna({c: 0 for c in COLS[7:]})
```
**Why it matters here**: Under log(|x| + 0.01) a filled 0 becomes a real low value. A cluster with many gaps would shift its scores without warning. The same line is in 10.1.
**Fix**: Log non-finite counts per feature per cluster, and assert they stay below a small fraction.

### Documentation & schema fidelity

#### Major

##### F30. The decision log says RC averaging is opt-in; the code and report 3 make it the default

`.living/decisions.md:267-271`
```markdown
### [2026-09-03] Test-time reverse-complement averaging stays opt-in
**Decision**: Ship it behind `--rc-average`, off by default.
# A/scripts/2.15.perfold_from_config.py:33  TEST-TIME REVERSE-COMPLEMENT AVERAGING IS ON BY DEFAULT (since 2026-09-03).
```
**Why it matters here**: RC averaging shifts every upstream number by +0.002 to +0.016, and the decision record is what says which tables are comparable. Report 3:379 tells readers to use `--rc-average`, a flag 2.15 no longer has.
**Fix**: Add a dated decision for the switch to default-on, mark the 09-03 entry superseded, and fix report 3:379.

##### F31. Report 5 lacks the required Headline figure, Goals, Open questions and Regenerate sections

`A/report5_multicell_p300_transfer.qmd:62, 263-298`
```text
ERROR: Required section missing: '## Goals'
ERROR: Required section missing: '## Headline figure' — the report must lead with a headline figure
warning: No '### Regenerate' subsection inside '## Methods'.
```
**Why it matters here**: Report 5 fails the lab lint, so it was rendered with lint off. Its figures (`scripts/5.1`) cannot be rebuilt from the document, and it is the only report without an environment or lock line. (doc-fidelity and reproducibility-env both flagged this.)
**Fix**: Add the four sections (headline = fig22; Regenerate = `pixi run -e multimodal python scripts/5.1...` plus render), then render with lint on.

##### F32. Report 3 headings still state claims that the callouts under them withdraw

`A/report3_design_decisions.qmd:1005, 468, 572`
```markdown
## DNase beats ATAC as the accessibility input, and it is the first H3K27ac model to clear the floor on transfer
> **Superseded (F-016, 2026-09-14):** a DNase-input model must be compared with a DNase-only floor; ...
```
**Why it matters here**: Headings are the findings in this report style and appear in the TOC. A skim reader gets the withdrawn claim, which also contradicts the report's own summary.
**Fix**: Retitle L1005, L572 and L468 to their corrected readings.

#### Minor

##### F33. Report 6 says no change moves ABC by more than 0.011, but its headline figure shows oracle rows at +0.019 and +0.041

`A/report6_improving_transfer.qmd:2, 66, 76`
```markdown
**Figure 1 | Several changes raise correlation with p300 by 0.1 or more, but none moves ABC by more than 0.011 AUPRC ...
... the ranking loss at weight 10 lowers it by 0.015 and gives the largest single-run ABC gain (+0.0087)
```
**Why it matters here**: The report counts the GATA oracle as one of its seven changes. "Largest single-run gain" is also false against its own figure: headpair is +0.0105 and the ATAC list +0.0100. (doc-fidelity, llm-failure-modes and figure-visual all flagged this.)
**Fix**: Write "no deployable change ... (oracle upper bounds reach +0.041)", and drop or restrict "largest".

##### F34. Report 2 still recommends the partial-correlation screen that F-031 found failing

`A/report2_evaluation_methodology.qmd:92, 172-174, 263`
```markdown
prediction knows about p300 that ATAC does not already supply. It is a candidate screen: it
orders models with large differences correctly, and it missed the one small gain it has been
```
**Why it matters here**: F-031 shows +0.06 to +0.22 partial rho with no AUPRC change, and report 6 shows +0.125 with no ABC change. A reader following report 2 would rank models with a screen that has now failed in both directions. (llm-failure-modes and doc-fidelity both flagged this.)
**Fix**: Add an Update (F-031) callout and limit the screen to exploration.

##### F35. Report 7 proposes the eQTL readout that F-066 has since built and failed

`A/report7_tf_evidence_multiome.qmd:639, 653`
```markdown
an eQTL readout restricted to expressed genes, where the activity term decides the ranking, is the alternative to build and calibrate with the ABC arm.
**Findings.** F-032 to F-041 (Part A), F-047 to F-056 (Part B) and F-057 to F-065 (Part C)
```
**Why it matters here**: A reader would queue work that has already been run and failed. (llm-failure-modes and doc-fidelity both flagged this.)
**Fix**: Replace the question with F-066's result, and extend the findings range to F-066.

##### F36. The F-063 and F-064 titles, Status lines and ledger still state the withdrawn readings

`.living/findings/predicting-regulatory-element-function-at-scale.md:2518, 2560, 2593`
```markdown
## F-064: The eQTL benchmark resolves a CRISPR-scale change: ... so F-063's null rules out the gain F-059 predicts
**Status:** supported (positive control for F-063; one comparison; RNA differs between the arms)
```
**Why it matters here**: Readers, the dashboard and the INDEX scan titles and Status lines. Other withdrawals in this file (F-010, F-019, F-025) carry the withdrawal in the title or Status line. (llm-failure-modes and doc-fidelity both flagged this.)
**Fix**: Retitle both findings to "inconclusive, see F-065/F-066" and update the Status and ledger cells.

##### F37. ANALYSIS_MANIFEST points at the superseded report and repeats revised findings

`analysis/ANALYSIS_MANIFEST.md:231-245`
```yaml
last_updated: 2026-08-24
  - Sequence adds ~2.4x less for H3K27ac than for p300 (+0.127 vs +0.304 ...). See F-001.
report: 2026_0824_H3K27ac_model/h3k27ac_model_report.html
```
**Why it matters here**: Reports 1-6 are registered nowhere, and the entry carries a deployment recommendation that F-016 undercuts.
**Fix**: List reports 1-7, update key_findings to their current status, and bump last_updated.

##### F38. Report 4's Regenerate runs `build_reports.py` on the superseded report

`A/report4_atac_to_dnase_converter.qmd:256`
```bash
python build_reports.py h3k27ac_model_report.qmd .
python render_report.py report4_atac_to_dnase_converter.qmd
```
**Why it matters here**: Reports 1-3 have diverged from the frozen source. If the script re-splits it into the current directory, it could overwrite them (not confirmed: the script was not in the mirror). Goals still give the converter's premise as "advantage surviving transfer", which F-016 withdrew.
**Fix**: Drop the line or say what it writes, and note the F-016 withdrawal in Goals.

##### F39. Reports 1, 2 and 4 have no Open questions section

`A/report1_data_characterisation.qmd:322`
```markdown
## Methods
# (no "Open questions / next steps" section; same in report2 and report4)
```
**Why it matters here**: The lab convention requires one, so open work can be traced into TODOLIST.
**Fix**: Add a short section to each (it can say "none, closed").

##### F40. No report uses the mandated arm label

`CLAUDE.md:307-308`
```markdown
- **Labels:** report every arm as `input -> target | trained on | scored on | ATAC-only at
  inference: y/n`.
```
**Why it matters here**: Deployable and oracle arms share forest plots and are told apart only in prose, which is the root of F33.
**Fix**: Apply the label in the arm tables and legends of reports 5-7 at least.

##### F41. INDEX.md and the findings front matter are a month stale

`.living/INDEX.md:3-9`
```markdown
Last audit: 2026-09-12
| findings/ | 12 findings across 4 topics | 2026-09-12 | ...
```
**Why it matters here**: CLAUDE.md says to trust INDEX, but tag lookup misses F-013 to F-066.
**Fix**: Rerun the INDEX generator and bump the findings `last_updated`.

### Code quality

#### Major

##### F42. "AUPRC" is three different estimators across the bootstraps, and two of them share one axis

`A/scripts/4.17.paired_auprc_bootstrap.py:63-70`
```python
return float((prec * yt).sum() / n_pos) if n_pos else np.nan          # 4.17 ABC: step AP, inner join
return statistic_aupr(np.asarray(y).astype(np.int64), np.asarray(s, dtype=float))   # 4.34 rE2G
p, r, _ = precision_recall_curve(yy, s); return auc(r, p)              # A2 4.2/8.0/8.3 scE2G
```
**Why it matters here**: fig39b draws ABC and rE2G AUPRCs as adjacent bars. fig16 and fig21/fig24 give absolute levels for the same arms on the same pairs that differ by about 0.011 (floor 0.457 vs 0.468), the size of the transfer effects. The output schemas also differ, so consumers rename columns by hand. (code-quality and figure-visual both flagged this.)
**Fix**: Use one shared bootstrap module with an explicit `estimator=` argument and one output schema, and name the estimator wherever absolute levels are plotted.

##### F43. Every rE2G AUPRC is computed by code imported from a live, unpinned checkout

`A/scripts/4.34.re2g_paired_bootstrap.py:29-31`
```python
sys.path.insert(0, "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/workflow/scripts/"
                   "model_training")
from training_functions import statistic_aupr  # noqa: E402
```
**Why it matters here**: A `git pull` there silently changes the metric behind about 15 rE2G tables, including the Corgi rE2G verdict.
**Fix**: Vendor `statistic_aupr` with a source-commit link, or pin and assert the checkout commit.

##### F44. fig39c compares ABC against one run and rE2G against the 3-seed mean in the same row

`A/scripts/5.4.plot_benchmark_round2.py:98-99`
```python
("pred. p300 x corrected\nH3K27ac - pred. p300", ("k27frag_pred_p300k27.ABC.Score", "p300pred_multicell_holdK562.ABC.Score"),
 ("ATAC_p300k27pred_powerlaw_all", "mse"))]
```
**Why it matters here**: Two points drawn as one contrast use different reference models, and the caption says panel c is single runs. Verified in the code.
**Fix**: Score both frameworks against the same reference, and keep the caption in step.

#### Minor

##### F45. fig32 parses per-seed AUPRCs from SLURM logs pinned to job IDs

`A/scripts/5.3.plot_tf_activity_report_figures.py:366-377`
```python
GA = groups(f"{P}/log/seeds_boot.46073940.txt")
GR = groups(f"{P}/log/re2g_seeds_boot.46073942.txt")
```
**Why it matters here**: A bootstrap rerun writes a new log that fig32 never reads, so fig32 can silently disagree with fig26.
**Fix**: Have 4.17/4.34 write the group-member AUPRCs to a `.groups.tsv`, and read that.

##### F46. Report 3 names 3.3 for a figure that 3.4 writes

`A/report3_design_decisions.qmd:233, 282`
```markdown
![](figures/fig5_residual_grid_both.png)
pixi run python scripts/3.3.plot_residual_figure.py
```
**Why it matters here**: Following report 3's own instructions does not rebuild its Figure 4, and the two near-duplicate scripts will drift.
**Fix**: Delete 3.3, and point both reports at 3.4.

##### F47. The numbers-manifest lint covers only reports 1-3

`A/scripts/3.7.build_numbers_manifest.py:59-454`
```python
d = tsv("window_tradeoff.tsv") ... d = tsv("coupling_panel_recomputed.tsv") ...
# (no *_auprc_pairs.tsv, p300_multicell_k562_per_fold.tsv, or A2/results/* entries)
```
**Why it matters here**: Report 4-7 numbers are hand-typed and never linted. The report 5 five-intervention table exists only in 5.1's stdout. The manifest caught three wrong numbers in reports 1-3 (learnings 2026-09-01).
**Fix**: Have 5.1 write its table to a TSV, and extend 3.7 to the bootstrap and A2 tables.

##### F48. Result values are hand-typed into plotting and analysis code

`A/scripts/4.18.plot_p300_benchmark.py:46-54`
```python
PAIRED = {
    "p300obs_k562":                ("+0.099", "[+0.079, +0.120]", ""),
    "p300pred_gm12878_multimodal": ("+0.009", "[-0.004, +0.022]", "n.s."),
```
**Why it matters here**: The same pattern appears in 5.6 (0.627 base, 0.021 nested line), in 10.2 (the 0.021 / 0.050 calibration scale), and in 5.1 (peak counts). A rerun moves none of these values. (llm-failure-modes and code-quality both flagged this.)
**Fix**: Read the values from the result tables at plot time, or register them as manifest constants.

##### F49. Shared helpers are reached through exec and argv hacks

`A2/scripts/10.3.eqtl_withingene_diagnostics.py:16-20`
```python
src = open("scripts/10.2.eqtl_withingene_gather.py").read()
for name in ["ap_steps", "weighted_ap", "unit_ap"]:
    exec(re.search(rf"^def {name}\(.*?(?=^\S)", src, re.S | re.M).group(0), ns)
```
**Why it matters here**: 10.3 breaks outside A2/ and would silently pick up a renamed function. Importing 9.0 runs argparse and `os.makedirs` at import time.
**Fix**: Move the replica, prep and AP helpers into an importable `sce2g_lib.py`.

##### F50. The paired within-fold t interval is reimplemented in about 15 scripts

`A/scripts/3.10.plot_rc_averaging.py:22, 41-45`
```python
TCRIT = tdist.ppf(0.975, df=4)
a = rc[rc["config"] == cfg].sort_values("fold")[metric].to_numpy(float)
d = a - b
```
**Why it matters here**: Arrays are paired by position with df hard-coded to 4. Every current table has 5 folds, so no number is wrong today, but a missing fold would mis-pair silently.
**Fix**: Add one `paired_fold_delta()` helper that pivots on fold and asserts matching fold sets.

##### F51. Training regimes are pairwise-excluded by hand in four places

`repo_scripts/train_multimodal_bpnet.py:823-841`
```python
if args.panel_mean:
    bad = [f for f, on in (("--count-offset-model", args.count_offset_model),
```
**Why it matters here**: main() is about 665 lines, and every new regime needs every flag pair re-checked; a missed combination trains silently.
**Fix**: Replace the flags with one `--regime` choice and per-regime data-assembly functions.

##### F52. The invalid within-gene permutation null is still runnable

`A2/scripts/8.3.nested_selection.py:53`
```python
ap.add_argument("--permute-seed", type=int, default=None, help="INVALID null, kept for the record")
```
**Why it matters here**: It makes it easy to rerun a null recorded as invalid on 2026-10-06; git already keeps the record.
**Fix**: Remove the branch.

### Reproducibility & environment

#### Major

##### F53. Regenerate blocks in four reports run the wrong pixi env

`A/report2_evaluation_methodology.qmd:394-408`
```bash
Environment: pixi env `multimodal`, pixi.lock sha256:069960d36312
pixi run python scripts/2.15.perfold_from_config.py config/<name>_configs.json \
```
**Why it matters here**: Plain `pixi run` uses the default env, which has no bpnetlite, so 2.15, 2.2, 2.4 and 4.1 fail on import. The same pattern is in h3k27ac_model_report, report 1 and report 3. Reports 6 and 7 already use `-e multimodal`.
**Fix**: Use `pixi run -e multimodal` throughout, or make `multimodal` the default environment.

##### F54. Report 7's figure script rebuilds a committed table from $SCRATCH on every run

`A/scripts/5.5.plot_tf_evidence_report_figures.py:168-183`
```python
if True:  # rebuilt every run from $SCRATCH; the committed copy is the record
        t = pd.read_csv(f"{RUNS}/{run}/K562_Xu_pl/{m}/model/model_coefficients.tsv", sep="\t")
pd.concat(rec).to_csv(cf, sep="\t", index=False)
```
**Why it matters here**: After the purge (about January 2027 for these runs), report 7's regenerate command crashes before fig43 is written. Until then, each run rewrites the source that the fig43b numbers cite. (reproducibility-env and code-quality both flagged this.)
**Fix**: Read the committed TSV by default, and rebuild only behind an explicit `--rebuild-coefficients` flag.

#### Minor

##### F55. Inputs the reports depend on live only on purgeable $SCRATCH

`A2/scripts/4.0.setup_featB.sh:20-29`
```bash
BASE=$SCRATCH/sce2g_runs/2026_1002_k562_xu
EXT=$SCRATCH/sce2g_external
# A/scripts/0.9.make_atac_fragment_channels.sh:31  BAM_DIR="${SCRATCH}/atac_pe"
```
**Why it matters here**: These live only on $SCRATCH: the baseline scE2G run and every Part B/C run, the pairs files, the feature library, the eQTL outputs (A2 scripts read `os.environ["SCRATCH"]` with no fallback), and the K562 PE ATAC BAMs (mtime 2026-08-25, purge-eligible about 2026-11-23). None has a DATA_MANIFEST entry.
**Fix**: Before about 2026-11-20, copy the baseline scE2G run and the three `.pe.bam` files to Oak or Elm (or record how to re-fetch them), and add DATA_MANIFEST rows with a purge note.

##### F56. Launchers that the reports name sit in gitignored `log/`

`A/report7_tf_evidence_multiome.qmd:484-486`
```bash
bash log/predp300_launch.sh            # F-057: 4.10 -> 4.8 -> driver -> 4.2
# .gitignore:85  **/log/
```
**Why it matters here**: The chains behind F-057 to F-059, eqtl_transfer_launch.sh, and the CLAUDE.md template launcher are not in git history.
**Fix**: Move the launchers into `scripts/` (or `git add -f` them), and update the report and CLAUDE.md paths.

##### F57. Several commands use bare `python` or `Rscript`

`A/report7_tf_evidence_multiome.qmd:372, 404`
```bash
python scripts/4.9.make_activity_pairs.py --out $SCRATCH/sce2g_external/activity_pairs.tsv.gz
Rscript scripts/5.0.export_atac_matrix.R ...
```
**Why it matters here**: Whatever interpreter is on PATH runs these, and neither pixi.toml contains R. Report 4:252-256 has the same issue.
**Fix**: Prefix each command with the env actually used, and add r-base and r-matrix to an A2 feature.

##### F58. The scE2G commit cited does not match the checkout that ran

`A/report7_tf_evidence_multiome.qmd:649`
```markdown
scE2G main (7cb2af7) with EngreitzLab/scE2G#118 applied locally, training workflow, `--rerun-triggers mtime`
# cluster: branch local-with-pr118, HEAD 729c648 on 4f57df9
```
**Why it matters here**: Every Part B/C AUPRC came from local, unpushed code in `scE2G_temp`, so nobody else can rebuild it exactly.
**Fix**: Record 729c648, and push the branch or commit a `git format-patch` copy into A2.

##### F59. The SLURM srun patch is unversioned, and only some drivers check for it

`A/scripts/4.33.submit_re2g_training.sh:32`
```bash
grep -q 'SLURM_CPUS_PER_TASK:-1' ... || { echo "ERROR: the srun patch is missing ..."
# no guard: A/scripts/4.5.submit_abc.sh, A/scripts/4.7.submit_crispr_comparison.sh, A2 log/sce2g_train.sh
```
**Why it matters here**: An env rebuild reverts the patch, and the ABC, CRISPR and scE2G drivers would then fail in a way that looks like a workflow bug.
**Fix**: Commit the diff with an apply-or-verify helper, and call the guard from 4.5, 4.7 and sce2g_train.sh.

##### F60. The regression-gate command points to a nonexistent path

`A/report2_evaluation_methodology.qmd:405`
```bash
pixi run python scripts/test_asymmetric_loss.py
# the file is at the repo root: scripts/test_asymmetric_loss.py
```
**Why it matters here**: Run from A/ as written, the gate that must pass before scoring fails with "file not found". Report 3:1164 repeats it.
**Fix**: Use `pixi run -e multimodal python ../scripts/test_asymmetric_loss.py`.

##### F61. fig25 and fig33 are PNG only

`A/scripts/5.2.plot_partial_corr_explainer.py:125-126`
```python
out = f"{OUT}/fig25_partial_corr_explainer.png"
fig.savefig(out); print("wrote", out)
```
**Why it matters here**: They are the only report figures without a vector copy.
**Fix**: Save through `nature_style.save_fig` (PNG + PDF) in 5.2 and 4.49.

##### F62. CLAUDE.md, ENVIRONMENTS_INSTALLATIONS.md and DATA_MANIFEST describe the old setup

`CLAUDE.md:421-426`
```markdown
Primary conda environments:
- `bpnet_37` - For BPNet training, prediction, and SHAP (requires `module load cuda/11.1.1 cudnn/8.1.1.33`)
```
**Why it matters here**: Someone setting up from these docs would build the wrong envs, and could not find the inputs for reports 4-7 (DATA_MANIFEST was last touched 2026-08-25 and has 13 entries).
**Fix**: Document the pixi envs, the A2 workspace, the Snakemake wrapper envs and patch, and add the datasets from report 6's table.

### Vault fidelity

#### Minor

##### F63. The analysis has no EngreitzLabVault pointer note

`EngreitzLabVault/Atlas/Projects/ (missing note)`
```markdown
# Atlas/Projects/ holds 53 notes; none for EP300_BPNet. Only mention, Atlas/Concepts/p300-collab.md:
"preliminary, not-yet-published lab science tied to [[Maya Sheth]]'s p300 BPNet project"
```
**Why it matters here**: Lab members searching the Vault for p300, H3K27ac or ABC-activity modelling won't find this work or its eight reports. The repo has no `vault-registration-skips.md` to explain why.
**Fix**: Register a lean project note (question, Oak path, GitHub repo, links to p300-collab, chrombpnet, ABC, scE2G; status active, maturity shaping) that references the eight reports by path, and link it from p300-collab.

### Figure integrity & design

#### Major

##### F64. Report 2's headline legend quotes 0.831 for a bar that is 0.814, and says "five models" for three

`A/report2_evaluation_methodology.qmd:73, 82`
```markdown
**Figure 1 | The same five models give different effect sizes ...**
falls from 0.494 to 0.380 and ATAC-only from 0.831 to 0.548, so the gap between modalities is
# results/fiveprime_stratified_fold_summary.tsv: atac5p atac 500 all 5 0.8136 ... 0.7941 0.8332
```
**Why it matters here**: 0.831 is close to the CI upper bound, and the derived "gap 0.34" is really 0.32. This is the headline of the scoring-rules report. Verified against the table.
**Fix**: Use 0.814 (gap 0.32), and say "three models".

##### F65. The fig15 caption contradicts panel c and gives the wrong tail size for panels b-d

`A/figures/fig15_error_strata.png` (caption `A/report3_design_decisions.qmd:690-693`)
```markdown
The 1% tails are n = 1,536 each ... peaks report enrichment for every mark, signal for CTCF only.
# panel c: EP300 signal 1.37x, H3K27me3 peaks 0.84x; panels b-d labelled "over-predicted 5%", n = 7,678
```
**Why it matters here**: Panel c exists to separate peak overlap from signal, and EP300 signal enriched at over-predicted elements bears on the CTCF reading.
**Fix**: Rewrite the panel c sentence to match the bars, and say panels b-d use 5% tails.

##### F66. fig39b draws AUPRC bars from a 0.45 baseline with one framework's floor and no intervals

`A/figures/fig39_benchmark_round2.png` (`A/scripts/5.4.plot_benchmark_round2.py:78-79`)
```text
y-axis 0.45-0.72, no break: ABC ATAC x H3K27ac (0.513) is ~0.06 tall, rE2G H3K27ac (0.620) ~0.17 tall.
Only reference line: "ABC ATAC-only floor" (0.468); the rE2G floor (0.544) is not drawn.
```
**Why it matters here**: A 1.2x difference looks like 3x, and the panel carries "observed marks do best without ATAC in ABC".
**Fix**: Use dots with paired intervals (or bars from 0), and draw each framework's own floor.

##### F67. fig41 draws AUROC bars from 0.3 when chance is 0.5

`A/figures/fig41_footprint_chip_validation.png`
```text
y-axis 0.3-0.9, dashed line at 0.5: homemade GATA1 (0.419, below chance) is drawn as a positive bar.
```
**Why it matters here**: The section's claim (homemade below chance, TOBIAS/PRINT above) is exactly the information the bar baseline hides.
**Fix**: Anchor the bars at 0.5, or use a dot plot with the chance line.

#### Minor

##### F68. fig5 has unexplained dashed lines, and its control is not "flat at zero"

`A/figures/fig5_residual_grid_both.png` (`A/scripts/3.4.plot_residual_grid_both.py:90`)
```python
ls=(0, (3, 1.5)) if hatched else "-")
# caption: "the ATAC control stays flat at zero across all ten strata" (GM12878 control -0.06 to +0.075)
```
**Why it matters here**: Readers cannot tell the objectives apart, and the negative-control claim is stronger than the picture.
**Fix**: Add "dashed = residual objective" to the legend, and say "within ±0.08 of zero".

##### F69. profile_comparison colours and panel letters do not match the captions

`A/figures/profile_comparison.png` (`A/scripts/0.6.compare_profiles.py:47`)
```python
COLORS = ["#2166AC", "#762A83", "#B2182B", ...]
# captions: "Teal, K562 ATAC; orange, K562 H3K27ac 5' ends"; panels a/b referenced, none drawn
```
**Why it matters here**: The caption describes a different image, and purple means "sequence + ATAC" under the report-wide convention.
**Fix**: Recolour to teal and orange as captioned, and add panel letters.

##### F70. Colours change meaning within and between reports

`A/figures/fig14_crispr_benchmark.png`, `fig25`, `fig19-24`, `fig40`, `fig46`, `fig48`, `fig50`
```text
report2:80 "red sequence only, blue ATAC only, purple sequence + ATAC" - fig14 red = GM12878 multimodal.
Multi-cell orange in fig25, green in fig19-24. Ochre = observed p300 ceiling in fig40/45/46, ABC control in fig50.
fig48: red/dark blue = H3K27ac/p300 in panel a, K562/GM12878 in panel b.
```
**Why it matters here**: Every figure has a legend, so none is wrong on its own. But colour meanings flip as a reader moves through the reports, and in fig50 the ceiling colour is reused for a control.
**Fix**: Fix one colour key per report (modality, arm class, mark, cell) and reuse it.

##### F71. The fragment-length figure does not show the bins report 3 cites it for

`A/figures/atac_fragment_lengths.png` (`A/scripts/0.10.plot_fragment_lengths.py:55`)
```python
return v[(v > 0) & (v <= max_len)]
# report3:338 "Bin edges sit in the troughs ... (Report 1, Fig. 5)"; the 620 bp edge is past the tri-nucleosome maximum
```
**Why it matters here**: The cited figure shades the superseded windows, and the "50.0% in neither" figure is really 50% of fragments up to 800 bp.
**Fix**: Shade the five current bins, justify or move the 620 bp edge, and note the 800 bp cap.

##### F72. fig4 and fig1 give different H3K27ac values for what read as the same models

`A/figures/fig4_p300_vs_h3k27ac.png`
```text
fig4a: ATAC-only 0.542, seq+ATAC 0.669 (+0.127); fig1b: 0.548, 0.685 (+0.137); same footer.
```
**Why it matters here**: Summaries quote +0.137 and +0.127 side by side as if they came from one model set.
**Fix**: Name the model generation in the fig4 caption, or rerun fig4 with the 5' models.

##### F73. fig18 caption ratios do not match the bars, and panel b mixes count sources

`A/figures/fig18_multicell_panel.png`
```text
GM12878/HepG2 = 2.58x (caption 2.4x); A549/MCF-7 = 6.6x (caption 4.6x);
panel_depth_audit.tsv: GM12878, K562 = tagAlign; A549, HepG2, MCF-7 = bam.
```
**Why it matters here**: The caption cannot be checked against the bars, and the depth bars count different things.
**Fix**: Correct the ratios, and count all cells the same way or label the source on each bar.

##### F74. fig11 plots DNase series that no caption mentions

`A/figures/fig11_profile_ceiling.png`
```text
Legend: "K562 DNase" and "GM12878 DNase" lines in both panels, stratum unstated.
```
**Why it matters here**: The two highest curves in the figure are unexplained.
**Fix**: Describe the DNase series and their stratum, and add "top quintile" to the solid-line legend entries.

##### F75. fig8 lacks a colourbar and a key for its grey shades

`A/figures/fig8_coupling_across_celltypes.png`
```text
Hexbin with no colourbar; three grey shades in panels d/e with no key; "−VEGF" unexplained.
```
**Why it matters here**: This is report 1's headline figure, and its shading cannot be decoded.
**Fix**: Add a colourbar, and explain the shades and "−VEGF" in the caption.

##### F76. The fig9 caption says p300 trends "the same way" when panel b trends the opposite way

`A/figures/fig9_atac_5prime_input.png` (`A/h3k27ac_model_report.qmd:242`)
```markdown
p300 is underpowered ... with point estimates trending the same way.
```
**Why it matters here**: In panel b, p300 ATAC-only is about -0.012 and multimodal +0.017, the reverse of the H3K27ac pattern. Report 3's copy of this caption already dropped the phrase.
**Fix**: Restrict the sentence to panel a, as report 3 does.

##### F77. Clipped labels, colliding annotations and cut-off whiskers in about ten figures

`A/figures/fig17_converter.png`, `fig12`, `fig5`, `fig19`, `fig15`, `fig9`, `fig14`, `fig22`, `fig36`, `fig47`, `fig39`
```text
fig17c x-label cut to "gain over the observed ATAC tracl"; fig36d legend sits on a data point;
fig47c "nested prediction" overlaps tick labels; fig14 and fig22 whiskers clipped at the axis limit.
```
**Why it matters here**: The bottom row of fig36d and the clipped intervals cannot be read.
**Fix**: Widen the margins, move legends outside the axes, and set limits that contain every interval.

##### F78. Reports 5-7 figures bypass nature_style and break the stated design rules

`A/scripts/5.1.plot_multicell_report_figures.py:29-57, 138`
```python
plt.rcParams.update({ "font.family": "sans-serif", ... "figure.dpi": 160, ...})
C = {"floor": "#bdbdbd", ... "multicell": "#238b45"}
ax.grid(axis="x", color="#eee", zorder=0)
```
**Why it matters here**: These figures break the stated design rules:
- gridlines in fig19;
- grey axes and ticks in fig25;
- grey footers and coloured text in fig26, fig31 and fig47;
- labels that are not sentence case across fig15-38;
- ColorBrewer and ad hoc palettes in 5.1, 5.3, 5.4, 4.56 and 4.59;
- DejaVu at 9 pt instead of Arial at 7 pt.

(code-quality, reproducibility-env and both figure reviewers flagged this.)
**Fix**: Route 5.x, 4.56 and 4.59 through `nature_style`, remove gridlines, make text black, apply sentence case, and map colours to the Nature palette.

##### F79. fig28 and fig26 give two values for one neighbour arm

`A/figures/fig28_rank_and_neighbour.png`
```text
fig28b +0.106 (baseline 0.291, crossclass refit) vs fig26a +0.124 (baseline 0.268, ATAC rank + square).
```
**Why it matters here**: Readers comparing the two figures see two gains and two baselines for one variant.
**Fix**: Use one partial-correlation definition, or state the covariates and baseline in each caption.

##### F80. fig48 says "percentile units" for a 0-1 rank fraction

`A/figures/fig48_cgi_activity_panel.png`
```text
y-label "(percentile units)", ticks -0.15 to 0.05; fig31 and fig45b use a 0-100 "percentile points" scale.
```
**Why it matters here**: The K562 deficit is 12.5 points; read on the other figures' scale, it looks 100x smaller.
**Fix**: Multiply by 100 and say "percentile points", or label it "rank fraction (0-1)".

##### F81. Missing axis titles, colour keys and readable feature names

`A/figures/fig35_nucleosome_profiles_k562.png`, `fig34`, `fig38`, `fig45`, `fig47`
```text
fig35 b/c no y title; fig34a colourbar unlabelled; fig45a three colours with no key and an undescribed "Xu ATAC" row;
fig47b raw IDs as ticks ("screen_composite", "Wt1_MA1627.2").
```
**Why it matters here**: Units and colour meanings cannot be read from the panels.
**Fix**: Add the titles, keys and readable names, and describe the Xu ATAC row.

##### F82. fig33's title promises four column blocks but draws two dividers

`A/figures/fig33_panel_tf_diversity.png`
```text
Title "K562 | current panel | greedy additions in order | rest"; dividers after K562 and MCF-7 only; no ATAC/DNase marker.
```
**Why it matters here**: The reader cannot see where the greedy picks end, or check the caption's assay explanation.
**Fix**: Add the third divider and an assay strip.

## What was checked but is fine

- **Statistics & causal inference**: the paired design is sound throughout. Every arm is scored on one pair set with shared resamples, seed groups average inside each draw, and the Part C nested CV with a stratified permutation null is well built.
- **Data pipeline & leakage**: rE2G is scored on LOCO predictions only, K562 is held out at the code level with incomplete panel entries refused, region sets are shared and guarded, and normalisation for downstream inference is label-free.
- **Bioinformatics**: the coordinate and assay plumbing is careful. Summit offsets, Tn5 shifts, binomial thinning, PE read-1 handling, N-masked motif hits, JASPAR symbol mapping, TOBIAS/PRINT insertion conventions and eQTL filters all check out.
- **LLM coding antipatterns**: no fabricated numbers. Headline numbers in every report match their result tables, and the revision history shows caveats kept and claims withdrawn rather than the reverse.
- **Documentation & schema fidelity**: every HTML matches its qmd (headings and every decimal), all image paths resolve, every cited script exists, and the superseded report carries a clear banner.
- **Code quality**: the plotting scripts fail fast (asserts on missing arms and rows), and the existing fold tables all have five folds, so F50 changes no current number.
- **Reproducibility & environment**: both pixi.lock files are committed and match the hashes the reports quote, bootstraps and subsampling are seeded, and big data is ignored and guarded by the 5 MB hook.
- **Vault fidelity**: skipped; there is no pointer note to check (F63).
- **Figure integrity & design**: no gridlines and black axes in reports 6-7, and nearly every plotted number checked (about 150) matches its caption and table.

## Finding tally

| Category | Major | Minor | Total |
| --- | ---: | ---: | ---: |
| Statistics & causal inference | 7 | 5 | 12 |
| Data pipeline & leakage | 2 | 4 | 6 |
| Bioinformatics | 2 | 4 | 6 |
| LLM coding antipatterns | 0 | 5 | 5 |
| Documentation & schema fidelity | 3 | 9 | 12 |
| Code quality | 3 | 8 | 11 |
| Reproducibility & environment | 2 | 8 | 10 |
| Vault fidelity | 0 | 1 | 1 |
| Figure integrity & design | 4 | 15 | 19 |
| **Total** | **23** | **59** | **82** |

## Notes

- **Compound remediations.** Some findings share a fix:
  - F1, F4, F5, F6, F8 and F33 all come down to how intervals are worded. One pass with a fixed vocabulary would clear most of them: "resolvable", "not resolved, interval includes X", "single run, provisional".
  - F7 and F42 share the bootstrap code: one shared module with a gene-level option and an explicit estimator.
  - F53, F55, F56 and F60 are one regenerate audit.
- **Corgi idea 1 (2026-10-13, not in these reports)** is exposed to F7 (its lower bounds are +0.0010 ABC and +0.0002 rE2G) and to F43 (its rE2G metric). Re-read it under a gene-level bootstrap before adopting it.
- **Dropped as low confidence.** Three findings were dropped, per the synthesis rules (low-confidence Minor):
  - a stated noise floor of 0.018 vs 0.0018 in report 2;
  - MCF-7 depth heterogeneity in the panel;
  - the unfiltered IgG control in 4.13.

  Each is worth a glance, not a fix.
