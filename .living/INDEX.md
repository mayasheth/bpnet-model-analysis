<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-09-08

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 17 sections | 2026-09-04 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 18 entries (large — read selectively) | 2026-09-06 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric, Paired-end H3K27ac targets use read 1 only, not both mates |
| learnings.md | 48 entries (large — read selectively) | 2026-09-06 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 9 findings across 4 topics | 2026-09-08 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-09-08 (heuristic)

## Tag clusters

- **h3k27ac** (22 entries) — D-1, D-2, D-3, D-4, D-5
- **5-prime** (11 entries) — L-18, L-28, L-29, D-5, D-6
- **atac** (9 entries) — L-29, L-33, D-4, D-6, D-8
- **ceiling** (9 entries) — L-14, L-17, L-18, L-24, D-2
- **architecture** (8 entries) — L-42, L-45, D-3, D-7, D-15
- **negative-result** (8 entries) — L-45, D-13, D-14, D-15, D-18

## Most recent (10)

- [2026-09-06] L-47: bpnet-gc-background hangs rather than reporting a shortfall when the foreground is GC-rich
- [2026-09-06] L-48: The ABC snakemake driver stalls after its last real rule; gate downstream work on files, not on the driver
- [2026-09-06] D-16: Make p300 the primary modelling target; train a GM12878 p300 model next
- [2026-09-06] D-17: Use a paired bootstrap for every CRISPR-benchmark comparison
- [2026-09-06] D-18: Close training-element composition: GC-matching the negatives changes nothing
- [2026-09-05] L-41: Why better H3K27ac does not help ABC: rank displacement of the functional elements
- [2026-09-05] L-42: Sequence already knows which elements are not acetylated; the additive trunk ignores it
- [2026-09-05] L-43: Correcting the previous entry: sequence does NOT see the CTCF signature, and my substitute control was also wrong
- [2026-09-05] L-44: The training negative pool was never GC-matched, and the filename says otherwise
- [2026-09-05] L-45: The gate and the asymmetric loss both fail, as the corrected sequence reading predicted

## By tag

- `h3k27ac`: L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14, L-15, L-16, L-17, L-18, L-22, L-24, L-27, L-30, L-31, D-1, D-2, D-3, D-4, D-5
- `5-prime`: L-8, L-9, L-10, L-11, L-13, L-17, L-18, L-28, L-29, D-5, D-6
- `atac`: L-12, L-22, L-25, L-28, L-29, L-33, D-4, D-6, D-8
- `ceiling`: L-2, L-7, L-9, L-10, L-14, L-17, L-18, L-24, D-2
- `architecture`: L-35, L-36, L-37, L-42, L-45, D-3, D-7, D-15
- `negative-result`: L-7, L-9, L-41, L-45, D-13, D-14, D-15, D-18
- `abc`: L-39, L-40, L-41, L-48, D-9, D-16, D-17
- `gm12878`: L-14, L-22, L-24, L-30, L-31, L-35, L-37
- `k562`: L-24, L-30, L-31, L-35, L-37, L-38
- `silent-failure`: L-6, L-15, L-19, L-20, L-21, L-47
- `accessibility`: L-29, L-37, D-6, D-7, D-8
- `bpnetlite`: L-1, L-3, L-8, L-11, D-3
- `chrombpnet`: L-28, L-29, L-33, L-47, D-6
- `evaluation`: L-2, L-3, L-34, D-4, D-10
- `gating`: L-39, L-40, L-42, L-45, D-15
- `loss-weighting`: L-1, L-11, L-13, L-16, D-3
- `negatives`: L-5, L-44, L-47, D-15, D-18
- `residual`: L-15, L-27, L-30, L-31, D-4
- `stratification`: L-17, L-18, L-31, L-35, L-36
- `transferability`: L-14, L-22, L-24, L-31, D-16
- `crispr-benchmark`: L-41, D-9, D-16, D-17
- `fragment-length`: L-12, L-33, L-37, D-8
- `multimodal`: L-27, L-29, L-36, D-7
- `p300`: L-3, L-29, D-13, D-16
- `receptive-field`: L-35, L-36, L-37, D-7
- `target-definition`: L-8, L-9, L-10, D-5
- `testing`: L-19, L-20, L-34, D-10
- `tooling`: L-4, L-6, L-20, L-32
- `controls`: L-43, L-46, D-13
- `ctcf`: L-42, L-43, L-46
- `diagnosis`: L-39, L-40, L-41
- `gc-matching`: L-44, L-47, D-18
- `hyperparameters`: L-1, L-13, L-16
- `loss-design`: L-42, L-45, D-15
- `mnll`: L-8, L-11, L-13
- `negative-control`: L-27, L-30, L-38
- `prediction-was-wrong`: L-10, L-14, L-24
- `profile-head`: L-11, L-36, D-3
- `self-correction`: L-42, L-43, L-46
- `sherlock`: L-4, L-21, L-48
- `slurm`: L-4, L-19, L-21
- `snakemake`: L-39, L-40, L-48
- `statistics`: L-2, L-23, D-17
- `telohaec`: L-26, L-28, D-5
- `training-composition`: L-44, D-15, D-18
- `validation`: L-19, L-28, L-33
- `variance`: L-16, L-18, L-23
- `window-selection`: L-7, D-1, D-2
- `deployment`: L-31, D-14
- `determinism`: L-34, D-10
- `fragment-extension`: L-8, L-9
- `guards`: L-19, L-20
- `inference`: L-38, D-11
- `methodology`: L-15, L-19
- `mtime`: L-39, L-40
- `nucleosome`: L-7, L-12
- `over-prediction`: L-42, L-43
- `overclaim`: L-16, L-17
- `paired-test`: L-27, L-29
- `panel`: L-25, D-5
- `qnorm`: L-41, D-9
- `read-length`: L-28, D-6
- `regression`: L-34, D-10
- `replication`: L-30, L-35
- `reproducibility`: L-16, L-32
- `reversal`: L-36, D-7
- `reverse-complement`: L-38, D-11
- `sequence-branch`: L-44, D-18
- `superset`: L-33, D-8
- `tn5`: L-33, D-8
- `tolerance`: L-34, D-10
- `top-quintile`: L-35, L-37
- `training`: L-1, L-5
- `training-objective`: L-27, L-30
- `activity`: D-9
- `adoption`: D-11
- `all-elements-artifact`: L-17
- `asymmetry`: L-24
- `augmentation`: L-38
- `bash`: L-4
- `bigwig`: L-8
- `bootstrap`: D-17
- `bpnet-gc-background`: L-47
- `buffering`: L-4
- `caching`: L-32
- `causal-inference`: L-46
- `cell-line-contamination`: L-26
- `channels`: L-12
- `checkpoints`: L-21
- `chrM`: L-28
- `churn`: D-12
- `comparability`: D-11
- `conda`: L-40
- `confidence-intervals`: L-23
- `confounding`: L-25
- `contamination`: D-2
- `coordinates`: L-33
- `count-target`: L-3
- `cpg-island`: L-42
- `cross-validation`: L-23
- `data-availability`: L-25
- `data-inspection`: L-12
- `data-provenance`: L-26
- `decisions`: L-6
- `denominator`: L-22
- `depth`: L-28
- `dnase`: L-25
- `documentation`: D-12
- `documentation-drift`: L-44
- `drift`: L-32
- `dry-run`: L-39
- `dynamic-range`: L-43
- `element-centric`: D-1
- `encode`: L-25
- `enrichment`: L-46
- `error-strata`: L-46
- `fragment-channels`: D-14
- `free-win`: L-38
- `generalization`: L-25
- `geo`: L-26
- `gpu`: L-34
- `hang`: L-47
- `index`: L-6
- `injection`: D-9
- `input-definition`: D-6
- `input-design`: D-8
- `known-behaviour`: L-48
- `latent-bug`: L-34
- `leakage`: D-9
- `learnability`: L-10
- `mechanism`: L-41
- `memory`: L-5
- `metric-choice`: D-4
- `model-free-baseline`: L-22
- `multi-task`: D-13
- `mycelium`: L-6
- `normalization`: L-14
- `numbers-manifest`: L-32
- `off-by-factor`: L-3
- `offset`: L-15
- `oom`: L-5
- `orchestration`: L-48
- `organisation`: D-12
- `orphaned-jobs`: L-40
- `owners`: L-21
- `paired-end`: D-5
- `paired-testing`: D-17
- `patching`: L-20
- `pooled-statistics`: L-46
- `positive-result`: D-16
- `prediction-before-measurement`: L-45
- `preemption`: L-21
- `premature-conclusion`: L-36
- `process`: L-40
- `profile-loss`: D-5
- `ranking`: L-41
- `reciprocal`: L-24
- `reliability`: L-2
- `replicates`: L-2
- `report`: L-32
- `reporting`: D-12
- `reporting-standard`: L-23
- `reuse`: L-39
- `sample-definition`: L-26
- `scaling`: L-5
- `scope`: L-25
- `sequence-vs-accessibility`: L-18
- `spearman-brown`: L-2
- `stall`: L-48
- `str-replace`: L-20
- `stranded`: L-3
- `submit-scripts`: L-4
- `superseding`: L-43
- `target-choice`: D-16
- `target-mismatch`: L-15
- `template-mismatch`: L-6
- `trade-off`: D-2
- `training-design`: D-1
- `transfer`: D-14

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
