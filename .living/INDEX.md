<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-09-03

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 13 sections | 2026-09-03 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 6 entries | 2026-08-30 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric, Paired-end H3K27ac targets use read 1 only, not both mates |
| learnings.md | 38 entries (large — read selectively) | 2026-09-03 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 3 findings across 4 topics | 2026-09-01 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-09-03 (heuristic)

## Tag clusters

- **h3k27ac** (22 entries) — D-1, D-2, D-3, D-4, D-5
- **5-prime** (11 entries) — L-18, L-28, L-29, D-5, D-6
- **ceiling** (9 entries) — L-14, L-17, L-18, L-24, D-2
- **atac** (8 entries) — L-28, L-29, L-33, D-4, D-6
- **gm12878** (7 entries) — L-24, L-30, L-31, L-35, L-37
- **k562** (6 entries) — L-30, L-31, L-35, L-37, L-38

## Most recent (10)

- [2026-09-03] L-36: The receptive-field conclusion inverted when the second input mode finished
- [2026-09-03] L-37: Both architecture changes that work act on the accessibility input
- [2026-09-03] L-38: Test-time reverse-complement averaging is free, always positive, and carries its own control
- [2026-09-02] L-34: Byte-identity is the wrong regression invariant for GPU inference
- [2026-09-02] L-35: A 4.2 kb receptive field buys dead-vs-active separation and nothing within active elements
- [2026-09-01] L-30: The residual-objective result reproduces in GM12878, so it is the objective and not K562
- [2026-09-01] L-31: Residual and multimodal transfer equivalently; the significant difference is in the stratum we distrust
- [2026-09-01] L-32: A generated numbers manifest caught three wrong figures that proofreading missed
- [2026-09-01] L-33: Measure a coordinate shift by correlation instead of trusting the documented one
- [2026-08-31] L-28: The 5' ATAC rebuild works, and the read-length artifact is confirmed quantitatively

## By tag

- `h3k27ac`: L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14, L-15, L-16, L-17, L-18, L-22, L-24, L-27, L-30, L-31, D-1, D-2, D-3, D-4, D-5
- `5-prime`: L-8, L-9, L-10, L-11, L-13, L-17, L-18, L-28, L-29, D-5, D-6
- `ceiling`: L-2, L-7, L-9, L-10, L-14, L-17, L-18, L-24, D-2
- `atac`: L-12, L-22, L-25, L-28, L-29, L-33, D-4, D-6
- `gm12878`: L-14, L-22, L-24, L-30, L-31, L-35, L-37
- `k562`: L-24, L-30, L-31, L-35, L-37, L-38
- `bpnetlite`: L-1, L-3, L-8, L-11, D-3
- `loss-weighting`: L-1, L-11, L-13, L-16, D-3
- `residual`: L-15, L-27, L-30, L-31, D-4
- `silent-failure`: L-6, L-15, L-19, L-20, L-21
- `stratification`: L-17, L-18, L-31, L-35, L-36
- `architecture`: L-35, L-36, L-37, D-3
- `chrombpnet`: L-28, L-29, L-33, D-6
- `evaluation`: L-2, L-3, L-34, D-4
- `target-definition`: L-8, L-9, L-10, D-5
- `tooling`: L-4, L-6, L-20, L-32
- `transferability`: L-14, L-22, L-24, L-31
- `accessibility`: L-29, L-37, D-6
- `fragment-length`: L-12, L-33, L-37
- `hyperparameters`: L-1, L-13, L-16
- `mnll`: L-8, L-11, L-13
- `multimodal`: L-27, L-29, L-36
- `negative-control`: L-27, L-30, L-38
- `prediction-was-wrong`: L-10, L-14, L-24
- `profile-head`: L-11, L-36, D-3
- `receptive-field`: L-35, L-36, L-37
- `slurm`: L-4, L-19, L-21
- `telohaec`: L-26, L-28, D-5
- `testing`: L-19, L-20, L-34
- `validation`: L-19, L-28, L-33
- `variance`: L-16, L-18, L-23
- `window-selection`: L-7, D-1, D-2
- `fragment-extension`: L-8, L-9
- `guards`: L-19, L-20
- `methodology`: L-15, L-19
- `negative-result`: L-7, L-9
- `nucleosome`: L-7, L-12
- `overclaim`: L-16, L-17
- `p300`: L-3, L-29
- `paired-test`: L-27, L-29
- `panel`: L-25, D-5
- `read-length`: L-28, D-6
- `replication`: L-30, L-35
- `reproducibility`: L-16, L-32
- `sherlock`: L-4, L-21
- `statistics`: L-2, L-23
- `top-quintile`: L-35, L-37
- `training`: L-1, L-5
- `training-objective`: L-27, L-30
- `all-elements-artifact`: L-17
- `asymmetry`: L-24
- `augmentation`: L-38
- `bash`: L-4
- `bigwig`: L-8
- `buffering`: L-4
- `caching`: L-32
- `cell-line-contamination`: L-26
- `channels`: L-12
- `checkpoints`: L-21
- `chrM`: L-28
- `confidence-intervals`: L-23
- `confounding`: L-25
- `contamination`: D-2
- `coordinates`: L-33
- `count-target`: L-3
- `cross-validation`: L-23
- `data-availability`: L-25
- `data-inspection`: L-12
- `data-provenance`: L-26
- `decisions`: L-6
- `denominator`: L-22
- `deployment`: L-31
- `depth`: L-28
- `determinism`: L-34
- `dnase`: L-25
- `drift`: L-32
- `element-centric`: D-1
- `encode`: L-25
- `free-win`: L-38
- `generalization`: L-25
- `geo`: L-26
- `gpu`: L-34
- `index`: L-6
- `inference`: L-38
- `input-definition`: D-6
- `latent-bug`: L-34
- `learnability`: L-10
- `memory`: L-5
- `metric-choice`: D-4
- `model-free-baseline`: L-22
- `mycelium`: L-6
- `negatives`: L-5
- `normalization`: L-14
- `numbers-manifest`: L-32
- `off-by-factor`: L-3
- `offset`: L-15
- `oom`: L-5
- `owners`: L-21
- `paired-end`: D-5
- `patching`: L-20
- `preemption`: L-21
- `premature-conclusion`: L-36
- `profile-loss`: D-5
- `reciprocal`: L-24
- `regression`: L-34
- `reliability`: L-2
- `replicates`: L-2
- `report`: L-32
- `reporting-standard`: L-23
- `reversal`: L-36
- `reverse-complement`: L-38
- `sample-definition`: L-26
- `scaling`: L-5
- `scope`: L-25
- `sequence-vs-accessibility`: L-18
- `spearman-brown`: L-2
- `str-replace`: L-20
- `stranded`: L-3
- `submit-scripts`: L-4
- `superset`: L-33
- `target-mismatch`: L-15
- `template-mismatch`: L-6
- `tn5`: L-33
- `tolerance`: L-34
- `trade-off`: D-2
- `training-design`: D-1

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
