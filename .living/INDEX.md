<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-27

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 11 sections | 2026-08-27 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 4 entries | 2026-08-25 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric |
| learnings.md | 23 entries | 2026-08-27 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 3 findings across 4 topics | 2026-08-27 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-27 (heuristic)

## Tag clusters

- **h3k27ac** (17 entries) — L-22, D-1, D-2, D-3, D-4
- **ceiling** (8 entries) — L-10, L-14, L-17, L-18, D-2
- **5-prime** (7 entries) — L-10, L-11, L-13, L-17, L-18
- **bpnetlite** (5 entries) — L-1, L-3, L-8, L-11, D-3
- **loss-weighting** (5 entries) — L-1, L-11, L-13, L-16, D-3
- **silent-failure** (5 entries) — L-6, L-15, L-19, L-20, L-21

## Most recent (10)

- [2026-08-27] L-22: The ATAC-only model transfers fine; only sequence collapses
- [2026-08-27] L-23: Confidence intervals confirm the main claims and kill one of them
- [2026-08-26] L-13: The 5' count_loss_weight optimum is 10, with a genuine interior peak
- [2026-08-26] L-14: GM12878 is an EASIER H3K27ac target than K562, which makes the transfer failure worse
- [2026-08-26] L-15: The offset model must be trained on the same target as the residual model
- [2026-08-26] L-16: Correction: the count_loss_weight "optimum" was within run-to-run noise
- [2026-08-26] L-17: Correction: the 5' target's gain on all-elements is mostly the ceiling rising
- [2026-08-26] L-18: The 5' switch is a real but modest gain — biggest for sequence-only
- [2026-08-26] L-19: A guard written to prevent a silent failure failed itself, because it was never tested
- [2026-08-26] L-20: A naive str.replace on a 2-character token corrupted a guard into always-fail

## By tag

- `h3k27ac`: L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14, L-15, L-16, L-17, L-18, L-22, D-1, D-2, D-3, D-4
- `ceiling`: L-2, L-7, L-9, L-10, L-14, L-17, L-18, D-2
- `5-prime`: L-8, L-9, L-10, L-11, L-13, L-17, L-18
- `bpnetlite`: L-1, L-3, L-8, L-11, D-3
- `loss-weighting`: L-1, L-11, L-13, L-16, D-3
- `silent-failure`: L-6, L-15, L-19, L-20, L-21
- `atac`: L-12, L-22, D-4
- `evaluation`: L-2, L-3, D-4
- `hyperparameters`: L-1, L-13, L-16
- `mnll`: L-8, L-11, L-13
- `slurm`: L-4, L-19, L-21
- `target-definition`: L-8, L-9, L-10
- `tooling`: L-4, L-6, L-20
- `variance`: L-16, L-18, L-23
- `window-selection`: L-7, D-1, D-2
- `fragment-extension`: L-8, L-9
- `gm12878`: L-14, L-22
- `guards`: L-19, L-20
- `methodology`: L-15, L-19
- `negative-result`: L-7, L-9
- `nucleosome`: L-7, L-12
- `overclaim`: L-16, L-17
- `prediction-was-wrong`: L-10, L-14
- `profile-head`: L-11, D-3
- `residual`: L-15, D-4
- `sherlock`: L-4, L-21
- `statistics`: L-2, L-23
- `stratification`: L-17, L-18
- `testing`: L-19, L-20
- `training`: L-1, L-5
- `transferability`: L-14, L-22
- `all-elements-artifact`: L-17
- `architecture`: D-3
- `bash`: L-4
- `bigwig`: L-8
- `buffering`: L-4
- `channels`: L-12
- `checkpoints`: L-21
- `confidence-intervals`: L-23
- `contamination`: D-2
- `count-target`: L-3
- `cross-validation`: L-23
- `data-inspection`: L-12
- `decisions`: L-6
- `denominator`: L-22
- `element-centric`: D-1
- `fragment-length`: L-12
- `index`: L-6
- `learnability`: L-10
- `memory`: L-5
- `metric-choice`: D-4
- `model-free-baseline`: L-22
- `mycelium`: L-6
- `negatives`: L-5
- `normalization`: L-14
- `off-by-factor`: L-3
- `offset`: L-15
- `oom`: L-5
- `owners`: L-21
- `p300`: L-3
- `patching`: L-20
- `preemption`: L-21
- `reliability`: L-2
- `replicates`: L-2
- `reporting-standard`: L-23
- `reproducibility`: L-16
- `scaling`: L-5
- `sequence-vs-accessibility`: L-18
- `spearman-brown`: L-2
- `str-replace`: L-20
- `stranded`: L-3
- `submit-scripts`: L-4
- `target-mismatch`: L-15
- `template-mismatch`: L-6
- `trade-off`: D-2
- `training-design`: D-1
- `validation`: L-19

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
