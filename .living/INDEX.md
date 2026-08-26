<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-26

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 10 sections | 2026-08-26 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 4 entries | 2026-08-25 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric |
| learnings.md | 17 entries | 2026-08-26 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 3 findings across 4 topics | 2026-08-26 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-26 (heuristic)

## Tag clusters

- **h3k27ac** (15 entries) — L-17, D-1, D-2, D-3, D-4
- **ceiling** (7 entries) — L-9, L-10, L-14, L-17, D-2
- **5-prime** (6 entries) — L-9, L-10, L-11, L-13, L-17
- **bpnetlite** (5 entries) — L-1, L-3, L-8, L-11, D-3
- **loss-weighting** (5 entries) — L-1, L-11, L-13, L-16, D-3
- **evaluation** (3 entries) — L-2, L-3, D-4

## Most recent (10)

- [2026-08-26] L-13: The 5' count_loss_weight optimum is 10, with a genuine interior peak
- [2026-08-26] L-14: GM12878 is an EASIER H3K27ac target than K562, which makes the transfer failure worse
- [2026-08-26] L-15: The offset model must be trained on the same target as the residual model
- [2026-08-26] L-16: Correction: the count_loss_weight "optimum" was within run-to-run noise
- [2026-08-26] L-17: Correction: the 5' target's gain on all-elements is mostly the ceiling rising
- [2026-08-25] L-6: Mycelium's decision-log template and its indexer disagree on heading level
- [2026-08-25] L-7: Excluding the nucleosome-free center makes the H3K27ac target worse, not better
- [2026-08-25] L-8: The IGV display track was the wrong target all along: use 5' ends
- [2026-08-25] L-9: 5'-end vs fragment-extended target: identical ceiling on active elements
- [2026-08-25] L-10: An identical ceiling does not mean identical learnability

## By tag

- `h3k27ac`: L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14, L-15, L-16, L-17, D-1, D-2, D-3, D-4
- `ceiling`: L-2, L-7, L-9, L-10, L-14, L-17, D-2
- `5-prime`: L-8, L-9, L-10, L-11, L-13, L-17
- `bpnetlite`: L-1, L-3, L-8, L-11, D-3
- `loss-weighting`: L-1, L-11, L-13, L-16, D-3
- `evaluation`: L-2, L-3, D-4
- `hyperparameters`: L-1, L-13, L-16
- `mnll`: L-8, L-11, L-13
- `target-definition`: L-8, L-9, L-10
- `window-selection`: L-7, D-1, D-2
- `atac`: L-12, D-4
- `fragment-extension`: L-8, L-9
- `negative-result`: L-7, L-9
- `nucleosome`: L-7, L-12
- `overclaim`: L-16, L-17
- `prediction-was-wrong`: L-10, L-14
- `profile-head`: L-11, D-3
- `residual`: L-15, D-4
- `silent-failure`: L-6, L-15
- `tooling`: L-4, L-6
- `training`: L-1, L-5
- `all-elements-artifact`: L-17
- `architecture`: D-3
- `bash`: L-4
- `bigwig`: L-8
- `buffering`: L-4
- `channels`: L-12
- `contamination`: D-2
- `count-target`: L-3
- `data-inspection`: L-12
- `decisions`: L-6
- `element-centric`: D-1
- `fragment-length`: L-12
- `gm12878`: L-14
- `index`: L-6
- `learnability`: L-10
- `memory`: L-5
- `methodology`: L-15
- `metric-choice`: D-4
- `mycelium`: L-6
- `negatives`: L-5
- `normalization`: L-14
- `off-by-factor`: L-3
- `offset`: L-15
- `oom`: L-5
- `p300`: L-3
- `reliability`: L-2
- `replicates`: L-2
- `reproducibility`: L-16
- `scaling`: L-5
- `sherlock`: L-4
- `slurm`: L-4
- `spearman-brown`: L-2
- `statistics`: L-2
- `stranded`: L-3
- `stratification`: L-17
- `submit-scripts`: L-4
- `target-mismatch`: L-15
- `template-mismatch`: L-6
- `trade-off`: D-2
- `training-design`: D-1
- `transferability`: L-14
- `variance`: L-16

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
