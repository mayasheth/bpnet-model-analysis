<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-25

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 9 sections | 2026-08-25 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 3 entries | 2026-08-25 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it |
| learnings.md | 8 entries | 2026-08-25 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 1 findings across 4 topics | 2026-08-25 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-25 (heuristic)

## Tag clusters

- **h3k27ac** (6 entries) — L-7, L-8, D-1, D-2, D-3
- **bpnetlite** (4 entries) — L-1, L-3, L-8, D-3
- **ceiling** (3 entries) — L-2, L-7, D-2
- **window-selection** (3 entries) — L-7, D-1, D-2
- **evaluation** (2 entries) — L-2, L-3
- **loss-weighting** (2 entries) — L-1, D-3

## Most recent (10)

- [2026-08-25] L-6: Mycelium's decision-log template and its indexer disagree on heading level
- [2026-08-25] L-7: Excluding the nucleosome-free center makes the H3K27ac target worse, not better
- [2026-08-25] L-8: The IGV display track was the wrong target all along: use 5' ends
- [2026-08-24] L-1: count_loss_weight must be calibrated to the actual loss magnitudes, not copied
- [2026-08-24] L-2: Inter-replicate r is not a model performance ceiling without two corrections
- [2026-08-24] L-3: bpnetlite's count target sums ALL channels, not one strand
- [2026-08-24] L-4: Three Sherlock/SLURM traps that cost a job each
- [2026-08-24] L-5: Peak count scales the negative pool, which can OOM by 50x
- [2026-08-24] D-1: Center H3K27ac training windows on candidate elements, not ChIP peaks
- [2026-08-24] D-2: Counting window is a trade-off between signal and neighbour contamination

## By tag

- `h3k27ac`: L-1, L-7, L-8, D-1, D-2, D-3
- `bpnetlite`: L-1, L-3, L-8, D-3
- `ceiling`: L-2, L-7, D-2
- `window-selection`: L-7, D-1, D-2
- `evaluation`: L-2, L-3
- `loss-weighting`: L-1, D-3
- `tooling`: L-4, L-6
- `training`: L-1, L-5
- `5-prime`: L-8
- `architecture`: D-3
- `bash`: L-4
- `bigwig`: L-8
- `buffering`: L-4
- `contamination`: D-2
- `count-target`: L-3
- `decisions`: L-6
- `element-centric`: D-1
- `fragment-extension`: L-8
- `hyperparameters`: L-1
- `index`: L-6
- `memory`: L-5
- `mnll`: L-8
- `mycelium`: L-6
- `negative-result`: L-7
- `negatives`: L-5
- `nucleosome`: L-7
- `off-by-factor`: L-3
- `oom`: L-5
- `p300`: L-3
- `profile-head`: D-3
- `reliability`: L-2
- `replicates`: L-2
- `scaling`: L-5
- `sherlock`: L-4
- `silent-failure`: L-6
- `slurm`: L-4
- `spearman-brown`: L-2
- `statistics`: L-2
- `stranded`: L-3
- `submit-scripts`: L-4
- `target-definition`: L-8
- `template-mismatch`: L-6
- `trade-off`: D-2
- `training-design`: D-1

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
