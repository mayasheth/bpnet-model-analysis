<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-25

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 7 sections | 2026-08-25 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 0 entries | 2026-08-25 | — |
| learnings.md | 5 entries | 2026-08-25 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 1 findings across 4 topics | 2026-08-25 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-25 (heuristic)

## Tag clusters

- **bpnetlite** (2 entries) — L-1, L-3
- **evaluation** (2 entries) — L-2, L-3
- **training** (2 entries) — L-1, L-5

## Most recent (10)

- [2026-08-24] L-1: count_loss_weight must be calibrated to the actual loss magnitudes, not copied
- [2026-08-24] L-2: Inter-replicate r is not a model performance ceiling without two corrections
- [2026-08-24] L-3: bpnetlite's count target sums ALL channels, not one strand
- [2026-08-24] L-4: Three Sherlock/SLURM traps that cost a job each
- [2026-08-24] L-5: Peak count scales the negative pool, which can OOM by 50x

## By tag

- `bpnetlite`: L-1, L-3
- `evaluation`: L-2, L-3
- `training`: L-1, L-5
- `bash`: L-4
- `buffering`: L-4
- `ceiling`: L-2
- `count-target`: L-3
- `h3k27ac`: L-1
- `hyperparameters`: L-1
- `loss-weighting`: L-1
- `memory`: L-5
- `negatives`: L-5
- `off-by-factor`: L-3
- `oom`: L-5
- `p300`: L-3
- `reliability`: L-2
- `replicates`: L-2
- `scaling`: L-5
- `sherlock`: L-4
- `slurm`: L-4
- `spearman-brown`: L-2
- `statistics`: L-2
- `stranded`: L-3
- `submit-scripts`: L-4
- `tooling`: L-4

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
