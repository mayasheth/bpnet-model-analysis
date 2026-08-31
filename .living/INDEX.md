<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-31

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 12 sections | 2026-08-31 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments, SLURM submit scripts (mandatory) |
| decisions.md | 6 entries | 2026-08-30 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric, Paired-end H3K27ac targets use read 1 only, not both mates |
| learnings.md | 29 entries (large — read selectively) | 2026-08-31 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each, Peak count scales the negative pool, which can OOM by 50x |
| findings/ | 3 findings across 4 topics | 2026-08-29 | predicting-regulatory-element-function-at-scale, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, how-enhancers-control-gene-expression |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-31 (heuristic)

## Tag clusters

- **h3k27ac** (20 entries) — D-1, D-2, D-3, D-4, D-5
- **5-prime** (11 entries) — L-18, L-28, L-29, D-5, D-6
- **ceiling** (9 entries) — L-14, L-17, L-18, L-24, D-2
- **atac** (7 entries) — L-25, L-28, L-29, D-4, D-6
- **bpnetlite** (5 entries) — L-1, L-3, L-8, L-11, D-3
- **loss-weighting** (5 entries) — L-1, L-11, L-13, L-16, D-3

## Most recent (10)

- [2026-08-31] L-28: The 5' ATAC rebuild works, and the read-length artifact is confirmed quantitatively
- [2026-08-31] L-29: The 5' ATAC input helps the ATAC-only model and does nothing for multimodal
- [2026-08-30] L-27: The residual objective helps only when the input is blind to accessibility
- [2026-08-30] D-6: Accessibility inputs should be 5'-end insertion counts (ChromBPNet convention), and ours are not
- [2026-08-29] L-25: ATAC-only is the right rule, and it cuts the panel from five new cell types to one
- [2026-08-29] L-26: TeloHAEC's recorded blocker was the wrong blocker — the real one is a second cell line in the directory
- [2026-08-29] D-5: Paired-end H3K27ac targets use read 1 only, not both mates
- [2026-08-27] L-22: The ATAC-only model transfers fine; only sequence collapses
- [2026-08-27] L-23: Confidence intervals confirm the main claims and kill one of them
- [2026-08-27] L-24: Transfer is strongly asymmetric: the K562->GM12878 collapse is largely a property of the target cell type, not of sequence

## By tag

- `h3k27ac`: L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14, L-15, L-16, L-17, L-18, L-22, L-24, L-27, D-1, D-2, D-3, D-4, D-5
- `5-prime`: L-8, L-9, L-10, L-11, L-13, L-17, L-18, L-28, L-29, D-5, D-6
- `ceiling`: L-2, L-7, L-9, L-10, L-14, L-17, L-18, L-24, D-2
- `atac`: L-12, L-22, L-25, L-28, L-29, D-4, D-6
- `bpnetlite`: L-1, L-3, L-8, L-11, D-3
- `loss-weighting`: L-1, L-11, L-13, L-16, D-3
- `silent-failure`: L-6, L-15, L-19, L-20, L-21
- `target-definition`: L-8, L-9, L-10, D-5
- `chrombpnet`: L-28, L-29, D-6
- `evaluation`: L-2, L-3, D-4
- `gm12878`: L-14, L-22, L-24
- `hyperparameters`: L-1, L-13, L-16
- `mnll`: L-8, L-11, L-13
- `prediction-was-wrong`: L-10, L-14, L-24
- `residual`: L-15, L-27, D-4
- `slurm`: L-4, L-19, L-21
- `telohaec`: L-26, L-28, D-5
- `tooling`: L-4, L-6, L-20
- `transferability`: L-14, L-22, L-24
- `variance`: L-16, L-18, L-23
- `window-selection`: L-7, D-1, D-2
- `accessibility`: L-29, D-6
- `fragment-extension`: L-8, L-9
- `guards`: L-19, L-20
- `methodology`: L-15, L-19
- `multimodal`: L-27, L-29
- `negative-result`: L-7, L-9
- `nucleosome`: L-7, L-12
- `overclaim`: L-16, L-17
- `p300`: L-3, L-29
- `paired-test`: L-27, L-29
- `panel`: L-25, D-5
- `profile-head`: L-11, D-3
- `read-length`: L-28, D-6
- `sherlock`: L-4, L-21
- `statistics`: L-2, L-23
- `stratification`: L-17, L-18
- `testing`: L-19, L-20
- `training`: L-1, L-5
- `validation`: L-19, L-28
- `all-elements-artifact`: L-17
- `architecture`: D-3
- `asymmetry`: L-24
- `bash`: L-4
- `bigwig`: L-8
- `buffering`: L-4
- `cell-line-contamination`: L-26
- `channels`: L-12
- `checkpoints`: L-21
- `chrM`: L-28
- `confidence-intervals`: L-23
- `confounding`: L-25
- `contamination`: D-2
- `count-target`: L-3
- `cross-validation`: L-23
- `data-availability`: L-25
- `data-inspection`: L-12
- `data-provenance`: L-26
- `decisions`: L-6
- `denominator`: L-22
- `depth`: L-28
- `dnase`: L-25
- `element-centric`: D-1
- `encode`: L-25
- `fragment-length`: L-12
- `generalization`: L-25
- `geo`: L-26
- `index`: L-6
- `input-definition`: D-6
- `k562`: L-24
- `learnability`: L-10
- `memory`: L-5
- `metric-choice`: D-4
- `model-free-baseline`: L-22
- `mycelium`: L-6
- `negative-control`: L-27
- `negatives`: L-5
- `normalization`: L-14
- `off-by-factor`: L-3
- `offset`: L-15
- `oom`: L-5
- `owners`: L-21
- `paired-end`: D-5
- `patching`: L-20
- `preemption`: L-21
- `profile-loss`: D-5
- `reciprocal`: L-24
- `reliability`: L-2
- `replicates`: L-2
- `reporting-standard`: L-23
- `reproducibility`: L-16
- `sample-definition`: L-26
- `scaling`: L-5
- `scope`: L-25
- `sequence-vs-accessibility`: L-18
- `spearman-brown`: L-2
- `str-replace`: L-20
- `stranded`: L-3
- `submit-scripts`: L-4
- `target-mismatch`: L-15
- `template-mismatch`: L-6
- `trade-off`: D-2
- `training-design`: D-1
- `training-objective`: L-27

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
