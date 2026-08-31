<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-31

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 11 sections | 2026-08-27 | Layout: analyses live at the repo root, not under `analysis/`, Reports, Figures, Environments |
| decisions.md | 6 entries | 2026-08-30 | Center H3K27ac training windows on candidate elements, not ChIP peaks, Counting window is a trade-off between signal and neighbour contamination, Keep the profile head, down-weighted, rather than removing it, Residual correlation beyond ATAC becomes the headline metric |
| learnings.md | 28 entries | 2026-08-31 | count_loss_weight must be calibrated to the actual loss magnitudes, not copied, Inter-replicate r is not a model performance ceiling without two corrections, bpnetlite's count target sums ALL channels, not one strand, Three Sherlock/SLURM traps that cost a job each |
| findings/ | 3 findings across 4 topics | 2026-08-29 | how-enhancers-control-gene-expression, linking-noncoding-variation-to-molecular-function, mapping-regulatory-perturbations-to-phenotype, predicting-regulatory-element-function-at-scale |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-31 (refresh_living_index.py, not the canonical mycelium indexer)

## Tag clusters

- **h3k27ac** (20 entries) — L-1, L-7, L-8, L-9, L-10, L-11, L-13, L-14
- **5-prime** (10 entries) — L-8, L-9, L-10, L-11, L-13, L-17, L-18, L-28
- **ceiling** (9 entries) — L-2, L-7, L-9, L-10, L-14, L-17, L-18, L-24
- **atac** (6 entries) — L-12, L-22, L-25, L-28, D-4, D-6
- **bpnetlite** (5 entries) — L-1, L-3, L-8, L-11, D-3
- **loss-weighting** (5 entries) — L-1, L-11, L-13, L-16, D-3
- **silent-failure** (5 entries) — L-6, L-15, L-19, L-20, L-21
- **target-definition** (4 entries) — L-8, L-9, L-10, D-5
- **evaluation** (3 entries) — L-2, L-3, D-4
- **gm12878** (3 entries) — L-14, L-22, L-24
- **hyperparameters** (3 entries) — L-1, L-13, L-16
- **mnll** (3 entries) — L-8, L-11, L-13
- **prediction-was-wrong** (3 entries) — L-10, L-14, L-24
- **residual** (3 entries) — L-15, L-27, D-4
- **slurm** (3 entries) — L-4, L-19, L-21
- **telohaec** (3 entries) — L-26, L-28, D-5
- **tooling** (3 entries) — L-4, L-6, L-20
- **transferability** (3 entries) — L-14, L-22, L-24

## Most recent entries

- *learnings.md*: The 5' ATAC rebuild works, and the read-length artifact is confirmed quantitatively
- *learnings.md*: The residual objective helps only when the input is blind to accessibility
- *learnings.md*: TeloHAEC's recorded blocker was the wrong blocker — the real one is a second cell line in the directory
- *decisions.md*: Accessibility inputs should be 5'-end insertion counts (ChromBPNet convention), and ours are not
- *decisions.md*: Paired-end H3K27ac targets use read 1 only, not both mates
- *decisions.md*: Residual correlation beyond ATAC becomes the headline metric
<!-- END KNOWLEDGE SUMMARY -->
