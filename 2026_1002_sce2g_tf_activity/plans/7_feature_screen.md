# Feature screen for the scE2G activity term (draft, 2026-10-05)

Find deployable element features that carry observed K562 p300 beyond Xu ATAC, using p300 over
all ~155k Xu candidate elements as the selection target. CRISPR labels are not used until one
pre-registered scE2G run at the end. Observed p300 as a feature is worth +0.056 to +0.060 CRISPR
AUPRC (F-056), so "p300 beyond ATAC" is the signal worth recovering.

## Rules

- Deployable only: sequence, and the Xu multiome (fragments, cell x peak matrix, TOBIAS on Xu).
  No ENCODE K562 assay, and no ENCODE cCRE class (built partly from K562 marks).
- No CRISPR label enters any step before the final run.
- Discovery on odd chromosomes, replication on even chromosomes and chrX.

## Feature library (`scripts/7.1`)

| Family | Features | Source |
|---|---|---|
| composition | GC, CpG o/e, CpG count, 16 dinucleotide frequencies, longest A/T run, poly(dA:dT) tracts >= 5 bp (midpoint +/- 250 bp) | hg38 |
| motifs | site count and max match score for each of 879 JASPAR 2024 motifs, sites within +/- 250 bp of the midpoint | TOBIAS BINDetect on Xu peaks (2.2) |
| element shape | element width, mappability (Umap k36) | Xu candidates, hg38_resources |
| fragments | share of insertions from fragments < 100 bp, 180-247 bp, mean log2 fragment length, share of insertions within +/- 50 bp of the midpoint | Xu fragments |
| cells | share of cells with the peak accessible, mean count per accessible cell | scE2G Kendall matrix (5.0 export, 7,821 cells) |
| neighbourhood | Xu elements within 50 kb, summed Xu ATAC within 50 kb (self excluded) | Xu candidates |
| footprints (reference) | TOBIAS `fp_any`, PRINT erythroid and CTCF footprint max | 2.4 tables |
| conservation, repeats | phyloP100way mean and share of bases > 2; share covered by any repeat and by each of 7 RepeatMasker classes; CpG-island share and log distance | UCSC hg38, `$SCRATCH/ucsc_hg38` |

## Screen (`scripts/7.2`)

Target: rank of log observed p300 (2.4 table), residualised on Z = Xu ATAC percentile and its
square, log distance to TSS, promoter flag, mappability, GC and CpG o/e: the 2.4 gate's controls
plus mappability, so the 0.15
bar below means what it meant there. GC is known already (F-055) and dominates without it in Z
(chr21-22 smoke test: partial rho -0.55 beyond ATAC alone; the erythroid PRINT footprint, null on
CRISPR, reaches 0.18 and would pass the bar). Every feature's value without GC and CpG in Z is
reported alongside. Mappability is technical: it was the smoke test's top hit (0.32), as expected
for short-read p300 ChIP in repeats.

1. Univariate: partial Spearman of each feature with the target given Z, on discovery and on
   replication. BH-FDR over the library on discovery; replication intervals from a 1 Mb block
   bootstrap, since neighbouring elements are not independent.
2. Multivariate: gradient-boosted trees and an elastic net on the residual, trained on one
   chromosome half and scored on the other (both directions). Held-out R2 of the residual is the
   share of "p300 beyond ATAC" the library recovers. Out-of-fold tree predictions give one
   composite feature per element.
3. Stability: lasso on 50 half-samples of discovery 1 Mb blocks; selection frequency per feature.

Calibration rows printed with every table: GC, CTCF motif count, erythroid motif count, TOBIAS and
PRINT footprints, whose CRISPR outcomes are known (F-050, F-051, F-055). Footprints reached partial
rho 0.07-0.10 against p300 and moved CRISPR AUPRC by nothing, so that is the floor a feature must
clear.

## Shortlist rule (fix before reading 7.2 output)

A feature goes forward if all hold: same sign and replication interval excluding 0, |partial
rho| >= 0.15 on replication, lasso selection frequency >= 0.8, and |rho| < 0.7 with every
feature already on the list. At most 5, plus the composite.

## Confirmation

One scE2G run through `4.8`: baseline n6 + GC (the F-055 arm, as the reference), n6 + GC +
shortlist, n6 + GC + composite. The composite is fitted beyond GC and CpG, so it is added on top
of GC, never instead of it. Pass rule written into `.living/decisions.md` before the run:
paired-bootstrap interval vs n6 + GC above 0.

## Open

- The composite is trained on K562 p300, so in another cell it is a transferred model with the
  same transfer question as the BPNet line.
- p300 is a proxy. GC points the same way on p300 and on CRISPR (negative on both, F-054), but
  the erythroid footprints pass on p300 and are null on CRISPR (F-050), so a p300 hit is a
  candidate, not a result.
