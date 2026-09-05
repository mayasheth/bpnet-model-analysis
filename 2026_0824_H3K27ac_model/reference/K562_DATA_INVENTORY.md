# K562 data inventory (from Maya's spreadsheet, 2026-09-05)

All under `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562` unless noted.
Peak-call type matters: IDR-thresholded is stricter than pseudoreplicated.

| assay | experiment | run | processed | peaks | peak type | notes |
|---|---|---|---|---|---|---|
| DNase-seq | ENCSR000EOT | pe, se | `ENCFF205FNC.filtered.sorted.bam`, `ENCFF860XAE.filtered.sorted.bam` | — | — | correct processing; SPOT1 0.55 / 0.52 |
| DNase-seq | ENCSR000EKS | se | `ENCFF987IUK.filtered.sorted.bam` | — | — | SPOT1 0.53 / 0.35 |
| H3K27ac | ENCSR000AKP | se | `ENCFF790GFL.se.filtered.sorted.bam`, `ENCFF817HMW.se.filtered.sorted.bam` | `ENCFF544LXB.bed.gz` | pseudoreplicated | the model target |
| ATAC-seq | ENCSR868FGK | pe | `ENCFF077FBI/ENCFF128WZG/ENCFF534DCE.tn5.sorted.tagAlign.gz` | — | — | PE BAMs `*.pe.bam`; fragment channels need these |
| CTCF | ENCSR000AKO | se | `ENCFF216XRV.filtered.sorted.bam`, `ENCFF678LBZ.filtered.sorted.bam` | `ENCFF519CXF.bed.gz` | **optimal IDR** | tests the over-prediction phenotype |
| EP300 | ENCSR000EGE | se | `ENCFF466WKF.filtered.sorted.bam`, `ENCFF163FSR.filtered.sorted.bam` | `ENCFF702XPO.bed.gz` | IDR | candidate activity term |
| H3K4me1 | ENCSR000AKS | se | `ENCFF524BOJ/ENCFF204MWI/ENCFF665JSC.filtered.sorted.bam` | `ENCFF135ZLM.bed.gz` | replicated | primed-enhancer mark |
| H3K27me3 | ENCSR000AKQ | se | `ENCFF549RYG/ENCFF483EPA/ENCFF351YGP.filtered.sorted.bam` | `ENCFF323WOT.bed.gz` | pseudoreplicated | repressive |
| IgG control | ENCSR000EHI | se | `ENCFF396DTD.sorted.bam` | — | — | ChIP background |

**Do not use** `ENCSR000EOT_archive` (`ENCFF325RTP_sorted.bam`, `ENCFF860XAE_sorted.bam`):
R1 is paired-end but was processed as single-end, and R1 depth >> R2.

Signal bigwigs are on the ENCODE portal or mitra; see the spreadsheet for the fold-change and
signal URLs. `ENCFF702XPO`/`ENCFF519CXF` are the two used so far (`4.11`).
