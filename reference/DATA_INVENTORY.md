# ENCODE / SRA data inventory

From Maya's spreadsheet `2023_0701 Maya's sequences - ENCODE data.tsv` (recorded 2026-09-05).
Supersedes the K562-only version of this file.

**Processing convention.** DNase-seq and ChIP-seq single-ended runs: download UNFILTERED
alignments and filter/sort/index manually. ATAC-seq: download FILTERED alignments and convert
to tagAlign per ENCODE standards.

Root for ENCODE cell types: `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/<cell>`.

## Assay availability — what constrains the panel

| cell type | DNase | ATAC | H3K27ac | EP300 | CTCF | H3K4me1 | H3K27me3 |
|---|---|---|---|---|---|---|---|
| K562     | yes (x2) | **yes** | yes | **yes** | yes | yes | yes |
| GM12878  | yes | **yes** | yes | **yes** | yes | yes | yes |
| TeloHAEC | — | **yes** (4 conditions) | yes (4 conditions) | — | — | — | — |
| WTC11    | yes | — | Mint-ChIP | — | yes | Mint-ChIP | Mint-ChIP |
| HCT116   | yes | — | yes | — | yes | yes | yes |
| Jurkat   | yes | — | yes (GEO) | — | yes (GEO) | yes (GEO) | yes (GEO) |
| H1       | yes | — | yes | — | yes | yes | yes |
| H9       | yes | — | yes | — | yes | yes | yes |
| THP-1    | yes (bigwig only) | — | yes (GEO) | — | — | — | — |

**Consequences.** ATAC exists in only three cell types (K562, GM12878, TeloHAEC), which is what
forces the ATAC-only panel rule. **EP300 exists in exactly two (K562, GM12878)** — enough to
test p300 transferability the same way H3K27ac was tested. CTCF, H3K4me1 and H3K27me3 are
available in six cell types, so annotation-based analyses generalise further than modelling
does.

## K562  `Data/ENCODE/K562`

| assay | experiment | run | processed | fold-change bigwig | peaks | peak type |
|---|---|---|---|---|---|---|
| DNase | ENCSR000EOT | pe,se | `ENCFF205FNC`, `ENCFF860XAE`.filtered.sorted.bam | ENCFF414OGC | — | — |
| DNase | ENCSR000EKS | se | `ENCFF987IUK.filtered.sorted.bam` | — | — | — |
| H3K27ac | ENCSR000AKP | se | `ENCFF790GFL.se`, `ENCFF817HMW.se`.filtered.sorted.bam | — (signal on mitra) | `ENCFF544LXB` | pseudoreplicated |
| ATAC | ENCSR868FGK | pe | `ENCFF077FBI`, `ENCFF128WZG`, `ENCFF534DCE`.tn5.sorted.tagAlign.gz | — | — | — |
| CTCF | ENCSR000AKO | se | `ENCFF216XRV`, `ENCFF678LBZ`.filtered.sorted.bam | ENCFF405AYC | `ENCFF519CXF` | optimal IDR |
| EP300 | ENCSR000EGE | se | `ENCFF466WKF`, `ENCFF163FSR`.filtered.sorted.bam | ENCFF636VVR | `ENCFF702XPO` | IDR |
| H3K4me1 | ENCSR000AKS | se | `ENCFF524BOJ`, `ENCFF204MWI`, `ENCFF665JSC`.filtered.sorted.bam | ENCFF607SUJ | `ENCFF135ZLM` | replicated |
| H3K27me3 | ENCSR000AKQ | se | `ENCFF549RYG`, `ENCFF483EPA`, `ENCFF351YGP`.filtered.sorted.bam | ENCFF242ENK | `ENCFF323WOT` | pseudoreplicated |
| IgG control | ENCSR000EHI | se | `ENCFF396DTD.sorted.bam` | — | — | — |

ATAC PE BAMs (`*.pe.bam`) are on `$SCRATCH/atac_pe`; fragment channels need them.

**DO NOT USE** `ENCSR000EOT_archive` (`ENCFF325RTP_sorted.bam`, `ENCFF860XAE_sorted.bam`):
R1 is paired-end but was processed as single-end, and R1 depth >> R2.

## GM12878  `Data/ENCODE/GM12878`

| assay | experiment | run | processed | bigwig | peaks | peak type |
|---|---|---|---|---|---|---|
| DNase | ENCSR000EMT | se | `ENCFF467CXY_sorted.bam`, `ENCFF940NSD_sorted.bam` | — | — | — |
| H3K27ac | ENCSR000AKC | se | `ENCFF645BAL`, `ENCFF865OOP`.filtered.sorted.bam | — | `ENCFF023LTU` | pseudoreplicated |
| ATAC | ENCSR637XSC | pe | `ATAC/ENCFF981FXV`, `ENCFF962FMH`, `ENCFF440GRZ`.tn5.sorted.tagAlign.gz | — | — | — |
| EP300 | ENCSR000DZG | se | `EP300/ENCFF515HYM`, `EP300/ENCFF215GSQ`.filtered.sorted.bam | ENCFF545BXW | `EP300/ENCFF926AKK.bed.gz` | IDR |
| CTCF | ENCSR000AKB | se | `ENCFF067RMO`, `ENCFF551LHV`.filtered.sorted.bam | ENCFF734CUT | `ENCFF797SDL` | IDR |
| H3K4me1 | ENCSR000AKF | se | `ENCFF757IRH`, `ENCFF579EEP`.filtered.sorted.bam | ENCFF564KBE | `ENCFF321BVG` | pseudoreplicated |
| H3K27me3 | ENCSR000DRX | se | `ENCFF867JWR`, `ENCFF539NLB`.filtered.sorted.bam | ENCFF486WAD | `ENCFF695ETB` | replicated |

GM12878 ATAC `.bam` files in `ATAC/` are the paired-end originals (verified PE,
coordinate-sorted); the fragment channels were built from them.

## TeloHAEC  `/oak/stanford/groups/engreitz/Projects/E2G/endothelial_cells/TeloHAEC_GSE210489_GSE210491/<condition>`

Four conditions: `TeloHAEC_ctrl`, `TeloHAEC_IL1b`, `TeloHAEC_TNFa`, `TeloHAEC_VEGF_ctrl`
(the last is "no VEGF"). ATAC and H3K27ac only. All paired-end, Engreitz-processed.

| condition | ATAC reps | H3K27ac reps |
|---|---|---|
| ctrl | 6 (`SRR20809416/419/420/430/431/433`) | 4 (`SRR20810532/533/544/545`) |
| +IL1b | 3 (`SRR20809411/414/415`) | 2 (`SRR20810530/531`) |
| +TNFa | 3 (`SRR20809409/410/412`) | 2 (`SRR20810528/529`) |
| no VEGF | 3 (`SRR20809407/408/413`) | 2 (`SRR20810526/527`) |

**`TeloHAEC_ctrl/ATAC` also contains 3 EA.hy926 files** (`SRR20809434/435/436`) under the same
sample name. Always use explicit accession lists, never a directory glob.

Caveats for modelling: 35-36 bp reads against K562/GM12878's 94-95 bp, shallower libraries,
and ATAC-derived elements. Element derivation was shown not to matter (p = 0.83), but read
length and depth remain confounded.

## Other cell types — H3K27ac and annotations, no ATAC

**WTC11** `Data/ENCODE/WTC11` — DNase ENCSR785ZUI (`ENCFF492WXQ`, `ENCFF715YXX`.sorted.bam);
H3K27ac **Mint-ChIP** ENCSR146DPQ (`ENCFF738QRT`, `ENCFF696GDY`.sorted.bam, peaks
`ENCFF655PNM`); CTCF ENCSR992XTY (peaks `ENCFF112GJQ`, IDR); H3K4me1 ENCSR473GYW
(`ENCFF199JTJ`); H3K27me3 ENCSR418YEV (`ENCFF801DJK`). H3K27ac is Mint-ChIP, excluded from
the panel.

**HCT116** `Data/ENCODE/HCT116` — DNase ENCSR000ENM; H3K27ac ENCSR661KMA (`ENCFF943GHK`,
`ENCFF977FPK`, peaks `ENCFF899XEF`); CTCF ENCSR240PRQ (`ENCFF803RIY`, IDR); H3K4me1
ENCSR161MXP (`ENCFF240LRP`); H3K27me3 ENCSR810BDB (`ENCFF294LZM`). H3K27ac replicates differ
in run type, so no inter-replicate ceiling is computable.

**H1** `Data/ENCODE/H1` — DNase ENCSR000EMU; H3K27ac ENCSR000ANP (`ENCFF120QMN`,
`ENCFF104RJG`, **no peaks**); CTCF ENCSR000AMF (`ENCFF692RPA`, IDR); H3K4me1 ENCSR000ANA
(`ENCFF984DGO`); H3K27me3 ENCSR000ALU (`ENCFF305KNA`).

**H9** `Data/ENCODE/H9` — DNase ENCSR275ICP; H3K27ac ENCSR876RGF (`ENCFF825BCF`,
`ENCFF709WUL`, **no peaks**); CTCF ENCFF963CHU (3 reps, `ENCFF101UJJ`, conservative IDR);
H3K4me1 ENCSR276HBK (`ENCFF188TGA`); H3K27me3 ENCSR792GCH (`ENCFF680AKW`). H3K27ac and the
histone marks are Bing Ren / Roadmap. H9 elements are much wider (mean 798 bp vs ~570-600).

**Jurkat** — DNase ENCSR000EOS in `Data/ENCODE/Jurkat`; everything else is GEO, processed by
jgalante under `Users/jgalante/making_jurkat/...`. H3K27ac GSE155555 (2 reps, peaks
`Users/sheth/Data/SRA/Jurkat/H3K27ac/macs2_peaks.filtered.bed.gz`, top 54k); CTCF GSE130140
(1 rep); H3K27me3 GSE85601 (1 rep); H3K4me1 GSE119439 (`Data/SRA/Jurkat/H3K4me1/
SRR7782877.filtered.sorted.dedup.bam`, top 100k).

**THP-1** — DNase ENCSR896ADX, **bigwig only** (filtered alignments came from
regulome.altius.org as `AG81591.filtered.cram`). H3K27ac GSE201352 under
`Projects/E2G/THP1/THP1_PRJNA830917/THP1_macrophages_0000min_R{1,2}/H3K27ac/`.
No local DNase alignments, so THP-1 is not modellable without reprocessing.

## Deployment constraint — what may be used as a MODEL INPUT

The application target is any cell type with ATAC-seq and nothing else. So:

- **ATAC** may be an input. It is the one assay assumed present.
- **p300, CTCF, H3K4me1, H3K27me3 may NOT be inputs**, however many cell types we hold them
  in. They are unavailable in the target, so a model depending on them cannot deploy.
- They MAY be targets, auxiliary training tasks, or analysis annotations, since all of those
  need the data only where we train and evaluate.

This is why a multi-head model (predict H3K27ac and p300 from ATAC + sequence, read only the
H3K27ac head at inference) is preferable to feeding p300 in: it captures p300's information
at training time and needs nothing extra at deployment.
