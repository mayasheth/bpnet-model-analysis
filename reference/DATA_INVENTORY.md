# ENCODE / SRA data inventory

From Maya's spreadsheet `2023_0701 Maya's sequences - ENCODE data.tsv` (recorded 2026-09-05).
Supersedes the K562-only version of this file.

**Processing convention.** DNase-seq and ChIP-seq single-ended runs: download UNFILTERED
alignments and filter/sort/index manually. ATAC-seq: download FILTERED alignments and convert
to tagAlign per ENCODE standards.

Root for ENCODE cell types: `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/<cell>`.

## Assay availability, what constrains the panel

| cell type | DNase | ATAC | H3K27ac | EP300 | CTCF | H3K4me1 | H3K27me3 |
|---|---|---|---|---|---|---|---|
| K562     | yes (x2) | **yes** | yes | **yes** | yes | yes | yes |
| GM12878  | yes | **yes** | yes | **yes** | yes | yes | yes |
| A549     |, | **yes** |, | **yes** |, |, |, |
| HepG2    |, | **yes** |, | **yes** |, |, |, |
| MCF-7    |, | **yes** |, | **yes** |, |, |, |
| TeloHAEC |, | **yes** (4 conditions) | yes (4 conditions) |, |, |, |, |
| WTC11    | yes |, | Mint-ChIP |, | yes | Mint-ChIP | Mint-ChIP |
| HCT116   | yes |, | yes |, | yes | yes | yes |
| Jurkat   | yes |, | yes (GEO) |, | yes (GEO) | yes (GEO) | yes (GEO) |
| H1       | yes |, | yes |, | yes | yes | yes |
| H9       | yes |, | yes |, | yes | yes | yes |
| THP-1    | yes (1 rep, BAM) |, | yes (GEO, 2 reps pe) |, |, |, |, |

**Consequences.** ATAC exists in only three cell types here (K562, GM12878, TeloHAEC), which is
what forces the ATAC-only panel rule. **That rule is in tension with the results rather than a
neutral choice**: DNase beats ATAC as a model input (F-010), the advantage is base-resolution
structure in K562 (F-013), and a DNase-input model is the only one whose predictions clear the
ABC benchmark floor (F-018). A DNase-input panel is internally consistent and reaches three
cell types here (K562, GM12878, THP-1).

**CORRECTED 2026-09-15: DNase is available for FEWER cell types than ATAC, not more.** An
earlier version of this paragraph claimed the opposite and used it to argue a DNase-input panel
was the more deployable choice. That is backwards and it contradicts the project's own premise:
the application target has ATAC and nothing else, which is the entire reason an ATAC-to-DNase
converter was attempted (and closed, F-014). ATAC is the more widely available assay, DNase the
scarcer one. **Consequence for how F-010 and F-018 are read: a DNase-input model is LESS
deployable than an ATAC-input one, so those findings say what information helps H3K27ac
prediction, not that a deployment route has been found.** Locally the two assays reach three
cell types each, so neither is broader in this inventory; the asymmetry is about what can be
obtained for a new target cell type.
**A converter training cell type needs BOTH assays, and only K562 and GM12878 have both**, which
is what blocks any claim about converter portability. **CORRECTED 2026-09-18: EP300 is available in far more than two cell types.** This
file previously said "EP300 exists in exactly two (K562, GM12878)", which was true
of what had been DOWNLOADED and false of what exists. A portal survey on 2026-09-18
found 55 released human EP300 ChIP-seq experiments, and **six cell lines have both
EP300 and ATAC**: K562, GM12878, A549, HepG2, MCF-7 and SK-N-SH, plus eight tissues
(sigmoid colon, transverse colon, ovary, stomach, upper lobe of left lung, esophagus
muscularis mucosa, gastroesophageal sphincter, tibial nerve). A549, HepG2 and MCF-7
were downloaded on 2026-09-18 for the multi-cell-type p300 panel. SK-N-SH was left
out because its only ATAC experiment is retinoic-acid treated while its untreated
EP300 experiments have no matching ATAC, so cell type and treatment would confound.
**The lesson is that "we have data for two cell types" was being read as "data
exists for two cell types" for weeks.** Check the portal before concluding an assay
is unavailable. CTCF, H3K4me1 and H3K27me3 are
available in six cell types, so annotation-based analyses generalise further than modelling
does.

## K562  `Data/ENCODE/K562`

| assay | experiment | run | processed | fold-change bigwig | peaks | peak type |
|---|---|---|---|---|---|---|
| DNase | ENCSR000EOT | pe,se | `ENCFF205FNC`, `ENCFF860XAE`.filtered.sorted.bam | ENCFF414OGC |, |, |
| DNase | ENCSR000EKS | se | `ENCFF987IUK.filtered.sorted.bam` |, |, |, |
| H3K27ac | ENCSR000AKP | se | `ENCFF790GFL.se`, `ENCFF817HMW.se`.filtered.sorted.bam |, (signal on mitra) | `ENCFF544LXB` | pseudoreplicated |
| ATAC | ENCSR868FGK | pe | `ENCFF077FBI`, `ENCFF128WZG`, `ENCFF534DCE`.tn5.sorted.tagAlign.gz |, |, |, |
| CTCF | ENCSR000AKO | se | `ENCFF216XRV`, `ENCFF678LBZ`.filtered.sorted.bam | ENCFF405AYC | `ENCFF519CXF` | optimal IDR |
| EP300 | ENCSR000EGE | se | `ENCFF466WKF`, `ENCFF163FSR`.filtered.sorted.bam | ENCFF636VVR | `ENCFF702XPO` | IDR |
| H3K4me1 | ENCSR000AKS | se | `ENCFF524BOJ`, `ENCFF204MWI`, `ENCFF665JSC`.filtered.sorted.bam | ENCFF607SUJ | `ENCFF135ZLM` | replicated |
| H3K27me3 | ENCSR000AKQ | se | `ENCFF549RYG`, `ENCFF483EPA`, `ENCFF351YGP`.filtered.sorted.bam | ENCFF242ENK | `ENCFF323WOT` | pseudoreplicated |
| IgG control | ENCSR000EHI | se | `ENCFF396DTD.sorted.bam` |, |, |, |

ATAC PE BAMs (`*.pe.bam`) are on `$SCRATCH/atac_pe`; fragment channels need them.

**DO NOT USE** `ENCSR000EOT_archive` (`ENCFF325RTP_sorted.bam`, `ENCFF860XAE_sorted.bam`):
R1 is paired-end but was processed as single-end, and R1 depth >> R2.

## GM12878  `Data/ENCODE/GM12878`

| assay | experiment | run | processed | bigwig | peaks | peak type |
|---|---|---|---|---|---|---|
| DNase | ENCSR000EMT | se | `ENCFF467CXY_sorted.bam`, `ENCFF940NSD_sorted.bam` |, |, |, |
| H3K27ac | ENCSR000AKC | se | `ENCFF645BAL`, `ENCFF865OOP`.filtered.sorted.bam |, | `ENCFF023LTU` | pseudoreplicated |
| ATAC | ENCSR637XSC | pe | `ATAC/ENCFF981FXV`, `ENCFF962FMH`, `ENCFF440GRZ`.tn5.sorted.tagAlign.gz |, |, |, |
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

## Multi-cell-type p300 panel, added 2026-09-18

Downloaded for the leave-one-out p300 transfer experiment (F-024 made p300 the default
activity target). Manifest with sizes and md5-verified accessions:
`2026_0824_H3K27ac_model/reference/p300_panel_manifest.tsv`. Fetched by
`scripts/0.43.fetch_p300_panel.sh`, which verifies md5 rather than existence, so a preempted
transfer cannot leave a truncated BAM that fails later somewhere confusing.

**ATAC is matched across the whole panel and that was deliberate**: paired-ended, 100 bp,
Snyder lab in all five cell types, the same experiment family as K562's ENCSR868FGK. A549 also
has a single-ended 51 bp ATAC experiment (ENCSR220ASC, Reddy) and it was NOT used, because
run type and read length are already a known confound here from TeloHAEC.

**EP300 is NOT matched and this is the panel's main caveat.** Single-ended everywhere, but
36 bp (K562, GM12878, HepG2), 50 bp (MCF-7) and 51 bp (A549), across three labs. Read length
should matter little for a 5-prime-end target, since the 5-prime position is the alignment
start regardless of read length, but antibody and protocol differ by lab and that is a real
batch effect in a model trained across cell types. Any cross-cell-type p300 result should be
checked against the possibility that the model is reading lab rather than biology.

| cell | assay | experiment | run | read | lab | files (`Data/ENCODE/<cell>/<assay>/`) | peaks |
|---|---|---|---|---|---|---|---|
| A549  | EP300 | ENCSR686BQM | se | 51 | Reddy, Duke | `ENCFF639MJZ`, `ENCFF981BEX`, `ENCFF371YBD`.bam (unfiltered) | `ENCFF143OQP` IDR |
| A549  | ATAC  | ENCSR032RGS | pe | 100 | Snyder | `ENCFF607DTB`, `ENCFF701BDT`, `ENCFF616DYV`.bam (filtered) | `ENCFF429FOV` IDR |
| HepG2 | EP300 | ENCSR000EDV | se | 36 | Snyder | `ENCFF922TSG`, `ENCFF713NWW`.bam (unfiltered) | `ENCFF488UHZ` IDR |
| HepG2 | ATAC  | ENCSR291GJU | pe | 100 | Snyder | `ENCFF990VCP`, `ENCFF624SON`, `ENCFF926KFU`.bam (filtered) | `ENCFF915FZC` IDR |
| MCF-7 | EP300 | ENCSR000BTR | se | 50 | Myers, HAIB | `ENCFF490PAU`, `ENCFF075FEU`.bam (unfiltered) | `ENCFF290DJX` IDR |
| MCF-7 | ATAC  | ENCSR422SUG | pe | 100 | Snyder | `ENCFF607OSL`, `ENCFF772EFK`.bam (filtered) | `ENCFF882OVP` IDR |

**BAM choice.** ChIP gets `unfiltered alignments` and is filtered locally, ATAC gets the
pipeline's filtered `alignments`, per the processing convention at the top of this file.
Where a replicate had several BAMs from different ENCODE reprocessings (HepG2 and MCF-7 EP300
each had a 2016/2019 and a 2020-12-26 version), the most recent was taken; no
`preferred_default` flag was set on any of them.

**Not yet processed.** These are raw downloads. Still to do: filter/sort/index the ChIP BAMs,
convert ATAC to tagAlign, build 5-prime bigwigs and per-cell-type element sets, and compare
depth and signal-to-noise against K562 and GM12878 BEFORE training on them.

## Other cell types, H3K27ac and annotations, no ATAC

**WTC11** `Data/ENCODE/WTC11`, DNase ENCSR785ZUI (`ENCFF492WXQ`, `ENCFF715YXX`.sorted.bam);
H3K27ac **Mint-ChIP** ENCSR146DPQ (`ENCFF738QRT`, `ENCFF696GDY`.sorted.bam, peaks
`ENCFF655PNM`); CTCF ENCSR992XTY (peaks `ENCFF112GJQ`, IDR); H3K4me1 ENCSR473GYW
(`ENCFF199JTJ`); H3K27me3 ENCSR418YEV (`ENCFF801DJK`). H3K27ac is Mint-ChIP, excluded from
the panel.

**HCT116** `Data/ENCODE/HCT116`, DNase ENCSR000ENM; H3K27ac ENCSR661KMA (`ENCFF943GHK`,
`ENCFF977FPK`, peaks `ENCFF899XEF`); CTCF ENCSR240PRQ (`ENCFF803RIY`, IDR); H3K4me1
ENCSR161MXP (`ENCFF240LRP`); H3K27me3 ENCSR810BDB (`ENCFF294LZM`). H3K27ac replicates differ
in run type, so no inter-replicate ceiling is computable.

**H1** `Data/ENCODE/H1`, DNase ENCSR000EMU; H3K27ac ENCSR000ANP (`ENCFF120QMN`,
`ENCFF104RJG`, **no peaks**); CTCF ENCSR000AMF (`ENCFF692RPA`, IDR); H3K4me1 ENCSR000ANA
(`ENCFF984DGO`); H3K27me3 ENCSR000ALU (`ENCFF305KNA`).

**H9** `Data/ENCODE/H9`, DNase ENCSR275ICP; H3K27ac ENCSR876RGF (`ENCFF825BCF`,
`ENCFF709WUL`, **no peaks**); CTCF ENCFF963CHU (3 reps, `ENCFF101UJJ`, conservative IDR);
H3K4me1 ENCSR276HBK (`ENCFF188TGA`); H3K27me3 ENCSR792GCH (`ENCFF680AKW`). H3K27ac and the
histone marks are Bing Ren / Roadmap. H9 elements are much wider (mean 798 bp vs ~570-600).

**Jurkat**, DNase ENCSR000EOS in `Data/ENCODE/Jurkat`; everything else is GEO, processed by
jgalante under `Users/jgalante/making_jurkat/...`. H3K27ac GSE155555 (2 reps, peaks
`Users/sheth/Data/SRA/Jurkat/H3K27ac/macs2_peaks.filtered.bed.gz`, top 54k); CTCF GSE130140
(1 rep); H3K27me3 GSE85601 (1 rep); H3K4me1 GSE119439 (`Data/SRA/Jurkat/H3K4me1/
SRR7782877.filtered.sorted.dedup.bam`, top 100k).

**THP-1** - DNase ENCSR896ADX, **alignments ARE on Oak** (corrected 2026-09-13; an earlier
entry here said bigwig only and that THP-1 was not modellable, which was wrong and had ruled
it out as a third cell type):
`Users/sheth/Data/ENCODE/THP1/DNase/AG81591.filtered.bam` (5.8 GB, indexed; the
`.filtered.cram` from regulome.altius.org sits beside it). **One replicate only**, so no
DNase inter-replicate ceiling is computable. That does not block the H3K27ac panel, where the
ceiling that matters is H3K27ac's, but it does block using THP-1 DNase as a CONVERTER TARGET,
since 0.25's shape ceiling needs two replicates.
H3K27ac GSE201352, two replicates, both paired-end, same run type, so the H3K27ac ceiling IS
computable:
`Projects/E2G/THP1/THP1_PRJNA830917/THP1_macrophages_0000min_R1/H3K27ac/SRR18899252.pe.filtered.sorted.dedup.bam`
and `.../R2/H3K27ac/SRR18899253.pe.filtered.sorted.dedup.bam`, with MACS2 peaks alongside R2.
**Paired-end, so read 1 only** when building 5-prime targets.
No ATAC, so THP-1 can be a DNase-input panel member but not a converter training cell type.

## Deployment constraint, what may be used as a MODEL INPUT

The application target is any cell type with ATAC-seq and nothing else. So:

- **ATAC** may be an input. It is the one assay assumed present.
- **p300, CTCF, H3K4me1, H3K27me3 may NOT be inputs**, however many cell types we hold them
  in. They are unavailable in the target, so a model depending on them cannot deploy.
- They MAY be targets, auxiliary training tasks, or analysis annotations, since all of those
  need the data only where we train and evaluate.

This is why a multi-head model (predict H3K27ac and p300 from ATAC + sequence, read only the
H3K27ac head at inference) is preferable to feeding p300 in: it captures p300's information
at training time and needs nothing extra at deployment.
