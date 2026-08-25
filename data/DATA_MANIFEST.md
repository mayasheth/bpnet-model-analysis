# Data Manifest

Per the lab data convention (`.living/conventions/engreitz-lab/data-conventions.md`) this
repo tracks **manifests and absolute Oak paths, not bytes**. `data/raw/`,
`data/processed/`, and `data/metadata/` were created by mycelium init and are intentionally
empty — the real files live in `reference/` (shared inputs) and each analysis's own `data/`
subdirectory, both under `$OAK/Users/sheth/EP300_BPNet/`.

`$OAK` = `/oak/stanford/groups/engreitz`.

---

### encsr000ege-p300-chipseq

```yaml
name: encsr000ege-p300-chipseq
type: genomic
source: ENCODE experiment ENCSR000EGE — p300/EP300 ChIP-seq in K562
format: narrowPeak + FASTA
size: 1.6 MB (narrowPeak), 9.2 MB (FASTA)
raw_path: reference/ENCSR000EGE_peaks_inliers.narrowPeak
processed_path: reference/ENCSR000EGE_peaks_inliers.fa
status: processed
access_restrictions: none
tags: [p300, chipseq, k562, encode, training-labels]
```

Training peaks for every p300 model in the repo. `_inliers` denotes the outlier-filtered
peak set. Models train these against GC-matched negatives drawn from
`genomewide_gc_stride_1000_flank_size_1057.gc.bed`.

---

### k562-dnase-candidate-elements

```yaml
name: k562-dnase-candidate-elements
type: genomic
source: K562 DNase-derived candidate regulatory elements
format: narrowPeak, bgzipped BED + tabix index, FASTA + fai
size: 9.6 MB (narrowPeak), 95 MB (FASTA)
raw_path: reference/K562_DNase_candidate_elements.narrowPeak
processed_path: reference/K562_DNase_candidate_elements.fa
metadata_path: reference/K562_DNase_candidate_elements.bed.bgz.tbi
status: processed
access_restrictions: none
tags: [candidate-elements, dnase, k562, prediction-set]
```

The element set predictions and SHAP are computed over — the denominator for "all
elements" performance numbers. Note the p300+/p300- window mismatch flagged under
`p300-k562-bpnet-v1` in `analysis/ANALYSIS_MANIFEST.md`.

---

### gm12878-candidate-elements

```yaml
name: gm12878-candidate-elements
type: genomic
source: GM12878 candidate regulatory elements
format: narrowPeak
size: 9.8 MB
raw_path: 2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak
status: processed
access_restrictions: none
tags: [candidate-elements, gm12878, transferability]
```

GM12878 evaluation element set for the cross-cell-type transferability analysis.

---

### gm12878-p300-chipseq

```yaml
name: gm12878-p300-chipseq
type: genomic
source: ENCODE — EP300 ChIP-seq in GM12878 (ENCFF960OFK plus / ENCFF941MGK minus)
format: BigWig (stranded)
size: 157 MB (plus), 88 MB (minus); merged EP300_plus.bw / EP300_minus.bw ~89 MB each
raw_path: 2026_0606_GM12878_transferability/data/ENCFF960OFK_plus.bw
processed_path: 2026_0606_GM12878_transferability/data/EP300_plus.bw
status: processed
access_restrictions: none
tags: [p300, chipseq, gm12878, encode, stranded, transferability]
```

Stranded p300 signal in GM12878, used both as the GM12878-trained model's target and as
the evaluation truth for K562→GM12878 transfer.

---

### k562-atac-bigwig

```yaml
name: k562-atac-bigwig
type: genomic
source: Merged from 3 K562 ATAC-seq tagAlign replicates
format: BigWig
size: 2.0 GB
raw_path: 2026_0529_multimodal_p300_model/data/atac.bw
status: processed
access_restrictions: none
tags: [atac, accessibility, k562, multimodal-input]
```

Base-resolution accessibility input to the multimodal model. Normalized log1p + z-score at
training time (see `2026_0529_multimodal_p300_model/scripts/train_multimodal_bpnet.py`).

---

### gm12878-atac-bigwig

```yaml
name: gm12878-atac-bigwig
type: genomic
source: GM12878 ATAC-seq
format: BigWig
size: 2.0 GB
raw_path: 2026_0606_GM12878_transferability/data/atac.bw
status: processed
access_restrictions: none
tags: [atac, accessibility, gm12878, multimodal-input, transferability]
```

Accessibility input for GM12878-trained ATAC-only and multimodal models, and for evaluating
K562 multimodal models on GM12878.

---

### motif-compendium-human

```yaml
name: motif-compendium-human
type: other
source: Human Motif Compendium Database; JASPAR2024 CORE vertebrates non-redundant also present
format: MEME
size: 642 KB (compendium), 469 KB (JASPAR2024)
raw_path: reference/MotifCompendium-Database-Human.meme.txt
status: validated
access_restrictions: none
tags: [motifs, meme, pwm, fimo]
```

PWM database for FIMO scanning. `reference/MotifCompendium-test.meme.txt` is a small subset
for pipeline smoke tests. `reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt`
is an alternative compendium.

---

### curated-motifs-finemo

```yaml
name: curated-motifs-finemo
type: other
source: Hand-curated from MoDISCo CWMs in this repo
format: NPZ + Parquet + TSV, two versions
raw_path: reference/curated_motif_data_v2.tsv
processed_path: reference/curated_motif_data_for_finemo_v2.motifs.npz
status: validated
known_issues:
  - v1 (2025-08-22) and v2 (2025-11-20) coexist. v2 is current — 8 active motifs after excluding FOS_JUN and the redundant GATA_TAL1_8BP.
  - Identity of REPEAT_G, NF2L_NFE, and ELF still unconfirmed (see HANDOVER_finemo_tf_analysis.md).
access_restrictions: none
tags: [motifs, finemo, curated, modisco-derived]
```

The curated motif set FiNeMo calls hits against. **Use the `_v2` files** — v1 is retained
only for reproducing older figures.

---

### k562-tf-chromatin-annotations

```yaml
name: k562-tf-chromatin-annotations
type: genomic
source: Sibling project — $OAK/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/
format: TSV
raw_path: $OAK/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv
status: processed
known_issues:
  - Lives outside this repo. Changes there silently affect p300+/p300- calls here.
access_restrictions: none
tags: [annotations, chromatin, p300-status, external-dependency]
```

Supplies the p300+/p300- status used throughout the FIMO and FiNeMo analyses. **External
dependency** — not versioned with this repo.

---

### hg38-cv-folds

```yaml
name: hg38-cv-folds
type: other
source: Standard hg38 chromosome-based 5-fold split
format: JSON
raw_path: reference/hg38_five_folds.json
metadata_path: reference/fold0.json … fold4.json
status: validated
access_restrictions: none
tags: [cross-validation, hg38, splits, reproducibility]
```

Chromosome-holdout fold definitions shared across all models, so CV numbers are comparable
between v1, v2, v3, multimodal, and the GM12878 models.

---

### gc-matched-negatives

```yaml
name: gc-matched-negatives
type: genomic
source: Genome-wide GC-content-binned background regions, hg38
format: BED
size: 89 MB (each of two versions)
raw_path: reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
status: processed
known_issues:
  - reference/hg38_gc_stride_1000_flank_size_1057.gc.bed (2025-03-25) is byte-identical in size to the 2025-10-31 genomewide_ file; confirm which one training actually reads before reusing.
access_restrictions: none
tags: [negatives, gc-matched, hg38, training]
```

Negative set for peak-vs-background training. Stride 1000, flank 1057.

---

### gata-chipseq

```yaml
name: gata-chipseq
type: genomic
source: GATA1 and GATA2 ChIP-seq in K562, downloaded via each project's download.sh
format: BigWig (stranded plus/minus) + peaks
raw_path: [K562_GATA1_BPNet/data/, K562_GATA2_BPNet/data/]
status: processed
access_restrictions: none
tags: [gata1, gata2, chipseq, k562, comparison-models]
```

Training signal and peaks for the GATA1/GATA2 comparison models — each directory holds
`plus.bigWig`, `minus.bigWig`, and `peaks_inliers.bed.gz` (GATA2 also has
`peaks_inliers.narrowPeak`). Provenance is in each project's `download.sh`.

---

### k562-h3k27ac-chipseq

```yaml
name: k562-h3k27ac-chipseq
type: genomic
source: ENCODE ENCSR000AKP — H3K27ac ChIP-seq in K562, 2 replicates
format: BAM (36 bp single-end, filtered/sorted/indexed) + derived coverage BigWig
size: 267 MB + 518 MB (BAMs), 164 MB (BigWig)
raw_path: $OAK/Users/sheth/Data/ENCODE/K562/ENCFF790GFL.se.filtered.sorted.bam
processed_path: $OAK/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw
status: processed
known_issues:
  - The BigWig lives under Data/share/IGV/, a visualization share directory, but is raw counts with 250 bp fragment extension (bam_to_bigWig.sh -r SINGLE) and is therefore valid as a training target. Do not assume files in that directory are normalized for display.
  - Unstranded, so n_outputs=1; the shared trainer currently requires both plus and minus BigWigs.
  - Reads are 36 bp with fragment extension to a fixed 250 bp. Extension length is an untested parameter.
access_restrictions: none
tags: [h3k27ac, chipseq, k562, encode, histone, training-target]
```

Second replicate is `ENCFF817HMW.se.filtered.sorted.bam`. The two were merged into the
coverage BigWig by `$OAK/Users/sheth/Data/scripts/bam_to_bigWig.sh`; the invocation is
recorded in `$OAK/Users/sheth/Data/ENCODE/log.sh`.

Per-replicate BigWigs were built 2026-08-24 by
`2026_0824_H3K27ac_model/scripts/0.2.make_replicate_bigwigs.sh` at
`2026_0824_H3K27ac_model/data/h3k27ac_rep{1,2}.bw` (47 MB / 120 MB), using settings
identical to the merged target so the ceiling is comparable. Results in
`results/replicate_ceiling_by_window.tsv`. Note the raw `r(rep1, rep2)` there is
single-replicate reliability, not a model ceiling — see `.living/learnings.md` for the
two corrections needed.
