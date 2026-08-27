# Todo List

Master list of future work items. Each item can have a detailed writeup
in a separate `.md` file in this directory.

## Items

### H3K27ac model — `h3k27ac-model.md`

Detailed writeup: [`h3k27ac-model.md`](h3k27ac-model.md). Ordered into waves by
dependency; everything inside a wave can run at the same time.

**Gate — DONE 2026-08-25.** 5′ and fragment-extended ceilings are *identical* on the
top quintile (0.760 vs 0.761 at ±500). Switched to 5′ on principle — MNLL well-posed, no
neighbour bleed — not for an accuracy gain. See `.living/learnings.md`.

- [x] **G1. Inter-replicate ceiling on the 5′ target** — decides whether 5′ ends or the
      250 bp fragment track is the model target. Every later training decision inherits
      this. ~10 min CPU. `0.3.replicate_ceiling_by_window.py` with the new
      `h3k27ac_rep{1,2}_5p_{plus,minus}.bw`. The 5′ track is sparse (max 27, ~6 reads per
      kb at background), so the ceiling could fall; if it does, test unextended 36 bp
      read coverage as a middle option.

**Wave 1 — independent of G1, all parallel**

- [x] **W1a. Residual EVALUATION** — done (jobs 40799753, 40803766). Residual metric
      is now the headline; it reverses the p300/H3K27ac ranking (F-002). Residual
      *training* is now item W2e below, promoted by that result.
- [ ] **W1a-bis. Residual / difference-from-ATAC model + evaluation** — the highest-value
      item and the one that matches the actual goal. Stage 1 already trained. ~4 GPU-h.
- [ ] **W1b. GM12878 cross-cell-type transfer** — separates model capacity from
      "H3K27ac is not sequence-determinable". Inference only on existing K562 models.
      ~0.5 GPU-h. Data all present (`ENCFF645BAL`, `ENCFF865OOP`).
- [x] **W1c-i. PE ATAC BAMs downloaded** to `$SCRATCH/atac_pe` (24.5 GB, 165M/153M/229M
      paired reads). Channels building: job 40808663.
- [ ] **W1c. Download paired-end ATAC BAMs → fragment-stratified channels** — unblocks
      the nucleosome-positioning input. `K562/log.sh:9` has the curl commands
      (commented). **Download to `$SCRATCH`, not Oak — Oak is at 97%.**
- [ ] **W1d. Implement nucleosome-resolution (binned) profile target** — code only, no
      compute. Touches `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared
      with the p300 models, so it needs the same backward-compatibility care and
      regression test as the unstranded change.

**Wave 2 — needs its Wave 1 / gate input**

- [~] **W2a. Re-sweep `count_loss_weight` on the 5′ target** — RUNNING (jobs 40808673/6/8,
      weights 10/100/1000). Expect the optimum to fall far below 1000 now that MNLL
      sees real read counts. Was: (needs G1) — expect it
      to fall a long way from 1000 once MNLL sees real read counts. 2–3 GPU jobs.
- [ ] **W2b. Train with fragment-size ATAC channels** (needs W1c) — 5–10 GPU jobs.
      **Revisit the split first.** The fragment distribution is a clean nucleosomal
      ladder (sub-nucleosomal mode ~42 bp, mono ~205 bp, di ~400 bp, tri ~600 bp; see
      `figures/atac_fragment_lengths.pdf`), but the current two channels capture only
      30% (sub 1–99) and 20% (mono 180–247) of fragments — **50% fall in neither**,
      including both higher-order peaks and the 100–180 bp trough.
      Preferred design: **[all, sub, mono, di]** — keep an all-fragments channel so the
      input is a strict superset of the current flat-ATAC baseline and cannot regress,
      plus a di-nucleosomal channel (~350–450 bp) for higher-order structure. Coverage
      fractions per replicate are in `results/atac_fragment_length_summary.tsv`.
      Architecture support landed 2026-08-25: `MultiModalBPNet(n_acc_channels=N)`; the
      forward pass already slices `X[:, 4:, :]` generically, so only the `acc_conv`
      in-channels needed changing. Still to do: multi-BigWig accessibility input in
      `extract_windows` / `--accessibility-bw`, which currently accepts one track.
- [ ] **W2c. Validate the binned profile target** (needs W1d) — 5–10 GPU jobs, ~4–8 GPU-h.
- [ ] **W2d. Re-check the ±500 vs ±1000 window on the 5′ target** (needs G1) — less
      smearing means less neighbour bleed, so the contamination penalty may differ.

- [ ] **W2e. Residual TRAINING — promoted, high value.** An independently-trained
      sequence model captures almost none of the ATAC residual (r = 0.100): it is
      *redundant* with accessibility, not complementary. So do not expect a sequence
      model to find the complement on its own — train it explicitly on
      `observed − atac_pred`. Needs a per-region count offset in the trainer, applied to
      the count loss only. ~4 GPU-h once implemented.
- [ ] **W2f. GM12878 transfer evaluation** — GM12878 fragment-extended target building
      (job 40808659); ATAC (`2026_0606_GM12878_transferability/data/atac.bw`) and
      elements already exist. Must match the target definition the models were trained
      on, so this runs against the fragment target until the 5′ retrain lands.

**Wave 3 — only once the target and weight are settled, so it is not redone**

- [ ] **W3a. Full retrain grid on the final target + weight** — 15–30 GPU jobs.
- [ ] **W3b. Wider receptive field (`--n-layers 10`)** — best single architecture bet, and
      the most expensive: `in_window` 2114 → 5186, ~3.1× FLOPs/epoch, 2–4 h per job,
      ~30–35 GB RSS. 30–60 GPU-h. Deliberately last of the training items.

**Wave 4**

- [ ] **W4. Motif syntax (SHAP / MoDISCo / FiNeMo)** on the final chosen model. Per F-001,
      expect less sequence signal to attribute here than for p300.

### Multi-cell-type generalization (promoted — supersedes the single GM12878 pair)

Inventory built from `2023_0701 Maya's sequences.xlsx` ("ENCODE data" sheet) into
`2026_0824_H3K27ac_model/results/celltype_inventory.tsv`, under two standing rules:

1. **Mint-ChIP is excluded** — different assay chemistry, not comparable to standard ChIP.
   This drops **WTC11** from the H3K27ac panel entirely (ENCSR146DPQ is its only H3K27ac
   experiment), and also its H3K27me3 and H3K4me1.
2. **Each (experiment, processing) pair is a separate SAMPLE.** Replicates may only be
   pooled, or compared as an inter-replicate ceiling, within one experiment and one
   processing. Two consequences below.

After exclusions: **11 biosamples with H3K27ac**, of which 6 have ATAC (K562, GM12878, and
the 4 TeloHAEC conditions) and 7 have DNase (K562, GM12878, H1, H9, HCT116, Jurkat, THP-1).
The structure matters more than the count.

- [~] **Model-free ATAC-vs-H3K27ac correlation across the panel — cheap, do first.**
      Correlate summed ATAC against summed H3K27ac over the same element windows, no model
      involved. This is the model-free reference for the ATAC-only baseline: it separates
      "how much accessibility can explain in principle" from "how well a network learns
      it". If accessibility explains less of H3K27ac in other cell types, then the
      ATAC-only model's strength in K562 is partly a property of K562 rather than of the
      assay pair — which would also mean the sequence margin is not directly comparable
      across cell types. Needs no GPU and no training; it is bigwig reads plus a
      correlation. `scripts/0.12.atac_vs_h3k27ac.py`, accumulating into
      `results/atac_vs_h3k27ac_by_celltype.tsv` one row-set per label.
      **K562 and GM12878 reference points running now (job 41010962).** Remaining cell
      types are blocked only on candidate elements (below).
- [x] **UNBLOCKED: candidate element sets built for all 11 panel cell types.** Converted
      from ENCODE_rE2G `Neighborhoods/EnhancerList.bed` by
      `scripts/0.13.make_candidate_elements.py` into
      `reference/celltype_elements/<label>_<assay>_candidate_elements.narrowPeak`, each with
      a `.provenance.txt` sidecar. EnhancerList.bed is BED4 with no summit column, so the
      conversion writes `summit = width // 2` to preserve element-centred windows.
      Verified: EnhancerList.bed is byte-identical between the `*_H3K27ac_megamap` and
      `*_megamap` variants, so that choice does not matter.
      Sources: `2025_0227_validation_new_inputs` (K562, GM12878, HCT116, Jurkat),
      `2025_0219_ESC` (H1, H9), `Projects/E2G/THP1/ENCODE_rE2G_results` (THP-1),
      `Projects/E2G/endothelial_cells/.../ENCODE_rE2G_results` (TeloHAEC x4).
      WTC11 is present in both rE2G dirs but stays excluded on the Mint-ChIP rule.

- [ ] **Element derivation is NOT uniform across the panel — manage as a confound.**
      DNase/DHS-derived (`dhs_*` model dirs): K562, GM12878, HCT116, Jurkat, H1, H9, THP-1.
      **ATAC-derived** (`atac_h3k27ac_powerlaw`): all four TeloHAEC conditions. So TeloHAEC
      differs from the rest in how its elements were called, not only in cell type. Either
      restrict cross-cell-type claims to the DNase-derived subset, or derive DNase-based
      TeloHAEC elements, or state the confound wherever TeloHAEC is compared to the others.
      Also note H9's elements are much wider (mean 798 bp vs ~570-600 for the rest), which
      will change its neighbour-contamination profile at any given window.

- [ ] **Decide whether to re-derive the K562 analysis on its rE2G EnhancerList.** GM12878's
      in-use element set IS its EnhancerList (154,224 / 154,224 exact). K562's is not — it
      is a separate derivation. The two K562 sets are largely the same loci with different
      boundaries: of 150,528 in-use elements, 144,194 (96%) overlap an EnhancerList element
      and 142,536 (95%) at 50% reciprocal overlap, with only 6,334 (4%) having no overlap.
      So this is a boundary difference, not a different set of loci — but it does mean K562
      is currently the odd one out in the panel, and every window is centred slightly
      differently there. Cheap to resolve by re-running the K562 evaluations on
      `reference/celltype_elements/K562_DNase_candidate_elements.narrowPeak`.
      (An earlier exact-interval count of 25,550 shared was an artifact of comparing files
      that were not lexicographically sorted; the bedtools overlap figures above supersede
      it.)
- [ ] **Train a DNase-input model — this is the enabling step, not an optional variant.**
      Only 3 distinct cell types have ATAC (K562, GM12878, TeloHAEC — the other 3 ATAC
      entries are TeloHAEC conditions); **7 have DNase** (K562, GM12878, H1, H9, HCT116,
      Jurkat, THP-1 — WTC11 is out on the Mint-ChIP rule). Applying an ATAC-trained model to a
      DNase cell type shifts the *input* domain, which confounds the transfer test with a
      change of assay. So DNase is the common denominator across the panel, and without a
      DNase model most of this data is unusable for generalization. Note the p300 project
      already has a `models/dnase` directory that was never trained (blocked on generating
      `data/dnase.bw`), so this was always intended.
- [ ] **Cross-cell-type matrix: train on each, evaluate on all others.** Currently every
      conclusion assumes K562 is representative, and it may well not be — it is a cancer
      line with an atypical chromatin landscape. If sequence transfers better *out of*
      GM12878 or an untransformed line than out of K562, the story changes from "H3K27ac
      has little sequence signal" to "K562 is a poor training cell type". This is the test
      that decides which.
- [x] **Checked: the existing K562 and GM12878 ceilings are unaffected by rule 2.**
      ENCSR000AKP and ENCSR000AKC are each a single experiment, `run_type = se`, with two
      replicates processed identically (`.filtered.sorted.bam`), so every ceiling and
      transfer number reported so far is computed within one experiment and one processing.
      No correction needed.
- [ ] **TeloHAEC ±IL1b / ±TNFa / ±VEGF is a different and sharper test.** Those four rows
      are *conditions of one cell type*, not four cell types, so they are weak evidence
      for cross-cell-type transfer — but they are strong evidence for something else:
      whether the model tracks condition-specific H3K27ac changes *within* a cell type,
      with no input-domain shift and no element-definition change. Worth running early
      because it is cheap and the interpretation is clean.
- [ ] **Split the TeloHAEC rows before use — each pools 2 experiments.** All four
      TeloHAEC rows list `GSE210489,GSE210491`, so under rule 2 each row is at least two
      samples, not one. The ctrl row has 4 replicate files (SRR20810532/533/544/545),
      likely 2 per experiment, but the sheet does not say which SRR belongs to which GSE —
      that needs SRA metadata before any ceiling is computed. Computing a single
      4-replicate ceiling across both experiments would fold batch effect into what is
      reported as biological reproducibility. Upside: once split, the same-condition
      cross-experiment pair is a clean estimate of experiment-level batch effect.
- [ ] **HCT116 cannot yield an inter-replicate ceiling as catalogued.** ENCSR661KMA is
      `run_type = "se, pe"`, so its two replicates differ in run type and are separate
      samples under rule 2. Their correlation would conflate a run-type difference with
      biological noise and understate the ceiling. Either use one run type only (and then
      there is a single replicate, so no ceiling), or find additional same-processing
      replicates. Same caution applies to HCT116's DNase, CTCF, H3K27me3 and H3K4me1 rows,
      which are all `"se, pe"`.
- [ ] **Watch two confounds when comparing absolute numbers across the panel.** Run type
      is mixed (K562 is 36 bp SE; THP-1, TeloHAEC, WTC11 are PE), and provenance is mixed
      (ENCODE versus GEO/SRA with different processing, some deduplicated). Within-cell-type
      ceilings are unaffected; cross-cell-type absolute correlations are not, which is
      another reason to report every transfer number as a fraction of that cell type's own
      inter-replicate ceiling rather than raw.

### Revisit with proper statistical power

- [ ] **Re-run the `count_loss_weight` sweep with >=3 folds and >=2 seeds.** The current
      choice of 10 rests on single-fold runs, and measured run-to-run variance is ~0.018
      (same config gave 0.4962 and 0.4785 on separate runs). Only the extremes are
      resolvable at that noise level: 10 clearly beats 1 and 1000, but 3 / 10 / 100 are
      indistinguishable. Weight 10 sits in a flat region so nothing is currently wrong —
      revisit if the weight is ever suspected of mattering, or after the profile target is
      binned (which changes MNLL's magnitude and so moves the optimum).
- [ ] **Quantify run-to-run variance properly.** One config x 5 seeds would give an error
      bar for every comparison in this project. Several conclusions so far rest on
      differences of 0.01-0.05 between single runs, which may or may not clear it.
- [ ] **Make submit scripts refuse to overwrite a completed fold directory** (or add a run
      tag to the path). The 5-prime grid silently destroyed the weight sweep's clw10 run
      because both wrote to `models/sequence5p_hw500_clw10/fold0`.

### Housekeeping

- [ ] Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history
      is confirmed good.
- [ ] `$OAK/Users/sheth/Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a
      **1.3 GB leftover** intermediate from `bam_to_bigWig.sh` (its cleanup did not run).
      Oak is at 97%, so worth removing — owner's call, not deleted automatically.
- [ ] `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      `bedGraphToBigWig`, which requires chrom+start order. It has worked so far because
      `genomecov` follows BAM header order, but that is not guaranteed to match
      `chrom.sizes`. `0.5.make_5prime_bigwigs.sh` sorts explicitly; consider backporting.

### Report (do not build yet)

- [ ] **Write up this analysis path** via `/engreitzlab-report` → `/analysis-report`
      (Quarto → self-contained HTML), once the story settles. Deliberately deferred: the
      narrative has already reversed twice (flanking-window idea rejected; residual metric
      flipping the p300/H3K27ac ranking), so wait for the residual-training and transfer
      results before fixing a storyline. Figures that would carry it, already generated:
      `h3k27ac_signal_vs_distance`, `profile_comparison`, `h3k27ac_model_comparison`,
      `h3k27ac_stratified_eval`, `residual_evaluation`, `h3k27ac_replicate_ceiling`.
