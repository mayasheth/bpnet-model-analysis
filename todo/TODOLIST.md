# Todo List

Open work only. Completed items live in the git history and in `.living/`; the analysis
narrative is in `2026_0824_H3K27ac_model/h3k27ac_model_report.html`.

## Settled (do not redo)

Target = 5′ ends, ±500 bp window, `count_loss_weight = 10`. Panel is **ATAC-only** (DNase
and ATAC are not interchangeable inputs). PE H3K27ac targets use **read 1 only**.
Accessibility inputs should be ChromBPNet-style 5′ insertion counts; tracks are built and
validated but **no model uses them yet**. Residual-objective training helps only a
sequence-blind input and costs a multimodal one — replicated in K562 and GM12878, so it is
a property of the objective and needs no further per-cell-type testing.

## Highest value

- [ ] **W4. Motif syntax (SHAP / TF-MoDISco / FiNeMo).** Not started; the actual scientific
      goal. Run on the **residual-trained** model — its attributions are forced onto
      accessibility-independent signal, which multimodal attributions cannot separate. Use
      the ±500 bp window (zero neighbour contamination). Expect less signal than p300.

- [ ] **TeloHAEC training + transfer.** The only new cell type available under the ATAC-only
      rule. Tracks, elements and model-free coupling (0.33–0.37 top quintile) are all ready.
      3 modes × 5 folds, then the four-evaluation transfer set against K562 and GM12878.
- [ ] **TeloHAEC ±IL1b / ±TNFa / −VEGF.** Same genome and cell line, different regulatory
      state, no input-domain shift — a sharp and cheap test of whether the model tracks
      condition-specific change. Inference only if trained on ctrl.

## Open questions

- [ ] **Repeat the residual comparison for p300.** Resolved only for H3K27ac; p300's
      residual *r* is already 0.654 so the headroom may be smaller.
- [ ] **Predict the composite activity metric directly.** Downstream uses
      `geomean(accessibility, H3K27ac)`, best as `geomean(DHS, H3K27ac)`. Is the composite
      easier to predict than predicting H3K27ac and combining afterwards?
      **Design caution:** if the accessibility input is the same assay as the accessibility
      term in the target, half the target is readable off the input and the comparison is
      circular. Predict `geomean(DHS, H3K27ac)` from sequence + ATAC, and baseline against
      the two-step route scored on the *same* composite.
- [ ] **[Q for Maya] Is a mostly-accessibility model useful for the intended application,**
      or does the goal require the sequence component? Decides whether the architecture list
      below is worth the GPU time.

## Architecture, in expected order of value

- [ ] **Wider receptive field** (`--n-layers 10`, ~4.2 kb): H3K27ac is still at 28% of
      maximum at ±2 kb against a ~1.1 kb receptive field. ~30–60 GPU-h, the most expensive
      item — defer until motif work shows the sequence component justifies it.
- [ ] **Nucleosome-resolution binned profile target** (~50–150 bp). Code only. Touches
      `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared with p300, so it
      needs the same backward-compatibility regression test as the unstranded change.
- [ ] **Fragment-size ATAC channels** `[all, sub, mono, di]`. Sub and mono channels built;
      the two cover only 50% of fragments. A superset of the current flat input, so it
      cannot regress.
- [ ] **Switch the accessibility input to 5′ counts, as a set.** ~35 fold-jobs. Within a
      cell type this moves only the ATAC-only model; across cell types it removes a
      read-length confound (TeloHAEC 36 bp vs K562 95 bp). Mixing the two inputs is invalid.
- [ ] **Re-check ±500 vs ±1000** on the 5′ target once a retrain happens anyway.

## Panel data caveats

- [ ] **Element derivation is not uniform** (`reference/ELEMENT_DERIVATION.md`): K562 and
      GM12878 DNase-derived, TeloHAEC ATAC-derived. An ATAC-derived K562 set exists in
      `K562_ATAC_ChromBPNet/data/`; adopting it makes two of three consistent at the cost of
      re-deriving every K562 number.
- [ ] **HCT116** has ENCODE ATAC but its H3K27ac replicates differ in run type, so no
      inter-replicate ceiling is computable. Trainable, not normalisable.
- [ ] **H1, H9, Jurkat, THP-1 have no ENCODE ATAC at all** (all 559 released experiments
      checked). Reaching them means GEO/SRA fastqs through
      `Data/scripts/sra_paired_fastq_to_bam.sh` — a separate decision.
- [ ] **TeloHAEC_ctrl/ATAC holds 3 EA.hy926 files** (`SRR20809434/435/436`) under the same
      sample name. Always use explicit accession lists, never a directory glob.

## Statistical power

- [ ] **Re-run the `count_loss_weight` sweep with ≥3 folds and ≥2 seeds.** The current pick
      (10) came from single folds; 3/10/100 are within noise of each other.
- [ ] **Quantify run-to-run variance properly** — one config × 5 seeds.
- [ ] **Make submit scripts refuse to overwrite a completed fold directory.** A grid once
      silently overwrote a sweep result that shared an `OUT_DIR`.

## Housekeeping

- [ ] Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history is
      confirmed good.
- [ ] `Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a stray directory.
- [ ] `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      merging; works today only because inputs happen to be sorted.
