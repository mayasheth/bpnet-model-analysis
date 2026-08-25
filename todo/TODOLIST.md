# Todo List

Master list of future work items. Each item can have a detailed writeup
in a separate `.md` file in this directory.

## Items

### H3K27ac model — `h3k27ac-model.md`

Detailed writeup: [`h3k27ac-model.md`](h3k27ac-model.md). Ordered into waves by
dependency; everything inside a wave can run at the same time.

**Gate — run before committing GPU time to anything expensive**

- [ ] **G1. Inter-replicate ceiling on the 5′ target** — decides whether 5′ ends or the
      250 bp fragment track is the model target. Every later training decision inherits
      this. ~10 min CPU. `0.3.replicate_ceiling_by_window.py` with the new
      `h3k27ac_rep{1,2}_5p_{plus,minus}.bw`. The 5′ track is sparse (max 27, ~6 reads per
      kb at background), so the ceiling could fall; if it does, test unextended 36 bp
      read coverage as a middle option.

**Wave 1 — independent of G1, all parallel**

- [ ] **W1a. Residual / difference-from-ATAC model + evaluation** — the highest-value
      item and the one that matches the actual goal. Stage 1 already trained. ~4 GPU-h.
- [ ] **W1b. GM12878 cross-cell-type transfer** — separates model capacity from
      "H3K27ac is not sequence-determinable". Inference only on existing K562 models.
      ~0.5 GPU-h. Data all present (`ENCFF645BAL`, `ENCFF865OOP`).
- [ ] **W1c. Download paired-end ATAC BAMs → fragment-stratified channels** — unblocks
      the nucleosome-positioning input. `K562/log.sh:9` has the curl commands
      (commented). **Download to `$SCRATCH`, not Oak — Oak is at 97%.**
- [ ] **W1d. Implement nucleosome-resolution (binned) profile target** — code only, no
      compute. Touches `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared
      with the p300 models, so it needs the same backward-compatibility care and
      regression test as the unstranded change.

**Wave 2 — needs its Wave 1 / gate input**

- [ ] **W2a. Re-sweep `count_loss_weight` on the chosen target** (needs G1) — expect it
      to fall a long way from 1000 once MNLL sees real read counts. 2–3 GPU jobs.
- [ ] **W2b. Train with fragment-size ATAC channels** (needs W1c) — 5–10 GPU jobs.
- [ ] **W2c. Validate the binned profile target** (needs W1d) — 5–10 GPU jobs, ~4–8 GPU-h.
- [ ] **W2d. Re-check the ±500 vs ±1000 window on the 5′ target** (needs G1) — less
      smearing means less neighbour bleed, so the contamination penalty may differ.

**Wave 3 — only once the target and weight are settled, so it is not redone**

- [ ] **W3a. Full retrain grid on the final target + weight** — 15–30 GPU jobs.
- [ ] **W3b. Wider receptive field (`--n-layers 10`)** — best single architecture bet, and
      the most expensive: `in_window` 2114 → 5186, ~3.1× FLOPs/epoch, 2–4 h per job,
      ~30–35 GB RSS. 30–60 GPU-h. Deliberately last of the training items.

**Wave 4**

- [ ] **W4. Motif syntax (SHAP / MoDISCo / FiNeMo)** on the final chosen model. Per F-001,
      expect less sequence signal to attribute here than for p300.

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
