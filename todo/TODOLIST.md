# Todo List

Open work only. Completed items live in the git history and in `.living/`; the analysis
narrative is in `2026_0824_H3K27ac_model/h3k27ac_model_report.html`.

## Settled (do not redo)

Target = 5' ends, ±500 bp window, `count_loss_weight = 10`. Panel is **ATAC-only** (DNase
and ATAC are not interchangeable inputs). PE H3K27ac targets use **read 1 only**.
Accessibility inputs should be ChromBPNet-style 5' insertion counts; tracks are built and
validated but **no model uses them yet**. Residual-objective training helps only a
sequence-blind input and costs a multimodal one, replicated in K562 and GM12878, so it is
a property of the objective and needs no further per-cell-type testing.

## Highest value

- [ ] **W4. Motif syntax (SHAP / TF-MoDISco / FiNeMo).** Not started; the actual scientific
      goal. Run on the **residual-trained** model, its attributions are forced onto
      accessibility-independent signal, which multimodal attributions cannot separate. Use
      the ±500 bp window (zero neighbour contamination). Expect less signal than p300.

- [ ] **TeloHAEC training + transfer.** The only new cell type available under the ATAC-only
      rule. Tracks, elements and model-free coupling (0.33-0.37 top quintile) are all ready.
      3 modes x 5 folds, then the four-evaluation transfer set against K562 and GM12878.
- [ ] **TeloHAEC ±IL1b / ±TNFa / −VEGF.** Same genome and cell line, different regulatory
      state, no input-domain shift, a sharp and cheap test of whether the model tracks
      condition-specific change. Inference only if trained on ctrl.

## Interrogating the ABC negative result, the main thread

The CRISPR benchmark says predicted H3K27ac adds nothing detectable over ATAC alone (best
predicted arm 0.482 vs floor 0.457, CI overlapping; deployment-scenario arms at or below the
floor). `4.9.diagnose_abc_gap.py` localised why: the prediction is accurate (Spearman 0.82 on
CRISPR-tested regions, better than genome-wide) but **ranks the functional elements too low**,
and since ABC's qnorm removes scale by construction, rank is the only channel available. All
arms catch the same positives, so the deficit is in suppressing negatives.

`4.10`/`4.12` then localised the error to a specific population: every model that sees ATAC
over-predicts H3K27ac at accessible-but-unacetylated elements by 7-8x its own median, while
sequence-only elevates them 1.6x and observed H3K27ac not at all. Those elements are GC 0.59,
CpG o/e 0.55, 2.4x promoter-enriched, ATAC/K27ac ratio 13.25 vs 2.20 typical, a CpG-island /
CTCF phenotype. **The signal is in sequence; the additive trunk gives it no way to veto
accessibility.**

- [x] **Gate x loss factorial. CLOSED 2026-09-05** (`1.19`, decision 2026-09-05). All 15
      fold-jobs finished and were scored; this entry was left marked RUNNING by mistake and the
      marker was cleared 2026-09-14 after re-deriving the paired deltas from
      `gatefac_k562_per_fold.tsv` with `2.38`. In-cell K562, paired against the matched ungated
      baseline on top-quintile Pearson: GATE -0.0001 (*p*=0.99), ASYM +0.0001 (*p*=0.99),
      GATE_ASYM -0.0032 (*p*=0.51). Transferred to GM12878: GATE -0.0016 (*p*=0.82), ASYM
      -0.0078 (*p*=0.11), **GATE_ASYM -0.0128 (*p*=0.039), i.e. resolvably worse**. Both
      interventions shared one premise, that the sequence branch can identify the
      accessible-but-unacetylated elements, and that premise is measured false (1.55x against
      1.44x where truth separates the tails 1.00x against 31.3x). The flags and
      `test_asymmetric_loss.py` stay in the code so the negative result stays reproducible.
      **Not scored, and the one thing that would change the reading:** nobody ran the `4.12`
      fold-elevation comparison on the gate arms, so the closure rests on aggregate *r* rather
      than on the failure population the gate was built for. Reopen only if some other result
      revives the premise; on its own it is not worth two `4.1` prediction passes.
- [ ] **Indicator-channel control.** GC and CpG-density tracks as extra accessibility
      channels, no gate, no loss change. If hand-supplied class information does as well as a
      learned gate, prefer it, simpler and interpretable. Both are sequence-derived, so
      unlike fragment channels they cost nothing in transferability. Needs `0.27` to build the
      tracks first.
- [x] **CTCF/EP300 enrichment in the error strata, PEAK OVERLAPS DONE** (`4.11`).
      Over-predicted 5%: CTCF 36.9% vs 17.4% typical (2.1x, confirming the CTCF hypothesis),
      but EP300 16.4% vs 8.8%, ENRICHED, contradicting my prediction of depletion, and
      H3K27ac peaks 30.6% vs 15.8% despite low H3K27ac RPM. Under-predicted 1%: EP300 57.7%
      (6.6x), H3K4me1 86.8% (3.1x), H3K27me3 0.0%, CTCF 9.3% (depleted), canonical active
      enhancers the model misses. CTCF enrichment is NON-MONOTONIC (36.9% at 5%, 23.4% at 1%),
      unexplained; the most extreme over-predictions may be a different population.
- [x] **Quantitative signal, not just peak overlaps. DONE 2026-09-05** (`4.13`, F-007);
      marker cleared 2026-09-14. Per-element RPKM from the local BAMs is in
      `error_strata_rpkm.tsv` and it changed the reading: the over-predicted 5% carries CTCF
      4.69 against 1.79 typical but H3K27ac 1.17 at 0.94x typical with **IgG 0.39 below
      background 0.48**, so the H3K27ac peaks called there at 2x background are threshold
      artefacts and the stratum really is unacetylated. The same control demoted the
      under-predicted stratum's p300 enrichment from 6.6x (peaks) to about 2.3x (signal minus
      input), which is what closed p300 as an auxiliary target (F-006). Peaks alone would have
      overstated both premises.
- [ ] **Characterise the 12 threshold-level false negatives** individually once `4.11` lands -
      observed H3K27ac 39x the genome median where the model predicts near-background.

## Sequence-gated multimodal for p300, diagnose before training

Naming: the architecture is **sequence-gated multimodal**; `gate` is the flag and filename
token. The asymmetric loss is named separately on purpose, because the factorial exists to
attribute any gain to the architecture or to the loss.

- [ ] **Confirm p300 has the same failure mode BEFORE training p300 with the gate.** The
      premise is that p300 suffers the same accessibility domination, and it may not: p300's
      residual *r* is 0.654 against H3K27ac's ~0.55, so sequence adds MORE beyond accessibility
      there. That is consistent with the problem being milder, or with it being equally severe
      but more fixable, different expected gains, same number.
      Cheap test, reusing existing machinery: predict p300 on the ABC candidate regions with
      the existing `2026_0529_multimodal_p300_model` models via `4.1`, then run the `4.12`
      fold-elevation comparison against observed p300 (ENCSR000EGE,
      `ENCFF466WKF`/`ENCFF163FSR`). Confirmed if p300's sequence-only arm gets high-ATAC /
      low-p300 elements right while its multimodal arm does not.
- [ ] **Then train p300 with the gate**, 15 fold-jobs, only if the diagnostic confirms it.
      `multimodal_bpnet.py` is already shared, so no further code change is needed, the p300
      submit scripts just gain the two flags.

## Which model transfers best, do not assume the in-cell-type winner

- [x] **Transfer matrix. CLOSED 2026-09-05** (`2.24`, F-005). Finished and written up; this
      entry was left marked RUNNING by mistake and the marker was cleared 2026-09-14 after
      re-deriving every delta from `txmatrix_*_per_fold.tsv` with `2.38`, which reproduced
      F-005 exactly. Paired against narrow+flat transferred, top-quintile Pearson: wide
      +0.0114 (*p*=0.14) K562->GM12878 and +0.0053 (*p*=0.53) GM12878->K562; **fragment
      channels -0.0011 (*p*=0.86) and -0.0029 (*p*=0.75), exactly null in both directions**
      despite gaining +0.0055 (*p*=0.032) and +0.0156 (*p*=0.002) in-cell. Fragment-length
      structure is a property of a specific ATAC library, which no in-cell test can detect.
      Wide moves `overall_pearson` (+0.0079, *p*=0.009) and `residual_pearson` (+0.0495,
      *p*=0.0009) transferred while leaving top-quintile unresolved, so it survives transfer
      on the mechanistic metric only. A transferred multimodal beats the target's own
      ATAC-only floor by +0.0306 [+0.0130, +0.0482] (*p*=0.008) into K562 but only
      +0.0132 [-0.0141, +0.0405] (*p*=0.25) into GM12878.
- [ ] **TeloHAEC as a third cell type, with conditions.** The only new cell type under the
      ATAC-only rule. Tracks, elements and model-free coupling are ready. Caveats to carry:
      ATAC-derived elements (though derivation was shown not to matter, p = 0.83), 36 bp reads
      vs 95 bp, and shallower libraries, so read length and depth remain confounded even
      though element derivation does not. Fragment channels there would need PE BAMs.
- [ ] **Train on multiple cell types to optimise transferability.** The most promising route
      to a deployable model, since it directly optimises what deployment needs rather than
      in-cell-type fit. Do it after the transfer matrix says which architecture to carry.

## Using p300 without needing p300 at deployment

The application target has ATAC and nothing else, so **p300 can never be a model input**, it
is unavailable in the target cell type, whatever we hold for training. It can be a target, an
auxiliary task, or an annotation. Ordered so the cheap prerequisite gates the expensive work:

- [ ] **First: do p300-target models predict the elements the H3K27ac model misses?** The
      under-predicted stratum is 6.6x EP300-peak-enriched, 3.1x H3K4me1, CTCF-depleted, zero
      H3K27me3, canonical active enhancers the H3K27ac model treats as background. If the
      existing p300 models in `2026_0529_multimodal_p300_model` also miss them, neither idea
      below helps and both are dropped. Reuses `4.1` (predict p300 on the ABC regions) and
      `4.12` (fold-elevation vs observed p300, ENCSR000EGE).
- [ ] **Then: multi-head model**, shared trunk, two counts heads, predicting H3K27ac and
      p300 from ATAC + sequence with a weighted loss. Deployment-valid because p300 is needed
      only at TRAINING time; inference reads the H3K27ac head. Directly supervises the failing
      population rather than adding capacity. Needs `extract_windows` to carry a second signal
      track and the loss to sum two count terms, more invasive than the gate, and
      `multimodal_bpnet.py` is shared with p300, so it needs the `test_asymmetric_loss.py`
      treatment.
      Trainable in K562 and GM12878, the same two cell types all the transfer work uses.
- [ ] **Fallback only: stack predicted p300 as an input channel.** Information-theoretically
      redundant, a predicted p300 track is a deterministic function of ATAC + sequence, which
      the H3K27ac model already sees, so it can only act as an inductive bias, and it costs
      two training runs and two models at inference to get what multi-head gets in one. Try
      only if multi-head underperforms.

## p300 as the activity term instead of H3K27ac

- [ ] **ABC + CRISPR with p300 predictions, and with observed p300.** p300 is a coactivator,
      one step closer to the sequence-specified event, and is arguably what ABC's activity
      term is really proxying. Observed p300 gives both a ceiling and, importantly, a qnorm
      REFERENCE, so the no-reference problem dissolves: build a p300 qnorm reference from the
      observed arm rather than running without qnorm and having to scale-match the geomean by
      hand. Data: EP300 ENCSR000EGE, `ENCFF466WKF/ENCFF163FSR.filtered.sorted.bam`, peaks
      `ENCFF702XPO.bed.gz`. p300 models already exist in `2026_0529_multimodal_p300_model`
      (residual r 0.654, higher than H3K27ac's).
- [ ] **Test removing qnorm** as its own arm (`use_qnorm: False`; `activity_base_no_qnorm` is
      already in every EnhancerList). Expectation is that it hurts, the model emits log1p
      counts over ±500 bp while observed H3K27ac is read counts over the element, so without
      qnorm the geomean multiplies incommensurate magnitudes, and ABC's thresholds are
      calibrated on qnorm'd values. Cheap enough to settle rather than argue.

## Downstream utility, decides whether the correlation metrics are the right target

- [ ] **Plug predicted H3K27ac into ABC and benchmark it.** The end-to-end test of whether
      this model is useful, and the concrete form of the standing question below. ABC's
      activity term is `geomean(accessibility, H3K27ac)`, so substitute the prediction for
      the observed mark, run ABC, and score against the CRISPR benchmark.
      **Three arms, or the result is uninterpretable:** observed H3K27ac (upper bound),
      predicted H3K27ac, accessibility alone (floor). If predicted lands near observed, the
      model is useful even at a top-quintile *r* of 0.69; if it lands near accessibility
      alone, it adds nothing downstream and the correlation gains we have been chasing do
      not matter.
      **Leakage trap:** ABC runs genome-wide but each fold model has seen four fifths of the
      genome. Assemble the genome-wide prediction from the 5 fold models, each applied only
      to its own held-out chromosomes, or the benchmark is contaminated.
      Repos: `~/Documents/ABC-Enhancer-Gene-Prediction`, CRISPR benchmark via
      `~/Documents/DC_TAP_Paper`.

- [x] **ATAC to DNase converter. BUILT AND BENCHMARKED 2026-09-12.** The cheap gate passed
      and the whole line of work ran. What it settled, and what it did not:

      **The assay is settled.** DNase beats ATAC as a model input, +0.037 in-cell K562 and
      +0.086 in-cell GM12878 top-quintile Pearson, and it is the only input change in the
      project that survives transfer, clearing GM12878's own ATAC-only floor by
      +0.056 [+0.014, +0.099] (F-010). As the ABC activity term it beats ATAC by
      **+0.0709 [+0.0476, +0.0938]** (F-012), a larger downstream effect than anything the
      modelling has produced. DNase also has a base-resolution profile where H3K27ac has
      none, 1 bp top-quintile ceiling 0.848 against 0.21.

      **The converter reproduces the profile and nothing downstream yet.** 1 bp top-quintile
      shape against observed DNase is 0.792 against a ceiling of 0.834 in K562 (95%), where
      the raw ATAC track manages 0.292; sequence does the work, since the ATAC-only arm
      reaches 0.563 (F-011). Transferred to GM12878 it reaches 0.490 against that cell type's
      0.684 ceiling (72%) and keeps 91% of its sequence margin. But in the CRISPR benchmark
      neither converter clears real ATAC resolvably (+0.021 [-0.005, +0.044] K562-trained,
      +0.014 [-0.012, +0.036] GM12878-trained) and real DNase still beats the better one by
      +0.050 [+0.030, +0.071]. Training cell type does not matter, unlike p300 (F-009).

      **That benchmark could not test the converter's actual claim.** ABC sums the activity
      track over each region and discards the profile, so F-012 scored the converter's counts,
      which are its weak half. The confound the original gate worried about was handled the
      other way round from what was planned: the converter trains on the ATAC-derived element
      set (D-19) so it never trains on regions chosen by the signal it predicts.

- [x] **H3K27ac under four accessibility inputs. DONE 2026-09-12, F-013.** Top-quintile
      Pearson: real DNase 0.728, DNase smoothed at 250 bp 0.704, ATAC 0.690, converted DNase
      0.680. Paired against ATAC: real DNase +0.0379 [+0.0213, +0.0546], smoothed
      +0.0141 [-0.0057, +0.0338] (not resolvable), **converted -0.0099 [-0.0175, -0.0023],
      i.e. resolvably worse than the track it would replace**.

      **The control answered: shape, not magnitude.** Destroying structure finer than 250 bp
      costs real DNase +0.0239 [+0.0027, +0.0451] of its +0.0379, and the remainder no longer
      separates from ATAC. So base-resolution structure is where the advantage lives, which is
      the branch that keeps a converter conceptually alive.

      **But this converter is worse than raw ATAC.** Reproducing the DNase profile at 95% of
      ceiling (F-011) while hurting downstream means the painted track carries a defect that
      outweighs its shape fidelity. Named candidate: the 1.54x magnitude inflation of the
      lowest observed-signal quintile (D-22), plus a 34% genome-wide total. Shape fidelity is
      necessary and demonstrably not sufficient.

      **Also measured:** `profile_pearson` is 0.063-0.064 across all four arms, so the
      H3K27ac profile head learns nothing whatever accessibility it is given.

- [x] **Magnitude fixes for the painted track. BOTH FAILED 2026-09-14, F-014.** Quantile
      mapping to another cell type's DNase distribution (`0.37`) matched the marginal exactly
      and made dynamic range WORSE, 18x to 9x against real DNase's 39x, costing 0.12 of 1 bp
      top-quintile shape. A power transform with gamma fitted out-of-cell (`0.39`, `0.40`)
      raised range to 24x for a shape cost of 0.006, and moved the downstream delta from
      -0.0099 to -0.0105, i.e. nothing. Magnitude was not what was costing the converter.

- [x] **The converter is finished as a model input. STOPPED 2026-09-14** (decision
      2026-09-14, F-014). Resolvably worse than plain ATAC in all four directions of the 2x2:
      -0.0105 in-cell K562, -0.0579 K562->GM12878, -0.0206 in-cell GM12878. Do not attempt a
      third transform. The machinery is kept (`4.20`, `4.22`, `0.36`, `0.37`, `0.39`, `0.40`)
      because a later attempt would reuse it and the shape-destroyed control generalises.

- [x] **Filter-width sweep for the shape control. SUPERSEDED 2026-09-14, F-014.** The 2x2
      answered the question the sweep was meant to tighten, and answered it differently than
      expected: whether DNase wins on shape or magnitude is CELL-TYPE DEPENDENT. Smoothing at
      250 bp costs real DNase +0.0239 (*p*=0.035) in K562 and +0.0677 (*p*=0.0002) transferred
      to GM12878, but **+0.0011 (*p*=0.72) in-cell GM12878**. Shape matters where DNase is deep
      enough to have reliable shape: GM12878's DNase is 53.8M reads against K562's 301.1M, with
      a profile ceiling of 0.686 against 0.848. A width sweep within K562 would have tightened
      a claim that does not generalise; the informative axis is depth, not width.

- [ ] **DNase-input panel, third cell type: THP-1. TRACKS READY 2026-09-14.** Now the
      successor to the converter line, because both surviving DNase results need REAL DNase:
      DNase as an input (F-010, F-014) and DNase as an ABC activity term (+0.0709, F-012). The
      deployment story is "use DNase where it exists", not "synthesise it where it does not",
      and DNase exists in far more cell types than ATAC.

      `0.38` built THP-1's DNase 5-prime input, pooled and per-replicate stranded H3K27ac
      targets (read 1 only; both replicates are paired-end), DNase-derived candidate elements,
      and the H3K27ac inter-replicate ceiling: 0.974 raw, 0.993 Spearman-Brown corrected on all
      elements. DNase alignments are at
      `Users/sheth/Data/ENCODE/THP1/DNase/AG81591.filtered.bam`; the inventory previously said
      THP-1 had none and was not modellable, which was wrong and had ruled it out.

      **What THP-1 can and cannot do.** One DNase replicate, so no DNase shape ceiling and no
      use as a converter target. No ATAC, so it cannot train a converter and cannot join an
      ATAC-input comparison. It is a DNase-input panel member, which is what is now wanted.

      **SUBMITTED 2026-09-14** (`1.30`, jobs 43442172-84): `multimodal` and the
      accessibility-only floor, 5 folds each. `1.30` is `1.22` with exactly three
      substitutions, THP-1's DNase input, H3K27ac target and element set, so all three cell
      types share every hyperparameter and can sit in one matrix. Pre-flight: 183,778
      elements over 25 contigs, all present in both tracks, 2.3% all-zero H3K27ac windows,
      125,187 elements in fold0's training chromosomes.

      **Depth cuts both ways here, so neither cell type is the easier target.** THP-1's
      H3K27ac is DEEPER than K562's (mean 208.2 per 1 kb window against 51.1) while its DNase
      is SHALLOWER (281.4 against 910.2), being one replicate. First floor model in: THP-1
      accessibility-only fold0 reaches validation count *r* 0.772 where K562's ATAC-only
      equivalent reaches 0.836.

      **The panel needed two floors that did not exist.** Every accessibility-only model in
      the project used ATAC, and THP-1 has no ATAC, so the three cells would have had floors
      built on different assays in a panel whose whole point is the DNase input. K562 and
      GM12878 DNase-only arms submitted 2026-09-14 (`1.22 atac` / `1.23 atac`, jobs
      43445448-533) into `atac5p_dnase_hw500_clw10` and
      `gm12878_atac5p_dnase_hw500_clw10`. No new script was needed; both already took
      `MODE atac` and nobody had run it.

      Then score the three-way matrix: nine transfer cells plus each target's own DNase-only
      floor. The config has to be built after the models exist, since `2.15` reads each
      checkpoint's geometry from `model.trimming`.

- [ ] **Sanity-check the GM12878->K562 collapse.** Every DNase-family input transferred in
      that direction lands far below ATAC (-0.133 real DNase, -0.297 smoothed), and the plain
      ATAC-input model is the best transferred model there. The in-cell GM12878 arms are strong
      (0.663 top-quintile), which argues against the GM12878 models simply being weak, but that
      is an inference rather than a test. Same asymmetry and same direction as p300 (F-009), so
      this is now the second finding blocked on the same gap.

- [x] **Multi-task arm: DNase profile head, H3K27ac counts head. DONE 2026-09-14, F-015
      established.** Kept, but it is a small architectural win and not progress on the real
      gap. Against an **epoch-matched** baseline: overall Pearson **+0.0026 [+0.0006,
      +0.0046]** (*p*=0.022), `incremental_r2` 0.005 [0.002, 0.007] and `incremental_r2_topq`
      0.015 [0.001, 0.030] both clearing zero. **Top-quintile Pearson does NOT clear against
      that control: +0.0076 [-0.0078, +0.0230], *p*=0.24** (0.701 against 0.693 and the
      early-stopped baseline's 0.690). The +0.0109 (*p*=0.043) first reported was against the
      early-stopped baseline and about a third of it was the epoch difference.

      **The auxiliary task itself is decisively learnable:** validation profile *r*
      0.530-0.561 in all five folds against 0.060-0.068, and the head drops H3K27ac shape to
      0.004. Needs no DNase at inference.

      **The epoch confound was real and is excluded** (`1.29`, jobs 43440391-96). Given the
      same 100-epoch budget the baseline moves +0.0000 [-0.0008, +0.0009] (*p*=0.92) and its
      best checkpoint still lands at epoch 33-47. Longer training was a CONSEQUENCE of the
      different loss, not the cause of the gain. **Carry this forward as a method rule: any
      intervention that changes the loss also changes the stopping rule, so an epoch-matched
      arm is the only way to attribute the result.** Three earlier arms in this project
      changed the loss (`1.19` asym, the residual objective, every `clw` sweep point) and
      none of them had one.

      **Only open follow-up, and it is optional:** `--profile-loss-weight 0.0561` is the
      measured depth ratio and the smallest defensible weight, so the ceiling on this effect
      is unknown. A sweep upward is five fold-jobs per weight. Not worth it unless something
      else revives interest, since even a doubling of the effect stays an order of magnitude
      below swapping the accessibility INPUT to real DNase (+0.038, F-010).

      Original scoping, for the record: Same trunk, same counts objective, profile head's target swapped from
      H3K27ac to DNase. Trainer change is in (`44e0e3c`): `extract_windows` takes an optional
      second signal pair, the dataset carries it as a fourth/fifth slot, and the loss takes
      `y_profile`. Gated by `scripts/test_profile_target.py`;
      `scripts/test_asymmetric_loss.py` still passes unchanged, and all ten existing callers
      of `extract_windows` keep the 4-tuple. Submit with `1.28`, score with `2.39`
      (`config/multitask_k562_configs.json`), paired against
      `multimodal5p_accs5p_hw500_clw10`.

      **`--profile-loss-weight` was NOT in the original plan and the arm is uninterpretable
      without it.** MNLL is `-sum(y * log_softmax)` plus a y-only term, so it scales with the
      target's read depth, and `loss = profile + clw * count`. K562 DNase carries **17.8x**
      the reads of K562 H3K27ac over these 1 kb windows (mean per-window total 910.2 against
      51.1). Swapping the target in unweighted would have grown the profile term about 18x
      and divided the effective count weight by the same factor, so a loss on counts would
      have measured "we trained it to care less about counts" rather than anything about the
      auxiliary task. The arm runs at the measured ratio, 0.0561, which starts the profile
      term at the baseline's magnitude. Default stays 1.0, so no earlier run moves.

      **Read only the counts columns.** `profile_*` in `2.39`'s table scores the profile head
      against H3K27ac for both arms, so for this arm it measures how well a DNase-trained head
      happens to predict H3K27ac shape. A drop there is expected and says nothing.

      **If it comes out null, the weight is the first thing to question, not the last.** 0.0561
      matches gradient magnitude to the baseline, which is the clean isolation of the change,
      but it is also the smallest defensible weight. A sweep upward is the follow-up, and it
      is cheap: no new tracks, no painting, five fold-jobs per weight.

      **Its premise held up and is now measured twice.** `profile_pearson` sits at 0.063-0.064
      across ALL FOUR accessibility inputs in K562 (F-013), so the H3K27ac profile head learns
      nothing whatever it is fed. It is dead weight with a target whose 1 bp inter-replicate
      ceiling is 0.21.

      **Run it in K562, not GM12878, and F-014 is the reason.** A DNase profile head is only a
      learnable auxiliary task where DNase has reliable base-resolution shape. K562's DNase
      profile ceiling is 0.848; GM12878's is 0.686 on a 5.6x shallower library, and F-014
      showed that in GM12878 destroying all sub-250 bp DNase structure costs the model nothing
      (*p*=0.72). Training the auxiliary head against GM12878 DNase would be fitting the same
      noise the H3K27ac head already fits.

      **Note what it does NOT inherit from the converter's failure.** This uses DNase as a
      TARGET at training time and needs no DNase at inference, so it is unaffected by the
      painted track being unusable (F-014) and is deployment-legal by the same argument that
      makes a multi-head model preferable to feeding p300 in (DATA_INVENTORY, deployment
      constraint).

      **Why.** The H3K27ac profile head currently trains against a target whose 1 bp
      inter-replicate ceiling is 0.21, i.e. mostly noise. D-3 kept it anyway, down-weighted,
      and that entry's own consequences note flags that `fconv` still consumes capacity for a
      near-zero gradient. DNase's 1 bp ceiling is 0.848, so this gives the head a learnable
      task for the first time. Hypothesis: a learnable base-resolution task makes the shared
      trunk better at H3K27ac COUNTS, which is the quantity of interest.

      **Why not joint prediction of both signals.** A second counts objective competes for the
      same gradient, and the p300 auxiliary head is the precedent: dropped 2026-09-05 after
      the p300 models predicted 1.83x enrichment at H3K27ac's failure elements against a
      2.10x INPUT control, i.e. below the control. What makes DNase different is specific.
      p300 failed because it could not identify the right elements; DNase would supply
      base-resolution structure that H3K27ac simply does not have.

      **Why it waits.** If DNase's input advantage is magnitude rather than shape the case
      weakens, though not fatally: DNase as an auxiliary TARGET teaches the trunk a
      sequence-to-structure mapping, a different channel from DNase as an INPUT. The control
      informs the prior rather than deciding it.

      **Implementation.** Needs a trainer change, not just a config: the profile and counts
      heads currently read one target. Cheapest correct route is a second pair of signal
      bigwig arguments consumed only by the profile loss. No painting and no second model, so
      materially cheaper than the converter path.

- [ ] **Converter count accuracy.** Considered and DEFERRED 2026-09-12, not closed. The
      converter already predicts DNase counts at top-quintile r 0.878 and still bought nothing
      in ABC, and the part of DNase counts carrying the downstream benefit is by construction
      the part that differs from ATAC counts. Revisit only if the H3K27ac arms say counts are
      the binding constraint. If they do, quantile-map the painted track to a reference
      accessibility distribution first, which is deployment-legal and fixes the residual
      magnitude distortion.

      **Caveat to carry into reading those arms.** The painted track compresses dynamic range:
      painted/observed ratio by observed-DNase quintile runs 2.73, 0.88, 0.73, 0.71, 0.73.
      Q2 to Q5 is near-uniform and training's own normalisation absorbs it, but the 1.54x
      residual inflation of the lowest quintile after thresholding (D-22) is the named
      alternative explanation if the converted arm underperforms.

      This does not contradict the ATAC-only panel rule. That rule exists because ATAC and
      DNase are not interchangeable as *inputs*; a converter is the principled way to get a
      DNase-like input everywhere without mixing assays across cell types.
      Also unblocks the composite-metric item below, whose best form needs DHS.

## Adopt-or-not decisions, all measured

- [ ] **Turn on test-time RC averaging as the default?** Free (one extra forward pass, no
      retraining) and positive for every model: +0.0162 sequence only, +0.0081 sequence +
      ATAC, +0.0016 and non-significant for ATAC only, which is the control (report Fig. 13).
      Currently opt-in via `2.15 --rc-average` so existing numbers stay comparable. Adopting
      it means re-scoring the report's tables as a set.
- [x] **Do the wider receptive field and the fragment channels combine? ANSWERED, they are
      SUB-ADDITIVE, and it does not matter for deployment.** Closed 2026-09-14; the model
      (`multimodal5p_fragchan_wide_hw500_clw10`) was trained and scored earlier and the answer
      was sitting unread in `txmatrix_gm12878_to_k562_per_fold.tsv` and
      `rc_fragwide_k562_per_fold.tsv`. In-cell K562 top-quintile, paired: wide+frag over
      wide-alone **+0.0065 (*p*=0.053)** against fragment channels' +0.0156 on their own, and
      over frag-alone +0.0188 (*p*=0.015) against wide's +0.0280 on its own. So each addition
      keeps roughly half its solo gain and they do not sum. **Deployment is decided by transfer,
      not by this:** wide+frag over wide-alone transferred is +0.0045 (*p*=0.52), and F-005
      already showed fragment channels are null transferred in both directions, so the
      combination adds nothing a deployed model could use. Prefer wide+flat.

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
      or does the goal require the sequence component? The ABC benchmark above is the
      empirical form of this question and would settle it without needing an answer in
      advance. Also decides whether the +0.043 all-elements receptive-field gain matters:
      it is real, but it is entirely dead-vs-active separation, so it pays off only for an
      on/off task.

## Architecture, in expected order of value

- [~] **Wider receptive field, ADOPT for multimodal; the earlier "do not pursue" was wrong.**
      `n_layers` 10 (~4.2 kb vs ~1.1 kb), 5 folds, both cell types, paired within fold on the
      intersection of valid regions.
      Top quintile: sequence-only −0.006 (K562) and +0.010 (GM12878), both p=0.53;
      **multimodal +0.027 (p=0.006) and +0.014 (p=0.025)**, all five folds rising in both.
      Accessibility residual 0.502 -> 0.547 and 0.397 -> 0.469.
      Results in `wide_{k562,gm12878}_*.tsv`; report Fig. 10.
      - [x] **ATAC-only wide arm. DONE, and it reproduces the multimodal gain.** Closed
            2026-09-14; the 10 fold-jobs finished and were scored into
            `rc_wide_{k562,gm12878}_per_fold.tsv`, where the marker was left stale. Top-quintile,
            paired within fold: ATAC-only widening gains **+0.0168 [+0.0048, +0.0288]**
            (*p*=0.018) in K562 and **+0.0144 [+0.0071, +0.0217]** (*p*=0.005) in GM12878,
            against the multimodal gains of +0.0280 and +0.0162 in the same tables. In GM12878
            the ATAC-only arm captures essentially the whole effect, so the extra context is
            being read as accessibility neighbourhood rather than as sequence context; K562
            leaves room for a sequence contribution (+0.0280 against +0.0168) but the intervals
            overlap. Conclusion as pre-registered: widen the deployed model.
            **Caveat:** both arms come from the RC-averaged tables, which is internally
            consistent since every arm in those tables is RC-averaged, but the numbers are not
            directly comparable to the non-RC `wide_*` tables.
      - [ ] Decide whether to re-run the transfer and deployment comparisons at
            `n_layers` 10, since those used the narrow multimodal model.
- [ ] **Decide what to do with the profile head, one decision, three options.** The 1 bp
      profile task is close to unlearnable: measured inter-replicate ceiling is 0.21 (K562)
      and 0.18 (GM12878) on the top quintile, rising to 0.72 and 0.70 at 50 bp binning
      (`profile_ceiling_binsize_*.tsv`, report Fig. 11). The options are (a) drop the profile
      loss entirely, (b) keep 1 bp, (c) bin to ~50 bp.

      **Test by hyperparameter before touching the architecture.** `loss = profile_loss +
      w · count_loss`, so a large `w` already approximates dropping the profile term, and the
      recorded sweep points AGAINST dropping it: `w` = 10 gave 0.496 while 100 gave 0.467 and
      1000 gave 0.464 (header of `1.11.submit_training_5prime_accs5p.sh`). Down-weighting the
      profile loss made the COUNTS worse. A head can be a poor predictor and still be a
      useful auxiliary task on the shared trunk, and those numbers are single-fold, so this
      is weak evidence pointing the opposite way from removal rather than a settled answer.
      The properly powered `count_loss_weight` sweep already listed under Statistical power
      answers it at zero architecture risk, run it jointly with bin size, since the two
      interact: binning changes how learnable the profile term is and therefore its optimal
      weight.

      **Gate on the real numbers.** Everything above about what the head achieves rests on a
      256-element smoke test (`profile_pearson` ~ 0.053 against a 0.21 ceiling). Job for the
      full 5-fold measurement is queued; use `prof_residual_grid_*.tsv` when it lands.

      **Cost of actual removal** (option a): it changes the model class, so every existing
      checkpoint becomes non-comparable and the whole grid needs retraining. It touches
      `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared with p300, so it
      needs a backward-compatibility regression, use the `2.18` pattern, exact on region
      counts and tolerant on metrics.
      **Payoff if it holds:** removes a term that is largely fitting Poisson noise, frees
      trunk capacity, and roughly halves the output tensor, buying a larger batch or a wider
      window at the same memory.
- [x] **Fragment-size ATAC channels, DONE and positive.** Top quintile 0.690 -> 0.703,
      paired +0.0135 [+0.0075, +0.0195], p = 0.0034, every fold rising; report Fig. 12.
      Strict superset of the flat input (bins sum to it exactly), so the gain is added
      information. Results in `fragchan_k562_per_fold.tsv`.
      - [ ] **Do the two accessibility-side gains combine?** `n_layers` 10 x fragment
            channels. Both may be reading the same neighbourhood structure, so +0.027 and
            +0.0135 may not sum. 5 fold-jobs, and it decides what the deployed model is.
      - [ ] **GM12878 replication** needs its paired-end BAMs downloaded, fragment length
            lives in TLEN, which the per-read tagAligns discard.
      - [ ] **ATAC-only fragment arm** (`1.15 atac`), asks whether fragment structure alone
            beats flat ATAC, with no sequence involved.
- [ ] **Switch the accessibility input to 5' counts, as a set.** ~35 fold-jobs. Within a
      cell type this moves only the ATAC-only model; across cell types it removes a
      read-length confound (TeloHAEC 36 bp vs K562 95 bp). Mixing the two inputs is invalid.
- [ ] **Re-check ±500 vs ±1000** on the 5' target once a retrain happens anyway.
- [ ] **Train on the ATAC-derived K562 element set**, easy (one path change in `1.11`),
      and now load-bearing for the ABC work rather than merely tidy.
      Verified path, 153,545 regions:
      `ENCODE_rE2G/results/2025_0226_ATAC_powerlaw_models/ATAC_H3K27ac_powerlaw/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed`
      (element set with class labels: the sibling `Neighborhoods/EnhancerList.bed`).
      **This file is byte-identical (md5 7d5995ce…) to the candidate regions of the July ABC
      run**, so the ABC experiment scores models on the ATAC-derived set while they were
      trained on the DNase-derived one (150,528 elements). The windows are convolutional so
      it is not fatal, but it is a genuine train/score element mismatch inside the ABC
      comparison, and training on this set removes it.
      Also the cheapest test of whether element derivation matters at all; if it does not,
      the panel-wide inconsistency stops being a caveat worth carrying.

## Panel data caveats

- [ ] **Element derivation is not uniform** (`reference/ELEMENT_DERIVATION.md`): K562 and
      GM12878 DNase-derived, TeloHAEC ATAC-derived. An ATAC-derived K562 set exists in
      `K562_ATAC_ChromBPNet/data/`; adopting it makes two of three consistent at the cost of
      re-deriving every K562 number.
- [ ] **HCT116** has ENCODE ATAC but its H3K27ac replicates differ in run type, so no
      inter-replicate ceiling is computable. Trainable, not normalisable.
- [ ] **H1, H9, Jurkat, THP-1 have no ENCODE ATAC at all** (all 559 released experiments
      checked). Reaching them means GEO/SRA fastqs through
      `Data/scripts/sra_paired_fastq_to_bam.sh`, a separate decision.
- [ ] **TeloHAEC_ctrl/ATAC holds 3 EA.hy926 files** (`SRR20809434/435/436`) under the same
      sample name. Always use explicit accession lists, never a directory glob.

## Statistical power

- [ ] **Re-run the `count_loss_weight` sweep with ≥3 folds and ≥2 seeds.** The current pick
      (10) came from single folds; 3/10/100 are within noise of each other.
- [ ] **Quantify run-to-run variance properly**, one config x 5 seeds.
- [ ] **Make submit scripts refuse to overwrite a completed fold directory.** A grid once
      silently overwrote a sweep result that shared an `OUT_DIR`.

## Housekeeping

- [ ] Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history is
      confirmed good.
- [ ] `Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a stray directory.
- [ ] `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      merging; works today only because inputs happen to be sorted.

## Carried over from the old root TODO.md (2026-09-12)

Open items that were stranded in `TODO.md` and `todo/h3k27ac-model.md`. Figure status and the
p300-model CV key values moved to `reference/PUBLICATION_FIGURES.md`; completed sections were
dropped, since git history and `.living/` hold them.

- [ ] **Revisit the p300+ definition across all figures.** `EP300_peak_overlap` (from
      `finemo_peaks_all_chr.chromatin_annotations.tsv`) flags any 1 kb window overlapping a
      p300 peak call, and **~40% of "p300+" elements by that definition have
      `true_logcounts = 0`**, most likely because the peak overlaps only the window edge,
      outside the 500 bp window where reads are counted. Consider replacing it with the top
      20% of elements by observed p300 counts. Figures that would change: `model_comparison`
      (1d), `k562_multimodal_comparison` (2 K562), `training_region_comparison` (S1b),
      `transferability_bar` (2d), every `*_p300plus.pdf`, and the `mean_predictions` scatters.
- [ ] **Two BioRender architecture schematics never started**: Fig 1a (BPNet) and Fig 2a
      (multimodal).
- [ ] **ATAC ChromBPNet true CV performance.** Fig 1d still carries the manuscript's
      Pearson r = 0.70 as a placeholder, with no p300+ number at all.
- [ ] **Section 3 and later publication figures.** MoDISCo motif logos, FiNeMo hits, and the
      motif spacing/pair experiments. Two pieces are already done: the individual-motif
      insertion violin plot (`scripts/plot_individual_motif_insertions.py`) and the motif-pair
      heatmaps for max log2FC and synergy (`scripts/plot_motif_pair_heatmaps.py`). The logos
      and hit panels are not.
- [ ] **Compare the GM12878 MoDISCo motifs to the K562 set.** Raised in the GM12878
      handover and never done. 26 GM12878 motifs exist
      (`modisco/max_seqlets_250k_30_10_0/logos/`). A shared motif vocabulary across the two
      cell types is the sequence-side counterpart to the transfer results.
- [ ] **Confirm the identity of `REPEAT_G`, `NF2L_NFE` and `ELF`** in the v2 FiNeMo motif set.
      Raised in the FiNeMo TF-analysis handover and unresolved.
- [ ] **hashFrag chromosome splits** to reduce train/test sequence-similarity leakage. Never
      assessed; every result in the project rests on plain chromosome holdout.
- [ ] **Update `P300_INTERACTORS`** in `plot_finemo_composite_figure.py` after a
      literature/BioGRID review. The composite figure's p300-interaction star depends on it.
- [ ] **SHAP on the K562 ATAC multimodal model.** Never submitted. Overlaps W4 above, but W4
      specifies the residual-trained model, so this is a separate run.
- [ ] **Reweight the accessible-but-unacetylated quadrant.** Those elements are already
      positives with near-zero targets, so up-weight rather than relabel. Aims directly at
      the over-predicted tail. Survivor of the training-composition family after GC-matching
      was refuted (D-18).
- [ ] **Sweep `negative_ratio` (0.1) and `--max-negatives` (50,000).** Both inherited from the
      p300 setup and never swept for H3K27ac.
- [ ] **Bin the profile target to nucleosome scale (~50-150 bp).** H3K27ac has no meaningful
      1 bp structure, so the profile head is fitting noise. Cheap to run and expensive to
      develop: touches `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, which the p300
      models share. **Read alongside the DNase multi-task arm above, which attacks the same
      dead profile head from the other direction**, by giving it a learnable target instead of
      coarsening an unlearnable one. Doing both would confound them.
- [ ] **Re-check the counting-window choice on the 5-prime target.** Less smearing means less
      neighbour bleed-through, so the contamination penalty at +/-1000 may differ from what
      the fragment-extended target showed.

## Where the project goes next (moved from TODO.md 2026-09-12; written 2026-09-07)

Every model so far fits one cell type's accessibility-to-signal mapping. Accessibility is the
part that generalises, so a transferred model decays toward the ATAC-only floor by
construction. The one arm that transferred well (K562-trained p300, +0.207 over the target's
floor, 83% of its local advantage) is also the one trained on 2.3x more reads in peaks, so
portable sequence grammar is learnable here and the models that failed to learn it were the
thinner ones. Order below reflects that.

### Running now

- [x] **Depth-subsampling causal test. REFUTED 2026-09-08**: the depth-matched K562 model
      still transfers at +0.202 against the full-depth +0.207 (paired -0.005, *p*=0.40), and
      peak-set geometry is excluded too. The asymmetry survives equalising volume, so a third
      cell type is now required rather than merely desirable. Original plan: Subsample K562 EP300 to GM12878's 30.0M mapped reads
      and 21,068 peaks, retrain the multimodal p300 model, re-run the transfer 2x2 (`2.27`).
      If the subsampled K562 model stops transferring, training-signal volume explains the
      asymmetry and the answer is more sequencing. If it still transfers, something else about
      K562 as a training cell type does, and a third cell type becomes urgent.
      Five folds, no new tooling beyond the subsampling and bigwig rebuild.

### Top section

- [ ] **Multi-cell-type joint training, leave-one-cell-type-out.** K562 + GM12878 + TeloHAEC,
      shared sequence trunk, per-cell-type accessibility input, held-out cell type as the
      evaluation. Measures the deployment case directly instead of inferring it from a two-way
      transfer. Needs a cell-type-invariant target scale, which the quantile item supplies.
- [ ] **Quantile-normalised targets within cell type, plus a ranking loss** in place of log1p
      MSE on counts. ABC consumes an ordering; the measured downstream failure is compressed
      dynamic range on the strongly acetylated elements; nothing tried so far targets spread.
      Also makes targets commensurable across cell types, which joint training requires.
- [ ] **Pretrained sequence trunk** (Borzoi or Enformer embeddings, frozen first, then
      fine-tuned) replacing the 8-layer dilated stack. Those representations were fit across
      hundreds of cell types, so cross-cell-type grammar is imported rather than learned from
      two experiments.
- [ ] **DNase-input panel across six or more cell types** instead of ATAC across three. The
      confound that ruled DNase out was about mixing assays inside one comparison; a
      DNase-only panel is internally consistent and doubles the available cell types.

### Lower confidence, cheap once the above exist

- [ ] Gradient-reversal on a cell-type classifier reading the sequence embedding, penalising
      anything cell-type-identifiable in the sequence representation. Meaningless before joint
      training exists.
- [ ] Multiplicative two-tower factorisation: sequence predicts regulatory potential,
      accessibility predicts availability, output is their product. The additive trunk lets
      accessibility dominate. The gate tested sequence modulating accessibility and failed;
      the reverse conditioning is untested.
- [ ] Multi-task across marks available in six cell types (H3K4me1, H3K27me3, CTCF) rather
      than p300, which exists in two. The p300 multi-head died on p300's learnability at
      specific elements, not on the multi-task idea.

### Do not reopen without new evidence

- Further single-cell-type architecture variants. Wide receptive field, fragment channels,
  the sequence gate, the asymmetric loss and GC-matched negatives have all either failed
  outright or failed to travel.

### Structural blocker

- [ ] **A downstream benchmark in a second cell type.** CRISPR data exists only for K562, so
      the terminal metric can only test the transfer direction that already fails. MPRA or
      eQTL data in a second cell type would be worth more right now than another model.
