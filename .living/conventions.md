# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

## Layout: analyses live at the repo root, not under `analysis/`

Mycelium was initialized here on **2026-08-24**, years into the project. The canonical
mycelium directories (`analysis/`, `data/`, `algorithms/`, `reference_material/`) were
created but are **intentionally empty except for their manifests**, no files were moved,
because doing so would break every absolute path in the workflow logs
(`*/0.0.log.sh`), the SLURM submit scripts, and `CLAUDE.md`.

The real layout is:

| Mycelium expects | This repo uses |
|---|---|
| `analysis/<name>/` | dated dirs at root: `2025_0517_official_EP300_K562_model/`, `2026_0606_GM12878_transferability/`, … and model dirs `K562_*_BPNet/`, `K562_*_ChromBPNet/` |
| `data/raw`, `data/processed` | `reference/` for shared inputs; each analysis's own `data/` for project-specific signal |
| `algorithms/` | `scripts/` at root (shared) and each analysis's `scripts/` |
| `reference_material/` | `external/` for vendored third-party code; `reference/` for data |

**The manifests are the index into that real layout**, every entry carries a `path:` field
pointing at the actual location. Read the manifest, not the directory tree.

Naming convention for analyses: `YYYY_MMDD_short_description/`, where the date is when the
analysis started. Model directories predating that convention use `<CELLTYPE>_<TARGET>_<ARCH>/`.

## Reports

`/engreitzlab-report` (-> `/analysis-report`, Quarto -> self-contained HTML) supersedes
mycelium's `report-generator` pack for anything written up here. The `report:` field in
`analysis/ANALYSIS_MANIFEST.md` is `null` for every entry because no analysis has been
routed through that path yet, manuscript figures are tracked in `README.md` instead.

## Figures

`scripts/STYLE_GUIDELINES.md` is the binding style reference for this project and takes
precedence over `.living/conventions/engreitz-lab/figure-conventions.md` where they differ:
sentence case axis titles, "p300" never "P300", no gridlines, black axes, diverging
`managua` / sequential `PuBu`, p300 status `#792374` (p300+) and `#49bcbc` (p300-).

## Environments

This repo predates the lab pixi standard and uses a mix: `pixi.toml` at root (Python 3.9)
plus three conda envs, `bpnet_37` (training/prediction/SHAP, needs
`module load cuda/11.1.1 cudnn/8.1.1.33`), `tfmodisco` (FIMO, MoDISCo, inference), and
`analysis` (FiNeMo formatting, downstream plots). The multimodal project has its own pixi
env named `multimodal`. Do not assume a single environment works across stages.

## SLURM submit scripts (mandatory)

Every new submit script in this repo must:

1. **`export PYTHONUNBUFFERED=1`.** Python block-buffers stdout to a file, so a running
   job shows an empty log and looks hung. Without this, a 35-minute job is
   undiagnosable while alive.
2. **Guard empty array expansion.** `set -u` plus `"${ARR[@]}"` on an empty array is an
   unbound-variable error in bash < 4.4. Use `${ARR[@]+"${ARR[@]}"}` for any optional
   argument array.
3. **Put the varied hyperparameter in the output path.** A sweep whose runs share an
   output directory silently overwrites itself. See
   `2026_0824_H3K27ac_model/scripts/1.1.submit_training.sh`, where the counting window
   and count-loss weight both appear in `OUT_DIR`.
4. **Bound memory explicitly when a set size scales with the data.** Anything of the
   form `n_items * k` needs a cap; see `--max-negatives`.

## Never run heavy work on a login node

Anything beyond a few seconds goes through `sbatch`. A long python process started over
`ssh` on a login node is killed when the calling shell detaches, it exits 0 and writes
no output, which reads as success. Login-node work is limited to inspecting files and
small tabulations; the one window-profiling script that does run there is explicitly
memory-bounded (264 MB) and documented as such.

## Model evaluation must be stratified by signal level

`reference/K562_DNase_candidate_elements.narrowPeak` holds 150,528 elements and most
carry little signal, so any correlation computed over all of them is dominated by the
dead-vs-active contrast. Report the top signal quintile alongside the overall number -
this repo's equivalent of the all-vs-p300+ split already used for p300. See
`2026_0824_H3K27ac_model/scripts/2.2.evaluate_stratified.py`. Concretely: ATAC-only
predicts H3K27ac at 0.746 over all elements but only 0.543 on the top quintile.

Ceilings derived from replicate agreement need both the Spearman-Brown and sqrt
corrections before they bound model performance, see `.living/learnings.md`.

## Decision-log entries use `###`, not `##`

Mycelium's `decision-log-entry.md` template says `## [YYYY-MM-DD] Title`, but
`generate_index.py` only counts `### ` headings in `decisions.md`. Entries written per
the template are silently absent from `.living/INDEX.md`. Use `### ` for decisions here.
See `.living/learnings.md` (2026-08-25).

## Never take a model target from a display directory

Tracks under `$OAK/Users/sheth/Data/share/IGV/` are built for visualization. They may be
raw-count and unnormalized, and so look usable, while carrying processing choices that
are wrong for modelling. `ENCSR000AKP_coverage.bw` extends 36 bp reads to a fixed 250 bp
fragment, which flattens nucleosome structure and breaks the multinomial assumption
behind bpnetlite's profile loss.

Derive training targets from the source BAMs inside the analysis project, so the
processing is explicit, versioned, and reviewable. For ChIP-seq in this repo that means
`bedtools genomecov -5 -dz` (single-base 5' ends, stranded), matching
`scripts/0.3.make_training_bw.sh`. If a display track is used as a stopgap, record the
exact command that produced it and treat replacing it as an open item.

## H3K27ac target is 5-prime ends, full stop

`2026_0824_H3K27ac_model/data/h3k27ac_5p_{plus,minus}.bw`. Nothing new is trained or
evaluated on the 250 bp fragment-extended track, it breaks MNLL's read-count assumption
and bleeds signal between neighbouring elements. See
`2026_0824_H3K27ac_model/results/TARGET_PROVENANCE.md` for which existing results used
which target, and which conclusions survive.

Enforcement: `scripts/train_multimodal_bpnet.py` writes `training_target.json` into every
model directory and refuses a `--count-offset-model` whose signal tracks or output window
differ from the current run. This exists because a fragment-trained offset model was
silently used against a 5-prime target on 2026-08-26 and trained to a plausible-looking
but meaningless result.

## Sample definition for multi-dataset comparisons

Two rules when assembling panels from the lab data inventory
(`2023_0701 Maya-s sequences.xlsx`, "ENCODE data" sheet):

1. **No Mint-ChIP.** Different assay chemistry from standard ChIP-seq and not comparable.
   Excludes WTC11 H3K27ac/H3K27me3/H3K4me1 (ENCSR146DPQ, ENCSR418YEV, ENCSR473GYW).
2. **A sample is one (experiment, processing) pair.** Replicates may be pooled, or
   compared as an inter-replicate ceiling, only within a single experiment AND a single
   processing. Rows in the sheet sometimes pool two experiments (all four TeloHAEC
   H3K27ac rows list GSE210489 and GSE210491) or mix run types across replicates
   (HCT116 rows are `se, pe`); both must be split before use. Pooling across either
   folds batch or assay differences into what gets reported as biological
   reproducibility, which understates the ceiling.

Recorded because the ceiling is the denominator for every performance number in the
H3K27ac work, so contaminating it silently rescales all of them.

## Refreshing the .living index

`.living/INDEX.md` is GENERATED and drifts silently. It sat at 24 learnings / 4 decisions
while the sources held 28 / 6, which defeats its purpose as the first thing a fresh session
reads. Regenerate it after adding entries:

```bash
module load python/3.12.1
MYC=/oak/stanford/groups/engreitz/Users/sheth/EngreitzLabAgents/mycelium-upstream
python3 $MYC/skills/core/scripts/generate_index.py \
    --living-dir /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.living \
    --summary-heuristic
```

Use `--summary-heuristic` (tag-based clustering, no LLM), that is what the existing index
was built with, so the format stays consistent. Requires Python 3.11+, hence the
`module load`; over a non-login `ssh host cmd` this must be wrapped in `bash -lc`.

Do NOT hand-write INDEX.md or reimplement the generator. A stand-in was written once during
a session where a `find` for the upstream copy timed out before reaching the
EngreitzLabAgents checkout; it produced a subtly different tag-cluster format. The canonical
script is on Oak at the path above.


## Never end a sentence with a p-value in the report

`render_report.py`'s numeric-claim regex captures a trailing period into the token, so
`p = 0.53.` is parsed as the number `0.53.`, which parses as nothing and is reported as an
untraceable claim even when `0.53` is registered in `numbers.json`. This has cost a
render-fix cycle three times.

Put a word or a comma after every p-value: `p = 0.53 in both strata`, `p = 0.090, consistent
with no effect`, `at p = 0.086, and`. The same applies to any figure or CI value that would
otherwise sit at the end of a sentence.


## Completion markers and the queue can disagree

`training_complete.json` is written when training finishes, but the process does not always
exit. One wide-receptive-field job wrote its marker and final checkpoint at 06:01 and was
still holding a GPU at 10:39, 4.5 hours, with the log ending at "Model saved" and no
"Done:" line. On `scancel` it ignored SIGTERM and slurmd gave up after 927 s
("JOB NOT ENDING WITH SIGNALS"), so it was a real teardown hang, most likely dataloader
workers not joining.

Two consequences:

- A family showing 5/5 markers does not mean the queue is clear. Check `squeue` separately
  before concluding a stage is finished, and before assuming a GPU shortage is external.
- A job whose log ends at the last thing the Python script printed, with the wrapper's own
  trailing `echo` missing, is hung rather than working. That absent final line is the
  cheapest signal, so keep an `echo` after the python call in every submit script.


## Match the partition to the job, not to habit

**Use the lab's own `engreitz` partition for anything that does not need a GPU.** Nine nodes
at 24+ cores and 192 GB+, 7-day limit, no GPUs. It does not compete with the general GPU
queues, so it sidesteps a fairshare depleted by training. Put it first and let `normal` and
`owners` follow as fallbacks: `#SBATCH -p engreitz,normal,owners`.

Three separate costs today came from putting the wrong job on the wrong partition.

- **Long CPU jobs belong on `normal`, not `owners`.** The GM12878 fragment-channel build
  (~3.5 h, a full pass over three 9 GB BAMs) was preempted after 40 minutes. It requeues
  automatically, but its expensive intermediate lives in `$SCRATCH` under a cleanup trap, so
  a preemption loses everything. `normal` is not preemptible; use it whenever a single
  uncheckpointed step runs longer than about an hour. `owners` remains right for training,
  which is short per fold and cheap to lose.
- **Pure inference does not need a GPU.** The genome-wide ABC predictions sat behind a
  depleted GPU fairshare while hundreds of CPU nodes were idle. `4.1` selects its device via
  `torch.cuda.is_available()`, so the CPU twin (`4.4`) needed no code change and 153,545
  regions finish comfortably within a few hours on 16 threads. Reserve GPU requests for
  training.
- **Submitting a large batch costs priority for everything after it.** Twenty training jobs
  submitted inside an hour drained fairshare, so the work that had actually been prioritised
  queued behind exploratory follow-ups. Submit the thing that unblocks downstream work first,
  and `scontrol hold` the speculative batch rather than racing it.


## Snakemake on Sherlock: wrapper env + profile, matched to the major version

Full table in `CLAUDE.md`. The short form: drive pipelines with
`$OAK/Users/sheth/.conda/envs/run_snakemake` (7.32.4) or `run_snakemake9` (9.6.0), always
pass `--profile` so each rule instance is submitted as its own job, and submit the driver
itself as a small long-lived job on a non-preemptible partition.

The trap worth remembering: `~/.config/snakemake/slurm` and `slurm_long` are **v7-style**,
driving submission through a `cluster:` command string. Snakemake 9 uses the executor plugin
and ignores `cluster:` entirely. v9 profiles live **inside their repo**, not in
`~/.config/snakemake`; copy the working one at
`$OAK/Users/sheth/seq-processing-snakemake/.snakemake_profile/slurm/config.yaml`.

Also: `--use-conda` needs the wrapper env's `bin` on `PATH`, because it shells out to mamba
from `/usr/bin/bash`. Calling snakemake by absolute path is not enough.

## Ask before running a broad search

Before a repo-wide grep, a `find` over Oak, or any hunt for where something lives, ask.
Maya usually knows the path or can narrow it to one directory. A `find` for an element set
timed out and returned nothing; a `grep -rl` across four repos missed a snakemake9 profile
that one question located immediately, in a hidden directory the patterns never covered.
An empty or timed-out search is not evidence of absence. Targeted `ls` of a named directory
does not need asking.

## The TODO list lives in `todo/TODOLIST.md`

`todo/TODOLIST.md` is the live list of open work, and is the file to update when work starts,
finishes or is deferred. It is the mycelium `todo/` directory.

It is the only one. Two files used to look like it: `TODO.md` at the repo root and
`todo/h3k27ac-model.md`. Both were folded into TODOLIST and deleted on 2026-09-12, and their
publication-figure status and p300-model CV values moved to `reference/PUBLICATION_FIGURES.md`.
Do not add a second list. The failure mode those two demonstrated is worth naming: a stale
list is worse than no list, because it describes finished work as unstarted and a reader
plans from it.

The list holds OPEN work. Anything finished belongs in the git history, in `decisions.md` and
in `findings/`. Two exceptions earn a `- [x]` entry that stays visible: a result worth stating
so nobody reopens it blind, and a deferral, which must carry the reason and the condition that
would reopen it. A deferral recorded as a closure is how a cheap idea gets re-litigated three
sessions later.

Entries quote the number and its interval rather than a verdict, and cross-reference finding
IDs, decision dates and script numbers instead of restating content. The TODO is an index of
what is open, not a second copy of the analysis.

This convention is duplicated in `CLAUDE.md` deliberately: CLAUDE.md is loaded into every
session's context automatically and this file is not, so the copy that gets read is the one
that matters.
