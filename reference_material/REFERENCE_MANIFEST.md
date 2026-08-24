# Reference Manifest

Reference material for this project lives in **`reference/`** (data inputs, also catalogued
in `data/DATA_MANIFEST.md`) and **`external/`** (vendored third-party code), not in
`reference_material/`, which mycelium init created empty. Nothing was moved.

---

### bpnet-refactor

```yaml
name: bpnet-refactor
type: code
location: $OAK/Users/sheth/bpnet-refactor
local_patch: external/bpnet-refactor-patch/
status: active
tags: [bpnet, training, external-codebase, vendored-patch]
```

The BPNet training/prediction codebase all sequence-only models in this repo use. It lives
**outside this repo**; local modifications are vendored as a patch in
`external/bpnet-refactor-patch/` (added 2026-07-08) containing `modified_files.diff`,
`new_files/`, `UPSTREAM_COMMIT.txt`, and a README. See `README.md` § "External codebases".

Models predict a **profile head** (base-resolution ChIP-seq shape) and a **counts head**
(total log-scale reads) — counts is the primary output used in analyses.

---

### style-guidelines

```yaml
name: style-guidelines
type: documentation
location: scripts/STYLE_GUIDELINES.md
status: active
tags: [figures, plotting, house-style]
```

Plot formatting rules: sentence case for axis titles, "p300" not "P300", no gridlines,
black axes. Palettes — diverging `managua`, sequential `PuBu`, p300 status `#792374`
(p300+) / `#49bcbc` (p300-). Complements the lab-wide
`.living/conventions/engreitz-lab/figure-conventions.md`.

---

### workflow-logs

```yaml
name: workflow-logs
type: documentation
location: [2025_0517_official_EP300_K562_model/scripts/0.0.log.sh, 2026_0606_GM12878_transferability/0.0.log.sh, log.sh]
status: active
tags: [provenance, commands, reproducibility]
```

Append-only shell logs holding the literal command history — model downloads, prediction
and SHAP submissions, MoDISCo parameterizations, FiNeMo calls, `conda activate` /
`module load` lines. These are the ground truth for how any given output was produced.

---

### multimodal-architecture-writeup

```yaml
name: multimodal-architecture-writeup
type: documentation
location: 2026_0529_multimodal_p300_model/multimodal_bpnet_architecture.html
status: active
tags: [architecture, multimodal, design-decisions]
```

Interactive write-up of every parameter and design decision in the multimodal
architecture — the reference for why the model is shaped the way it is.

---

### motif-compendium-and-inputs

```yaml
name: motif-compendium-and-inputs
type: data
location: reference/
status: active
tags: [see-data-manifest]
```

Shared data inputs (peaks, candidate-element FASTAs, motif databases, CV folds, GC
negatives). Fully catalogued in `data/DATA_MANIFEST.md` rather than duplicated here.
