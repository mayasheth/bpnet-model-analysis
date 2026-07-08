# bpnet-refactor local patch

This directory vendors the local modifications made on top of
[kundajelab/bpnet-refactor](https://github.com/kundajelab/bpnet-refactor)
(BPNet reference implementation used for model training, prediction, and SHAP
in this project). It does not contain a full copy of the reference repo —
only what diverges from the pinned upstream commit.

## Contents

- `UPSTREAM_COMMIT.txt` — the upstream commit the local checkout is based on.
- `modified_files.diff` — unified diff of edits to existing upstream files:
  - `bpnet/cli/argparsers.py`
  - `bpnet/cli/gc/get_gc_matched_negatives.py`
  - `bpnet/cli/shap_scores.py`
  - `bpnet/generators/generators.py`
  - `bpnet/utils/mean_predictions.py`
- `new_files/` — new scripts that don't exist upstream, added under this
  project's local checkout. All four are actively used by this repo's
  pipeline scripts (`0.0.log.sh`, `3.1.submit_mean_shap_one_fold.sh`, etc.):
  - `bpnet/cli/shap_split.py`
  - `bpnet/utils/mean_shap_plus_peaks.py`
  - `bpnet/utils/merge_shap_across_chrom.py`
  - `bpnet/utils/run_importance_hdf5_to_bigwig.py`

Excluded: a `mine/` scratch directory and a stray test file in the local
checkout that aren't referenced anywhere in this project's scripts.

## Reconstructing the environment

```bash
git clone https://github.com/kundajelab/bpnet-refactor.git
cd bpnet-refactor
git checkout $(cat ../external/bpnet-refactor-patch/UPSTREAM_COMMIT.txt)
git apply ../external/bpnet-refactor-patch/modified_files.diff
cp -r ../external/bpnet-refactor-patch/new_files/* .
```

The live working copy used to run this project's pipeline is at
`/oak/stanford/groups/engreitz/Users/sheth/bpnet-refactor`.
