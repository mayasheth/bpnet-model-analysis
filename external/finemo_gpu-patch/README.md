# finemo_gpu local patch

Local edits on top of [austintwang/finemo_gpu](https://github.com/austintwang/finemo_gpu), the FiNeMo version used for motif hit calling in this project. Only what differs from the pinned upstream commit is kept here, not a full copy of the repo.

## Contents

- `UPSTREAM_COMMIT.txt`: the upstream commit the local checkout was based on.
- `modified_files.diff`: edits to existing upstream files.
  - `src/finemo/data_io.py`: adds `load_curated_modisco_motifs(prefix)`, which reads a curated motif set saved as `<prefix>.motifs.parquet` plus `<prefix>.motifs.npz` (written by `scripts/extract_modisco_motifs.py`).
  - `src/finemo/main.py`: `call-hits` and `report` accept either a modisco `.h5`/`.hdf5` or a curated-motif prefix for `--modisco-h5`; with a prefix, recall against seqlets is skipped.
  - `src/finemo/evaluation.py`: motif logos are saved as SVG as well as PNG.
- `new_files/environment_v2.yml`: conda spec for the `finemo_gpu` environment that `scripts/5.1.submit_finemo.sh` activates.

## Rebuilding

```bash
git clone https://github.com/austintwang/finemo_gpu.git
cd finemo_gpu
git checkout $(cat ../EP300_BPNet/external/finemo_gpu-patch/UPSTREAM_COMMIT.txt)
git apply ../EP300_BPNet/external/finemo_gpu-patch/modified_files.diff
conda env create -n finemo_gpu -f ../EP300_BPNet/external/finemo_gpu-patch/new_files/environment_v2.yml
conda activate finemo_gpu && pip install -e .
```

The original working copy was at `/oak/stanford/groups/engreitz/Users/sheth/finemo_gpu` and was removed during offboarding (2026-11).
