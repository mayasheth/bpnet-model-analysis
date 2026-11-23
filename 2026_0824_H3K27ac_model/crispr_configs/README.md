# CRISPR benchmark configs for the EP300 runs

Copies of the configs used to benchmark the EP300 ABC variants against CRISPR data with [EngreitzLab/CRISPR_comparison](https://github.com/EngreitzLab/CRISPR_comparison). The originals sit untracked in `/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison` (branch `ms_dev`), which is where `scripts/4.6.setup_crispr_comparison.py` writes them and where `4.7` runs the pipeline. Copied 2026-11-23 so they are kept in git.

- `config/`: one `config_p300_<run>.yml` per benchmark run (19).
- `pred_config/`: the matching `pred_config_p300_<run>.tsv` (predictor columns, aggregation, colors).
- `cell_type_mapping/`: the per-sample cell-type mapping files the configs point to.

Paths inside the configs are absolute and still point at the CRISPR_comparison_v3 clone. Results are in that clone under `workflow/results/2026_09*_p300_*` and `2026_10*_p300_*`; `4.8`, `4.18`, `4.35`, `4.41` and `5.1` read them from there.

To rerun: the configs read ABC predictions from `ABC_working/ABC-Enhancer-Gene-Prediction/results/<run>/<sample>/Predictions/`, which were deleted during offboarding (2026-11). Rerun ABC with the configs in `../abc_configs/` first.
