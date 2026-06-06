# TODO

## Model improvements

- [ ] Consider switching chromosome fold splits to reduce train/test leakage: evaluate hashFrag (https://github.com/de-Boer-Lab/hashFrag) for sequence-similarity-aware splits
- [ ] Train a p300 BPNet-like model that takes DNA sequence AND base-pair-level chromatin accessibility as input (multimodal architecture)
  - Architecture implemented: `scripts/multimodal_bpnet.py` (MultiModalBPNet), `scripts/train_multimodal_bpnet.py`, `scripts/shap_multimodal_bpnet.py`
  - Model dir: `2025_0529_multimodal_p300_model/`; pixi env: `multimodal` (see `pixi.toml`)
  - Blocked on: ATAC-seq and DNase-seq BAM/tagAlign file paths (needed to generate BigWig inputs)
  - Future improvement: bias-correct accessibility input signal

## Evaluation

- [ ] Compute inter-replicate p300 ChIP-seq Pearson/Spearman correlation as an empirical ceiling to compare against model predictive performance
