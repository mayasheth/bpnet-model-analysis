# INVALID — do not use these numbers

Epoch logs from the first residual-training run (jobs 40876817-22, 2026-08-26).

The offset model was `models/atac_hw500_clw1000`, trained on the **250 bp
fragment-extended** target, while these models fit the **5-prime** target. Fragment counts
are ~250x larger, so offsets arrived near 8.15 in log space against a target near 2.9.
Training converged cleanly (valCountPearson 0.813-0.854) because correlation is
shift-invariant, but the learned component is the residual PLUS a scale correction, and
there was no matched ATAC-only baseline on the 5-prime target to compare against.

Kept only as the record of the mistake. The model checkpoints were deleted so they cannot
be evaluated by accident. See `.living/learnings.md` (2026-08-26) and
`results/TARGET_PROVENANCE.md`.
