#!/usr/bin/env python3
"""Checks for attach_cell_conditioning (FiLM on a per-cell vector), both directions.

1. Zero-initialised: an attached model fed the extra row-index channel matches the plain model.
2. Live: once the FiLM output layer is non-zero, two cells' rows give different predictions,
   and the same row gives the same prediction.
3. Dropout acts only in training mode.
4. The trainer rejects a table missing a training cell and a table with non-finite values.

Run: pixi run -e multimodal python scripts/test_cell_conditioning.py
"""
import copy
import os
import subprocess
import sys
import tempfile

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multimodal_bpnet import MultiModalBPNet, attach_cell_conditioning  # noqa: E402

torch.manual_seed(0)
N, L = 6, 2114
base = MultiModalBPNet(n_acc_channels=3, n_outputs=2, verbose=False).eval()
X = torch.randn(N, 7, L)
X[:, :4] = torch.nn.functional.one_hot(torch.randint(0, 4, (N, L)), 4).permute(0, 2, 1).float()

table = np.array([[1.0, -0.5, 0.2], [-1.0, 0.5, 0.0], [0.3, 0.3, -0.3]], dtype=np.float32)
cond = attach_cell_conditioning(copy.deepcopy(base), table, ["A", "B", "K562"], dropout=0.5).eval()
rows = torch.tensor([0, 1, 2, 0, 1, 2], dtype=torch.float32)
Xc = torch.cat([X, rows.reshape(N, 1, 1).expand(N, 1, L)], dim=1)

with torch.no_grad():
    p0, c0 = base(X)
    p1, c1 = cond(Xc)
assert torch.equal(p0, p1) and torch.equal(c0, c1), "zero-initialised FiLM changed the output"
print("1. zero-init identical: ok")

with torch.no_grad():
    torch.nn.init.normal_(cond.film[2].weight, std=0.05)
    p2, c2 = cond(Xc)
assert not torch.allclose(c2[0], c2[1]) or not torch.allclose(c2[0], c2[2]), "rows have no effect"
Xs = Xc.clone()
Xs[:, :-1] = Xc[0:1, :-1]                      # same DNA/ATAC in every row, rows still 0,1,2,...
with torch.no_grad():
    _, c3 = cond(Xs)
assert torch.allclose(c3[0], c3[3]) and torch.allclose(c3[1], c3[4]), "same row, different output"
assert not torch.allclose(c3[0], c3[1]), "different rows, same output"
print("2. live FiLM: rows change predictions, equal rows agree: ok")

with torch.no_grad():
    a, b = cond(Xs)[1], cond(Xs)[1]
assert torch.equal(a, b), "eval mode is not deterministic"
cond.train()
torch.manual_seed(1)
with torch.no_grad():
    outs = torch.stack([cond(Xs)[1] for _ in range(8)])
assert outs.std(0).max() > 0, "dropout never fired in training mode"
print("3. dropout only in training: ok")

# 4. Trainer guards, run up to the point where the table is read; both must exit non-zero
# with the expected message.
here = os.path.dirname(os.path.abspath(__file__))
panel = os.environ.get("TEST_PANEL_JSON")
if panel:
    for bad, msg in ((["cell\tf1", "nonexistent\t0.0"], "has no row for"),
                     (None, "non-finite")):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            if bad is None:
                import json
                cells = list(json.load(open(panel)))
                f.write("cell\tf1\n" + "".join(f"{c}\tnan\n" for c in cells))
            else:
                f.write("\n".join(bad) + "\n")
        r = subprocess.run([sys.executable, os.path.join(here, "train_multimodal_bpnet.py"),
                            "--cell-types-json", panel, "--holdout-cell", "K562", "--panel-mean",
                            "--peaks", "unused", "--signal-plus-bw", "unused",
                            "--mode", "multimodal", "--fold", os.environ["TEST_FOLD_JSON"],
                            "--output-dir", tempfile.mkdtemp(), "--max-negatives", "50",
                            "--genome", os.environ["TEST_GENOME"],
                            "--negatives", os.environ["TEST_NEGATIVES"],
                            "--cell-cond-tsv", f.name, "--max-epochs", "0"],
                           capture_output=True, text=True)
        assert r.returncode != 0 and msg in (r.stdout + r.stderr), \
            f"guard '{msg}' did not fire: rc={r.returncode}\n{(r.stdout + r.stderr)[-2000:]}"
    print("4. trainer guards reject a missing row and non-finite values: ok")
else:
    print("4. skipped (set TEST_PANEL_JSON, TEST_FOLD_JSON, TEST_GENOME, TEST_NEGATIVES)")
print("all cell-conditioning checks passed")
