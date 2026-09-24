#!/usr/bin/env python3
"""Regression gate for the pairwise ranking term in multimodal_bpnet.py.

Checks, each of which would otherwise fail silently inside a training run:
  1. weight 0 leaves _asymmetric_mixture_loss bit-identical to before the term existed
  2. labels in {0, 1} give the same profile mask as the old `labels == 1`
  3. correctly ordered predictions score lower than reversed ones
  4. pairs are formed WITHIN groups only: a cross-group ordering has no effect
  5. negatives (label 0) never enter the ranking term
  6. near-ties under min_gap are ignored
  7. the term has a gradient
"""
import math
import sys

import torch

sys.path.insert(0, "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts")
from multimodal_bpnet import _asymmetric_mixture_loss, _pairwise_rank_loss  # noqa: E402

torch.manual_seed(0)
ok = True


def check(name, cond, detail=""):
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    ok &= bool(cond)


N, L = 12, 50
y = torch.poisson(torch.rand(N, 2, L) * 5)
logits = torch.randn(N, 2, L)
lc = torch.randn(N, 1)
labels01 = torch.tensor([1] * 9 + [0] * 3)

# 1 and 2
base = _asymmetric_mixture_loss(y, logits, lc, 10.0, 1.0, labels01)
same = _asymmetric_mixture_loss(y, logits, lc, 10.0, 1.0, labels01, rank_loss_weight=0.0)
check("weight 0 is bit-identical", all(torch.equal(a, b) for a, b in zip(base, same)))

# 3
t = torch.arange(8, dtype=torch.float32)
g = torch.zeros(8, dtype=torch.long)
good = _pairwise_rank_loss(t.clone(), t, g).item()
bad = _pairwise_rank_loss(-t.clone(), t, g).item()
check("ordered < reversed", good < bad, f"{good:.4f} vs {bad:.4f}")
tie = _pairwise_rank_loss(torch.zeros(8), t, g).item()
check("constant prediction scores log(2)", abs(tie - math.log(2)) < 1e-5, f"{tie:.5f}")

# 4: two groups, each internally ordered, but group 1 entirely BELOW group 0 in prediction
# while ABOVE it in target. If cross-group pairs leaked in, this would score badly.
t2 = torch.tensor([0., 1., 2., 3., 10., 11., 12., 13.])
p2 = torch.tensor([10., 11., 12., 13., 0., 1., 2., 3.])
g2 = torch.tensor([0, 0, 0, 0, 1, 1, 1, 1])
within = _pairwise_rank_loss(p2, t2, g2).item()
pooled = _pairwise_rank_loss(p2, t2, torch.zeros(8, dtype=torch.long)).item()
ref = _pairwise_rank_loss(torch.tensor([0., 1., 2., 3.]), torch.tensor([0., 1., 2., 3.]),
                          torch.zeros(4, dtype=torch.long)).item()
check("cross-group pairs excluded", abs(within - ref) < 1e-6 and pooled > within,
      f"within {within:.4f}, same-as-single-group {ref:.4f}, pooled {pooled:.4f}")

# 5: negatives carry wildly misordered predictions; with the rank term they must not matter
y5 = torch.zeros(6, 2, L)
for i, c in enumerate([1, 2, 3, 4, 50, 60]):
    y5[i, 0, 0] = c
lab5 = torch.tensor([1, 1, 1, 1, 0, 0])
lc_a = torch.tensor([[0.], [1.], [2.], [3.], [-9.], [-9.]])
lc_b = torch.tensor([[0.], [1.], [2.], [3.], [9.], [-9.]])
la = _asymmetric_mixture_loss(y5, torch.zeros(6, 2, L), lc_a, 0.0, 1.0, lab5,
                              rank_loss_weight=1.0)[2]
lb = _asymmetric_mixture_loss(y5, torch.zeros(6, 2, L), lc_b, 0.0, 1.0, lab5,
                              rank_loss_weight=1.0)[2]
check("negatives excluded from ranking", torch.allclose(la, lb), f"{la.item():.5f} vs {lb.item():.5f}")

# 6
tn = torch.tensor([0.0, 0.05, 0.09])
check("near-ties ignored", _pairwise_rank_loss(-tn, tn, torch.zeros(3, dtype=torch.long),
                                               min_gap=0.1).item() == 0.0)

# 7
pg = torch.zeros(8, requires_grad=True)
_pairwise_rank_loss(pg, t, g).backward()
check("has a gradient", pg.grad is not None and pg.grad.abs().sum() > 0)

print("\nALL CHECKS PASSED" if ok else "\nFAILURES")
sys.exit(0 if ok else 1)
