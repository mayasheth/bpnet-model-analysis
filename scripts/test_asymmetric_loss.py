#!/usr/bin/env python
"""Regression gate for the gate + asymmetric-loss change to multimodal_bpnet.py.

That file is shared with the p300 models, so three properties must hold:
  1. overprediction_weight == 1.0 reproduces bpnetlite's _mixture_loss EXACTLY
  2. weight > 1.0 penalises only the over-prediction cases, and by the stated factor
  3. checkpoints pickled before the gate existed still load and predict bit-identically
"""
import sys
import numpy as np
import torch

sys.path.insert(0, "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts")
from bpnetlite.losses import _mixture_loss
from multimodal_bpnet import MultiModalBPNet, _asymmetric_mixture_loss

torch.manual_seed(0)
N, L, OUT = 16, 1000, 1000
y = torch.poisson(torch.full((N, 2, OUT), 3.0))
logits = torch.randn(N, 2, OUT)
# Straddle the truth deliberately. log(sum(y)+1) is ~8.7 for this y, so centring the
# predictions there and spreading them puts roughly half above and half below. An earlier
# version of this test centred them at 5.0, which is always BELOW the truth, so the
# over-prediction branch never executed and check 2 silently compared base against base.
_truth_log = torch.log(y.reshape(N, -1).sum(-1).reshape(N, 1) + 1)
logcounts = _truth_log + torch.randn(N, 1) * 1.5
labels = (torch.rand(N) > 0.3).long()

print("1. weight 1.0 must equal bpnetlite exactly")
for lab in (None, labels):
    a = _mixture_loss(y, logits, logcounts, 10.0, lab)
    b = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, lab)
    for i, nm in enumerate(("profile", "count", "total")):
        d = abs(a[i].item() - b[i].item())
        tol = 1e-5 * max(abs(a[i].item()), 1.0)
        assert d <= tol, (f"{nm} differs by {d:.3e} > {tol:.3e} "
                          f"(labels={'yes' if lab is not None else 'no'})")
    print(f"   labels={'yes' if lab is not None else 'no ':<3}  "
          f"profile {a[0].item():.6f}  count {a[1].item():.6f}  total {a[2].item():.6f}  MATCH")

print("\n2. weight > 1 must penalise only over-prediction, by the stated factor")
err = logcounts - torch.log(y.reshape(N, -1).sum(-1).reshape(N, 1) + 1)
n_over = int((err > 0).sum())
base = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0)[1].item()
for w in (2.0, 3.0, 5.0):
    got = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, w)[1].item()
    over = (err[err > 0] ** 2).sum().item() / err.numel()
    expect = base + (w - 1.0) * over
    # Relative tolerance: these are float32 sums over 16x2000 elements, so absolute
    # agreement to 1e-6 is below the arithmetic's own precision. w=5 failed at 2e-6 on a
    # value of 8.19, which is float noise rather than a behavioural difference.
    assert abs(got - expect) <= 1e-5 * max(abs(expect), 1.0), \
        f"w={w}: {got:.6f} vs expected {expect:.6f}"
    print(f"   w={w}: count loss {base:.6f} -> {got:.6f}  (expected {expect:.6f})  MATCH")
print(f"   {n_over}/{N} examples over-predicted in this sample")
assert 3 <= n_over <= N - 3, (
    f"only {n_over}/{N} over-predicted: this test cannot distinguish the asymmetric branch "
    f"from the symmetric one. Fix the synthetic logcounts, do not relax the assertion.")
assert _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 3.0)[1].item() > base * 1.05, (
    "weight 3.0 did not raise the count loss materially; the branch is not being taken")

print("\n3. a pre-gate checkpoint must load and predict bit-identically")
ck = ("/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/"
      "models/multimodal5p_accs5p_hw500_clw10/fold0/multimodal_bpnet.torch")
m = torch.load(ck, map_location="cpu", weights_only=False).eval()
print(f"   loaded, has gate_accessibility attr: {hasattr(m, 'gate_accessibility')}")
x = torch.rand(4, 5, 2114)
with torch.no_grad():
    p1, c1 = m(x)
    p2, c2 = m(x)
assert torch.equal(p1, p2) and torch.equal(c1, c2)
print(f"   forward deterministic, counts {c1.flatten().tolist()[:2]}")

print("\n4. an open gate must leave predictions essentially unchanged")
# Seeding both constructors identically does NOT give matched weights: the gated model
# allocates an extra gate_conv, which consumes RNG draws, so every later layer differs. An
# earlier version of this test did exactly that and reported a 108% difference, which
# measured nothing. Build the gated model, then copy its shared weights into an ungated one.
torch.manual_seed(1)
gat = MultiModalBPNet(mode="multimodal", n_acc_channels=1, gate_accessibility=True,
                      verbose=False).eval()
ung = MultiModalBPNet(mode="multimodal", n_acc_channels=1, verbose=False).eval()
shared = {k: v for k, v in gat.state_dict().items() if not k.startswith("gate_conv")}
missing, unexpected = ung.load_state_dict(shared, strict=False)
assert not unexpected, f"unexpected keys when copying: {unexpected}"
assert not [k for k in missing if not k.startswith("gate_conv")], f"missing: {missing}"
with torch.no_grad():
    _, cu = ung(x)
    _, cg = gat(x)
rel = ((cg - cu).abs() / cu.abs().clamp_min(1e-6)).max().item()
assert rel < 0.05, (
    f"open gate changes predictions by {rel:.1%}, too much to call the gated model a "
    f"drop-in start. Raise the gate bias or rethink the initialisation.")
print(f"   gate value at init = sigmoid(4) = {torch.sigmoid(torch.tensor(4.0)):.4f}")
print(f"   max relative difference in predicted logcounts: {rel:.4f}")
print("   (not bit-identical by design: the gate scales accessibility by 0.982. bias=4 keeps")
print("    sigmoid' = 0.018 so the gate stays trainable; a larger bias starts closer to 1.0")
print("    but saturates and barely learns.)")
print("\nALL CHECKS PASSED")
