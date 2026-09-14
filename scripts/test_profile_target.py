#!/usr/bin/env python
"""Regression gate for the separate profile-head target in multimodal_bpnet.py.

The change lets the profile head train against one signal while the counts head trains
against another (DNase shape, H3K27ac counts). multimodal_bpnet.py and
train_multimodal_bpnet.py are shared with the p300 models and with every existing H3K27ac
arm, so the properties that must hold are:

  1. y_profile=None is bit-identical to before, i.e. this cannot move a published number
  2. y_profile=y is identical to y_profile=None -- the feature is a no-op when the two
     targets are the same signal, which is the only way to be sure the plumbing is not
     quietly reshaping something
  3. a DIFFERENT y_profile moves the profile term and leaves the count term untouched.
     This is the actual separation claim and it is the one that would fail silently:
     a mis-plumbed target would still train and still report a plausible loss.
  4. a profile target whose channel count disagrees with n_outputs raises, rather than
     broadcasting into a meaningless comparison
  5. the four dataset layouts put the right tensor in the right slot, and fit()'s
     positional reads recover them. (X, offset, y, label) and (X, y, label, y_profile)
     are both length 4, so arity cannot disambiguate and the flag must.
  6. reverse-complement augmentation flips each target within its own block
"""
import sys
import numpy as np
import torch

sys.path.insert(0, "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts")
from bpnetlite.losses import _mixture_loss
from multimodal_bpnet import _asymmetric_mixture_loss
from train_multimodal_bpnet import MultiModalPeakNegativeSampler

torch.manual_seed(0)
N, OUT = 16, 1000
y = torch.poisson(torch.full((N, 2, OUT), 3.0))
y_alt = torch.poisson(torch.full((N, 2, OUT), 7.0))     # the "DNase" stand-in
logits = torch.randn(N, 2, OUT)
logcounts = torch.log(y.reshape(N, -1).sum(-1).reshape(N, 1) + 1) + torch.randn(N, 1) * 1.5
labels = (torch.rand(N) > 0.3).long()

print("1. y_profile=None must still equal bpnetlite exactly")
for lab in (None, labels):
    a = _mixture_loss(y, logits, logcounts, 10.0, lab)
    b = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, lab)
    for i, nm in enumerate(("profile", "count", "total")):
        d = abs(a[i].item() - b[i].item())
        assert d <= 1e-5 * max(abs(a[i].item()), 1.0), f"{nm} differs by {d:.3e}"
    print(f"   labels={'yes' if lab is not None else 'no ':<3}  MATCH")

print("\n2. y_profile=y must be identical to y_profile=None")
for lab in (None, labels):
    a = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, lab)
    b = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, lab, y_profile=y)
    for i, nm in enumerate(("profile", "count", "total")):
        assert abs(a[i].item() - b[i].item()) <= 1e-6 * max(abs(a[i].item()), 1.0), \
            f"{nm}: {a[i].item():.8f} vs {b[i].item():.8f}"
    print(f"   labels={'yes' if lab is not None else 'no ':<3}  "
          f"profile {b[0].item():.6f}  count {b[1].item():.6f}  MATCH")

print("\n3. a different y_profile must move the profile term and NOT the count term")
base_p, base_c, _ = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels)
alt_p, alt_c, _ = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels,
                                           y_profile=y_alt)
assert alt_c.item() == base_c.item(), (
    f"count loss moved from {base_c.item():.8f} to {alt_c.item():.8f}; the counts head is "
    f"reading the profile target, which is the whole failure this test exists to catch")
assert abs(alt_p.item() - base_p.item()) > 1e-3 * abs(base_p.item()), (
    f"profile loss did not move ({base_p.item():.6f} -> {alt_p.item():.6f}); the profile "
    f"head is NOT reading y_profile and the auxiliary task would be a no-op")
print(f"   profile {base_p.item():.4f} -> {alt_p.item():.4f}   MOVED")
print(f"   count   {base_c.item():.8f} -> {alt_c.item():.8f}   UNCHANGED (bit-identical)")

print("\n4. a channel-count mismatch must raise")
try:
    _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels,
                             y_profile=y_alt[:, :1, :])
except ValueError as e:
    assert "n_outputs" in str(e), f"raised, but not the intended error: {e}"
    print(f"   raised ValueError as intended")
else:
    raise AssertionError("a 1-channel profile target against a 2-channel head did not "
                         "raise; it would have been compared after a silent reshape")

print("\n5. the four dataset layouts")
NP, NN, LIN, LOUT = 8, 4, 2114, 1000
rng = np.random.RandomState(0)
mk = lambda n, c, L: rng.rand(n, c, L).astype(np.float32)
pk = dict(peak_seqs=mk(NP, 4, LIN), peak_accs=mk(NP, 1, LIN), peak_signals=mk(NP, 2, LOUT),
          neg_seqs=mk(NN, 4, LIN), neg_accs=mk(NN, 1, LIN), neg_signals=mk(NN, 2, LOUT))
prof = dict(peak_profiles=mk(NP, 2, LOUT), neg_profiles=mk(NN, 2, LOUT))
offs = dict(peak_offsets=rng.rand(NP).astype(np.float32),
            neg_offsets=rng.rand(NN).astype(np.float32))
common = dict(in_window=LIN, out_window=LOUT, max_jitter=0, mode="multimodal",
              random_state=1, negative_ratio=0.1)

for name, extra, want_len, has_prof in (
        ("plain", {}, 3, False),
        ("residual", offs, 4, False),
        ("profile target", prof, 4, True),
        ("both", {**offs, **prof}, 5, True)):
    ds = MultiModalPeakNegativeSampler(**pk, **common, **extra)
    assert ds.use_profile_target == has_prof
    item = ds[0]
    assert len(item) == want_len, f"{name}: arity {len(item)}, expected {want_len}"
    # Exactly fit()'s positional reads.
    if has_prof:
        X, yy, lab, yp = item[0], item[-3], item[-2], item[-1]
        off = item[1] if len(item) > 4 else None
    else:
        X, yy, lab = item[0], item[-2], item[-1]
        off = item[1] if len(item) > 3 else None
        yp = None
    assert X.shape == (5, LIN), f"{name}: X {tuple(X.shape)}"
    assert yy.shape == (2, LOUT), f"{name}: y {tuple(yy.shape)}"
    assert lab in (0, 1), f"{name}: label {lab!r} is not a label"
    assert (off is not None) == bool(extra.get("peak_offsets") is not None), \
        f"{name}: offset slot wrong"
    if off is not None:
        assert off.ndim == 0, f"{name}: offset should be a scalar, got {tuple(off.shape)}"
    assert (yp is not None) == has_prof, f"{name}: profile slot wrong"
    if yp is not None:
        assert yp.shape == (2, LOUT), f"{name}: y_profile {tuple(yp.shape)}"
        assert not torch.equal(yp, yy), \
            f"{name}: y_profile equals y, so the slots are crossed"
    print(f"   {name:<15} arity {len(item)}  X {tuple(X.shape)}  y {tuple(yy.shape)}"
          f"  offset {'yes' if off is not None else 'no '}"
          f"  y_profile {'yes' if yp is not None else 'no'}")

print("\n6. reverse complement must flip each target inside its own block")
# One deterministic region, RC forced by driving the dataset's own RNG.
ds = MultiModalPeakNegativeSampler(**pk, **common, **prof, reverse_complement=True)
seen_flipped = False
for k in range(200):
    it = ds[k]
    yy, yp = it[-3], it[-1]
    i_pk = None
    for i in range(NP):
        if torch.equal(yy, torch.from_numpy(pk["peak_signals"][i])):
            i_pk = i
            break
        if torch.equal(yy, torch.flip(torch.from_numpy(pk["peak_signals"][i]), [0, 1])):
            i_pk, seen_flipped = i, True
            # The paired profile target must be flipped THE SAME WAY and taken from the
            # same region. A single flip over a concatenated tensor would instead hand
            # back the two blocks in swapped order.
            assert torch.equal(
                yp, torch.flip(torch.from_numpy(prof["peak_profiles"][i]), [0, 1])), \
                "y was flipped but y_profile was not flipped to match"
            break
    if i_pk is not None and seen_flipped:
        break
assert seen_flipped, ("never drew a reverse-complemented peak in 200 draws, so this check "
                      "proved nothing. Fix the sampling, do not drop the assertion.")
print("   a flipped y came back with a flipped y_profile from the same region")

print("\n7. profile_loss_weight must scale only the total, and 1.0 must be a no-op")
p0, c0, t0 = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels,
                                      y_profile=y_alt)
p1, c1, t1 = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels,
                                      y_profile=y_alt, profile_loss_weight=1.0)
assert (p0.item(), c0.item(), t0.item()) == (p1.item(), c1.item(), t1.item()), \
    "profile_loss_weight=1.0 is not a no-op"
for w in (0.0561, 0.5, 2.0):
    pw, cw, tw = _asymmetric_mixture_loss(y, logits, logcounts, 10.0, 1.0, labels,
                                          y_profile=y_alt, profile_loss_weight=w)
    # The reported profile loss stays UNWEIGHTED so it remains comparable to every
    # existing run's logs; only its contribution to the total is scaled.
    assert pw.item() == p0.item(), f"w={w}: reported profile loss was scaled too"
    assert cw.item() == c0.item(), f"w={w}: count loss moved"
    expect = w * p0.item() + 10.0 * c0.item()
    assert abs(tw.item() - expect) <= 1e-4 * max(abs(expect), 1.0), \
        f"w={w}: total {tw.item():.6f} vs expected {expect:.6f}"
    print(f"   w={w:<7} total {t0.item():>12.2f} -> {tw.item():>12.2f}  "
          f"(expected {expect:>12.2f})  MATCH")

print("\nALL CHECKS PASSED")
