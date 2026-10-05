"""Reweight ABC's activity term in scE2G output without the (2-3 digit) numerators.

scE2G's ABC (powerlaw) per gene g: the self-promoter pair has ABC = 1 exactly; every other pair has
ABC_e = num_e / D_g with D_g = sum of numerators over ALL the gene's candidate pairs, the
self-promoter included (checked 2026-10-04 on CEP104: D constant to 4 digits, self-promoter 1.0).
So the non-self ABC values are shares of D_g, and the self-promoter's share is 1 - sum(non-self ABC)
(0 when the gene's promoter is not a candidate). Multiplying each element's activity by w_e gives

    ABC'_e = ABC_e w_e T / (sum_nonself ABC w + (1 - sum_nonself ABC) w_self)    (non-self pairs)
    ABC'_self = 1

ARC-E2G is ABC x a per-pair Kendall factor (exp(21.5 x Kendall) in the Xu run), so
ARC'_e = ARC_e x ABC'_e / ABC_e = ARC_e w_e T / denominator, and the self-promoter's ARC is unchanged.
T = the gene's total share (1 up to rounding), so that w = 1 reproduces ABC exactly.
"""
import numpy as np
import pandas as pd


def reweight(gene, is_self, abc, arc, w):
    """gene, is_self, abc, arc, w: aligned 1-D arrays over a gene's full candidate set."""
    gene = np.asarray(gene)
    is_self = np.asarray(is_self, bool)
    abc, arc, w = (np.asarray(x, float) for x in (abc, arc, w))
    ns = ~is_self
    s_ns = pd.Series(np.where(ns, abc, 0.0)).groupby(gene).transform("sum").values
    sw_ns = pd.Series(np.where(ns, abc * w, 0.0)).groupby(gene).transform("sum").values
    w_self = pd.Series(np.where(is_self, w, np.nan)).groupby(gene).transform("max").values
    share_self = np.clip(1.0 - s_ns, 0.0, None)
    # dividing by the total share makes w = 1 exact even where rounded ABC shares do not sum to 1
    den = (sw_ns + share_self * np.nan_to_num(w_self, nan=0.0)) / (s_ns + share_self)
    f = np.divide(w, den, out=np.ones_like(den), where=den > 0)
    abc_new = np.where(ns, abc * f, abc)
    arc_new = np.where(ns, arc * f, arc)
    return abc_new, arc_new


def check(gene, is_self, abc, arc):
    """Assert the model of scE2G's ABC holds and that w = 1 reproduces ABC and ARC."""
    is_self = np.asarray(is_self, bool)
    assert np.allclose(np.asarray(abc)[is_self], 1.0), "self-promoter ABC is not 1"
    s_ns = pd.Series(np.where(~is_self, abc, 0.0)).groupby(np.asarray(gene)).sum()
    assert (s_ns <= 1 + 1e-3).all(), f"non-self ABC sums exceed 1 (max {s_ns.max():.4f})"
    a1, r1 = reweight(gene, is_self, abc, arc, np.ones(len(abc)))
    assert np.allclose(a1, abc, rtol=1e-9) and np.allclose(r1, arc, rtol=1e-9), "w = 1 does not reproduce ABC / ARC"
    return s_ns
