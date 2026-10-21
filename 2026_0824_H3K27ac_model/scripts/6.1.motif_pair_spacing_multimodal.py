#!/usr/bin/env python
"""Motif-pair syntax (identity, signed spacing, orientation) learned by the Corgi p300 model.

A cooperator-insertion assay in the style of the Bolero preprint, run on the multi-cell
panel-mean p300 model (1.42). An anchor motif sits at the centre of a dinucleotide-shuffled
2,114 bp window and a partner is inserted at signed edge-to-edge spacing d in both
orientations; the readout is the model's count head.

WHY SHUFFLED DNA WITH REAL ACCESSIBILITY. The model reads 4 DNA channels plus 3 accessibility
channels, and accessibility dominates its predictions. Shuffling the DNA while keeping the
element's own accessibility channels for the chosen cell fixes the accessibility contribution
and leaves a pure sequence contrast, so what is measured is what the sequence branch learned,
at a realistic accessibility operating point.

COOPERATIVITY. For each background b,

    coop(b) = f(anchor + partner) - f(anchor) - f(partner) + f(neither)

with f the fold-averaged predicted log counts. Every term is evaluated on the SAME background
with the SAME accessibility channels, so the panel-mean count offset -- a per-window constant
that 4.1 adds to get deployment-scale counts -- cancels identically and is not added here.
f(partner) is taken at the partner's ACTUAL absolute position in the pair, not at the window
centre, because position inside the window changes the prediction.

CELL CONTEXT. The same shuffled DNA is read under each cell's own accessibility channels at the
same genomic windows, so context enters only through accessibility. Channel construction
follows 4.1 exactly: channel 0 is the target cell's log1p ATAC standardised with statistics
computed over the ABC candidate windows (4.26's recipe; --stage prep recomputes them for every
cell so all five contexts are normalised the same way, which is what the K562 deployment arm
did and is NOT the stored per-cell training statistic for the four panel cells), channel 1
subtracts the panel-mean ATAC z, channel 2 is the standardised log1p panel-mean p300. The panel
is leave-one-out: for K562, which the model never saw, the mean is over all four training
cells, exactly as the deployment arm does; for a training cell it is over the other three,
exactly as training did. Channels 1 and 2 use each FOLD's own panel_mean_stats.json, as 4.1
does, so the accessibility input is rebuilt per (context, fold).

FOLDS. All five folds of the seed-0 model are averaged. 4.1's held-out-fold-per-chromosome rule
does not apply: the DNA is shuffled and carries no genomic identity, and the real accessibility
channels are identical across the four terms of the contrast.

NO REVERSE-COMPLEMENT AVERAGING. 4.1 averages the two strands, but the reverse complement of a
designed window maps both motifs to the other orientation and flips the sign of d, so averaging
it in would mix design cells that this assay is built to separate.

Stages:
  prep     choose backgrounds, compute per-cell accessibility normalisation (CPU)
  run      predict one chunk of backgrounds, write partial sums (GPU)
  combine  pool the partials, t-test and BH-FDR, write the result tables (CPU)
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
import torch
from scipy import stats

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
A = f"{D}/2026_0824_H3K27ac_model"
sys.path.insert(0, f"{D}/scripts")
from train_multimodal_bpnet import (  # noqa: E402
    extract_windows, normalize_accessibility, panel_mean_tracks)

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
MEME = f"{D}/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
MODEL = f"{A}/models/p300_multicell_panelmean_holdK562_hw500_clw10"
PANEL_CFG = f"{A}/config/p300_panel.json"
REG = ("/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/"
       "results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/"
       "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
SIGNAL_BOUNDS = f"{A}/data/h3k27ac_5p_plus.bw"   # values discarded; bounds check only
IN_W, OUT_W = 2114, 1000
CONTEXTS = ["K562", "A549", "GM12878", "HepG2", "MCF-7"]

# label -> (JASPAR matrix id, column slice or None). The consensus is the max-probability base
# at each kept column.
ANCHORS = {
    "AP1_FOSJUN": ("MA0099.4", None),
    "GATA1_TAL1": ("MA0140.3", None),
    "GATA1": ("MA0035.5", None),
    # MA0035.5's max-probability consensus is CTAATCT, i.e. AGATTAG on the other strand, one
    # base off the canonical WGATAA core (position 2 is A 0.692 / T 0.286). MA0140.3's first
    # six columns are the same half-site read cleanly (TTATCT = AGATAA), so this fourth anchor
    # carries the positive control's GATA site with that base taken from a JASPAR matrix rather
    # than written by hand. Drop it with --anchors to get exactly the three asked for.
    "GATA1_WGATAA": ("MA0140.3", (0, 6)),
}
PARTNERS = {
    "FOXA1": ("MA0148.5", None),
    "TEAD1": ("MA0090.4", None),
    "ETS1": ("MA0098.4", None),
    "SPI1": ("MA0080.7", None),
    "NFIA": ("MA0670.2", None),
    "CEBPB": ("MA0466.4", None),
    "GATA1": ("MA0035.5", None),
    "EBOX_TAL1_TCF3": ("MA0091.2", None),
    "RUNX1": ("MA0002.3", None),
    "IRF4": ("MA1419.2", None),
    "EGR1": ("MA0162.5", None),
    "MEF2C": ("MA0497.2", None),
    "CTCF": ("MA0139.2", None),
}

_LUT = np.zeros((256, 4), dtype=np.float32)
for _i, _c in enumerate("ACGT"):
    _LUT[ord(_c), _i] = 1.0
    _LUT[ord(_c.lower()), _i] = 1.0
_COMP = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def read_meme_consensus():
    """Probability matrix of every motif in the JASPAR meme file, by matrix id."""
    lines = open(MEME).read().splitlines()
    out, i = {}, 0
    while i < len(lines):
        if not lines[i].startswith("MOTIF "):
            i += 1
            continue
        mid = lines[i].split()[1]
        j = i + 1
        while j < len(lines) and not lines[j].startswith("letter-probability"):
            j += 1
        w = int(re.search(r"w=\s*(\d+)", lines[j]).group(1))
        rows = [[float(x) for x in lines[j + 1 + k].split()] for k in range(w)]
        out.setdefault(mid, rows)          # first occurrence of an id wins
        i = j + w + 1
    return out


def consensus(mats, spec):
    mid, sl = spec
    if mid not in mats:
        raise SystemExit(f"error: {mid} not in {MEME}")
    rows = mats[mid]
    if sl is not None:
        rows = rows[sl[0]:sl[1]]
    return "".join("ACGT"[int(np.argmax(r))] for r in rows)


def dinuc_shuffle(seq, rng):
    """Dinucleotide-frequency-preserving shuffle (Altschul-Erikson).

    Same algorithm as motif_exp_utils.dinuc_shuffle, reimplemented because that module imports
    tensorflow at the top level: the `multimodal` env has no tensorflow, and pulling a second
    deep-learning runtime into a torch job would be wrong even if it did.
    """
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.int8)
    chars, tokens = np.unique(arr, return_inverse=True)
    nxt = []
    for t in range(len(chars)):
        inds = np.where(tokens[:-1] == t)[0] + 1
        if len(inds) > 1:
            order = np.arange(len(inds))
            order[:-1] = rng.permutation(len(inds) - 1)   # last stays put: Eulerian validity
            inds = inds[order]
        nxt.append(inds)
    counters = [0] * len(chars)
    ind = 0
    result = np.empty_like(tokens)
    result[0] = tokens[ind]
    for j in range(1, len(tokens)):
        t = tokens[ind]
        ind = nxt[t][counters[t]]
        counters[t] += 1
        result[j] = tokens[ind]
    return chars[result].astype(np.int8).tobytes().decode("ascii")


def build_design(anchor_names, partner_names, spacings, mats):
    """Every sequence the assay needs, de-duplicated, plus the row index arrays.

    Positions depend only on motif lengths, never on the background, so this is built once and
    reused for every background, context and fold. `partner alone` is keyed on the partner's
    ABSOLUTE start, so the overlapping position ranges that different anchor lengths produce
    share sequences instead of duplicating them (~2,400 instead of ~8,400 sequences).
    """
    seqs = {}
    for n in anchor_names:
        if n not in ANCHORS:
            raise SystemExit(f"error: unknown anchor {n}")
        seqs[n] = consensus(mats, ANCHORS[n])
    for n in partner_names:
        if n not in PARTNERS:
            raise SystemExit(f"error: unknown partner {n}")
        seqs[n] = consensus(mats, PARTNERS[n])
    centre = IN_W // 2

    specs, index = [], {}

    def spec_of(ins):
        key = tuple(sorted(ins))
        if key not in index:
            index[key] = len(specs)
            specs.append(key)
        return index[key]

    none_ix = spec_of(())
    anchor_ix, anchor_pos = {}, {}
    for an in anchor_names:
        s = seqs[an]
        a_s = centre - len(s) // 2
        anchor_pos[an] = (a_s, a_s + len(s))
        anchor_ix[an] = spec_of(((s, a_s),))

    rec = []
    for an in anchor_names:
        a_s, a_e = anchor_pos[an]
        aseq = seqs[an]
        for pn in partner_names:
            base = seqs[pn]
            for orient in ("+", "-"):
                pseq = base if orient == "+" else revcomp(base)
                for d in spacings:
                    # d >= 0: partner downstream of the anchor with an edge-to-edge gap of d.
                    # d < 0: partner upstream with a gap of |d|. d = 0 is abutting downstream.
                    p_s = a_e + d if d >= 0 else a_s + d - len(pseq)
                    if p_s < 0 or p_s + len(pseq) > IN_W:
                        raise SystemExit(f"error: {an}/{pn} d={d} falls outside the window")
                    if not (p_s >= a_e or p_s + len(pseq) <= a_s):
                        raise SystemExit(f"error: {an}/{pn} d={d} overlaps the anchor")
                    rec.append((an, pn, orient, d,
                                spec_of(((aseq, a_s), (pseq, p_s))),
                                spec_of(((pseq, p_s),)),
                                anchor_ix[an]))
    rows = pd.DataFrame(rec, columns=["anchor", "partner", "orient", "d",
                                      "pair_ix", "palone_ix", "anchor_ix"])
    return specs, rows, {"none_ix": none_ix, "anchor_ix": anchor_ix, "motifs": seqs}


def materialize(specs, base_seq):
    """(n_specs, IN_W) uint8 ASCII matrix: the background with each spec's inserts written in."""
    base = np.frombuffer(base_seq.encode("ascii"), dtype=np.uint8)
    if len(base) != IN_W:
        raise SystemExit(f"error: background is {len(base)} bp, expected {IN_W}")
    toks = np.repeat(base[None, :], len(specs), axis=0)
    for si, ins in enumerate(specs):
        for mseq, st in ins:
            toks[si, st:st + len(mseq)] = np.frombuffer(mseq.encode("ascii"), dtype=np.uint8)
    return toks


# --------------------------------------------------------------------------- prep

def stage_prep(a):
    panel = json.load(open(PANEL_CFG))
    bw0 = pyBigWig.open(panel["K562"]["accessibility_bw"])
    chroms = dict(bw0.chroms())
    bw0.close()
    reg = pd.read_csv(REG, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chr", "start", "end"])
    reg = reg[reg["chr"].isin(chroms)].reset_index(drop=True)
    print(f"{len(reg):,} ABC candidate regions", flush=True)

    rng = np.random.default_rng(a.seed)
    pick = np.sort(rng.choice(len(reg), size=min(a.stat_sample, len(reg)), replace=False))
    sub = reg.iloc[pick].reset_index(drop=True)
    sub["summit"] = ((sub["end"] - sub["start"]) // 2).astype(int)
    sub["centre"] = sub["start"] + sub["summit"]
    half = IN_W // 2

    accstats, per_region = {}, {}
    for cell in CONTEXTS:
        bw = pyBigWig.open(panel[cell]["accessibility_bw"])
        cl = dict(bw.chroms())
        n = tot = totsq = 0
        means = np.full(len(sub), np.nan)
        for i, row in sub.iterrows():
            s, e = int(row["centre"]) - half, int(row["centre"]) + half
            if s < 0 or e > cl.get(row["chr"], 0):
                continue
            v = bw.values(row["chr"], s, e, numpy=True)
            if v is None or len(v) != IN_W:
                continue
            lv = np.log1p(np.clip(np.nan_to_num(v, nan=0.0), 0, None))
            n += lv.size
            tot += float(lv.sum())
            totsq += float((lv ** 2).sum())
            means[i] = float(lv.mean())
        bw.close()
        mean = tot / n
        std = float(np.sqrt(max(totsq / n - mean ** 2, 0.0)))
        accstats[cell] = {"acc_mean": mean, "acc_std": std, "n_windows": int(n // IN_W)}
        per_region[cell] = means
        print(f"  {cell:8s} log1p ATAC over {n // IN_W:,} windows: "
              f"mean {mean:.4f} std {std:.4f}", flush=True)

    # BACKGROUND RULE: middle tercile of the window-mean log1p K562 ATAC over the sampled ABC
    # candidate regions. Accessible but not the most accessible, and chosen with no reference to
    # observed or predicted p300, so nothing about the label enters the background set.
    k = per_region["K562"]
    ok = np.isfinite(k)
    for cell in CONTEXTS:
        ok &= np.isfinite(per_region[cell])
    lo, hi = np.quantile(k[ok], [1 / 3, 2 / 3])
    pool = np.where(ok & (k > lo) & (k <= hi))[0]
    print(f"\nmiddle K562 ATAC tercile: {lo:.4f} < mean log1p <= {hi:.4f}, "
          f"{len(pool):,} candidates", flush=True)

    order = rng.permutation(pool)[:4 * a.n_bg]
    cand = sub.iloc[order].reset_index(drop=True)
    cand["atac_mean_log1p"] = k[order]

    # Keep only windows every cell's tracks can serve: panel_mean_tracks refuses a region any
    # panel cell cannot serve, and a row dropped in one cell but not another would pair
    # predictions with the wrong window.
    valid = np.ones(len(cand), dtype=bool)
    for cell in CONTEXTS:
        spec = panel[cell]
        tracks = ",".join([spec["signal_plus_bw"]]
                          + ([spec["signal_minus_bw"]] if spec.get("signal_minus_bw") else [])
                          + [spec["accessibility_bw"]])
        _, _, _, v = extract_windows(cand, None, spec["signal_plus_bw"],
                                     spec.get("signal_minus_bw"), tracks,
                                     IN_W, OUT_W, 0, is_peak=True)
        valid &= v
        print(f"  {cell:8s} valid {int(v.sum()):,}/{len(cand):,}", flush=True)
    _, _, _, vseq = extract_windows(cand, GEN, SIGNAL_BOUNDS, None, None,
                                    IN_W, OUT_W, 0, is_peak=True)
    valid &= vseq
    print(f"  sequence valid {int(vseq.sum()):,}/{len(cand):,}; "
          f"valid everywhere {int(valid.sum()):,}", flush=True)

    bg = cand[valid].reset_index(drop=True).iloc[:a.n_bg]
    if len(bg) < a.n_bg:
        raise SystemExit(f"error: only {len(bg)} valid backgrounds for --n-bg {a.n_bg}; "
                         f"raise --stat-sample")
    bg[["chr", "start", "end", "summit", "atac_mean_log1p"]].to_csv(
        a.bg_tsv, sep="\t", index=False)
    json.dump({"acc_stats_prediction_regions": accstats,
               "rule": ("middle tercile of window-mean log1p K562 ATAC over a "
                        f"{a.stat_sample}-region sample of the ABC candidate set (seed "
                        f"{a.seed}); label-free, no reference to p300"),
               "tercile_bounds": [float(lo), float(hi)], "n_backgrounds": int(len(bg)),
               "stat_sample": int(a.stat_sample), "regions": REG},
              open(a.stats_json, "w"), indent=2)
    print(f"\nwrote {a.bg_tsv} ({len(bg)} backgrounds) and {a.stats_json}", flush=True)


# --------------------------------------------------------------------------- run

def stage_run(a):
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.backends.cudnn.benchmark = True
    if not a.tf32:
        # TF32 keeps ~10 mantissa bits, so a log count near 4.5 carries ~5e-3 of rounding --
        # the same size as the cooperativity being measured. Off by default.
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.cuda.matmul.allow_tf32 = False
    print(f"device {dev}, tf32 {a.tf32}", flush=True)

    panel_full = json.load(open(PANEL_CFG))
    accstats = json.load(open(a.stats_json))["acc_stats_prediction_regions"]
    contexts = a.contexts.split(",")
    anchors = a.anchors.split(",")
    partners = a.partners.split(",")
    folds = a.folds.split(",")
    spacings = list(range(a.dmin, a.dmax + 1))
    for c in contexts:
        if c not in panel_full:
            raise SystemExit(f"error: unknown context {c}")

    mats = read_meme_consensus()
    specs, rows, meta = build_design(anchors, partners, spacings, mats)
    print(f"{len(rows):,} (anchor, partner, orient, d) rows; {len(specs):,} unique sequences "
          f"per background", flush=True)
    for n, s in meta["motifs"].items():
        print(f"  {n:16s} {s}", flush=True)

    bg = pd.read_csv(a.bg_tsv, sep="\t")
    start = a.chunk_start
    stop = min(a.chunk_start + a.chunk_size, len(bg))
    chunk = bg.iloc[start:stop].reset_index(drop=True)
    if not len(chunk):
        raise SystemExit(f"error: empty chunk at {start} of {len(bg)} backgrounds")
    print(f"backgrounds {start}..{stop - 1} of {len(bg)}", flush=True)

    models, pstats = {}, {}
    for f in folds:
        md = f"{MODEL}/fold{f}"
        if not os.path.exists(f"{md}/training_complete.json"):
            raise SystemExit(f"error: {md} has no training_complete.json")
        m = torch.load(f"{md}/multimodal_bpnet.torch", map_location="cpu", weights_only=False)
        if m.mode != "multimodal" or m.n_acc_channels != 3:
            raise SystemExit(f"error: fold{f} is mode {m.mode}, {m.n_acc_channels} acc channels")
        if m.trimming * 2 + OUT_W != IN_W:
            raise SystemExit(f"error: fold{f} trimming {m.trimming} does not give {IN_W}")
        models[f] = m.to(dev).eval()
        pstats[f] = json.load(open(f"{md}/panel_mean_stats.json"))
    nf = len(models)
    training_cells = pstats[folds[0]]["training_cells"]

    # Accessibility per (context, fold): channel 0 is fold-independent, but channel 1's
    # subtrahend and channel 2's standardisation come from each fold's own panel_mean_stats,
    # exactly as 4.1 reads them per fold.
    acc = {}
    for ctx in contexts:
        _, _, raw, valid = extract_windows(chunk, None, SIGNAL_BOUNDS, None,
                                           panel_full[ctx]["accessibility_bw"],
                                           IN_W, OUT_W, 0, is_peak=True)
        if not valid.all():
            raise SystemExit(f"error: {int((~valid).sum())} windows invalid in {ctx}'s ATAC")
        x0 = normalize_accessibility(raw, mean=accstats[ctx]["acc_mean"],
                                     std=accstats[ctx]["acc_std"])[0]
        others = {c: panel_full[c] for c in training_cells if c != ctx}
        for f in folds:
            pst = pstats[f]
            zm, pl, om = panel_mean_tracks(chunk, 0, True, others, pst["acc_stats"],
                                           pst["p300_depth"], pst["dref"], IN_W, OUT_W)
            pl = (pl - pst["p300_channel_mean"]) / pst["p300_channel_std"]
            acc[(ctx, f)] = np.concatenate([x0, x0 - zm, pl], axis=1).astype(np.float32)
        print(f"  {ctx:8s} panel mean over {', '.join(others)}; fold{folds[-1]} channel sds "
              f"{np.round(acc[(ctx, folds[-1])].std(axis=(0, 2)), 3)}; "
              f"offset mean {np.log1p(om).mean():.3f}", flush=True)

    genome = pyfaidx.Fasta(GEN)
    half = IN_W // 2
    nr, nc = len(rows), len(contexts)
    pair_ix = rows["pair_ix"].to_numpy()
    palone_ix = rows["palone_ix"].to_numpy()
    anc_ix = rows["anchor_ix"].to_numpy()
    none_ix = meta["none_ix"]

    S = np.zeros((nc, nr))                  # sum of coop over backgrounds
    SS = np.zeros((nc, nr))                 # sum of squares
    MA = np.zeros((nc, len(anchors)))       # anchor alone minus neither
    MP = np.zeros((nc, nr))                 # partner alone at this position minus neither
    NB = 0
    dp_ix = None
    SD = SSD = None
    if "K562" in contexts and "HepG2" in contexts:
        dp_ix = (contexts.index("K562"), contexts.index("HepG2"))
        SD, SSD = np.zeros(nr), np.zeros(nr)

    t0, nseq = time.time(), 0
    for bi in range(len(chunk)):
        row = chunk.iloc[bi]
        gi = start + bi
        centre = int(row["start"]) + int(row["summit"])
        real = str(genome[row["chr"]][centre - half:centre + half]).upper()
        shuf = dinuc_shuffle(real, np.random.RandomState(a.seed * 1000003 + gi))
        toks = materialize(specs, shuf)

        pred = np.zeros((nc, len(specs)))
        with torch.inference_mode():
            for ci, ctx in enumerate(contexts):
                for f in folds:
                    m = models[f]
                    at = torch.from_numpy(acc[(ctx, f)][bi]).to(dev)
                    for i in range(0, len(specs), a.batch):
                        oh = np.ascontiguousarray(
                            _LUT[toks[i:i + a.batch]].transpose(0, 2, 1))
                        xb = torch.from_numpy(oh).to(dev)
                        xb = torch.cat(
                            [xb, at.unsqueeze(0).expand(xb.shape[0], -1, -1)], dim=1)
                        _, lc = m(xb)
                        pred[ci, i:i + a.batch] += lc.squeeze(-1).double().cpu().numpy()
                    nseq += len(specs)
                pred[ci] /= nf

        base = pred[:, none_ix][:, None]
        coop = pred[:, pair_ix] - pred[:, anc_ix] - pred[:, palone_ix] + base
        S += coop
        SS += coop ** 2
        MP += pred[:, palone_ix] - base
        for ai, an in enumerate(anchors):
            MA[:, ai] += pred[:, meta["anchor_ix"][an]] - pred[:, none_ix]
        if dp_ix:
            dv = coop[dp_ix[0]] - coop[dp_ix[1]]
            SD += dv
            SSD += dv ** 2
        NB += 1
        if NB % a.report == 0 or NB == len(chunk):
            el = time.time() - t0
            print(f"  {NB}/{len(chunk)} backgrounds, {nseq:,} forwards, "
                  f"{nseq / el:,.0f}/s, {el / 60:.1f} min elapsed, "
                  f"{el / NB * (len(chunk) - NB) / 60:.1f} min left", flush=True)
    genome.close()

    os.makedirs(os.path.dirname(a.part), exist_ok=True)
    out = dict(S=S, SS=SS, MA=MA, MP=MP, n=np.array([NB]),
               contexts=np.array(contexts), anchors=np.array(anchors),
               partners=np.array(partners), spacings=np.array(spacings))
    if dp_ix:
        out.update(SD=SD, SSD=SSD)
    np.savez_compressed(a.part, **out)
    print(f"wrote {a.part}: {NB} backgrounds, {nseq:,} forwards in "
          f"{(time.time() - t0) / 60:.1f} min", flush=True)


# --------------------------------------------------------------------------- combine

def _tstat(S, SS, n):
    """One-sample t-test from running sums. S, SS shaped (..., n_rows)."""
    mean = S / n
    var = np.clip((SS - n * mean ** 2) / max(n - 1, 1), 0, None)
    se = np.sqrt(var / n)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, mean / se, 0.0)
    return mean, se, 2 * stats.t.sf(np.abs(t), df=max(n - 1, 1))


def _bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def stage_combine(a):
    parts = sorted(f for f in os.listdir(a.part_dir) if f.endswith(".npz"))
    if not parts:
        raise SystemExit(f"error: no partials in {a.part_dir}")
    acc = None
    for f in parts:
        z = np.load(os.path.join(a.part_dir, f), allow_pickle=False)
        if acc is None:
            acc = {k: z[k].copy() for k in z.files}
        else:
            for k in ("S", "SS", "MA", "MP", "n", "SD", "SSD"):
                if k in z.files:
                    acc[k] = acc[k] + z[k]
            for k in ("contexts", "anchors", "partners", "spacings"):
                if not np.array_equal(acc[k], z[k]):
                    raise SystemExit(f"error: {f} disagrees with the first partial on {k}")
        print(f"  {f}: n={int(z['n'][0])}", flush=True)
    n = int(acc["n"][0])
    contexts = [str(c) for c in acc["contexts"]]
    anchors = [str(c) for c in acc["anchors"]]
    partners = [str(c) for c in acc["partners"]]
    spacings = [int(c) for c in acc["spacings"]]
    print(f"{len(parts)} partials, {n} backgrounds, contexts {contexts}", flush=True)

    _, rows, _ = build_design(anchors, partners, spacings, read_meme_consensus())
    if len(rows) != acc["S"].shape[1]:
        raise SystemExit(f"error: design has {len(rows)} rows, partials have "
                         f"{acc['S'].shape[1]}; the flags must match the run")

    mean, se, p = _tstat(acc["S"], acc["SS"], n)
    frames = []
    for ci, ctx in enumerate(contexts):
        f = rows[["anchor", "partner", "orient", "d"]].copy()
        f["context"] = ctx
        f["mean_coop"] = mean[ci]
        f["se"] = se[ci]
        f["p"] = p[ci]
        frames.append(f)
    res = pd.concat(frames, ignore_index=True)
    res["fdr"] = _bh(res["p"].to_numpy())
    res["n_bg"] = n
    res = res[["anchor", "partner", "orient", "d", "context",
               "mean_coop", "se", "p", "fdr", "n_bg"]]
    res.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    sz = os.path.getsize(a.out) / 1e6
    print(f"wrote {a.out} ({len(res):,} rows, {sz:.2f} MB)", flush=True)
    if sz > 5:
        print("WARNING: over the repo's 5 MB git limit", flush=True)

    # Marginal single-motif effects, so the cooperativity can be read against them.
    mrec = []
    for ci, ctx in enumerate(contexts):
        for ai, an in enumerate(anchors):
            mrec.append((an, an, "anchor_alone", "+", 0, ctx, acc["MA"][ci, ai] / n))
    for ci, ctx in enumerate(contexts):
        for j in range(len(rows)):
            r = rows.iloc[j]
            mrec.append((r["anchor"], r["partner"], "partner_alone", r["orient"], r["d"],
                         ctx, acc["MP"][ci, j] / n))
    marg = pd.DataFrame(mrec, columns=["anchor", "motif", "role", "orient", "d", "context",
                                       "mean_effect"])
    mpath = a.out.replace(".tsv", ".marginals.tsv")
    marg.to_csv(mpath, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {mpath} ({len(marg):,} rows, "
          f"{os.path.getsize(mpath) / 1e6:.2f} MB)", flush=True)

    if "SD" in acc:
        dm, dse, dpv = _tstat(acc["SD"][None, :], acc["SSD"][None, :], n)
        dif = rows[["anchor", "partner", "orient", "d"]].copy()
        dif["mean_diff_K562_minus_HepG2"] = dm[0]
        dif["se"] = dse[0]
        dif["p"] = dpv[0]
        dif["fdr"] = _bh(dpv[0])
        dif["n_bg"] = n
        dpath = a.out.replace(".tsv", ".context_diff.tsv")
        dif.to_csv(dpath, sep="\t", index=False, float_format="%.6g")
        print(f"wrote {dpath} ({len(dif):,} rows)", flush=True)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, choices=["prep", "run", "combine"])
    ap.add_argument("--bg-tsv", default=f"{A}/results/2026_1020_motif_pair_backgrounds.tsv")
    ap.add_argument("--stats-json", default=f"{A}/results/2026_1020_motif_pair_accstats.json")
    ap.add_argument("--out", default=f"{A}/results/2026_1020_motif_pair_spacing.tsv")
    ap.add_argument("--part", default=None, help="run: partial .npz to write")
    ap.add_argument("--part-dir", default=None, help="combine: directory of partial .npz files")
    ap.add_argument("--n-bg", type=int, default=500)
    ap.add_argument("--stat-sample", type=int, default=20000,
                    help="prep: ABC regions sampled for the accessibility statistics and terciles")
    ap.add_argument("--chunk-start", type=int, default=0)
    ap.add_argument("--chunk-size", type=int, default=50)
    ap.add_argument("--contexts", default=",".join(CONTEXTS))
    ap.add_argument("--anchors", default=",".join(ANCHORS))
    ap.add_argument("--partners", default=",".join(PARTNERS))
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--dmin", type=int, default=-40)
    ap.add_argument("--dmax", type=int, default=40)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--report", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tf32", action="store_true",
                    help="allow TF32 convolutions: faster, ~5e-3 of rounding on each log count")
    return ap.parse_args()


def main():
    a = parse_args()
    if a.stage == "prep":
        stage_prep(a)
    elif a.stage == "run":
        if not a.part:
            raise SystemExit("error: --stage run needs --part")
        stage_run(a)
    else:
        if not a.part_dir:
            raise SystemExit("error: --stage combine needs --part-dir")
        stage_combine(a)


if __name__ == "__main__":
    main()
