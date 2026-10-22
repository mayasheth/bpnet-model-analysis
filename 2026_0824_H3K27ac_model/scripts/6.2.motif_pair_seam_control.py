#!/usr/bin/env python
"""Seam control for the motif-pair insertion assay (6.1, F-068).

F-068's leading pairs (AP1-IRF4 d = -2 and +2, AP1-ETS1 d = 0) are spikes at |d| <= 2, about
ten times the same curve at d = 3-8. Two abutting consensus sequences can form a new k-mer at
their junction, so a seam spike may be a junction artefact rather than motif-pair syntax. This
script reruns the selected seam rows with two controls on the same backgrounds, contexts and
folds, design fixed in decisions [2026-10-21] "Seam control for F-068".

SELECTION. From 6.1's K562 table: for each (anchor, partner, orientation) curve, the |d| <= 2
row with the largest |coop|; kept if spike = coop - median(coop at 3 <= |d| <= 8) exceeds 0.01
in absolute value. The selection is recomputed from the table at every stage, so all stages
agree by construction.

ARMS, per background b and selected row:
  real      pair - anchor - partner + neither, exactly as 6.1 (and gated against 6.1, below)
  scramble  the partner replaced by a mononucleotide permutation of itself at the same
            position, K fresh permutations per background (never the partner or its reverse
            complement); coop averaged over the K. Asks whether any composition-matched
            sequence abutting the anchor gives the spike.
  jitter    the 2 junction-side bases of the anchor and of the partner redrawn per background
            from their PWM columns (oriented as inserted), used in the pair AND in each motif's
            alone term. Asks whether the spike needs those exact junction bases.

GATE. With 25 backgrounds per block, every block's real coop is compared with the partial sums
6.1 wrote for the same 25 backgrounds (6.1 used 25 per array task). A mismatch above --gate-tol
stops the run: it means this script does not reproduce the assay it is controlling.
--gate-shift N compares against block + N instead; it exists to show the gate fails.

Stages:
  run      GPU, a range of backgrounds, one gated partial per 25-background block
  combine  CPU, pools the partials and applies the pre-registered verdicts
"""
import argparse
import importlib.util
import os
import time

import numpy as np
import pandas as pd
import pyfaidx
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "mp61", os.path.join(HERE, "6.1.motif_pair_spacing_multimodal.py"))
mp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mp)

A = mp.A
BLOCK = 25
JUNC = 2
_COMP_COL = [3, 2, 1, 0]          # A C G T -> T G C A


def select_rows(ref_tsv, min_spike):
    t = pd.read_csv(ref_tsv, sep="\t")
    key = ["anchor", "partner", "orient"]
    flank = (t[(t.d.abs() >= 3) & (t.d.abs() <= 8)]
             .groupby(key + ["context"]).mean_coop.median().rename("flank_med"))
    k = t[(t.context == "K562") & (t.d.abs() <= 2)]
    seam = k.loc[k.groupby(key).mean_coop.apply(lambda s: s.abs().idxmax())]
    seam = seam.merge(flank.reset_index(), on=key + ["context"])
    seam["spike"] = seam.mean_coop - seam.flank_med
    sel = seam[seam.spike.abs() > min_spike][key + ["d"]].reset_index(drop=True)
    long = t.merge(sel, on=key + ["d"]).merge(flank.reset_index(), on=key + ["context"])
    long["spike"] = long.mean_coop - long.flank_med
    long = long.rename(columns={"mean_coop": "ref_coop", "se": "ref_se"})
    return sel, long[key + ["d", "context", "ref_coop", "ref_se", "flank_med", "spike"]]


def motif_cols(mats, spec):
    mid, sl = spec
    rows = np.asarray(mats[mid], float)
    return rows if sl is None else rows[sl[0]:sl[1]]


def orient_cols(cols, orient):
    return cols if orient == "+" else cols[::-1][:, _COMP_COL]


def placements(sel, mats):
    """Positions, sequences and junction columns of every selected row (6.1's geometry)."""
    centre = mp.IN_W // 2
    out = []
    for r in sel.itertuples(index=False):
        aseq = mp.consensus(mats, mp.ANCHORS[r.anchor])
        acol = motif_cols(mats, mp.ANCHORS[r.anchor])
        a_s = centre - len(aseq) // 2
        a_e = a_s + len(aseq)
        base = mp.consensus(mats, mp.PARTNERS[r.partner])
        pseq = base if r.orient == "+" else mp.revcomp(base)
        pcol = orient_cols(motif_cols(mats, mp.PARTNERS[r.partner]), r.orient)
        if "".join("ACGT"[i] for i in pcol.argmax(1)) != pseq:
            raise SystemExit(f"error: oriented PWM of {r.partner} {r.orient} != {pseq}")
        p_s = a_e + r.d if r.d >= 0 else a_s + r.d - len(pseq)
        if r.d >= 0:
            aj, pj = list(range(len(aseq) - JUNC, len(aseq))), list(range(JUNC))
        else:
            aj, pj = list(range(JUNC)), list(range(len(pseq) - JUNC, len(pseq)))
        out.append(dict(aseq=aseq, acol=acol, a_s=a_s, pseq=pseq, pcol=pcol, p_s=p_s,
                        aj=aj, pj=pj))
    return out


def scramble(seq, rng):
    bad = {seq, mp.revcomp(seq)}
    for _ in range(1000):
        s = "".join(rng.permutation(list(seq)))
        if s not in bad:
            return s
    raise SystemExit(f"error: cannot scramble {seq}")


def jitter(seq, cols, idx, rng):
    s = list(seq)
    for i in idx:
        p = cols[i] / cols[i].sum()
        s[i] = "ACGT"[rng.choice(4, p=p)]
    s = "".join(s)
    return s, sum(a != b for a, b in zip(s, seq))


def ref_block(ref_parts, bi, contexts, sel):
    """6.1's mean real coop for the selected rows over one 25-background block."""
    z = np.load(f"{ref_parts}/part_{bi:03d}.npz", allow_pickle=False)
    anchors = [str(x) for x in z["anchors"]]
    partners = [str(x) for x in z["partners"]]
    spacings = [int(x) for x in z["spacings"]]
    _, rows, _ = mp.build_design(anchors, partners, spacings, mp.read_meme_consensus())
    ix = {k: i for i, k in enumerate(zip(rows.anchor, rows.partner, rows.orient, rows.d))}
    cix = [[str(c) for c in z["contexts"]].index(c) for c in contexts]
    j = [ix[k] for k in zip(sel.anchor, sel.partner, sel.orient, sel.d)]
    return z["S"][np.ix_(cix, j)] / int(z["n"][0]), int(z["n"][0])


def stage_run(a):
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cuda.matmul.allow_tf32 = False
    contexts = a.contexts.split(",")
    folds = a.folds.split(",")
    mats = mp.read_meme_consensus()
    sel, _ = select_rows(a.ref_tsv, a.min_spike)
    if a.expect_n and len(sel) != a.expect_n:
        raise SystemExit(f"error: selected {len(sel)} rows, pre-registered {a.expect_n}")
    pl = placements(sel, mats)
    nr, nc = len(sel), len(contexts)
    print(f"device {dev}; {nr} seam rows; contexts {contexts}; K = {a.k}", flush=True)

    panel_full = mp.json.load(open(mp.PANEL_CFG))
    accstats = mp.json.load(open(a.stats_json))["acc_stats_prediction_regions"]
    models, pstats = mp.load_models(folds, dev)
    nf = len(models)
    bg = pd.read_csv(a.bg_tsv, sep="\t")
    genome = pyfaidx.Fasta(mp.GEN)
    half = mp.IN_W // 2
    os.makedirs(a.part_dir, exist_ok=True)

    if a.start % BLOCK:
        raise SystemExit(f"error: --start must be a multiple of {BLOCK}")
    for b0 in range(a.start, min(a.stop, len(bg)), BLOCK):
        bi_ref = b0 // BLOCK
        part = f"{a.part_dir}/part_{bi_ref:03d}.npz"
        if os.path.exists(part):
            print(f"already done: {part}", flush=True)
            continue
        chunk = bg.iloc[b0:b0 + BLOCK].reset_index(drop=True)
        acc = mp.build_acc(chunk, contexts, folds, pstats, panel_full, accstats)
        acc_s = {k: np.zeros((nc, nr)) for k in
                 ("S_real", "SS_real", "S_scr", "SS_scr", "S_d1", "SS_d1", "S_jit", "SS_jit")}
        nchg = np.zeros(nr)
        anychg = np.zeros(nr)
        t0 = time.time()
        for bi in range(len(chunk)):
            row = chunk.iloc[bi]
            gi = b0 + bi
            c0 = int(row["start"]) + int(row["summit"])
            real = str(genome[row["chr"]][c0 - half:c0 + half]).upper()
            shuf = mp.dinuc_shuffle(real, np.random.RandomState(a.seed * 1000003 + gi))
            rng = np.random.RandomState(a.ctrl_seed * 1000003 + gi)

            specs, index = [], {}

            def sid(*ins):
                key = tuple(sorted(ins))
                if key not in index:
                    index[key] = len(specs)
                    specs.append(key)
                return index[key]

            none = sid()
            terms = []
            for ri, p in enumerate(pl):
                A_, P_ = (p["aseq"], p["a_s"]), (p["pseq"], p["p_s"])
                real_t = (sid(A_, P_), sid(A_), sid(P_))
                scr_t = []
                for _ in range(a.k):
                    Pk = (scramble(p["pseq"], rng), p["p_s"])
                    scr_t.append((sid(A_, Pk), sid(Pk)))
                aj, na = jitter(p["aseq"], p["acol"], p["aj"], rng)
                pj, npj = jitter(p["pseq"], p["pcol"], p["pj"], rng)
                Aj, Pj = (aj, p["a_s"]), (pj, p["p_s"])
                jit_t = (sid(Aj, Pj), sid(Aj), sid(Pj))
                nchg[ri] += na + npj
                anychg[ri] += (na + npj) > 0
                terms.append((real_t, scr_t, jit_t))

            toks = mp.materialize(specs, shuf)
            pred = np.zeros((nc, len(specs)))
            with torch.inference_mode():
                for ci, ctx in enumerate(contexts):
                    for f in folds:
                        at = torch.from_numpy(acc[(ctx, f)][bi]).to(dev)
                        for i in range(0, len(specs), a.batch):
                            oh = np.ascontiguousarray(
                                mp._LUT[toks[i:i + a.batch]].transpose(0, 2, 1))
                            xb = torch.from_numpy(oh).to(dev)
                            xb = torch.cat(
                                [xb, at.unsqueeze(0).expand(xb.shape[0], -1, -1)], dim=1)
                            _, lc = models[f](xb)
                            pred[ci, i:i + a.batch] += lc.squeeze(-1).double().cpu().numpy()
                    pred[ci] /= nf

            nb = pred[:, none]
            for ri, (rt, st, jt) in enumerate(terms):
                cr = pred[:, rt[0]] - pred[:, rt[1]] - pred[:, rt[2]] + nb
                cs = np.mean([pred[:, pp] - pred[:, rt[1]] - pred[:, pa] + nb
                              for pp, pa in st], axis=0)
                cj = pred[:, jt[0]] - pred[:, jt[1]] - pred[:, jt[2]] + nb
                for nm, v in (("real", cr), ("scr", cs), ("d1", cr - cs), ("jit", cj)):
                    acc_s[f"S_{nm}"][:, ri] += v
                    acc_s[f"SS_{nm}"][:, ri] += v ** 2
        n = len(chunk)

        ref, nref = ref_block(a.ref_parts, bi_ref + a.gate_shift, contexts, sel)
        if nref != n:
            raise SystemExit(f"error: 6.1 block {bi_ref + a.gate_shift} has {nref} "
                             f"backgrounds, this block {n}")
        diff = np.abs(acc_s["S_real"] / n - ref).max()
        print(f"block {bi_ref}: {len(specs)} sequences/background, "
              f"{(time.time() - t0) / 60:.1f} min; gate max |real - 6.1| = {diff:.2e}",
              flush=True)
        if not diff <= a.gate_tol:
            raise SystemExit(f"GATE FAILED: block {bi_ref} real coop differs from 6.1 block "
                             f"{bi_ref + a.gate_shift} by {diff:.3e} > {a.gate_tol}")
        np.savez_compressed(part, n=np.array([n]), nchg=nchg, anychg=anychg,
                            contexts=np.array(contexts),
                            rows=np.array([f"{r.anchor}|{r.partner}|{r.orient}|{r.d}"
                                           for r in sel.itertuples(index=False)]),
                            gate=np.array([diff]), **acc_s)
        print(f"wrote {part}", flush=True)
    genome.close()


def stage_combine(a):
    sel, ref = select_rows(a.ref_tsv, a.min_spike)
    keys = [f"{r.anchor}|{r.partner}|{r.orient}|{r.d}" for r in sel.itertuples(index=False)]
    parts = sorted(f for f in os.listdir(a.part_dir) if f.endswith(".npz"))
    if not parts:
        raise SystemExit(f"error: no partials in {a.part_dir}")
    tot = None
    for f in parts:
        z = np.load(os.path.join(a.part_dir, f), allow_pickle=False)
        if list(z["rows"]) != keys:
            raise SystemExit(f"error: {f} rows differ from the selection")
        if tot is None:
            tot = {k: z[k].copy() for k in z.files}
        else:
            for k in z.files:
                if k.startswith(("S_", "SS_")) or k in ("n", "nchg", "anychg"):
                    tot[k] = tot[k] + z[k]
            tot["gate"] = np.maximum(tot["gate"], z["gate"])
    n = int(tot["n"][0])
    contexts = [str(c) for c in tot["contexts"]]
    print(f"{len(parts)} partials, {n} backgrounds, worst gate {tot['gate'][0]:.2e}",
          flush=True)
    tcrit = mp.stats.t.ppf(0.975, n - 1)

    frames = []
    for ci, ctx in enumerate(contexts):
        f = sel.copy()
        f["context"] = ctx
        for nm in ("real", "scr", "d1", "jit"):
            m, se, _ = mp._tstat(tot[f"S_{nm}"][ci], tot[f"SS_{nm}"][ci], n)
            f[f"{nm}_mean"], f[f"{nm}_se"] = m, se
        f["d1_lo"] = f.d1_mean - tcrit * f.d1_se
        f["d1_hi"] = f.d1_mean + tcrit * f.d1_se
        frames.append(f)
    res = pd.concat(frames, ignore_index=True).merge(
        ref, on=["anchor", "partner", "orient", "d", "context"], how="left")
    res["frac_junction_changed"] = np.tile(tot["anychg"] / n, len(contexts))
    res["mean_bases_changed"] = np.tile(tot["nchg"] / n, len(contexts))
    sg = np.sign(res.spike)
    res["arm1_pass"] = (sg * res.d1_lo > 0) & (sg * res.d1_hi > 0) & \
        (res.d1_mean.abs() >= 0.5 * res.spike.abs())
    res["jit_retention"] = (res.jit_mean - res.flank_med) / res.spike
    res["arm2_informative"] = res.frac_junction_changed >= 0.25
    res["arm2_pass"] = res.jit_retention >= 0.5
    res["verdict"] = np.select(
        [~res.arm1_pass, res.arm2_informative & ~res.arm2_pass,
         res.arm2_informative & res.arm2_pass],
        ["not_partner_specific", "junction_kmer", "syntax"], "undetermined")
    res["n_bg"] = n
    res.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {a.out} ({len(res)} rows)", flush=True)
    show = ["anchor", "partner", "orient", "d", "context", "spike", "real_mean", "scr_mean",
            "d1_lo", "d1_hi", "jit_retention", "frac_junction_changed", "verdict"]
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(res[show].round(4).to_string(index=False))
        print(res.groupby("context").verdict.value_counts().to_string())


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, choices=["run", "combine"])
    ap.add_argument("--ref-tsv", default=f"{A}/results/2026_1020_motif_pair_spacing.tsv")
    ap.add_argument("--ref-parts", default=f"{A}/results/parts_2026_1020_motif_pair_spacing")
    ap.add_argument("--bg-tsv", default=f"{A}/results/2026_1020_motif_pair_backgrounds.tsv")
    ap.add_argument("--stats-json", default=f"{A}/results/2026_1020_motif_pair_accstats.json")
    ap.add_argument("--part-dir", required=True)
    ap.add_argument("--out", default=None, help="combine: result table")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=500)
    ap.add_argument("--contexts", default="K562,HepG2")
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--min-spike", type=float, default=0.01)
    ap.add_argument("--expect-n", type=int, default=31, help="0 skips the check")
    ap.add_argument("--k", type=int, default=5, help="scrambles per background")
    ap.add_argument("--seed", type=int, default=0, help="6.1's background shuffle seed")
    ap.add_argument("--ctrl-seed", type=int, default=1)
    ap.add_argument("--gate-tol", type=float, default=1e-4)
    ap.add_argument("--gate-shift", type=int, default=0)
    ap.add_argument("--batch", type=int, default=512)
    return ap.parse_args()


def main():
    a = parse_args()
    if a.stage == "run":
        stage_run(a)
    else:
        if not a.out:
            raise SystemExit("error: --stage combine needs --out")
        stage_combine(a)


if __name__ == "__main__":
    main()
