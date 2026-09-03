#!/usr/bin/env python
"""Per-fold metrics with CIs and paired tests, driven by a 2.4-style config JSON.

Generalises 2.8 / 2.11, which were copies with hardcoded cell-type paths. Any comparison
expressible as a 2.4 config -- in-cell grids, cross-cell-type transfer, input-definition
swaps -- can now be scored per fold without another copy.

Emits, for every compare entry against the config's baseline:
  overall_pearson       r(observed, full prediction). THE metric when the question is
                        "how well would this predict H3K27ac in a cell type where I have
                        none" -- it scores the quantity you actually want to produce.
  overall_pearson_topq  the same on the top signal quintile (project reporting standard).
  residual_pearson      r(observed - atac_pred, model_pred - atac_pred), the mechanistic
                        readout of what the model adds beyond the baseline.
  incremental_r2        R2(model) - R2(baseline) against the observed signal.
  profile_pearson       mean over elements of the correlation between the observed base-
                        resolution profile and the predicted probabilities, via bpnetlite's
                        calculate_performance_measures with the same arguments the training
                        loop uses, so it is comparable to the Validation Profile Pearson
                        column in the training logs. Every metric above uses only the counts
                        head; this is the first thing here that scores the profile head at
                        all.
  profile_jsd           Jensen-Shannon distance between observed and predicted profiles.
                        LOWER IS BETTER, unlike every other column.

Residual-objective models are handled via `"residual": true` in the config entry: their
forward() emits the residual, so the baseline prediction is added back before scoring.

MODELS MAY HAVE DIFFERENT RECEPTIVE FIELDS. Each entry's input window is read from its own
saved model (`model.trimming`), not assumed. Windows are extracted once at the LARGEST
in-window across entries and cropped centrally per model. That is also what makes a
wide-vs-narrow comparison fair: extracting at the largest window means every model is
scored on the same regions -- the intersection of what each could accept -- so the
receptive field is not confounded with which regions near chromosome ends were scorable.

ENTRIES MAY HAVE DIFFERENT ACCESSIBILITY INPUTS. An entry may override
`accessibility_bw` (e.g. a comma-separated list of fragment-size channels); it defaults to
the baseline's. Distinct inputs are extracted separately and asserted to yield the same
valid-region mask, so rows stay aligned across entries.

Usage: 2.15.perfold_from_config.py CONFIG_JSON OUT_PREFIX ELEMENTS [--pair A B]...
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, ttest_rel, t as tdist

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks, normalize_accessibility
from bpnetlite.performance import calculate_performance_measures

GEN = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
FOLDS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/hg38_five_folds.json"
HW = 500
OUT_W = 2 * HW
TC = tdist.ppf(0.975, df=4)
dev = "cuda" if torch.cuda.is_available() else "cpu"

ap = argparse.ArgumentParser()
ap.add_argument("config"); ap.add_argument("out_prefix"); ap.add_argument("elements")
ap.add_argument("--pair", nargs=2, action="append", default=[])
ap.add_argument("--rc-average", action="store_true",
                help="Average each prediction with its reverse-complement. Off by "
                     "default so stored numbers stay comparable; turn on to measure "
                     "what test-time RC averaging is worth.")
a = ap.parse_args()
spec = json.load(open(a.config))
entries = [spec["baseline"]] + spec["compare"]


def model_path(cfg, fold):
    return f'{cfg["model_dir"]}/fold{fold}/multimodal_bpnet.torch'


# Resolve each entry's input window from its own fold0 checkpoint. Deriving it from
# model.trimming rather than a config key means the geometry can never disagree with the
# weights being scored.
for cfg in entries:
    m0 = torch.load(model_path(cfg, 0), map_location="cpu", weights_only=False)
    cfg["_in_window"] = OUT_W + 2 * m0.trimming
    cfg["_trimming"] = m0.trimming
    # Pre-2.15 configs (consumed by 2.8/2.11) omit this entirely and hardcoded the path in
    # the script, so say what is missing instead of raising KeyError from a setdefault.
    if "accessibility_bw" not in cfg:
        base_acc = spec["baseline"].get("accessibility_bw")
        if base_acc is None:
            raise SystemExit(
                f"error: entry {cfg['label']!r} has no 'accessibility_bw' and neither does "
                f"the baseline. This config predates 2.15, which needs an explicit "
                f"accessibility track per entry. Add one -- and check whether the models "
                f"were trained on atac.bw or atac_5p.bw, because the two are not "
                f"interchangeable (results/TARGET_PROVENANCE.md).")
        cfg["accessibility_bw"] = base_acc
    del m0
IN_W_MAX = max(c["_in_window"] for c in entries)
for cfg in entries:
    print(f'{cfg["label"]:<34} trimming={cfg["_trimming"]:>5} '
          f'in_window={cfg["_in_window"]:>5}')
print(f'extracting at in_window={IN_W_MAX} (regions valid for every entry)\n')
ACC_SPECS = sorted({c["accessibility_bw"] for c in entries})


def crop(arr, in_w):
    """Centre-crop an (N, C, IN_W_MAX) array to (N, C, in_w)."""
    if arr.shape[2] == in_w:
        return arr
    off = (arr.shape[2] - in_w) // 2
    return arr[:, :, off:off + in_w]


def predict(cfg, fold, seqs_max, accs_by_spec):
    md, mode = cfg["model_dir"], cfg["mode"]
    marker = f"{md}/fold{fold}/training_complete.json"
    if not os.path.exists(marker):
        raise SystemExit(f"error: {marker} missing; refusing to score an unfinished fold.")
    in_w = cfg["_in_window"]
    seqs = crop(seqs_max, in_w)
    accs_raw = crop(accs_by_spec[cfg["accessibility_bw"]], in_w)
    x = accs_raw
    if mode in ("multimodal", "atac"):
        st = json.load(open(f"{md}/fold{fold}/acc_normalization_stats.json"))
        x = normalize_accessibility(accs_raw, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, x], axis=1) if mode == "multimodal"
         else seqs if mode == "sequence" else x).astype(np.float32)
    m = torch.load(f"{md}/fold{fold}/multimodal_bpnet.torch", map_location="cpu",
                   weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = mode
    m = m.to(dev).eval()

    n_seq = 4 if mode in ("multimodal", "sequence") else 0

    def rc_input(xb):
        """Reverse-complement a batch, respecting the channel layout.

        One-hot channels are ACGT, so reversing the channel axis maps A<->T and C<->G;
        reversing the length axis completes the reverse complement. Accessibility channels
        are strand-agnostic coverage, so they are only reversed along length. This mirrors
        exactly what the training-time augmentation in ChIPSeqDataset does.
        """
        if n_seq:
            seq_part = torch.flip(xb[:, :n_seq], dims=[1, 2])
            if xb.shape[1] > n_seq:
                acc_part = torch.flip(xb[:, n_seq:], dims=[2])
                return torch.cat([seq_part, acc_part], dim=1)
            return seq_part
        return torch.flip(xb, dims=[2])

    def rc_profile(pr):
        """Undo the RC transform on a profile: reverse positions and swap strands."""
        return torch.flip(pr, dims=[1, 2])

    lcs, profs = [], []
    with torch.no_grad():
        for i in range(0, len(X), 256):
            xb = torch.from_numpy(X[i:i + 256]).to(dev)
            pr, lc = m(xb)
            if a.rc_average:
                pr_rc, lc_rc = m(rc_input(xb))
                lc = (lc + lc_rc) / 2
                # Average in probability space over the flattened strand+position
                # multinomial, which is the space the loss is defined in, then return to
                # log space. Log-probabilities are valid logits for the metric call.
                sh = pr.shape
                p1 = torch.softmax(pr.reshape(sh[0], -1), dim=-1)
                p2 = torch.softmax(rc_profile(pr_rc).reshape(sh[0], -1), dim=-1)
                pr = torch.log(((p1 + p2) / 2).clamp_min(1e-12)).reshape(sh)
            lcs.append(lc.squeeze(-1).cpu().numpy())
            profs.append(pr.cpu().numpy())
    m.to("cpu"); del X
    return np.concatenate(lcs), np.concatenate(profs)


def profile_metrics(logits, sigs, logcounts, top):
    """Profile metrics via bpnetlite, called exactly as the training loop calls it."""
    t = lambda z: torch.from_numpy(np.ascontiguousarray(z)).float()
    out = {}
    for suffix, mask in (("", slice(None)), ("_topq", top)):
        msr = calculate_performance_measures(
            t(logits[mask]), t(sigs[mask]), t(logcounts[mask]).reshape(-1, 1),
            kernel_sigma=7, kernel_width=81,
            measures=["profile_pearson", "profile_jsd"])
        out["profile_pearson" + suffix] = float(
            np.nan_to_num(np.asarray(msr["profile_pearson"], dtype=float)).mean())
        if not suffix:
            out["profile_jsd"] = float(
                np.nan_to_num(np.asarray(msr["profile_jsd"], dtype=float)).mean())
    return out


rows = []
folds_json = json.load(open(FOLDS))
for fold in range(5):
    els = load_peaks(a.elements, folds_json[str(fold)]["val"])
    b = spec["baseline"]
    accs_by_spec, seqs, sigs, ref_valid = {}, None, None, None
    for acc_spec in ACC_SPECS:
        sq, sg, ac, valid = extract_windows(
            els, GEN, b["signal_plus_bw"], b.get("signal_minus_bw"),
            acc_spec, IN_W_MAX, OUT_W, 0, is_peak=True)
        if ref_valid is None:
            seqs, sigs, ref_valid = sq, sg, valid
        else:
            # Rows must line up across accessibility inputs or every metric silently
            # compares different elements.
            assert np.array_equal(valid, ref_valid), (
                f"fold{fold}: accessibility input {acc_spec} kept a different region set "
                f"({valid.sum()} vs {ref_valid.sum()})")
            del sq, sg
        accs_by_spec[acc_spec] = ac
    obs = np.log1p(sigs.sum(axis=(1, 2)))
    base, base_prof = predict(b, fold, seqs, accs_by_spec)
    true_resid = obs - base
    top = obs >= np.quantile(obs, 0.8)
    r2b, r2bt = pearsonr(obs, base)[0] ** 2, pearsonr(obs[top], base[top])[0] ** 2
    for cfg in spec["compare"]:
        raw, prof = predict(cfg, fold, seqs, accs_by_spec)
        full = raw + base if cfg.get("residual") else raw
        mres = raw if cfg.get("residual") else raw - base
        row = {"fold": fold, "config": cfg["label"], "n": len(obs),
               "overall_pearson": pearsonr(obs, full)[0],
               "overall_pearson_topq": pearsonr(obs[top], full[top])[0],
               "residual_pearson": pearsonr(true_resid, mres)[0],
               "incremental_r2": pearsonr(obs, full)[0] ** 2 - r2b,
               "incremental_r2_topq": pearsonr(obs[top], full[top])[0] ** 2 - r2bt}
        row.update(profile_metrics(prof, sigs, full, top))
        rows.append(row); del prof
    row = {"fold": fold, "config": b["label"], "n": len(obs),
           "overall_pearson": pearsonr(obs, base)[0],
           "overall_pearson_topq": pearsonr(obs[top], base[top])[0],
           "residual_pearson": np.nan, "incremental_r2": 0.0,
           "incremental_r2_topq": 0.0}
    row.update(profile_metrics(base_prof, sigs, base, top))
    rows.append(row); del base_prof, sigs
    del seqs, accs_by_spec
    print(f"fold{fold}: n={len(obs):,}", flush=True)

df = pd.DataFrame(rows)
p1 = f"{P}/results/{a.out_prefix}per_fold.tsv"
df.round(4).to_csv(p1, sep="\t", index=False)
print("\nWrote", p1)

METRICS = ["overall_pearson", "overall_pearson_topq", "residual_pearson",
           "incremental_r2", "incremental_r2_topq",
           "profile_pearson", "profile_pearson_topq", "profile_jsd"]
srows = []
for label, g in df.groupby("config", sort=False):
    r = {"config": label}
    for m in METRICS:
        v = g[m].to_numpy(float)
        if np.isnan(v).all():
            r[m] = "-"
            continue
        mu = v.mean(); half = TC * v.std(ddof=1) / np.sqrt(len(v))
        r[m] = f"{mu:.3f} [{mu-half:.3f}, {mu+half:.3f}]"
    srows.append(r)
summ = pd.DataFrame(srows)
p2 = f"{P}/results/{a.out_prefix}fold_summary.tsv"
summ.to_csv(p2, sep="\t", index=False)
print("Wrote", p2)
print(summ.to_string(index=False))

if a.pair:
    print("\n--- paired differences (A - B), within fold ---")
    piv = {m: df.pivot(index="fold", columns="config", values=m) for m in METRICS}
    for A, B in a.pair:
        for m in ["overall_pearson", "overall_pearson_topq",
                  "profile_pearson", "profile_pearson_topq"]:
            d = (piv[m][A] - piv[m][B]).to_numpy()
            mu = d.mean(); half = TC * d.std(ddof=1) / np.sqrt(len(d))
            p = ttest_rel(piv[m][A], piv[m][B]).pvalue
            print(f"  {m:<22} {A} - {B}: {mu:+.4f} [{mu-half:+.4f}, {mu+half:+.4f}]  p={p:.4f}")
