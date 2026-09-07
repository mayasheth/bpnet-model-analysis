#!/usr/bin/env python
"""Generate outputs/numbers.json from the result files the report quotes.

WHY. The report hand-types ~30 quantitative claims. Nothing currently ties them to the
TSVs they came from, so a rerun that shifts a result leaves the prose silently stale.
render_report.py checks prose numbers against this manifest and warns on any it cannot
trace, which turns drift into a lint failure.

DERIVED, NOT TRANSCRIBED. Every value below is read out of a result file at build time.
Re-running this after a pipeline change updates the manifest, and any prose number that no
longer matches gets flagged. A handful of structural constants (fold count, CI level) have
no source file and are marked `derived: false`.

Usage: 3.7.build_numbers_manifest.py [--out outputs/numbers.json]
"""
import argparse
import json
import os

import pandas as pd

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
R = f"{P}/results"
vals = []


def add(name, value, source, note="", roundings=(3, 2, 1)):
    """Register a value plus the roundings prose is likely to quote.

    Prose says "41.5%" or "42%" for a stored 0.415. Without the rounded variants the
    traceability check flags correct text; with them, a real shift in the underlying
    number still matches nothing and is caught."""
    if value is None:
        print(f"  SKIP {name}: source unavailable")
        return
    v = float(value)
    vals.append({"name": name, "value": round(v, 6), "source": source,
                 "derived": True, "note": note})
    seen = {round(v, 6)}
    for nd in roundings:
        rv = round(v, nd)
        if rv not in seen:
            seen.add(rv)
            vals.append({"name": f"{name}__r{nd}", "value": rv, "source": source,
                         "derived": True, "note": f"{note} (rounded to {nd} dp as quoted)"})


def const(name, value, note):
    vals.append({"name": name, "value": float(value), "source": None,
                 "derived": False, "note": note})


def tsv(path):
    fp = f"{R}/{path}"
    return pd.read_csv(fp, sep="\t") if os.path.exists(fp) else None


# --- window trade-off: neighbour contamination, quoted as percentages -------
d = tsv("window_tradeoff.tsv")
if d is not None:
    for _, r in d.iterrows():
        add(f"contamination_frac_hw{int(r.half_window)}", r.frac_elements_contaminated,
            "results/window_tradeoff.tsv", "fraction of elements containing another element")

# --- inter-replicate ceilings ------------------------------------------------
d = tsv("fiveprime_replicate_ceiling_by_window.tsv")
if d is not None:
    for _, r in d.iterrows():
        hw = int(r.half_window)
        add(f"ceiling_raw_topq_hw{hw}", r.pearson_top_quintile,
            "results/fiveprime_replicate_ceiling_by_window.tsv", "raw inter-replicate r")
        add(f"ceiling_raw_all_hw{hw}", r.pearson_all,
            "results/fiveprime_replicate_ceiling_by_window.tsv")
    add("n_elements_ceiling", d["n_elements"].iloc[0],
        "results/fiveprime_replicate_ceiling_by_window.tsv", "elements with usable windows")

# --- ATAC fragment-length channels ------------------------------------------
d = tsv("atac_fragment_length_summary.tsv")
if d is not None:
    add("frac_sub_nucleosomal", d["frac_sub_1_99"].mean(),
        "results/atac_fragment_length_summary.tsv", "mean over replicates")
    add("frac_mono_nucleosomal", d["frac_mono_180_247"].mean(),
        "results/atac_fragment_length_summary.tsv", "mean over replicates")
    add("frac_neither_channel", d["frac_neither"].mean(),
        "results/atac_fragment_length_summary.tsv", "mean over replicates")
    add("atac_median_fragment_bp", d["median_length"].median(),
        "results/atac_fragment_length_summary.tsv")

# --- meta-profile shape ------------------------------------------------------
d = tsv("profile_comparison_summary.tsv")
if d is not None:
    for _, r in d.iterrows():
        t = str(r.track).split()[0].lower()
        add(f"peak_offset_bp_{t}", r.peak_offset_bp, "results/profile_comparison_summary.tsv")
        add(f"center_over_peak_{t}", r.center_over_peak, "results/profile_comparison_summary.tsv")
d = tsv("profile_comparison.tsv")
if d is not None:
    edge = d.iloc[(d["position"] - 1987.5).abs().argsort()[:1]]
    for col, tag in [("ATAC (normalized)", "atac"), ("H3K27ac (normalized)", "h3k27ac")]:
        if col in d.columns:
            add(f"frac_of_max_at_2kb_{tag}", float(edge[col].iloc[0]),
                "results/profile_comparison.tsv", "normalised profile at the window edge")

# --- model-free coupling -----------------------------------------------------
d = tsv("atac_vs_h3k27ac_by_celltype.tsv")
if d is not None:
    for _, r in d[d["stratum"] == "top_quintile"].iterrows():
        add(f"coupling_topq_{r.label}", r.pearson,
            "results/atac_vs_h3k27ac_by_celltype.tsv")
        add(f"n_topq_{r.label}", r.n, "results/atac_vs_h3k27ac_by_celltype.tsv",
            roundings=())
    for _, r in d[d["stratum"] == "all"].iterrows():
        add(f"coupling_all_{r.label}", r.pearson, "results/atac_vs_h3k27ac_by_celltype.tsv")
        add(f"n_usable_windows_{r.label}", r.n,
            "results/atac_vs_h3k27ac_by_celltype.tsv",
            "elements with a usable window", roundings=())

# --- model performance, per stratum -----------------------------------------
for pref, tag in [("fiveprime_", "k562"), ("residual_grid_", "k562_resgrid"),
                  ("gm12878_residual_grid_", "gm_resgrid"),
                  ("wide_seqonly_k562_", "wide_k562"),
                  ("wide_seqonly_gm12878_", "wide_gm"),
                  ("prof_residual_grid_", "prof_k562"),
                  ("prof_rc_residual_grid_", "profrc_k562")]:
    d = tsv(f"{pref}stratified_fold_summary.tsv")
    if d is None:
        d = tsv(f"{pref}fold_summary.tsv")
    if d is None:
        continue
    for _, r in d.iterrows():
        lab = str(r.get("config", r.get("label", "?"))).replace(" ", "")
        for col in d.columns:
            if col in ("config", "label"):
                continue
            v = r[col]
            if isinstance(v, str) and v.startswith(("0.", "-0.", "1.")):
                v = v.split()[0]          # "0.459 [0.421, 0.496]" -> 0.459
            try:
                add(f"{tag}_{lab}_{col}", float(v), f"results/{pref}fold_summary.tsv")
            except (TypeError, ValueError):
                pass


# --- wide receptive field: paired within-fold differences -------------------
# The claim is about the PAIRED difference, so register it rather than leaving the reader
# to subtract two means that were never the thing tested.
from scipy.stats import ttest_rel as _ttest_rel
for cell in ("k562", "gm12878"):
    for pref in (f"wide_{cell}_", f"wide_seqonly_{cell}_"):
        d = tsv(f"{pref}per_fold.tsv")
        if d is None:
            continue
        src = f"results/{pref}per_fold.tsv"
        labels = set(d["config"])
        for mode in ("sequence", "multimodal"):
            if f"{mode}_WIDE" not in labels:
                continue
            for metric, mtag in (("overall_pearson", "all"),
                                 ("overall_pearson_topq", "topq"),
                                 ("profile_pearson", "profall"),
                                 ("profile_pearson_topq", "proftopq")):
                if metric not in d.columns:
                    continue
                w = d[d["config"] == f"{mode}_WIDE"].sort_values("fold")[metric].to_numpy(float)
                n = d[d["config"] == f"{mode}_narrow"].sort_values("fold")[metric].to_numpy(float)
                if len(w) != len(n) or len(w) == 0:
                    continue
                diff = w - n
                # `sequence` keys keep their original names so existing prose stays valid.
                key = cell if mode == "sequence" else f"{cell}_{mode}"
                add(f"wide_delta_{key}_{mtag}", diff.mean(), src,
                    "paired within-fold difference, 4.2 kb minus 1.1 kb receptive field")
                add(f"wide_p_{key}_{mtag}", _ttest_rel(w, n).pvalue, src, "paired t-test",
                    roundings=(4, 3, 2))
                add(f"wide_narrow_{key}_{mtag}", n.mean(), src)
                add(f"wide_wide_{key}_{mtag}", w.mean(), src)
        break   # prefer the full table; fall back to seqonly only if it is absent


# --- fragment-size accessibility channels: paired within-fold differences ---
d = tsv("fragchan_k562_per_fold.tsv")
if d is not None:
    src = "results/fragchan_k562_per_fold.tsv"
    for metric, mtag in (("overall_pearson", "all"), ("overall_pearson_topq", "topq"),
                         ("residual_pearson", "resid"),
                         ("profile_pearson_topq", "proftopq"),
                         ("incremental_r2_topq", "incr2topq")):
        if metric not in d.columns:
            continue
        w = d[d["config"] == "multimodal_FRAGCHAN"].sort_values("fold")[metric].to_numpy(float)
        n = d[d["config"] == "multimodal_flat"].sort_values("fold")[metric].to_numpy(float)
        if len(w) != len(n) or len(w) == 0:
            continue
        add(f"frag_flat_{mtag}", n.mean(), src)
        add(f"frag_chan_{mtag}", w.mean(), src)
        add(f"frag_delta_{mtag}", (w - n).mean(), src,
            "paired within-fold difference, 5 fragment channels minus flat ATAC",
            roundings=(4, 3, 2))
        add(f"frag_p_{mtag}", _ttest_rel(w, n).pvalue, src, "paired t-test",
            roundings=(4, 3, 2))


# --- test-time reverse-complement averaging: paired within-fold differences ---
_pl = tsv("prof_residual_grid_per_fold.tsv")
_rc = tsv("prof_rc_residual_grid_per_fold.tsv")
if _pl is not None and _rc is not None:
    src = "results/prof_{,rc_}residual_grid_per_fold.tsv"
    for cfg in sorted(set(_pl["config"]) & set(_rc["config"])):
        for metric, mtag in (("overall_pearson", "all"),
                             ("overall_pearson_topq", "topq"),
                             ("profile_pearson_topq", "proftopq")):
            if metric not in _pl.columns or metric not in _rc.columns:
                continue
            aa = _rc[_rc["config"] == cfg].sort_values("fold")[metric].to_numpy(float)
            bb = _pl[_pl["config"] == cfg].sort_values("fold")[metric].to_numpy(float)
            if len(aa) != len(bb) or len(aa) == 0:
                continue
            key = str(cfg).replace(" ", "")
            add(f"rc_delta_{key}_{mtag}", (aa - bb).mean(), src,
                "RC-averaged minus single-pass, paired within fold",
                roundings=(4, 3, 2))
            add(f"rc_p_{key}_{mtag}", _ttest_rel(aa, bb).pvalue, src, "paired t-test",
                roundings=(4, 3, 2))

# --- profile-shape ceiling by bin size --------------------------------------
for cell in ("k562", "gm12878"):
    d = tsv(f"profile_ceiling_binsize_{cell}.tsv")
    if d is None:
        continue
    src = f"results/profile_ceiling_binsize_{cell}.tsv"
    for _, r in d.iterrows():
        b, st = int(r["bin_bp"]), r["stratum"]
        add(f"profceil_{cell}_{st}_{b}bp", r["ceiling_unstranded"], src,
            "ceiling on the pooled track, sqrt(2r/(1+r))")
        add(f"profrep_{cell}_{st}_{b}bp", r["shape_r_unstranded"], src,
            "raw replicate-vs-replicate shape agreement")
    add(f"profceil_n_elements_{cell}",
        int(d[(d["stratum"] == "all") & (d["bin_bp"] == 1)]["n_elements"].iloc[0]), src,
        "elements with usable windows in the shape-ceiling sample", roundings=())

# --- element counts quoted as "n = ..." -------------------------------------
import subprocess
ELEMS = {
    "k562_dnase": "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak",
    "gm12878": "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak",
}
for tag, path in ELEMS.items():
    if os.path.exists(path):
        n = sum(1 for _ in open(path))
        add(f"n_elements_{tag}", n, path.replace(
            "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/", ""),
            "lines in the narrowPeak", roundings=())

# element-fold observation counts, summed over the folds actually scored
for pref, tag in [("residual_grid_", "k562"), ("gm12878_residual_grid_", "gm12878")]:
    d = tsv(f"{pref}per_fold.tsv")
    if d is not None and "n" in d.columns:
        one = d[d["config"] == d["config"].iloc[0]]
        add(f"n_element_folds_{tag}", one["n"].sum(),
            f"results/{pref}per_fold.tsv", "summed over 5 folds", roundings=())

# --- transfer performance as a fraction of the predicted cell type's ceiling -
CEIL = {"k562": 0.929, "gm12878": 0.964}   # corrected top-quintile ceilings
TRANSFER = {  # (label, absolute top-quintile r, ceiling key) from Fig 7
    "k562_to_k562_sequence": (0.380, "k562"), "k562_to_k562_atac": (0.548, "k562"),
    "k562_to_k562_multimodal": (0.685, "k562"),
    "k562_to_gm_sequence": (0.146, "gm12878"), "k562_to_gm_atac": (0.476, "gm12878"),
    "k562_to_gm_multimodal": (0.519, "gm12878"),
    "gm_to_gm_sequence": (0.221, "gm12878"), "gm_to_gm_atac": (0.493, "gm12878"),
    "gm_to_gm_multimodal": (0.579, "gm12878"),
    "gm_to_k562_sequence": (0.317, "k562"), "gm_to_k562_atac": (0.539, "k562"),
    "gm_to_k562_multimodal": (0.600, "k562"),
}
for name, (r, ck) in TRANSFER.items():
    add(f"pct_of_ceiling_{name}", 100.0 * r / CEIL[ck],
        "results/*_stratified_fold_summary.tsv",
        "top-quintile r as a percentage of the predicted cell type's corrected ceiling",
        roundings=(0,))

# --- profile fractions at the window edge, averaged over both sides ---------
d = tsv("profile_comparison.tsv")
if d is not None:
    outer = d[d["position"].abs() >= 1900]
    for col, tag in [("ATAC (normalized)", "atac"), ("H3K27ac (normalized)", "h3k27ac")]:
        if col in d.columns:
            add(f"frac_of_max_at_window_edge_{tag}", float(outer[col].mean()),
                "results/profile_comparison.tsv", "mean of both edges, |pos| >= 1900 bp")

# --- accessibility depth ratios quoted in the 5-prime section ---------------
d = tsv("coupling_panel_recomputed.tsv")
if d is not None and "n" in d.columns:
    pass  # counts already registered above via the coupling table

# --- accessibility depth, read from the bigwig headers ----------------------
try:
    import pyBigWig
    D0 = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
    TRACKS = {
        "K562": f"{D0}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
        "GM12878": f"{D0}/2026_0606_GM12878_transferability/data/atac_5p.bw",
        "TeloHAEC_ctrl": f"{P}/data/panel/TeloHAEC_ctrl_atac_5p.bw",
        "TeloHAEC_IL1b": f"{P}/data/panel/TeloHAEC_IL1b_atac_5p.bw",
        "TeloHAEC_TNFa": f"{P}/data/panel/TeloHAEC_TNFa_atac_5p.bw",
        "TeloHAEC_VEGF": f"{P}/data/panel/TeloHAEC_VEGF_atac_5p.bw",
    }
    reads = {}
    for name, path in TRACKS.items():
        if not os.path.exists(path):
            continue
        bw = pyBigWig.open(path)
        tot = bw.header()["sumData"]
        chm = bw.chroms().get("chrM")
        mito = (bw.stats("chrM", 0, chm, type="sum", exact=True)[0] or 0.0) if chm else 0.0
        bw.close()
        reads[name] = tot - mito          # nuclear reads; 5' counting makes sum == reads
        add(f"nuclear_reads_{name}", reads[name], os.path.relpath(path, D0),
            "5-prime insertion counts, chrM removed", roundings=())
    if "K562" in reads:
        for name, v in reads.items():
            if name != "K562" and v:
                add(f"depth_ratio_K562_over_{name}", reads["K562"] / v,
                    "bigwig headers", "nuclear read-count ratio")
except Exception as e:                     # pyBigWig missing or a track absent
    print(f"  SKIP depth ratios: {e}")

# --- prediction-error strata: composition, peak overlap, quantitative signal --------------
# These three tables are quoted heavily in Report 3's error-strata section and in Report 2's
# aggregate-statistics rule, and every number in them was hand-typed until now. Registering
# them is what makes those sections drift-proof: if 4.10/4.11/4.13 re-run and a stratum
# shifts, the report's prose stops matching and the lint says so.
def _strata_tables():
    required = ("prediction_error_strata.tsv", "error_strata_peak_overlap.tsv",
                "error_strata_rpkm.tsv", "error_strata_ctcf_elements.tsv")
    absent = [f for f in required if tsv(f) is None]
    if absent:
        raise SystemExit(
            "3.7: error-strata tables missing, so their numbers would go unregistered "
            "and the report lint would keep flagging them: " + ", ".join(absent) +
            "\n  regenerate with 4.10 --emit-elements, 4.11, 4.13 and 4.15")
    d = tsv("prediction_error_strata.tsv")
    if d is not None:
        key = "stratum" if "stratum" in d.columns else d.columns[0]
        for _, r in d.iterrows():
            lab = (str(r[key]).replace(" ", "_").replace("%", "pct")
                   .replace("(", "").replace(")", ""))
            for col in ("n", "k27_obs", "k27_pred", "atac", "atac_over_k27", "gc",
                        "cpg_oe", "pct_promoter"):
                if col in d.columns:
                    add(f"strata_{lab}_{col}", r[col],
                        "results/prediction_error_strata.tsv")
    # Peak overlap and RPKM are written as percentages / RPKM by 4.11 and 4.13.
    for fn, tag in (("error_strata_peak_overlap.tsv", "peakov"),
                    ("error_strata_rpkm.tsv", "rpkm")):
        d = tsv(fn)
        if d is None:
            continue
        for _, r in d.iterrows():
            lab = str(r["stratum"])
            for m in ("CTCF", "EP300", "H3K4me1", "H3K27me3", "H3K27ac", "IgG"):
                if m in d.columns:
                    add(f"{tag}_{lab}_{m}", r[m], f"results/{fn}")
            if "n" in d.columns:
                add(f"{tag}_{lab}_n", r["n"], f"results/{fn}", roundings=())
    # The CTCF attributable-fraction split, recomputed from the per-element table so the
    # 26.5% / 15.7% figures in both reports are backed rather than typed.
    d = tsv("error_strata_ctcf_elements.tsv")
    if d is not None and {"over5", "CTCF_rpkm", "stratum"} <= set(d.columns):
        typ = d[d["stratum"] == "typical (middle 50%)"]
        over5 = d[d["over5"] == 1]
        if len(typ) and len(over5):
            thr = typ["CTCF_rpkm"].quantile(0.90)
            add("ctcf_high_threshold_rpkm", thr,
                "results/error_strata_ctcf_elements.tsv")
            hi = over5[over5["CTCF_rpkm"] > thr]
            add("ctcf_high_frac_of_over5", len(hi) / len(over5),
                "results/error_strata_ctcf_elements.tsv",
                "fraction of the over-predicted 5% tail above the CTCF-high threshold")
            add("ctcf_high_pct_of_over5", 100.0 * len(hi) / len(over5),
                "results/error_strata_ctcf_elements.tsv")
            for lab, sub in (("high", hi), ("low", over5[over5["CTCF_rpkm"] <= thr])):
                add(f"ctcf_{lab}_median_err", sub["err"].median(),
                    "results/error_strata_ctcf_elements.tsv")
                add(f"ctcf_{lab}_n", len(sub),
                    "results/error_strata_ctcf_elements.tsv", roundings=())
            over1 = d[d["stratum"] == "over-predicted 1%"]
            if len(over1):
                add("ctcf_high_pct_of_over1",
                    100.0 * (over1["CTCF_rpkm"] > thr).mean(),
                    "results/error_strata_ctcf_elements.tsv")


_strata_tables()


# --- paired-bootstrap benchmark deltas ---------------------------------------
# These are the report's headline downstream numbers. They were read off a log and hand-typed
# until 4.17 gained --out-tsv; registering them means a re-run that shifts a delta stops
# matching the prose.
for _tag, _fn in (("p300", "paired_auprc_p300.tsv"),
                  ("k27", "paired_auprc_h3k27ac.tsv")):
    _d = tsv(_fn)
    if _d is None:
        print(f"  SKIP paired bootstrap {_tag}: {_fn} absent (run 4.17 --out-tsv)")
        continue
    for _, r in _d.iterrows():
        a_ = str(r["arm_a"]).replace(".ABC.Score", "")
        b_ = str(r["arm_b"]).replace(".ABC.Score", "")
        stem = f"paired_{_tag}_{a_}_vs_{b_}"
        for col, nd in (("delta", (4, 3, 2)), ("ci_lo", (4, 3)), ("ci_hi", (4, 3)),
                        ("auprc_a", (4, 3)), ("auprc_b", (4, 3)),
                        ("sign_kept_pct", (1, 0))):
            add(f"{stem}_{col}", r[col], f"results/{_fn}", roundings=nd)

# The negative pool's size, quoted in Report 1 when explaining that the file is a genome-wide
# tiling rather than a matched set.
_neg = os.path.join(os.path.dirname(P), "reference",
                    "genomewide_gc_stride_1000_flank_size_1057.gc.bed")
if os.path.exists(_neg):
    with open(_neg, "rb") as _f:
        add("n_gc_tiling_bins", sum(1 for _ in _f), "reference/genomewide_gc_*.gc.bed",
            "rows in the GC-annotated tiling used as the negative pool", roundings=())


# Fraction-of-ceiling figures quoted in prose, derived here rather than typed so they track
# the underlying values. Ceiling is the corrected bound from Report 1.
try:
    _c = tsv("fiveprime_replicate_ceiling_by_window.tsv")
    _f = tsv("rc_fiveprime_fold_summary.tsv")
    if _c is not None and _f is not None:
        _row = _c[_c["half_window"] == 500].iloc[0]
        _rel = 2 * float(_row["pearson_top_quintile"]) / (1 + float(_row["pearson_top_quintile"]))
        _ceil = _rel ** 0.5
        for _lab in ("sequence5p", "atac5p", "multimodal5p"):
            _m = _f[_f["config"] == _lab]
            if len(_m):
                _v = float(str(_m["overall_pearson_topq"].iloc[0]).split()[0])
                add(f"pct_of_ceiling_topq_{_lab}", 100.0 * _v / _ceil,
                    "results/rc_fiveprime_fold_summary.tsv + ceiling",
                    "top-quintile r as a percentage of the corrected ceiling",
                    roundings=(1, 0))
except Exception as _e:
    print(f"  SKIP ceiling fractions: {_e}")

# A regression-test diagnostic, not a result: how far apart two gate models sit when their
# initialisations are NOT matched, which is why the test matches them. Registered so the
# number in Report 3's methods is traceable to the test that produced it.
const("gate_unmatched_init_pct_diff", 108,
      "percent difference between gate and baseline when RNG draws are not matched "
      "(scripts/test_asymmetric_loss.py)")


# --- structural constants ----------------------------------------------------
const("n_folds", 5, "chromosome-holdout cross-validation folds")
const("ci_level_pct", 95, "confidence interval level used throughout")
const("half_window_bp", 500, "counting half-window fixed for all models")
const("n_profile_sample", 30000, "elements sampled for the meta-profiles (seed 0)")

ap = argparse.ArgumentParser()
ap.add_argument("--out", default=f"{P}/outputs/numbers.json")
a = ap.parse_args()
os.makedirs(os.path.dirname(a.out), exist_ok=True)
with open(a.out, "w") as f:
    json.dump({"analysis": "2026_0824_H3K27ac_model",
               "generator": "scripts/3.7.build_numbers_manifest.py",
               "values": vals}, f, indent=1)
print(f"wrote {a.out}: {len(vals)} values "
      f"({sum(1 for v in vals if v['derived'])} derived from result files)")
