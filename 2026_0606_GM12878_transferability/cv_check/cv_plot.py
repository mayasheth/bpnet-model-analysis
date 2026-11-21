"""Draft bar charts for Fig. C.2 and Fig. 4.4c,d using held-out (CV) predictions.

Bars: CV Pearson r. Open circles: current mean-of-folds value, for comparison.
Bar order: descending CV r on all elements, shared between panels (as in the current figures).
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import pandas as pd

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/cv_check"
res = pd.read_csv(f"{D}/cv_vs_mean_accuracy.tsv", sep="\t")

REPL = {"K562": {"all": 0.8762, "p300+": 0.7458, "n": (150519, 27176)},
        "GM12878": {"all": 0.8812, "p300+": 0.8345, "n": (154224, 20396)}}

COLORS = {
    ("seq+atac", True): "#4f2d5e", ("seq+atac", False): "#7a5d82",
    ("seq", True): "#5b7ec2", ("seq", False): "#869fd0",
    ("atac", True): "#b8673f", ("atac", False): "#cc8b69",
}


def kind(label):
    if "sequence + ATAC" in label:
        return "seq+atac"
    return "atac" if "ATAC only" in label else "seq"


def draw(fig_id, cell, title_suffix, out):
    sub = res[res.figure == fig_id]
    allv = sub[sub.subset == "all"].set_index("model")
    pk = sub[sub.subset == "p300+"].set_index("model")
    order = allv.sort_values("r_cv", ascending=False).index.tolist()
    labels = [f"{cell} ChIP-seq replicates"] + order
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 3.6), gridspec_kw={"wspace": 0.08})
    for ax, tab, subset, title, n in [
        (axes[0], allv, "all", "All elements", REPL[cell]["n"][0]),
        (axes[1], pk, "p300+", "p300$^{+}$ elements", REPL[cell]["n"][1]),
    ]:
        vals = [REPL[cell][subset]] + [tab.loc[m, "r_cv"] for m in order]
        means = [None] + [tab.loc[m, "r_mean"] for m in order]
        cols = ["#c5ccce"] + [COLORS[(kind(m), m.startswith(cell))] for m in order]
        y = list(range(len(labels)))[::-1]
        ax.barh(y, vals, color=cols, height=0.72)
        for yi, v, mv in zip(y, vals, means):
            txt = f"{v:.3f}" if mv is None else f"{v:.3f} ({mv:.3f})"
            ax.text(max(v, mv or 0) + 0.025, yi, txt, va="center", fontsize=8, color="black")
            if mv is not None:
                ax.plot(mv, yi, marker="o", mfc="none", mec="black", ms=5, mew=0.9)
        ax.set_xlim(0, 1.1)
        ax.set_title(title, fontsize=11, color="black")
        ax.set_xlabel(f"Pearson's $r$ (log p300 counts, {cell})", color="black")
        ax.text(1.08, -0.6, f"$n$ = {n:,}", ha="right", va="bottom", fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors="black")
        ax.set_yticks(y)
        ax.set_yticklabels(labels if ax is axes[0] else [], fontsize=9, color="black")
    fig.suptitle(f"{title_suffix}: bars = held-out (CV) predictions; circles and parentheses = current mean-of-folds",
                 fontsize=9, y=1.02)
    fig.savefig(f"{D}/{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{D}/{out}.png", bbox_inches="tight", dpi=200)


draw("C.2", "K562", "Fig. C.2 draft", "figC2_cv")
draw("4.4c,d", "GM12878", "Fig. 4.4c,d draft", "fig4_4cd_cv")
print("plots written")
