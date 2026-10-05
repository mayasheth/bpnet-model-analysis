#!/usr/bin/env python3
"""Move every TF-motif section of report 6 into report 7 (2026-10-04, Maya: all motif work in one
report). Sections move by title, unchanged except for figure numbers, so legends and numbers
cannot drift. Figure references are remapped through placeholders; a reference to a figure that
now lives in the other report becomes "report N Fig. X". Every reference must map, or the script
stops.

Writes:  report6_improving_transfer.qmd   (motif sections removed, title / summary / goals updated)
         report7_partA.inc.md             (the moved sections + the TF terms, for build_report7.py)
Usage:   build_report_tf_split.py [--suffix .new]   (suffix: write beside the originals for review)
"""
import argparse
import re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--suffix", default="")
a = ap.parse_args()
R6 = Path("report6_improving_transfer.qmd")
text = R6.read_text()

MOVE = [
    "The transferred model under-ranks GATA elements, and fixing only those would recover half the gap",
    "The simplest TF-list estimate from ATAC keeps the cell-specific signal best",
    "Given the true TF list, the model removes a third of the GATA error; the cross-cell 'pair' loss removes two thirds",
    "Three training runs per model put the noise at 0.002 AUPRC in ABC and 0.005 in ENCODE-rE2G",
    "Among 19 ENCODE samples with p300, only stomach shows K562-like GATA activity",
]
R6_KEEP = {1: 1, 2: 2, 3: 3, 9: 4, 10: 5, 11: 6, 12: 7}   # old report 6 fig -> new report 6 fig
R6_TO_R7 = {4: 2, 5: 3, 6: 4, 7: 5, 8: 6}                # old report 6 fig -> report 7 fig
assert set(R6_KEEP) | set(R6_TO_R7) == set(range(1, 13)) and not set(R6_KEEP) & set(R6_TO_R7)

heads = list(re.finditer(r"^## (.+?)\s*$", text, re.M))
spans = {}
for i, h in enumerate(heads):
    end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
    spans[h.group(1).strip()] = (h.start(), end)
for t in MOVE:
    assert t in spans, f"section not found: {t!r}"

REF = re.compile(r"\b(Figs?\.?|Figures?)(\s+)(\d+)([a-z]?)(?:(\s*-\s*)(\d+)([a-z]?))?")


def remap(s, local, cross, cross_name):
    def one(n):
        n = int(n)
        if n in local:
            return None, local[n]
        if n in cross:
            return cross_name, cross[n]
        raise SystemExit(f"unmapped figure reference {n} in: {s[:80]!r}")

    def sub(m):
        word, sp, n1, l1, dash, n2, l2 = m.groups()
        r1 = one(n1)
        if n2 is None:
            where, k = r1
            return f"{where + ' ' if where else ''}{word}{sp}@@{k}@@{l1}"
        r2 = one(n2)
        assert r1[0] == r2[0], f"range spans reports: {m.group(0)}"
        where = r1[0]
        return f"{where + ' ' if where else ''}{word}{sp}@@{r1[1]}@@{l1}{dash}@@{r2[1]}@@{l2}"

    out = REF.sub(sub, s)
    return re.sub(r"@@(\d+)@@", r"\1", out)


moved = "\n".join(remap(text[slice(*spans[t])].rstrip("\n") + "\n", R6_TO_R7, R6_KEEP, "report 6") for t in MOVE)
terms = [l for l in text.splitlines() if re.match(r"- \*(Partial correlation given ATAC|TF-activity track|True TF list|ATAC TF list|Seed replicate)\*", l)]
assert len(terms) == 5, terms
Path("report7_partA.inc.md" + a.suffix).write_text(
    "TF-activity terms (from report 6):\n\n" + "\n".join(remap(l, R6_TO_R7, R6_KEEP, "report 6") for l in terms)
    + "\n\n@@TERMS_END@@\n\n" + moved)

# ---------------------------------------------------------------- new report 6
keep = text
for t in sorted(MOVE, key=lambda t: -spans[t][0]):
    s, e = spans[t]
    keep = keep[:s] + keep[e:]
OLD_TITLE = re.search(r'^title: ".*"$', keep, re.M).group(0)
keep = keep.replace(OLD_TITLE, 'title: "Seven changes to the transferred p300 model raise correlation with p300 by up to 0.12 '
                               'but move ABC by at most 0.011 AUPRC, against a 0.042 gap to the in-cell model"', 1)
s_sum, e_sum = re.search(r"^## Summary\n", keep, re.M).end(), re.search(r"^## Goals\n", keep, re.M).start()
NEW_SUMMARY = """
- **No change closes the transfer gap.** Every variant moves ABC by at most 0.011 AUPRC against a 0.042 gap between the transferred and the in-cell model (Fig. 1). The variants that give the model TF information (the TF-activity track and its training schemes), the GATA error analysis behind them, the TF-list estimators and the panel's TF coverage are in report 7, which collects all the TF-motif work; the best of them, the ATAC-derived TF list, gains +0.006 in ABC averaged over three runs (report 7, Fig. 5).
- **Single training runs cannot separate these variants.** The same model trained three times varies with a standard deviation of 0.002 AUPRC in ABC and 0.005 in ENCODE-rE2G (report 7, Fig. 5).
- **Correlation gains do not become AUPRC gains, and the reverse also happens.** The 'pair' loss raises correlation beyond ATAC by +0.125 and neighbour re-ranking by +0.10 to +0.12 over 153,545 elements, with no resolvable ABC change; the ranking loss at weight 10 lowers it by 0.015 and gives the largest single-run ABC gain (+0.0087) (Fig. 1a, Fig. 3). Only the oracle re-rankings move both. Only the benchmark decides.
- **Fragment length adds no deployable signal.** Fragment-size input channels lose to the flat input (Fig. 7a). Fragment features predict H3K27ac at the element level, but through where accessibility lies rather than nucleosome length (Figs 5-6), and correcting transferred H3K27ac with them does not reach the benchmark (Fig. 7c).
- **In ABC, the ATAC term costs observed marks AUPRC.** Observed p300 alone scores 0.591 against 0.566 as ATAC x p300, and H3K27ac 0.531 against 0.513; multiplying p300 by H3K27ac adds nothing in ABC, while in ENCODE-rE2G the combination is the best observed activity term (0.672 against 0.660 for p300, not resolved) (Fig. 7b).

"""
old_summary = keep[s_sum:e_sum]
keep = keep[:s_sum] + NEW_SUMMARY + keep[e_sum:]
# the old summary is replaced wholesale, so only the remaining sections are remapped
head_end = keep.index("## Summary")
s_sum2 = keep.index("## Goals")
body = remap(keep[:head_end], R6_KEEP, R6_TO_R7, "report 7") + keep[head_end:s_sum2] + remap(keep[s_sum2:], R6_KEEP, R6_TO_R7, "report 7")
body = body.replace("This report asks why the transferred model falls short and tests seven kinds of change meant to close the gap, all evaluated on the K562 CRISPR benchmark.",
                    "This report asks why the transferred model falls short and tests seven kinds of change meant to close the gap, all evaluated on the K562 CRISPR benchmark. The changes that give the model TF information, and the analyses behind them, are reported in report 7.", 1)
assert "reported in report 7" in body, "goals anchor missed"
for t in MOVE:
    assert f"## {t}" not in body
Path(str(R6) + a.suffix).write_text(body)
print(f"moved {len(MOVE)} sections ({len(moved):,} chars); report 6 {len(text):,} -> {len(body):,} chars")
for m in sorted(set(re.findall(r"report [67] Figs?\.? \d+[a-z]?", body + moved))):
    print("  cross-ref:", m)
