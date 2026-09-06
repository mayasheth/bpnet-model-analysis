#!/usr/bin/env python3
"""Split h3k27ac_model_report.qmd into three reports ordered by churn rate.

Rationale is recorded in .living/decisions.md [2026-09-03]. Reports 1 and 2 are
citable references; report 3 is the workbench and is expected to change weekly.

Sections that move unchanged are pulled from the source by title rather than
retyped, so their legends and numbers cannot drift during the split. Figure
numbers restart per report, remapped through a placeholder pass so that a chain
of substitutions cannot collide.
"""
import re, sys
from pathlib import Path

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2])
text = SRC.read_text()

CSS = re.search(r"```\{=html\}\n<style>.*?</style>\n```", text, re.S).group(0)

def sections(t):
    """title -> body, for every H2 in the source."""
    heads = list(re.finditer(r"^##\s+(.+?)\s*$", t, re.M))
    out = {}
    for i, h in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(t)
        out[h.group(1).strip()] = t[h.end():end].strip("\n")
    return out

OLD = sections(text)

def pull(title, mapping, new_title=None, fixups=()):
    """Lift a section from the source, remapping its figure numbers.

    Renumbering applies ONLY to pulled text. New prose is written with this
    report's final figure numbers already, so running the remapper over it
    would map them a second time.
    """
    if title not in OLD:
        raise SystemExit(f"source section not found: {title!r}")
    body = OLD[title]
    for pat, rep in fixups:
        body = re.sub(pat, rep, body)
    return f"## {new_title or title}\n\n{renumber(body, mapping)}\n"

NUM = re.compile(r"\b(Figs?\.?|Figures?)(\s+)(\d+(?:\s*,\s*\d+)*)")

def renumber(body, mapping):
    """Remap figure numbers via placeholders so substitutions cannot chain."""
    def sub(m):
        nums = [n.strip() for n in m.group(3).split(",")]
        mapped = []
        for n in nums:
            if int(n) not in mapping:
                raise SystemExit(f"figure {n} referenced but not in this report's map:\n"
                                 f"  ...{m.string[max(0, m.start()-90):m.end()+40]}...")
            mapped.append(f"@@{mapping[int(n)]}@@")
        return f"{m.group(1)}{m.group(2)}{', '.join(mapped)}"
    body = NUM.sub(sub, body)
    body = re.sub(r"@@(\d+)@@", r"\1", body)
    # Cross-report references are written "Report 1, Fig.~2" so NUM cannot see
    # them; restore the space once this report's own numbers are settled.
    return re.sub(r"\b(Figs?\.?|Figures?)~", r"\1 ", body)

def front(title, date, status="draft"):
    return (f'---\ntitle: "{title}"\ndate: "{date}"\nstatus: "{status}"\n'
            "format:\n  html:\n    embed-resources: true\n    toc: true\n"
            "    code-fold: true\n---\n\n" + CSS + "\n")

def write(name, title, date, parts):
    body = "\n".join(parts)
    body = re.sub(r"\b(Figs?\.?|Figures?)~", r"\1 ", body)
    p = OUT / name
    p.write_text(front(title, date) + "\n" + body)
    figs = len(re.findall(r"^!\[\]\(", body, re.M))
    print(f"wrote {p.name}: {len(body.splitlines())} lines, {figs} images, "
          f"{len(re.findall(r'^## ', body, re.M))} sections")

import report_text as R
R.build(pull, write)
