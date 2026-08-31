#!/usr/bin/env python3
"""Refresh .living/INDEX.md from the source files.

NOT the canonical mycelium indexer. Mycelium's `generate_index.py` is not installed on
Sherlock or on the local machine, and INDEX.md had drifted badly (reporting 24 learnings and
4 decisions against actual counts of 28 and 6), which defeats its purpose as the discovery
surface for a fresh session. This reproduces the existing format from the source files so
the counts are true. If the real indexer is ever run it will simply overwrite this.

Preserves the <!-- BEGIN/END --> marker structure so the canonical tool can still find its
sections.
"""
import os, re, sys, datetime, collections

R = sys.argv[1] if len(sys.argv) > 1 else "."
L = os.path.join(R, ".living")


def mtime(p):
    return datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat()


def entries(path, pat):
    """Return list of (id_index, title) for headings matching pat."""
    out = []
    for line in open(path):
        m = re.match(pat, line)
        if m:
            out.append(m.group(1).strip())
    return out


learn = entries(f"{L}/learnings.md", r"^### \[[\d-]+\]\s*(.+)$")
dec = entries(f"{L}/decisions.md", r"^### \[[\d-]+\]\s*(.+)$")
conv = entries(f"{L}/conventions.md", r"^##\s+(.+)$")

fdir = f"{L}/findings"
ffiles = sorted(f for f in os.listdir(fdir) if f.endswith(".md"))
nfind = 0
for f in ffiles:
    nfind += sum(1 for ln in open(os.path.join(fdir, f)) if ln.startswith("## F-"))

today = datetime.date.today().isoformat()


def topics(items, n=5):
    return ", ".join(items[:n]) if items else "-"


rows = [
    ("conventions.md", f"{len(conv)} sections", mtime(f"{L}/conventions.md"), topics(conv, 4)),
    ("decisions.md", f"{len(dec)} entries", mtime(f"{L}/decisions.md"), topics(dec, 4)),
    ("learnings.md", f"{len(learn)} entries", mtime(f"{L}/learnings.md"), topics(learn, 4)),
    ("findings/", f"{nfind} findings across {len(ffiles)} topics",
     max(mtime(os.path.join(fdir, f)) for f in ffiles),
     ", ".join(f[:-3] for f in ffiles)),
]

quick = ["<!-- BEGIN QUICK REFERENCE -->", "# .living/ Index", f"Last audit: {today}", "",
         "| File | Entries | Last updated | Key topics |",
         "|------|---------|--------------|------------|"]
for r_ in rows:
    quick.append("| %s | %s | %s | %s |" % r_)
quick += ["", "## Local skills", "See `.living/skills/` for project-specific skill packs.",
          "<!-- END QUICK REFERENCE -->"]

# tag clusters, with stable L-n / D-n ids assigned in file order
tags = collections.defaultdict(list)
for label, path in (("L", f"{L}/learnings.md"), ("D", f"{L}/decisions.md")):
    idx = 0
    cur = None
    for line in open(path):
        if re.match(r"^### \[[\d-]+\]", line):
            idx += 1
            cur = f"{label}-{idx}"
        m = re.match(r"^\*\*Tags\*\*:\s*(.+)$", line)
        if m and cur:
            for t in re.split(r"[,\s]+", m.group(1).strip()):
                t = t.strip().strip(".")
                if t:
                    tags[t].append(cur)

summary = ["<!-- BEGIN KNOWLEDGE SUMMARY -->",
           f"Last summarized: {today} (refresh_living_index.py, not the canonical mycelium indexer)",
           "", "## Tag clusters", ""]
for t, ids in sorted(tags.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:18]:
    if len(ids) < 2:
        continue
    summary.append(f"- **{t}** ({len(ids)} entries) — " + ", ".join(ids[:8]))
summary += ["", "## Most recent entries", ""]
for label, items in (("learnings.md", learn), ("decisions.md", dec)):
    for t in items[-3:][::-1]:
        summary.append(f"- *{label}*: {t}")
summary.append("<!-- END KNOWLEDGE SUMMARY -->")

out = "\n".join(quick) + "\n\n" + "\n".join(summary) + "\n"
open(f"{L}/INDEX.md", "w").write(out)
print(f"refreshed {L}/INDEX.md")
print(f"  conventions {len(conv)} | decisions {len(dec)} | learnings {len(learn)} | findings {nfind}")
