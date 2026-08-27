#!/usr/bin/env python3
"""Dependency-tracked propagation for analysis-report .qmd files.

A report is two overlapping DAGs:

  - Semantic DAG  (figure -> legend -> Outcome -> Findings -> headline),
    auto-derived from the .qmd section structure.
  - Computational DAG (one figure's analysis feeds another's data), which lives
    in the analysis scripts, not the .qmd. Declared into the state file as
    per-figure `depends_on` edges + a `command` to regenerate.

Changing one node should propagate to its dependents WITHOUT regenerating the
whole report. Two phases:

  1. Dirty-mark (this script, deterministic, free): a node whose own content or
     file changed since the last `commit` is a source change; every transitive
     dependent is marked needs-update.
  2. Diff-gate (run by Claude, per needs-update node): the worklist hands Claude
     each stale node in topological order with the diffs of what changed
     upstream. Claude decides unchanged (stop) or regenerate. Code nodes that
     carry a `command` become RECOMPUTE actions instead.

Subcommands:
  graph   <report.qmd>   print the derived DAG (inspection)
  status  <report.qmd>   print the needs-update worklist in topo order with diffs
  commit  <report.qmd>   snapshot every node's fingerprint into the state file

State lives next to the report at <report>.report-state.json. The graph itself
is re-derived from the .qmd every run, so the markdown stays the single source
of truth; only fingerprints, snapshots, and the non-derivable computational /
cross-section edges are persisted.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

STATE_SUFFIX = ".report-state.json"
STATE_VERSION = 1

# H2 titles that map to fixed special sections (matched case-insensitively on a
# prefix so "Open questions / next steps" and "Open questions" both hit).
SPECIAL_TITLE_PREFIXES = {
    "headline figure": "headline",
    "summary": "summary",
    "goals": "goals",
    "open questions": "open_questions",
    "methods": "methods",
    "references": "skip",
    "data sources": "skip",
    "quantity definitions": "skip",
    "definitions": "skip",
}

IMAGE_RE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
# A link to a local .html companion — an interactive artifact (e.g. an IGV report).
ARTIFACT_LINK_RE = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+\.html)\)")
LEGEND_RE = re.compile(r"\*\*Figure\s+(\w+)\s*\|", re.IGNORECASE)
QA_RE = re.compile(r"^\*{0,2}([QA]):\*{0,2}\s*(.+)$", re.MULTILINE)
METHOD_RE = re.compile(r"\*\*Methods?\.?\*\*", re.IGNORECASE)
DETAILS_RE = re.compile(r"<details[^>]*>(.*?)</details>", re.DOTALL | re.IGNORECASE)
SUMMARY_SPLIT_RE = re.compile(r"<summary[^>]*>(.*?)</summary>(.*)", re.DOTALL | re.IGNORECASE)
HTML_COMMENT_RE = re.compile(r"<!--.*?-->", re.DOTALL)


@dataclass
class Node:
    id: str
    kind: str  # "code" (figure file) | "llm" (generated prose)
    type: str  # headline_fig | headline_legend | goals | findings | legend |
    #            section_goal | section_methods | section_outcome | figure |
    #            open_questions
    label: str  # human-facing, e.g. "Figure 3 legend"
    inputs: list[str] = field(default_factory=list)
    text: str | None = None  # prose content (llm nodes)
    path: str | None = None  # figure path relative to report dir (code nodes)


@dataclass
class Graph:
    report_path: Path
    nodes: dict[str, Node]

    def topo_order(self) -> list[str]:
        """Kahn's algorithm; inputs point from dependency -> dependent."""
        indegree = {nid: 0 for nid in self.nodes}
        dependents: dict[str, list[str]] = {nid: [] for nid in self.nodes}
        for nid, node in self.nodes.items():
            for dep in node.inputs:
                if dep in self.nodes:
                    indegree[nid] += 1
                    dependents[dep].append(nid)
        ready = sorted(nid for nid, d in indegree.items() if d == 0)
        order: list[str] = []
        while ready:
            nid = ready.pop(0)
            order.append(nid)
            for child in dependents[nid]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    ready.append(child)
                    ready.sort()
        if len(order) != len(self.nodes):
            unresolved = [nid for nid in self.nodes if nid not in order]
            raise ValueError(f"dependency cycle involving: {unresolved}")
        return order


# --------------------------------------------------------------------------- #
# Parsing                                                                      #
# --------------------------------------------------------------------------- #


def strip_comments(text: str) -> str:
    return HTML_COMMENT_RE.sub("", text)


def split_sections(body: str) -> list[tuple[str, str]]:
    """Split markdown body into (h2_title, section_body) pairs, in order."""
    sections: list[tuple[str, str]] = []
    current_title: str | None = None
    current_lines: list[str] = []
    for line in body.splitlines():
        if line.startswith("## "):
            if current_title is not None:
                sections.append((current_title, "\n".join(current_lines)))
            current_title = line[3:].strip()
            current_lines = []
        elif current_title is not None:
            current_lines.append(line)
    if current_title is not None:
        sections.append((current_title, "\n".join(current_lines)))
    return sections


def classify_title(title: str) -> str | None:
    low = title.strip().lower()
    for prefix, kind in SPECIAL_TITLE_PREFIXES.items():
        if low.startswith(prefix):
            return kind
    return None


def extract_paragraphs(section_body: str) -> list[str]:
    """Comment-stripped, blank-line-delimited paragraphs."""
    cleaned = strip_comments(section_body)
    paras = [p.strip() for p in re.split(r"\n\s*\n", cleaned)]
    return [p for p in paras if p]


def first_image(section_body: str) -> str | None:
    match = IMAGE_RE.search(strip_comments(section_body))
    return match.group(1).strip() if match else None


def details_blocks(section_body: str) -> list[tuple[str, str]]:
    """Return ``(summary_label, body)`` for every <details> block in a section."""
    blocks: list[tuple[str, str]] = []
    for inner in DETAILS_RE.findall(section_body):
        m = SUMMARY_SPLIT_RE.search(inner)
        if m:
            blocks.append((m.group(1).strip(), m.group(2).strip()))
        else:
            blocks.append(("", inner.strip()))
    return blocks


def extract_legend(section_body: str) -> str:
    """The legend lives in a <details> block whose summary mentions 'legend'."""
    blocks = details_blocks(section_body)
    for label, body in blocks:
        if "legend" in label.lower():
            return body
    return blocks[0][1] if blocks else ""


def extract_method(section_body: str) -> str:
    """Visible Method bullets plus the collapsed full-methods detail."""
    visible = ""
    m = METHOD_RE.search(section_body)
    if m:
        # From '**Method.**' up to the first <details> (or end of section).
        tail = section_body[m.start():]
        visible = tail.split("<details", 1)[0].strip()
    detail = ""
    for label, body in details_blocks(section_body):
        if "method" in label.lower():
            detail = body
            break
    return (visible + "\n\n" + detail).strip()


def extract_qa(section_body: str) -> str:
    """The Q: and A: lines (without the figure/legend/method below them)."""
    lines = []
    for tag, text in QA_RE.findall(section_body):
        lines.append(f"{tag}: {text.strip()}")
    return "\n".join(lines)


def figure_number(section_body: str, fallback: str) -> str:
    m = LEGEND_RE.search(section_body)
    return m.group(1) if m else fallback


def parse_report(report_path: Path) -> Graph:
    raw = report_path.read_text()
    # Drop YAML front matter.
    body = raw
    if raw.startswith("---"):
        end = raw.find("\n---", 3)
        if end != -1:
            body = raw[end + 4 :]

    nodes: dict[str, Node] = {}
    sections = split_sections(body)
    analysis_index = 0

    for title, section_body in sections:
        kind = classify_title(title)
        clean = strip_comments(section_body)
        image = first_image(section_body)

        if kind == "headline":
            if image:
                nodes["headline_fig"] = Node(
                    "headline_fig", "code", "headline_fig", "Headline figure",
                    path=image,
                )
            nodes["headline_legend"] = Node(
                "headline_legend", "llm", "headline_legend", "Headline legend",
                inputs=(["headline_fig"] if image else []) + ["summary"],
                text=extract_legend(clean),
            )
        elif kind == "summary":
            nodes["summary"] = Node(
                "summary", "llm", "summary", "Summary",
                text=clean.strip(),
            )  # inputs (all section_qa) filled after all sections parsed
        elif kind == "goals":
            nodes["goals"] = Node(
                "goals", "llm", "goals", "Goals", text=clean.strip(),
            )
        elif kind == "open_questions":
            nodes["open_questions"] = Node(
                "open_questions", "llm", "open_questions", "Open questions",
                inputs=["summary"], text=clean.strip(),
            )
        elif kind in ("methods", "skip"):
            continue  # static / framing; not part of the propagation graph
        else:
            # An analysis slide. Its source artifact is either a static figure
            # or a local interactive-artifact HTML (e.g. an IGV report).
            artifact = image
            is_interactive = False
            if artifact is None:
                m = ARTIFACT_LINK_RE.search(clean)
                if m:
                    artifact, is_interactive = m.group(1).strip(), True
                else:
                    continue  # a prose framing section — not a slide

            analysis_index += 1
            k = analysis_index
            base = f"sec{k}"
            label = f"Interactive {k}" if is_interactive else f"Figure {figure_number(clean, str(k))}"
            nodes[f"{base}_fig"] = Node(
                f"{base}_fig", "code", "artifact" if is_interactive else "figure",
                label, path=artifact,
            )
            fig_in = [f"{base}_fig"]
            nodes[f"{base}_qa"] = Node(
                f"{base}_qa", "llm", "section_qa",
                f"{label} Q/A ({title[:32]})", inputs=fig_in, text=extract_qa(clean),
            )
            nodes[f"{base}_method"] = Node(
                f"{base}_method", "llm", "section_method", f"{label} method",
                inputs=fig_in, text=extract_method(clean),
            )
            # Static figures carry a folded Nature legend; interactive artifacts
            # are self-describing and have no separate legend node.
            if not is_interactive:
                nodes[f"{base}_legend"] = Node(
                    f"{base}_legend", "llm", "legend", f"{label} legend",
                    inputs=fig_in, text=extract_legend(clean),
                )

    # Summary depends on every slide's Q/A (the per-figure answers it synthesizes).
    qa_ids = [nid for nid, n in nodes.items() if n.type == "section_qa"]
    if "summary" in nodes:
        nodes["summary"].inputs = sorted(qa_ids)
    # Drop edges to nodes that don't exist (e.g. headline when Summary is omitted).
    for node in nodes.values():
        node.inputs = [dep for dep in node.inputs if dep in nodes]

    return Graph(report_path, nodes)


# --------------------------------------------------------------------------- #
# State + fingerprints                                                         #
# --------------------------------------------------------------------------- #


def state_path_for(report_path: Path) -> Path:
    return report_path.with_name(report_path.name + STATE_SUFFIX)


def load_state(report_path: Path) -> dict:
    path = state_path_for(report_path)
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def empty_state() -> dict:
    return {
        "version": STATE_VERSION,
        "nodes": {},  # id -> {fingerprint, snapshot}
        "computational_edges": {},  # fig_node_id -> [upstream_node_id, ...]
        "commands": {},  # fig_node_id -> shell command to regenerate
        "extra_semantic_edges": {},  # node_id -> [upstream_node_id, ...]
    }


def apply_declared_edges(graph: Graph, state: dict) -> None:
    """Merge non-derivable edges (computational + cross-section) onto the graph."""
    for node_id, ups in state.get("computational_edges", {}).items():
        if node_id in graph.nodes:
            for up in ups:
                if up in graph.nodes and up not in graph.nodes[node_id].inputs:
                    graph.nodes[node_id].inputs.append(up)
    for node_id, ups in state.get("extra_semantic_edges", {}).items():
        if node_id in graph.nodes:
            for up in ups:
                if up in graph.nodes and up not in graph.nodes[node_id].inputs:
                    graph.nodes[node_id].inputs.append(up)


def fingerprint(graph: Graph, node: Node) -> str:
    if node.kind == "code":
        target = (graph.report_path.parent / node.path).resolve() if node.path else None
        if not target or not target.exists():
            return f"MISSING:{node.path}"
        return "sha256:" + hashlib.sha256(target.read_bytes()).hexdigest()
    return "sha256:" + hashlib.sha256((node.text or "").encode()).hexdigest()


def snapshot(node: Node) -> str:
    """Human-readable content stored in state to diff against next run."""
    if node.kind == "code":
        return f"[figure file] {node.path}"
    return node.text or ""


# --------------------------------------------------------------------------- #
# Propagation                                                                  #
# --------------------------------------------------------------------------- #


def compute_worklist(graph: Graph, state: dict) -> dict:
    """Return changed sources + topo-ordered needs-update nodes."""
    stored = state.get("nodes", {})
    order = graph.topo_order()

    current_fp = {nid: fingerprint(graph, graph.nodes[nid]) for nid in graph.nodes}
    has_baseline = bool(stored)

    # A source change: the node's own file/content moved since last commit, or
    # it is brand new (no baseline for it).
    changed_sources = set()
    for nid in graph.nodes:
        prior = stored.get(nid)
        if prior is None or prior.get("fingerprint") != current_fp[nid]:
            changed_sources.add(nid)

    # Propagate: a node needs update if any input changed or itself needs update.
    needs_update: set[str] = set()
    for nid in order:
        if any(dep in changed_sources or dep in needs_update
               for dep in graph.nodes[nid].inputs):
            needs_update.add(nid)

    return {
        "has_baseline": has_baseline,
        "order": order,
        "current_fp": current_fp,
        "changed_sources": changed_sources,
        "needs_update": needs_update,
    }


def diff_text(old: str, new: str) -> str:
    lines = list(difflib.unified_diff(
        old.splitlines(), new.splitlines(),
        fromfile="before", tofile="after", lineterm="",
    ))
    return "\n".join(lines) if lines else "(no textual diff)"


# --------------------------------------------------------------------------- #
# Subcommands                                                                  #
# --------------------------------------------------------------------------- #


def cmd_graph(report_path: Path) -> int:
    state = load_state(report_path)
    graph = parse_report(report_path)
    apply_declared_edges(graph, state)
    order = graph.topo_order()
    print(f"# Dependency graph for {report_path.name}\n")
    print(f"{len(graph.nodes)} nodes (topological order):\n")
    commands = state.get("commands", {})
    for nid in order:
        node = graph.nodes[nid]
        deps = ", ".join(node.inputs) if node.inputs else "—(root)"
        tag = "CODE" if node.kind == "code" else "llm "
        print(f"  [{tag}] {nid:<16} ({node.label})")
        print(f"          depends on: {deps}")
        if node.kind == "code" and nid in commands:
            print(f"          regenerate: {commands[nid]}")
    return 0


def cmd_status(report_path: Path) -> int:
    state = load_state(report_path)
    graph = parse_report(report_path)
    apply_declared_edges(graph, state)
    work = compute_worklist(graph, state)
    stored = state.get("nodes", {})
    commands = state.get("commands", {})

    if not work["has_baseline"]:
        print("No baseline state found.")
        print(f"Run:  python report_graph.py commit {report_path}")
        print("to record the current report as the propagation baseline.")
        return 0

    needs = [nid for nid in work["order"] if nid in work["needs_update"]]
    # Source changes the user already made (not themselves needing regeneration).
    pure_sources = sorted(work["changed_sources"] - work["needs_update"])

    if pure_sources:
        print("Changed since last commit (source edits — already up to date):")
        for nid in pure_sources:
            print(f"  • {graph.nodes[nid].label}  [{nid}]")
        print()

    if not needs:
        print("Nothing downstream needs updating. Report is in sync.")
        return 0

    print(f"{len(needs)} node(s) need updating, in propagation order:\n")
    for nid in needs:
        node = graph.nodes[nid]
        if node.kind == "code":
            if nid in commands:
                action = f"RECOMPUTE — run: {commands[nid]}"
            else:
                action = "MANUAL RERUN — regenerate this figure, no command declared"
        else:
            action = "DIFF-GATE — decide unchanged vs regenerate (see upstream changes)"
        print(f"── {node.label}  [{nid}]")
        print(f"   action: {action}")

        changed_inputs = [dep for dep in node.inputs
                          if dep in work["changed_sources"] or dep in work["needs_update"]]
        for dep in changed_inputs:
            dep_node = graph.nodes[dep]
            if dep in work["needs_update"]:
                print(f"   ← {dep_node.label} [{dep}] (awaiting its update above)")
            elif dep_node.kind == "code":
                print(f"   ← {dep_node.label} [{dep}]: figure file changed")
            else:
                old = stored.get(dep, {}).get("snapshot", "")
                new = dep_node.text or ""
                print(f"   ← {dep_node.label} [{dep}] changed:")
                for line in diff_text(old, new).splitlines():
                    print(f"        {line}")
        print()

    print("After updating the prose / rerunning figures, record the new baseline:")
    print(f"  python report_graph.py commit {report_path}")
    return 0


def cmd_commit(report_path: Path) -> int:
    state = load_state(report_path) or empty_state()
    state.setdefault("version", STATE_VERSION)
    for key in ("computational_edges", "commands", "extra_semantic_edges"):
        state.setdefault(key, {})
    graph = parse_report(report_path)
    apply_declared_edges(graph, state)

    new_nodes: dict[str, dict] = {}
    for nid, node in graph.nodes.items():
        new_nodes[nid] = {
            "fingerprint": fingerprint(graph, node),
            "snapshot": snapshot(node),
        }
    state["nodes"] = new_nodes

    out = state_path_for(report_path)
    out.write_text(json.dumps(state, indent=2) + "\n")
    print(f"Baseline recorded for {len(new_nodes)} nodes → {out.name}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("graph", "status", "commit"):
        sp = sub.add_parser(name)
        sp.add_argument("report", type=Path)
    args = parser.parse_args(argv)

    report_path = args.report
    if not report_path.exists():
        print(f"error: report not found: {report_path}", file=sys.stderr)
        return 2

    return {"graph": cmd_graph, "status": cmd_status, "commit": cmd_commit}[args.command](report_path)


if __name__ == "__main__":
    raise SystemExit(main())
