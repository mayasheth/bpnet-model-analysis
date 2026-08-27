#!/usr/bin/env python3
"""Lint, stale-figure check, and render an analysis report to self-contained HTML.

Usage
-----
    python render_report.py path/to/report.qmd          # lint + render
    python render_report.py path/to/report.qmd --no-render  # lint only
    python render_report.py path/to/report.qmd --no-lint    # render only
    python render_report.py path/to/report.qmd --strict     # treat warnings as errors

What it checks
--------------
1. Broken image links — every ``![](path)`` resolves to a file on disk
2. Empty required sections — Goals, Methods (Summary is optional)
3. Slide structure — each figure slide carries Q:/A: lines, a <details> legend,
   and a Method block (warning; errors under --strict)
4. Placeholder leftovers — ``TODO``, ``<placeholder>``-style angle-bracketed tokens
5. Stale prose — any data figure with mtime newer than the report's own mtime
   (concept-art schematics are exempt)
6. Regenerate environment — a Regenerate block that runs analysis code but pins
   no environment (lockfile / env name / hash) warns (errors under --strict)
7. Number traceability — statistical claims in prose (n=, p=, FDR<, r=, %) are
   checked against mycelium's register_value numbers.json when present
   (untraceable claims warn); otherwise a soft nudge to register them
8. Style (soft, never fails) — over-long Q:/A: lines, too many Method bullets,
   over-long Summary, hype words banned by lander_voice.md

Diagnostics come in three tiers: errors (always fail), warnings (fail under
``--strict``), and style notes (never fail).

Then runs ``quarto render`` (or falls back to ``pandoc --standalone
--embed-resources`` if quarto isn't available).
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


# Sections that must exist and must contain non-whitespace content beyond
# the heading and HTML comments. The "Headline figure" section is required
# but checked separately because its body is an image rather than prose.
# "Summary" is intentionally NOT required — it is cross-cutting synthesis that
# is dropped when the work is a single thread.
REQUIRED_SECTIONS = ("Goals", "Methods")
HEADLINE_SECTION = "Headline figure"

# H2 titles that are framing/prose sections, not figure "slides". Any other H2
# that contains an image is treated as an analysis slide and checked for the
# Q:/A:/legend/Method structure.
NON_SLIDE_TITLES = {
    "headline figure", "summary", "goals", "open questions / next steps",
    "open questions", "methods", "references", "data sources",
    "quantity definitions", "definitions",
}

# Figures whose path contains one of these markers are hand-authored concept
# cartoons that don't regenerate with the data — exempt from the stale check.
CONCEPT_ART_MARKERS = ("schematic", "cartoon", "concept")

# Hype words the Lander voice bans (lander_voice.md). Soft warning only.
HYPE_WORDS = (
    "novel", "powerful", "robust", "seamless", "cutting-edge",
    "state-of-the-art", "striking", "elegant", "remarkable",
)

# Regex for markdown image links: ![alt](path)
IMAGE_RE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")

# A link to a local .html companion (an interactive artifact, e.g. an IGV report).
# Negative lookbehind on '!' so it doesn't match image links.
ARTIFACT_LINK_RE = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+\.html)\)")

# Regex for ATX headings (## Section)
HEADING_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*$", re.MULTILINE)

# Anything that looks like an unfilled placeholder.
# Matches:
#   - lines containing TODO, FIXME, XXX (case-insensitive, word-boundary)
#   - <PLACEHOLDER_LIKE_THIS> tokens (uppercase + underscores, in angle brackets)
#   - ``YYYY-MM-DD`` literal in places it shouldn't be (rendered) — checked
#     separately because it's expected inside the regenerate command example
PLACEHOLDER_TODO_RE = re.compile(r"\b(TODO|FIXME|XXX)\b", re.IGNORECASE)
PLACEHOLDER_ANGLE_RE = re.compile(r"<([A-Z][A-Z0-9_]{2,}|[a-z][a-z0-9_]*\s+[a-z0-9_\s]+)>")
PLACEHOLDER_DATE_RE = re.compile(r"\bYYYY-MM-DD\b")

# Status -> (color, label) for the banner. Picked for accessibility.
STATUS_STYLES = {
    "draft": ("#fff3cd", "#664d03", "Draft — work in progress"),
    "in-review": ("#cfe2ff", "#084298", "In review — feedback wanted"),
    "final": ("#d1e7dd", "#0f5132", "Final"),
    "superseded": ("#e2e3e5", "#41464b", "Superseded — see newer report"),
}


@dataclass
class LintResult:
    """Accumulated diagnostics from a lint pass.

    Three tiers (the F2 policy):
    - ``errors``        — always fail.
    - ``warnings``      — structural gaps; become errors under ``--strict``
                          (the finalize path).
    - ``soft_warnings`` — style heuristics (prose length, hype words); printed
                          but NEVER fail, even under ``--strict``.
    """

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    soft_warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def report(self, strict: bool) -> int:
        """Print diagnostics and return the exit code (0 if clean)."""
        for e in self.errors:
            print(f"ERROR: {e}", file=sys.stderr)
        for w in self.warnings:
            print(f"warning: {w}", file=sys.stderr)
        for s in self.soft_warnings:
            print(f"style: {s}", file=sys.stderr)
        if self.errors:
            return 1
        if strict and self.warnings:
            return 1
        return 0


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def split_front_matter(text: str) -> tuple[str, str]:
    """Return ``(yaml_block, body)``. ``yaml_block`` is empty if no front matter."""
    if not text.startswith("---\n"):
        return "", text
    end = text.find("\n---\n", 4)
    if end == -1:
        return "", text
    return text[4:end], text[end + 5:]


def strip_html_comments(text: str) -> str:
    """Remove ``<!-- ... -->`` blocks for content checks."""
    return re.sub(r"<!--.*?-->", "", text, flags=re.DOTALL)


def strip_code_blocks(text: str) -> str:
    """Remove fenced code blocks for placeholder/heading checks.

    Code blocks legitimately contain TODO, <placeholder>, and YYYY-MM-DD as
    examples; we don't want to flag them.
    """
    return re.sub(r"```.*?```", "", text, flags=re.DOTALL)


def find_section_body(text: str, section_name: str) -> str | None:
    """Return the body of a section heading (any level), or None if missing.

    Body runs from the heading to the next heading of the same or higher level.
    HTML comments and code blocks are stripped from the returned body.
    """
    pattern = re.compile(
        rf"^(#{{1,6}})\s+{re.escape(section_name)}\s*$",
        re.MULTILINE,
    )
    m = pattern.search(text)
    if not m:
        return None
    level = len(m.group(1))
    start = m.end()
    # Find next heading of same-or-higher level
    next_heading = re.compile(rf"^#{{1,{level}}}\s+\S", re.MULTILINE)
    nxt = next_heading.search(text, pos=start)
    body = text[start:nxt.start()] if nxt else text[start:]
    return strip_html_comments(body).strip()


# ---------------------------------------------------------------------------
# Lint checks
# ---------------------------------------------------------------------------


def check_required_sections(body_no_comments: str, result: LintResult) -> None:
    for section in REQUIRED_SECTIONS:
        body = find_section_body(body_no_comments, section)
        if body is None:
            result.errors.append(f"Required section missing: '## {section}'")
        elif not body:
            result.errors.append(
                f"Section '{section}' is empty — fill it in or remove the heading"
            )


def check_headline_figure(body_no_comments: str, result: LintResult) -> None:
    """The headline figure must be the first H2 section and contain an image.

    Two checks:
    1. A '## Headline figure' section exists and contains at least one image.
    2. The first H2 section in the document IS '## Headline figure', so a
       skim reader sees the headline image before any other section.
    """
    body = find_section_body(body_no_comments, HEADLINE_SECTION)
    if body is None:
        result.errors.append(
            f"Required section missing: '## {HEADLINE_SECTION}' "
            "— the report must lead with a headline figure"
        )
        return
    if not IMAGE_RE.search(body):
        result.errors.append(
            f"Section '## {HEADLINE_SECTION}' has no image link "
            "— add the headline figure (![](path))"
        )

    # Check that this is the first H2-or-deeper heading in the body.
    first_heading = HEADING_RE.search(body_no_comments)
    if first_heading and first_heading.group(2).strip() != HEADLINE_SECTION:
        result.errors.append(
            f"'## {HEADLINE_SECTION}' must be the first section in the report "
            f"(found '## {first_heading.group(2).strip()}' first). "
            "Move the headline figure to the top so a skim reader sees it before any prose."
        )


def check_regenerate_in_methods(body_no_comments: str, result: LintResult) -> None:
    """The Regenerate command lives at the bottom, inside the Methods section.

    Warns if a top-level '## Regenerate' heading exists (old-style report) or
    if there is no '### Regenerate' subsection inside Methods.
    """
    # Old-style top-level Regenerate section
    top_level = re.search(r"^##\s+Regenerate\s*$", body_no_comments, re.MULTILINE)
    if top_level:
        result.errors.append(
            "'## Regenerate' should be a '### Regenerate' subsection inside '## Methods', "
            "not a top-level section. Move it to the bottom of the report."
        )
    methods_body = find_section_body(body_no_comments, "Methods")
    if methods_body is not None:
        if not re.search(r"^###\s+Regenerate\s*$", methods_body, re.MULTILINE):
            result.warnings.append(
                "No '### Regenerate' subsection inside '## Methods'. "
                "Add the regenerate command(s) so collaborators can rebuild the report."
            )
        else:
            check_regenerate_environment(methods_body, result)


# Commands that run analysis code — their presence means a bare command was
# recorded but the environment that ran it may not have been.
RUN_COMMAND_RE = re.compile(
    r"\b(pixi run|mise exec|Rscript|jupyter|conda run|uv run|poetry run|python\d?)\b"
)
# Evidence that the environment is pinned, not just the command. A lockfile
# reference, an explicit env name, or a content hash all satisfy this.
ENV_FINGERPRINT_RE = re.compile(
    r"pixi\.lock|renv\.lock|conda-lock|environment\.ya?ml|requirements\.txt"
    r"|\bsha256\b|\bEnvironment:\s|\benv\b.*@|pixi\.toml",
    re.IGNORECASE,
)


def check_regenerate_environment(regenerate_body: str, result: LintResult) -> None:
    """A regenerate command is only reproducible if the environment that ran it
    is pinned too. If the Regenerate block runs analysis code but records no
    environment fingerprint (lockfile, env name, or content hash), warn.

    Captures the environment, not just the command — a bare `pixi run …` with no
    `pixi.lock` reference reproduces a different result once the env drifts.
    """
    # Only relevant if there is an actual run command to reproduce.
    if not RUN_COMMAND_RE.search(regenerate_body):
        return
    if not ENV_FINGERPRINT_RE.search(regenerate_body):
        result.warnings.append(
            "Regenerate records a command but no environment fingerprint. "
            "Add an 'Environment:' line pinning the env + lockfile so the command "
            "reproduces the same result later, e.g.:\n"
            "    Environment: pixi env `default`, pixi.lock @ sha256 "
            "`$(shasum -a 256 pixi.lock | cut -c1-12)`"
        )


def iter_h2_sections(body_no_comments: str):
    """Yield ``(title, section_body)`` for every H2 heading, in order."""
    heads = list(re.finditer(r"^##\s+(.+?)\s*$", body_no_comments, re.MULTILINE))
    for i, h in enumerate(heads):
        start = h.end()
        end = heads[i + 1].start() if i + 1 < len(heads) else len(body_no_comments)
        yield h.group(1).strip(), body_no_comments[start:end]


SENTENCE_ABBREVIATIONS = (
    "Fig.", "Figs.", "Eq.", "Eqs.", "vs.", "e.g.", "i.e.", "et al.",
    "approx.", "ca.", "cf.", "No.", "Suppl.", "ref.", "Dr.",
)


def count_sentences(text: str) -> int:
    """Rough sentence count — terminal . ! ? not inside a number or abbreviation."""
    stripped = re.sub(r"\b\d+\.\d+\b", "", text)  # don't count decimals
    for abbr in SENTENCE_ABBREVIATIONS:
        stripped = stripped.replace(abbr, abbr.replace(".", ""))
    return len(re.findall(r"[.!?](?:\s|$)", stripped)) or (1 if text.strip() else 0)


def check_slides(body_no_comments: str, result: LintResult) -> None:
    """Each analysis slide (an H2 with a figure, not a framing section) must
    carry Q:/A: lines, a collapsed <details> legend, and a Method block.

    Missing structure → warning (errors under --strict). Prose-length issues →
    soft warning (never fails).
    """
    for title, sect in iter_h2_sections(body_no_comments):
        if title.lower() in NON_SLIDE_TITLES:
            continue
        has_image = bool(IMAGE_RE.search(sect))
        is_interactive = bool(ARTIFACT_LINK_RE.search(sect))
        if not has_image and not is_interactive:
            continue  # a prose section, not a slide

        q = re.search(r"^\*{0,2}Q:\*{0,2}\s*(.+)$", sect, re.MULTILINE)
        a = re.search(r"^\*{0,2}A:\*{0,2}\s*(.+)$", sect, re.MULTILINE)
        if not q:
            result.warnings.append(f"Slide '{title}': missing a '**Q:**' question line")
        if not a:
            result.warnings.append(f"Slide '{title}': missing an '**A:**' answer line")
        # A figure slide needs its legend folded into <details>; an interactive
        # artifact slide (no static figure) is self-describing and exempt.
        if has_image and "<details" not in sect.lower():
            result.warnings.append(
                f"Slide '{title}': no collapsed <details> legend — "
                "put the full Nature legend in a <details> block"
            )
        if not re.search(r"\*\*Method", sect):
            result.warnings.append(
                f"Slide '{title}': no '**Method.**' block — add 1–3 method bullets"
            )

        # Soft (never-fail) style checks.
        if q and count_sentences(q.group(1)) > 1:
            result.soft_warnings.append(
                f"Slide '{title}': Q: should be one sentence")
        if a and count_sentences(a.group(1)) > 1:
            result.soft_warnings.append(
                f"Slide '{title}': A: should be one sentence")
        bullets = re.findall(r"^\s*[-*]\s+\S", sect, re.MULTILINE)
        # Count only bullets before the first <details> (the visible Method bullets).
        visible = sect.split("<details")[0]
        visible_bullets = re.findall(r"^\s*[-*]\s+\S", visible, re.MULTILINE)
        if len(visible_bullets) > 3:
            result.soft_warnings.append(
                f"Slide '{title}': {len(visible_bullets)} visible Method bullets "
                "(keep to 1–3; move detail into the <details> block)")


def check_interactive_artifacts(body: str, report_dir: Path, result: LintResult) -> None:
    """Every local interactive-artifact link (e.g. an IGV report HTML) must resolve.

    These companion files are how the report links out to an interactive view
    while staying self-contained itself. A broken link ships a dead button.
    """
    for match in ARTIFACT_LINK_RE.finditer(body):
        link = match.group(1).strip()
        if link.startswith(("http://", "https://", "data:", "#")):
            continue
        target = (report_dir / link).resolve()
        if not target.exists():
            result.warnings.append(
                f"Broken interactive-artifact link: '{link}' — no file at {target}. "
                "Build it (e.g. igv build_igv_report.py) or fix the path."
            )


def check_summary(body_no_comments: str, result: LintResult) -> None:
    """Summary is optional, but if present should be a short synthesis."""
    body = find_section_body(body_no_comments, "Summary")
    if body is None:
        return
    bullets = re.findall(r"^\s*[-*]\s+\S", body, re.MULTILINE)
    if len(bullets) > 5:
        result.soft_warnings.append(
            f"Summary has {len(bullets)} bullets — keep cross-cutting synthesis to ~5")


# --- Number traceability (bridge to mycelium's register_value manifest) -----
#
# scitexlintr keeps a LaTeX report honest by wrapping every reportable number in
# a \SciVal macro checked against a manifest. A Quarto report has no macros —
# numbers are typed straight into prose — so we adapt the two scitexlintr rules
# that survive without wrappers:
#   - handwritten-numeric-claim (manifest-free): flag hand-typed statistical
#     claims (n = 48, p = 1e-8, FDR < 0.05, r = 0.72, 96.5%) wherever they sit.
#   - unsourced-numeric-token (manifest-backed): if mycelium's numbers.json is
#     present, a flagged claim whose value IS in the manifest is traceable and
#     stays silent; one that is NOT is an untraceable number to register/cite.
# The claim detector bounds what we inspect (high precision); the manifest tells
# us which of those claims are backed.

# Statistical-claim patterns. Each captures the numeric literal in group 1.
NUMERIC_CLAIM_RES = (
    re.compile(r"\b[nN]\s*=\s*([\d,]+(?:\.\d+)?)"),                       # n = 48
    re.compile(r"\b[pPqQ]\s*[=<>]\s*([\d.]+(?:[eE][-+]?\d+)?)"),          # p = 1e-8, q < 0.01
    re.compile(r"\b(?:FDR|alpha|α)\s*[=<>]\s*([\d.]+(?:[eE][-+]?\d+)?)"), # FDR < 0.05
    re.compile(r"\b(?:r|ρ|rho|R\^?2|R²|AUC|AUROC)\s*=\s*(-?[\d.]+)"),     # r = 0.72
    re.compile(r"([\d.]+)\s*%"),                                          # 96.5%
    re.compile(r"([\d.]+)\s*(?:-?\s*fold\b|×|-fold)"),                    # 3.2-fold
)


def find_numbers_manifest(report_dir: Path, explicit: Path | None) -> Path | None:
    """Locate mycelium's register_value output for this analysis.

    Search order: an explicit ``--numbers`` path, then the conventional
    locations relative to the report — ``outputs/numbers.json`` and the
    enriched ``reports/.manifest.json`` — walking up so an
    ``analysis/<name>/reports/report.qmd`` finds ``analysis/<name>/outputs/``.
    """
    if explicit is not None:
        return explicit if explicit.exists() else None
    candidates = []
    for base in (report_dir, *list(report_dir.parents)[:3]):
        candidates += [base / "outputs" / "numbers.json",
                       base / ".manifest.json",
                       base / "reports" / ".manifest.json"]
    for c in candidates:
        if c.exists():
            return c
    return None


def load_manifest_values(manifest_path: Path) -> set[float]:
    """Return the numeric values registered in a numbers.json fragment or an
    enriched .manifest.json. Non-numeric (str) values are ignored — they can't
    collide with a numeric claim."""
    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return set()
    entries = data.get("values") or data.get("numbers") or []
    out: set[float] = set()
    for e in entries:
        if not isinstance(e, dict):
            continue
        v = e.get("value")
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            out.add(float(v))
        elif isinstance(v, str):
            try:
                out.add(float(v.replace(",", "")))
            except ValueError:
                pass
    return out


def token_matches_manifest(token: str, values: set[float]) -> bool:
    """A prose number is traceable if it (or its percent-scaled forms) matches a
    registered value. Matching is deliberately generous — a false 'traceable'
    is quieter than a false 'untraceable' nag."""
    try:
        x = float(token.replace(",", ""))
    except ValueError:
        return False
    for v in values:
        for cand in (x, x / 100.0, x * 100.0):
            if v == 0:
                if abs(cand) < 1e-12:
                    return True
            elif abs(cand - v) <= 1e-6 * abs(v):
                return True
    return False


def check_number_traceability(
    body_no_comments: str,
    report_dir: Path,
    result: LintResult,
    numbers_path: Path | None = None,
) -> None:
    """Flag quantitative claims in prose. Manifest-backed when numbers.json is
    found (untraceable claims → warning), manifest-free otherwise (all claims →
    soft nudge to register them)."""
    prose = strip_code_blocks(body_no_comments)
    claims: list[str] = []
    seen: set[str] = set()
    for rx in NUMERIC_CLAIM_RES:
        for m in rx.finditer(prose):
            tok = m.group(1)
            ctx = m.group(0).strip()
            if ctx not in seen:
                seen.add(ctx)
                claims.append(tok if tok == ctx else ctx)
    if not claims:
        return

    manifest = find_numbers_manifest(report_dir, numbers_path)
    if manifest is None:
        result.soft_warnings.append(
            "Hand-typed numeric claims found but no register_value manifest "
            "(numbers.json) to check them against: "
            + ", ".join(sorted(seen))
            + ". Register these with register_value() so the report can't drift, "
            "or pass --numbers PATH if the manifest lives elsewhere."
        )
        return

    values = load_manifest_values(manifest)
    untraceable = []
    for rx in NUMERIC_CLAIM_RES:
        for m in rx.finditer(prose):
            if not token_matches_manifest(m.group(1), values):
                untraceable.append(m.group(0).strip())
    untraceable = sorted(set(untraceable))
    if untraceable:
        result.warnings.append(
            f"Untraceable numeric claims (not in {manifest.name}): "
            + ", ".join(untraceable)
            + ". Register each with register_value() at its compute site, or "
            "cite its source — an unbacked number silently drifts when the "
            "analysis re-runs."
        )


def check_hype_words(body_no_comments: str, result: LintResult) -> None:
    cleaned = strip_code_blocks(body_no_comments).lower()
    found = sorted({w for w in HYPE_WORDS if re.search(rf"\b{w}\b", cleaned)})
    if found:
        result.soft_warnings.append(
            "Hype words the Lander voice bans (lander_voice.md): "
            + ", ".join(found))


def check_image_links(
    body: str,
    report_dir: Path,
    result: LintResult,
) -> list[tuple[str, Path]]:
    """Verify every image link resolves to a file. Returns (link, abspath) for resolved ones."""
    resolved: list[tuple[str, Path]] = []
    for match in IMAGE_RE.finditer(body):
        link = match.group(1).strip()
        # Skip URLs (http/https) and data URIs
        if link.startswith(("http://", "https://", "data:")):
            continue
        target = (report_dir / link).resolve()
        if not target.exists():
            result.errors.append(
                f"Broken image link: '{link}' — no file at {target}"
            )
        else:
            resolved.append((link, target))
    return resolved


def check_placeholders(body: str, result: LintResult) -> None:
    """Flag TODO/FIXME, angle-bracketed placeholders, YYYY-MM-DD outside front matter."""
    # Strip code blocks first — placeholders inside them are intentional examples
    cleaned = strip_code_blocks(body)

    todo_lines = [
        (i, line)
        for i, line in enumerate(cleaned.splitlines(), 1)
        if PLACEHOLDER_TODO_RE.search(line)
    ]
    for lineno, line in todo_lines:
        result.warnings.append(f"Line {lineno}: leftover TODO/FIXME — '{line.strip()[:80]}'")

    for match in PLACEHOLDER_ANGLE_RE.finditer(cleaned):
        token = match.group(0)
        result.warnings.append(f"Possible unfilled placeholder: '{token}'")

    if PLACEHOLDER_DATE_RE.search(cleaned):
        result.warnings.append(
            "Literal 'YYYY-MM-DD' appears outside the front matter or code block "
            "— probably a placeholder you meant to fill in"
        )


def check_stale_figures(
    report_path: Path,
    images: list[tuple[str, Path]],
    result: LintResult,
) -> None:
    """Warn about figures newer than the report itself."""
    if not images:
        return
    report_mtime = report_path.stat().st_mtime
    stale = []
    for link, target in images:
        if any(marker in link.lower() for marker in CONCEPT_ART_MARKERS):
            continue  # hand-authored cartoon; doesn't track the data
        fig_mtime = target.stat().st_mtime
        # 60-second grace period to absorb "saved figure then immediately edited prose"
        if fig_mtime > report_mtime + 60:
            stale.append((link, fig_mtime - report_mtime))
    if stale:
        msgs = []
        for link, delta in stale:
            mins = int(delta // 60)
            hrs = mins // 60
            if hrs > 24:
                age = f"{hrs // 24}d ago"
            elif hrs > 0:
                age = f"{hrs}h ago"
            else:
                age = f"{mins}m ago"
            msgs.append(f"  - {link} (regenerated {age} after the report was last edited)")
        result.warnings.append(
            "Figures are newer than the report prose — legends may be stale:\n"
            + "\n".join(msgs)
        )


def check_status_field(yaml: str, result: LintResult) -> str | None:
    """Return the status value if present, None otherwise."""
    m = re.search(r"^status:\s*[\"']?([a-z-]+)[\"']?\s*$", yaml, re.MULTILINE)
    if not m:
        result.warnings.append(
            "No 'status:' field in YAML — consider adding "
            "draft | in-review | final | superseded"
        )
        return None
    status = m.group(1)
    if status not in STATUS_STYLES:
        result.errors.append(
            f"Unknown status: '{status}'. "
            f"Use one of: {', '.join(STATUS_STYLES)}"
        )
        return None
    return status


# ---------------------------------------------------------------------------
# Lint orchestrator
# ---------------------------------------------------------------------------


def lint_report(report_path: Path, numbers_path: Path | None = None) -> LintResult:
    result = LintResult()
    if not report_path.exists():
        result.errors.append(f"Report file not found: {report_path}")
        return result

    text = report_path.read_text()
    yaml, body = split_front_matter(text)
    body_no_comments = strip_html_comments(body)

    check_status_field(yaml, result)
    check_required_sections(body_no_comments, result)
    check_headline_figure(body_no_comments, result)
    check_regenerate_in_methods(body_no_comments, result)
    check_slides(body_no_comments, result)
    check_summary(body_no_comments, result)
    check_number_traceability(body_no_comments, report_path.parent, result, numbers_path)
    check_hype_words(body_no_comments, result)
    check_interactive_artifacts(body_no_comments, report_path.parent, result)
    images = check_image_links(body_no_comments, report_path.parent, result)
    check_placeholders(body, result)
    check_stale_figures(report_path, images, result)

    return result


# ---------------------------------------------------------------------------
# Status banner injection
# ---------------------------------------------------------------------------


def status_banner_html(status: str) -> str:
    bg, fg, label = STATUS_STYLES[status]
    return (
        f'<div style="background:{bg};color:{fg};padding:8px 14px;'
        f'border-radius:4px;margin:0 0 16px 0;font-weight:600;'
        f'font-family:system-ui,sans-serif;">'
        f"Status: {label}"
        f"</div>"
    )


DETAILS_CSS = """
<style>
details.report-collapse, details {
  margin: 0.4em 0 1.1em 0;
  border-left: 2px solid #d0d7de;
  padding: 0.1em 0 0.1em 0.9em;
}
details > summary {
  cursor: pointer;
  color: #57606a;
  font-size: 0.9em;
  font-weight: 600;
  list-style: none;
  user-select: none;
}
details > summary::before { content: "\\25B8  "; color: #8c959f; }
details[open] > summary::before { content: "\\25BE  "; }
details > summary:hover { color: #0969da; }
details[open] > summary { margin-bottom: 0.5em; }
</style>
"""


def inject_details_css(html_path: Path) -> None:
    """Style raw <details> blocks so progressive disclosure looks intentional."""
    html = html_path.read_text()
    # Use function replacements — DETAILS_CSS contains backslash escapes
    # (e.g. \25B8) that re would otherwise misread as group references.
    new_html, n = re.subn(r"</head>", lambda m: DETAILS_CSS + m.group(0), html, count=1)
    if n == 0:
        new_html, n = re.subn(
            r"<body[^>]*>", lambda m: m.group(0) + "\n" + DETAILS_CSS, html, count=1)
    if n == 0:
        new_html = DETAILS_CSS + html
    html_path.write_text(new_html)


def inject_status_banner(html_path: Path, status: str) -> None:
    """Inject a status banner at the top of the rendered HTML body."""
    html = html_path.read_text()
    banner = status_banner_html(status)
    # Insert after the opening <body...> tag. Both Quarto and pandoc emit one.
    # Function replacement so banner text is inserted literally.
    new_html, n = re.subn(
        r"<body[^>]*>",
        lambda m: m.group(0) + "\n" + banner,
        html,
        count=1,
    )
    if n == 0:
        # Fallback: prepend (shouldn't happen with standard renderers)
        new_html = banner + html
    html_path.write_text(new_html)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def render(report_path: Path) -> Path:
    """Render report.qmd → report.html. Returns the output path."""
    html_path = report_path.with_suffix(".html")
    if shutil.which("quarto"):
        subprocess.run(
            ["quarto", "render", str(report_path), "--to", "html"],
            check=True,
        )
    elif shutil.which("pandoc"):
        print(
            "note: quarto not found; falling back to pandoc. Install quarto via "
            "`mise install quarto` for full features.",
            file=sys.stderr,
        )
        # --embed-resources was introduced in pandoc 2.19; older pandoc (e.g. 2.10,
        # which is what ships with the anaconda install on Sherlock) spells the same
        # thing --self-contained. Pick the flag the installed pandoc actually accepts,
        # otherwise the render dies with "Unknown option --embed-resources".
        ver = subprocess.run(["pandoc", "--version"], capture_output=True,
                             text=True).stdout.split()
        embed_flag = "--embed-resources"
        try:
            major, minor = (int(x) for x in ver[1].split(".")[:2])
            if (major, minor) < (2, 19):
                embed_flag = "--self-contained"
                print(f"note: pandoc {ver[1]} predates --embed-resources; "
                      "using --self-contained.", file=sys.stderr)
        except (IndexError, ValueError):
            pass
        subprocess.run(
            [
                "pandoc",
                str(report_path),
                "-o",
                str(html_path),
                "--standalone",
                embed_flag,
                "--toc",
            ],
            check=True,
        )
    else:
        raise SystemExit(
            "ERROR: neither quarto nor pandoc found in PATH. "
            "Install quarto (`mise install quarto`) or pandoc."
        )

    if not html_path.exists():
        raise SystemExit(f"ERROR: render claimed success but {html_path} doesn't exist")
    return html_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("report", type=Path, help="Path to the .qmd report")
    parser.add_argument("--no-lint", action="store_true", help="Skip lint checks")
    parser.add_argument("--no-render", action="store_true", help="Lint only; do not render")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat warnings as errors (fail render on any warning)",
    )
    parser.add_argument(
        "--numbers",
        type=Path,
        default=None,
        help="Path to a register_value numbers.json / .manifest.json for the "
        "number-traceability check (auto-discovered near the report if omitted)",
    )
    args = parser.parse_args(argv)

    report_path = args.report.resolve()
    numbers_path = args.numbers.resolve() if args.numbers else None

    if not args.no_lint:
        result = lint_report(report_path, numbers_path)
        rc = result.report(strict=args.strict)
        if rc != 0:
            print(
                "\nLint failed. Fix the errors above, or pass --no-lint to skip.",
                file=sys.stderr,
            )
            return rc
        if not result.errors and not result.warnings and not result.soft_warnings:
            print("Lint clean.")

    if args.no_render:
        return 0

    # Read status before rendering so we know whether to inject the banner
    yaml, _ = split_front_matter(report_path.read_text())
    status_match = re.search(r"^status:\s*[\"']?([a-z-]+)[\"']?\s*$", yaml, re.MULTILINE)
    status = status_match.group(1) if status_match else None

    html_path = render(report_path)

    inject_details_css(html_path)
    if status and status in STATUS_STYLES:
        inject_status_banner(html_path, status)

    print(f"\nRendered: {html_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
