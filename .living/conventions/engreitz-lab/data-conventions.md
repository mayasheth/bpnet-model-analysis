# Data Conventions — Engreitz Lab

Mycelium's default rule is "`data/raw/` is immutable — never modify original
files." In this lab the raw bytes usually **don't live in the repo at all** —
they live on Sherlock (`$GROUP_HOME` / `$SCRATCH`) or in Google Drive. The repo
stores **pointers and manifests, not bytes.** This replaces mycelium's
byte-in-repo data model with a path-in-repo one.

## What goes where

| Thing | Where it lives | What the repo stores |
|---|---|---|
| Raw sequencing / large matrices | Sherlock `$GROUP_HOME/<lab path>` (non-sensitive) or `$SCRATCH` (sensitive) | The **path**, not the file |
| Sensitive / Stanford High-Risk data | `$SCRATCH` only — never `$GROUP_HOME`, never Git/GitHub | A pointer + a `sensitivity` note; **never the data** |
| Shared reference data / project docs | Google Drive | The Drive URL |
| Processed outputs small enough to version | `data/processed/` in the repo | The bytes (only if small and non-sensitive) |
| Every dataset | — | A row in `data/DATA_MANIFEST.md` + metadata in `data/metadata/<dataset>/` |

## data/DATA_MANIFEST.md rows point outward

Each dataset entry records the **canonical location**, not a repo path, when the
bytes live off-repo. Minimum per row:

- **name** — dataset display name (matches the Vault `dataset` note if one exists)
- **location** — one of:
  - `sherlock:$GROUP_HOME/<...>` (non-sensitive lab data)
  - `sherlock:$SCRATCH/<user>/<...>` (sensitive)
  - `drive:https://drive.google.com/...`
  - `repo:data/processed/<...>` (only for small, non-sensitive, versioned files)
- **owner** — the person who owns/produced it (link to their Vault `person` note)
- **sensitivity** — `normal` | `private` | `sensitive`
- **access** — how another lab member gets it (DTN transfer, Drive share, etc.)

## Rules

- **Never commit raw or sensitive bytes.** Large/sensitive files are gitignored;
  the manifest documents how to fetch them (DTN `rsync`, Drive link).
- **Never `chmod`/`setfacl` on Oak from Claude.** Emit the command and ask the
  user to run it (repo-level rule in `EngreitzLabAgents/CLAUDE.md`).
- **Sherlock transfers go through the DTN**
  (`engreitz@dtn.sherlock.stanford.edu`), never the login node.
- **A dataset that has a Vault `dataset` note** should link to it by
  `[[wikilink]]` in its metadata so the analysis's Vault pointer note (see
  `analysis-conventions.md`) can reference the same entity.
