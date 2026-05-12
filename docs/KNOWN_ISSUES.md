# Known issues

This document tracks issues that are known and **not yet resolved** in the
publication-ready 1.0.0-uksh release. Resolved items from the
pre-cleanup `WORKFLOW_ISSUES_REPORT.md` have been removed.

For the full pre-cleanup audit (with all 7 issues, including the 5 that
were addressed during the publication-cleanup pass), see the snapshot at
`Chapter3-Archive/2026-05-11_UKSH-TalosSV_pre-cleanup_snapshot/WORKFLOW_ISSUES_REPORT.md`
(internal workspace path, not part of the released repo).

---

## Moderate — Large reference files: directory layout not documented

The Nextflow configs (`nextflow/talos.config`, `nextflow/annotation.config`)
expect a `large_files/` directory containing reference data (Ensembl GFF3,
HPO/Phenio DB, MANE summary, ClinVar archive, etc.) but the repo does not
ship a manifest of what should live there, where to download each file
from, and what file name is expected. New users have to read the configs
line by line.

**Workaround**: see `nextflow/talos.config` parameters under
`params.large_files = ...`. Each commented-out line above a `params.xxx`
declaration names the upstream source.

**Fix idea**: add a `docs/REFERENCE_FILES.md` with a one-row-per-file table
(parameter name, expected filename, source URL, approximate size, last
verified date).

---

## Minor — `uv.lock` still references cpg-utils

`cpg-utils` is now declared as an optional `[cpg]` extra in
`pyproject.toml`, but `uv.lock` was generated when it was still a hard
runtime dependency. The lockfile therefore still pins `cpg-utils` and its
transitive dependencies.

**Impact**: none at runtime - `pip install .` without `[cpg]` still
succeeds and `cpg-utils` is simply not imported. The wasted resolution
cost is tiny.

**Fix idea**: regenerate `uv.lock` with `uv lock` after the cpg-utils
extra was moved out of `dependencies`. Deferred so that the cleanup pass
does not silently rev other transitive versions.

---

## Minor — Two upstream tests fail collection on v8.2.0

Tests under `test/_pending_v8_adaption/` were copied from
upstream `main` and target post-v8.2.0 module names (snake_case). They
need their imports rewritten before they will collect. See
`test/_pending_v8_adaption/README.md`.

This is tracked here so that future contributors know to look there
when picking up test coverage work.
