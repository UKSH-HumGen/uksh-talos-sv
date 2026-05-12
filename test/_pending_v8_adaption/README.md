# Pending tests (need v8.2.0-uksh API adaption)

These test files were copied verbatim from
[populationgenomics/talos@main](https://github.com/populationgenomics/talos)
(snapshot used during the cleanup). They target a **post-v8.2.0** API
where the module file names have been switched to `snake_case`
(`run_hail_filtering.py`, `download_panelapp.py`, etc.).

This fork is pinned to upstream **v8.2.0**, where those same modules
still use `CamelCase` (`RunHailFiltering.py`, `DownloadPanelApp.py`,
...). Importing the tests as-is therefore fails at collection time
with `ModuleNotFoundError`.

The two upstream tests whose API still matches v8.2.0
(`test_utils.py`, `test_moi_tests.py`) live in the parent
`test/` directory and are part of the active suite.

## What to do here

For each file in this directory, before re-enabling it:

1. Rewrite the imports to the v8.2.0 module names, e.g.
   `from talos.run_hail_filtering import ...`
   becomes
   `from talos.RunHailFiltering import ...`.
2. Re-check the called function/class names against the v8.2.0 source.
   Some helpers were renamed between v8.2 and current `main`.
3. Run `pytest -q test/the_test.py` and fix any remaining drift.
4. Move the file back up to `test/` and add a line to `CHANGELOG.md`.

Pull-requests that restore a file are very welcome.
