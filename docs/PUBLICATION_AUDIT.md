# GitHub publication audit — 2026-10-08

## Scope and destination

Reviewed local `main` at `3de57df` before publication. The working tree was clean and no remote existed. Created the private repository https://github.com/BioinformaticaUDESC/PocketToDynamics and configured `origin`, after the owner selected this account. Existing commits are preserved; no history rewriting or force push is used.

## Content review

- Enumerated tracked paths and all reachable historical objects: 114 file blobs at the initial revision.
- Largest historical blob: 309,339 bytes, a receptor PDB fixture. No historical MD trajectory/checkpoint files, environment directories, private keys or credential-file paths were found by the filename checks.
- Scanned all initial historical blobs for private-key headers, GitHub token formats, AWS access-key identifiers and quoted credential assignments. No matches were found. Pattern-based review cannot guarantee detection of every possible secret.
- Versioned content consists of source, tests, lightweight receptor/ligand fixtures, configuration and documentation. Existing machine paths are configuration references; environment contents and heavy-disk artifacts are excluded. Ordinary existing author metadata remains in the preserved Git history.
- Strengthened `.gitignore` for runtime state, inputs, heavy MD formats, archives, environments, caches, local configuration and common credential files. `config/environment.yml` remains versioned as a dependency specification.
- Historical ignored runtime state was not assumed safe merely because it is ignored now: reachable historical blobs were reviewed separately.

## Scientific evidence and checks

The earlier conversation snapshot described M15 as under development and M16–M21 as planned. The audited checkout contains implemented M01–M21 and the recorded audit in `docs/AUDIT_2026-10-08.md`. The manuscript therefore follows the checkout, not the older snapshot.

Existing production audit: 12 PASS, 13 SKIP; no candidates promoted. Existing fixture audit: 24 PASS, 1 SKIP, with an explicit fixture branch at M10C. The portable `docs/paper/validation_summary.json` records module states and hashes of the local audit records without copying local command paths or runtime outputs.

Fresh checks during this publication task:

- `tests/test_core_regressions.py`: PASS.
- `tests/test_pipeline_failures.py`: PASS (invalid candidate rejection messages are expected negative-test output).
- `tests/test_m15_trajectory.py`: PASS against the existing trajectory.
- `tests/test_analysis_modules.py 16` through `21`: all PASS against existing artifacts, including manifest integrity verification.

These checks were run using the existing `pocket2dynamics` Python environment. No full MD rerun was performed for manuscript preparation. The 20 ps fixture demonstrates software execution and numerical consistency only; it establishes neither convergence nor biological efficacy or binding affinity.

## Release limitations

The repository is private, not a public software release. No software license or DOI has been assigned by this work. Cross-machine paths and external dependencies require configuration. The evidence bundle and trajectories remain outside GitHub. A future submission requires confirmed author/affiliation information, journal formatting, licensing and an appropriate data-access arrangement.
