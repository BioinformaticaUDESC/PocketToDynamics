# PocketToDynamics

Local molecular screening and dynamics pipeline. Source and lightweight state live in
`~/PocketToDynamics`; heavy outputs live in `/mnt/labio-fast-disk/PocketToDynamics/dynamics`.
This repository produces software/QC artifacts, not a manuscript.

## Run

Activate `pocket2dynamics` (Python 3.11), then:

```bash
cd ~/PocketToDynamics
bash scripts/run_pipeline.sh --profile production
```

Without activation, set `P2D_PYTHON` to that environment's Python executable. The runner
adds its bin directory to PATH. Dependencies are in `config/environment.yml`;
GROMACS and AmberTools paths are in `config/dynamics.yaml`.
The validated stack uses GROMACS 2026.3, AmberTools 25/MMPBSA.py 14.0,
MDAnalysis 2.10.0, NumPy, SciPy, ParmEd and Matplotlib. MDAnalysis reads GRO/XTC,
because its TPR parser does not yet support this GROMACS TPR format.

For a complete computational fixture audit:

```bash
bash scripts/run_pipeline.sh --profile fixture
```

This requires exactly the committed `test_receptor.pdb` and `test_ligand.sdf` in
`input/receptors` and `input/ligands`, with matching hashes. It does not replace user
inputs. M01–M10B exercise the regular pipeline. At M10C, the explicitly selected
fixture supplies its recorded ethanol pose despite rejection by M09, so that
M11–M21 can be tested. Candidate thresholds are never relaxed. The fixture rebuilds
parameterization, minimization, equilibration and 20 ps production, then validates
each module before proceeding. Earlier fixture production is archived under
`dynamics/audit_history`; M11–M13 preserve previous attempts as `.previous-*`.
Expect this audit to require several minutes of real computation.

## Module registry

| Modules | Output |
|---|---|
| M01–M04 | Input hashes, normalized structures, pockets, job matrix |
| M05A/B/C | Engine/preparation QC and docking with fingerprinted resume |
| M06–M09 | Parsed scores, ranking, descriptors, finite-value candidate gate |
| M10A/B/C | Dynamics jobs, dependency preflight, selected pose |
| M11–M14 | AMBER/GAFF2 parameters, solvent/minimization, equilibration, segmented MD |
| M15 | Whole/nojump/centered/fitted trajectories; frame, atom and hash QC |
| M16 | RMSD, RMSF, Rg, SASA, H-bonds, distances and contacts |
| M17 | MM/GBSA (igb=5, mbondi2), selected frames and residue decomposition |
| M18 | Unweighted Cartesian C-alpha PCA with reconstruction QC |
| M19 | PC1/PC2 occupancy landscape and watershed representative structures |
| M20 | Per-residue contact/H-bond occupancy and ranked residues |
| M21 | Verified archive, source/environment snapshots, figures, tables, software report |

Each M15–M21 script supports `--profile`, `--manifest`, `--output-root`, and `--log`.
The default profile is production; fixture paths remain separate. Example:

```bash
python scripts/15_preprocess_trajectory.py --profile fixture
python tests/test_m15_trajectory.py
python scripts/16_analyze_trajectory.py --profile fixture
python tests/test_analysis_modules.py 16
```

## Integrity and status

`manifest/STATUS.tsv` reports all 25 production module entries; `STATUS_fixture.tsv`
reports the fixture branch. Missing logs are UNKNOWN, not omitted. The runner marks
unexecuted modules NOT_RUN on each invocation, stops on nonzero exits or QC failures,
and records commands, exit codes, duration and log hashes in `manifest/audit_*.json`.
Per-run stdout/stderr are in `logs/audit_*`. SKIP means no eligible upstream jobs and
must not be interpreted as a completed scientific calculation.

M05C resume requires identical input files, box, docking parameters, task,
executable contents/version and output pose hash. Previous attempts are retained.
M14 validates source coordinates/checkpoint/topology, configuration and executable
signature before checkpoint reuse; completed segment outputs are hashed. Legacy
checkpoints without provenance are rejected: choose a new output root. The explicit
M14 `--overwrite` option is destructive and should only be used intentionally.
M15 checks segment hashes and ordered intervals, all frame times/coordinates/boxes,
full-system atom counts and trajectory coverage. Downstream modules verify artifacts
before reading them. M15–M21 recompute rather than silently resuming incomplete work.

Negative regression checks:

```bash
python tests/test_core_regressions.py
python tests/test_pipeline_failures.py
```

## Interpretation and units

The 20 ps adenylate-kinase/ethanol fixture validates software only. No convergence,
binding affinity or biological effect is established. Production simulation duration
and sampling are explicit configuration choices that still need scientific assessment.

Structural tables use ps and nm; SASA uses nm². Protein and ligand RMSD use the first
frame in the common backbone fit, with no ligand self-fit. RMSF uses the mean C-alpha
position. Contacts use heavy atoms and minimum-image distances on the unrotated
centered trajectory. M16 H-bonds use GROMACS's angular convention; M20 uses a bonded
donor, charge-based acceptor/hydrogen heuristic and DHA >150° (its convention is
recorded separately). These counts are not assumed identical. Zero contacts or
H-bonds are valid observations.

M17 reports a single-trajectory endpoint estimate in kcal/mol, without entropy.
The frame ledger is explicit; its ordinary SEM does not correct serial correlation.
M18 uses common fitted coordinates in nm and covariance normalization N−1. M19 uses
`-RT ln(P/Pmax)` in kJ/mol, with no pseudocounts and infinite energy for unsampled bins.
Discrete density basins are sampling-dependent and are not validated kinetic states.
M20 persistence means fraction of sampled frames, not a residence-time estimate.

## Bundle

M21 writes `reproducibility.tar.gz` under the appropriate heavy `bundles` directory.
The archive contains repository sources, inputs, intermediate dynamics artifacts,
commands and environment versions, with an inventory mapping original paths to
relative archive paths. Every entry is read back and hash-checked after compression.
After extraction, run `python bundle/verify_bundle.py` from the extracted root.
Source paths intentionally retain the documented home/heavy-disk layout; a move to
another machine requires explicit path configuration and dependency installation.
External binaries and conda environments are described, not embedded. The outer
M21 manifest records archive completion without a circular self-hash.

Method references: [GROMACS analysis documentation](https://manual.gromacs.org/current/reference-manual/analysis.html),
[GROMACS H-bond definitions](https://manual.gromacs.org/current/reference-manual/analysis/hydrogen-bonds.html),
[MDAnalysis documentation](https://docs.mdanalysis.org/2.10.0/), and the installed
AmberTools `MMPBSA.py --input-file-help` / `ante-MMPBSA.py --help`.
