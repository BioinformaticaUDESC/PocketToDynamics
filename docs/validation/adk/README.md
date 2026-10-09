# Adenylate kinase Ap5A pilot

This is a known-site pose-recovery pilot, separate from the M01–M21 ethanol integration fixture. It is not an enrichment benchmark, affinity measurement or validation of automated pocket discovery. The ligand is the experimentally co-crystallized inhibitor Ap5A (CCD AP5) in 1AKE. The apo structure 4AKE is a conformational comparison, not a second crystallographic ligand reference.

## Prespecified design

- One ligand; chain A of 1AKE and 4AKE; blank/A alternative conformers; all receptor heteroatoms removed.
- Chemical graph from the RCSB CCD. One terminal phosphate hydroxyl proton is removed per phosphorus, explicitly assuming formal charge −5. Protonation/metal-state uncertainty remains untested.
- Crystallographic heavy-atom coordinates are reserved as the reference. Docking starts from an independent RDKit ETKDGv3 conformer (seed 20261008), relaxed with UFF. No crystal ligand self-fitting is allowed in RMSD evaluation.
- One global C-alpha Kabsch alignment maps 4AKE into the 1AKE frame using 214 matched chain-A residues. The reported 4AKE ligand RMSD is therefore descriptive and alignment-dependent.
- Search box: crystallographic ligand bounding box plus 6 Å per side, with each dimension at least 20 Å; same box in the aligned receptor frame for both receptors. This uses prior site knowledge.
- AutoDock Vina 1.2.5; exhaustiveness 8; up to 9 modes; energy range 4 kcal/mol; seeds 20261007, 20261008, 20261009; two CPU threads per search.
- Six searches run concurrently on a 16-logical-CPU machine. The first two-worker scheduling attempt was interrupted before completion and retained locally under `attempts`; no completed result was selected or discarded to improve recovery.
- Project pose-recovery criterion: symmetry-aware 57-heavy-atom RMSD ≤2 Å in the fixed receptor frame. Report top-ranked recovery and best-returned recovery separately for 1AKE. Do not apply a success/failure threshold to the apo comparison.
- All results, including failures to recover the reference, are retained. No MD is launched from these outputs and no candidate-gate thresholds are changed.

## Reproduction

Use the `pocket2dynamics` environment described in `config/environment.yml`, with Vina installed at `/usr/bin/vina`. The recorded run used RDKit 2026.09.1 and Meeko 0.8.0. A fresh output directory is required; existing calculations are not overwritten.

```bash
python scripts/validation/adk/fetch.py --root work/adk_pilot
python scripts/validation/adk/prepare.py --root work/adk_pilot
python scripts/validation/adk/dock.py --root work/adk_pilot
python scripts/validation/adk/analyze.py --root work/adk_pilot
```

The source downloader enforces the SHA-256 hashes of the RCSB snapshots used here. If a public structure changes, review the revision explicitly instead of silently accepting new inputs. Preparation refuses to change inputs after a docking directory exists. Docking refuses to overwrite an attempt. Analysis verifies output hashes and chemical graphs before RMSD calculation and checks that a rigidly translated ligand retains its 10 Å displacement (no hidden self-fit).

The figure for the older ethanol fixture is generated separately:

```bash
python scripts/validation/adk/fixture_figures.py \
  --heavy-root /path/to/PocketToDynamics/dynamics \
  --output work/adk_pilot/report
```

This command reads existing manifests and checks artifact/input hashes. It does not generate new MD data. Tables and figure-source data are published in `docs/validation/adk`; larger raw execution directories remain in local `work/`.

## Primary structural sources

- [1AKE](https://www.rcsb.org/structure/1AKE): Müller and Schulz, 1992, DOI [10.1016/0022-2836(92)90582-5](https://doi.org/10.1016/0022-2836(92)90582-5).
- [4AKE](https://www.rcsb.org/structure/4AKE): Müller et al., 1996, DOI [10.1016/S0969-2126(96)00018-4](https://doi.org/10.1016/S0969-2126(96)00018-4).
- [AP5 chemical dictionary](https://www.rcsb.org/ligand/AP5).

Limitations include a single large flexible ligand, one assumed charge state, rigid receptors, excluded waters/ions, one initial conformer and three search seeds at fixed exhaustiveness. Search seeds are algorithmic repeats, not biological or MD replicates. Vina scores must not be interpreted as experimental binding free energies.

Published text tables and SVG figures use normalized line endings/whitespace. Source-artifact hashes in the fixture evidence refer to the original verified files; numerical values are unchanged.
