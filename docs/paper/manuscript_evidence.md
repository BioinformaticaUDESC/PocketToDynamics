# Evidence ledger for the PocketToDynamics technical manuscript

Baseline inspected: `3de57df8a2980207111a8cb71d145eea0aecb382` (8 October 2026). Paths below are relative to the repository root. This ledger distinguishes source-supported methods, recorded execution, and fresh checks; it is not a claim of fresh end-to-end computation.

| Manuscript claims | Primary repository evidence |
|---|---|
| M01–M21 implemented; 25 executable entries | `scripts/p2d_steps.py`; `scripts/run_pipeline.py`; actual module source files |
| Production 12 PASS/13 SKIP; fixture 24 PASS/1 SKIP | `manifest/STATUS.tsv`; `manifest/STATUS_fixture.tsv`; `manifest/audit_production.json`; `manifest/audit_fixture.json` |
| Recorded end-to-end fixture and M17 recovery history | `docs/AUDIT_2026-10-08.md`; `manifest/audit_fixture.json`; `logs/audit_*` |
| 4AKE identity; two 214-residue chains | `tests/fixtures/receptors/test_receptor.pdb` HEADER/TITLE/DBREF records |
| Explicit fixture, exact input identities, candidate-gate separation | `scripts/run_pipeline.py`; `scripts/10c_extract_selected_pose.py`; `README.md` |
| ConCavity score >=0.10; 8 Å C-alpha graph; minimum 3 residues; centroid | `scripts/03_map_pockets.py` |
| Docking configuration, box, seed, preparation engines | `config/screening.yaml`; `scripts/05_prepare_inputs.py`; `scripts/05c_run_screening.py` |
| Per-receptor/ligand aggregation across pocket best scores; rank weights | `scripts/07_rank_candidates.py`; `config/ranking.yaml` |
| Descriptors and gate checks | `scripts/08_chemoinformatics.py`; `scripts/09_candidate_gate.py`; `config/candidate_gate.yaml` |
| Chemical reference/pose reconstruction and ethanol atom identity | `scripts/10c_extract_selected_pose.py`; `tests/test_m10c_pose_extraction.py`; `README.md` |
| Parameterization, solvent, minimization and simulation configuration | `config/dynamics.yaml`; `scripts/11_parameterize_system.py` through `scripts/14_production_md.py` |
| 64,587 atoms, minimization, 20 ps NVT/NPT, 20 ps production | `docs/AUDIT_2026-10-08.md`; per-stage fixture logs and external fixture manifests |
| Preprocessing and integrity | `scripts/15_preprocess_trajectory.py`; `tests/test_m15_trajectory.py` |
| Structural metrics, first-frame RMSD and common fitting, contact/PBC convention | `scripts/16_analyze_trajectory.py`; `scripts/p2d_analysis.py`; `tests/test_analysis_modules.py` |
| MM/GBSA igb=5/mbondi2, no entropy, three-frame ledger | `scripts/17_binding_energetics.py`; `config/dynamics.yaml`; `tests/test_analysis_modules.py` |
| Unweighted C-alpha PCA, N−1 normalization, SVD and reconstruction | `scripts/18_essential_dynamics.py`; `tests/test_analysis_modules.py` |
| Histogram landscape formula, empty-bin infinity, eight-neighbor basins | `scripts/19_free_energy_landscape.py`; `tests/test_analysis_modules.py` |
| Contact and hydrogen-bond occupancies and heuristic definitions | `scripts/20_interaction_persistence.py`; `tests/test_analysis_modules.py`; `README.md` |
| Hash-verified archive and environment/source snapshots | `scripts/21_reproducibility_bundle.py`; `tests/test_analysis_modules.py`; `README.md` |
| Negative regression cases | `tests/test_core_regressions.py`; `tests/test_pipeline_failures.py`; `tests/test_m14_resume_integration.py` |
| Dependency specification and installation limits | `config/environment.yml`; `config/dynamics.yaml`; `README.md` |
| Sequence retrieval, representative selection, reference-guided alignment, entropy and ligand-contact mapping | `scripts/validation/evolution/`; `docs/validation/evolution/README.md`; archived UniProt TSV, FASTA, alignment, results and checksums |
| Evolutionary conservation figure and reported descriptive means | `docs/paper/figures/evolutionary_conservation.png`; `docs/validation/evolution/summary.json`; `adk_conservation.tsv` |

The sanitized `validation_summary.json` prepared alongside this manuscript summarizes existing audit evidence and fresh read-only QC without distributing heavy trajectories. `docs/PUBLICATION_AUDIT.md` records the publication checks.

## Fresh verification during publication preparation

Publication preparation included successful fresh runs of `tests/test_core_regressions.py` and `tests/test_pipeline_failures.py` using the installed Python environment, followed by seven passing read-only checks covering M15 and M16–M21 against existing artifacts. These verify selected invariants and manifests/artifact hashes. The full docking/parameterization/MD computation was not rerun. Do not describe those rechecks as new simulations or independent replicates.

## Availability

Publication preparation created the private repository `https://github.com/BioinformaticaUDESC/PocketToDynamics`. Public availability, a software license and an archival DOI were not established. The manuscript intentionally makes no assertion of a public release or public heavy-data deposit.

## Verified external reference metadata

The following primary-paper metadata was verified using online records:

- Vina: Eberhardt et al., 2021, doi:10.1021/acs.jcim.1c00203; https://pubmed.ncbi.nlm.nih.gov/34278794/
- ff19SB: Tian et al., 2020, doi:10.1021/acs.jctc.9b00591; https://pubmed.ncbi.nlm.nih.gov/31714766/
- MMPBSA.py: Miller et al., 2012, doi:10.1021/ct300418h; https://pubmed.ncbi.nlm.nih.gov/26605738/

These references identify methods; they do not independently validate this workflow. Complete third-party software/method citations remain an editorial prerequisite for journal submission.

## Deliberately excluded claims

No efficacy, affinity, enrichment, pose accuracy, experimental validation, equilibrium convergence, kinetics, superior performance, universal reproducibility, or independent replication is claimed. The configured 100 ns production duration is not reported as an executed result. Three fixture histogram basins are not labeled metastable states. The uncorrected MM/GBSA SEM is not treated as convergence evidence. Earlier M15-in-development wording is superseded by inspected implementation and recorded fixture evidence for M15–M21.

## Ap5A pilot and publication figures

A separate six-search known-site pilot generated 54 Ap5A poses. All three 1AKE top-ranked poses met the prespecified 2 Å project criterion. The aligned 4AKE comparison is descriptive and is not evaluated as a crystallographic recovery test. No new MD was performed.

- `docs/validation/adk/RESULTS.md`, `summary.csv`, `poses.csv`: all measured results.
- `docs/validation/adk/chemistry_qc.json`: eight carbon stereocenters, graph identity and no-fit/symmetry checks.
- `docs/validation/adk/protocol.json`, `inputs/preparation.json`, `run_records.json`: settings, source hashes, commands and execution history.
- `docs/validation/adk/inputs/` and `poses/`: small prepared inputs and all output poses; raw MD data remain excluded.
- `scripts/validation/adk/`: downloadable source snapshots, preparation, docking, graph-aware RMSD and figure generation.
- `docs/paper/figures/`: four PNG/SVG figures. The ethanol panel was drawn from the already verified 21-frame fixture; its reduced tables and source hashes are in `docs/validation/adk/`.

The pilot was performed separately from the M01–M21 production runner and used known-site information, a larger ligand-envelope box, isolated chain A and an assumed Ap5A charge state of −5. It does not establish pocket-discovery performance, candidate-gate calibration, binding affinity or biological efficacy.

## Added ADP validation and comparative structural analysis

- ADP redocking protocol, formal-charge assumption, receptor filtering, box construction, Vina version/settings and source SHA-256 values: `docs/validation/adk_multi/protocol.json`; `docs/validation/adk_multi/inputs/preparation.json`.
- Source coordinate/chemical-component files, prepared inputs, all 54 poses, six Vina logs, per-pose and run-level results, and hashes: `docs/validation/adk_multi/`.
- Six top-ranked site-specific pose RMSDs (0.6215–0.6370 Å), scores, and pose ranks: `summary.csv`; pose-level file `poses.csv`.
- Chain-A 1AKE/4AKE/8RJ4 global C-alpha displacement and crystallographic contact sets: `structure_summary.json`; `contact_residues.tsv`.
- ADP summary and structural-context figures are generated by `scripts/validation/adk_multi/analyze.py` and `structure_analysis.py`; native/docked coordinate projection uses `make_pose_figure.py`.
- Figure inputs are crystallographic references and independently generated ADP starting conformers; crystallographic pose coordinates are used only for evaluation and box definition. Both docking sites are in one enzyme sequence, and each search removes the crystallographic ligand, Mg2+, chloride and waters from the receptor while assuming ADP charge −3. This is known-site redocking and not independent-target validation.

## Expanded reference audit

The bibliography was expanded to cite every external program or method discussed in the manuscript, including Vina, Meeko, Open Babel, RDKit, ff19SB, GAFF, AM1-BCC, ParmEd, TIP3P, GROMACS, Amber MMPBSA.py, MDAnalysis, sequence-alignment resources, PDB structures, and published workflow comparators. Four additional papers on molecular-dynamics workflows, FAIR workflow recommendations, simulation reproducibility, and docking controls were added to the Discussion. DOI/title/journal metadata were cross-checked against publisher pages, PubMed, Zenodo, official software citation guidance, or the source project. Numeric citations were checked against the bibliography; no listed reference is uncited. The current draft has a 330-word abstract and 7,419 main-text words by the documented count, and includes nine figures and six tables.

The manuscript is still a preparation draft rather than a completed submission package: repository license/public DOI and author-supplied disclosures remain unresolved. The verified topic deadline and scope are tracked in the editorial notes, not inferred from the manuscript evidence.
