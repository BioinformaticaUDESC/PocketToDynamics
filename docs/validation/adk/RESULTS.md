# Ap5A validation results

Six Vina searches completed successfully, yielding 54 poses. **1AKE top-ranked pose recovery passed in 3/3 searches** under the prespecified project criterion of symmetry-aware heavy-atom RMSD ≤2 Å without ligand fitting. This is a one-ligand known-site pilot, not a screening-enrichment or affinity benchmark.

| Receptor | Seed | Top score (kcal/mol) | Top RMSD (Å) | Best RMSD (Å) | Best-pose rank |
|---|---|---|---|---|---|
| 1AKE | 20261007 | -13.441 | 1.764 | 1.457 | 2 |
| 1AKE | 20261008 | -13.615 | 1.482 | 1.482 | 1 |
| 1AKE | 20261009 | -13.789 | 1.762 | 1.720 | 2 |
| 4AKE | 20261007 | -7.658 | 9.626 | 6.324 | 9 |
| 4AKE | 20261008 | -7.252 | 8.950 | 5.505 | 3 |
| 4AKE | 20261009 | -8.517 | 5.947 | 5.947 | 1 |

The 4AKE comparison uses an apo receptor globally aligned to 1AKE (214 C-alpha atoms; protein RMSD 7.131 Å). Its ligand deviations are descriptive; 4AKE has no experimental Ap5A reference pose. The structures were modeled as isolated chain A. No MD was performed for Ap5A.

Inputs used an explicitly assumed −5 phosphate charge state, eight verified carbon stereocenters, and an independent ETKDGv3/UFF starting conformer. Three seeds (20261007–20261009), exhaustiveness 8, nine modes, a 4 kcal/mol energy window, and the same known-site box were used for both receptors. No parameters were tuned after examining completed results. Terminal phosphate resonance symmetry was included in the no-fit metric. All poses and scores are retained.

![All searches](../../paper/figures/docking_summary.png)

![Pose projections for the prespecified first seed](../../paper/figures/pose_comparison.png)

The additional fixture figure describes the existing 21 frames over 20 ps of the ethanol integration test, not a new simulation. The PCA projection accounts for 47.57% of sampled variance; ten of 25 landscape bins were occupied. Neither statistic demonstrates thermodynamic convergence.

![Existing ethanol fixture](../../paper/figures/fixture_analysis.png)

The manuscript and repository contain methods, limitations, source hashes, commands, tables and reproducible analysis scripts. Further validation would require additional ligands/targets, reference controls, protonation and receptor-ensemble sensitivity, and independently designed MD sampling.
