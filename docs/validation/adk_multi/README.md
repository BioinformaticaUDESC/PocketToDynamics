# ADP known-site redocking validation

This supplemental dataset records six AutoDock Vina 1.2.5 runs in two crystallographic ADP sites (residues 301 and 302) of chain A in PDB 8RJ4. Each ligand input was generated from an independent ETKDGv3 conformer; crystal coordinates were used only as the evaluation reference and to define the known-site box. The docking receptor contains chain-A protein atoms only. Crystallographic ADP, Mg2+, chloride and water were omitted. The modeled ADP charge state was assumed to be −3; these choices limit the interpretation.

Each site box used the crystallographic ligand heavy-atom bounding box plus 6 Å on each side (minimum dimension 20 Å). Vina settings: exhaustiveness 16, nine output modes, energy range 4 kcal/mol, two CPUs, seeds 20261010–20261012. RMSD was symmetry-aware, no-fit heavy-atom RMSD to the crystallographic ADP in the corresponding site. The 2 Å threshold was specified before interpreting results.

The coordinate model and ligand dictionary were retrieved from RCSB PDB (8RJ4, resolution 2.11 Å; ADP CCD). SHA-256 hashes for all files are in `SHA256SUMS`. `inputs/` contains source and prepared inputs; `poses/` contains all 54 output poses and run logs; CSV/JSON files contain pose-level, run-level and structural summaries. `scripts/validation/adk_multi/` contains preparation, docking, analysis and contact-analysis scripts.

Results: top-ranked pose recovery was 6/6; RMSD range 0.622–0.637 Å; mean 0.628 Å (sample SD 0.006 Å). These are repeated searches of two known sites in a single protein, not six independent targets or a prospective docking benchmark. No affinity or biological activity is inferred.

The PDB coordinates and CCD chemical component data are attributed to the wwPDB/RCSB PDB; consult their current data-use terms. The repository license status remains unresolved.
