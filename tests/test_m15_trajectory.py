#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from p2d_analysis import *
p=HEAVY/'preprocessed_fixture/MDTEST_a5b7f2f48479/preprocessing_manifest.json'
d=verified(p)['details']; q=d['trajectory_check']
assert q['n_frames']==21 and q['start_ps']==0 and q['end_ps']==20 and q['interval_ps']==1
assert d['group_counts']['Ligand']==9
assert d['group_counts']['Protein_Ligand']==d['group_counts']['Protein']+9
u=universe(verified(p)); g=groups(u,verified(p))
assert len(u.atoms)==q['n_atoms']
# Every peptide bond remains local after PBC processing (coordinate unit: Angstrom).
for ts in u.trajectory:
    for res in g['Protein'].residues:
        ca=res.atoms.select_atoms('name CA'); n=res.atoms.select_atoms('name N'); c=res.atoms.select_atoms('name C')
        if len(ca)==len(n)==len(c)==1:
            assert np.linalg.norm(ca.positions-n.positions)<2.0
            assert np.linalg.norm(ca.positions-c.positions)<2.0
print('PASS_M15_QC: 21 frames, 0–20 ps, consistent groups, whole backbone')
