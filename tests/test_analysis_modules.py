#!/usr/bin/env python3
"""Numerical fixture QC, run after each module before advancing."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from p2d_analysis import *
m=int(sys.argv[1]); out=HEAVY/(DIRS[m]+'_fixture')/'MDTEST_a5b7f2f48479'
d=verified(out/(NAMES[m]+'_manifest.json'))['details']
if m==16:
    rows=read_table(out/'metrics.tsv'); a=np.array([[float(v) for v in r.values()] for r in rows])
    assert len(rows)==21 and np.isfinite(a).all() and (a>=0).all()
    assert float(rows[0]['protein_rmsd_nm'])<1e-6 and float(rows[0]['ligand_rmsd_nm'])<1e-6
    assert all(float(r['protein_sasa_nm2'])>0 for r in rows)
    assert all(float(r['protein_rg_nm'])>0 for r in rows)
    assert len(read_table(out/'rmsf.tsv'))==428
elif m==17:
    assert d['n_frames']==3 and np.isfinite(d['delta_total_kcal_mol'])
    assert [r['frame_0based'] for r in d['sampling']]==[0,10,20]
    assert (out/'decomposition.csv').stat().st_size>100
elif m==18:
    z=np.load(out/'pca.npz'); v=z['eigenvalues']; modes=z['modes']; score=z['scores']
    assert len(score)==21 and np.all(v>=0) and np.all(np.diff(v)<=1e-10)
    assert np.allclose(modes@modes.T,np.eye(len(modes)),atol=1e-6)
    assert np.allclose(score.mean(axis=0),0,atol=1e-7)
    assert np.isclose(np.var(score,axis=0,ddof=1).sum(),v.sum())
elif m==19:
    z=np.load(out/'landscape.npz'); count=z['counts']; energy=z['free_energy_kj_mol']; mask=count>0
    assert count.sum()==21 and np.isclose(energy[mask].min(),0)
    assert np.isinf(energy[~mask]).all()
    rows=read_table(out/'basins.tsv'); assert rows
    assert sum(int(r['frame_count']) for r in rows)==21
    assert all(Path(r['structure']).is_file() for r in rows)
elif m==20:
    rows=read_table(out/'residue_contacts.tsv'); assert len(rows)==428
    assert all(0<=float(r['occupancy'])<=1 for r in rows)
    z=np.load(out/'contacts.npz'); assert z['contact_matrix'].shape==(21,428)
    assert np.allclose(z['contact_matrix'].mean(axis=0),[float(r['occupancy']) for r in rows])
elif m==21:
    import tarfile
    with tarfile.open(out/'reproducibility.tar.gz') as t:
        names=t.getnames(); assert 'bundle/inventory.json' in names and 'bundle/report.md' in names
print(f'PASS_M{m}_NUMERICAL_FIXTURE_QC')
