"""Verify docking artifacts, compute no-fit symmetry RMSD, and plot all six searches."""
from pathlib import Path
import json,csv,re,copy,hashlib
import numpy as np
from rdkit import Chem
from rmsd import matching_graph,fixed_frame_rmsd,self_check
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);args=parser.parse_args();root=args.root.resolve();root.mkdir(parents=True,exist_ok=True);
out=root/'report';out.mkdir(exist_ok=True)
protocol=json.loads((root/'protocol.json').read_text())
for name,digest in protocol['input_sha256'].items():
 assert hashlib.sha256((root/'prepared'/name).read_bytes()).hexdigest()==digest
expected={f'{rec}_{seed}' for rec in protocol['receptors'] for seed in protocol['seeds']}
assert {p.name for p in (root/'docking').iterdir()}==expected
ref=Chem.SDMolSupplier(str(root/'prepared/reference.sdf'),removeHs=True)[0];assert ref is not None
self_check(ref)
rows=[];mols={};summaries=[]
for dest in sorted((root/'docking').iterdir()):
 run=json.loads((dest/'run.json').read_text());assert run['returncode']==run['export_returncode']==0
 for name,digest in run['sha256'].items():assert hashlib.sha256((dest/name).read_bytes()).hexdigest()==digest
 energies=[float(x) for x in re.findall(r'REMARK VINA RESULT:\s+([-\d.]+)',(dest/'poses.pdbqt').read_text())]
 poses=list(Chem.SDMolSupplier(str(dest/'poses.sdf'),removeHs=True));assert 0<len(poses)==len(energies)<=9
 assert np.isfinite(energies).all() and all(a<=b for a,b in zip(energies,energies[1:]))
 local=[]
 for rank,(m,e) in enumerate(zip(poses,energies),1):
  assert m is not None and m.GetNumAtoms()==57
  assert np.isfinite(m.GetConformer().GetPositions()).all()
  # Reject graph/stereo mismatch; symmetry-aware metric never aligns the ligand.
  assert Chem.MolToSmiles(matching_graph(m))==Chem.MolToSmiles(matching_graph(ref)),(dest,rank)
  rms=fixed_frame_rmsd(m,ref)
  row={'receptor':run['receptor'],'seed':run['seed'],'rank':rank,'vina_score_kcal_mol':e,'symmetry_rmsd_angstrom':rms};rows.append(row);local.append(row);mols[(run['receptor'],run['seed'],rank)]=m
 best=min(local,key=lambda x:x['symmetry_rmsd_angstrom'])
 summaries.append({'receptor':run['receptor'],'seed':run['seed'],'n_poses':len(local),'top_score_kcal_mol':local[0]['vina_score_kcal_mol'],'top_rmsd_angstrom':local[0]['symmetry_rmsd_angstrom'],'best_rmsd_angstrom':best['symmetry_rmsd_angstrom'],'best_rmsd_rank':best['rank'],'top_pose_recovered':local[0]['symmetry_rmsd_angstrom']<=2 if run['receptor']=='1AKE' else None,'any_pose_recovered':best['symmetry_rmsd_angstrom']<=2 if run['receptor']=='1AKE' else None,'runtime_seconds':run['runtime_seconds']})
for name,records in [('poses.csv',rows),('summary.csv',summaries)]:
 with (out/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
(out/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(3,1,figsize=(8.5,8),layout='constrained');colors={'1AKE':'#007c91','4AKE':'#b74816'}
for rec in ['1AKE','4AKE']:
 v=[s for s in summaries if s['receptor']==rec];x=np.arange(3)+(-.08 if rec=='1AKE' else .08)
 for ax,key in zip(axes,['top_score_kcal_mol','top_rmsd_angstrom','best_rmsd_angstrom']):ax.plot(x,[s[key] for s in v],'o-',label=rec,color=colors[rec],lw=1.2);ax.set_xticks(range(3),['20261007','20261008','20261009'],rotation=30);ax.set_xlabel('Search seed')
axes[0].set_ylabel('Vina score (kcal/mol)');axes[0].set_title('A  Top-ranked score');axes[1].set_ylabel('Heavy-atom RMSD (Å)');axes[1].set_title('B  Top-ranked pose');axes[2].set_ylabel('Heavy-atom RMSD (Å)');axes[2].set_title('C  Best returned pose')
for ax in axes[1:]:ax.axhline(2,color='#777777',linestyle='--',lw=1,label='2 Å criterion (1AKE only)')
axes[0].legend();axes[2].legend(fontsize=8)
fig.suptitle('Ap5A pilot: closed-state redocking and open-state comparison',fontsize=13)
for ext in ['png','svg']:fig.savefig(out/f'docking_summary.{ext}',dpi=220)
plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
for i,rec in enumerate(['1AKE','4AKE']):
 m=mols[(rec,20261007,1)];pts=m.GetConformer().GetPositions();rp=ref.GetConformer().GetPositions();pr=[]
 for l in (root/f'prepared/{rec}.pdb').read_text().splitlines():
  if l.startswith('ATOM') and l[12:16].strip()=='CA':pr.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 pr=np.array(pr)
 for j,(a,b) in enumerate([(0,1),(0,2)]):
  ax=axes[i,j];ax.plot(pr[:,a],pr[:,b],color='#cccccc',lw=.6,label='C-alpha trace',zorder=0)
  for mol,p,color,label in [(ref,rp,'#222222','Crystal reference'),(m,pts,colors[rec],'Top-ranked docking pose')]:
   for bond in mol.GetBonds():z=p[[bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]];ax.plot(z[:,a],z[:,b],color=color,lw=1.3)
   ax.scatter(p[:,a],p[:,b],s=9,c=color,label=label)
  ax.set_aspect('equal');ax.set_xlabel('XYZ'[a]+' (Å)');ax.set_ylabel('XYZ'[b]+' (Å)');ax.set_title(f'{rec} | seed 20261007 | '+('XY' if j==0 else 'XZ'))
axes[0,0].legend(fontsize=8)
fig.suptitle('Pose projections in the fixed receptor frame (no ligand superposition)',fontsize=12)
for ext in ['png','svg']:fig.savefig(out/f'pose_comparison.{ext}',dpi=220)
plt.close(fig)
print(json.dumps(summaries,indent=2));print('ANALYSIS_PASS')
