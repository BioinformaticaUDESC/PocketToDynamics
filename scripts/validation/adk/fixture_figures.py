"""Render descriptive figures from hash-verified existing ethanol fixture data."""
from pathlib import Path
import argparse,csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--heavy-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
job='MDTEST_a5b7f2f48479';evidence={}
for folder,name in [('analyzed_fixture','analysis'),('pca_fixture','pca'),('landscape_fixture','landscape')]:
 file=a.heavy_root/folder/job/(name+'_manifest.json');d=json.loads(file.read_text());assert d['status']=='PASS'
 for label in ['artifact_hashes','input_hashes']:
  for fn,digest in d.get(label,{}).items():assert hashlib.sha256(Path(fn).read_bytes()).hexdigest()==digest
 evidence[name]={'manifest_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'artifact_hashes':{Path(k).name:v for k,v in d['artifact_hashes'].items()}}
f=a.heavy_root/'analyzed_fixture'/job/'metrics.tsv';rows=list(csv.DictReader(f.open(),delimiter='\t'));v={k:np.array([float(r[k]) for r in rows]) for k in rows[0]};assert len(rows)==21
pc=np.load(a.heavy_root/'pca_fixture'/job/'pca.npz');ls=np.load(a.heavy_root/'landscape_fixture'/job/'landscape.npz');assert ls['counts'].sum()==21
(out/'fixture_metrics.tsv').write_bytes(f.read_bytes())
np.savetxt(out/'fixture_pca_projection.csv',np.column_stack([pc['times_ps'],pc['scores'][:,:2]]),delimiter=',',header='time_ps,PC1_nm,PC2_nm',comments='')
landrows=[]
for i in range(5):
 for j in range(5):landrows.append({'pc1_low_nm':ls['x_edges'][i],'pc1_high_nm':ls['x_edges'][i+1],'pc2_low_nm':ls['y_edges'][j],'pc2_high_nm':ls['y_edges'][j+1],'count':ls['counts'][i,j],'relative_energy_kj_mol':ls['free_energy_kj_mol'][i,j]})
with (out/'fixture_landscape.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=landrows[0]);w.writeheader();w.writerows(landrows)
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(10,8),layout='constrained');t=v['time_ps']
ax=axs[0,0];ax.plot(t,v['protein_rmsd_nm'],color='#007c91',label='Protein');ax.plot(t,v['ligand_rmsd_nm'],color='#b74816',label='Ethanol');ax.set(xlabel='Time (ps)',ylabel='RMSD (nm)',title='A  RMSD after common backbone fitting');ax.legend()
ax=axs[0,1];ax.plot(t,v['protein_rg_nm'],color='#007c91');ax.set(xlabel='Time (ps)',ylabel='Radius of gyration (nm)',title='B  Protein compactness');ax.ticklabel_format(useOffset=False)
ax=axs[1,0];sc=ax.scatter(pc['scores'][:,0],pc['scores'][:,1],c=pc['times_ps'],cmap='viridis',s=36);ax.set(xlabel='PC1 (nm)',ylabel='PC2 (nm)',title=f"C  PCA projection ({100*pc['variance_fraction'][:2].sum():.1f}% variance)");fig.colorbar(sc,ax=ax,label='Time (ps)')
ax=axs[1,1];z=np.ma.masked_invalid(ls['free_energy_kj_mol'].T);im=ax.pcolormesh(ls['x_edges'],ls['y_edges'],z,cmap='magma',shading='flat');ax.set(xlabel='PC1 (nm)',ylabel='PC2 (nm)',title='D  Empirical occupancy landscape');fig.colorbar(im,ax=ax,label='Relative −RT ln(P/Pmax) (kJ/mol)')
fig.suptitle('Ethanol software fixture: 21 frames over 20 ps; no convergence claim',fontsize=13)
for ext in ['png','svg']:fig.savefig(out/f'fixture_analysis.{ext}',dpi=220)
summary={'n_frames':21,'duration_ps':20,'protein_rmsd_range_nm':[float(v['protein_rmsd_nm'].min()),float(v['protein_rmsd_nm'].max())],'ligand_rmsd_range_nm':[float(v['ligand_rmsd_nm'].min()),float(v['ligand_rmsd_nm'].max())],'protein_rg_range_nm':[float(v['protein_rg_nm'].min()),float(v['protein_rg_nm'].max())],'pc1_pc2_variance_fraction':float(pc['variance_fraction'][:2].sum()),'occupied_landscape_bins':int(np.sum(ls['counts']>0)),'evidence':evidence}
(out/'fixture_figure_evidence.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='evidence'},indent=2))
