#!/usr/bin/env python3
"""Empirical PC1/PC2 occupancy landscape; discrete histogram watershed basins."""
from p2d_analysis import *

def histogram_basins(counts):
    occupied=[tuple(v) for v in np.argwhere(counts>0)]; roots={}
    def ascend(cell):
        trail=[]
        while cell not in roots:
            trail.append(cell)
            neighbors=[(a,b) for a in range(max(0,cell[0]-1),min(counts.shape[0],cell[0]+2)) for b in range(max(0,cell[1]-1),min(counts.shape[1],cell[1]+2)) if counts[a,b]>0]
            best=min(neighbors,key=lambda p:(-counts[p],p))
            if (-counts[best],best)>=(-counts[cell],cell):
                roots[cell]=cell;break
            cell=best
        root=roots[cell]
        for p in trail: roots[p]=root
        return root
    for cell in occupied: ascend(cell)
    peaks=sorted(set(roots.values()),key=lambda p:(-counts[p],p)); labels=np.full(counts.shape,-1,dtype=int)
    for cell,root in roots.items(): labels[cell]=peaks.index(root)
    return labels,peaks

def worker(row,out,cfg,args):
    prior=verified(row['manifest']); pre=preprocessing(prior); pca_path=Path(row['manifest']).parent/'pca.npz'; z=np.load(pca_path)
    scores=z['scores'][:,:2]; frames=z['frames']; temperature=float(json.loads(Path(pre['details']['production_manifest']).read_text())['temperature_k'])
    if not np.isfinite(temperature) or temperature<=0: raise ValueError('Invalid temperature')
    n=int(cfg['landscape']['fixture_bins' if pre['details']['profile']=='fixture' else 'bins'])
    if n<2: raise ValueError('At least two histogram bins required')
    counts,xe,ye=np.histogram2d(scores[:,0],scores[:,1],bins=n)
    energy=np.full_like(counts,np.inf); occupied=counts>0
    energy[occupied]=-0.008314462618*temperature*np.log(counts[occupied]/counts.max())
    labels,peaks=histogram_basins(counts)
    ix=np.clip(np.searchsorted(xe,scores[:,0],side='right')-1,0,n-1);iy=np.clip(np.searchsorted(ye,scores[:,1],side='right')-1,0,n-1)
    basin_ids=labels[ix,iy]; u=universe(pre); gr=groups(u,pre); rows=[]
    # Keep every watershed basin; no silent truncation of occupancy.
    for basin,peak in enumerate(peaks):
        members=np.where(basin_ids==basin)[0]; center=np.array([(xe[peak[0]]+xe[peak[0]+1])/2,(ye[peak[1]]+ye[peak[1]+1])/2])
        representative=members[np.argmin(np.sum((scores[members]-center)**2,axis=1))];frame=int(frames[representative]);u.trajectory[frame]
        path=out/f'basin_{basin+1:03d}.gro';gr['Protein_Ligand'].write(str(path))
        rows.append({'basin':basin+1,'frame_count':len(members),'occupancy':len(members)/len(scores),'representative_frame_0based':frame,'time_ps':float(u.trajectory.ts.time),'structure':str(path)})
    table(out/'basins.tsv',rows,list(rows[0]));table(out/'frame_basins.tsv',[{'frame_0based':int(f),'basin':int(b)+1} for f,b in zip(frames,basin_ids)],['frame_0based','basin'])
    np.savez_compressed(out/'landscape.npz',counts=counts,free_energy_kj_mol=energy,x_edges=xe,y_edges=ye,basin_labels=labels)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,5));plot=ax.pcolormesh(xe,ye,np.ma.masked_invalid(energy.T),shading='auto');fig.colorbar(plot,ax=ax,label='Relative -RT ln P (kJ/mol)');ax.set(xlabel='PC1 (nm)',ylabel='PC2 (nm)');fig.tight_layout();fig.savefig(out/'landscape.png',dpi=150);plt.close(fig)
    artifact_manifest(out,19,row['dynamics_job_id'],[row['manifest'],pca_path],
        [out/f for f in ['basins.tsv','frame_basins.tsv','landscape.npz','landscape.png']]+[Path(r['structure']) for r in rows],
        {'preprocessing_manifest':prior['details']['preprocessing_manifest'],'pca_manifest':row['manifest'],'temperature_k':temperature,'n_frames':len(scores),'n_basins':len(rows),
         'method':'Equal-width histogram -RT ln(P/Pmax); empty bins +infinity; 8-neighbor deterministic watershed of empirical density',
         'interpretation':'Sampling-dependent occupancy landscape, not validated thermodynamic free energy or kinetic metastability; no pseudocounts.'})

if __name__=='__main__': run_module(19,worker)
