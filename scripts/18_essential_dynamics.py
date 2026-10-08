#!/usr/bin/env python3
"""Unweighted Cartesian C-alpha PCA on the common M15 backbone fit."""
from p2d_analysis import *

def worker(row,out,cfg,args):
    prior=verified(row['manifest']); pre=preprocessing(prior); u=universe(pre); ca=groups(u,pre)['C_alpha']
    stride=int(cfg['pca']['stride']); limit=int(cfg['pca']['max_frames'])
    if stride<1 or limit<3: raise ValueError('Invalid PCA sampling')
    frames=np.arange(0,len(u.trajectory),stride)
    if len(frames)>limit: raise ValueError('PCA exceeds max_frames; explicitly increase stride or memory budget')
    coordinates=[]; times=[]
    for i in frames:
        ts=u.trajectory[int(i)]; coordinates.append(ca.positions.copy().ravel()/10); times.append(float(ts.time))
    x=np.asarray(coordinates,dtype=np.float64)
    if len(x)<3 or not np.isfinite(x).all(): raise ValueError('Invalid PCA coordinates')
    mean=x.mean(axis=0); centered=x-mean
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=2): left,s,modes=np.linalg.svd(centered,full_matrices=False)
    # Remove the centering null mode; retain the full nonzero sample covariance spectrum.
    keep=s>max(s[0]*1e-10,1e-12); s=s[keep]; modes=modes[keep]; scores=centered@modes.T
    if len(s)<2: raise ValueError('PCA rank below two')
    eigenvalues=s*s/(len(x)-1); fraction=eigenvalues/eigenvalues.sum()
    if not np.allclose(scores@modes,centered,atol=1e-7): raise ValueError('PCA reconstruction error')
    np.savez_compressed(out/'pca.npz',mean=mean,modes=modes,scores=scores,eigenvalues=eigenvalues,variance_fraction=fraction,
                        frames=frames,times_ps=times,atom_indices=ca.indices)
    table(out/'eigenvalues.tsv',[{'component':i+1,'eigenvalue_nm2':v,'variance_fraction':f} for i,(v,f) in enumerate(zip(eigenvalues,fraction))],['component','eigenvalue_nm2','variance_fraction'])
    table(out/'projection.tsv',[{'frame_0based':int(frame),'time_ps':t,'PC1_nm':float(v[0]),'PC2_nm':float(v[1])} for frame,t,v in zip(frames,times,scores)],['frame_0based','time_ps','PC1_nm','PC2_nm'])
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,5)); plot=ax.scatter(scores[:,0],scores[:,1],c=times,cmap='viridis'); fig.colorbar(plot,ax=ax,label='Time (ps)');ax.set(xlabel='PC1 (nm)',ylabel='PC2 (nm)'); fig.tight_layout();fig.savefig(out/'pca.png',dpi=150);plt.close(fig)
    artifact_manifest(out,18,row['dynamics_job_id'],[row['manifest'],prior['details']['preprocessing_manifest']],
        [out/f for f in ['pca.npz','eigenvalues.tsv','projection.tsv','pca.png']],
        {'preprocessing_manifest':prior['details']['preprocessing_manifest'],'energetics_manifest':row['manifest'],'n_frames':len(x),'n_atoms':len(ca),'rank':len(s),
         'method':'Unweighted Cartesian C-alpha PCA; common M15 backbone fit; covariance normalization N-1; coordinates nm',
         'variance_PC1_PC2':float(fraction[:2].sum()),'stride':stride})

if __name__=='__main__': run_module(18,worker)
