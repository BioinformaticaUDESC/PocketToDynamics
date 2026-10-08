#!/usr/bin/env python3
"""Structural metrics on common fitted coordinates; periodic contacts on unrotated box."""
from p2d_analysis import *
from MDAnalysis.lib.distances import distance_array

def worker(row,out,cfg,args):
    pre=verified(row['manifest']); d=pre['details']; u=universe(pre); gr=groups(u,pre)
    u.trajectory[0]; prot=gr['Protein']; lig=gr['Ligand']; ca=gr['C_alpha']
    p0=prot.positions.copy(); l0=lig.positions.copy(); mean=np.zeros((len(ca),3)); m2=mean.copy()
    records=[]
    for i,ts in enumerate(u.trajectory):
        delta=ca.positions-mean; mean+=delta/(i+1); m2+=delta*(ca.positions-mean)
        records.append({'time_ps':float(ts.time),'protein_rmsd_nm':float(np.sqrt(np.mean(np.sum((prot.positions-p0)**2,axis=1)))/10),
                        'ligand_rmsd_nm':float(np.sqrt(np.mean(np.sum((lig.positions-l0)**2,axis=1)))/10),
                        'protein_rg_nm':float(prot.radius_of_gyration()/10)})
    table(out/'rmsf.tsv',[{'resindex':int(r),'resid':int(rid),'resname':rn,'rmsf_nm':float(v)} for r,rid,rn,v in zip(ca.resindices,ca.resids,ca.resnames,np.sqrt(m2.sum(axis=1)/len(records))/10)],['resindex','resid','resname','rmsf_nm'])
    centered=universe(pre,False); cg=groups(centered,pre)
    p=cg['Protein'].select_atoms('not name H*'); l=cg['Ligand'].select_atoms('not name H*')
    cutoff=float(cfg['analysis']['contact_cutoff_nm'])*10
    for i,ts in enumerate(centered.trajectory):
        distances=distance_array(p.positions,l.positions,box=ts.dimensions)
        records[i]['minimum_heavy_distance_nm']=float(distances.min()/10)
        records[i]['heavy_atom_contact_pairs']=int((distances<cutoff).sum())
    gmx=cfg['trajectory_preprocessing']['gmx_executable']
    base=['-s',d['reference_tpr'],'-f',d['centered_xtc'],'-n',d['index']]
    command([gmx,'sasa',*base,'-surface','group "Protein"','-output','group "Protein"','-o',out/'sasa.xvg'],out,'sasa')
    command([gmx,'hbond-legacy',*base,'-num',out/'hbonds.xvg','-r','0.35','-a','30','-nthreads','2'],out,'hbonds','1\n2\n')
    for name,file in [('protein_sasa_nm2','sasa.xvg'),('protein_ligand_hbonds','hbonds.xvg')]:
        arr=xvg(out/file)
        if len(arr)!=len(records) or not np.allclose(arr[:,0],[r['time_ps'] for r in records],atol=1e-4): raise ValueError('Metric frame alignment mismatch')
        for r,v in zip(records,arr[:,1]): r[name]=float(v)
    if not np.isfinite([[r[k] for k in records[0]] for r in records]).all(): raise ValueError('Nonfinite metrics')
    table(out/'metrics.tsv',records,list(records[0]))
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(2,2,figsize=(10,7))
    for a,k in zip(ax.flat,['protein_rmsd_nm','ligand_rmsd_nm','protein_rg_nm','protein_sasa_nm2']):
        a.plot([r['time_ps'] for r in records],[r[k] for r in records]); a.set(xlabel='Time (ps)',ylabel=k)
    fig.tight_layout(); fig.savefig(out/'metrics.png',dpi=150); plt.close(fig)
    artifact_manifest(out,16,row['dynamics_job_id'],[row['manifest']],
        [out/f for f in ['metrics.tsv','rmsf.tsv','sasa.xvg','hbonds.xvg','metrics.png']],
        {'preprocessing_manifest':row['manifest'],'n_frames':len(records),'rmsd_reference':'First fitted frame; common backbone fit; no ligand self-fit',
         'rmsf_reference':'Mean C-alpha position in common backbone fit','contact_cutoff_nm':cutoff/10,'hbond_cutoffs':{'donor_acceptor_nm':0.35,'deviation_from_linear_degrees':30},
         'summary':{k:float(np.mean([r[k] for r in records])) for k in records[0] if k!='time_ps'}})

if __name__=='__main__': run_module(16,worker)
