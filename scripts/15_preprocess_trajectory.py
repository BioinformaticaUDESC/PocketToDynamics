#!/usr/bin/env python3
"""Canonical trajectory with atom identity, temporal continuity and hashes QC."""
from p2d_analysis import *
import MDAnalysis as mda

GROUPS=['System','Protein','Ligand','Protein_Ligand','Backbone','C_alpha','Solvent_Ions']

def trajectory_qc(tpr, xtc, expected_end):
    u=mda.Universe(str(tpr),str(xtc)); times=[]
    for ts in u.trajectory:
        if not np.isfinite(ts.positions).all(): raise ValueError('Nonfinite coordinates')
        if ts.dimensions is None or not np.isfinite(ts.dimensions).all() or (ts.dimensions[:3]<=0).any(): raise ValueError('Invalid periodic box')
        times.append(float(ts.time))
    delta=np.diff(times)
    if len(times)<3 or (delta<=0).any(): raise ValueError('Insufficient or nonmonotonic frames')
    if not np.allclose(delta, np.median(delta),rtol=1e-4,atol=1e-4): raise ValueError('Trajectory gap/irregular sampling')
    if not np.isclose(times[0],0,atol=1e-3) or not np.isclose(times[-1],expected_end,atol=1e-3): raise ValueError('Trajectory does not cover production interval')
    return {'n_frames':len(times),'start_ps':times[0],'end_ps':times[-1],'interval_ps':float(np.median(delta)),'n_atoms':len(u.atoms)}

def worker(row,out,cfg,args):
    pcfg=cfg['trajectory_preprocessing']; gmx=pcfg['gmx_executable']
    if pcfg['output_group']!='System': raise ValueError('Full System output required to preserve topology indexing')
    production_path=Path(row['final_gro']).parent.parent/'production_manifest.json'
    production=json.loads(production_path.read_text()); base=production_path.parent
    if production['status']!='PASS': raise ValueError('Production not PASS')
    inputs=[production_path,ROOT/'config/dynamics.yaml']; xtcs=[]; last=0.0
    for s in production['segments']:
        if not np.isclose(s['start_ps'],last): raise ValueError('Segment ordering/continuity mismatch')
        last=s['end_ps']; directory=base/s['segment']
        marker=directory/'segment_complete.json'; meta=json.loads(marker.read_text())
        if meta['status']!='PASS': raise ValueError('Incomplete segment')
        for ext,h in meta['sha256'].items():
            if sha256(directory/('md.'+ext))!=h: raise ValueError(f'Segment hash mismatch: {ext}')
        xtcs.append(directory/'md.xtc'); inputs.extend([marker,directory/'md.xtc',directory/'md.tpr'])
    if not np.isclose(last,production['target_duration_ps']): raise ValueError('Incomplete production duration')
    tpr=base/production['segments'][0]['segment']/'md.tpr'
    gro=tpr.with_suffix('.gro')
    u=mda.Universe(str(gro))
    protein=u.select_atoms('protein'); ligand=u.select_atoms('resname '+' '.join(pcfg['ligand_resnames']))
    if not len(protein) or not len(ligand): raise ValueError(f'Empty protein/ligand; observed {sorted(set(u.atoms.resnames))}')
    if np.intersect1d(protein.indices,ligand.indices).size: raise ValueError('Overlapping groups')
    gr={'System':u.atoms,'Protein':protein,'Ligand':ligand,'Protein_Ligand':protein+ligand,
        'Backbone':protein.select_atoms('name N CA C'),'C_alpha':protein.select_atoms('name CA'),
        'Solvent_Ions':u.atoms[np.setdiff1d(u.atoms.indices,np.r_[protein.indices,ligand.indices])]}
    if len(gr['C_alpha'])<3: raise ValueError('Insufficient C-alpha atoms')
    ndx=out/'analysis.ndx'
    with ndx.open('w') as f:
        for name in GROUPS:
            f.write(f'[ {name} ]\n')
            ids=gr[name].indices+1
            for start in range(0,len(ids),15): f.write(' '.join(map(str,ids[start:start+15]))+'\n')
    raw=out/'production_raw.xtc'
    command([gmx,'trjcat','-f',*xtcs,'-o',raw],out,'trjcat')
    qc=trajectory_qc(gro,raw,last)
    current=raw
    # Make each molecule whole before nojump; then image whole molecules about complex.
    steps=[('whole',['-pbc','whole'],'0\n')]
    if pcfg.get('pbc_nojump',True): steps.append(('nojump',['-pbc','nojump'],'0\n'))
    steps.extend([('centered',['-center','-pbc','mol','-ur','compact'],str(GROUPS.index(pcfg['center_group']))+'\n0\n'),
                  ('fit',['-fit','rot+trans'],str(GROUPS.index(pcfg['fit_group']))+'\n0\n')])
    outputs=[ndx,raw]
    for name,opts,stdin in steps:
        dest=out/f'production_{name}.xtc'
        command([gmx,'trjconv','-s',tpr,'-f',current,'-n',ndx,'-o',dest,*opts],out,name,stdin)
        current=dest; outputs.append(dest)
    finalqc=trajectory_qc(gro,current,last)
    if finalqc!=qc: raise ValueError('Frame count/time/atom count changed during processing')
    reference=out/'reference_fit.gro'
    command([gmx,'trjconv','-s',tpr,'-f',current,'-n',ndx,'-o',reference,'-dump','0'],out,'reference','0\n')
    check=command([gmx,'check','-f',current],out,'gmx_check')
    (out/'gmx_check.txt').write_text(check)
    details={'preprocessing_manifest':str(out/'preprocessing_manifest.json'),'reference_tpr':str(tpr),
             'production_manifest':str(production_path),'fit_xtc':str(current),'centered_xtc':str(out/'production_centered.xtc'),
             'reference_fit':str(reference),'index':str(ndx),'trajectory_check':qc,
             'indices':{k:v.indices.tolist() for k,v in gr.items()},'group_counts':{k:len(v) for k,v in gr.items()},
             'profile':production['profile'],'processing':pcfg,'gromacs_version':command([gmx,'--version'],out,'version')}
    artifact_manifest(out,15,row['dynamics_job_id'],inputs,outputs+[reference,tpr,out/'gmx_check.txt'],details)

if __name__=='__main__': run_module(15,worker)
