#!/usr/bin/env python3
"""Single trajectory MM/GBSA with AmberTools, explicit frame ledger and residue decomposition."""
from p2d_analysis import *
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors
import parmed


def worker(row,out,cfg,args):
    prior=verified(row['manifest']); pre=preprocessing(prior); d=pre['details']; ecfg=cfg['energetics']
    if int(ecfg['igb'])!=5: raise ValueError('This validated MMGBSA implementation requires igb=5 with mbondi2')
    fixture=d['profile']=='fixture'
    root=Path(ecfg['minimized_fixture_root'] if fixture else ecfg['minimized_root'])
    minpath=root/row['dynamics_job_id']/'minimization_manifest.json'
    mini=json.loads(minpath.read_text()); top=Path(mini['files']['solvated_prmtop'])
    if mini['status']!='PASS' or sha256(top)!=mini['sha256']['solvated_prmtop']: raise ValueError('Original Amber topology hash mismatch')
    prod=json.loads(Path(d['production_manifest']).read_text())
    if sha256(prod['files']['topology'])!=sha256(mini['files']['gromacs_top']): raise ValueError('Production topology differs from minimized Amber-derived topology')
    u=universe(pre,False); gr=groups(u,pre); parm=parmed.load_file(str(top))
    if len(parm.atoms)!=len(u.atoms) or [a.name for a in parm.atoms]!=list(u.atoms.names): raise ValueError('Amber/trajectory atom identity mismatch')
    amber=Path(ecfg['amber_root']); env=os.environ.copy(); env['AMBERHOME']=str(amber); env['PATH']=str(amber/'bin')+os.pathsep+env['PATH']
    ligand_resnames=sorted(set(gr['Ligand'].resnames)); solvent=sorted(set(gr['Solvent_Ions'].resnames))
    command([amber/'bin/ante-MMPBSA.py','-p',top,'-c',out/'complex.prmtop','-r',out/'receptor.prmtop','-l',out/'ligand.prmtop',
             '-s',':'+','.join(solvent),'-n',':'+','.join(ligand_resnames),'--radii','mbondi2'],out,'ante',env=env)
    dry=parmed.load_file(str(out/'complex.prmtop')); complex_atoms=gr['Protein_Ligand'].sort()
    if [a.name for a in dry.atoms]!=list(complex_atoms.names) or [a.residue.name for a in dry.atoms]!=list(complex_atoms.resnames): raise ValueError('Stripped topology ordering mismatch')
    start=int(ecfg['start_frame'])-1; stride=int(ecfg['fixture_interval' if fixture else 'interval']); limit=int(ecfg['fixture_max_frames' if fixture else 'max_frames'])
    if start<0 or stride<1 or limit<2: raise ValueError('Invalid sampling configuration')
    frames=list(range(start,len(u.trajectory),stride))[:limit]
    if len(frames)<2: raise ValueError('At least two energetic frames required')
    ledger=[]
    with mda.Writer(str(out/'sampled.nc'),n_atoms=len(complex_atoms)) as w:
        for frame in frames:
            ts=u.trajectory[frame]
            # Rigid translation only: keep ligand in the nearest periodic image to receptor.
            displacement=(gr['Ligand'].center_of_geometry()-gr['Protein'].center_of_geometry()).reshape(1,3)
            shift=minimize_vectors(displacement,ts.dimensions)[0]-displacement[0]
            gr['Ligand'].positions+=shift
            w.write(complex_atoms)
            ledger.append({'frame_0based':frame,'time_ps':float(ts.time),'ligand_shift_angstrom':','.join(map(str,shift))})
    table(out/'sampled_frames.tsv',ledger,['frame_0based','time_ps','ligand_shift_angstrom'])
    inp=out/'mmpbsa.in'
    inp.write_text(f"Audited single-trajectory MMGBSA\n&general\n startframe=1, endframe={len(frames)}, interval=1, keep_files=0, verbose=1,\n/\n&gb\n igb={int(ecfg['igb'])}, saltcon={float(ecfg['saltcon'])},\n/\n&decomp\n idecomp=2, dec_verbose=0,\n/\n")
    command([amber/'bin/MMPBSA.py','-O','-i',inp,'-cp',out/'complex.prmtop','-rp',out/'receptor.prmtop','-lp',out/'ligand.prmtop',
             '-y',out/'sampled.nc','-o',out/'FINAL_RESULTS_MMPBSA.dat','-do',out/'FINAL_DECOMP_MMPBSA.dat',
             '-eo',out/'energies.csv','-deo',out/'decomposition.csv'],out,'mmpbsa',env=env)
    result=(out/'FINAL_RESULTS_MMPBSA.dat').read_text()
    matches=re.findall(r'^DELTA TOTAL\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)',result,re.M)
    if len(matches)!=1 or not np.isfinite(np.array(matches,dtype=float)).all(): raise ValueError('Missing/nonfinite MMGBSA binding estimate')
    decomposition=(out/'FINAL_DECOMP_MMPBSA.dat').read_text()
    if len(decomposition)<100 or re.search(r'\b(?:nan|inf)\b',decomposition,re.I): raise ValueError('Invalid residue decomposition')
    files=[out/f for f in ['complex.prmtop','receptor.prmtop','ligand.prmtop','sampled.nc','sampled_frames.tsv','mmpbsa.in','FINAL_RESULTS_MMPBSA.dat','FINAL_DECOMP_MMPBSA.dat','energies.csv','decomposition.csv']]
    artifact_manifest(out,17,row['dynamics_job_id'],[row['manifest'],minpath,top,mini['files']['gromacs_top']],files,
        {'preprocessing_manifest':prior['details']['preprocessing_manifest'],'analysis_manifest':row['manifest'],'method':'MM/GBSA igb=5 mbondi2; single trajectory; no entropy',
         'n_frames':len(frames),'delta_total_kcal_mol':float(matches[0][0]),'sd_kcal_mol':float(matches[0][1]),'naive_sem_kcal_mol':float(matches[0][2]),
         'sampling':ledger,'interpretation':'Endpoint estimate only. Correlated-frame SEM is not convergence evidence; fixture is software QC only.'})

if __name__=='__main__': run_module(17,worker)
