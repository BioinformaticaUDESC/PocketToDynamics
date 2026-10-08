#!/usr/bin/env python3
"""Per-residue heavy-atom contact and hydrogen-bond persistence."""
from p2d_analysis import *
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import parmed

def worker(row,out,cfg,args):
    prior=verified(row['manifest']); pre=preprocessing(prior); u=universe(pre,False); gr=groups(u,pre)
    p=gr['Protein'].select_atoms('not name H*'); l=gr['Ligand'].select_atoms('not name H*')
    residues=gr['Protein'].residues; residue_map={int(r.ix):i for i,r in enumerate(residues)}
    atom_res=np.array([residue_map[int(x)] for x in p.resindices]); cutoff=float(cfg['analysis']['contact_cutoff_nm'])*10
    contacts=np.zeros((len(u.trajectory),len(residues)),dtype=bool);minimum=np.full(contacts.shape,np.inf);times=[]
    for ts in u.trajectory:
        dist=distance_array(p.positions,l.positions,box=ts.dimensions).min(axis=1)
        np.minimum.at(minimum[ts.frame],atom_res,dist/10)
        contacts[ts.frame]=minimum[ts.frame]<cutoff/10;times.append(float(ts.time))
    prod=json.loads(Path(pre['details']['production_manifest']).read_text()); topology=Path(prod['files']['topology'])
    parm=parmed.load_file(str(topology))
    if len(parm.atoms)!=len(u.atoms) or [a.name for a in parm.atoms]!=list(u.atoms.names): raise ValueError('H-bond topology identity mismatch')
    u.add_TopologyAttr('charges',[a.charge for a in parm.atoms]);u.add_TopologyAttr('masses',[a.mass for a in parm.atoms])
    u.add_TopologyAttr('bonds',[(b.atom1.idx,b.atom2.idx) for b in parm.bonds])
    ligand='resname '+' '.join(sorted(set(l.resnames))); selected='protein or ('+ligand+')'
    hb=HydrogenBondAnalysis(u,between=['protein',ligand],d_a_cutoff=3.5,d_h_a_angle_cutoff=150,update_selections=False)
    hb.hydrogens_sel=hb.guess_hydrogens(selected);hb.acceptors_sel=hb.guess_acceptors(selected)
    hb.run(); hmat=np.zeros_like(contacts); events=[]
    for frame,donor,hydrogen,acceptor,dist,angle in hb.results.hbonds:
        atoms=u.atoms[[int(donor),int(acceptor)]]
        for atom in atoms:
            if int(atom.resindex) in residue_map: hmat[int(frame),residue_map[int(atom.resindex)]]=True
        events.append({'frame_0based':int(frame),'donor_atom_0based':int(donor),'hydrogen_atom_0based':int(hydrogen),'acceptor_atom_0based':int(acceptor),'distance_angstrom':float(dist),'angle_degrees':float(angle)})
    table(out/'hbond_events.tsv',events,['frame_0based','donor_atom_0based','hydrogen_atom_0based','acceptor_atom_0based','distance_angstrom','angle_degrees'])
    rows=[]
    for i,res in enumerate(residues):
        v=contacts[:,i];runs=np.diff(np.r_[False,v,False].astype(int));starts=np.where(runs==1)[0];stops=np.where(runs==-1)[0]
        rows.append({'resindex':int(res.ix),'resid':int(res.resid),'resname':res.resname,'occupancy':float(v.mean()),'hbond_occupancy':float(hmat[:,i].mean()),
                     'contact_frames':int(v.sum()),'max_consecutive_frames':int(max(stops-starts,default=0)),'minimum_distance_nm':float(minimum[:,i].min())})
    table(out/'residue_contacts.tsv',rows,list(rows[0]));table(out/'ranked_residues.tsv',sorted(rows,key=lambda r:(-r['occupancy'],-r['hbond_occupancy'],r['resindex'])),list(rows[0]))
    np.savez_compressed(out/'contacts.npz',contact_matrix=contacts,hbond_matrix=hmat,minimum_distance_nm=minimum,times_ps=times,resindices=residues.resindices)
    artifact_manifest(out,20,row['dynamics_job_id'],[row['manifest'],topology],
        [out/f for f in ['residue_contacts.tsv','ranked_residues.tsv','contacts.npz','hbond_events.tsv']],
        {'preprocessing_manifest':prior['details']['preprocessing_manifest'],'landscape_manifest':row['manifest'],'n_frames':len(times),'contact_cutoff_nm':cutoff/10,
         'hbond_definition':'Bonded donor/hydrogen from topology; charge-based acceptor/hydrogen guess; DA <=3.5 Angstrom, DHA >150 degrees; differs from GROMACS angular convention',
         'residue_identity':'Zero-based topology resindex disambiguates chains; resid is topology numbering',
         'interpretation':'Occupancies and consecutive sampled frames; not residence times or independent observations.'})

if __name__=='__main__': run_module(20,worker)
