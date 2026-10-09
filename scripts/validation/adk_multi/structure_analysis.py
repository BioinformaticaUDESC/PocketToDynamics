#!/usr/bin/env python3
"""Compare structure-defined contacts and C-alpha conformations for 1AKE, 4AKE and 8RJ4."""
from pathlib import Path
import argparse,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem

def protein_atoms(path,chain='A'):
    atoms=[];ca={};seq={}
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21]!=chain or line[16] not in (' ','A'):continue
        el=line[76:78].strip().upper()
        if el in ('H','D'):continue
        key=(int(line[22:26]),line[26].strip());coord=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        atoms.append((key,line[17:20].strip(),line[12:16].strip(),coord))
        if line[12:16].strip()=='CA':ca[key]=coord;seq[key]=line[17:20].strip()
    return atoms,ca,seq

def ligand_atoms(pdb,lig_resid=None,sdf=None,site=None):
    if sdf is not None:return Chem.SDMolSupplier(str(sdf),removeHs=True)[0].GetConformer().GetPositions()
    coords=[]
    for line in Path(pdb).read_text().splitlines():
        if line.startswith('HETATM') and line[21]=='A' and int(line[22:26])==lig_resid and line[17:20]==site and line[76:78].strip().upper() not in ('H','D'):
            coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    return np.asarray(coords,float)

def contact_map(atoms,lig,cutoff=5.):
    residues={}
    for key,res,atom,xyz in atoms:residues.setdefault(key,(res,[]))[1].append(xyz)
    result={}
    for key,(res,coords) in residues.items():
        d=np.linalg.norm(np.asarray(coords)[:,None,:]-lig[None,:,:],axis=2)
        result[key[0]]={'resname':res,'min_distance_angstrom':float(d.min()),'atom_pairs_within_cutoff':int(np.count_nonzero(d<=cutoff)),'in_contact':bool(d.min()<=cutoff)}
    return result

def superpose_rmsd(mobile,reference):
    keys=sorted(set(mobile)&set(reference));X=np.array([mobile[k] for k in keys]);Y=np.array([reference[k] for k in keys]);xc=X.mean(0);yc=Y.mean(0)
    u,_,vt=np.linalg.svd((X-xc).T@(Y-yc));r=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;t=yc-xc@r;disp=np.linalg.norm(X@r+t-Y,axis=1)
    return keys,disp,float(np.sqrt(np.mean(disp**2)))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--reference-pdb',type=Path,required=True);p.add_argument('--apo-pdb',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    atoms1,ca1,seq1=protein_atoms(a.reference_pdb);atoms4,ca4,seq4=protein_atoms(a.apo_pdb);atoms8,ca8,seq8=protein_atoms(a.root/'raw/8RJ4.pdb')
    if len(ca1)!=214 or len(ca4)!=214 or len(ca8)!=214:raise SystemExit('Unexpected CA count')
    if [seq1[k] for k in sorted(seq1)] != [seq8[k] for k in sorted(seq8)]:raise SystemExit('8RJ4/1AKE chain-A sequences differ')
    ligands={'Ap5A 1AKE':ligand_atoms(None,sdf=a.reference_pdb.parent/'reference.sdf'),'ADP 301 8RJ4':ligand_atoms(a.root/'raw/8RJ4.pdb',301,site='ADP'),'ADP 302 8RJ4':ligand_atoms(a.root/'raw/8RJ4.pdb',302,site='ADP')}
    protein={'Ap5A 1AKE':atoms1,'ADP 301 8RJ4':atoms8,'ADP 302 8RJ4':atoms8};contacts={};summaries=[];cutoff=5.0
    for name,lig in ligands.items():
        cm=contact_map(protein[name],lig,cutoff);contacts[name]=cm;chosen=[(k,v) for k,v in cm.items() if v['in_contact']]
        summaries.append({'complex':name,'ligand_heavy_atoms':int(len(lig)),'protein_residues_within_5A':len(chosen),'heavy_atom_pairs_within_5A':sum(v['atom_pairs_within_cutoff'] for _,v in chosen),'residue_numbers_within_5A':','.join(str(k) for k,_ in chosen)})
    # Compare common-sequence structures after one global C-alpha Kabsch fit.
    s4,disp4,rms4=superpose_rmsd(ca4,ca1);s8,disp8,rms8=superpose_rmsd(ca8,ca1)
    struct={'n_ca_atoms':214,'rmsd_4AKE_to_1AKE_angstrom':rms4,'median_ca_displacement_4AKE_angstrom':float(np.median(disp4)),'p95_ca_displacement_4AKE_angstrom':float(np.percentile(disp4,95)),'rmsd_8RJ4_to_1AKE_angstrom':rms8,'median_ca_displacement_8RJ4_angstrom':float(np.median(disp8)),'p95_ca_displacement_8RJ4_angstrom':float(np.percentile(disp8,95)),'alignment':'Single global Kabsch fit using all 214 matched chain-A C-alpha atoms; per-residue displacements are descriptive.'}
    (a.out/'structure_summary.json').write_text(json.dumps({'contacts':summaries,'structure_comparison':struct},indent=2)+'\n')
    with (a.out/'contact_residues.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['complex','residue_number','residue_name','minimum_heavy_atom_distance_angstrom','heavy_atom_pairs_within_5A'])
        for name,cm in contacts.items():
            for pos,v in sorted(cm.items()):
                if v['in_contact']:w.writerow([name,pos,v['resname'],f"{v['min_distance_angstrom']:.3f}",v['atom_pairs_within_cutoff']])
    fig,(ax,bx)=plt.subplots(2,1,figsize=(10,6.3),gridspec_kw={'height_ratios':[2,1]},layout='constrained')
    pos4=[k[0] for k in s4];pos8=[k[0] for k in s8]
    ax.plot(pos4,disp4,color='#9a5f75',lw=1.1,label=f'4AKE apo vs 1AKE holo; global RMSD {rms4:.2f} Å')
    ax.plot(pos8,disp8,color='#366f9c',lw=1.1,label=f'8RJ4 ADP complex vs 1AKE holo; global RMSD {rms8:.2f} Å')
    for label,color in [('Ap5A 1AKE','#d36a35'),('ADP 301 8RJ4','#087e8b'),('ADP 302 8RJ4','#7b61a8')]:
        xs=[k for k,v in contacts[label].items() if v['in_contact']]
        y=np.interp(xs,pos8,disp8) if '8RJ4' in label else np.interp(xs,pos4,disp4)
        ax.scatter(xs,y,s=16,color=color,label=f'{label} contact residues',alpha=.8)
    ax.set(xlabel='E. coli adenylate kinase reference residue',ylabel='C-alpha displacement after global fit (Å)',xlim=(1,214));ax.legend(frameon=False,fontsize=8,ncol=2)
    names=list(ligands);mat=np.array([[int(contacts[n][i]['in_contact']) for i in range(1,215)] for n in names])
    bx.imshow(mat,aspect='auto',interpolation='nearest',cmap=matplotlib.colors.ListedColormap(['#f1f3f5','#d36a35']),extent=[.5,214.5,2.5,-.5],vmin=0,vmax=1)
    bx.set_yticks(range(3),names);bx.set_xlabel('E. coli adenylate kinase reference residue');bx.set_title('Residues with any protein heavy atom within 5 Å of ligand')
    fig.suptitle('Ligand-contact footprints and conformational differences in adenylate kinase',fontsize=13)
    for ext in ('png','svg'):fig.savefig(a.out/f'structure_context.{ext}',dpi=220)
    plt.close(fig);print(json.dumps({'contacts':summaries,'structure_comparison':struct},indent=2))
if __name__=='__main__':main()
