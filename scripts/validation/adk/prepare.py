"""Prepare chain-A receptors and an independent Ap5A conformer for a known-site pilot."""
from pathlib import Path
import json,hashlib,subprocess,sys
import gemmi,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);args=parser.parse_args();root=args.root.resolve();root.mkdir(parents=True,exist_ok=True);
raw=root/'raw';out=root/'prepared'
if (root/'docking').exists():
 raise SystemExit('Refusing to change prepared inputs after docking. Use a fresh --root.')
out.mkdir(exist_ok=True)
def pdb_atoms(name,record,resname=None):
 lines=[]
 for l in (raw/name).read_text().splitlines():
  if l.startswith(record) and l[21]=='A' and l[16] in (' ','A') and (resname is None or l[17:20]==resname):lines.append(l[:16]+' '+l[17:])
 return lines
rec1=pdb_atoms('1AKE.pdb','ATOM  ');rec4=pdb_atoms('4AKE.pdb','ATOM  ')
lig=pdb_atoms('1AKE.pdb','HETATM','AP5');xyz={l[12:16].strip():[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lig}
b=gemmi.cif.read(str(raw/'AP5.cif')).sole_block(); atoms=list(b.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol','_chem_comp_atom.charge']));rw=Chem.RWMol();idx={}
for row in atoms:
 name=gemmi.cif.as_string(row[0]);a=Chem.Atom(row[1]);a.SetFormalCharge(int(row[2]));a.SetProp('ccd_name',name);idx[name]=rw.AddAtom(a)
for row in b.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']):
 rw.AddBond(idx[gemmi.cif.as_string(row[0])],idx[gemmi.cif.as_string(row[1])],{'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE}[row[2]])
m=rw.GetMol();Chem.SanitizeMol(m)
# Explicit assumed phosphate state: remove one terminal OH proton per phosphorus.
remove=[];deprot=[]
for a in m.GetAtoms():
 if a.GetSymbol()=='O' and any(n.GetSymbol()=='P' for n in a.GetNeighbors()):
  hs=[n for n in a.GetNeighbors() if n.GetSymbol()=='H']
  if hs:
   assert len(hs)==1;remove.append(hs[0].GetIdx());deprot.append(a.GetProp('ccd_name'));a.SetFormalCharge(-1)
rw=Chem.RWMol(m)
for i in sorted(remove,reverse=True):rw.RemoveAtom(i)
m=Chem.RemoveHs(rw.GetMol());Chem.SanitizeMol(m);assert len(remove)==5
assert set(xyz)=={a.GetProp('ccd_name') for a in m.GetAtoms()}
c=Chem.Conformer(m.GetNumAtoms())
for a in m.GetAtoms():c.SetAtomPosition(a.GetIdx(),xyz[a.GetProp('ccd_name')])
m.AddConformer(c);Chem.AssignAtomChiralTagsFromStructure(m,replaceExistingTags=True)
for a in m.GetAtoms():
 if a.GetSymbol()=='P':a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
Chem.AssignStereochemistry(m,cleanIt=True,force=True)
m.SetProp('_Name','AP5_reference_chain_A_altloc_A_charge_minus5')
w=Chem.SDWriter(str(out/'reference.sdf'));w.write(m);w.close()
h=Chem.AddHs(m);h.RemoveAllConformers();p=AllChem.ETKDGv3();p.randomSeed=20261008;p.numThreads=1
assert AllChem.EmbedMolecule(h,p)==0
uff=AllChem.UFFOptimizeMolecule(h,maxIters=2000)
h.SetProp('_Name','AP5_independent_ETKDG_charge_minus5');w=Chem.SDWriter(str(out/'ligand.sdf'));w.write(h);w.close()
def ca(lines):return {int(l[22:26]):np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in lines if l[12:16].strip()=='CA'}
a=ca(rec1);d=ca(rec4);keys=sorted(set(a)&set(d));assert len(keys)==214
X=np.array([d[k] for k in keys]);Y=np.array([a[k] for k in keys]);xc=X.mean(0);yc=Y.mean(0);u,_,vt=np.linalg.svd((X-xc).T@(Y-yc));R=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;t=yc-xc@R
aligned=[]
for l in rec4:
 q=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])@R+t;aligned.append(l[:30]+''.join(f'{v:8.3f}' for v in q)+l[54:])
for name,lines in [('1AKE',rec1),('4AKE',aligned)]: (out/f'{name}.pdb').write_text('\n'.join(lines)+'\nTER\nEND\n')
coords=np.array(list(xyz.values()));center=(coords.min(0)+coords.max(0))/2;size=np.maximum(coords.max(0)-coords.min(0)+12,20)
metadata={'ligand':'Ap5A / CCD AP5','heavy_atoms':m.GetNumAtoms(),'formula':rdMolDescriptors.CalcMolFormula(h),'formal_charge':Chem.GetFormalCharge(h),'deprotonated_atom_names':deprot,'charge_state_note':'Assumed minus-five phosphate state; no pH titration or metal/protonation ensemble validation','altloc':'A or blank','chain':'A','water_policy':'All HETATM removed from receptors; no waters/ions','ligand_start':'Independent ETKDGv3 conformation; UFF relaxation, not crystal conformation','embedding_seed':20261008,'uff_returncode':uff,'stereo_smiles':Chem.MolToSmiles(m),'box_center_angstrom':center.tolist(),'box_size_angstrom':size.tolist(),'alignment':'4AKE to 1AKE, all 214 matched chain-A C-alpha atoms, rigid Kabsch; no ligand fit','alignment_ca_rmsd_angstrom':float(np.sqrt(np.mean(np.sum((X@R+t-Y)**2,axis=1)))),'rotation':R.tolist(),'translation':t.tolist(),'raw_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in raw.iterdir()},'sources':{'1AKE.pdb':'https://files.rcsb.org/download/1AKE.pdb','4AKE.pdb':'https://files.rcsb.org/download/4AKE.pdb','AP5.cif':'https://files.rcsb.org/ligands/download/AP5.cif'}}
(out/'preparation.json').write_text(json.dumps(metadata,indent=2)+'\n');print(json.dumps(metadata,indent=2))
for name in ['1AKE','4AKE']:
 cmd=[sys.executable,str(Path(sys.executable).parent/'mk_prepare_receptor.py'),'--read_pdb',str(out/f'{name}.pdb'),'-o',str(out/name),'-p',str(out/f'{name}.pdbqt'),'--compute_charges']
 r=subprocess.run(cmd,capture_output=True,text=True);(out/f'{name}_prep.log').write_text(r.stdout+r.stderr);print(name,r.returncode);assert r.returncode==0,r.stderr
cmd=[sys.executable,str(Path(sys.executable).parent/'mk_prepare_ligand.py'),'-i',str(out/'ligand.sdf'),'-o',str(out/'ligand.pdbqt')];r=subprocess.run(cmd,capture_output=True,text=True);(out/'ligand_prep.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
print('PREPARATION_PASS')
