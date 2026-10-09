#!/usr/bin/env python3
"""Prepare ADP reference and independent conformer for 8RJ4 chain-A redocking."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys
import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve();raw=root/'raw';out=root/'prepared'
    if (root/'docking').exists():raise SystemExit('Refusing to alter prepared inputs after docking; choose a fresh root')
    out.mkdir(parents=True,exist_ok=True)
    block=gemmi.cif.read(str(raw/'ADP.cif')).sole_block();rw=Chem.RWMol();ids={}
    for row in block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol','_chem_comp_atom.charge']):
        name=gemmi.cif.as_string(row[0]);atom=Chem.Atom(gemmi.cif.as_string(row[1]));atom.SetFormalCharge(int(gemmi.cif.as_string(row[2])));atom.SetProp('ccd_name',name);ids[name]=rw.AddAtom(atom)
    for row in block.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']):
        order=gemmi.cif.as_string(row[2]);rw.AddBond(ids[gemmi.cif.as_string(row[0])],ids[gemmi.cif.as_string(row[1])],{'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE}[order])
    m=rw.GetMol();Chem.SanitizeMol(m)
    # ADP is prepared in an explicitly assumed -3 phosphate state; remove one
    # alpha-phosphate and two beta-phosphate hydroxyl protons from the CCD graph.
    deprotonated={'HOA2':'O2A','HOB2':'O2B','HOB3':'O3B'}
    for hname,oname in deprotonated.items():
        h=next(x for x in m.GetAtoms() if x.GetProp('ccd_name')==hname);o=next(x for x in m.GetAtoms() if x.GetProp('ccd_name')==oname);o.SetFormalCharge(-1)
    rw=Chem.RWMol(m)
    for name in deprotonated:rw.RemoveAtom(next(x.GetIdx() for x in rw.GetAtoms() if x.GetProp('ccd_name')==name))
    m=Chem.RemoveHs(rw.GetMol());Chem.SanitizeMol(m);assert Chem.GetFormalCharge(m)==-3
    pdb=(raw/'8RJ4.pdb').read_text().splitlines();protein=[l for l in pdb if l.startswith('ATOM  ') and l[21]=='A' and l[16] in (' ','A')]
    if not protein:raise SystemExit('No chain-A protein atoms in 8RJ4')
    (out/'receptor.pdb').write_text('\n'.join(protein+['TER','END'])+'\n')
    site_meta={}
    for resid in (301,302):
        lines=[l for l in pdb if l.startswith('HETATM') and l[21]=='A' and int(l[22:26])==resid and l[17:20]=='ADP' and l[16] in (' ','A')]
        xyz={l[12:16].strip():[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lines}
        names={x.GetProp('ccd_name') for x in m.GetAtoms()}
        if set(xyz)!=names:raise SystemExit(f'ADP {resid} atom-name mismatch; missing {names-set(xyz)}, extra {set(xyz)-names}')
        ref=Chem.Mol(m);conf=Chem.Conformer(ref.GetNumAtoms())
        for atom in ref.GetAtoms():conf.SetAtomPosition(atom.GetIdx(),xyz[atom.GetProp('ccd_name')])
        ref.AddConformer(conf);Chem.AssignAtomChiralTagsFromStructure(ref,replaceExistingTags=True);Chem.AssignStereochemistry(ref,cleanIt=True,force=True)
        ref.SetProp('_Name',f'8RJ4_chain_A_ADP_{resid}_CCD_ADP_charge_minus3')
        Chem.SDWriter(str(out/f'reference_{resid}.sdf')).write(ref)
        start=Chem.AddHs(ref);start.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=20261010;params.numThreads=1
        if AllChem.EmbedMolecule(start,params)!=0:raise SystemExit(f'ETKDG failed for ADP {resid}')
        uff=AllChem.UFFOptimizeMolecule(start,maxIters=2000);start.SetProp('_Name',f'ADP_independent_ETKDG_charge_minus3_seed20261010')
        w=Chem.SDWriter(str(out/f'ligand_{resid}.sdf'));w.write(start);w.close()
        coords=np.array(list(xyz.values()));center=(coords.min(0)+coords.max(0))/2;size=np.maximum(coords.max(0)-coords.min(0)+12,20)
        site_meta[str(resid)]={'center_angstrom':center.tolist(),'size_angstrom':size.tolist(),'reference_heavy_atoms':ref.GetNumAtoms(),'uff_returncode':uff}
    cmd=[sys.executable,str(Path(sys.executable).parent/'mk_prepare_receptor.py'),'--read_pdb',str(out/'receptor.pdb'),'-o',str(out/'receptor'),'-p',str(out/'receptor.pdbqt'),'--compute_charges']
    result=subprocess.run(cmd,capture_output=True,text=True);(out/'receptor_prep.log').write_text(result.stdout+result.stderr)
    if result.returncode:raise SystemExit(result.stderr)
    for resid in (301,302):
        cmd=[sys.executable,str(Path(sys.executable).parent/'mk_prepare_ligand.py'),'-i',str(out/f'ligand_{resid}.sdf'),'-o',str(out/f'ligand_{resid}.pdbqt')]
        result=subprocess.run(cmd,capture_output=True,text=True);(out/f'ligand_{resid}_prep.log').write_text(result.stdout+result.stderr)
        if result.returncode:raise SystemExit(result.stderr)
    def sequence_from_pdb(path,chain):
        lines=path.read_text().splitlines();seen={}
        for line in lines:
            if line.startswith('ATOM  ') and line[21]==chain and line[12:16].strip()=='CA' and line[16] in (' ','A'):
                seen[(int(line[22:26]),line[26])]=gemmi.find_tabulated_residue(line[17:20].strip()).one_letter_code
        return ''.join(seen[k] for k in sorted(seen)),len(seen)
    seq8,n8=sequence_from_pdb(raw/'8RJ4.pdb','A')
    repo=Path(__file__).resolve().parents[3]
    seq1,n1=sequence_from_pdb(repo/'docs/validation/adk/inputs/1AKE.pdb','A')
    if seq8!=seq1:raise SystemExit('8RJ4 chain-A sequence differs from the 1AKE reference; do not treat as an identical construct')
    meta={'pdb_entry':'8RJ4','resolution_angstrom':2.11,'chain':'A','chain_sequence_length':n8,'sequence_comparison':'Chain A sequence is identical to chain A in 1AKE across all 214 C-alpha positions; RCSB marks the expressed construct engineered.','ligand':'ADP, CCD component ADP','ligand_sites':[301,302],'ligand_formal_charge':-3,'charge_assumption':'One alpha- and two beta-phosphate hydroxyls deprotonated; Mg2+ is omitted from the docking receptor; charge/metal-state effects are not tested.','starting_conformer':'Independent ETKDGv3 conformer; UFF relaxation; seed 20261010','receptor_policy':'Chain A protein ATOM records only; all HETATM including ADP, Mg2+, chloride and waters omitted.','box_policy':'Each site uses its crystallographic ADP heavy-atom bounding box plus 6 A per side, with a 20 A minimum dimension.','sites':site_meta,'sources':{'8RJ4.pdb':'https://files.rcsb.org/download/8RJ4.pdb','ADP.cif':'https://files.rcsb.org/ligands/download/ADP.cif'},'source_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in raw.iterdir()}}
    (out/'preparation.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
if __name__=='__main__':main()
