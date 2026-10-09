"""No-fit ligand RMSD with carbon stereochemistry and phosphate resonance symmetry."""
import numpy as np
from rdkit import Chem

def matching_graph(mol):
    mol=Chem.RemoveHs(Chem.Mol(mol))
    # P=O / P-O(-) terminal oxygens are resonance-equivalent. Use a charge-
    # separated single-bond representation ONLY for atom matching, not docking.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()!=15:
            continue
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if atom.HasProp('_CIPCode'):
            atom.ClearProp('_CIPCode')
        for bond in atom.GetBonds():
            other=bond.GetOtherAtom(atom)
            if other.GetAtomicNum()==8 and other.GetDegree()==1 and bond.GetBondType()==Chem.BondType.DOUBLE:
                bond.SetBondType(Chem.BondType.SINGLE)
                other.SetFormalCharge(-1)
                atom.SetFormalCharge(atom.GetFormalCharge()+1)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    return mol

def fixed_frame_rmsd(probe, reference):
    p=matching_graph(probe);r=matching_graph(reference)
    if Chem.MolToSmiles(p)!=Chem.MolToSmiles(r):
        raise ValueError('Ligand graph or carbon stereochemistry mismatch')
    matches=r.GetSubstructMatches(p,uniquify=False,useChirality=True,maxMatches=100000)
    if not matches or len(matches)>=100000:
        raise ValueError('Missing or truncated symmetry mapping')
    x=probe.GetConformer().GetPositions();y=reference.GetConformer().GetPositions()
    return min(float(np.sqrt(np.mean(np.sum((x-y[list(match)])**2,axis=1)))) for match in matches)

def self_check(reference):
    assert fixed_frame_rmsd(reference,reference)<1e-8
    shifted=Chem.Mol(reference)
    for i,pos in enumerate(reference.GetConformer().GetPositions()):
        shifted.GetConformer().SetAtomPosition(i,pos+np.array([10.,0.,0.]))
    assert abs(fixed_frame_rmsd(shifted,reference)-10.)<1e-6
    swapped=Chem.Mol(reference)
    for atom in reference.GetAtoms():
        if atom.GetAtomicNum()!=15:continue
        oxy=[n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum()==8 and n.GetDegree()==1]
        assert len(oxy)==2
        x=reference.GetConformer().GetPositions()
        swapped.GetConformer().SetAtomPosition(oxy[0],x[oxy[1]])
        swapped.GetConformer().SetAtomPosition(oxy[1],x[oxy[0]])
    assert fixed_frame_rmsd(swapped,reference)<1e-8

    inverted=Chem.Mol(reference)
    for atom in inverted.GetAtoms():
        if atom.GetAtomicNum()==6 and atom.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
            atom.InvertChirality()
            break
    try:
        fixed_frame_rmsd(inverted,reference)
    except ValueError:
        pass
    else:
        raise AssertionError('A different carbon stereoisomer must be rejected')
