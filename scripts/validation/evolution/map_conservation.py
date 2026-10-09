#!/usr/bin/env python3
"""Project MSA entropy and coverage onto a protein structure and ligand site."""
from pathlib import Path
import argparse, csv, json, math
from Bio import Align, SeqIO
from Bio.PDB import PDBParser, Polypeptide
from rdkit import Chem

AA=set('ACDEFGHIKLMNPQRSTVWY')

def seq_from_chain(chain):
    residues=[]; seq=[]
    for r in chain:
        if r.id[0] != ' ': continue
        aa=Polypeptide.protein_letters_3to1_extended.get(r.resname.upper(), 'X')
        if aa=='X' or 'CA' not in r: continue
        residues.append(r); seq.append(aa)
    return ''.join(seq),residues

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--alignment',type=Path,required=True);p.add_argument('--structure',type=Path,required=True);p.add_argument('--ligand',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--chain',default='A');p.add_argument('--reference',default='P69441');p.add_argument('--contact-cutoff',type=float,default=5.0)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    records=list(SeqIO.parse(a.alignment,'fasta'));ref=[r for r in records if r.id.split('|')[0]==a.reference]
    if len(ref)!=1:raise SystemExit(f'Expected exactly one reference {a.reference}; found {len(ref)}')
    alen=len(ref[0].seq)
    if any(len(r.seq)!=alen for r in records):raise SystemExit('Alignment rows differ in length')
    chain=PDBParser(QUIET=True).get_structure('adk',str(a.structure))[0][a.chain];pdbseq,pdbres=seq_from_chain(chain);query=str(ref[0].seq).upper().replace('-','')
    aligner=Align.PairwiseAligner();aligner.mode='global';aligner.match_score=2;aligner.mismatch_score=-1;aligner.open_gap_score=-8;aligner.extend_gap_score=-0.5
    aln=aligner.align(query,pdbseq)[0];q_to_p={}
    for (q0,q1),(s0,s1) in zip(*aln.aligned):
        for qpos,spos in zip(range(q0,q1),range(s0,s1)):q_to_p[qpos]=spos
    if len(q_to_p)<0.9*min(len(query),len(pdbseq)):raise SystemExit('Reference and structure chain do not align sufficiently')
    qpos=0;col_to_q={}
    for col,aa in enumerate(str(ref[0].seq)):
        if aa!='-':col_to_q[col]=qpos;qpos+=1
    ligand=Chem.SDMolSupplier(str(a.ligand),removeHs=True)[0]
    if ligand is None:raise SystemExit('Could not read ligand SDF')
    ligxyz=ligand.GetConformer().GetPositions();lig_heavy=[i for i,x in enumerate(ligand.GetAtoms()) if x.GetAtomicNum()>1];ligxyz=ligxyz[lig_heavy]
    contacts=set()
    for i,r in enumerate(pdbres):
        if i not in set(q_to_p.values()):continue
        atoms=[x for x in r if x.element not in ('H','D')]
        if atoms and min(sum((float(x.coord[k])-float(y[k]))**2 for k in range(3)) for x in atoms for y in ligxyz)**0.5<=a.contact_cutoff:contacts.add(r.id)
    rows=[]
    for col in range(alen):
        aas=[str(rec.seq[col]).upper() for rec in records];nongap=[x for x in aas if x in AA];coverage=len(nongap)/len(records)
        if nongap:
            counts={x:nongap.count(x) for x in set(nongap)};entropy=-sum((n/len(nongap))*math.log2(n/len(nongap)) for n in counts.values());norm=entropy/math.log2(20)
        else:entropy=norm=float('nan')
        q=col_to_q.get(col);pi=q_to_p.get(q) if q is not None else None;residue=pdbres[pi] if pi is not None and pi<len(pdbres) else None
        rows.append({'alignment_column':col+1,'reference_position':q+1 if q is not None else '', 'reference_aa':str(ref[0].seq[col]),'entropy_bits':entropy,'normalized_entropy':norm,'conservation':1-norm if math.isfinite(norm) else '', 'coverage':coverage,'n_sequences':len(records),'pdb_residue':f"{residue.id[1]}{residue.id[2].strip()}" if residue else '', 'pdb_aa':Polypeptide.protein_letters_3to1_extended.get(residue.resname.upper(),'') if residue else '', 'ap5_contact':bool(residue and residue.id in contacts)})
    with (a.out/'adk_conservation.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
    mapped=[r for r in rows if r['pdb_residue'] and r['conservation']!=''];site=[r for r in mapped if r['ap5_contact']];non=[r for r in mapped if not r['ap5_contact']]
    result={'reference':'E. coli K-12 adenylate kinase P69441','chain':a.chain,'n_sequences':len(records),'alignment_columns':alen,'structure_chain_residues':len(pdbres),'mapped_reference_positions':len(q_to_p),'site_cutoff_angstrom':a.contact_cutoff,'site_residue_count':len(site),'site_mean_conservation':sum(float(x['conservation']) for x in site)/len(site) if site else None,'noncontact_residue_count':len(non),'noncontact_mean_conservation':sum(float(x['conservation']) for x in non)/len(non) if non else None,'interpretation':'Descriptive, single protein family projection; no phylogenetic model, independence correction, or ligand-affinity claim.'}
    (a.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
