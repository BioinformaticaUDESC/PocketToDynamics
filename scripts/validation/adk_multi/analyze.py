#!/usr/bin/env python3
"""Verify ADP docking artifacts, compute no-fit symmetry RMSD and summarize runs."""
from pathlib import Path
import argparse,csv,hashlib,json,re
import sys
import numpy as np
from rdkit import Chem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from adk.rmsd import matching_graph,fixed_frame_rmsd

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve();output=a.output.resolve();output.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((root/'protocol.json').read_text())
    for name,digest in protocol['input_sha256'].items():
        if hashlib.sha256((root/'prepared'/name).read_bytes()).hexdigest()!=digest:raise SystemExit(f'Input hash mismatch: {name}')
    expected={f'{r}_{s}' for r in protocol['sites'] for s in protocol['seeds']};actual={q.name for q in (root/'docking').iterdir()}
    if expected!=actual:raise SystemExit(f'Unexpected docking directories: expected {expected}, got {actual}')
    rows=[];summaries=[]
    for site in protocol['sites']:
        ref=Chem.SDMolSupplier(str(root/'prepared'/f'reference_{site}.sdf'),removeHs=True)[0]
        if ref is None:raise SystemExit(f'Cannot read ADP reference at site {site}')
        for seed in protocol['seeds']:
            dest=root/'docking'/f'{site}_{seed}';run=json.loads((dest/'run.json').read_text())
            if run['returncode']!=0 or run.get('export_returncode')!=0:raise SystemExit(f'Incomplete run {site}/{seed}')
            for name,digest in run['sha256'].items():
                if hashlib.sha256((dest/name).read_bytes()).hexdigest()!=digest:raise SystemExit(f'Output hash mismatch {site}/{seed}/{name}')
            energies=[float(x) for x in re.findall(r'REMARK VINA RESULT:\s+([-\d.]+)',(dest/'poses.pdbqt').read_text())]
            poses=list(Chem.SDMolSupplier(str(dest/'poses.sdf'),removeHs=True))
            if len(poses)!=len(energies) or not poses:raise SystemExit(f'Pose count mismatch {site}/{seed}')
            rms=[]
            ref_smiles=Chem.MolToSmiles(matching_graph(ref))
            for rank,(mol,energy) in enumerate(zip(poses,energies),1):
                if mol is None or mol.GetNumAtoms()!=ref.GetNumAtoms():raise SystemExit(f'Bad pose {site}/{seed}/{rank}')
                if Chem.MolToSmiles(matching_graph(mol))!=ref_smiles:raise SystemExit(f'Graph/stereochemistry mismatch {site}/{seed}/{rank}')
                dist=fixed_frame_rmsd(mol,ref);rms.append(dist)
                rows.append({'site_residue':site,'seed':seed,'rank':rank,'vina_score_kcal_mol':energy,'no_fit_symmetry_rmsd_angstrom':dist})
            summaries.append({'site_residue':site,'seed':seed,'n_poses':len(poses),'top_score_kcal_mol':energies[0],'top_rmsd_angstrom':rms[0],'best_rmsd_angstrom':min(rms),'best_rmsd_rank':int(np.argmin(rms)+1),'top_pose_recovered_le_2A':bool(rms[0]<=2.0),'any_pose_recovered_le_2A':bool(min(rms)<=2.0),'runtime_seconds':run['runtime_seconds']})
    for filename,data in [('poses.csv',rows),('summary.csv',summaries)]:
        with (output/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    (output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    colors={301:'#087e8b',302:'#d96b27'};fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for ax,key,title,ylabel in [(axes[0],'top_score_kcal_mol','Top ranked docking score','Vina score (kcal/mol)'),(axes[1],'top_rmsd_angstrom','Top ranked pose recovery','No fit symmetry RMSD (Å)')]:
        for site in protocol['sites']:
            v=[x for x in summaries if x['site_residue']==site];ax.plot(range(1,4),[x[key] for x in v],'o-',color=colors[site],label=f'ADP {site}',lw=1.6)
        ax.set_xticks([1,2,3],protocol['seeds']);ax.set_xlabel('Vina search seed');ax.set_ylabel(ylabel);ax.set_title(title)
    axes[1].axhline(2,color='#555',ls='--',lw=1,label='2 Å project criterion');axes[0].legend(frameon=False);axes[1].legend(frameon=False,fontsize=8)
    fig.suptitle('ADP known-site redocking in 8RJ4 chain A',fontsize=13)
    for ext in ('png','svg'):fig.savefig(output/f'adp_redocking.{ext}',dpi=220)
    plt.close(fig);print(json.dumps(summaries,indent=2));print('ADP_ANALYSIS_PASS')
if __name__=='__main__':main()
