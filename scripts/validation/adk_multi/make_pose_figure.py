#!/usr/bin/env python3
"""Plot no-fit projections of native and top-ranked ADP poses for 8RJ4."""
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,default=Path('docs/validation/adk_multi'))
p.add_argument('--output',type=Path,default=Path('docs/paper/figures/adp_pose_overlay.png'))
a=p.parse_args(); root=a.root.resolve(); out=a.output.resolve(); out.parent.mkdir(parents=True,exist_ok=True)
fig,axes=plt.subplots(1,2,figsize=(10,4.6))
colors={'native':'#2878a8','docked':'#e87522'}
for ax,site in zip(axes,(301,302)):
    ref=Chem.SDMolSupplier(str(root/'inputs'/f'reference_{site}.sdf'),removeHs=True)[0]
    dock=Chem.SDMolSupplier(str(root/'poses'/f'{site}_20261010_poses.sdf'),removeHs=True)[0]
    if ref is None or dock is None: raise SystemExit(f'Could not read site {site} SDF')
    bonds=lambda m:{(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType())) for b in m.GetBonds()}
    if (ref.GetNumAtoms()!=dock.GetNumAtoms() or [x.GetAtomicNum() for x in ref.GetAtoms()] != [x.GetAtomicNum() for x in dock.GetAtoms()] or bonds(ref)!=bonds(dock)):
        raise SystemExit(f'Ligand graph/order mismatch at site {site}')
    xyz=np.asarray(ref.GetConformer().GetPositions());xyz2=np.asarray(dock.GetConformer().GetPositions())
    # Rotate both structures using axes derived from the crystal reference only.
    centered=xyz-xyz.mean(0);_,_,vt=np.linalg.svd(centered,full_matrices=False);basis=vt.T
    xy=centered@basis[:,:2];xy2=(xyz2-xyz.mean(0))@basis[:,:2]
    for mol,coords,key,style in ((ref,xy,'native','-'),(dock,xy2,'docked','--')):
        for bond in mol.GetBonds():
            i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            ax.plot(coords[[i,j],0],coords[[i,j],1],linestyle=style,color=colors[key],linewidth=1.6,alpha=.82,zorder=1)
        for atom in mol.GetAtoms():
            i=atom.GetIdx();elem=atom.GetSymbol()
            ax.scatter(coords[i,0],coords[i,1],s=34 if elem=='O' else 26,marker='o' if key=='native' else 's',
                       facecolor=colors[key] if elem=='C' else ('#c43c35' if elem=='O' else '#5369b2'),
                       edgecolor='white',linewidth=.45,zorder=2 if key=='native' else 3)
    ax.set_title(f'ADP site {site}',weight='bold');ax.set_xlabel('Reference principal axis 1 (Å)');ax.set_ylabel('Reference principal axis 2 (Å)')
    ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.16)
handles=[plt.Line2D([0],[0],color=colors['native'],lw=2,label='Crystallographic ADP'),
         plt.Line2D([0],[0],color=colors['docked'],lw=2,ls='--',label='Top-ranked pose, seed 20261010')]
fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.035))
fig.suptitle('Crystallographic and docked ADP coordinate projections',y=.99,fontsize=14)
fig.tight_layout(rect=[0,.08,1,.94]);fig.savefig(out,dpi=300,bbox_inches='tight',facecolor='white');fig.savefig(out.with_suffix('.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
print(out)
