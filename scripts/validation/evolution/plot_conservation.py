#!/usr/bin/env python3
"""Create a descriptive conservation figure from a mapped residue table."""
import argparse,csv
from pathlib import Path
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 rows=[r for r in csv.DictReader(a.input.open(),delimiter='\t') if r['pdb_residue'] and r['conservation']]
 pos=[int(r['reference_position']) for r in rows];cons=[float(r['conservation']) for r in rows];site=[r['ap5_contact']=='True' for r in rows]
 fig,(ax,bx)=plt.subplots(2,1,figsize=(10,6.5),gridspec_kw={'height_ratios':[2,1]},constrained_layout=True)
 ax.plot(pos,cons,color='#9aa6ac',lw=1,zorder=1);ax.scatter([x for x,s in zip(pos,site) if not s],[y for y,s in zip(cons,site) if not s],s=19,color='#32869a',label='Outside 5 Å contact shell',zorder=2);ax.scatter([x for x,s in zip(pos,site) if s],[y for y,s in zip(cons,site) if s],s=26,color='#d36a35',label='Within 5 Å of Ap5A',zorder=3)
 ax.set(xlabel='E. coli adenylate kinase reference position',ylabel='Conservation (1 − normalized entropy)',ylim=(-.03,1.03),title='Evolutionary sequence conservation mapped to the Ap5A-bound structure');ax.legend(frameon=False,loc='lower right')
 v1=[y for y,s in zip(cons,site) if s];v2=[y for y,s in zip(cons,site) if not s]
 bx.boxplot([v1,v2],tick_labels=[f'Ap5A contact\nn={len(v1)}',f'Other mapped residues\nn={len(v2)}'],showfliers=False,patch_artist=True,boxprops={'facecolor':'#dcebef'},medianprops={'color':'#273f48'});bx.set(ylabel='Conservation',ylim=(-.03,1.03));bx.text(.99,.03,'Descriptive comparison; residues are not independent evolutionary observations.',transform=bx.transAxes,ha='right',fontsize=8,color='#555')
 a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=220);fig.savefig(a.output.with_suffix('.svg'));plt.close(fig)
if __name__=='__main__':main()
