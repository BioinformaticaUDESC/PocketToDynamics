#!/usr/bin/env python3
"""Build a deterministic reference-guided MSA from reviewed UniProt homologs."""
import argparse,json
from pathlib import Path
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--reference',default='P69441');p.add_argument('--output',type=Path,required=True);p.add_argument('--min-identity',type=float,default=.25);p.add_argument('--min-coverage',type=float,default=.75);a=p.parse_args()
 records=list(SeqIO.parse(a.input,'fasta')); refs=[r for r in records if r.id.split('|')[0]==a.reference]
 if len(refs)!=1:raise SystemExit(f'Expected exactly one {a.reference}, found {len(refs)}')
 ref=str(refs[0].seq).upper(); aligner=Align.PairwiseAligner();aligner.mode='global';aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-0.5
 accepted=[]; rejected=[]; maxins=[0]*(len(ref)+1)
 for record in records:
  seq=str(record.seq).upper();aln=aligner.align(ref,seq)[0];x,y=map(str,aln)
  pos=0;ins=['']*(len(ref)+1);mapped=['-']*len(ref);ident=paired=0
  for aa,bb in zip(x,y):
   if aa=='-':
    if bb!='-':ins[pos]+=bb
   else:
    if bb!='-':mapped[pos]=bb;paired+=1;ident+=aa==bb
    pos+=1
  coverage=paired/len(ref);identity=ident/paired if paired else 0
  if record.id.split('|')[0]==a.reference or (coverage>=a.min_coverage and identity>=a.min_identity):
   accepted.append((record,ins,mapped,identity,coverage));maxins=[max(m,len(s)) for m,s in zip(maxins,ins)]
  else:rejected.append((record.id,identity,coverage))
 if len(accepted)<2:raise SystemExit('Fewer than two sequences passed homology coverage filters')
 def render(ins,mapped):
  out=[]
  for i in range(len(ref)+1):
   out.extend(ins[i]+'-'*(maxins[i]-len(ins[i])))
   if i<len(ref):out.append(mapped[i])
  return ''.join(out)
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w') as f:
  for record,ins,mapped,identity,coverage in accepted:f.write(f'>{record.id}|pairwise_identity={identity:.4f}|reference_coverage={coverage:.4f}\n{render(ins,mapped)}\n')
 summary={'reference':a.reference,'input_sequences':len(records),'accepted_sequences':len(accepted),'rejected_sequences':len(rejected),'alignment_columns':len(render(accepted[0][1],accepted[0][2])),'method':'reference-guided global pairwise alignments; BLOSUM62; gap open -10; gap extension -0.5','minimum_reference_coverage':a.min_coverage,'minimum_pairwise_identity':a.min_identity,'identity_definition':'identical amino acids / paired non-gap aligned positions'}
 a.output.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
