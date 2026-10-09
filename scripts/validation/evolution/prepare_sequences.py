#!/usr/bin/env python3
"""Select one reviewed adenylate-kinase sequence per species from UniProt TSV."""
import argparse,csv,json
from pathlib import Path
from Bio import SeqIO

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--reference',default='P69441');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 with a.input.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 needed={'Entry','Protein names','Organism','Length','Sequence'}
 if not rows or not needed.issubset(rows[0]):raise SystemExit('Input does not match expected UniProt TSV fields')
 ref=[r for r in rows if r['Entry']==a.reference]
 if len(ref)!=1:raise SystemExit(f'Expected one reference {a.reference}')
 target_length=int(ref[0]['Length']);by_species={}
 for r in rows:
  species=r['Organism'].split(' (')[0]
  rank=(abs(int(r['Length'])-target_length),r['Entry'])
  if species not in by_species or rank<by_species[species][0]:by_species[species]=(rank,r)
 selected=[ref[0]];seen={ref[0]['Sequence']};species_kept={ref[0]['Organism'].split(' (')[0]}
 for species,(rank,r) in sorted(by_species.items()):
  if species in species_kept:continue
  species_kept.add(species)
  if r['Sequence'] in seen:continue
  selected.append(r);seen.add(r['Sequence'])
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w') as f:
  for r in selected:f.write(f">{r['Entry']}|{r['Organism'].replace(' ','_')}"+('|reference' if r['Entry']==a.reference else '')+f"\n{r['Sequence']}\n")
 meta={'query_records':len(rows),'species_representatives':len(by_species),'exact_unique_sequences':len(selected),'reference':a.reference,'selection':'One reviewed entry per species minimizing absolute length difference from reference; accession is the deterministic tie break; exact duplicate sequences removed.','input_tsv':str(a.input)}
 a.output.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
if __name__=='__main__':main()
