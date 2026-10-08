#!/usr/bin/env python3
"""Report all expected modules; never silently omit missing logs or mix fixture states."""
from pathlib import Path
import argparse,csv
from p2d_steps import STEPS,log_name
ROOT=Path.home()/'PocketToDynamics'

def main():
    p=argparse.ArgumentParser();p.add_argument('--profile',choices=['production','fixture'],default='production');a=p.parse_args()
    rows=[]
    for module,_,filename in STEPS:
        log=ROOT/'logs'/log_name(module,filename,a.profile)
        lines=log.read_text(errors='replace').splitlines() if log.exists() else []
        message=lines[0].strip() if lines else 'MISSING_LOG'
        status=message.split('_')[0]
        if status not in {'PASS','SKIP','FAIL','WARN'}: status='UNKNOWN'
        rows.append({'module':module,'status':status,'message':message,'log_file':str(log.relative_to(ROOT))})
    out=ROOT/'manifest'/('STATUS_fixture.tsv' if a.profile=='fixture' else 'STATUS.tsv');out.parent.mkdir(exist_ok=True)
    with out.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    for r in rows: print(r['module'],r['status'],r['message'])
    if any(r['status'] in {'FAIL','UNKNOWN'} for r in rows): raise SystemExit(1)

if __name__=='__main__':main()
