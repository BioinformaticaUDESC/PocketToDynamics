"""Shared audited I/O for M15–M21. Paths remain in ~/PocketToDynamics."""
from pathlib import Path
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import argparse
import numpy as np
import yaml

ROOT = Path.home() / 'PocketToDynamics'
HEAVY = Path('/mnt/labio-fast-disk/PocketToDynamics/dynamics')
NAMES = {15:'preprocessing',16:'analysis',17:'energetics',18:'pca',19:'landscape',20:'interactions',21:'bundle'}
DIRS = {15:'preprocessed',16:'analyzed',17:'energetics',18:'pca',19:'landscape',20:'interactions',21:'bundles'}
LOGS = {15:'trajectory_preprocessing',16:'trajectory_analysis',17:'energetics',18:'pca',19:'landscape',20:'interactions',21:'bundle'}

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def write_json(path, data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n'); tmp.replace(path)

def table(path, rows, fields):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    tmp.replace(path)

def read_table(path):
    with Path(path).open() as f: return list(csv.DictReader(f,delimiter='\t'))

def command(args, out, name, stdin=None, env=None):
    args=list(map(str,args))
    p=subprocess.run(args,input=stdin,text=True,capture_output=True,cwd=out,env=env)
    (out/(name+'.stdout')).write_text(p.stdout or '')
    (out/(name+'.stderr')).write_text(p.stderr or '')
    with (out/'commands.jsonl').open('a') as f:
        f.write(json.dumps({'argv':args,'stdin':stdin,'returncode':p.returncode})+'\n')
    if p.returncode: raise RuntimeError(f'{name}: exit {p.returncode}: {p.stderr[-2000:]}')
    return p.stdout+'\n'+p.stderr

def verified(path):
    path=Path(path); data=json.loads(path.read_text())
    if data.get('status')!='PASS': raise ValueError(f'Upstream not PASS: {path}')
    for p,h in data.get('artifact_hashes',{}).items():
        if not Path(p).is_file() or sha256(p)!=h: raise ValueError(f'Artifact hash mismatch: {p}')
    if not data.get('artifact_hashes'): raise ValueError(f'No artifact hashes: {path}')
    return data

def artifact_manifest(out, module, job, inputs, files, details):
    hashes={str(Path(p).resolve()):sha256(p) for p in files}
    data={'schema_version':1,'module':module,'dynamics_job_id':job,'status':'PASS',
          'input_hashes':{str(Path(p).resolve()):sha256(p) for p in inputs},
          'artifact_hashes':hashes,'details':details,
          'software':{'python':sys.version,'numpy':np.__version__},
          'config_sha256':sha256(ROOT/'config/dynamics.yaml')}
    write_json(out/(NAMES[module]+'_manifest.json'),data)
    return data

def run_module(module, worker):
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--log',type=Path)
    parser.add_argument('--profile',choices=['production','fixture'],default='production')
    args=parser.parse_args()
    suffix='_fixture' if args.profile=='fixture' else ''
    prev=module-1
    default_input=(HEAVY/('production'+suffix)/'production_manifest.tsv' if module==15 else HEAVY/(DIRS[prev]+suffix)/(NAMES[prev]+'_manifest.tsv'))
    source=args.manifest or default_input
    out=args.output_root or HEAVY/(DIRS[module]+suffix)
    log=args.log or ROOT/'logs'/f'M{module}_{LOGS[module]}{suffix}.log'
    out=out.expanduser().resolve(); out.mkdir(parents=True,exist_ok=True)
    log.parent.mkdir(parents=True,exist_ok=True)
    log.write_text(f'FAIL_M{module}_INCOMPLETE\n')
    results=[]; failure=False
    try:
        if not source.exists(): raise FileNotFoundError(f'Upstream manifest missing: {source}')
        rows=read_table(source)
        cfg=yaml.safe_load((ROOT/'config/dynamics.yaml').read_text())
        seen=set()
        for row in rows:
            job=row['dynamics_job_id']
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',job) or job in {'.','..'} or job in seen: raise ValueError('Invalid/duplicate job id')
            seen.add(job); jobdir=out/job; jobdir.mkdir(exist_ok=True)
            marker=jobdir/(NAMES[module]+'_manifest.json')
            write_json(marker,{'status':'FAIL','reason':'INCOMPLETE'})
            try:
                if row.get('status')!='PASS': raise ValueError('Upstream job not PASS')
                worker(row,jobdir,cfg,args)
                verified(marker)
                results.append({'dynamics_job_id':job,'manifest':str(marker),'status':'PASS','error':''})
                print(f'[M{module}] {job} PASS',flush=True)
            except Exception as exc:
                failure=True
                write_json(marker,{'status':'FAIL','reason':str(exc)})
                results.append({'dynamics_job_id':job,'manifest':str(marker),'status':'FAIL','error':str(exc).replace('\n',' ')})
                print(f'[M{module}] {job} FAIL {exc}',flush=True)
    except Exception as exc:
        failure=True; print(f'[M{module}] FAIL {exc}',flush=True)
    table(out/(NAMES[module]+'_manifest.tsv'),results,['dynamics_job_id','manifest','status','error'])
    status=('FAIL' if failure else 'PASS' if results else 'SKIP')+f'_M{module}_{LOGS[module].upper()}'
    log.write_text(status+f'\nSystems: {len(results)}\nManifest: {out/(NAMES[module]+"_manifest.tsv")}\n')
    print(status,flush=True)
    if failure: raise SystemExit(1)

def preprocessing(data):
    return verified(data['details']['preprocessing_manifest'])

def universe(pre, fitted=True):
    import MDAnalysis as mda
    d=pre['details']; return mda.Universe(d['reference_fit'],d['fit_xtc' if fitted else 'centered_xtc'])

def groups(u, pre):
    return {name:u.atoms[np.asarray(idx,dtype=int)] for name,idx in pre['details']['indices'].items()}

def xvg(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith(('#','@')): rows.append([float(v) for v in line.split()])
    data=np.array(rows)
    if data.ndim!=2 or not len(data) or not np.isfinite(data).all(): raise ValueError(f'Invalid numeric XVG: {path}')
    return data
