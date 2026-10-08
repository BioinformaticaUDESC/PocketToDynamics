#!/usr/bin/env python3
"""Self-contained, checksummed software/QC bundle, not a manuscript."""
from p2d_analysis import *
import tarfile
import gzip
import io

def worker(row,out,cfg,args):
    prior=verified(row['manifest']); pre=preprocessing(prior); fixture=pre['details']['profile']=='fixture'; suffix='_fixture' if fixture else ''
    job=row['dynamics_job_id']; manifests=[]; current=Path(row['manifest'])
    links={20:'landscape_manifest',19:'pca_manifest',18:'energetics_manifest',17:'analysis_manifest',16:'preprocessing_manifest'}
    for m in range(20,14,-1):
        data=verified(current)
        if data['module']!=m or data['dynamics_job_id']!=job: raise ValueError('Broken module provenance chain')
        manifests.append((current,data))
        if m in links: current=Path(data['details'][links[m]])
    files=set()
    def include(path):
        path=Path(path)
        if path.is_file() and not path.is_symlink() and not path.name.startswith(('#','_MMPBSA','.')) and '__pycache__' not in path.parts:
            files.add(path.resolve())
    tracked=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    for p in tracked: include(ROOT/p)
    for path,data in manifests:
        include(path)
        for p in data['artifact_hashes']: include(p)
        for p,h in data['input_hashes'].items():
            # Config snapshots can change after upstream run; inputs used for provenance cannot.
            if p!=str(ROOT/'config/dynamics.yaml') and sha256(p)!=h: raise ValueError(f'Upstream input changed: {p}')
            include(p)
    for stage in ['parameterized','minimized','equilibrated','production']+[DIRS[m] for m in range(15,21)]:
        for p in (HEAVY/(stage+suffix)/job).rglob('*'): include(p)
    for folder in ['input','work','manifest','logs']:
        for p in (ROOT/folder).rglob('*'):
            if 'history' not in p.parts and not p.name.startswith('M21_') and p.name!='STATUS.tsv' and p.name!='STATUS_fixture.tsv': include(p)
    def relative(path):
        if path.is_relative_to(ROOT): return 'repository/'+str(path.relative_to(ROOT))
        if path.is_relative_to(HEAVY): return 'heavy/'+str(path.relative_to(HEAVY))
        raise ValueError(f'External artifact needs explicit bundle mapping: {path}')
    inventory=[{'original_path':str(p),'archive_path':relative(p),'sha256':sha256(p),'bytes':p.stat().st_size} for p in sorted(files)]
    versions={'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'git_status':subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True),
              'python':sys.version,'pip_freeze':subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True),
              'gromacs':pre['details']['gromacs_version'],
              'amber_mmpbsa':command([Path(cfg['energetics']['amber_root'])/'bin/MMPBSA.py','-v'],out,'amber_version')}
    write_json(out/'inventory.json',inventory);write_json(out/'versions.json',versions)
    description='Fixture: software validation only (20 ps); no convergence or binding claim.' if fixture else 'Production: inspect sampling and convergence before interpreting estimates.'
    report=f'# PocketToDynamics software and QC bundle\n\nJob: `{job}`\n\n{description}\n\nSource commit: `{versions["git_commit"]}`\n\n'
    report+='Modules M15–M20 have verified artifact hashes and PASS status. Coordinates, topologies, selected energetic frames, figures, tables, commands, configs and source are included.\n\n'
    for _,data in reversed(manifests):
        report+=f'## M{data["module"]}\n\n```json\n{json.dumps({k:v for k,v in data["details"].items() if k not in {"indices","gromacs_version","sampling"}},indent=2)}\n```\n\n'
    report+='## Verification and replay\n\nExtract the archive, then run `python bundle/verify_bundle.py` from its extracted root. The verifier checks every inventory entry. `repository/` contains the source and `heavy/` contains heavy artifacts. Original absolute paths are recorded in inventory.json. Replaying uses the documented ~/PocketToDynamics and /mnt/labio-fast-disk layout; relocate with an explicit path migration, not a blind textual replacement. See repository/README.md for commands. Binaries and conda environments are not embedded; versions and dependency specifications are included.\n'
    (out/'report.md').write_text(report)
    verifier='''from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
for entry in json.loads((root/'bundle/inventory.json').read_text()):
    path=root/entry['archive_path']
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''): h.update(chunk)
    if h.hexdigest()!=entry['sha256']: raise SystemExit('FAIL '+str(path))
print('PASS_BUNDLE_ALL_HASHES')
'''
    (out/'verify_bundle.py').write_text(verifier)
    archive=out/'reproducibility.tar.gz'; tmp=archive.with_suffix('.tmp')
    with tmp.open('wb') as stream,gzip.GzipFile(filename='',fileobj=stream,mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w|') as tar:
        entries=[(p,relative(p)) for p in sorted(files)]+[(out/f,'bundle/'+f) for f in ['inventory.json','versions.json','report.md','verify_bundle.py']]
        for path,name in entries:
            info=tar.gettarinfo(str(path),arcname=name);info.mtime=0;info.uid=info.gid=0;info.uname=info.gname=''
            with path.open('rb') as f: tar.addfile(info,f)
    tmp.replace(archive)
    with tarfile.open(archive) as tar:
        for entry in inventory:
            h=hashlib.sha256(); stream=tar.extractfile(entry['archive_path'])
            for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
            if h.hexdigest()!=entry['sha256']: raise ValueError('Archive verification failed')
    artifact_manifest(out,21,job,[row['manifest']],[archive,out/'inventory.json',out/'versions.json',out/'report.md',out/'verify_bundle.py'],
        {'preprocessing_manifest':prior['details']['preprocessing_manifest'],'interaction_manifest':row['manifest'],'files_packaged':len(inventory),
         'archive_sha256':sha256(archive),'profile':pre['details']['profile'],'archive_verified':True,'interpretation':description})

if __name__=='__main__': run_module(21,worker)
