"""Run six independent Vina searches; preserve every returned pose and command."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,json,time,hashlib,sys,platform
import argparse
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);args=parser.parse_args();root=args.root.resolve();root.mkdir(parents=True,exist_ok=True);
prep=root/'prepared';out=root/'docking'
if out.exists() and any(out.iterdir()):
 raise SystemExit('Refusing to overwrite an existing docking attempt. Use a fresh --root.')
out.mkdir(exist_ok=True)
p=json.loads((prep/'preparation.json').read_text())
plan={'receptors':['1AKE','4AKE'],'seeds':[20261007,20261008,20261009],'exhaustiveness':8,'num_modes':9,'energy_range':4,'cpu_per_run':2,'parallel_runs':6,'criterion':'Project criterion: symmetry-aware heavy-atom RMSD <= 2.0 angstrom without ligand superposition; top-ranked and best-returned evaluated separately. For 4AKE RMSD is descriptive only.','vina_version':subprocess.check_output(['/usr/bin/vina','--version'],text=True).strip(),'python':sys.version,'platform':platform.platform(),'input_sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in prep.iterdir() if q.suffix in ('.sdf','.pdbqt','.json')}}
(root/'protocol.json').write_text(json.dumps(plan,indent=2)+'\n')
def run(job):
 rec,seed=job;dest=out/f'{rec}_{seed}';dest.mkdir(exist_ok=False)
 cmd=['/usr/bin/vina','--receptor',str(prep/f'{rec}.pdbqt'),'--ligand',str(prep/'ligand.pdbqt'),'--out',str(dest/'poses.pdbqt'),'--cpu','2','--seed',str(seed),'--exhaustiveness','8','--num_modes','9','--energy_range','4']
 for key,values in [('center',p['box_center_angstrom']),('size',p['box_size_angstrom'])]:
  for axis,val in zip('xyz',values):cmd.extend([f'--{key}_{axis}',str(val)])
 start=time.monotonic()
 with (dest/'vina.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 info={'receptor':rec,'seed':seed,'command':cmd,'returncode':r.returncode,'runtime_seconds':time.monotonic()-start}
 if r.returncode==0:
  exp=[sys.executable,str(Path(sys.executable).parent/'mk_export.py'),str(dest/'poses.pdbqt'),'-s',str(dest/'poses.sdf')]
  r2=subprocess.run(exp,capture_output=True,text=True);(dest/'export.log').write_text(r2.stdout+r2.stderr);info['export_returncode']=r2.returncode
 info['sha256']={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in dest.iterdir()};(dest/'run.json').write_text(json.dumps(info,indent=2)+'\n');print(rec,seed,info,flush=True);return info
with ThreadPoolExecutor(max_workers=6) as pool:
 results=list(pool.map(run,[(r,s) for s in plan['seeds'] for r in plan['receptors']]))
assert all(x['returncode']==0 and x.get('export_returncode')==0 for x in results)
print('DOCKING_EXECUTION_PASS',flush=True)
