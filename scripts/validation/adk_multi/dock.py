#!/usr/bin/env python3
"""Run fixed-seed Vina redocking for two crystallographic ADP sites in 8RJ4."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse,hashlib,json,platform,subprocess,sys,time

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve();prep=root/'prepared';out=root/'docking'
    if out.exists() and any(out.iterdir()):raise SystemExit('Refusing to overwrite a docking attempt; choose a new root')
    out.mkdir(parents=True,exist_ok=True);meta=json.loads((prep/'preparation.json').read_text());seeds=[20261010,20261011,20261012]
    protocol={'entry':'8RJ4','chain':'A','ligand':'ADP','sites':[301,302],'seeds':seeds,'vina_version':subprocess.check_output(['/usr/bin/vina','--version'],text=True).strip(),'exhaustiveness':16,'num_modes':9,'energy_range_kcal_mol':4,'cpu_per_run':2,'box_policy':meta['box_policy'],'site_specific_threshold':'Top pose no-fit symmetry-aware heavy-atom RMSD <= 2.0 A, evaluated per crystallographic reference site.','input_sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in prep.iterdir() if q.is_file() and q.suffix in ('.sdf','.pdbqt','.json')}}
    (root/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    def run(job):
        resid,seed=job;dest=out/f'{resid}_{seed}';dest.mkdir(exist_ok=False)
        cmd=['/usr/bin/vina','--receptor',str(prep/'receptor.pdbqt'),'--ligand',str(prep/f'ligand_{resid}.pdbqt'),'--out',str(dest/'poses.pdbqt'),'--cpu','2','--seed',str(seed),'--exhaustiveness','16','--num_modes','9','--energy_range','4']
        for key,vals in [('center',meta['sites'][str(resid)]['center_angstrom']),('size',meta['sites'][str(resid)]['size_angstrom'])]:
            for axis,val in zip('xyz',vals):cmd.extend([f'--{key}_{axis}',str(val)])
        t=time.monotonic()
        with (dest/'vina.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
        recorded_cmd=[str(Path(x).relative_to(root)) if isinstance(x,str) and x.startswith(str(root)+'/') else x for x in cmd]
        info={'site_residue':resid,'seed':seed,'command':recorded_cmd,'command_working_directory':'.','returncode':r.returncode,'runtime_seconds':time.monotonic()-t}
        if not r.returncode:
            exp=[sys.executable,str(Path(sys.executable).parent/'mk_export.py'),str(dest/'poses.pdbqt'),'-s',str(dest/'poses.sdf')]
            r2=subprocess.run(exp,capture_output=True,text=True);(dest/'export.log').write_text(r2.stdout+r2.stderr);info['export_returncode']=r2.returncode
        info['sha256']={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in dest.iterdir() if q.is_file()};(dest/'run.json').write_text(json.dumps(info,indent=2)+'\n')
        if r.returncode or info.get('export_returncode')!=0:raise RuntimeError(f'Vina failed {resid}/{seed}; inspect {dest}')
        return info
    with ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(run,[(r,s) for r in (301,302) for s in seeds]))
    print(json.dumps(results,indent=2));print('ADP_DOCKING_EXECUTION_PASS')
if __name__=='__main__':main()
