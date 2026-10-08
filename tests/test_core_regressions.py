import importlib.util,sys,json,tempfile
from pathlib import Path
sys.path.insert(0,'/home/juan/PocketToDynamics/scripts')
def module(name):
    spec=importlib.util.spec_from_file_location('tested','/home/juan/PocketToDynamics/scripts/'+name+'.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
m=module('05c_run_screening')
with tempfile.TemporaryDirectory() as d:
    p=Path(d)/'poses';p.write_text('pose');f=Path(d)/'manifest.json'
    f.write_text(json.dumps({'task_id':'T','status':'PASS','returncode':0,'fingerprint':'current','output_sha256':m.sha256(p)}))
    assert m.valid_existing_task(f,p,'T','current')
    assert not m.valid_existing_task(f,p,'T','old')
    p.write_text('corrupted');assert not m.valid_existing_task(f,p,'T','current')
    f.write_text('{broken');assert not m.valid_existing_task(f,p,'T','current')
g=module('09_candidate_gate')
for value in ['nan','inf','-inf','',None]:
    try:g.as_float(value)
    except (ValueError,TypeError):pass
    else:raise AssertionError(value)
for value in ['1.2','0','-1']:
    try:g.as_int(value)
    except ValueError:pass
    else:raise AssertionError(value)
print('PASS_CORE_HASH_NUMERIC_REGRESSIONS')
