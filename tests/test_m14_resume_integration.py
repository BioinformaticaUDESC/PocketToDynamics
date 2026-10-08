#!/usr/bin/env python3
"""Real fixture-only checkpoint resume and completed-segment reuse regression."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path.home()/'PocketToDynamics';HEAVY=Path('/mnt/labio-fast-disk/PocketToDynamics/dynamics')
job=HEAVY/'production_fixture/MDTEST_a5b7f2f48479'
command=[sys.executable,str(ROOT/'scripts/14_production_md.py'),'--profile','fixture','--manifest',str(HEAVY/'equilibrated_fixture/equilibration_manifest.tsv'),'--output-root',str(HEAVY/'production_fixture'),'--log',str(ROOT/'logs/M14_resume_test.log')]
first=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
assert first.returncode==0,first.stdout+first.stderr
assert first.stdout.count('SKIP_COMPLETED')==2,first.stdout
marker=job/'segment_0002/segment_complete.json';before=marker.read_text()
(job/'resume_test_original_marker.json').write_text(before)
marker.unlink()  # Remove only the completion flag; all simulation files remain intact.
second=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
assert second.returncode==0,second.stdout+second.stderr
assert 'RESUME' in second.stdout,second.stdout
result=json.loads(marker.read_text());assert result['mode']=='RESUME' and result['end_ps']==20
for ext,expected in result['sha256'].items():
    assert hashlib.sha256((marker.parent/('md.'+ext)).read_bytes()).hexdigest()==expected
(job/'resume_test_evidence.json').write_text(json.dumps({'status':'PASS','skip_returncode':first.returncode,'resume_returncode':second.returncode,'skip_stdout':first.stdout,'resume_stdout':second.stdout},indent=2)+'\n')
print('PASS_M14_REAL_SKIP_AND_CHECKPOINT_RESUME')
