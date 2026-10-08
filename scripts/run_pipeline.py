#!/usr/bin/env python3
"""Fail-closed M01–M21 orchestration with a separately authorized test-fixture branch."""
from p2d_analysis import *
from p2d_steps import STEPS,log_name
from datetime import datetime,timezone
import time

QC={'M04':'test_m01_m04_integration.py','M06':'test_m06_parsing.py','M07':'test_m07_ranking.py','M08':'test_m08_chemoinformatics.py','M09':'test_m09_candidate_gate.py',
    'M10A':'test_m10a_dynamics_jobs.py','M10C':'test_m10c_pose_extraction.py','M11':'test_m11_parameterization.py','M12':'test_m12_minimization.py','M13':'test_m13_equilibration.py','M14':'test_m14_production_md.py','M15':'test_m15_trajectory.py'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--profile',choices=['production','fixture'],default='production');p.add_argument('--resume-audit',choices=[s[0] for s in STEPS]);a=p.parse_args()
    os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH']
    fixture=a.profile=='fixture'; suffix='_fixture' if fixture else ''
    if fixture:
        for folder,name in [('receptors','test_receptor.pdb'),('ligands','test_ligand.sdf')]:
            path=ROOT/'input'/folder/name
            if not path.exists() or sha256(path)!=sha256(ROOT/'tests/fixtures'/folder/name): raise SystemExit('Fixture audit requires unchanged committed test inputs in input/. Existing user inputs are never replaced.')
            if len(list((ROOT/'input'/folder).glob('*')))!=1: raise SystemExit('Fixture audit requires isolated test inputs')
    runid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'); auditdir=ROOT/'logs'/('audit_'+a.profile+'_'+runid);auditdir.mkdir(parents=True)
    record={'profile':a.profile,'started_utc':runid,'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'steps':[],
            'fixture_branch':'Committed ethanol fixture bypasses candidate promotion only at M10C; production thresholds unchanged.' if fixture else None}
    auditfile=ROOT/'manifest'/('audit_'+a.profile+'.json');cfg=yaml.safe_load((ROOT/'config/dynamics.yaml').read_text());previous=None;failed=False
    start_index=0
    if a.resume_audit:
        previous_record=json.loads(auditfile.read_text())
        start_index=[s[0] for s in STEPS].index(a.resume_audit)
        completed=previous_record['steps'][:start_index]
        if len(completed)!=start_index or any(s['status'] not in {'PASS','SKIP','WARN'} or s.get('qc_returncode',0) for s in completed):
            raise SystemExit('Cannot resume audit: upstream steps did not pass')
        for step,(_,_,filename) in zip(completed,STEPS):
            path=ROOT/'logs'/log_name(step['module'],filename,a.profile)
            if sha256(path)!=step['log_sha256']: raise SystemExit('Cannot resume audit: upstream log changed')
        if fixture and start_index>=[s[0] for s in STEPS].index('M16'):
            for m in range(15,int(re.match(r'M(\d+)',a.resume_audit)[1])):
                for row in read_table(HEAVY/(DIRS[m]+'_fixture')/(NAMES[m]+'_manifest.tsv')):
                    verified(row['manifest'])
        record=previous_record
        record.setdefault('attempt_history',[]).append({'steps':record['steps'][start_index:],'status':record.get('status')})
        record['steps']=completed;record.pop('status',None)
        record.setdefault('resumptions',[]).append({'from':a.resume_audit,'time_utc':runid,'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()})
        previous=completed[-1]['status'] if completed else None
    # Existing fixture production has checkpoints; rebuilding M11–M13 invalidates them.
    if fixture and not a.resume_audit:
        old=HEAVY/'production_fixture'
        if old.exists():
            history=HEAVY/'audit_history'/runid;history.mkdir(parents=True,exist_ok=True);old.rename(history/'production_fixture')
    for label,_,filename in STEPS[start_index:]:
        (ROOT/'logs'/log_name(label,filename,a.profile)).write_text(f'FAIL_{label}_NOT_RUN\n')
    for module,script,filename in STEPS[start_index:]:
        logfile=ROOT/'logs'/log_name(module,filename,a.profile)
        logfile.write_text(f'FAIL_{module}_INCOMPLETE\n')
        cmd=[sys.executable,str(ROOT/'scripts'/script)]; number=int(re.match(r'M(\d+)',module)[1]);params=[]
        if fixture and module=='M10C': params=['--manifest',str(ROOT/'tests/fixtures/dynamics/dynamics_job_manifest.tsv'),'--output-root',str(ROOT/'work/dynamics/test_fixture'),'--log',str(logfile)]
        if fixture and number in {11,12,13,14}:
            source={11:ROOT/'work/dynamics/test_fixture/pose_extraction_manifest.tsv',12:HEAVY/'parameterized_fixture/parameterization_manifest.tsv',13:HEAVY/'minimized_fixture/minimization_manifest.tsv',14:HEAVY/'equilibrated_fixture/equilibration_manifest.tsv'}[number]
            dest={11:'parameterized_fixture',12:'minimized_fixture',13:'equilibrated_fixture',14:'production_fixture'}[number]
            params=['--manifest',str(source),'--output-root',str(HEAVY/dest),'--log',str(logfile)]
            if number==14:params+=['--profile','fixture']
        if number>=15:params=['--profile',a.profile]
        start=time.monotonic(); print(f'[{a.profile}] {module} {script}',flush=True)
        if not fixture and number>=11 and previous=='SKIP':
            logfile.write_text(f'SKIP_{module}_NO_UPSTREAM_JOBS\n')
            if number>=15:
                table(HEAVY/DIRS[number]/(NAMES[number]+'_manifest.tsv'),[],['dynamics_job_id','manifest','status','error'])
            rc=0
        else:
            with (auditdir/(module+'.txt')).open('w') as stream:
                proc=subprocess.run(cmd+params,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
            rc=proc.returncode
        first=logfile.read_text().splitlines()[0];status=first.split('_')[0]
        if rc or status not in {'PASS','SKIP','WARN'}:
            logfile.write_text(f'FAIL_{module}_EXECUTION\nExit: {rc}\nDetails: {auditdir/(module+".txt")}\n');status='FAIL';failed=True
        step={'module':module,'command':cmd+params,'returncode':rc,'status':status,'runtime_seconds':round(time.monotonic()-start,3),'log_sha256':sha256(logfile)}
        if fixture and not failed:
            qc=([sys.executable,str(ROOT/'tests'/QC[module])] if module in QC else [sys.executable,str(ROOT/'tests/test_analysis_modules.py'),str(number)] if number>=16 else None)
            if qc:
                with (auditdir/(module+'_qc.txt')).open('w') as stream: result=subprocess.run(qc,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
                step['qc_returncode']=result.returncode
                if result.returncode:
                    failed=True;logfile.write_text(f'FAIL_{module}_QC\n');step['status']='FAIL'
        record['steps'].append(step);write_json(auditfile,record);print(f'[{a.profile}] {module} {step["status"]}',flush=True)
        if failed:break
        previous=status
    record['status']='FAIL' if failed else 'PASS';record['completed_modules']=len(record['steps']);record['completed_git_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();write_json(auditfile,record)
    subprocess.run([sys.executable,str(ROOT/'scripts/update_status.py'),'--profile',a.profile],cwd=ROOT)
    print(f'{record["status"]}_PIPELINE_{a.profile.upper()} audit={auditfile}',flush=True)
    if failed:raise SystemExit(1)

if __name__=='__main__':main()
