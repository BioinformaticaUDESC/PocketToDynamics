#!/usr/bin/env python3
"""Adversarial fixtures: corrupt outputs, invalid candidates, missing manifests, trajectory gaps."""
import sys,importlib.util,tempfile,json,subprocess
from pathlib import Path
ROOT=Path.home()/'PocketToDynamics';sys.path.insert(0,str(ROOT/'scripts'))
from p2d_analysis import *

def load(name):
    spec=importlib.util.spec_from_file_location('under_test',ROOT/'scripts'/name);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

with tempfile.TemporaryDirectory() as dirname:
    tmp=Path(dirname)
    artifact=tmp/'artifact';artifact.write_text('valid');manifest=tmp/'manifest.json'
    write_json(manifest,{'status':'PASS','artifact_hashes':{str(artifact):sha256(artifact)}})
    verified(manifest);artifact.write_text('tampered')
    try:verified(manifest)
    except ValueError:pass
    else:raise AssertionError('tampered output accepted')
    artifact.write_text('valid'); source=tmp/'source';source.write_text('old')
    write_json(manifest,{'status':'PASS','artifact_hashes':{str(artifact):sha256(artifact)},'input_hashes':{str(source):sha256(source)}})
    source.write_text('changed')
    try:verified(manifest)
    except ValueError:pass
    else:raise AssertionError('stale upstream accepted')
    # Any-mode must never promote NaN/Inf or malformed rank, even if other metrics pass.
    gate=load('09_candidate_gate.py')
    for mode in ['any','all']:
        gate.CONFIG=tmp/'gate.yaml';gate.CONFIG.write_text(yaml.safe_dump({'candidate_gate':{'logic':{'mode':mode},'criteria':{'max_integrated_rank':20,'max_best_score_kcal_mol':-7,'min_hit_fraction':0.25,'min_docking_ligand_efficiency':0.2}}}))
        gate.INPUT=tmp/'input.tsv';gate.OUT_DIR=tmp;gate.OUT_ALL=tmp/'all.tsv';gate.OUT_PROMOTED=tmp/'promoted.tsv';gate.SUMMARY=tmp/'summary.json';gate.MODULE_LOG=tmp/'module.log'
        base={'target_id':'T','ligand_id':'L','integrated_rank':'1','best_pocket':'P','best_score':'-10','hit_fraction':'1','docking_ligand_efficiency':'1'}
        for key,bad in [('best_score','nan'),('best_score','-inf'),('hit_fraction','1.2'),('integrated_rank','1.5')]:
            row=dict(base);row[key]=bad;table(gate.INPUT,[row],list(row))
            try:gate.main()
            except SystemExit as exc:assert exc.code==1
            else:raise AssertionError('Invalid candidate did not fail')
            assert not read_table(gate.OUT_PROMOTED)
    # Missing explicit upstream is FAIL, never SKIP; stale PASS marker is replaced.
    missing=subprocess.run([sys.executable,str(ROOT/'scripts/16_analyze_trajectory.py'),'--manifest',str(tmp/'missing.tsv'),'--output-root',str(tmp/'output'),'--log',str(tmp/'missing.log')],capture_output=True,text=True)
    assert missing.returncode==1 and (tmp/'missing.log').read_text().startswith('FAIL_')
    # Segment corruption must not be reported as completed.
    production=load('14_production_md.py');segment=tmp/'segment';segment.mkdir()
    for ext in ['gro','cpt','edr','xtc','log','tpr']:(segment/('md.'+ext)).write_text(ext)
    write_json(segment/'segment_complete.json',{'status':'PASS','sha256':{ext:sha256(segment/('md.'+ext)) for ext in ['gro','cpt','edr','xtc','log','tpr']}})
    assert production.segment_is_complete(segment)
    (segment/'md.xtc').write_text('corrupt');assert not production.segment_is_complete(segment)
    # Synthetic trajectory is valid only if every expected time is present.
    import MDAnalysis as mda
    trajectory=load('15_preprocess_trajectory.py')
    u=mda.Universe.empty(3,trajectory=True);u.add_TopologyAttr('names',['CA','N','C']);u.add_TopologyAttr('resnames',['ALA']);u.add_TopologyAttr('resids',[1])
    u.atoms.positions=[[0,0,0],[1,0,0],[0,1,0]];u.dimensions=[20,20,20,90,90,90];u.atoms.write(str(tmp/'reference.gro'))
    for times,valid in [([0,1,2],True),([0,1,3],False),([1,2,3],False)]:
        with mda.Writer(str(tmp/'test.xtc'),n_atoms=3) as writer:
            for t in times:u.trajectory.ts.time=t;writer.write(u.atoms)
        try:trajectory.trajectory_qc(tmp/'reference.gro',tmp/'test.xtc',times[-1])
        except ValueError:assert not valid
        else:assert valid
print('PASS_PIPELINE_FAILURE_REGRESSIONS')
