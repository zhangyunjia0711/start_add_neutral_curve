"""Logged, cached START scan executor for the complete-case calculation."""
import configparser, datetime, hashlib, json, math, os, pathlib, re, subprocess, time
ROOT=pathlib.Path(__file__).resolve().parent.parent
CASE=ROOT.parent
PREVIOUS=CASE/'neutral_beta_envelope_20261022'
FIXED=CASE/'neutral_xomega_20261022'
SOLVER=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE/bin/START_LST')
ENV=dict(os.environ,LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
COLUMNS=['outer','inner','iloc','x','Re','alpha_r','alpha_i','beta_r','beta_i','omega_r','omega_i','growth','branch','status']
DEFAULT_OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_ncv','8']

def sha(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for v in iter(lambda:f.read(1048576),b''):h.update(v)
    return h.hexdigest()

def snapshot():
    original=json.loads((FIXED/'input_manifest.json').read_text());checks=[]
    for m in original:
        for p,e in [(m['config'],m['config_sha256'])]+[(v['path'],v['sha256']) for v in m['files']]:
            actual=sha(p);assert actual==e,p;checks.append(dict(path=p,sha256=actual))
    actual=sha(SOLVER);assert actual=='807169bb401e558d9509c11f75304e5694d1d8af6925660fca251e8bc29f469a'
    checks.append(dict(path=str(SOLVER),sha256=actual));return checks

def table(path):
    rows=[];zone=0
    if not path.exists():return rows
    for line in path.read_text().splitlines():
        m=re.search(r'ZONE\s+T="row_(\d+)"',line)
        if m:zone=int(m[1]);continue
        try:
            v=[float(x.replace('D','E')) for x in line.split()]
            if len(v)==14:rows.append(dict(zip(COLUMNS,v),row=zone))
        except ValueError:pass
    return rows

def run(spec):
    d=ROOT/'runs'/spec['name'];d.mkdir(parents=True,exist_ok=True)
    if (d/'summary.json').exists():
        s=json.loads((d/'summary.json').read_text());assert s['spec']==spec,'cache mismatch';return s
    assert not (d/'out/neutral_run.log').exists(),'Incomplete existing run '+str(d)
    mode=spec.get('mode','M1');cfg=configparser.ConfigParser();cfg.read(CASE/mode/'config.ini')
    cfg['filenames']['dir']=str(CASE/mode/'in')+'/'
    seed=spec['seed'];cfg['domain']['iloc']=str(spec.get('seed_station',spec['station']))
    cfg['lst']['profile_coordinates']=spec.get('coordinates','global')
    cfg['lst']['omega']='%.15g, 0'%seed['omega_r'];cfg['lst']['beta']='%.15g, 0'%seed['beta_r']
    cfg['lst']['alpha']='%.15g, %.15g'%(seed['alpha_r']+spec.get('seed_shift',1e-3),seed['alpha_i'])
    cfg['neutral']={k:str(v) for k,v in dict(enabled='T',track_initial_target='T',task=spec.get('task',31),outer_parameter=spec.get('outer_parameter','omega'),inner_parameter=spec.get('inner_parameter','beta'),outer_min=spec.get('outer_min',spec.get('omega_min')),outer_max=spec.get('outer_max',spec.get('omega_max')),outer_points=spec.get('outer_points',spec.get('omega_points')),inner_min=spec.get('inner_min',spec.get('beta_min')),inner_max=spec.get('inner_max',spec.get('beta_max')),scan_points=spec.get('inner_points',spec.get('beta_points')),growth_tol=spec.get('growth_tol',1e-8),parameter_tol=spec.get('parameter_tol',1e-7),max_iterations=70,eigen_error_tol=spec.get('eigen_tol',1e-7),max_eigen_jump=spec.get('jump',2),continuation_step=spec.get('step',.08),output_prefix='out/neutral').items()}
    cfg['neutral']['adaptive_steps']='T' if spec.get('adaptive',False) else 'F'
    cfg['neutral']['minimum_step']=str(spec.get('minimum_step',1e-5))
    cfg['neutral']['vector_tracking']='T' if spec.get('vector_tracking',False) else 'F'
    cfg['neutral']['minimum_overlap']=str(spec.get('minimum_overlap',.95))
    cfg['neutral']['sample_only']='T' if spec.get('sample_only',False) else 'F'
    for key in ['reverse_outer','reverse_inner','anchor_row_seed']:
        cfg['neutral'][key]='T' if spec.get(key,False) else 'F'
    cfg['neutral']['retry_eigenpairs']=str(spec.get('retry_eigenpairs',0))
    cfg['neutral']['station_predictor']='T' if spec.get('station_predictor',False) else 'F'
    with (d/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(spec.get('solver',SOLVER)),'-f','config.ini']+spec.get('options',DEFAULT_OPTIONS)
    binary_hash=sha(cmd[0])
    if 'solver_sha256' in spec:assert spec['solver_sha256']==binary_hash,'Solver binary changed'
    started=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat();t=time.time()
    with (d/'solver.log').open('w') as out:
        out.write('COMMAND '+json.dumps(cmd)+'\n');out.flush()
        try:p=subprocess.run(cmd,cwd=d,env=ENV,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,timeout=spec.get('timeout',2400));code=p.returncode
        except subprocess.TimeoutExpired:code=124
    path=d/'out/neutral_run.log';log=path.read_text() if path.exists() else ''
    m=re.search(r'evaluations=\s*(\d+)\s+solves=\s*(\d+)\s+rejected=\s*(\d+)\s+result=(.+)',log)
    peaks=table(d/'out/neutral_ridge.plt');field=table(d/'out/neutral_field.plt')
    solver_log=(d/'solver.log').read_text(errors='replace')
    errors=len(re.findall(r'PETSC ERROR|SLEPC ERROR',solver_log,re.I))
    unused=len(re.findall(r'Option left:|options you set that were not used',solver_log,re.I))
    missing_coordinates=spec.get('coordinates','global')=='wall_local' and not (d/'wall_local_profile_metadata.csv').exists()
    s=dict(spec=spec,solver_sha256=binary_hash,command=cmd,started=started,seconds=round(time.time()-t,3),exit=code,terminal=m.group(4).strip() if m else 'MISSING',evaluations=int(m.group(1)) if m else None,solves=int(m.group(2)) if m else None,rejected=int(m.group(3)) if m else None,peaks=peaks,field_records=len(field),failed_field_records=sum(v['status']!=0 for v in field),petsc_errors=errors,unused_options=unused)
    s['coordinate_metadata_missing']=missing_coordinates
    (d/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
    print(json.dumps(dict(name=spec['name'],seconds=s['seconds'],terminal=s['terminal'],rejected=s['rejected'],peaks=[{k:v[k] for k in ['iloc','omega_r','beta_r','alpha_r','alpha_i','growth']} for v in peaks])),flush=True)
    return s
