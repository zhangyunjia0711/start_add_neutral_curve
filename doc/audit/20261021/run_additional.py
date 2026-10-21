import os,pathlib,subprocess,json,re,time
root=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE');audit=root.parent/'neutral_curve_audit_20261021'
base=(root/'lst_test/lst_test/neutral_curve.ini').read_text()
env=dict(os.environ,LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''))
summary=[]
for name,changes,mpi,expected in [('re37',{'outer_parameter':'Re','outer_min':79000,'outer_max':79200,'outer_points':3},False,0),('station_stride',{'task':31,'outer_parameter':'iloc','outer_min':269,'outer_max':271,'outer_points':2},False,0),('mpi_rejection',{},True,1)]:
 case=audit/'real_cases'/name;case.mkdir(exist_ok=True)
 if not (case/'in').exists(): (case/'in').symlink_to(root/'lst_test/lst_test/in',target_is_directory=True)
 cfg=base
 for key,value in changes.items(): cfg=re.sub(r'(?m)^\s*'+re.escape(key)+r' = .*$',f'  {key} = {value}',cfg)
 (case/'config.ini').write_text(cfg)
 cmd=[str(root/'bin/START_LST'),'-f','config.ini']
 if mpi: cmd=['/usr/local/mpich-4.1.2/bin/mpiexec','-n','2']+cmd
 with (case/'solver.log').open('w') as out:
  out.write('COMMAND '+json.dumps(cmd)+'\n');out.flush()
  result=subprocess.run(cmd,cwd=case,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,timeout=90)
 rows=[];p=case/'out/neutral_roots.plt'
 if p.exists():
  for line in p.read_text().splitlines():
   try:
    if len(line.split())==14: rows.append([float(v) for v in line.split()])
   except ValueError: pass
 passed=result.returncode==expected and (mpi or len(rows)>0) and all(abs(r[11])<=1e-7 for r in rows)
 if mpi: passed &= 'Neutral curves require a single MPI process.' in (case/'solver.log').read_text()
 record=dict(case=name,exit=result.returncode,roots=len(rows),passed=passed);print(json.dumps(record),flush=True);summary.append(record)
(audit/'additional_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if not all(r['passed'] for r in summary): raise SystemExit(1)
