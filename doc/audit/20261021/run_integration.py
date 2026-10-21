import os,pathlib,subprocess,time,json,re,sys
root=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE');audit=root.parent/'neutral_curve_audit_20261021'
base=(root/'lst_test/lst_test/config.ini').read_text()
env=dict(os.environ,LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''))
standard=dict(enabled='T',task=37,outer_parameter='beta',inner_parameter='alpha',outer_min=-3.25,outer_max=-3.15,outer_points=3,inner_min=2.8,inner_max=4.5,scan_points=31,growth_tol='1e-7',parameter_tol='1e-9',eigen_error_tol='1e-7',max_eigen_jump=0.5,continuation_step=0.02,output_prefix='out/neutral')
cases=[]
for task in (31,32,33):
 cases.append((f'task{task}',dict(task=task,outer_parameter='iloc',outer_min=269,outer_max=271,outer_points=3),None,None,0))
cases += [('task37',{},None,None,0),('task36',dict(task=36,outer_min=-3.3,outer_max=-3.1,ds=0.02,ds_min='1e-5',ds_max=0.03,max_points=100,direction=1),None,None,0),('spatial37',dict(outer_min=-3.2059,outer_max=-3.2059,outer_points=1,inner_parameter='omega',inner_min=-0.35,inner_max=0.30), 'spatial',None,0),('interactive',dict(enabled='F',task=37,outer_min=-3.2059,outer_max=-3.2059,outer_points=1),None,'37 1 0\n0 0 0\n-1 0 0\n',0),('invalid_axes',dict(inner_parameter='omega'),None,None,1),('strict_jump',dict(max_eigen_jump='1e-12',outer_min=-3.2059,outer_max=-3.2059,outer_points=1,inner_min=3.49,inner_max=3.51,scan_points=3),None,None,1),('no_roots',dict(outer_min=-3.2059,outer_max=-3.2059,outer_points=1,inner_min=3.45,inner_max=3.55,scan_points=5),None,None,1)]
summary=[]
for name,updates,mode,stdin,expected in cases:
 case=audit/'real_cases'/name;case.mkdir(parents=True,exist_ok=True)
 if not (case/'in').exists(): (case/'in').symlink_to(root/'lst_test/lst_test/in',target_is_directory=True)
 options=dict(standard,**updates);cfg=base
 if mode=='spatial':
  cfg=re.sub(r'alpha = .*','alpha = 3.5d0, -0.04d0',cfg)
  cfg=re.sub(r'omega = .*','omega = -0.03933479054543d0, 0.0d0',cfg)
  cfg=cfg.replace('LST_TEMPORAL_MODE','LST_SPATIAL_MODE')
 cfg+='\n[neutral]\n'+''.join(f' {k} = {v}\n' for k,v in options.items())
 (case/'config.ini').write_text(cfg)
 cmd=[str(root/'bin/START_LST'),'-f','config.ini'];t=time.time()
 with (case/'solver.log').open('w') as out:
  out.write('COMMAND '+json.dumps(cmd)+'\n');out.flush()
  result=subprocess.run(cmd,input=stdin or '',text=True,cwd=case,env=env,stdout=out,stderr=subprocess.STDOUT,timeout=120)
 rows=[];f=case/'out/neutral_roots.plt'
 if f.exists():
  for line in f.read_text().splitlines():
   fields=line.split()
   if len(fields)==14:
    try: rows.append([float(v) for v in fields])
    except ValueError: pass
 record=dict(case=name,exit=result.returncode,expected=expected,roots=len(rows),max_abs_growth=max([abs(r[11]) for r in rows],default=None),seconds=round(time.time()-t,3))
 record['passed']=result.returncode==expected and (expected!=0 or bool(rows)) and all(abs(r[11])<=1e-7 and r[13]==0 for r in rows)
 if name=='interactive': record['passed'] &= bool(list((case/'out').glob('LST_eigenfunction*')))
 log=case/'out/neutral_run.log'
 if log.exists(): record['termination']=log.read_text().splitlines()[-1].strip()
 print(json.dumps(record),flush=True);summary.append(record)
 (audit/'integration_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 if not record['passed']: print('LOG TAIL\n'+'\n'.join((case/'solver.log').read_text().splitlines()[-20:]),flush=True)
if not all(r['passed'] for r in summary): sys.exit(1)
