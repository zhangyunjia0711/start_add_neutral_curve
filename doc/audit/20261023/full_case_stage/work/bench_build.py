"""Validate independent solver build against archived original-binary benchmark."""
import configparser, datetime, json, pathlib, re, statistics, subprocess, time
from full_driver import ROOT, PREVIOUS, CASE, ENV, run, table
from bench_speed import inspect, regular, selected_residuals, FLOAT
NEW=ROOT/'work/solver_source/bin/START_LST'
TIGHT=['-pep_tol','1e-12','-pep_max_it','1000']
FAST=TIGHT+['-st_pc_factor_mat_solver_type','superlu','-pep_ncv','8']

def first_selected(log):
    pairs=[]
    for line in log.read_text(errors='replace').splitlines():
        if 'ncov=' in line:pairs=[]
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if m:pairs.append(dict(index=int(m[1]),alpha_r=float(m[2]),alpha_i=float(m[3]),error=float(m[4])))
        if 'alpha------------>' in line:
            v=[float(x) for x in re.findall(FLOAT,line.split('alpha------------>',1)[1])]
            assert len(v)==2 and pairs
            selected=min(pairs,key=lambda p:abs(complex(p['alpha_r']-v[0],p['alpha_i']-v[1])))
            return dict(selected=selected,pairs=pairs)
    raise AssertionError('No initial selected eigenpair')

def manual_copy(name,source,flag,enabled='T'):
    d=ROOT/'runs'/name;d.mkdir(exist_ok=True)
    cfg=configparser.ConfigParser();cfg.read(source)
    cfg['neutral']['track_initial_target']=flag;cfg['neutral']['enabled']=enabled
    with (d/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(NEW),'-f','config.ini']+FAST+['-pep_nev','3','-options_left'];t=time.time()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        p=subprocess.run(cmd,cwd=d,env=ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=120)
    a=selected_residuals(d/'solver.log');first=first_selected(d/'solver.log')
    s=dict(name=name,flag=flag,neutral_enabled=enabled,exit=p.returncode,seconds=time.time()-t,audit=a,first=first)
    (d/'validation.json').write_text(json.dumps(s,indent=2)+'\n');print('TARGET CONTROL '+json.dumps(s),flush=True)
    return s

def main():
    reference=inspect(json.loads((ROOT/'runs/bench_tight_r1/summary.json').read_text()))
    seed=next(v for v in json.loads((PREVIOUS/'outputs/seeds.json').read_text()) if v['station']==321)
    groups={'default':TIGHT,'serial_ncv8':FAST};results=[]
    for repeat in [1,2]:
        for group in (list(groups) if repeat==1 else list(reversed(groups))):
            s=run(dict(name='buildbench_%s_r%d'%(group,repeat),solver=str(NEW),mode='M1',station=321,omega_min=20.,omega_max=20.,omega_points=1,beta_min=-12.,beta_max=12.,beta_points=25,step=.16,seed=seed,options=groups[group]+['-options_left'],timeout=180))
            a=inspect(s);a['group']=group;a['repeat']=repeat;p=a['point'];q=reference['point'];rr=reference['regular']
            a['comparison']=dict(delta_peak_growth=abs(p['growth']-q['growth']),delta_peak_beta=abs(p['beta_r']-q['beta_r']),delta_peak_alpha=abs(complex(p['alpha_r']-q['alpha_r'],p['alpha_i']-q['alpha_i'])),regular_count=len(a['regular']),max_regular_growth_difference=max(abs(v['growth']-rr[b]['growth']) for b,v in a['regular'].items()),max_regular_alpha_difference=max(abs(complex(v['alpha_r']-rr[b]['alpha_r'],v['alpha_i']-rr[b]['alpha_i'])) for b,v in a['regular'].items()))
            c=a['comparison'];a['accepted']=a['valid'] and c['regular_count']==25 and c['delta_peak_growth']<1e-8 and c['delta_peak_alpha']<1e-6 and c['delta_peak_beta']<1e-3 and c['max_regular_growth_difference']<1e-8 and c['max_regular_alpha_difference']<1e-8
            del a['regular'];results.append(a)
    exact=dict(iloc=321.,x=.6143816857,omega_r=22.61994307483593,omega_i=0.,beta_r=10.57775824698773,beta_i=0.,alpha_r=44.02144804943075,alpha_i=5.812488106864006e-8)
    spec=dict(name='buildbench_initial_target_true',solver=str(NEW),mode='M1',station=321,task=33,omega_min=exact['omega_r'],omega_max=exact['omega_r'],omega_points=1,beta_min=exact['beta_r']-.001,beta_max=exact['beta_r']+.001,beta_points=3,step=.16,seed=exact,options=FAST+['-pep_nev','3','-options_left'],timeout=120)
    s=run(spec);audit=selected_residuals(ROOT/'runs'/spec['name']/'solver.log');first=first_selected(ROOT/'runs'/spec['name']/'solver.log')
    target=complex(exact['alpha_r']+.001,exact['alpha_i']);nearest=min(first['pairs'],key=lambda p:abs(complex(p['alpha_r'],p['alpha_i'])-target))
    control=manual_copy('buildbench_initial_target_false',ROOT/'runs'/spec['name']/'config.ini','F')
    interactive=manual_copy('buildbench_non_neutral_flag_true',ROOT/'runs'/spec['name']/'config.ini','T','F')
    checks=dict(opt_in_nearest=first['selected']==nearest,opt_in_target_error=abs(complex(first['selected']['alpha_r']-exact['alpha_r'],first['selected']['alpha_i']-exact['alpha_i'])),opt_in_error=first['selected']['error'],opt_in_audit=audit,default_growth_selection=control['first']['selected']==max(control['first']['pairs'],key=lambda p:-p['alpha_i']),non_neutral_flag_ignored=interactive['first']['selected']==max(interactive['first']['pairs'],key=lambda p:-p['alpha_i']))
    target_accepted=checks['opt_in_nearest'] and checks['opt_in_target_error']<1e-8 and checks['opt_in_error']<1e-7 and audit['petsc_error_lines']==0 and not audit['unconsumed_benchmark_options'] and checks['default_growth_selection'] and checks['non_neutral_flag_ignored'] and control['audit']['petsc_error_lines']==0 and interactive['audit']['petsc_error_lines']==0
    times={name:[v['seconds'] for v in results if v['group']==name] for name in groups};med={name:statistics.median(v) for name,v in times.items()}
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),binary=str(NEW),results=results,times_seconds=times,median_seconds=med,speedup_new_default_vs_serial_ncv8=med['default']/med['serial_ncv8'],all_benchmarks_accepted=all(v['accepted'] for v in results),initial_target=dict(checks=checks,accepted=target_accepted,first=first,false_control=control,non_neutral_control=interactive),recommended_options=FAST)
    (ROOT/'outputs/solver_build_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('BUILD VALIDATION SUMMARY '+json.dumps({k:result[k] for k in ['times_seconds','median_seconds','speedup_new_default_vs_serial_ncv8','all_benchmarks_accepted']}),flush=True)
    print('INITIAL TARGET CHECKS '+json.dumps(dict(accepted=target_accepted,checks=checks)),flush=True)
    assert result['all_benchmarks_accepted'] and target_accepted
    cmd=[str(ROOT/'work/solver_source/bin/test_slepc_adapter'),'-st_pc_factor_mat_solver_type','superlu','-options_left']
    print('SERIAL ADAPTER TEST COMMAND '+json.dumps(cmd),flush=True)
    p=subprocess.run(cmd,env=ENV,capture_output=True,text=True,timeout=60);(ROOT/'outputs/serial_adapter_test.log').write_text(p.stdout+p.stderr)
    assert p.returncode==0 and 'PETSC ERROR' not in p.stdout+p.stderr and 'Option left:' not in p.stdout+p.stderr and 'PASS:' in p.stdout
    print('SERIAL ADAPTER TEST PASS',flush=True)

if __name__=='__main__':main()