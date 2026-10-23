"""Focused real-case and resource regression for the immutable v7 snapshot."""
import datetime, json, os, pathlib, re, shlex, subprocess, time
from full_driver import ROOT, PREVIOUS, run, table, sha, snapshot, ENV

V5 = ROOT/'outputs/bin/START_LST_wall_local_v5'
V6 = ROOT/'outputs/bin/START_LST_wall_local_v7'
V5HASH = '528d6d6a1b7aab34c73a70b09721295089ba293b2ef2c9bec173aba8e589c266'
V6HASH = '7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
COMMON = ['-pep_tol','1e-12','-pep_max_it','1000','-st_pc_factor_mat_solver_type','superlu']

def checked_case(spec):
    summary = run(spec)
    directory = ROOT/'runs'/spec['name']
    rows = table(directory/'out/neutral_field.plt')
    log = (directory/'solver.log').read_text(errors='replace')
    assert summary['petsc_errors']==0 and summary['unused_options']==0, summary
    assert all(abs(r['growth']+r['alpha_i'])<=1e-12 for r in rows if r['status']==0)
    return dict(summary=summary, field=rows, retry_calls=log.count('Candidate retry consuming owned reference:'),
                retry_preparations=log.count('Candidate retry prepared:'), native_index6=len(re.findall(r'MAC best index=\s*6\b', log)))

def compile_fixture():
    source = ROOT/'work/solver_source'; build = source/'build'
    env = dict(ENV, PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1', PETSC_ARCH='arch-linux-c-debug',
               SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',
               PKG_CONFIG_PATH='/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/lib/pkgconfig:/usr/local/lib/pkgconfig')
    env['PATH']='/usr/local/mpich-4.1.2/bin:'+env.get('PATH','')
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','PETSc','slepc','fortran_stdlib'],env=env,text=True))
    libs=shlex.split(subprocess.check_output(['pkg-config','--libs','PETSc','slepc','fortran_stdlib'],env=env,text=True))
    objects=sorted(p for p in (build/'src/lst/CMakeFiles/START_LST.dir').rglob('*.o') if p.name!='main.f90.o')
    exe=ROOT/'work/test_candidate_retry_v7'
    command=['/usr/local/mpich-4.1.2/bin/mpif90','-cpp','-fopenmp','-std=legacy','-I'+str(build/'src/lst')]+flags+[str(ROOT/'work/test_candidate_retry_v7.f90')]+list(map(str,objects))+libs+['-llapack','-lblas','-Wl,-rpath,/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/lib','-o',str(exe)]
    print('FIXTURE COMPILE '+json.dumps(command),flush=True)
    subprocess.run(command,env=env,check=True)
    command=[str(exe),'-pep_tol','1e-12','-st_pc_factor_mat_solver_type','superlu','-malloc_debug','-malloc_dump','-options_left','-memory_cycles','50']
    started=time.time(); process=subprocess.run(command,env=env,capture_output=True,text=True,timeout=120)
    text=process.stdout+process.stderr
    (ROOT/'outputs/candidate_retry_v7_fixture.log').write_text(text)
    result=dict(command=command,exit=process.returncode,seconds=time.time()-started,passed='CANDIDATE RETRY PASS:' in text,
                petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),
                final_malloc_dump_bytes=sum(map(int,re.findall(r'\[ 0\] (\d+) bytes',text))),
                memory_lines=[line for line in text.splitlines() if 'MEMORY CYCLES/' in line or 'RETRY MEMORY BEFORE/AFTER/DELTA:' in line])
    print('FIXTURE '+json.dumps(result),flush=True)
    assert result['exit']==0 and result['passed'] and result['petsc_errors']==0 and result['unused_options']==0
    assert result['final_malloc_dump_bytes']==0
    assert len(result['memory_lines'])==3 and all(abs(float(line.split()[-1]))<1e4 for line in result['memory_lines'])
    # Explicitly unsupported EPS must fail, before it can silently use PEP-only settings.
    reject_command=command+['-reject_eps']
    process=subprocess.run(reject_command,env=env,capture_output=True,text=True,timeout=120)
    text=process.stdout+process.stderr
    (ROOT/'outputs/candidate_retry_v7_eps_rejection.log').write_text(text)
    result['eps_rejection']=dict(command=reject_command,exit=process.returncode,
                                 explicit='Candidate retry rejected: EPS is unsupported' in text,
                                 unexpected_segv='Caught signal number 11' in text)
    assert process.returncode!=0 and result['eps_rejection']['explicit']
    assert not result['eps_rejection']['unexpected_segv']
    return result

def main():
    assert sha(V5)==V5HASH and sha(V6)==V6HASH
    snapshot()
    fixture=compile_fixture()
    seed=dict(omega_r=-1e-8,beta_r=-7.,alpha_r=51.76909408658165,alpha_i=-.5594665817891874)
    base=dict(mode='M1',coordinates='global',station=120,seed_station=119,seed=seed,task=33,
              outer_parameter='iloc',outer_min=119,outer_max=120,outer_points=2,
              inner_parameter='omega',inner_min=-1e-8,inner_max=1e-8,inner_points=3,
              adaptive=True,vector_tracking=True,sample_only=True,minimum_overlap=.95,step=.08,timeout=120)
    records={}
    for tag,binary,digest,nev,retry in [('v5_single',V5,V5HASH,1,0),('v5_all8',V5,V5HASH,8,0),('v7_retry8',V6,V6HASH,1,8)]:
        spec=dict(base,name='candidate_'+tag+'_m3_i119_120',solver=str(binary),solver_sha256=digest,
                  retry_eigenpairs=retry,options=COMMON+['-pep_nev',str(nev),'-pep_ncv','24' if nev==8 else '8'])
        records[tag]=checked_case(spec)
    expected=[r for r in records['v5_all8']['field'] if int(r['iloc'])==120 and r['status']==0]
    recovered=[r for r in records['v7_retry8']['field'] if int(r['iloc'])==120 and r['status']==0]
    assert len(expected)==3 and len(recovered)==3
    max_delta=max(abs(complex(a['alpha_r'],a['alpha_i'])-complex(b['alpha_r'],b['alpha_i'])) for a,b in zip(expected,recovered))
    single_failed=sum(int(r['iloc'])==120 and r['status']!=0 for r in records['v5_single']['field'])
    assert single_failed==3, single_failed
    assert max_delta<1e-8 and records['v7_retry8']['retry_calls']>0
    assert all(r['status']==0 for r in records['v7_retry8']['field'])
    assert records['v7_retry8']['summary']['rejected']==0
    # A small ordinary beta scan checks that retry support does not alter a
    # successful path when no expansion is needed. Same finite sample points.
    seed=next(r for r in json.loads((PREVIOUS/'outputs/seeds.json').read_text()) if r['station']==321)
    seed={k:seed[k] for k in ['omega_r','beta_r','alpha_r','alpha_i']}
    base=dict(mode='M1',coordinates='global',station=321,seed=seed,task=33,
              omega_min=20.,omega_max=20.,omega_points=1,beta_min=10.,beta_max=11.,beta_points=5,
              adaptive=True,vector_tracking=True,sample_only=True,minimum_overlap=.95,step=.08,timeout=120)
    beta={}
    for tag,binary,digest,nev,retry in [('v5_all8',V5,V5HASH,8,0),('v7_retry8',V6,V6HASH,1,8)]:
        beta[tag]=checked_case(dict(base,name='candidate_'+tag+'_i321_w20_beta10_11',solver=str(binary),
                                   solver_sha256=digest,retry_eigenpairs=retry,
                                   options=COMMON+['-pep_nev',str(nev),'-pep_ncv','24' if nev==8 else '8']))
    assert len(beta['v5_all8']['field'])==len(beta['v7_retry8']['field'])==5
    assert all(r['status']==0 for result in beta.values() for r in result['field'])
    max_beta_delta=max(abs(complex(a['alpha_r'],a['alpha_i'])-complex(b['alpha_r'],b['alpha_i'])) for a,b in zip(beta['v5_all8']['field'],beta['v7_retry8']['field']))
    assert max_beta_delta<1e-7
    frozen_v6=json.loads((ROOT/'outputs/candidate_retry_v6_regression.json').read_text())
    old_beta=frozen_v6['beta_scan']['v6_retry8']['field']
    max_delta_v6=max(abs(complex(a['alpha_r'],a['alpha_i'])-complex(b['alpha_r'],b['alpha_i'])) for a,b in zip(old_beta,beta['v7_retry8']['field']))
    assert max_delta_v6<1e-7
    output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
                solver=str(V6),solver_sha256=V6HASH,fixture=fixture,m3=records,beta_scan=beta,
                comparison=dict(m3_single_failed=single_failed,m3_max_delta_alpha=max_delta,beta_max_delta_alpha=max_beta_delta,
                                v6_v7_beta_max_delta_alpha=max_delta_v6),
                timing_scope='Sequential short cases under shared host load; actual wall times, no speed optimum claim.',passed=True)
    path=ROOT/'outputs/candidate_retry_v7_regression.json'
    assert not path.exists();path.write_text(json.dumps(output,indent=2)+'\n')
    snapshot();assert sha(V5)==V5HASH and sha(V6)==V6HASH
    print('CANDIDATE REGRESSION '+json.dumps(output['comparison']),flush=True)

if __name__=='__main__':main()
