"""Build only new repository test targets; keep the frozen production solver intact."""
import datetime, json, os, pathlib, subprocess, time
from full_driver import ROOT, ENV, sha

source=ROOT/'work/solver_source';build=source/'build'
expected='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
assert sha(source/'bin/START_LST')==expected and sha(ROOT/'outputs/bin/START_LST_wall_local_v7')==expected
env=dict(ENV,PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PETSC_ARCH='arch-linux-c-debug',SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1')
env['PATH']='/usr/local/mpich-4.1.2/bin:'+env.get('PATH','')
commands=[['cmake','-S',str(source),'-B',str(build),'-DBUILD_TESTING=ON'],
          ['cmake','--build',str(build),'--target','test_slepc_candidate_retry','test_raw_vector_resources','--parallel','4'],
          ['ctest','--test-dir',str(build),'--output-on-failure','-R','slepc_candidate_retry|raw_vector_resources']]
records=[]
for command in commands:
    print('REPO TEST COMMAND '+json.dumps(command),flush=True);started=time.time()
    process=subprocess.run(command,env=env,capture_output=True,text=True)
    tag='configure' if len(records)==0 else 'build' if len(records)==1 else 'ctest'
    path=ROOT/f'outputs/v7_repo_{tag}.log'
    if path.exists():path=ROOT/f'outputs/v7_repo_{tag}_{int(time.time())}.log'
    path.write_text(process.stdout+process.stderr)
    print((process.stdout+process.stderr)[-3500:],flush=True)
    records.append(dict(command=command,exit=process.returncode,seconds=time.time()-started,log=str(path)))
    if process.returncode:break
result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),commands=records,
            passed=len(records)==3 and all(r['exit']==0 for r in records),solver_sha256=sha(source/'bin/START_LST'),
            frozen_solver_sha256=sha(ROOT/'outputs/bin/START_LST_wall_local_v7'))
assert result['solver_sha256']==result['frozen_solver_sha256']==expected
path=ROOT/'outputs/v7_repo_test_validation.json'
if path.exists():path=ROOT/f'outputs/v7_repo_test_validation_{int(time.time())}.json'
path.write_text(json.dumps(result,indent=2)+'\n');print('REPO TEST RESULT '+json.dumps(result),flush=True)
if not result['passed']:raise SystemExit(1)
