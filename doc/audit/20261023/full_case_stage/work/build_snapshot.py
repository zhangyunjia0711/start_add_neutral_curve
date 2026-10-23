"""Build, test and freeze a named immutable solver snapshot."""
import argparse,datetime,json,os,shutil,subprocess,time
from full_driver import ROOT,sha,snapshot
p=argparse.ArgumentParser();p.add_argument('--tag',required=True);a=p.parse_args()
assert a.tag.isalnum()
source=ROOT/'work/solver_source'
env=dict(os.environ,PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PETSC_ARCH='arch-linux-c-debug',
 SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''))
env['PATH']='/usr/local/mpich-4.1.2/bin:'+env.get('PATH','')
commands=[['cmake','-S',str(source),'-B',str(source/'build'),'-DBUILD_TESTING=ON'],
 ['cmake','--build',str(source/'build'),'--target','START_LST','test_wall_local_profiles','test_slepc_adapter','test_neutral_numerics','test_slepc_candidate_retry','test_raw_vector_resources','--parallel','4'],
 ['ctest','--test-dir',str(source/'build'),'--output-on-failure','-R','wall_local|neutral_numerics|slepc_adapter_safety|slepc_candidate_retry|raw_vector_resources']]
snapshot();results=[]
for command in commands:
 print('BUILD COMMAND '+json.dumps(command),flush=True);started=time.time();process=subprocess.run(command,env=env)
 results.append(dict(command=command,exit=process.returncode,seconds=time.time()-started))
 if process.returncode:break
output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),commands=results,passed=len(results)==len(commands) and all(r['exit']==0 for r in results))
if output['passed']:
 target=ROOT/f'outputs/bin/START_LST_wall_local_{a.tag}'
 if target.exists():assert sha(target)==sha(source/'bin/START_LST'),'Immutable snapshot mismatch'
 else:shutil.copy2(source/'bin/START_LST',target)
 output.update(solver=str(target),sha256=sha(target),source_hashes={str(q.relative_to(source)):sha(q) for q in source.glob('src/**/*.f90')},
  git_head=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
  git_status=subprocess.check_output(['git','-C',str(source),'status','--short'],text=True))
 (ROOT/f'outputs/source_{a.tag}_tracked.patch').write_text(subprocess.check_output(['git','-C',str(source),'diff'],text=True))
path=ROOT/f'outputs/{a.tag}_build_validation.json'
if path.exists():path=ROOT/f'outputs/{a.tag}_build_validation_{int(time.time())}.json'
path.write_text(json.dumps(output,indent=2)+'\n');print('SNAPSHOT BUILD '+json.dumps({k:v for k,v in output.items() if k not in ['source_hashes','git_status']}),flush=True)
snapshot()
if not output['passed']:raise SystemExit(1)
