import datetime,hashlib,json,pathlib,shutil,subprocess
from full_driver import ROOT,snapshot
S=ROOT/'work/solver_source';target=ROOT/'outputs/bin/START_LST_vectors_v4';assert not target.exists()
assert json.loads((ROOT/'outputs/raw_vector_regression_attempt3.json').read_text())['passed']
shutil.copy2(S/'bin/START_LST',target)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(target)==sha(S/'bin/START_LST')
subprocess.run(['git','-C',str(S),'diff','--check'],check=True)
patch=subprocess.run(['git','-C',str(S),'diff'],capture_output=True,text=True,check=True).stdout
out=ROOT/'outputs/source_v4_snapshot.patch';assert not out.exists();out.write_text(patch)
paths=['src/public/module_slepc_adapter.f90','src/lst/lst_eqn.f90','src/lst/lst_solver_abstract.f90','src/lst/module_lst.f90','src/lst/module_neutral_curve.f90']
result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),solver=str(target),sha256=sha(target),source_hashes={p:sha(S/p) for p in paths},git_head=subprocess.run(['git','-C',str(S),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip(),git_branch=subprocess.run(['git','-C',str(S),'branch','--show-current'],capture_output=True,text=True,check=True).stdout.strip(),git_status=subprocess.run(['git','-C',str(S),'status','--short'],capture_output=True,text=True,check=True).stdout,raw_regression='outputs/raw_vector_regression_attempt3.json',original_integrity=snapshot())
out=ROOT/'outputs/solver_v4_freeze.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print('FROZEN V4 '+json.dumps(result),flush=True)
source=(ROOT/'work/validate_vectors.py').read_text().replace('START_LST_vectors_v3','START_LST_vectors_v4').replace('veccheck_','v4veccheck_').replace('vector_validation.json','vector_validation_v4.json')
source=source.replace('snapshot()\n',"assert all(c['common_samples']==25 and c['max_delta_alpha']<=1e-8 and c['fine_rejected']==0 and c['adaptive_rejected']==0 and c['petsc_errors']==0 for c in checks)\nassert all(r['rejected']==0 and r['petsc_errors']==0 and r['unused_options']==0 and r['terminal'] in ['converged','no bracketed root'] for r in results)\nsnapshot()\n",1) if False else source
source+='\nassert all(c["common_samples"]==25 and c["max_delta_alpha"]<=1e-8 and c["fine_rejected"]==0 and c["adaptive_rejected"]==0 and c["petsc_errors"]==0 for c in checks)\nassert all(r["rejected"]==0 and r["petsc_errors"]==0 and r["unused_options"]==0 and r["terminal"] in ["converged","no bracketed root"] for r in results)\nprint("V4 VALIDATION ALL CHECKS PASS",flush=True)\n'
out=ROOT/'work/validate_vectors_v4.py';assert not out.exists();out.write_text(source)