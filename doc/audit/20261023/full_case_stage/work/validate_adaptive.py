import concurrent.futures,json,shutil
from full_driver import ROOT,run,table,sha,snapshot
binary=ROOT/'work/solver_source/bin/START_LST'
target=ROOT/'outputs/bin/START_LST_adaptive_v2'
if not target.exists():shutil.copy2(binary,target)
assert sha(binary)==sha(target)
axial=json.loads((ROOT/'outputs/axial_seeds.json').read_text())['records']
trial=json.loads((ROOT.parent/'neutral_beta_envelope_20261022/outputs/seeds.json').read_text())
seeds=[('nose',119,min((r for r in axial if r['family']=='M3' and r['status']==0),key=lambda r:abs(r['iloc']-119)),0,-12,-2),('neutral',321,next(r for r in trial if r['station']==321),20,-12,12)]
specs=[]
for label,i,seed,om,lo,hi in seeds:
 for step,adaptive in [(.08,False),(.4,True)]:
  specs.append(dict(name=f'adaptcheck_{label}_{"adaptive" if adaptive else "fine"}',station=i,seed={k:seed[k] for k in ['omega_r','beta_r','alpha_r','alpha_i']},solver=str(target),solver_sha256=sha(target),task=31,omega_min=om,omega_max=om,omega_points=1,beta_min=lo,beta_max=hi,beta_points=25,step=step,adaptive=adaptive,options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu']))
snapshot()
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(run,specs))
checks=[]
for label,*_ in seeds:
 fine=next(r for r in results if r['spec']['name']==f'adaptcheck_{label}_fine');adaptive=next(r for r in results if r['spec']['name']==f'adaptcheck_{label}_adaptive')
 def grid(s):
  spec=s['spec'];d={}
  for r in table(ROOT/'runs'/spec['name']/'out/neutral_field.plt'):
   frac=(r['beta_r']-spec['beta_min'])*24/(spec['beta_max']-spec['beta_min'])
   if abs(frac-round(frac))<1e-8 and r['status']==0:d[round(frac)]=complex(r['alpha_r'],r['alpha_i'])
  return d
 a,b=grid(fine),grid(adaptive);common=set(a)&set(b)
 checks.append(dict(label=label,common_samples=len(common),max_delta_alpha=max((abs(a[k]-b[k]) for k in common),default=None),fine_rejected=fine['rejected'],adaptive_rejected=adaptive['rejected'],petsc_errors=fine['petsc_errors']+adaptive['petsc_errors'],fine_peak=fine['peaks'],adaptive_peak=adaptive['peaks']))
output=dict(solver=str(target),sha256=sha(target),checks=checks,results=results)
(ROOT/'outputs/adaptive_validation.json').write_text(json.dumps(output,indent=2)+'\n');print('ADAPTIVE CHECKS '+json.dumps(checks),flush=True)
snapshot()
