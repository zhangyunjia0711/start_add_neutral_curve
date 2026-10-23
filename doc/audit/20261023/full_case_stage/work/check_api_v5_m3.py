"""Frozen-v5 global M3 continuation across one discrete station, no long scan."""
import json,pathlib,re
from full_driver import ROOT,run,table,sha
BINARY=ROOT/'outputs/bin/START_LST_wall_local_v5';EXPECTED='528d6d6a1b7aab34c73a70b09721295089ba293b2ef2c9bec173aba8e589c266'
assert sha(BINARY)==EXPECTED
seed=dict(omega_r=-1e-8,beta_r=-7.,alpha_r=51.76909408658165,alpha_i=-.5594665817891874)
spec=dict(name='api_v5_m3_global_i119_120',solver=str(BINARY),solver_sha256=EXPECTED,mode='M1',coordinates='global',station=120,seed_station=119,seed=seed,task=33,outer_parameter='iloc',outer_min=119,outer_max=120,outer_points=2,inner_parameter='omega',inner_min=-1e-8,inner_max=1e-8,inner_points=3,adaptive=True,vector_tracking=True,sample_only=True,minimum_overlap=.95,step=.08,options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','8','-pep_ncv','24','-st_pc_factor_mat_solver_type','superlu'],timeout=120)
s=run(spec);d=ROOT/'runs'/spec['name'];field=table(d/'out/neutral_field.plt');text=(d/'solver.log').read_text(errors='replace')
expected={-1e-8:complex(51.037936370808886,-.9777827733592761),0.:complex(51.037936458351986,-.977782755768114),1e-8:complex(51.03793654589409,-.9777827381797919)}
station120=[r for r in field if int(r['iloc'])==120 and r['status']==0]
errors=[abs(complex(r['alpha_r'],r['alpha_i'])-expected[min(expected,key=lambda om:abs(om-r['omega_r']))]) for r in station120]
result=dict(spec=spec,run_summary=s,field=field,station120_count=len(station120),max_delta_alpha_vs_independent_nev8=max(errors,default=None),all_growth_bound=sum(abs(r['growth']+r['alpha_i'])>1e-12 for r in field)==0,mac_index6_log_count=len(re.findall(r'MAC best index=\s*6\b',text)),log_has_mac_selection='MAC accepted=' in text)
result['accepted']=len(field)==6 and len(station120)==3 and s['rejected']==0 and s['petsc_errors']==0 and s['unused_options']==0 and all(r['status']==0 for r in field) and result['max_delta_alpha_vs_independent_nev8']<1e-8 and result['all_growth_bound'] and result['log_has_mac_selection']
p=ROOT/'outputs/api_v5_m3_continuation.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print('V5 M3 CONTINUATION '+json.dumps({k:result[k] for k in ['station120_count','max_delta_alpha_vs_independent_nev8','all_growth_bound','mac_index6_log_count','accepted']}),flush=True)
assert result['accepted']