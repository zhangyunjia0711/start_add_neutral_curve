"""Small reproducible solver speed comparisons; original executable unchanged."""
import datetime, json, math, pathlib, re, statistics, subprocess
from full_driver import ROOT, PREVIOUS, run, table, DEFAULT_OPTIONS

FLOAT = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'

def selected_residuals(path):
    pairs=[]; selected=[]; unmatched=[]; unused=[]; petsc_errors=[]
    for line in path.read_text(errors='replace').splitlines():
        if 'ncov=' in line: pairs=[]
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if m:
            pairs.append((int(m[1]),complex(float(m[2].replace('D','E')),float(m[3].replace('D','E'))),float(m[4].replace('D','E'))))
        if 'alpha------------>' in line:
            values=re.findall(FLOAT,line.split('alpha------------>',1)[1])
            if len(values)==2 and pairs:
                actual=complex(*(float(v.replace('D','E')) for v in values))
                best=min(pairs,key=lambda p:abs(p[1]-actual))
                if abs(best[1]-actual)<1e-6:selected.append(best[2])
                else:unmatched.append(line)
            else:unmatched.append(line)
        if 'Option left:' in line:unused.append(line.strip())
        if 'PETSC ERROR:' in line:petsc_errors.append(line.strip())
    return dict(petsc_error_lines=len(petsc_errors),first_petsc_error_lines=petsc_errors[:8],selected_count=len(selected),max_selected_error=max(selected,default=None),selected_errors_over_tol=sum(v>1e-7 for v in selected),unmatched=unmatched,unused_option_lines=len(unused),unused_options=sorted(set(unused)),unconsumed_benchmark_options=[v for v in unused if any(k in v for k in ['pep_tol','pep_max_it','pep_ncv','st_pc_factor_mat_solver_type'])])

def continuity(field):
    issues=[]
    for omega in sorted(set(v['omega_r'] for v in field)):
        rows=sorted((v for v in field if v['omega_r']==omega and v['status']==0),key=lambda v:v['beta_r'])
        for a,b in zip(rows,rows[1:]):
            db=abs(b['beta_r']-a['beta_r']);da=abs(complex(b['alpha_r']-a['alpha_r'],b['alpha_i']-a['alpha_i']))
            if db<.02 and da>max(1e-4,30*db):issues.append(dict(beta_a=a['beta_r'],beta_b=b['beta_r'],delta_alpha=da,delta_beta=db))
    return issues

def regular(field):
    rows={}
    for v in field:
        beta=round(v['beta_r'])
        if -12<=beta<=12 and abs(v['beta_r']-beta)<1e-9 and v['status']==0:rows[beta]=v
    return rows

def inspect(s):
    d=ROOT/'runs'/s['spec']['name'];field=table(d/'out/neutral_field.plt');audit=selected_residuals(d/'solver.log')
    issues=continuity(field);good=s['exit'] in [0,1] and s['terminal'] in ['converged','no bracketed root'] and s['rejected']==0 and s['failed_field_records']==0 and len(s['peaks'])==1 and not issues and not audit['unconsumed_benchmark_options'] and not audit['unmatched'] and audit['petsc_error_lines']==0 and audit['selected_errors_over_tol']==0 and audit['selected_count']>0
    return dict(name=s['spec']['name'],seconds=s['seconds'],exit=s['exit'],terminal=s['terminal'],solves=s['solves'],evaluations=s['evaluations'],rejected=s['rejected'],point=s['peaks'][0] if len(s['peaks'])==1 else None,continuity_issues=issues,residual_audit=audit,valid=good,regular=regular(field),options=s['spec']['options'])

def host():
    cfg=pathlib.Path('/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/include/petscconf.h').read_text()
    return dict(cpu=subprocess.run(['lscpu'],capture_output=True,text=True).stdout,memory=subprocess.run(['free','-m'],capture_output=True,text=True).stdout,petsc_archs=[str(v) for v in pathlib.Path('/home/zaza/ProgramFiles/petsc-3.22.1').glob('arch*') if v.is_dir()],petsc_debug='#define PETSC_USE_DEBUG 1' in cfg,superlu_installed='#define PETSC_HAVE_SUPERLU 1' in cfg,superlu_dist_installed='#define PETSC_HAVE_SUPERLU_DIST 1' in cfg)

def main():
    seed=next(v for v in json.loads((PREVIOUS/'outputs/seeds.json').read_text()) if v['station']==321)
    groups={'tight':list(DEFAULT_OPTIONS),'superlu':list(DEFAULT_OPTIONS)+['-st_pc_factor_mat_solver_type','superlu'],'superlu_ncv8':list(DEFAULT_OPTIONS)+['-st_pc_factor_mat_solver_type','superlu','-pep_ncv','8'],'tight_ncv8':list(DEFAULT_OPTIONS)+['-pep_ncv','8']}
    results=[]
    for repeat in [1,2]:
        for name in (list(groups) if repeat==1 else list(reversed(groups))):
            spec=dict(name='bench_%s_r%d'%(name,repeat),mode='M1',station=321,omega_min=20.,omega_max=20.,omega_points=1,beta_min=-12.,beta_max=12.,beta_points=25,step=.16,seed=seed,options=groups[name]+['-options_left'],timeout=300)
            s=run(spec);a=inspect(s);a['group']=name;a['repeat']=repeat;results.append(a)
    reference=next(v for v in results if v['name']=='bench_tight_r1')
    reference_regular=dict(reference['regular'])
    for a in results:
        p=a['point'];q=reference['point']
        a['comparison']=dict(delta_peak_growth=abs(p['growth']-q['growth']),delta_peak_beta=abs(p['beta_r']-q['beta_r']),delta_peak_alpha=abs(complex(p['alpha_r']-q['alpha_r'],p['alpha_i']-q['alpha_i'])),regular_beta_count=len(a['regular']),max_regular_growth_difference=max((abs(v['growth']-reference_regular[b]['growth']) for b,v in a['regular'].items() if b in reference_regular),default=None),max_regular_alpha_difference=max((abs(complex(v['alpha_r']-reference_regular[b]['alpha_r'],v['alpha_i']-reference_regular[b]['alpha_i'])) for b,v in a['regular'].items() if b in reference_regular),default=None))
        c=a['comparison'];a['equivalent_to_reference']=a['valid'] and c['regular_beta_count']==25 and c['delta_peak_growth']<1e-7 and c['delta_peak_beta']<1e-3 and c['delta_peak_alpha']<1e-3 and c['max_regular_growth_difference']<1e-6 and c['max_regular_alpha_difference']<1e-6
        del a['regular']
    times={name:[a['seconds'] for a in results if a['group']==name] for name in groups}
    medians={name:statistics.median(v) for name,v in times.items()}
    eligible=[name for name in groups if all(a['equivalent_to_reference'] for a in results if a['group']==name)]
    recommendation=min(eligible,key=lambda name:medians[name]) if eligible else None
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),scope='M1 i321 omega20 beta[-12,12] N25 continuation .16, two sequential repeats per group; no recompilation',host=host(),option_evidence=dict(factor_solver='Outer PEP ST KSP prefix st_: SLEPc3.22 stsles.c404-405; PETSc factor options factimpl.c246',eps_prefix='PEPLINEAR internal EPS prefix pep_linear_: linear.c920-921; PEP dimensions/tol forwarded linear.c369-371',sources_root='/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/externalpackages/slepc-3.22.0',warning='Target-mode speed benchmark, not universal spectral equivalence proof. Timings include system load and I/O; selected pair residuals independently parsed. Serial SuperLU option is rejected when adapter hardcoded superlu_dist conflicts after setup; no performance recommendation from runs containing PETSc errors.'),results=results,times_seconds=times,median_seconds=medians,speedup_vs_tight={name:medians['tight']/v for name,v in medians.items()},recommendation=recommendation,recommended_options=groups[recommendation] if recommendation else None)
    path=ROOT/'outputs/speed_summary.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print('SPEED SUMMARY '+json.dumps({k:result[k] for k in ['times_seconds','median_seconds','speedup_vs_tight','recommendation','recommended_options']}),flush=True)

if __name__=='__main__':main()