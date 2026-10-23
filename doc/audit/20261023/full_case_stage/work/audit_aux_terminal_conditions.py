#!/usr/bin/env python3
"""Audit terminal J/parameter steps and actual solver acceptance, not physical degeneration."""
import json
import numpy as np
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import stamp

def collect_names(obj,names):
    if isinstance(obj,dict):
        for k,v in obj.items():
            if isinstance(v,str) and (k=='name' or k.endswith('_run')):
                if (ROOT/'runs'/v/'summary.json').exists():names.add(v)
            if k=='runs' and isinstance(v,list):
                names.update(n for n in v if isinstance(n,str) and (ROOT/'runs'/n/'summary.json').exists())
            collect_names(v,names)
    elif isinstance(obj,list):
        for v in obj:collect_names(v,names)

def main():
    source=ROOT/'outputs/aux_neutral_M1p_17_16_s64_20261023.json';d=json.loads(source.read_text())
    out=ROOT/'outputs/aux_M1p_terminal_numerics_20261023.json'
    if out.exists():raise RuntimeError('Existing audit preserved')
    records=d['records'][-3:];failure=d['failures'][-1];names=set()
    collect_names(records,names);collect_names(failure,names)
    summaries=[]
    for name in sorted(names):
        s=json.loads((ROOT/'runs'/name/'summary.json').read_text())
        summaries.append(dict(name=name,minimum_step=s['spec'].get('minimum_step'),
            evaluations=s['evaluations'],solves=s['solves'],rejected=s['rejected'],field_records=s['field_records'],
            failed_field_records=s['failed_field_records'],petsc_errors=s['petsc_errors'],exit=s['exit'],terminal=s['terminal']))
    def compact(t):
        return dict(name=t['name'],h=t['h'],center=t['center'],F=t['F'],J=t['J'],
            condition_number_2=float(np.linalg.cond(t['J'])),unlimited_newton_delta=t.get('unlimited_newton_delta'),
            trust_scale=t.get('trust_scale'),stencil_refinement_reason=t.get('stencil_refinement_reason'),
            line_search=[dict(factor=x.get('factor'),merit=x.get('merit'),status=x.get('status'),
                center=x.get('stencil',{}).get('center'),F=x.get('stencil',{}).get('F'),J=x.get('stencil',{}).get('J')) for x in t.get('line_search',[])],
            local_peak_reseed=t.get('local_peak_reseed'))
    report=dict(created=stamp(),source=str(source),source_sha256=sha(source),auxiliary_only=True,
        last_three=[dict(auxiliary_column=r['auxiliary_column'],theta=r['theta'],point=r['point'],
            verification_h=r['result']['verification']['h'],verification_F=r['result']['verification']['F'],
            verification_J=r['result']['verification']['J'],
            verification_J_condition_2=float(np.linalg.cond(r['result']['verification']['J'])),
            verification_correction=r['result']['verification_correction'],
            trajectory=[compact(t) for t in r['result']['trace']]) for r in records],
        failure=dict(auxiliary_column=failure['auxiliary_column'],theta=failure['theta'],reason=failure['reason'],
            prediction=failure['prediction'],transported_actual=failure['transported_actual'],predicted_actual=failure['predicted_actual'],
            trajectory=[compact(t) for t in failure['result']['trace']]),solver_summaries=summaries,
        findings=dict(actual_minimum_step_values=sorted(set(s['minimum_step'] for s in summaries)),
            number_of_run_summaries=len(summaries),sum_rejected_requested=sum(s['rejected'] or 0 for s in summaries),
            sum_failed_requested_fields=sum(s['failed_field_records'] for s in summaries),
            tiny_parameter_samples_all_accepted=all(s['failed_field_records']==0 and s['rejected']==0 for s in summaries),
            stopping_stage='Positive-curvature predictor followed by no unique interior beta maximum on actual refined window; no failed-point Newton delta or line-search was taken'),
        qualification='Numerical Jacobian and accepted-sample audit only. A large condition number or a predictor curvature sign is not evidence of physical degeneration or absence of native roots.')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps(report['findings']),flush=True)
    for r in report['last_three']:print(json.dumps({k:v for k,v in r.items() if k!='trajectory'}),flush=True)

if __name__=='__main__':main()
