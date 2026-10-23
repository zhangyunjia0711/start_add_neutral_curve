import json
import compensate as c
plan=json.loads((c.AREA/'compensation_plan.json').read_text());j=next(j for j in plan['jobs'] if j['name']=='i053_r00025_both')
r=c.beta_line(j['station'],j['start'],j['reference']['beta_r'],'append_protocol_fixture_bridge')
p=c.endpoint(r,j['cap'],j['reference']['beta_r']);info,q=c.vec(r,p)
fresh,reference,v=c.fresh(j['station'],p,'append_protocol_fixture_reference')
result=dict(actual_point=p,bridge=c.brief(r),bridge_vector=info,reference=c.brief(fresh),reference_vector=reference,
 mac=c.mac(q,v),passed=c.mac(q,v)>=1-1e-8,documentation='https://petsc.org/release/manualpages/Sys/PetscObjectViewFromOptions.html',
 note='Append preserves each Solve, including requested field endpoint and later restored seed state. Vector/log candidate counts must agree exactly.')
c.write('append_vector_protocol_validation',result);print(json.dumps({k:v for k,v in result.items() if k not in ('bridge','reference','bridge_vector','reference_vector')}))

