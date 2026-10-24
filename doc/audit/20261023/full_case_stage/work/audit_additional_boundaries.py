"""Bind each exported additional neutral point to its actual h/2 field file."""
import datetime,json
from full_driver import ROOT,sha,table,snapshot


def main():
    output=ROOT/'outputs';source=output/'neutral_variable_beta_candidates.json'
    d=json.loads(source.read_text());records=[];families={}
    solver=output/'bin/START_LST_wall_local_v10';solver_hash=sha(solver)
    for name,kind in [('tail841_lobe_lower','lower_local'),('tail841_low_lobe_upper','upper_local')]:
        rows=[r for r in d['unique_candidates_data'] if r['source']==name and r['frequency_boundary']==kind]
        for r in rows:
            run=ROOT/'runs'/r['verification_run'];summary_path=run/'summary.json'
            summary=json.loads(summary_path.read_text());field=run/'out/neutral_field.plt'
            assert summary['solver_sha256']==solver_hash==r['solver_sha256']
            assert not summary['petsc_errors'] and not summary['unused_options']
            assert not summary['failed_field_records'] and not summary['coordinate_metadata_missing']
            points=table(field);h=r['fd_h'];grid={}
            for q in points:
                ia=round((q['omega_r']-r['omega_r'])/h)
                ib=round((q['beta_r']-r['beta_r'])/h)
                if abs(ia)<=1 and abs(ib)<=1:
                    assert q['status']==0 and int(q['iloc'])==r['station']
                    assert abs(q['omega_r']-r['omega_r']-ia*h)<1e-10
                    assert abs(q['beta_r']-r['beta_r']-ib*h)<1e-10
                    grid[ia,ib]=q
            assert len(grid)==9
            q=grid[0,0]
            delta=max(abs(q[k]-r[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i','growth'])
            assert delta<1e-7
            g=lambda i,j:grid[i,j]['growth']
            gw=(g(1,0)-g(-1,0))/(2*h)
            gb=(g(0,1)-g(0,-1))/(2*h)
            k=(g(0,1)-2*g(0,0)+g(0,-1))/(h*h)
            assert abs(g(0,0))<1e-7 and k<0 and (gw>0)==(kind=='lower_local')
            assert abs(gb)<=max(1e-6,abs(k)*1e-5)
            assert abs(r['verification_delta_omega'])<1e-4 and abs(r['verification_delta_beta'])<1e-5
            assert abs(gw-r['g_omega'])<1e-8 and abs(k-r['g_beta_beta'])<1e-6
            records.append(dict(source=name,station=r['station'],omega=r['omega_r'],beta=r['beta_r'],
                growth=q['growth'],g_omega=gw,g_beta=gb,g_beta_beta=k,
                exported_point_field_delta=delta,summary=str(summary_path),summary_sha256=sha(summary_path),
                field=str(field),field_sha256=sha(field),solver_sha256=solver_hash,
                result_file=r['result_file'],result_sha256=sha(r['result_file'])))
        indices=sorted(r['station'] for r in rows)
        families[name]=dict(count=len(rows),station_range=[min(indices),max(indices)] if indices else None,
            missing_inside_range=[i for i in range(min(indices),max(indices)+1) if i not in indices] if indices else [])
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        candidate_sha256=sha(source),families=families,records=records,
        max_abs_growth=max(abs(r['growth']) for r in records),
        original_input_integrity=snapshot(),
        qualification='Actual local beta-maximum neutral candidates, not a complete beta/spectral or causal certification.')
    (output/'additional_boundaries_field_audit_20261024.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['created','families','max_abs_growth']}))


if __name__=='__main__':main()
