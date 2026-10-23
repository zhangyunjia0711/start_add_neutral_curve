#!/usr/bin/env python3
"""Refine both observed tail peaks independently to compare neutral frequencies."""
import concurrent.futures,datetime,json
from full_driver import ROOT,sha,snapshot
from fd_envelope import refine

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def main():
    snapshot();source=ROOT/'outputs/tail_beta_competition_20261023.json'
    data=json.loads(source.read_text());solver=ROOT/'outputs/bin/START_LST_wall_local_v7'
    assert sha(solver)=='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
    jobs=[]
    for station in [783,788,793]:
        curve=next(v for v in data['curves'] if v['station']==station and v['origin_label']=='old_beta6')
        for beta in [4.4,6]:
            actual=min(curve['local_sampled_peaks'],key=lambda v:abs(v['beta_r']-beta))
            jobs.append(dict(station=station,source='tail_beta4' if beta==4.4 else 'tail_beta6',
                name=f'tail_peak_fd_i{station}_b{4 if beta==4.4 else 6}_20261023',actual=actual))
    output=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),
        solver=str(solver),solver_sha256=sha(solver),results=[],
        qualification='Independent neutral refinements from both actual finite-window beta peaks; unresolved roots retained and no outer-envelope inference from growth alone')
    target=ROOT/'outputs/tail_beta_roots_20261023.json'
    def worker(job):
        p=job['actual']
        result=refine(job['station'],p,p['omega_r'],p['beta_r'],solver,
            coordinates='wall_local',label=job['name'],omega_bounds=(3,7),beta_bounds=(3,8))
        return dict(job=job,result=result)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(worker,jobs):
            output['results'].append(result);output['updated']=stamp()
            tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');tmp.replace(target)
            print(json.dumps(dict(station=result['job']['station'],source=result['job']['source'],
                status=result['result']['status'],point=result['result'].get('point'),reason=result['result'].get('reason'))),flush=True)
    snapshot()

if __name__=='__main__':main()
