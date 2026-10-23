"""Find actual positive growth witnesses at omega=0; never infer stability
from the absence of a positive sample. Beta guesses are interpolated only as
initial scan parameters; every output growth is independently solved.
"""
import argparse,concurrent.futures,datetime,hashlib,json,time
import numpy as np
from full_driver import ROOT,sha,run,table,snapshot
BIN=ROOT/"outputs/bin/START_LST_wall_local_v10"
AREA=ROOT/"outputs/zero_frequency_witness_v10_20261023"

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def write(p,d):
    q=p.with_suffix(p.suffix+".tmp");q.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");q.replace(p)

def roots():
    records=[]
    for pattern in ["neutral_x_cont_v7a_M1p_*_progress.json","neutral_M1p_*_M1p_upstream_progress.json",
                    "neutral_M3near*_M3near_*progress.json","local_peak_nose_M1p*progress.json"]:
        for p in (ROOT/"outputs").glob(pattern):
            d=json.loads(p.read_text())
            records.extend(r["point"] for r in d.get("records",[]) if r.get("status")=="candidate_converged")
    bystation={}
    for r in records:bystation[int(r["iloc"])]=r
    assert 32 in bystation and 841 in bystation
    return bystation

def solve(job):
    i=job["station"];path=AREA/f"stations/i{i:03d}.json"
    if path.exists():return json.loads(path.read_text())
    origin=job["origin"];beta=job["beta_guess"]
    spec=dict(station=i,mode="M1",seed={k:origin[k] for k in ["omega_r","beta_r","alpha_r","alpha_i"]},
        solver=str(BIN),solver_sha256=sha(BIN),coordinates="wall_local",task=33,
        omega_min=0,omega_max=0,omega_points=1,beta_min=beta-2,beta_max=beta+job['upper_pad'],beta_points=job['points'],
        step=.2,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        minimum_step=1e-6,sample_only=True,anchor_row_seed=True,retry_eigenpairs=8,
        seed_shift=1e-4,eigen_tol=1e-7,
        options=["-pep_tol","1e-12","-pep_max_it","1000","-pep_nev","1","-pep_ncv","8",
                 "-st_pc_factor_mat_solver_type","superlu"],timeout=120)
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec["name"]=f"zero_witness_v10_i{i:03d}_{token}"
    try:
        summary=run(spec)
        raw=table(ROOT/"runs"/spec["name"]/"out/neutral_field.plt")
        usable=[r for r in raw if r["status"]==0 and abs(r["omega_r"])<1e-12]
        okay=not(summary["petsc_errors"] or summary["unused_options"] or summary["coordinate_metadata_missing"])
        positive=[r for r in usable if r["growth"]>1e-7] if okay else []
        best=max(positive,key=lambda r:r["growth"]) if positive else None
        result=dict(station=i,status="positive_growth_witness" if best else "unknown_no_positive_witness",
            created=stamp(),guess=beta,run=spec["name"],summary=summary,samples=raw,witness=best,
            qualification="Positive is evidence G(omega=0)>0 for a followed pair; otherwise unknown, not stable.")
    except Exception as e:
        result=dict(station=i,status="unknown_driver_error",created=stamp(),reason=repr(e),run=spec["name"])
    write(path,result);return result

def main():
    global AREA
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--stations",type=int,nargs="+")
    p.add_argument("--phase",default=AREA.name)
    p.add_argument("--upper-pad",type=float,default=0)
    p.add_argument("--points",type=int,default=3)
    a=p.parse_args();AREA=ROOT/'outputs'/a.phase
    snapshot();AREA.mkdir(exist_ok=True);(AREA/"stations").mkdir(exist_ok=True)
    nominal_path=ROOT/"outputs/axial_wall_local_v10_near_homotopy.json"
    nominal=json.loads(nominal_path.read_text())["records"]
    bystation={int(r["iloc"]):r for r in nominal}
    neutral=roots();keys=sorted(neutral)
    jobs=[];known=[]
    for i in a.stations or range(32,842):
        if bystation[i]["growth"]>1e-7:
            known.append(dict(station=i,status="positive_growth_witness",witness=bystation[i],
                run="inherited_nominal_actual_endpoint",qualification="Actual omega=0 beta=-7 native endpoint."))
        else:
            b=float(np.interp(i,keys,[neutral[k]["beta_r"] for k in keys]))
            jobs.append(dict(station=i,origin=bystation[i],beta_guess=b,upper_pad=a.upper_pad,points=a.points))
    state=dict(created=stamp(),solver_sha256=sha(BIN),nominal_path=str(nominal_path),
        nominal_sha256=sha(nominal_path),scheduled=len(jobs),inherited_positive=len(known),
        finished=False,results=known)
    write(AREA/"progress.json",state);start=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(solve,j) for j in jobs]
        for f in concurrent.futures.as_completed(futures):
            r=f.result()
            state["results"].append({k:r[k] for k in ["station","status","run"]})
            state.update(updated=stamp(),elapsed_seconds=time.monotonic()-start)
            write(AREA/"progress.json",state)
            print("ZERO_WITNESS "+json.dumps({k:r[k] for k in ["station","status","run"]}),flush=True)
    state.update(finished=True,updated=stamp(),elapsed_seconds=time.monotonic()-start)
    write(AREA/"progress.json",state);snapshot()
    print("ZERO_WITNESS_FINISHED "+json.dumps(dict(total=len(state["results"]),
        positive=sum(r["status"]=="positive_growth_witness" for r in state["results"]))),flush=True)

if __name__=="__main__":main()

