"""Verify omega=0 positive witnesses by following actual M1 high neutral roots.
An earlier nominal fixed-beta x path can reach a different pair at the same
parameters; a negative value on that path is not a stability conclusion.
"""
import concurrent.futures,datetime,hashlib,json,time
from full_driver import ROOT,sha,run,table,snapshot
from trace_envelope_x import transfer
from zero_frequency_witness import roots
BIN=ROOT/"outputs/bin/START_LST_wall_local_v10"
AREA=ROOT/"outputs/zero_frequency_neutral_witness_v10_20261023"

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def write(p,d):
    q=p.with_suffix(p.suffix+".tmp");q.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");q.replace(p)

def solve(job):
    i=job["station"];path=AREA/f"stations/i{i:03d}.json"
    if path.exists():return json.loads(path.read_text())
    old=job["origin"];j=int(old["iloc"])
    try:
        actual,transport=transfer(j,i,old,BIN,f"zero_neutral_i{i:03d}",True)
        assert actual["omega_r"]>0
        spec=dict(station=i,mode="M1",seed={k:actual[k] for k in ["omega_r","beta_r","alpha_r","alpha_i"]},
            solver=str(BIN),solver_sha256=sha(BIN),coordinates="wall_local",task=33,
            omega_min=0,omega_max=actual["omega_r"],omega_points=2,
            beta_min=actual["beta_r"],beta_max=actual["beta_r"]+1e-8,beta_points=2,
            step=.1,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
            minimum_step=1e-6,sample_only=True,anchor_row_seed=True,retry_eigenpairs=8,
            reverse_outer=True,seed_shift=1e-4,eigen_tol=1e-7,
            options=["-pep_tol","1e-12","-pep_max_it","1000","-pep_nev","1","-pep_ncv","8",
                     "-st_pc_factor_mat_solver_type","superlu"],timeout=120)
        token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
        spec["name"]=f"zero_neutral_v10_i{i:03d}_{token}"
        summary=run(spec);raw=table(ROOT/"runs"/spec["name"]/"out/neutral_field.plt")
        usable=[r for r in raw if r["status"]==0 and abs(r["omega_r"])<1e-12]
        okay=not(summary["petsc_errors"] or summary["unused_options"] or summary["coordinate_metadata_missing"])
        positive=[r for r in usable if r["growth"]>1e-7] if okay else []
        result=dict(station=i,created=stamp(),origin_station=j,origin_point=old,
            transport_run=transport,transported_actual=actual,summary=summary,samples=raw,run=spec["name"],
            status="positive_growth_witness" if positive else "unknown_no_positive_witness",
            witness=max(positive,key=lambda r:r["growth"]) if positive else None)
    except Exception as e:
        result=dict(station=i,status="unknown_driver_error",created=stamp(),reason=repr(e))
    result["qualification"]="Actual omega=0 positive pair proves G>0 there only; negative/missing stays unknown."
    write(path,result);return result

if __name__=="__main__":
    snapshot();AREA.mkdir(exist_ok=True);(AREA/"stations").mkdir(exist_ok=True)
    neutral=roots();keys=sorted(k for k in neutral if neutral[k]["omega_r"]>0)
    nominal=ROOT/"outputs/axial_wall_local_v10_near_homotopy.json"
    known={int(r["iloc"]):dict(station=int(r["iloc"]),status="positive_growth_witness",witness=r,
        result_file=str(nominal)) for r in json.loads(nominal.read_text())["records"] if r["growth"]>1e-7}
    for name in ["zero_frequency_witness_v10_20261023","zero_frequency_expanded_v10_20261023"]:
        for p in (ROOT/"outputs"/name/"stations").glob("i*.json"):
            d=json.loads(p.read_text())
            if d["status"]=="positive_growth_witness":known[d["station"]]=dict(d,result_file=str(p))
    jobs=[dict(station=i,origin=neutral[min(keys,key=lambda k:abs(k-i))]) for i in range(32,842) if i not in known]
    state=dict(created=stamp(),finished=False,scheduled=len(jobs),inherited_positive=len(known),
        solver_sha256=sha(BIN),results=[dict(station=i,status="positive_growth_witness",
            result_file=r["result_file"]) for i,r in known.items() if i>=32],
        qualification="No positive omega interval is filled merely because zero and upper endpoint are known.")
    write(AREA/"plan.json",dict(jobs=jobs,known=known,created=stamp()))
    start=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(solve,j) for j in jobs]
        for f in concurrent.futures.as_completed(futures):
            r=f.result();state["results"].append(dict(station=r["station"],status=r["status"],
                result_file=str(AREA/f"stations/i{r['station']:03d}.json")))
            state.update(updated=stamp(),elapsed_seconds=time.monotonic()-start);write(AREA/"progress.json",state)
            print("ZERO_NEUTRAL "+json.dumps({k:r[k] for k in ["station","status"]}),flush=True)
    state.update(finished=True,updated=stamp(),elapsed_seconds=time.monotonic()-start)
    write(AREA/"progress.json",state);snapshot()
    print("ZERO_NEUTRAL_FINISHED "+json.dumps(dict(total=len(state["results"]),
        positive=sum(r["status"]=="positive_growth_witness" for r in state["results"]))),flush=True)

