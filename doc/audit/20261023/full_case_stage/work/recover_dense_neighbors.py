"""Recover failed coarse-origin densification through actual adjacent roots.
No growth values are interpolated. A new accepted native endpoint is mandatory.
"""
import argparse, concurrent.futures, datetime, json, time
from full_driver import ROOT,sha
from export_neutral_results import collect
from trace_envelope_x import transfer
from fd_envelope import refine

AREA=ROOT/"outputs/dense_native_neighbor_v10_20261023"
FIRST=ROOT/"outputs/dense_native_upper_v10_20261023"
SOLVER=ROOT/"outputs/bin/START_LST_wall_local_v10"

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def write(p,d):
    q=p.with_suffix(p.suffix+".tmp")
    q.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");q.replace(p)

def solve(job):
    station=job["station"];results=[];failures=[]
    for old in job["origins"]:
        label=f"dense_neighbor_v10_i{station:03d}_from{old['station']:03d}"
        origin={k:old[k] for k in ["omega_r","beta_r","alpha_r","alpha_i"]}
        entry=dict(station=station,source=old["source"],name=label,
            origin_station=old["station"],origin_point=origin,
            origin_result_file=old["result_file"],origin_verification_run=old["verification_run"])
        try:
            actual,transport=transfer(old["station"],station,origin,SOLVER,label,True)
            entry["transport"]=transport
            result=refine(station,actual,actual["omega_r"],actual["beta_r"],SOLVER,
                label=label,max_iterations=25,h=.0025,omega_bounds=(0,80),
                beta_bounds=(old["beta_r"]-.5,old["beta_r"]+.5))
            results.append(dict(job=entry,result=result))
            if result["status"]=="candidate_converged" and result["verification"]["J"][0][0]<0:
                break
            failures.append(dict(job=entry,reason=result.get("reason"),status=result["status"]))
        except Exception as e:
            failures.append(dict(job=entry,reason=repr(e)))
    accepted=[r["result"]["point"] for r in results
        if r["result"].get("status")=="candidate_converged"
        and r["result"]["verification"]["J"][0][0]<0]
    data=dict(station=station,created=stamp(),finished=True,
        status="actual_upper_candidate" if accepted else "unknown",
        highest_known=max(accepted,key=lambda r:r["omega_r"]) if accepted else None,
        results=results,failures=failures,qualification="Actual native endpoint; finite domain competition pending.")
    path=AREA/f"stations/i{station:03d}_wave{job['wave']}.json"
    write(path,data)
    return dict(station=station,status=data["status"],path=str(path),
        origin_stations=[r["station"] for r in job["origins"]],failures=len(failures))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--max-waves",type=int,default=10)
    a=p.parse_args();AREA.mkdir(exist_ok=True);(AREA/"stations").mkdir(exist_ok=True)
    progress_path=AREA/"progress.json"
    state=json.loads(progress_path.read_text()) if progress_path.exists() else dict(
        created=stamp(),solver_sha256=sha(SOLVER),waves=[],attempted_origins={})
    start=time.monotonic()
    for unused in range(a.max_waves):
        rows,inputs,failures=collect()
        upper={}
        for r in rows:
            if r["omega_r"]>=0 and r["frequency_boundary"]=="upper_local":
                if r["station"] not in upper or r["omega_r"]>upper[r["station"]]["omega_r"]:
                    upper[r["station"]]=r
        unknown=[i for i in range(11,33) if i not in upper]
        for path in FIRST.glob("stations/i*.json"):
            d=json.loads(path.read_text())
            if d["status"]=="unknown" and d["station"] not in upper:
                if d["station"] not in unknown:unknown.append(d["station"])
        wave=len(state["waves"])+1;jobs=[]
        for station in sorted(unknown):
            tried=state["attempted_origins"].get(str(station),[])
            origins=[upper[i] for i in (station-1,station+1) if i in upper and i not in tried]
            origins.sort(key=lambda r:r["omega_r"],reverse=True)
            if origins:jobs.append(dict(station=station,origins=origins,wave=wave))
        if not jobs:break
        write(AREA/f"wave{wave:02d}_plan.json",dict(created=stamp(),jobs=jobs,
            inputs=inputs,unknown_before=unknown,
            reason="Coarse-origin beta window can exclude the root. An actual adjacent root shortens transport."))
        wave_result=dict(wave=wave,jobs=len(jobs),results=[])
        state["waves"].append(wave_result);state["finished"]=False
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
            future_jobs={pool.submit(solve,j):j for j in jobs}
            for future in concurrent.futures.as_completed(future_jobs):
                j=future_jobs[future]
                try:r=future.result()
                except Exception as e:r=dict(station=j["station"],status="driver_error",reason=repr(e))
                wave_result["results"].append(r)
                state["attempted_origins"].setdefault(str(j["station"]),[]).extend(
                    old["station"] for old in j["origins"])
                state.update(updated=stamp(),elapsed_seconds=time.monotonic()-start)
                write(progress_path,state);print("DENSE_NEIGHBOR "+json.dumps(r),flush=True)
    state.update(finished=True,updated=stamp(),elapsed_seconds=time.monotonic()-start)
    write(progress_path,state)
    print("DENSE_NEIGHBOR_FINISHED "+json.dumps(dict(waves=len(state["waves"]),
        actual=sum(r["status"]=="actual_upper_candidate" for w in state["waves"] for r in w["results"]))),flush=True)

if __name__=="__main__":main()

