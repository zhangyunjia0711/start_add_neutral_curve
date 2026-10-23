"""Densify actual positive upper neutral candidates at native x stations.
All origins are frozen verified roots. Each target requires true-iloc MAC
transport and an independent g=g_beta=0 solve with half-step verification.
The result remains a highest-known candidate, not a complete spectrum.
"""
import argparse, concurrent.futures, datetime, json, shutil, time
from pathlib import Path
from full_driver import ROOT,sha,snapshot
from trace_envelope_x import transfer
from fd_envelope import refine

AREA=ROOT/"outputs/dense_native_upper_v10_20261023"
SOLVER=ROOT/"outputs/bin/START_LST_wall_local_v10"

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def write(path,d):
    temp=path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");temp.replace(path)

def plan():
    AREA.mkdir(exist_ok=True);(AREA/"stations").mkdir(exist_ok=True)
    frozen=AREA/"candidate_inputs.json"
    if not frozen.exists():
        shutil.copyfile(ROOT/"outputs/neutral_variable_beta_candidates.json",frozen)
    source=json.loads(frozen.read_text())
    candidates=[r for r in source["unique_candidates_data"]
        if r["omega_r"]>=0 and r["frequency_boundary"]=="upper_local"]
    upper={}
    for r in candidates:
        station=r["station"]
        if station not in upper or r["omega_r"]>upper[station]["omega_r"]:
            upper[station]=r
    jobs=[]
    for station in range(33,842):
        if station in upper:continue
        left=max((i for i in upper if i<station),default=None)
        right=min((i for i in upper if i>station),default=None)
        origins=[upper[i] for i in (left,right) if i is not None and abs(i-station)<=20]
        origins.sort(key=lambda r:abs(r["station"]-station))
        if len(origins)==2 and abs(origins[0]["beta_r"]-origins[1]["beta_r"])<=1:
            origins=origins[:1]
        jobs.append(dict(station=station,origins=origins))
    data=dict(created=stamp(),solver=str(SOLVER),solver_sha256=sha(SOLVER),
        candidate_input_path=str(frozen),candidate_input_sha256=sha(frozen),
        known_upper_stations=len(upper),jobs=jobs,
        qualification="Actual native-station upper candidates; finite beta/family competition is a separate audit.",
        numerical_policy=dict(minimum_MAC=.95,transport_jump=4,normal_pep_tol=1e-12,
            verification_pep_tol=1e-13,fd_h=.0025,beta_window_halfwidth=.5,
            omega_position_tol=1e-4,beta_position_tol=1e-5,
            no_positive_abs_frequency_substitution=True))
    write(AREA/"plan.json",data);return data

def solve(job):
    station=job["station"];path=AREA/f"stations/i{station:03d}.json"
    if path.exists():
        return json.loads(path.read_text())
    start=time.monotonic();results=[];failures=[]
    for number,old in enumerate(job["origins"]):
        label=f"dense_upper_v10_i{station:03d}_from{old['station']:03d}_p{number}"
        origin={k:old[k] for k in ["omega_r","beta_r","alpha_r","alpha_i"]}
        entry=dict(station=station,source=old["source"],name=label,
            origin_station=old["station"],origin_point=origin,
            origin_result_file=old["result_file"],origin_verification_run=old["verification_run"])
        try:
            actual,transport=transfer(old["station"],station,origin,SOLVER,label,True)
            entry["transport"]=transport
            result=refine(station,actual,actual["omega_r"],actual["beta_r"],SOLVER,
                label=label,max_iterations=25,h=.0025,
                omega_bounds=(0,80),beta_bounds=(old["beta_r"]-.5,old["beta_r"]+.5))
            if result["status"]!="candidate_converged":
                failures.append(dict(job=entry,reason=result.get("reason"),status=result["status"]))
            elif result["verification"]["J"][0][0]>=0:
                failures.append(dict(job=entry,reason="Verified root is not an upper local frequency boundary."))
            results.append(dict(job=entry,result=result))
        except Exception as e:
            failures.append(dict(job=entry,reason=repr(e)))
    accepted=[r["result"]["point"] for r in results
        if r["result"].get("status")=="candidate_converged"
        and r["result"]["verification"]["J"][0][0]<0]
    data=dict(station=station,created=stamp(),finished=True,seconds=time.monotonic()-start,
        status="actual_upper_candidate" if accepted else "unknown",
        highest_known=max(accepted,key=lambda r:r["omega_r"]) if accepted else None,
        results=results,failures=failures,qualification="Separate domain audit pending.")
    write(path,data);return data

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--limit",type=int)
    a=p.parse_args();snapshot();spec=plan();jobs=spec["jobs"]
    if a.limit:jobs=jobs[:a.limit]
    status=dict(created=stamp(),plan_sha256=sha(AREA/"plan.json"),
        total_jobs=len(spec["jobs"]),scheduled_jobs=len(jobs),workers=a.workers,results=[],finished=False)
    write(AREA/"progress.json",status);start=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(solve,j):j for j in jobs}
        for future in concurrent.futures.as_completed(futures):
            job=futures[future]
            try:
                d=future.result()
                short={k:d[k] for k in ["station","status","seconds","highest_known"]}
                short["failures"]=len(d["failures"])
            except Exception as e:
                short=dict(station=job["station"],status="driver_error",reason=repr(e))
            status["results"].append(short)
            status.update(updated=stamp(),elapsed_seconds=time.monotonic()-start)
            write(AREA/"progress.json",status)
            print("DENSE_NATIVE "+json.dumps(short),flush=True)
    status.update(finished=True,updated=stamp(),elapsed_seconds=time.monotonic()-start)
    write(AREA/"progress.json",status);snapshot()
    print("DENSE_FINISHED "+json.dumps(dict(stations=len(status["results"]),
        actual=sum(r["status"]=="actual_upper_candidate" for r in status["results"]),
        elapsed_seconds=status["elapsed_seconds"])),flush=True)

if __name__=="__main__":main()

