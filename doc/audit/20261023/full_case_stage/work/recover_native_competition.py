"""Recover competing lobes at detected apparent upper-candidate jumps.
Having any neutral root at a station is not evidence that its competitors
have been computed. All new points use actual adjacent-root transport.
"""
import concurrent.futures,json,threading,datetime
from full_driver import ROOT,sha,snapshot
from trace_envelope_x import transfer
from fd_envelope import refine
AREA=ROOT/"outputs/native_competition_v10_20261023"
BIN=ROOT/"outputs/bin/START_LST_wall_local_v10"
LOCK=threading.Lock()
CHAINS=[dict(label="M1high_nose_cross",source="M1high_competition",origin=103,targets=list(range(102,84,-1)),beta_range=(0,10)),
        dict(label="M1high_i206",source="M1high_competition",origin=205,targets=[206,207,208],beta_range=(15,25)),
        dict(label="M1high_i319",source="M1high_competition",origin=318,targets=[319,320,321],beta_range=(5,15)),
        dict(label="tail_beta4_i777",source="tail_beta4",origin=776,targets=[777,778],beta_range=(3,5))]

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def write(p,d):
    q=p.with_suffix(p.suffix+".tmp");q.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");q.replace(p)

def chain(spec,candidates):
    path=AREA/(spec["label"]+".json")
    if path.exists() and json.loads(path.read_text()).get("finished"):
        return json.loads(path.read_text())
    choices=[r for r in candidates if r["station"]==spec["origin"]
        and spec["beta_range"][0]<r["beta_r"]<spec["beta_range"][1]]
    assert choices,spec
    old=max(choices,key=lambda r:r["omega_r"]);origin=spec["origin"]
    report=dict(spec=spec,created=stamp(),initial=old,results=[],failures=[],finished=False)
    for i in spec["targets"]:
        label=spec["label"]+f"_to{i:03d}"
        job=dict(station=i,source=spec["source"],name=label,origin_station=origin,origin_point=old)
        try:
            actual,transport=transfer(origin,i,old,BIN,label,True);job["transport_run"]=transport
            result=refine(i,actual,actual["omega_r"],actual["beta_r"],BIN,
                label=label,max_iterations=25,h=.0025,omega_bounds=(0,80),
                beta_bounds=(old["beta_r"]-.5,old["beta_r"]+.5))
            report["results"].append(dict(job=job,result=result))
            if result["status"]!="candidate_converged":
                report["failures"].append(dict(station=i,reason=result.get("reason")));break
            assert result["verification"]["J"][0][0]<0
            old=result["point"];origin=i
        except Exception as e:
            report["failures"].append(dict(station=i,reason=repr(e)));break
        write(path,report)
    report.update(finished=True,updated=stamp());write(path,report)
    print("COMPETING_CHAIN "+json.dumps(dict(label=spec["label"],accepted=len(report["results"]),
        failures=report["failures"])),flush=True)
    return report

if __name__=="__main__":
    snapshot();AREA.mkdir(exist_ok=True)
    source=ROOT/"outputs/neutral_variable_beta_candidates.json"
    candidates=json.loads(source.read_text())["unique_candidates_data"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        outputs=list(pool.map(lambda s:chain(s,candidates),CHAINS))
    final=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),
        solver_sha256=sha(BIN),chains=outputs,
        results=[r for o in outputs for r in o["results"]],
        qualification="Actual competing stationary upper roots recovered; no source label proves permanent mode identity.")
    write(AREA/"results.json",final);snapshot()

