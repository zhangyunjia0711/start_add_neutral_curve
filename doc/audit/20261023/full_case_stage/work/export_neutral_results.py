#!/usr/bin/env python3
"""Collect verified wall-local neutral candidates; unknown coverage stays explicit.
This is an evidence export, not a claim of complete spectral enumeration.
"""
import argparse, csv, datetime, json, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from full_driver import ROOT, sha
from centered_survey import coarse_stations

TRACE_PATTERNS = [
    "neutral_x_cont_v7a_*_progress.json",
    "neutral_M3near_x_v7_*_progress.json",
    "neutral_M3near_signed_v7_*_progress.json",
    "neutral_M3nearc_x_v7_*_progress.json",
    "neutral_M1p_upstream_v10_*_progress.json",
    "neutral_M1p_nose_v10_*_progress.json",
    "neutral_M1c_nose_v10_*_progress.json",
    "neutral_tail_beta4_v10_*_progress.json",
    "local_peak_tail_*_hsmall_progress.json",
    "local_peak_nose_*_progress.json",
    "tail841_lower_x_v10_20261024_*_progress.json",
    "tail841_low_lobe_upper_x_v10_20261024_*_progress.json",
    "local_peak_tail841_*_progress.json",
]
EXTRA_PATTERNS = [
    "fd_short_v7_actual_origin_progress.json",
    "M1c_i32_resume_seed.json",
    "native_competition_v10_20261023/results.json",
    "local_mode_audit_M3near_seed_results.json",
    "M3nearc_i119_seed_results_20261023.json",
    "tail_beta_roots_20261023.json",
    "tail841_all_boundaries_20261024.json",
    "dense_tail841_*_20261024/stations/i*.json",
    "tail_crossover_native_i773.json",
    "primary_crossover_native_i372.json",
    "dense_native_upper_v10_20261023/stations/i*.json",
    "dense_native_neighbor_v10_20261023/stations/i*.json",
    "nose_local_candidates_20261023_1616.json",
    "upper_boundary_139_v10_20261023/boundary_new_roots.json",
    "nose_M1c_i11_fine_refine_20261023.json",
    "nose_*_fine_refine_20261023.json",
    "upper_boundary_audit_v10_20261023/i650_competition_refine_v10_actual.json",
    "upper_boundary_audit_v10_20261023/i650_competition_refine_v10_neighbor.json",
]
FIELDS = ["station","x_mm","omega_r","beta_r","alpha_r","alpha_i","growth",
          "g_omega","g_beta","g_beta_beta","frequency_boundary","fd_h",
          "verification_delta_omega","verification_delta_beta","source",
          "solver_sha256","verification_run","result_file","record_number"]
OUT = ROOT/"outputs"

def dumpcsv(path, rows, fields):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def valid_result(result, source, path, index):
    if result.get("status") != "candidate_converged":
        return None
    if result.get("coordinates","wall_local") != "wall_local":
        return None
    v=result.get("verification")
    point=result.get("point")
    if not v or not point or not v.get("J"):
        return None
    # Keep the point and derivatives belonging to the same final verification.
    verified=v.get("center",point)
    if any(abs(verified[k]-point[k])>1e-7 for k in
           ["omega_r","beta_r","alpha_r","alpha_i"]):
        raise ValueError("Point/verification mismatch: "+str(path))
    curvature=v["J"][1][1]
    if abs(point["growth"])>1e-7 or curvature>=0:
        raise ValueError("Not a verified neutral beta maximum: "+str(path))
    gw=v["J"][0][0]
    boundary="upper_local" if gw<0 else "lower_local" if gw>0 else "degenerate"
    correction=result.get("verification_correction",[None,None])
    return dict(point, station=int(round(point.get("iloc",result.get("station")))),
        x_mm=35*point["x"],g_omega=gw,g_beta=v["F"][1],
        g_beta_beta=curvature,frequency_boundary=boundary,fd_h=v["h"],
        verification_delta_omega=correction[0],verification_delta_beta=correction[1],
        source=source,solver_sha256=result.get("solver_sha256"),
        verification_run=v.get("name"),result_file=str(path),record_number=index)

def collect():
    rows=[];inputs=[];failures=[]
    paths=set()
    for pattern in TRACE_PATTERNS+EXTRA_PATTERNS:
        paths.update(OUT.glob(pattern))
    for path in sorted(paths):
        d=json.loads(path.read_text())
        before=len(rows)
        if "records" in d:
            source=d.get("source",d.get("spec",{}).get("seed_provenance",{}).get("source"))
            if not source: continue
            for i,r in enumerate(d["records"]):
                result=dict(r.get("result",{}))
                if not result.get("solver_sha256"):
                    result["solver_sha256"]=d.get("solver_sha256",d.get("spec",{}).get("solver_sha256"))
                row=valid_result(result,source,path,i)
                if row:rows.append(row)
            failures.extend(dict(source=source,result_file=str(path),
                station=r.get("station"),reason=r.get("reason",r.get("status","unknown")))
                for r in d.get("failures",[]))
        elif "results" in d:
            for i,r in enumerate(d["results"]):
                source=r.get("job",{}).get("source","unspecified")
                row=valid_result(r.get("result",{}),source,path,i)
                if row:rows.append(row)
        elif "result" in d:
            source="M2p" if "i650" in path.name else "unspecified"
            row=valid_result(d["result"],source,path,0)
            if row:rows.append(row)
        inputs.append(dict(path=str(path),sha256=sha(path),accepted=len(rows)-before,
                           finished=d.get("finished"),updated=d.get("updated")))
    return rows,inputs,failures

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",default="neutral_variable_beta_candidates")
    a=p.parse_args();base=OUT/a.output
    rows,inputs,failures=collect()
    rows.sort(key=lambda r:(r["station"],r["omega_r"],r["beta_r"],r["source"]))
    dumpcsv(base.with_suffix(".csv"),rows,FIELDS)
    # Strict numerical deduplication keeps the best available FD verification.
    unique=[]
    for r in rows:
        match=next((u for u in reversed(unique[-30:])
            if u["station"]==r["station"]
            and abs(u["omega_r"]-r["omega_r"])<1e-7
            and abs(u["beta_r"]-r["beta_r"])<1e-5
            and abs(u["alpha_r"]-r["alpha_r"])<1e-5),None)
        if match:
            match["provenance"].append({k:r[k] for k in
                ["source","result_file","record_number","verification_run"]})
        else:
            unique.append(dict(r,provenance=[{k:r[k] for k in
                ["source","result_file","record_number","verification_run"]}]))
    grid=json.loads((OUT/"seed_catalog.json").read_text())["grids"]["main_841"]["x_nondim"]
    upper=[];coverage=[]
    coarse=set(coarse_stations(grid))
    for station,x in enumerate(grid,1):
        eligible=[r for r in unique if r["station"]==station
            and r["omega_r"]>=0 and r["frequency_boundary"]=="upper_local"]
        signed=[r for r in unique if r["station"]==station]
        best=max(eligible,key=lambda r:r["omega_r"]) if eligible else None
        if best:upper.append(best)
        coverage.append(dict(station=station,x_mm=35*x,representative=station in coarse,
            status="highest_verified_candidate_audit_pending" if best else
                   "unknown_no_verified_positive_upper_candidate",
            candidate_count=len(signed),positive_upper_count=len(eligible),
            omega_r=best["omega_r"] if best else None,beta_r=best["beta_r"] if best else None,
            source=best["source"] if best else None))
    dumpcsv(base.with_name(base.name+"_highest_known_upper.csv"),upper,FIELDS)
    dumpcsv(base.with_name(base.name+"_coverage_841.csv"),coverage,list(coverage[0]))
    dumpcsv(base.with_name(base.name+"_coverage_139.csv"),
            [r for r in coverage if r["representative"]],list(coverage[0]))
    positive=[r for r in unique if r["omega_r"]>=0]
    fig,axs=plt.subplots(3,1,figsize=(11,10),constrained_layout=True)
    # Plot accepted samples rather than connect uncalculated x stations.
    scatter=axs[0].scatter([r["x_mm"] for r in positive],[r["omega_r"] for r in positive],
        c=[r["beta_r"] for r in positive],cmap="viridis",s=8,alpha=.65)
    axs[0].scatter([r["x_mm"] for r in upper],[r["omega_r"] for r in upper],
        s=12,facecolors="none",edgecolors="#202020",linewidths=.45,
        label="Highest verified positive upper candidate")
    fig.colorbar(scatter,ax=axs[0],label=r"$\beta_*$ (dimensionless)",pad=.02)
    axs[0].set(xlim=(-35,1100),ylabel=r"$\omega$ (dimensionless)",
        title="75-degree blunt plate: variable-beta neutral candidates")
    axs[0].legend(loc="upper right",fontsize=8,frameon=False)
    for boundary,marker,color in [("upper_local","o","#2166ac"),("lower_local","x","#b35806")]:
        q=[r for r in unique if r["x_mm"]<=50 and r["frequency_boundary"]==boundary]
        axs[1].scatter([r["x_mm"] for r in q],[r["omega_r"] for r in q],
            s=10,marker=marker,c=color,label=boundary.replace("_"," "))
    axs[1].set(xlim=(-35,50),ylabel=r"Signed $\omega$ (dimensionless)",
        title="Nose detail; negative-frequency points remain signed")
    axs[1].legend(ncol=2,frameon=False,fontsize=8)
    axs[2].scatter([r["x_mm"] for r in upper],[r["beta_r"] for r in upper],
        s=12,c="#1b7837",label=r"$\beta_*$ of highest verified candidate")
    missing=[r for r in coverage if r["representative"] and not r["positive_upper_count"]]
    axs[2].scatter([r["x_mm"] for r in missing],np.zeros(len(missing)),
        s=16,c="#b2182b",marker="x",label="Representative station: upper candidate unknown")
    axs[2].set(xlim=(-35,1100),ylabel=r"$\beta_*$ (dimensionless)",xlabel="Global wall x (mm)")
    axs[2].legend(frameon=False,fontsize=8)
    for ax in axs:
        ax.grid(alpha=.2);ax.axvline(0,color=".6",linewidth=.7)
    fig.suptitle("Preliminary: multi-path beta-domain audit in progress; missing portions are not filled.",
        fontsize=10,color="#9b3229")
    fig.savefig(base.with_suffix(".png"),dpi=180)
    fig.savefig(base.with_suffix(".svg"))
    plt.close(fig)
    meta=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        qualification="Verified individual local beta maxima at g=0. Not a complete spectral or finite-domain envelope.",
        coordinate_model="Frozen wall-local tangent plane, zero explicit curvature and local parallel approximation.",
        units="x is global wall x in mm; omega and beta are dimensionless.",
        negative_frequency_policy="Actual signed roots retained; no absolute-value or beta-only symmetry substitution.",
        raw_verified_records=len(rows),unique_candidates=len(unique),upper_station_count=len(upper),
        representative_stations=len(coarse),
        representative_upper_count=sum(r["representative"] and r["positive_upper_count"]>0 for r in coverage),
        max_abs_growth=max(abs(r["growth"]) for r in rows),
        inputs=inputs,failures=failures,unique_candidates_data=unique,coverage=coverage)
    base.with_suffix(".json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:meta[k] for k in ["created","raw_verified_records","unique_candidates",
        "upper_station_count","representative_stations","representative_upper_count","max_abs_growth"]}),flush=True)

if __name__=="__main__":
    main()

