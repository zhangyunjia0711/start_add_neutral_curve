"""Normal-grid sensitivity of two earlier unclassified outer-localized roots.
An outer peak is not rejected merely because it is far from the wall.
"""
import concurrent.futures,datetime,importlib.util,json,shutil
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from full_driver import ROOT,sha,snapshot

spec=importlib.util.spec_from_file_location("flat_helpers",ROOT/"work/flat_neutral_normal_validation_20261023/validate_flat.py")
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
AREA=ROOT/"work/outer_grid_audit_20261023"
AREA.mkdir(exist_ok=True)
BIN=ROOT/"outputs/bin/START_LST_wall_local_v10"
OLD=ROOT/"work/local_mode_normal_validation_20261022/inputs"

def mapped_mac(base,other):
    n0=helper.normal(Path(base["input_dir"]),119,base["metadata"])
    n1=helper.normal(Path(other["input_dir"]),119,other["metadata"])
    q0=helper.vectors(Path(base["native_vector_file"]),base["jn"])[base["raw_index"]]
    q1=helper.vectors(Path(other["native_vector_file"]),other["jn"])[other["raw_index"]]
    keep=n0<=n1[-1]+1e-12;n=n0[keep];a=q0[keep];b=CubicSpline(n1,q1,axis=0,extrapolate=False)(n)
    w=np.empty(len(n));w[1:-1]=(n[2:]-n[:-2])/2;w[0]=(n[1]-n[0])/2;w[-1]=(n[-1]-n[-2])/2;w=w[:,None]
    return float(abs(np.sum(w*np.conj(a)*b))**2/(np.sum(w*abs(a)**2)*np.sum(w*abs(b)**2)))

def case(label):
    source=ROOT/"outputs/local_mode_audit.json"
    record=next(r for r in json.loads(source.read_text())["results"] if r["label"]==label)
    old=next(r for r in record["stencil"] if r["requested_omega"]==0)
    point=dict(omega_r=0.,beta_r=-7.,alpha_r=float(old["alpha_r"]),alpha_i=float(old["alpha_i"]))
    results=[]
    for name,inp,jn in [("original400",ROOT.parent/"M1/in",400),
        ("pchip799",OLD/"pchip799",799),("pchip1597",OLD/"pchip1597",1597),
        ("crop340",OLD/"crop340",340),("crop360",OLD/"crop360",360)]:
        r=helper.solve("exact_"+label+"_"+name,inp,jn,119,point,
            here=AREA,binary=BIN,sha=sha(BIN))
        results.append(r)
    base=results[0];comparisons=[]
    for r in results[1:]:
        item=dict(name=r["name"],passed=r["passed"])
        if base["passed"] and r["passed"]:
            item.update(delta_alpha=abs(complex(r["alpha_r"]-base["alpha_r"],r["alpha_i"]-base["alpha_i"])),
                delta_growth=r["growth"]-base["growth"],normal_integral_equal_primitive_MAC=mapped_mac(base,r))
        comparisons.append(item)
    return dict(original_candidate=old,source_file=str(source),source_sha256=sha(source),
        point=point,results=results,comparisons=comparisons)

if __name__=="__main__":
    snapshot()
    previous=ROOT/"outputs/outer_grid_audit_20261023.json"
    archive=ROOT/"outputs/outer_grid_initial_nearest_target_20261023.json"
    if previous.exists() and not archive.exists():shutil.copyfile(previous,archive)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        cases=list(pool.map(case,["local119_stable","local119_growing"]))
    out=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        scope="Two unclassified outer-localized roots, i119 omega=0 beta=-7; fixed original CFD data.",
        solver=str(BIN),solver_sha256=sha(BIN),cases=cases,
        qualification="Nearest targets are not guaranteed cross-grid mode identity. Mapped MAC and growth sensitivity reported; no automatic physical rejection or full-envelope inclusion.",
        neutral_reconverged=False,CFD_grid_independence=False)
    (ROOT/"outputs/outer_grid_audit_20261023.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    snapshot()
    print("OUTER_GRID_AUDIT "+json.dumps([{k:c[k] for k in ["point","comparisons"]} for c in cases]),flush=True)

