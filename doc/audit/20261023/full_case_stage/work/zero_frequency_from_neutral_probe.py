"""Follow the actual tail neutral pair to zero frequency at fixed beta."""
import json,hashlib
from full_driver import ROOT,run,table,sha
path=ROOT/"outputs/neutral_x_cont_v7a_M1p_downstream_progress.json"
point=next(r["point"] for r in json.loads(path.read_text())["records"] if r["station"]==841)
binary=ROOT/"outputs/bin/START_LST_wall_local_v10"
spec=dict(station=841,mode="M1",seed={k:point[k] for k in ["omega_r","beta_r","alpha_r","alpha_i"]},
    solver=str(binary),solver_sha256=sha(binary),coordinates="wall_local",task=33,
    omega_min=0,omega_max=point["omega_r"],omega_points=9,
    beta_min=point["beta_r"],beta_max=point["beta_r"]+1e-8,beta_points=2,
    step=.02,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
    minimum_step=1e-6,sample_only=True,anchor_row_seed=True,retry_eigenpairs=8,
    reverse_outer=True,seed_shift=1e-4,eigen_tol=1e-7,
    options=["-pep_tol","1e-12","-pep_max_it","1000","-pep_nev","1","-pep_ncv","8",
             "-st_pc_factor_mat_solver_type","superlu"],timeout=120)
spec["name"]="zero_from_neutral_i841_"+hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
summary=run(spec)
rows=table(ROOT/"runs"/spec["name"]/"out/neutral_field.plt")
out=dict(source_path=str(path),source_sha256=sha(path),point=point,summary=summary,rows=rows,
    qualification="Actual fixed-beta frequency path from an independently verified neutral root.")
(ROOT/"outputs/zero_from_neutral_i841_probe.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
print("ZERO_FROM_NEUTRAL",[(r["omega_r"],r["growth"],r["status"]) for r in rows],flush=True)

