"""Remove exact, owned, superseded attempts after preserving hashes and conclusions."""
import datetime,json,shutil
from pathlib import Path
from common import ROOT,sha

def main():
    info=json.loads((ROOT/'outputs/input_conversion.json').read_text());delivery=Path(info['source_case'])/'outputs/M3_beta_minus7_20261024'
    output=ROOT/'outputs/cleanup_receipt.json';assert not output.exists()
    candidates=[(ROOT/'work'/n,'Superseded j-index operator-homotopy attempt; conclusions and raw JSON/logs retained') for n in
        ['transfer_m3_seed.py','transfer_spatial_seed.py','finish_spatial_seed.py','inspect_plot3d_source.py','seed_homotopy_vectors','spatial_seed_homotopy_vectors','spatial_seed_v2_vectors','rough_frequency_negative','rough_frequency_positive','rough_M3_spatial_seed.npy']]
    candidates += [(delivery/n,'Completed native roundtrip comparison or empty native output; SHA preserved') for n in ['rough_roundtrip.dat','tecplot_stdout.log','tecplot_stderr.log']]
    removed=[]
    for p,reason in candidates:
        if not p.exists():continue
        resolved=p.resolve();allowed=ROOT.resolve() if str(resolved).startswith(str(ROOT.resolve())+'/') else delivery.resolve()
        resolved.relative_to(allowed);assert not p.is_symlink()
        files=sorted(q for q in p.rglob('*') if q.is_file()) if p.is_dir() else [p]
        if p.is_dir():assert all(q.suffix=='.npy' and not q.is_symlink() for q in files)
        for q in files:
            q.resolve().relative_to(allowed)
            removed.append(dict(path=str(q),bytes=q.stat().st_size,sha256=sha(q),reason=reason))
        if p.is_dir():shutil.rmtree(p)
        else:p.unlink()
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),removed=removed,
        removed_files=len(removed),removed_bytes=sum(r['bytes'] for r in removed),
        retained='All valid physical M3 seed/frequency evidence, v2/v3 vectors referenced by the closed contour, native roots/vectors, physical profiles/matrices, source files, original solver runs and logs.')
    output.write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(output,delivery/output.name)
    print(json.dumps({k:report[k] for k in ['removed_files','removed_bytes']}))

if __name__=='__main__':main()
