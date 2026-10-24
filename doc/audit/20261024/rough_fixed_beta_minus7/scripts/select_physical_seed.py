"""Freeze the spatially aligned M3 seed and stop its superseded index-based survey."""
import datetime,json,os,signal
from pathlib import Path
from common import ROOT,sha,snapshot
from profile_interpolation import evaluate,mac
import numpy as np

def main():
    stopped=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            args=(p/'cmdline').read_bytes().split(b'\0')
            if b'work/survey_rough_seed.py' in args and (p/'cwd').resolve()==ROOT.resolve():
                os.kill(int(p.name),signal.SIGTERM);stopped.append(int(p.name))
        except (FileNotFoundError,PermissionError,ProcessLookupError):pass
    source=ROOT/'outputs/physical_seed_probe_interleaved.json';d=json.loads(source.read_text())
    info=json.loads((ROOT/'outputs/input_conversion.json').read_text());station=info['start_station']
    assert d['point']['raw_MAC']>.999
    fresh=evaluate('M3',station,d['point'],0.,1e-7,'rough_physical_M3_omega0_independent')
    overlap=mac(np.load(d['vector_file']),fresh['q']);assert overlap>1-1e-7
    assert abs(complex(fresh['point']['alpha_r']-d['point']['alpha_r'],fresh['point']['alpha_i']-d['point']['alpha_i']))<1e-6
    vector=ROOT/'work/rough_physical_M3_seed.npy';np.save(vector,fresh['q'])
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),finished=True,
        station=station,point=dict(fresh['point'],xi=float(station),x_mm=info['start_x_mm']),vector_file=str(vector),vector_sha256=sha(vector),
        physical_normal_projection_source=str(source),source_sha256=sha(source),physical_shape_MAC=d['point']['raw_MAC'],
        independent_spatial_run=fresh['run'],independent_raw_MAC=overlap,
        superseded_survey_stopped_pids=stopped,
        qualification='M3 correspondence is supported by matching actual wall-normal eigenfunction positions and nearby eigenvalues, not by matching j-index on different normal grids. The index-homotopy acoustic-like branch is excluded from the delivered M3 curve.',original_inputs=snapshot())
    (ROOT/'outputs/rough_physical_M3_seed.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['point']))

if __name__=='__main__':main()
