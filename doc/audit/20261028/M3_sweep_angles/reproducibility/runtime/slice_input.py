"""Export using a bounded native stencil, verified against full-field exports.

Independent frozen START checks still read the original complete angle dataset.
"""
import configparser,json
from pathlib import Path
import h5py,numpy as np
from common import ROOT,CASES,sha
def template(case,station):
    source=ROOT/'cases'/str(case)
    if str(case) not in CASES:return source,station,None
    count=CASES[case]['stations'];start=max(1,min(station-4,count-8));stop=start+9;local=station-start+1
    dest=ROOT/'inputs/native_stencils'/str(case)/f'i{station:03d}';receipt=dest/'slice.json'
    if receipt.exists():return dest,local,json.loads(receipt.read_text())
    (dest/'in').mkdir(parents=True,exist_ok=True)
    files=[]
    for name,key in [('grid_vertical2d_N.h5','grid'),('flow_vertical2d_N.h5','baseflow')]:
        with h5py.File(source/'in'/name) as h:a=h[key][:,start-1:stop-1,:,:]
        path=dest/'in'/name
        with h5py.File(path,'w') as h:h.create_dataset(key,data=a)
        with h5py.File(path) as h:assert np.array_equal(a,h[key][...])
        files.append(dict(path=str(path),sha256=sha(path)))
    cfg=configparser.ConfigParser();cfg.read(source/'config.ini');cfg['domain'].update({'in':'9','iloc':str(local),'istart':'1','iend':'9'})
    with (dest/'config.ini').open('w') as f:cfg.write(f)
    r=dict(case=str(case),global_station=station,local_station=local,native_stencil_first=start,native_stencil_last=stop-1,source=str(source),files=files)
    receipt.write_text(json.dumps(r,indent=2)+'\n');return dest,local,r
