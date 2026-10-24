"""Project the prior M3 eigenfunction onto the new physical normal coordinates."""
import csv,json
from pathlib import Path
import h5py,numpy as np
from common import ROOT,sha
from profile_interpolation import evaluate,mac
from fixed_beta_pencil import OperatorLinearFamily
from spatial_pencil import solve

def metadata(path,i):
    with path.open() as f:return next(q for q in csv.DictReader(line for line in f if not line.startswith('#')) if int(q['iloc'])==i)

def main():
    info=json.loads((ROOT/'outputs/input_conversion.json').read_text());src=Path(info['old_seed_source']);old=json.loads(src.read_text());oldmode=src.parents[3]/'M3'
    fresh=evaluate(oldmode,119,dict(old['point'],alpha_i=0.),0.,1e-7,'old_M3_omega0_identity')
    state=fresh['point'];q=fresh['q'];native=OperatorLinearFamily().native(info['start_station'])
    oldmeta=metadata(ROOT/'runs'/fresh['run']/'wall_local_profile_metadata.csv',119)
    newmeta=metadata(Path(native.exports[0]['files'][0]['path']).parent.parent/'wall_local_profile_metadata.csv',info['start_station'])
    with h5py.File(oldmode/'in/grid_vertical2d_N.h5') as f:g=f['grid'][:,118,:,0]
    nold=(g-g[0])@np.array([float(oldmeta['n_x']),float(oldmeta['n_y'])])
    with h5py.File(ROOT/'M3/in/grid_vertical2d_N.h5') as f:g=f['grid'][:,info['start_station']-1,:,0]
    nnew=(g-g[0])@np.array([float(newmeta['n_x']),float(newmeta['n_y'])])
    P=native.A[0].tocsr();print('FIRST_ROWS',[(j,P.indices[P.indptr[j]:P.indptr[j+1]].tolist()) for j in range(7)])
    for layout in ['interleaved','variable_major']:
        fields=q.reshape((400,5)) if layout=='interleaved' else q.reshape((5,400)).T
        reconstructed=np.column_stack([np.interp(nnew,nold,fields[:,k].real)+1j*np.interp(nnew,nold,fields[:,k].imag) for k in range(5)])
        pred=reconstructed.ravel() if layout=='interleaved' else reconstructed.T.ravel();pred/=np.linalg.norm(pred)
        try:
            r,v=solve(native,0.,complex(state['alpha_r'],state['alpha_i']),pred,.65)
            path=ROOT/'work'/('rough_physical_seed_'+layout+'.npy');np.save(path,v)
            result=dict(layout=layout,point=r,vector_file=str(path),vector_sha256=sha(path),normal_grid_old=nold.tolist(),normal_grid_new=nnew.tolist(),
                old_component_norms=np.linalg.norm(fields,axis=0).tolist(),old_component_peak_j=np.argmax(abs(fields),axis=0).tolist())
            (ROOT/'outputs'/('physical_seed_probe_'+layout+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print('PHYSICAL_SEED',layout,r)
        except Exception as error:print('PHYSICAL_SEED_FAILED',layout,repr(error))

if __name__=='__main__':main()
