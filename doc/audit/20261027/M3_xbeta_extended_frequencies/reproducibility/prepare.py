"""Isolated extension of signed physical frequencies; frozen equations unchanged."""
import datetime,hashlib,json,shutil
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SRC=BASE/'neutral_xbeta_extra_frequencies_20261027'
ROOT=BASE/'neutral_xbeta_extended_frequencies_20261027'
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
FS=[-200,-150,-125,125,150,200]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not ROOT.exists(),'Do not overwrite an earlier run'
    for d in ['work','outputs','runs','matrix_exports','inputs']:(ROOT/d).mkdir(parents=True,exist_ok=True)
    for d in ['native_pencils','cases']:(ROOT/d).symlink_to(SRC/d,target_is_directory=True)
    sources=[]
    for name in ['common.py','run_adapter.py','core.py','solver.py','contour.py','contour3d.py','physical.py','finalize.py','physical_fold_bordered.py','tip_correction.py']:
        p=SRC/'work'/name;s=p.read_text();sources.append(dict(path=str(p),sha256=sha(p)))
        if name in ['contour3d.py','finalize.py']:
            assert 'choices=[-75,-25,25,75]' in s
            s=s.replace('choices=[-75,-25,25,75]','choices='+repr(FS))
        if name=='contour3d.py':
            anchor="        section=float(np.dot(z-z0,t0));distance=float(np.linalg.norm(z-z0))\n"
            assert s.count(anchor)==1
            s=s.replace(anchor,anchor+"        if half==0 and len(state['records'])>30 and distance<1e-7 and mac(q0,q)>1-1e-7:\n            state['closure']=dict(passed=True,scaled_distance=distance,raw_MAC=mac(q0,q),method='Direct corrected native-section seed return');state.update(status='closed_verified',finished=True);break\n")
        if name=='tip_correction.py':s=s.rstrip()+'\n'
        (ROOT/'work'/name).write_text(s)
    protected=json.loads((SRC/'outputs/input_manifest.json').read_text())
    protected += [dict(path=str(p),sha256=sha(p),role='prior_nine_frequency_curve') for p in sorted((WIN/'outputs/M3_xbeta_nine_frequencies_20261027').glob('*.csv'))]
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    shutil.copy2(SRC/'outputs/cases.json',ROOT/'outputs/cases.json')
    c=json.loads((SRC/'outputs/frequency_conversion.json').read_text());c['frequencies']=[dict(f_kHz=f,omega_nd=6.0727867021443425*f/50) for f in FS]
    (ROOT/'outputs/frequency_conversion.json').write_text(json.dumps(c,indent=2)+'\n')
    (ROOT/'outputs/bootstrap_receipt.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),root=str(ROOT),sources=sources,protected_files=len(protected),frequencies_kHz=FS,model='Frozen local wall basis, zero explicit curvature; unchanged from prior nine frequencies'),indent=2)+'\n')
    for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,ROOT/'work'/p.name)
    print('PREPARED',ROOT,len(protected),flush=True)
if __name__=='__main__':main()
