import ast,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'work'))
import solver
from scipy.sparse.linalg import norm
def main():
    old=ROOT.parent/'75°/neutral_m3_frequency_limits_20261027/work/solver.py'
    tree=ast.parse(old.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='export')
    scope=dict(solver.__dict__);exec(compile(ast.Module(body=[fn],type_ignores=[]),str(old),'exec'),scope);full=scope['export'];rows=[]
    for case,i in [('s30',119),('s45',241),('s60',400),('s60',2),('s75',119)]:
        a,r=solver.export(case,i,-7.,2.3,label='stencil_validation');b,s=full(case,i,-7.,2.3,label='full_validation')
        errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(a,b)]
        rows.append(dict(case=case,station=i,errors=errors,cropped=r,full=s));print('STENCIL',case,i,errors,flush=True)
        assert max(errors)<1e-11
    (ROOT/'outputs/stencil_validation.json').write_text(json.dumps(dict(passed=True,checks=rows),indent=2)+'\n')
if __name__=='__main__':main()
