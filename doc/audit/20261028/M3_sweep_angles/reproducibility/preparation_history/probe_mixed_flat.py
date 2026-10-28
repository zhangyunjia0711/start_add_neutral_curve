"""Diagnostic only: test the reason for conservative mixed-frame rejection."""
import json,sys
from pathlib import Path
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_angles_fpp_20261028')
sys.path.insert(0,str(ROOT/'work'))
import physical
from solver import Pencil
from scipy.sparse.linalg import norm
def main():
    s=Path(physical.__file__).read_text().replace(";assert flags[0]==flags[1]", ";self.mixed=flags[0]!=flags[1]")
    scope=dict(physical.__dict__);exec(compile(s,'diagnostic_physical_pair','exec'),scope);pair=scope['Pair']('s75',680,-6.072786702144341,'mixed_flat_diagnostic');rows=[]
    for t,i in [(0.,680),(1.,681)]:
        a=pair.pencil(t);b=Pencil('s75',i,-6.072786702144341)
        errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(a.spatial(-7.5),b.spatial(-7.5))]
        row=dict(station=i,theta=t,mixed_flags=pair.mixed,relative_errors=errors);rows.append(row);print('MIXED_FRAME',row,flush=True)
    (ROOT/'outputs/mixed_flat_diagnostic.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
