"""Restart a guide from an actually corrected native point away from its tip."""
import json,shutil,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import sha
def main():
    case='rough';f=-112.;tag=f'{case}_f{f}kHz_S';path=ROOT/'outputs'/f'{tag}_v3_guide.json';d=json.loads(path.read_text());assert d['status']=='running'
    choices=[r for r in d['records'] if abs(r['point']['xi']-round(r['point']['xi']))<1e-8 and abs(r['tangent'][0])>.1];assert choices;r=max(choices,key=lambda r:abs(r['tangent'][0]));point=dict(r['point']);qpath=Path(r['vector_file']);point.update(case=case,f_kHz=f,station=int(round(point['xi'])),side='lower',kind='native',vector_file=str(qpath),vector_sha256=sha(qpath),guide_vector_suffix='_native_restart',restart_reason='Initial Poincare closure at an x turning point left the piecewise native cell; use a verified regular native landing')
    seeds=ROOT/'outputs/neutral_seeds.json';s=json.loads(seeds.read_text());s['records']=[point if x['case']==case and x['f_kHz']==f else x for x in s['records']];seeds.write_text(json.dumps(s,indent=2)+'\n')
    previous=path.with_name(path.stem+'_attempt1.json');assert not previous.exists();path.rename(previous);log=path.with_suffix('.log');shutil.copy2(log,log.with_name(log.stem+'_attempt1.log'))
    script=ROOT/'work/contour3d.py';text=script.read_text();old="wave=ROOT/'work'/f'{tag}_v3_guide';wave.mkdir(exist_ok=True)";new="wave=ROOT/'work'/(f'{tag}_v3_guide'+seed.get('source_record',{}).get('guide_vector_suffix',''));wave.mkdir(exist_ok=True)";assert text.count(old)==1;text=text.replace(old,new);script.write_text(text)
    from component_guides import one
    result=one((case,f));(ROOT/'outputs/rough_secondary_112_restart_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
