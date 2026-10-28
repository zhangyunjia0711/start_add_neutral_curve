"""Inspect every existing 1.1m LST table; do not infer curves from maxima."""
import json,re,hashlib
from pathlib import Path
import numpy as np
ROOT=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
OUT=ROOT/'work/M3_sweep_angles_20261027'
def parse(p):
    zones=[];header=[];current=None
    for line in p.read_text(errors='replace').splitlines():
        if line.lower().startswith('zone'):
            current=dict(title=line,rows=[]);zones.append(current)
        elif current is None:header.append(line)
        else:
            try:row=[float(v.replace('D','E')) for v in line.split()]
            except ValueError:continue
            if len(row)==12:current['rows'].append(row)
    records=[]
    for zone in zones:
        a=np.array(zone.pop('rows'))
        if not len(a):continue
        unique,counts=np.unique(a[:,[9,10,11]],axis=0,return_counts=True)
        zone.update(count=len(a),growth_range=[float(a[:,3].min()),float(a[:,3].max())],beta_per_m_range=[float(a[:,4].min()),float(a[:,4].max())],physical_f_kHz_range=[float(a[:,6].min()),float(a[:,6].max())],x_mm_range=[float(a[:,0].min()*1000),float(a[:,0].max()*1000)],station_count=len(unique),samples_per_station=[int(counts.min()),int(counts.max())],first=a[0].tolist(),last=a[-1].tolist())
        records.append(zone)
    return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),header=header,zones=records)
def main():
    allrows=[]
    for angle in sorted(ROOT.iterdir()):
        fine=angle/'1.1m'
        if not fine.is_dir():continue
        files=sorted(fine.glob('en_eigenValue*.plt'))
        for p in files:
            r=parse(p);allrows.append(r)
            z=r['zones'];target=[v for v in z if re.search(r'fkHz=(0|50|100)"',v['title'])]
            print(str(p.relative_to(ROOT)),len(z),'zones',[(v['title'],v['count'],v['growth_range'],v['samples_per_station']) for v in target],flush=True)
    (OUT/'fine_result_audit.json').write_text(json.dumps(allrows,indent=2)+'\n')
    for r in allrows:
        if any(a in r['path'] for a in ['50°','55°','65°','70°']):
            z=next((v for v in r['zones'] if 'fkHz=50"' in v['title']),None)
            if z:print('EXAMPLE',r['path'],json.dumps(z),flush=True)
if __name__=='__main__':main()
