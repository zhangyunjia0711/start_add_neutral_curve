"""Compare source ASCII with native Tecplot binary read-back, value for value."""
import argparse,datetime,hashlib,json,math,re
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def read(path):
    zones=[];variables=[];in_variables=False
    for line in path.read_text(encoding='ascii').splitlines():
        s=line.strip()
        if s.startswith('VARIABLES'):
            in_variables=True;variables+=re.findall(r'"([^"]+)"',s);continue
        if in_variables and s.startswith('"'):
            variables+=re.findall(r'"([^"]+)"',s);continue
        in_variables=False
        if s.startswith('ZONE '):
            zones.append(dict(name=re.search(r'T\s*=\s*"([^"]+)"',s).group(1),rows=[],double=False))
        elif s.startswith('DT=') or s.startswith('DT ='):
            assert len(re.findall('DOUBLE',s))==len(variables);zones[-1]['double']=True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])',s):
            q=[float(t) for t in s.split()];assert len(q)==len(variables)
            assert all(math.isfinite(x) for x in q);zones[-1]['rows'].append(q)
    assert zones and all(z['double'] for z in zones)
    return variables,zones

def main():
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--roundtrip-dir',type=Path,required=True);a=p.parse_args()
    manifest=a.directory/'tecplot_export_manifest.json';m=json.loads(manifest.read_text());reports=[]
    for item in m['datasets']:
        source=a.directory/item['file'];binary=a.directory/item['binary_file']
        native=a.roundtrip_dir/(source.stem+'_roundtrip.dat')
        original_vars,original=read(source);actual_vars,actual=read(native)
        assert original_vars==actual_vars==item['variables']
        assert len(original)==len(actual)==len(item['zones'])
        for x,y,expected in zip(original,actual,item['zones']):
            assert x['name']==y['name']==expected['name']
            assert len(x['rows'])==len(y['rows'])==expected['points']
            assert x['rows']==y['rows'],'Binary conversion changed a value'
        magic=binary.read_bytes()[:8].decode('ascii');assert magic.startswith('#!TDV')
        item.update(binary_sha256=sha(binary),binary_bytes=binary.stat().st_size,binary_magic=magic,
            native_roundtrip_sha256=sha(native),native_roundtrip_exact=True)
        reports.append(dict(file=binary.name,zones=len(actual),points=sum(len(z['rows']) for z in actual),
            values=sum(len(z['rows'])*len(actual_vars) for z in actual),exact=True))
    for f in m['input_files']:assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256']
    m['native_validation']=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        tool='Tecplot 360 EX 2024 R1 and Preplot 2024.1.0.52134',results=reports,
        all_values_exact=True,source_inputs_unchanged=True)
    manifest.write_text(json.dumps(m,indent=2,allow_nan=False)+'\n')
    print(json.dumps(m['native_validation'],ensure_ascii=False))

if __name__=='__main__':main()
