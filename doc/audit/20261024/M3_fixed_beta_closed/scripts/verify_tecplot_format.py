"""Verify native Tecplot conversion against the original double-precision data."""
import argparse,datetime,hashlib,json,math,re
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    variables=[];zones=[];invars=False
    for line in p.read_text(encoding='ascii').splitlines():
        s=line.strip()
        if s.startswith('VARIABLES'):
            invars=True;variables+=re.findall(r'"([^"]+)"',s);continue
        if invars and s.startswith('"'):
            variables+=re.findall(r'"([^"]+)"',s);continue
        invars=False
        if s.startswith('ZONE '):
            zones.append(dict(name=re.search(r'T\s*=\s*"([^"]+)"',s).group(1),rows=[],double=False))
        elif s.startswith('DT=') or s.startswith('DT ='):
            assert len(re.findall('DOUBLE',s))==len(variables);zones[-1]['double']=True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])',s):
            row=[float(t) for t in s.split()];assert len(row)==len(variables) and all(math.isfinite(t) for t in row)
            zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones)
    return variables,zones

def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();d=a.directory
    source=d/'M3_closed_beta_minus7.dat';binary=d/'M3_closed_beta_minus7.plt';back=d/'M3_roundtrip.dat'
    va,za=read(source);vb,zb=read(back);assert va==vb and za==zb
    assert len(va)==9 and len(za)==1 and len(za[0]['rows'])==1571
    data=za[0]['rows'];assert all(r[2]==-7 for r in data)
    assert data[0][:5]==data[-1][:5] and data[0][6]==0 and data[-1][6]==2
    assert binary.read_bytes()[:8]==b'#!TDV112'
    log=d/'batch.log';text=log.read_text()
    assert 'ERROR' not in text.upper() and 'ACTION: WriteDataSet: M3_roundtrip.dat' in text
    assert (d/'M3_closed_beta_minus7.lay').is_file()
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        tools='Tecplot 360 EX 2024 R1; Preplot 2024.1.0.52134',binary_file=binary.name,
        binary_magic='#!TDV112',binary_bytes=binary.stat().st_size,binary_sha256=sha(binary),
        source_sha256=sha(source),roundtrip_sha256=sha(back),native_batch_log_sha256=sha(log),
        variables=va,zones=1,points=1571,double_values_exactly_preserved=14139,
        all_numeric_values_exact=True,closed_coordinates_exact=True,all_beta_strictly_minus7=True,
        macro_log_error_free=True,layout_generated=True)
    (d/'tecplot_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
