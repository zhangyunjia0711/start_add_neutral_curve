"""Create native Tecplot view/roundtrip macros, or verify exact numeric read-back."""
import argparse,datetime,hashlib,json,math,re
from pathlib import Path
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_omega0_xbeta_20261026')
NAMES=['M3_omega0_xbeta_comparison','bluntplate_M3_omega0_xbeta','rough_M3_omega0_xbeta']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):
    variables=[];zones=[];invars=False
    for line in path.read_text(encoding='ascii').splitlines():
        s=line.strip()
        if s.startswith('VARIABLES'):invars=True;variables+=re.findall(r'"([^"]+)"',s);continue
        if invars and s.startswith('"'):variables+=re.findall(r'"([^"]+)"',s);continue
        invars=False
        if s.startswith('ZONE '):zones.append(dict(name=re.search(r'T\s*=\s*"([^"]+)"',s).group(1),rows=[],double=False))
        elif s.startswith('DT=') or s.startswith('DT ='):
            assert len(re.findall('DOUBLE',s))==len(variables);zones[-1]['double']=True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])',s):
            row=[float(t) for t in s.split()];assert len(row)==len(variables) and all(math.isfinite(t) for t in row);zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones);return variables,zones

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    if not a.verify:
        full=['#!MC 1410']
        for name in NAMES:
            maps=[('bluntplate',1,'BLUE'),('rough',2,'RED')] if name.startswith('M3_') else [(name.split('_')[0],1,'BLUE' if name.startswith('bluntplate') else 'RED')]
            lines=['#!MC 1410',f'$!READDATASET \'"{name}.plt"\'','  READDATAOPTION = NEW','  RESETSTYLE = YES','  INITIALPLOTTYPE = XYLINE','$!PLOTTYPE = XYLINE','$!DELETELINEMAPS']
            for j,(label,zone,color) in enumerate(maps,1):
                lines+=['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f"  NAME = '{label}'",f'  ASSIGN{{ZONE={zone} XAXISVAR=1 YAXISVAR=2}}',f'  LINES{{SHOW=YES COLOR={color}}}','  SYMBOLS{SHOW=NO}']
            lines+=['$!ACTIVELINEMAPS = ['+','.join(str(j) for j in range(1,len(maps)+1))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"']
            (DEST/(name+'.mcr')).write_text('\n'.join(lines)+'\n',encoding='ascii')
            full += [f'$!INCLUDEMACRO "{name}.mcr"',f'$!WRITEDATASET "{name}_roundtrip.dat"','  BINARY = NO','  USEPOINTFORMAT = YES','  PRECISION = 17','  INCLUDETEXT = NO','  INCLUDEGEOM = NO','  INCLUDECUSTOMLABELS = NO']
        full+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(full)+'\n',encoding='ascii')
        return
    log=(DEST/'batch.log').read_text();assert 'ERROR' not in log.upper();receipts=[]
    for name in NAMES:
        source=DEST/(name+'.dat');back=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt')
        va,za=read(source);vb,zb=read(back);assert va==vb and za==zb,(name,'Native roundtrip mismatch')
        assert binary.read_bytes()[:8]==b'#!TDV112';assert (DEST/(name+'.lay')).exists()
        assert all(r[2]==0 for z in za for r in z['rows'])
        assert all(z['rows'][0][:2]!=z['rows'][-1][:2] for z in za),'Artificial closure'
        receipts.append(dict(file=binary.name,binary_magic='#!TDV112',sha256=sha(binary),source_sha256=sha(source),roundtrip_sha256=sha(back),zones=len(za),points=sum(len(z['rows']) for z in za),double_values=sum(len(z['rows'])*len(va) for z in za),all_values_exact=True,all_omega_zero=True,all_boundaries_open=True))
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),tools='Preplot and Tecplot 360 EX 2024 R1',batch_log_sha256=sha(DEST/'batch.log'),macro_log_error_free=True,files=receipts)
    (DEST/'tecplot_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
