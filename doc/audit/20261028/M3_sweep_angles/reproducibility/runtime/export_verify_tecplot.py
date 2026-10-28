"""Native PREPLOT export and independent exact-value binary read-back."""
import csv,datetime,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from read_plt import Reader
DEST=Path(sys.argv[1]) if len(sys.argv)>1 else Path(r'E:\work\calc\Instability of attachment line\Han2019_SweptBluntPlate\outputs\M3_sweep_angles_20261028')
PREPLOT=Path(r'E:\Program Files\Tecplot\Tecplot 360 EX 2024 R1\bin\preplot.exe')
VARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC','sweep_deg','normal_method']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    selection=json.loads((DEST/'selection_manifest.json').read_text(encoding='utf8'));records=selection['curves'];checks=[]
    names=['M3_sweep_angles_all_neutral_curves']+[f'M3_{a}deg_xbeta' for a in [30,45,50,55,60,65,70,73,75]]
    with (DEST/'tecplot_export.log').open('w',encoding='utf8') as log:
        for name in names:
            dat=DEST/(name+'.dat');plt=DEST/(name+'.plt')
            if not dat.exists():continue
            command=[str(PREPLOT),str(dat),str(plt)];r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT);assert r.returncode==0
            expected=records if name==names[0] else [q for q in records if q['angle']==int(name.split('_')[1][:-3])]
            binary=Reader(plt.read_bytes()).data();assert binary.variables==VARS and len(binary.zones)==len(expected)
            values=0
            for zone,rec in zip(binary.zones,expected):
                with (DEST/rec['output_csv']).open(newline='',encoding='utf8') as f:
                    reader=csv.DictReader(f);keys=reader.fieldnames;a=np.array([[float(row[k]) for k in keys] for row in reader])
                if rec['angle']==75:a=np.column_stack([a,np.full(len(a),75.),np.zeros(len(a))])
                b=np.column_stack(zone['arrays']);assert zone['name']==rec['zone'] and b.shape==a.shape and np.array_equal(a,b),(zone['name'],rec['zone'],a.shape,b.shape)
                values+=a.size
            checks.append(dict(file=plt.name,sha256=sha(plt),zones=len(expected),values=values,exact_double_readback=True,bytes=plt.stat().st_size,command=command))
    receipt=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),passed=True,exporter='Installed native Tecplot PREPLOT 2024.1',reader='Independent restricted ordered-zone reader, installed Tecplot Binary Data File Format guide',checks=checks)
    (DEST/'tecplot_validation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8');print('EXACT_TECPLOT_READBACK',[(q['file'],q['zones'],q['values']) for q in checks])
if __name__=='__main__':main()
