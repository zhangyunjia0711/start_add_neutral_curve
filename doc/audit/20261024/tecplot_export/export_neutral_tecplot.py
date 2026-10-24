"""Export unchanged neutral-curve evidence as Tecplot double-precision XY data."""
import argparse,csv,datetime,hashlib,json,math,re
from pathlib import Path

VARIABLES=['X_mm','Omega','Beta','Alpha_r','Alpha_i','Growth','Station','Evidence','Boundary','X_plus_35_mm']
GUIDE='https://tecplot.azureedge.net/products/360/2024r1m1/360-data-format.html'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quote(t):return '"'+str(t).replace('\\','\\\\').replace('"','\\"')+'"'
def point(r,evidence=0):
    return dict(X_mm=float(r.get('x_mm',35*float(r.get('x',0)))),Omega=float(r['omega_r']),
        Beta=float(r['beta_r']),Alpha_r=float(r['alpha_r']),Alpha_i=float(r['alpha_i']),
        Growth=float(r['growth']),Station=float(r.get('station',r.get('iloc'))),Evidence=float(evidence),
        Boundary=float({'upper_local':1,'lower_local':-1}.get(r.get('frequency_boundary'),0)),
        X_plus_35_mm=float(r.get('x_mm',35*float(r.get('x',0))))+35)

def zone(name,rows,role='neutral_curve',style='line'):
    assert rows
    return dict(name=name,rows=rows,role=role,style=style)

def write_dataset(path,variables,zones,source_hash):
    with path.open('w',encoding='ascii',newline='\n') as f:
        f.write('TITLE="75-degree blunt plate: verified local neutral candidates"\n')
        f.write('VARIABLES='+ ' '.join(quote(v) for v in variables)+'\n')
        for key,value in dict(Status='Stage result; nose stations 1-10 unknown; not a certified full-spectrum envelope',
            CoordinateModel='Frozen wall-local parallel model; zero explicit curvature',
            Units='X_mm in millimetres; Omega and Beta dimensionless',
            EvidenceCodes='0 actual neutral; 1 conjugate-derived neutral; 2 positive growth NOT neutral',
            BoundaryCodes='1 local upper; -1 local lower; 0 not classified',SourceSHA256=source_hash).items():
            f.write('DATASETAUXDATA '+key+'='+quote(value)+'\n')
        for z in zones:
            f.write('ZONE T='+quote(z['name'])+', I='+str(len(z['rows']))+', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT\n')
            f.write('DT=('+ ' '.join(['DOUBLE']*len(variables))+')\n')
            f.write('AUXDATA Role='+quote(z['role'])+'\nAUXDATA RecommendedStyle='+quote(z['style'])+'\n')
            for r in z['rows']:
                values=[r[v] for v in variables];assert all(math.isfinite(x) for x in values)
                f.write(' '.join(format(x,'.17e') for x in values)+'\n')
    # Independent lexical read-back of the written data, with exact float comparison.
    blocks=[];current=None
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('ZONE '):
            current=[];blocks.append(current)
        elif current is not None and line[:1] in '+-.0123456789':
            values=[float(x) for x in line.split()];assert len(values)==len(variables);current.append(values)
    expected=[[[r[v] for v in variables] for r in z['rows']] for z in zones]
    assert blocks==expected,'Export changed a numeric value'
    return dict(file=path.name,binary_file=path.with_suffix('.plt').name,variables=variables,
        zones=[dict(name=z['name'],points=len(z['rows']),role=z['role'],style=z['style']) for z in zones],
        total_points=sum(len(z['rows']) for z in zones),ascii_sha256=sha(path),numeric_roundtrip_exact=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    candidate=a.source_dir/'neutral_variable_beta_candidates.json';status_path=a.source_dir/'case_stage_status.json'
    d=json.loads(candidate.read_text());status=json.loads(status_path.read_text())
    assert status['candidate_input_sha256']==sha(candidate)
    rows=d['unique_candidates_data'];by={}
    for r in rows:
        if r['omega_r']>=0 and r['frequency_boundary']=='upper_local':
            if r['station'] not in by or r['omega_r']>by[r['station']]['omega_r']:by[r['station']]=r
    upper=[by[i] for i in sorted(by)];assert len(upper)==status['native_upper_count']
    lower=sorted([r for r in rows if r['source']=='tail841_lobe_lower'],key=lambda r:r['station'])
    lowupper=sorted([r for r in rows if r['source']=='tail841_low_lobe_upper'],key=lambda r:r['station'])
    derived_path=a.source_dir/'neutral_conjugate_derived.csv'
    derived=[r for r in csv.DictReader(derived_path.open()) if float(r['omega_r'])>=0]
    core=[zone('Highest_verified_positive_upper',[point(r) for r in upper]),
        zone('High_frequency_lobe_local_lower',[point(r) for r in lower]),
        zone('Low_frequency_lobe_local_upper',[point(r) for r in lowupper])]
    mainzones=core+[zone('Conjugate_derived_lower_candidates',[point(r,1) for r in derived],role='conjugate_derived_not_an_extra_solve',style='scatter'),
        zone('Positive_growth_test_i10_NOT_NEUTRAL',[point(status['native10_positive_cap_test']['point'],2)],role='not_neutral_positive_growth_test',style='scatter')]
    groups={}
    for r in rows:
        if r['omega_r']>=0:groups.setdefault((r['source'],r['frequency_boundary']),[]).append(point(r))
    mainzones.extend(zone('Other_candidates_'+source+'_'+kind,values,role='actual_local_neutral_candidates',style='scatter')
                     for (source,kind),values in sorted(groups.items()))
    datasets=[('neutral_variable_beta_stage',VARIABLES,mainzones)]
    nose_upper=[point(r) for r in upper if r['x_mm']<-30]
    nose_derived=[point(r,1) for r in derived if float(r['x_mm'])<-30]
    nosevars=['X_plus_35_mm','Omega']+[v for v in VARIABLES if v not in ['X_plus_35_mm','Omega']]
    datasets.append(('neutral_variable_beta_nose',nosevars,[zone('Nose_highest_verified_upper',nose_upper),
        zone('Nose_conjugate_derived_lower',nose_derived,role='conjugate_derived_not_an_extra_solve',style='scatter')]))
    datasets.append(('neutral_variable_beta_tail',VARIABLES,[zone(z['name'],[r for r in z['rows'] if r['X_mm']>=200]) for z in core]))
    segments=[];part=[]
    switch_starts={v['stations'][1] for v in status['candidate_winner_switches']}
    for r in upper:
        if part and (r['station'] in switch_starts or r['station']!=part[-1]['Station']+1):
            segments.append(part);part=[]
        part.append(point(r))
    if part:segments.append(part)
    betavars=['X_mm','Beta']+[v for v in VARIABLES if v not in ['X_mm','Beta']]
    datasets.append(('neutral_variable_beta_beta',betavars,[zone('Beta_winner_segment_'+str(i+1),q) for i,q in enumerate(segments)]))
    allgroups={}
    for r in rows:allgroups.setdefault((r['source'],r['frequency_boundary']),[]).append(point(r))
    datasets.append(('neutral_variable_beta_all_candidates',VARIABLES,[zone(source+'_'+kind,q,style='scatter') for (source,kind),q in sorted(allgroups.items())]))
    reports=[]
    for name,variables,zones in datasets:reports.append(write_dataset(a.output_dir/(name+'.dat'),variables,zones,sha(candidate)))
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        input_files=[dict(path=str(f),sha256=sha(f)) for f in [candidate,status_path,derived_path]],
        format_guide=GUIDE,datasets=reports,source_unique_candidate_count=len(rows),
        upper_native_count=len(upper),additional_lower_count=len(lower),additional_low_upper_count=len(lowupper),
        qualification='Format conversion only. Unknowns and finite-family qualifications remain unchanged.')
    (a.output_dir/'tecplot_export_manifest.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    # A portable macro also documents scatter-only clouds and the non-neutral witness.
    for report in reports:
        lines=['#!MC 1410',"$!READDATASET '\"|MACROFILEPATH|/"+report['binary_file']+"\"'",'  READDATAOPTION = NEW','  RESETSTYLE = YES','  INITIALPLOTTYPE = XYLINE',
            '$!PLOTTYPE = XYLINE','$!DELETELINEMAPS']
        for i,z in enumerate(report['zones'],1):
            color=['BLUE','ORANGE','PURPLE'][i-1] if report['file'].startswith('neutral_variable_beta_stage') and i<=3 else 'BLACK'
            # Tecplot's standard colour set has no ORANGE.
            if color=='ORANGE':color='RED'
            lines+=['$!CREATELINEMAP',f'$!LINEMAP [{i}]',f'  NAME = {quote(z["name"])}',f'  ASSIGN{{ZONE={i} XAXISVAR=1 YAXISVAR=2}}',
                '  LINES{SHOW='+('YES' if z['style']=='line' else 'NO')+' COLOR='+color+'}',
                '  SYMBOLS{SHOW='+('NO' if z['style']=='line' else 'YES')+' COLOR='+('RED' if z['role'].startswith('not_neutral') else 'BLACK')+'}']
        lines+=['$!ACTIVELINEMAPS = [1-'+str(len(report['zones']))+']','$!VIEW FIT']
        if 'nose' in report['file']:lines+=['$!XYLINEAXIS XDETAIL 1 { COORDSCALE = LOG }']
        (a.output_dir/(Path(report['file']).stem+'.mcr')).write_text('\n'.join(lines)+'\n',encoding='ascii')
    print(json.dumps(dict(files=len(reports),upper=len(upper),additional=[len(lower),len(lowupper)],
        all_candidates=len(rows),ascii_readback='exact')))

if __name__=='__main__':main()
