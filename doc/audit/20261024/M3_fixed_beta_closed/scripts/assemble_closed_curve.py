"""Publish the closed M3 curve from native roots and re-solved physical folds.

Operator-interpolation arc nodes are used only to determine ordering and are
not exported as physical curve points. The first point is repeated only after
the independent section-return/eigenvector closure test has passed.
"""
import csv,datetime,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,sha,snapshot
from profile_interpolation import mac

FIELDS=['order','x_mm','omega_r','beta_r','alpha_r','alpha_i','native_station','xi',
    'temporal_imag_residual','polynomial_residual','point_kind','source_file']

def writecsv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    O=ROOT/'outputs';arc_path=O/'M3b7_closed_arc_v2.json';arc=json.loads(arc_path.read_text())
    check_path=O/'closed_contour_validation.json';check=json.loads(check_path.read_text());assert check['closure']['passed']
    phys_path=O/'physical_turning_points.json';phys=json.loads(phys_path.read_text());assert phys['all_passed']
    native_manifest=O/'native_neutral_points.json';native=json.loads(native_manifest.read_text());assert native['finished']
    spatial=json.loads((O/'native_spatial_verification.json').read_text());assert spatial['all_passed']
    assert all(v['status']=='converged' for v in native['results'])
    native_data={};unique={};duplicates=[]
    for entry in native['results']:
        doc=json.loads(Path(entry['file']).read_text());point=doc['point'];native_data[(entry['station'],entry['ordinal'])]=doc
        hit=next((old for old in unique.get(entry['station'],[]) if abs(old['point']['alpha_r']-point['alpha_r'])<1e-6
            and abs(old['point']['omega_r']-point['omega_r'])<1e-6),None)
        if hit:
            overlap=mac(np.load(hit['vector_file']),np.load(doc['vector_file']));assert overlap>1-1e-7
            duplicates.append(dict(file=entry['file'],same_root_as=hit['source_file'],raw_MAC=overlap))
        else:
            doc['source_file']=entry['file'];unique.setdefault(entry['station'],[]).append(doc)
    counts={i:len(v) for i,v in unique.items()};assert set(counts)==set(range(10,712))
    for i,n in counts.items():assert n==(4 if 439<=i<=478 or 578<=i<=617 else 2),(i,n)
    bysegment={}
    for job in native['jobs']:bysegment.setdefault(job['arc_segment'],[]).append(job)
    oldfold={v['arc_segment'][0]:(j,v) for j,v in enumerate(check['folds'])};aug=[]
    for r in arc['records']:
        aug.append(r)
        if r['index'] in oldfold:
            j,f=oldfold[r['index']];aug.append(dict(index='fold',point=f['point'],physical_fold=j))
    sequence=[]
    def append(point,kind,vector_file,source_file,station=0):
        row=dict(x_mm=point['x_mm'],omega_r=point['omega_r'],beta_r=-7.,alpha_r=point['alpha_r'],alpha_i=0.,
            native_station=station,xi=point.get('xi',float(station)),temporal_imag_residual=abs(point['omega_i']),
            polynomial_residual=point['original_real_frequency_polynomial_residual'],point_kind=kind,
            source_file=str(source_file),vector_file=str(vector_file))
        assert row['temporal_imag_residual']<1e-9 and row['polynomial_residual']<1e-7
        if sequence and all(abs(sequence[-1][k]-row[k])<1e-6 for k in ['x_mm','omega_r','alpha_r']):return
        sequence.append(row)
    for j,(a,b) in enumerate(zip(aug,aug[1:])):
        if 'physical_fold' in a:
            f=phys['results'][a['physical_fold']];append(f['point'],1,f['vector_file'],O/('physical_fold_%d.json'%f['index']))
        jobs=bysegment.get(j,[]);ax,bx=a['point']['xi'],b['point']['xi']
        jobs=sorted(jobs,key=lambda z:(z['station']-ax)/(bx-ax))
        for job in jobs:
            doc=native_data[(job['station'],job['ordinal'])]
            append(doc['point'],0,doc['vector_file'],O/'native_points'/('i%03d_p%d.json'%(job['station'],job['ordinal'])),job['station'])
    assert len(sequence)>1000
    # Remove the already verified terminal copy and append the canonical native
    # start exactly once; the untouched independently corrected return is above.
    assert all(abs(sequence[-1][k]-sequence[0][k])<1e-6 for k in ['x_mm','omega_r','alpha_r'])
    sequence.pop();sequence.append(dict(sequence[0],point_kind=2))
    for i,r in enumerate(sequence):r['order']=i
    overlaps=[mac(np.load(a['vector_file']),np.load(b['vector_file'])) for a,b in zip(sequence,sequence[1:])]
    native_rows=[]
    for station,docs in sorted(unique.items()):
        for number,doc in enumerate(sorted(docs,key=lambda z:z['point']['omega_r'])):
            q=doc['point'];native_rows.append(dict(order=number,x_mm=q['x_mm'],omega_r=q['omega_r'],beta_r=-7.,
                alpha_r=q['alpha_r'],alpha_i=0.,native_station=station,xi=station,temporal_imag_residual=abs(q['omega_i']),
                polynomial_residual=q['original_real_frequency_polynomial_residual'],point_kind=0,source_file=doc['source_file']))
    assert len(native_rows)==1564
    writecsv(O/'M3_closed_beta_minus7_curve.csv',sequence,FIELDS)
    writecsv(O/'M3_closed_beta_minus7_native_points.csv',native_rows,FIELDS)
    plotrows=[r for r in sequence if r['point_kind']!=2]
    folds=[r for r in sequence if r['point_kind']==1]
    X=np.array([r['x_mm'] for r in sequence]);W=np.array([r['omega_r'] for r in sequence])
    fig,ax=plt.subplots(figsize=(11,6),constrained_layout=True)
    ax.plot(X,W,color='#14578c',lw=1.6,label=r'M3 neutral curve, $\beta=-7$')
    ax.scatter([r['x_mm'] for r in folds],[r['omega_r'] for r in folds],s=28,marker='D',
        facecolors='white',edgecolors='#9b3d89',linewidths=1,zorder=4,label='Solved turning points (interpolated local profiles)')
    ax.axhline(0,color='.65',lw=.7);ax.axvline(0,color='.65',lw=.7);ax.grid(alpha=.18)
    ax.set(xlim=(-45,625),ylim=(-9,10),xlabel='Global wall x (mm)',ylabel=r'Signed $\omega$ (dimensionless)',
        title=r'Closed M3 neutral curve: fixed real $\beta=-7$')
    ax.legend(loc='upper right',frameon=False,fontsize=9)
    fig.savefig(O/'M3_closed_beta_minus7.png',dpi=220);fig.savefig(O/'M3_closed_beta_minus7.svg');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True)
    nose=[r for r in sequence if r['x_mm']<-30]
    ax.plot([r['x_mm']+35 for r in nose],[r['omega_r'] for r in nose],'.-',ms=2.5,lw=1.2,color='#14578c')
    ax.set(xscale='log',xlabel='Global wall x + 35 mm (mm; logarithmic scale)',ylabel=r'Signed $\omega$',title='M3 closed curve: nose detail')
    ax.grid(alpha=.2);fig.savefig(O/'M3_closed_beta_minus7_nose.png',dpi=200);plt.close(fig)
    variables=['X_mm','Omega','Beta','Alpha_r','Alpha_i','NativeStation','Kind','GammaResidual','PolynomialResidual']
    keys=['x_mm','omega_r','beta_r','alpha_r','alpha_i','native_station','point_kind','temporal_imag_residual','polynomial_residual']
    with (O/'M3_closed_beta_minus7.dat').open('w',encoding='ascii') as f:
        f.write('TITLE="M3 closed neutral curve, beta=-7"\nVARIABLES='+ ' '.join('"'+v+'"' for v in variables)+'\n')
        f.write('DATASETAUXDATA Model="Frozen wall-local parallel model; zero explicit curvature"\n')
        f.write('DATASETAUXDATA PointKinds="0 native neutral root; 1 solved interpolated-flow turning point; 2 verified canonical start repeat"\n')
        f.write('DATASETAUXDATA Units="X_mm in mm; Omega and Beta dimensionless; Omega remains signed"\n')
        f.write('ZONE T="M3_beta_minus7_closed", I=%d, DATAPACKING=POINT\n'%len(sequence))
        f.write('DT=('+ ' '.join(['DOUBLE']*len(variables))+')\n')
        for row in sequence:f.write(' '.join('%.17e'%row[k] for k in keys)+'\n')
    summary=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),status='closed_M3_component_computed',
        fixed_beta=-7.,native_station_range=[min(counts),max(counts)],native_station_count=len(counts),native_unique_neutral_points=len(native_rows),
        closed_curve_rows=len(sequence),solved_interpolated_flow_turning_points=len(folds),
        curve_first_last_coordinates_exactly_equal=all(sequence[0][k]==sequence[-1][k] for k in ['x_mm','omega_r','beta_r','alpha_r','alpha_i']),
        terminal_point_kind=2,terminal_point_note='Coordinates repeat the independently verified start; Kind intentionally differs from the first row.',
        closure_validation=check['closure'],minimum_adjacent_display_point_raw_MAC=min(overlaps),
        maximum_neutral_temporal_imag=max(r['temporal_imag_residual'] for r in sequence),
        maximum_original_real_frequency_polynomial_residual=max(r['polynomial_residual'] for r in sequence),
        x_range_mm=[float(X.min()),float(X.max())],observed_frequency_range=[float(W.min()),float(W.max())],
        native_four_root_station_intervals=[[439,478],[578,617]],deduplicated_jobs=duplicates,
        native_spatial_independent_checks=len(spatial['results']),native_spatial_maximum_abs_growth=max(abs(r['actual_spatial_point']['growth']) for r in spatial['results']),
        native_spatial_minimum_raw_MAC=min(r['raw_MAC'] for r in spatial['results']),original_inputs=snapshot(),
        turning_points=[{k:r[k] for k in ['xi','x_mm','omega_r','alpha_r','temporal_imag_residual']} for r in folds],
        qualification='One connected M3 seed continuation, not enumeration of every eigenmode at beta=-7. Final points use native operators or re-solved interpolated local-flow profiles; operator-linear guide nodes are excluded.')
    (O/'M3_closed_beta_minus7_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k] for k in ['status','native_unique_neutral_points','closed_curve_rows','x_range_mm','minimum_adjacent_display_point_raw_MAC','maximum_neutral_temporal_imag','native_spatial_maximum_abs_growth']}))

if __name__=='__main__':main()
