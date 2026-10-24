"""Publish native rough-case neutral roots and re-solved physical-profile folds."""
import csv,datetime,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,sha,snapshot
from profile_interpolation import mac
STEM='rough_M3_closed_beta_minus7'
FIELDS=['order','x_mm','omega_r','beta_r','alpha_r','alpha_i','native_station','xi','temporal_imag_residual','polynomial_residual','point_kind','source_file']

def writecsv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    O=ROOT/'outputs';arc=json.loads((O/'M3b7_closed_arc_v3.json').read_text());assert arc['closure']['passed']
    check=json.loads((O/'closed_contour_validation.json').read_text());assert check['closure']['passed']
    phys=json.loads((O/'physical_turning_points.json').read_text());assert phys['all_passed']
    native=json.loads((O/'native_neutral_points.json').read_text());assert native['finished'] and all(e['status']=='converged' for e in native['results'])
    spatial=json.loads((O/'native_spatial_verification.json').read_text());assert spatial['all_passed']
    documents={};unique={};duplicates=[]
    for entry in native['results']:
        doc=json.loads(Path(entry['file']).read_text());p=doc['point'];documents[(entry['station'],entry['ordinal'])]=doc;doc['source_file']=entry['file']
        same=next((v for v in unique.get(entry['station'],[]) if abs(v['point']['alpha_r']-p['alpha_r'])<1e-6 and abs(v['point']['omega_r']-p['omega_r'])<1e-6),None)
        if same:
            overlap=mac(np.load(same['vector_file']),np.load(doc['vector_file']));assert overlap>1-1e-7
            duplicates.append(dict(source_file=entry['file'],same_as=same['source_file'],raw_MAC=overlap))
        else:unique.setdefault(entry['station'],[]).append(doc)
    counts={i:len(v) for i,v in unique.items()};assert set(counts)==set(range(min(counts),max(counts)+1))
    bysegment={}
    for job in native['jobs']:bysegment.setdefault(job['arc_segment'],[]).append(job)
    folded={f['arc_segment'][0]:(j,f) for j,f in enumerate(check['folds'])};aug=[]
    for r in arc['records']:
        aug.append(r)
        if r['index'] in folded:
            j,f=folded[r['index']];aug.append(dict(index='fold',point=f['point'],physical_fold=j))
    sequence=[]
    def append(point,kind,vector_file,source_file,station=0):
        row=dict(x_mm=point['x_mm'],omega_r=point['omega_r'],beta_r=-7.,alpha_r=point['alpha_r'],alpha_i=0.,native_station=station,
            xi=point.get('xi',float(station)),temporal_imag_residual=abs(point['omega_i']),
            polynomial_residual=point['original_real_frequency_polynomial_residual'],point_kind=kind,source_file=str(source_file),vector_file=str(vector_file))
        assert row['temporal_imag_residual']<1e-9 and row['polynomial_residual']<1e-7
        if sequence and all(abs(sequence[-1][k]-row[k])<1e-6 for k in ['x_mm','omega_r','alpha_r']):return
        sequence.append(row)
    for j,(a,b) in enumerate(zip(aug,aug[1:])):
        if 'physical_fold' in a:
            f=phys['results'][a['physical_fold']];append(f['point'],1,f['vector_file'],O/f'physical_fold_{f["index"]}.json')
        ax,bx=a['point']['xi'],b['point']['xi'];jobs=bysegment.get(j,[])
        jobs=sorted(jobs,key=lambda e:(e['station']-ax)/(bx-ax))
        for job in jobs:
            doc=documents[(job['station'],job['ordinal'])];append(doc['point'],0,doc['vector_file'],doc['source_file'],job['station'])
    if all(abs(sequence[-1][k]-sequence[0][k])<1e-6 for k in ['x_mm','omega_r','alpha_r']):sequence.pop()
    expected=sum(counts.values());assert sum(r['point_kind']==0 for r in sequence)==expected
    assert sum(r['point_kind']==1 for r in sequence)==len(phys['results'])
    sequence.append(dict(sequence[0],point_kind=2))
    for i,r in enumerate(sequence):r['order']=i
    overlaps=[mac(np.load(a['vector_file']),np.load(b['vector_file'])) for a,b in zip(sequence,sequence[1:])];assert min(overlaps)>.95
    original_rows=[]
    for station,docs in sorted(unique.items()):
        for j,doc in enumerate(sorted(docs,key=lambda d:d['point']['omega_r'])):
            p=doc['point'];original_rows.append(dict(order=j,x_mm=p['x_mm'],omega_r=p['omega_r'],beta_r=-7.,alpha_r=p['alpha_r'],alpha_i=0.,native_station=station,
                xi=station,temporal_imag_residual=abs(p['omega_i']),polynomial_residual=p['original_real_frequency_polynomial_residual'],point_kind=0,source_file=doc['source_file']))
    writecsv(O/(STEM+'_curve.csv'),sequence);writecsv(O/(STEM+'_native_points.csv'),original_rows)
    X=np.array([r['x_mm'] for r in sequence]);W=np.array([r['omega_r'] for r in sequence]);folds=[r for r in sequence if r['point_kind']==1]
    fig,ax=plt.subplots(figsize=(11,6),constrained_layout=True)
    ax.plot(X,W,color='#b05024',lw=1.6,label=r'Rough-case M3 neutral curve, $\beta=-7$')
    ax.scatter([r['x_mm'] for r in folds],[r['omega_r'] for r in folds],s=28,marker='D',facecolors='white',edgecolors='#733c86',linewidths=1,zorder=4,
        label='Solved turning points (interpolated local profiles)')
    ax.axhline(0,color='.65',lw=.7);ax.axvline(0,color='.65',lw=.7);ax.grid(alpha=.18)
    ax.set(xlim=(min(-45,X.min()-10),X.max()+20),ylim=(W.min()-.8,W.max()+1),xlabel='Global wall x (mm)',ylabel=r'Signed $\omega$ (dimensionless)',
        title=r'Rough case: closed M3 neutral curve, fixed real $\beta=-7$')
    ax.legend(loc='upper right',frameon=False,fontsize=9);fig.savefig(O/(STEM+'.png'),dpi=220);fig.savefig(O/(STEM+'.svg'));plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True);nose=[r for r in sequence if r['x_mm']<-30]
    ax.plot([r['x_mm']+35 for r in nose],[r['omega_r'] for r in nose],'.-',ms=2.5,lw=1.2,color='#b05024')
    ax.set(xscale='log',xlabel='Global wall x + 35 mm (mm; logarithmic scale)',ylabel=r'Signed $\omega$',title='Rough-case M3 curve: nose detail')
    ax.grid(alpha=.2);fig.savefig(O/(STEM+'_nose.png'),dpi=200);plt.close(fig)
    variables=['X_mm','Omega','Beta','Alpha_r','Alpha_i','NativeStation','Kind','GammaResidual','PolynomialResidual']
    keys=['x_mm','omega_r','beta_r','alpha_r','alpha_i','native_station','point_kind','temporal_imag_residual','polynomial_residual']
    with (O/(STEM+'.dat')).open('w',encoding='ascii') as f:
        f.write('TITLE="Rough-case M3 closed neutral curve, beta=-7"\nVARIABLES='+' '.join('"'+v+'"' for v in variables)+'\n')
        f.write('DATASETAUXDATA Model="Wall-local parallel, zero explicit curvature; supplied FPP baseflow"\n')
        f.write('DATASETAUXDATA PointKinds="0 native root; 1 solved interpolated-flow turning point; 2 verified start repeat"\n')
        f.write('DATASETAUXDATA Units="X_mm mm; Omega/Beta dimensionless; Omega remains signed"\n')
        f.write('ZONE T="rough_M3_beta_minus7_closed", I=%d, DATAPACKING=POINT\n'%len(sequence))
        f.write('DT=('+ ' '.join(['DOUBLE']*9)+')\n')
        for r in sequence:f.write(' '.join('%.17e'%r[k] for k in keys)+'\n')
    summary=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),status='closed_rough_M3_component_computed',fixed_beta=-7.,
        native_station_range=[min(counts),max(counts)],native_station_count=len(counts),native_unique_neutral_points=expected,native_root_counts=counts,
        closed_curve_rows=len(sequence),solved_interpolated_flow_turning_points=len(folds),deduplicated_jobs=duplicates,
        curve_first_last_coordinates_exactly_equal=all(sequence[0][k]==sequence[-1][k] for k in ['x_mm','omega_r','beta_r','alpha_r','alpha_i']),
        terminal_kind=2,closure_validation=check['closure'],minimum_adjacent_display_point_raw_MAC=min(overlaps),
        maximum_neutral_temporal_imag=max(r['temporal_imag_residual'] for r in sequence),maximum_original_real_frequency_polynomial_residual=max(r['polynomial_residual'] for r in sequence),
        x_range_mm=[float(X.min()),float(X.max())],observed_frequency_range=[float(W.min()),float(W.max())],
        native_spatial_independent_checks=len(spatial['results']),native_spatial_maximum_abs_growth=max(abs(r['actual_spatial_point']['growth']) for r in spatial['results']),
        native_spatial_minimum_raw_MAC=min(r['raw_MAC'] for r in spatial['results']),native_spatial_side_sign_crossings=sum(r['side_growth_sign_crossing'] for r in spatial['results']),
        original_inputs=snapshot(),turning_points=[{k:r[k] for k in ['xi','x_mm','omega_r','alpha_r','temporal_imag_residual']} for r in folds],
        qualification='One M3 component supported by physical-normal eigenfunction alignment; no all-spectrum or explicit-curvature/CFD-resolution completeness claim. Operator-guide nodes are excluded from final public points.')
    (O/(STEM+'_summary.json')).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k] for k in ['status','native_unique_neutral_points','closed_curve_rows','native_station_range','x_range_mm','minimum_adjacent_display_point_raw_MAC','native_spatial_maximum_abs_growth']}))

if __name__=='__main__':main()
