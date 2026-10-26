"""Independent spatial-polynomial sign and derivative checks of the final bands."""
import json
import numpy as np
from common import ROOT,CASES,snapshot
from solver import Pencil,mac

def main():
    report=[]
    for case,c in CASES.items():
        data=json.loads((ROOT/'outputs'/f'{case}_final.json').read_text());native=[q for q in data['points'] if q['kind']==0]
        indices=sorted({int(q['xi']) for q in native});chosen=sorted(set([indices[2],c['seed_station'],indices[-1]]+[indices[j] for j in np.linspace(0,len(indices)-1,6,dtype=int)]))
        checks=[]
        for i in chosen:
            pair=sorted([q for q in native if int(q['xi'])==i],key=lambda q:q['beta_r']);width=pair[1]['beta_r']-pair[0]['beta_r'];p=Pencil(case,i);boundaries=[]
            for j,row in enumerate(pair):
                b=row['beta_r'];q=np.load(row['vector_file']);alpha=complex(row['alpha_r'],row['alpha_i']);delta=min(1e-4*(1+abs(b)),width*.01)
                minus,_,rm=p.solve_alpha(b-delta,alpha,q,False);plus,_,rp=p.solve_alpha(b+delta,alpha,q,False)
                fd=(plus-minus)/(2*delta);_,_,center=p.solve_alpha(b,alpha,q,True)
                derivative=complex(center['alpha_beta_r'],center['alpha_beta_i']);err=abs(fd-derivative)/max(1.,abs(derivative));assert err<2e-5,err
                inside=-plus.imag if j==0 else -minus.imag;outside=-minus.imag if j==0 else -plus.imag
                assert inside>0 and outside<0,(case,i,j,inside,outside)
                boundaries.append(dict(beta=b,inside_growth=inside,outside_growth=outside,relative_derivative_difference=err))
            # Trace the physical eigenvector from the high-|beta| boundary into the band.
            row=pair[0];b=row['beta_r'];alpha=complex(row['alpha_r'],row['alpha_i']);q=np.load(row['vector_file']);overlaps=[];mid=(pair[0]['beta_r']+pair[1]['beta_r'])/2
            for beta in np.linspace(b,mid,41)[1:]:
                alpha,q,r=p.solve_alpha(float(beta),alpha,q,False);overlaps.append(r['raw_MAC'])
            assert -alpha.imag>0,(case,i,'mid-band growth',-alpha.imag)
            checks.append(dict(station=i,x_mm=c['x_mm'][i-1],boundaries=boundaries,midpoint_beta=mid,midpoint_growth=-alpha.imag,minimum_continuation_MAC=min(overlaps)))
            print('AUDIT',case,i,round(-alpha.imag,8),flush=True)
        report.append(dict(case=case,cross_sections=checks))
    (ROOT/'outputs/band_audit.json').write_text(json.dumps(dict(cases=report,inputs=snapshot()),indent=2)+'\n')
    print('AUDIT_PASS',flush=True)

if __name__=='__main__':main()
