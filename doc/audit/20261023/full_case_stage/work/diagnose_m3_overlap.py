"""Independent frozen-v4 multi-pair MAC diagnostic; no neutral scan or source edits."""
import configparser,csv,datetime,hashlib,json,math,pathlib,re,struct,subprocess,time
import numpy as np
from full_driver import ROOT,CASE,ENV
BIN=ROOT/'outputs/bin/START_LST_vectors_v4'
EXPECTED='64344bee5324a25deaa2c1d669092121e1be690639eaedb75405189df19a1d6c'
TARGET=complex(51.76909408658165,-.5594665817891874)
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'

def vectors(path):
    data=path.read_bytes();offset=0;arrays=[]
    while offset<len(data):
        assert len(data)-offset>=8,(path,offset,len(data))
        cls,n=struct.unpack_from('>ii',data,offset);offset+=8
        assert cls==1211214 and n==2000,(cls,n,offset)
        end=offset+16*n;assert end<=len(data)
        arrays.append(np.frombuffer(data[offset:end],dtype='>c16').astype(np.complex128));offset=end
    return arrays

def probe(station,omega,nev):
    name=f'm3mac_i{station}_w{omega:+.1e}_nev{nev}';d=ROOT/'runs'/name;d.mkdir(exist_ok=True)
    out=d/'diagnostic.json'
    if out.exists():
        meta=json.loads(out.read_text());return meta,vectors(d/'out/all_vectors.petsc')
    assert not (d/'solver.log').exists(),'incomplete existing diagnostic'
    cfg=configparser.ConfigParser();cfg.read(CASE/'M1/config.ini')
    cfg['filenames']['dir']=str(CASE/'M1/in')+'/'
    cfg['domain']['iloc']=str(station);cfg['lst']['profile_coordinates']='global'
    cfg['lst']['alpha']='%.16g, %.16g'%(TARGET.real+.001,TARGET.imag)
    cfg['lst']['omega']='%.16g, 0'%omega;cfg['lst']['beta']='-7, 0'
    cfg['neutral']={'enabled':'F'}
    with (d/'config.ini').open('w') as f:cfg.write(f)
    options=['-pep_tol','1e-12','-pep_max_it','1000','-st_pc_factor_mat_solver_type','superlu','-pep_nev',str(nev),'-pep_ncv','16','-pep_view_vectors','binary:out/all_vectors.petsc','-options_left']
    cmd=[str(BIN),'-f','config.ini']+options;t=time.time()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        p=subprocess.run(cmd,cwd=d,env=ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=120)
    text=(d/'solver.log').read_text(errors='replace');pairs=[]
    for line in text.splitlines():
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if m:pairs.append(dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    vv=vectors(d/'out/all_vectors.petsc')
    meta=dict(name=name,station=station,omega=omega,beta=-7,nev=nev,ncv=16,coordinates='global',source_mode='M1 main841',solver=str(BIN),solver_sha256=EXPECTED,seconds=time.time()-t,exit=p.returncode,petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),pairs=pairs,vector_file=str(d/'out/all_vectors.petsc'),vector_sha256=hashlib.sha256((d/'out/all_vectors.petsc').read_bytes()).hexdigest(),command=cmd)
    assert len(vv)==len(pairs) and len(vv)>=nev and p.returncode==0 and meta['petsc_errors']==0 and meta['unused_options']==0
    out.write_text(json.dumps(meta,indent=2)+'\n');print('MAC PROBE '+json.dumps({k:meta[k] for k in ['name','seconds','exit','petsc_errors','unused_options']}),flush=True)
    return meta,vv

def main():
    assert hashlib.sha256(BIN.read_bytes()).hexdigest()==EXPECTED
    ref,rv=probe(119,0.,8)
    ref_index=min(range(len(ref['pairs'])),key=lambda k:abs(complex(ref['pairs'][k]['alpha_r'],ref['pairs'][k]['alpha_i'])-TARGET))
    reference=rv[ref_index];ref_eigen=complex(ref['pairs'][ref_index]['alpha_r'],ref['pairs'][ref_index]['alpha_i'])
    assert abs(ref_eigen-TARGET)<1e-4 and ref['pairs'][ref_index]['residual']<=1e-7
    records=[];summaries=[]
    for station,omega,nev in [(120,0.,4),(120,0.,8),(120,-1e-8,8),(120,1e-8,8)]:
        meta,vv=probe(station,omega,nev)
        candidates=[]
        for pair,vec in zip(meta['pairs'],vv):
            alpha=complex(pair['alpha_r'],pair['alpha_i']);overlap=float(abs(np.vdot(reference,vec))**2/(np.vdot(reference,reference).real*np.vdot(vec,vec).real))
            record=dict(run=meta['name'],station=station,omega=omega,beta=-7,nev=nev,ncv=16,coordinates='global',**pair,growth=-pair['alpha_i'],delta_alpha=abs(alpha-ref_eigen),mac=overlap,residual_accepted=pair['residual']<=1e-7,jump_accepted=abs(alpha-ref_eigen)<=2.,mac_accepted=overlap>=.95)
            record['all_accepted']=record['residual_accepted'] and record['jump_accepted'] and record['mac_accepted'];records.append(record);candidates.append(record)
        nearest=min(candidates,key=lambda r:r['delta_alpha']);eligible=[r for r in candidates if r['residual_accepted'] and r['jump_accepted']]
        best=max(eligible,key=lambda r:r['mac']) if eligible else None
        summaries.append(dict(run=meta['name'],nearest=nearest,best_within_jump_and_residual=best,acceptable_candidates=sum(r['all_accepted'] for r in candidates),pairs=len(candidates)))
    output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),solver=str(BIN),sha256=EXPECTED,reference=dict(run=ref['name'],pair=ref['pairs'][ref_index],index=ref_index),thresholds=dict(jump=2.,residual=1e-7,mac=.95),summaries=summaries,candidates=records,scope='Independent single-station multi-pair observation; does not continue/close global or new wall-local curves')
    p=ROOT/'outputs/m3_overlap_diagnostic.json';assert not p.exists();p.write_text(json.dumps(output,indent=2)+'\n')
    with (ROOT/'outputs/m3_overlap_candidates.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    print('MAC DIAGNOSTIC '+json.dumps(summaries),flush=True)

if __name__=='__main__':main()