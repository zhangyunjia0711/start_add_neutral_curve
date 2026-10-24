"""Three actual single-point PEP solves through the original interactive omega_i control."""
import configparser,datetime,hashlib,json,pathlib,re,struct,subprocess,sys,time
import numpy as np
AREA=pathlib.Path(__file__).resolve().parent
RUN=AREA.parent.parent
COMP=RUN/'outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148'
sys.path.insert(0,str(COMP));import compensate as c
CUT=RUN/'outputs/native10_beta_cuts_v10_20261023_2233'
SOLVER=RUN/'outputs/bin/START_LST_wall_local_v10'
SHA='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(n,v):(AREA/n).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
assert sha(SOLVER)==SHA
before=c.a.driver.snapshot()
source=json.loads((CUT/'results/N10_w40_center_6e8ec28957.json').read_text())
point=source['point'];c.AREA=CUT;refinfo,refq=c.vec(source['run'],point)
d=AREA/'runs/interactive_omega_i_three_points';d.mkdir(parents=True,exist_ok=True)
assert not (d/'solver.log').exists(),'Never overwrite a previous physical solve'
cfg=configparser.ConfigParser();cfg.read(CUT/'runs/N10_w40_center_6e8ec28957/config.ini')
cfg['lst']['omega']='40, 0';cfg['lst']['beta']='%.17g, 0'%point['beta_r']
cfg['lst']['alpha']='%.17g, %.17g'%(point['alpha_r']+1e-6,point['alpha_i'])
cfg['neutral']['enabled']='F'
with (d/'config.ini').open('w') as f:cfg.write(f)
stdin='7 1 0.005\n7 1 0.005\n-1 0 0\n';(d/'stdin.txt').write_text(stdin)
options=json.loads((CUT/'runs/N10_w40_center_6e8ec28957/summary.json').read_text())['spec']['options']
cmd=[str(SOLVER),'-f','config.ini']+options
started=stamp();t=time.time()
with (d/'solver.log').open('w') as f:
    f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
    try:p=subprocess.run(cmd,input=stdin,text=True,cwd=d,env=c.a.driver.ENV,stdout=f,stderr=subprocess.STDOUT,timeout=30);code=p.returncode
    except subprocess.TimeoutExpired:code=124
text=(d/'solver.log').read_text(errors='replace')
groups=[];flat=[]
for line in text.splitlines():
    if 'ncov=' in line:groups.append([])
    m=c.PAIR.search(line)
    if m:
        assert groups
        item=dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E')))
        groups[-1].append(item);flat.append(item)
data=(d/'out/eigenvectors.petsc').read_bytes();qs=[];off=0
while off<len(data):
    cls,n=struct.unpack_from('>ii',data,off);off+=8;end=off+16*n
    assert cls==1211214 and 0<n<1000000 and end<=len(data)
    qs.append(np.frombuffer(data[off:end],dtype='>c16').astype(np.complex128));off=end
assert len(flat)==len(qs)
localfile=next((d/'out').glob('*LOCAL_eigenvalues.plt'))
local=[]
for l in localfile.read_text().splitlines():
    try:
        v=[float(x) for x in l.split()]
        if len(v)==13:local.append(v)
    except ValueError:pass
assert code==0 and len(groups)==3 and len(local)==3,(code,len(groups),len(local))
offset=0;points=[];previous=refq
for k,(g,row) in enumerate(zip(groups,local)):
    om=complex(row[5],row[6]);alpha=complex(row[7],row[8])
    assert abs(om-complex(40,k*.005))<1e-10
    matches=[j for j,z in enumerate(g) if abs(c.point_eigen(z)-alpha)<1e-7]
    assert len(matches)==1,(k,matches)
    j=matches[0];z=g[j];q=qs[offset+j]
    mac=c.mac(previous,q)
    item=dict(station=int(row[0]),omega_r=om.real,omega_i=om.imag,beta_r=point['beta_r'],printed_beta_r=row[9],beta_i=row[10],
        alpha_r=z['alpha_r'],alpha_i=z['alpha_i'],growth=-z['alpha_i'],polynomial_residual=z['residual'],
        raw_vector_index=offset+j,solve_group=k,raw_vector_file=str(d/'out/eigenvectors.petsc'),
        raw_vector_file_sha256=sha(d/'out/eigenvectors.petsc'),raw_vector_component_sha256=hashlib.sha256(q.astype('>c16').tobytes()).hexdigest(),
        scalar_output_delta=abs(c.point_eigen(z)-alpha),MAC_to_previous_actual=mac,MAC_to_real_reference=c.mac(refq,q),
        vector_norm=float(np.linalg.norm(q)),passed=z['residual']<1e-7 and mac>=.95)
    points.append(item);previous=q;offset+=len(g)
sourcepaths=[RUN/'work/solver_source/src/lst/cfgio_adapter.f90',RUN/'work/solver_source/src/lst/module_lst.f90',
    RUN/'work/solver_source/src/lst/module_neutral_curve.f90',RUN/'work/solver_source/src/lst/lst_eqn.f90',RUN/'work/full_driver.py',
    pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE/src/lst/cfgio_adapter.f90'),
    pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE/src/lst/module_lst.f90'),
    pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE/src/lst/lst_eqn.f90')]
sourcechecks=[dict(path=str(p),sha256=sha(p)) for p in sourcepaths]
after=c.a.driver.snapshot();assert before==after and sha(SOLVER)==SHA
report=dict(created=stamp(),started=started,seconds=round(time.time()-t,3),command=cmd,stdin=stdin,exit=code,
    actual_pep_solves=len(groups),solve_group_sizes=[len(x) for x in groups],points=points,all_passed=all(x['passed'] for x in points),
    reference_point=point,reference_vector=refinfo,reference_run=source['name'],
    petsc_errors=len(re.findall('PETSC ERROR|SLEPC ERROR',text,re.I)),
    config_sha256=sha(d/'config.ini'),local_output_sha256=sha(localfile),solver_log_sha256=sha(d/'solver.log'),solver_sha256=SHA,
    original_inputs_and_original_binary_unchanged=True,source_checks=sourcechecks,
    interface_facts=[
        'Original lst cfgio reads two omega components; module_lst constructs omega(1)+i*omega(2).',
        'Original interactive k=7 adds i*pk and invokes the spatial PEP without dropping omega_i.',
        'Complex omega enters D0 += omega*G0 in spatial matrix assembly.',
        'The neutral entrance rejects non-real omega or beta above 1e-12, even when omega is not a scan axis.',
        'If omega is scanned, make_state explicitly sets its imaginary part to zero.',
        'full_driver.py formats the seed omega as real,0 and must not be used unchanged for complex-frequency work.',
        'This probe uses the original interactive loop and its legacy growth ordering; all selected eigenpairs were checked independently with actual raw-vector MAC and polynomial residual.',
        'No built-in omega_i MAC-controlled adaptive continuation is exposed by the neutral module; an independent driver can bound steps and accept/reject exported eigenpairs using raw MAC.'
    ],
    qualification='Three neighboring spatial eigenpairs only. This establishes complex-frequency interface and local vector continuity, not global analytic branch topology, pinching, physical absolute instability or downstream causality.')
write('complex_omega_probe.json',report)
write('complex_omega_integrity.json',dict(created=stamp(),before=before,after=after,unchanged=True,solver_sha256=SHA,source_checks=sourcechecks))
import csv
with (AREA/'complex_omega_points.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(points[0]));w.writeheader();w.writerows(points)
print(json.dumps(dict(seconds=report['seconds'],actual_pep_solves=len(groups),all_passed=report['all_passed'],points=points),indent=2))

