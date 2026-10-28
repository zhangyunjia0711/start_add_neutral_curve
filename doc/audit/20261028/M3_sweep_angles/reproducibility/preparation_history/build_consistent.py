"""Optional consistent wall-normal derivatives; preserve legacy coordinate modes."""
import datetime,difflib,hashlib,json,os,shutil,subprocess
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate')
ROOT=BASE/'neutral_m3_sweep_consistent_20261028'
SOURCE=BASE/'75°/bluntplate/neutral_beta_full_20261022/work/solver_source'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not ROOT.exists(),'Preserve previous build'
    for name in ['work','outputs','cases','runs','matrix_exports','native_pencils','inputs']:(ROOT/name).mkdir(parents=True)
    dest=ROOT/'work/solver_source';dest.mkdir()
    for name in ['src','tests']:shutil.copytree(SOURCE/name,dest/name)
    shutil.copy2(SOURCE/'CMakeLists.txt',dest/'CMakeLists.txt')
    p=dest/'src/public/wall_local_coordinates.f90';old=p.read_text();s=old
    s=s.replace('logical, save :: configured=.false., enabled=.false.','logical, save :: configured=.false., enabled=.false., consistent_normal_derivatives=.false.')
    s=s.replace("case('wall_local')","case('wall_local','wall_local_consistent')",1)
    s=s.replace('      configured=.true.','      consistent_normal_derivatives=trim(selection)==\'wall_local_consistent\'\n      configured=.true.',1)
    s=s.replace('maxval(abs(dx))<=16*epsilon(flat_scale)*flat_scale','maxval(abs(dx))<=16*epsilon(flat_scale)*flat_scale.and..not.consistent_normal_derivatives')
    s=s.replace("global or wall_local.'","global, wall_local or wall_local_consistent.'")
    assert s!=old;p.write_text(s)
    (ROOT/'outputs/consistent_normal_derivatives.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='a/src/public/wall_local_coordinates.f90',tofile='b/src/public/wall_local_coordinates.f90')))
    test=dest/'tests/wall_local_profiles_test.f90';s=test.read_text();needle="   open(newunit=unit,file='wall_local_profile_metadata.csv',status='old')"
    insert='''   ! Consistent mode must ignore sentinel global derivatives even on an exact flat column.
   call configure_wall_local_coordinates('wall_local_consistent',.true.,.true.,.false.,.true.)
   q=qq(:,3,:); qn=dq0dy(:,3,:); qnn=dq0dyy(:,3,:); en=eta_y(3,:); enn=eta_yy(3,:)
   call wall_local_grid(3,n,en,enn)
   call wall_local_flow(3,q,qn,qnn)
   if (maxval(abs(qn(2,:)-1.0_real64))>1e-10_real64.or.maxval(abs(qnn(2,:)))>1e-9_real64) &
      error stop 'Consistent flat derivative did not use the projected normal profile.'
   reference_en=en; reference_enn=enn; reference_qn=qn; reference_qnn=qnn
   do i=1,in
      do j=1,jn
         xx(i,j)=xx(i,j)+1e-13_real64*(yy(i,j)-yy(i,1))
      enddo
   enddo
   call configure_wall_local_coordinates('wall_local_consistent',.true.,.true.,.false.,.true.)
   q=qq(:,3,:); qn=dq0dy(:,3,:); qnn=dq0dyy(:,3,:); en=eta_y(3,:); enn=eta_yy(3,:)
   call wall_local_grid(3,n,en,enn)
   call wall_local_flow(3,q,qn,qnn)
   if (maxval(abs(qn-reference_qn))>1e-11_real64.or.maxval(abs(qnn-reference_qnn))>1e-10_real64) &
      error stop 'Consistent derivatives changed under flat coordinate roundoff.'
'''
    assert needle in s;test.write_text(s.replace(needle,insert+needle))
    env=dict(os.environ,PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PETSC_ARCH='arch-linux-c-debug',SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PKG_CONFIG_PATH='/usr/local/lib/pkgconfig',LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    commands=[['cmake','-S',str(dest),'-B',str(dest/'build')],['cmake','--build',str(dest/'build'),'--target','START_LST','test_wall_local_profiles','test_neutral_numerics','-j','4'],['ctest','--test-dir',str(dest/'build'),'-R','wall_local|neutral_numerics','--output-on-failure']]
    receipts=[]
    for k,cmd in enumerate(commands):
        log=ROOT/'outputs'/f'build_step{k}.log'
        with log.open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
        receipts.append(dict(command=cmd,returncode=r.returncode,log=str(log),sha256=sha(log)));print('BUILD_STEP',k,r.returncode,flush=True)
        assert r.returncode==0,log
    binary=dest/'bin/START_LST';(ROOT/'outputs/build_receipt.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),source=str(SOURCE),solver=str(binary),solver_sha256=sha(binary),new_coordinate_mode='wall_local_consistent',legacy_modes_preserved=['global','wall_local'],commands=receipts,source_hashes=[dict(path=str(p.relative_to(dest)),sha256=sha(p)) for p in sorted((dest/'src').rglob('*')) if p.is_file()]),indent=2)+'\n');print('BUILT',binary,sha(binary),flush=True)
if __name__=='__main__':main()
