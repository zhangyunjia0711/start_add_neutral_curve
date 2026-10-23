"""Build and run focused raw-vector/normalization regression after the shared source is ready."""
import datetime, hashlib, json, os, pathlib, shlex, subprocess, sys
ROOT=pathlib.Path(__file__).resolve().parent.parent;S=ROOT/'work/solver_source';B=S/'build'
ENV=dict(os.environ,PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PETSC_ARCH='arch-linux-c-debug',SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
ENV['PATH']='/usr/local/mpich-4.1.2/bin:'+ENV.get('PATH','');ENV['PKG_CONFIG_PATH']='/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/lib/pkgconfig:/usr/local/lib/pkgconfig'

def cmd(args):
    print('COMMAND '+json.dumps([str(v) for v in args]),flush=True)
    subprocess.run([str(v) for v in args],env=ENV,check=True)

def main():
    # Configure/build is intentionally separate from the initial patch script:
    # do not reapply/overwrite the original v1 patch or any archived manifests.
    cmd(['cmake','--build',B,'--target','START_LST','test_slepc_adapter','test_neutral_numerics','--parallel','4'])
    cmd(['ctest','--test-dir',B,'--output-on-failure','-R','neutral_numerics|slepc_adapter_safety'])
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','PETSc','slepc','fortran_stdlib'],env=ENV,check=True,capture_output=True,text=True).stdout)
    libs=shlex.split(subprocess.run(['pkg-config','--libs','PETSc','slepc','fortran_stdlib'],env=ENV,check=True,capture_output=True,text=True).stdout)
    objects=sorted(v for v in (B/'src/lst/CMakeFiles/START_LST.dir').rglob('*.o') if v.name!='main.f90.o')
    assert objects
    exe=ROOT/'work/test_raw_vector_copy'
    args=['/usr/local/mpich-4.1.2/bin/mpif90','-cpp','-fopenmp','-std=legacy','-I'+str(B/'src/lst')]+flags+[str(ROOT/'work/test_raw_vector_copy.f90')]+[str(v) for v in objects]+libs+['-llapack','-lblas','-Wl,-rpath,/home/zaza/ProgramFiles/petsc-3.22.1/arch-linux-c-debug/lib','-o',str(exe)]
    cmd(args)
    args=[str(exe),'-st_pc_factor_mat_solver_type','superlu','-malloc_debug','-options_left']
    print('COMMAND '+json.dumps(args),flush=True)
    p=subprocess.run(args,env=ENV,capture_output=True,text=True,timeout=120)
    log=p.stdout+p.stderr
    tag=sys.argv[1] if len(sys.argv)>1 else 'attempt2'
    outfile=ROOT/('outputs/raw_vector_regression_'+tag+'.log');assert not outfile.exists();outfile.write_text(log)
    for line in log.splitlines():
        if 'MEMORY BEFORE/AFTER/DELTA' in line or 'MEMORY PASS' in line or 'adapter_safety_test PASS' in line or 'PETSC ERROR:' in line:print(line,flush=True)
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),command=args,exit=p.returncode,petsc_errors=log.count('PETSC ERROR'),unused_options=log.count('Option left:'),passed='RAW COPY AND NORMALIZATION MEMORY PASS' in log,solver_sha256=hashlib.sha256((S/'bin/START_LST').read_bytes()).hexdigest(),source_hashes={str(v.relative_to(S)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [S/'src/lst/module_neutral_curve.f90',S/'src/lst/lst_eqn.f90',S/'src/lst/lst_solver_abstract.f90',S/'src/public/module_slepc_adapter.f90']})
    outfile=ROOT/('outputs/raw_vector_regression_'+tag+'.json');assert not outfile.exists();outfile.write_text(json.dumps(result,indent=2)+'\n')
    print('RAW REGRESSION '+json.dumps(result),flush=True)
    assert p.returncode==0 and result['petsc_errors']==0 and result['unused_options']==0 and result['passed']

if __name__=='__main__':main()