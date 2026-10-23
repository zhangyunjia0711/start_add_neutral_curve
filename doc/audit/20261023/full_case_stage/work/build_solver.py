"""Narrow isolated solver patch and build; no changes to original checkout."""
import datetime, hashlib, json, os, pathlib, subprocess
ROOT=pathlib.Path(__file__).resolve().parent.parent
SOURCE=ROOT/'work/solver_source'
ORIGINAL=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE')
OUT=ROOT/'outputs'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args,env=None):
    print('COMMAND '+json.dumps([str(v) for v in args]),flush=True)
    subprocess.run([str(v) for v in args],check=True,env=env)
def patch():
    paths=['bin/START_LST','src/public/module_slepc_adapter.f90','src/lst/module_lst.f90']
    baseline={p:sha(ORIGINAL/p) for p in paths}
    (OUT/'solver_original_integrity_before.json').write_text(json.dumps(baseline,indent=2)+'\n')
    p=SOURCE/'src/public/module_slepc_adapter.f90';s=p.read_text()
    if 'subroutine SetDefaultDirectSolver' not in s:
        for kind in ['PEP','EPS']:
            a=s.index('   subroutine '+kind+'SolverCreate');b=s.index('   end subroutine '+kind+'SolverCreate',a)
            block=s[a:b]
            block=block.replace('      ST :: st\n','')
            block=block.replace('call '+kind+'GetST(this%solver, st, ierr)','call '+kind+'GetST(this%solver, this%st, ierr)')
            block=block.replace('      call STSetType(st, STSINVERT, ierr)\n','      call SetDefaultDirectSolver(this)\n')
            lines=[]
            for line in block.splitlines():
                lines.append(line)
                if line.lstrip().startswith('call ') and line.rstrip().endswith(', ierr)'):
                    lines.append('      call CheckSolverConfiguration(ierr)')
            s=s[:a]+'\n'.join(lines)+'\n'+s[b:]
        start=s.index('      call PEPGetST(this%solver, this%st, ierr)',s.index('   subroutine PEPSolverSolve'))
        stop=s.index('      call PEPSetFromOptions(this%solver, ierr)',start)
        s=s[:start]+'      ! Direct-solver defaults were installed once in Create; preserve runtime overrides.\n'+s[stop:]
        start=s.index('      call EPSGetST(this%solver, this%st, ierr)',s.index('   subroutine EPSSolverSolve'))
        stop=s.index('      call EPSSetFromOptions(this%solver, ierr)',start)
        s=s[:start]+'      ! Direct-solver defaults were installed once in Create; preserve runtime overrides.\n'+s[stop:]
        s=s.replace('      call PEPSetFromOptions(this%solver, ierr)\n','      call PEPSetFromOptions(this%solver, ierr)\n      call CheckSolverConfiguration(ierr)\n')
        s=s.replace('      call EPSSetFromOptions(this%solver, ierr)\n','      call EPSSetFromOptions(this%solver, ierr)\n      call CheckSolverConfiguration(ierr)\n')
        marker='   !> '+chr(27714)+chr(35299)+chr(29305)+chr(24449)+chr(20540)+'\n   subroutine PEPSolverSolve'
        assert marker in s
        helper='''   !> Set defaults only while the new ST/PC has not been set up. Runtime options
   !> may override these defaults later, including a serial SuperLU factorization.
   subroutine SetDefaultDirectSolver(this)
      class(slepc_type), intent(inout) :: this
      call STSetType(this%st, STSINVERT, ierr)
      call CheckSolverConfiguration(ierr)
      call STGetKSP(this%st, this%ksp, ierr)
      call CheckSolverConfiguration(ierr)
      call KSPGetPC(this%ksp, this%pc, ierr)
      call CheckSolverConfiguration(ierr)
      call KSPSetType(this%ksp, KSPPREONLY, ierr)
      call CheckSolverConfiguration(ierr)
      call PCSetType(this%pc, PCLU, ierr)
      call CheckSolverConfiguration(ierr)
      call PCFactorSetMatSolverType(this%pc, MATSOLVERSUPERLU_DIST, ierr)
      call CheckSolverConfiguration(ierr)
   end subroutine SetDefaultDirectSolver

   !> Configuration failures are fatal: do not overwrite an unchecked error code
   !> and continue with a partially configured eigensolver.
   subroutine CheckSolverConfiguration(error)
      PetscErrorCode, intent(in) :: error
      if (error /= 0) then
         call PetscErrorF(error, __LINE__, "module_slepc_adapter.f90")
         call MPIU_Abort(PETSC_COMM_WORLD, error)
      end if
   end subroutine CheckSolverConfiguration

'''
        s=s.replace(marker,helper+marker)
        p.write_text(s)
    s=p.read_text().replace('call PetscErrorF(error)', 'call PetscErrorF(error, __LINE__, "module_slepc_adapter.f90")');p.write_text(s)
    p=SOURCE/'src/lst/module_lst.f90';s=p.read_text()
    if 'track_initial_target' not in s:
        s=s.replace('        integer :: k, ki, input_status, curve_status\n','        integer :: k, ki, input_status, curve_status\n        logical :: track_initial_target\n')
        marker='        call lst%solve(ome, alf, bat)'
        i=s.index(marker)
        insertion='''        ! Opt-in target selection for the first eigenpair of a neutral batch.
        ! Traditional interactive/default selection keeps its original behavior.
        track_initial_target = .false.
        if (neutral_batch_enabled()) then
            if (cfg%has_key('neutral','track_initial_target')) &
                call cfg%get('neutral','track_initial_target',track_initial_target)
            if (track_initial_target) call lst_track_target(.true.)
        end if
'''
        s=s[:i]+insertion+s[i:];p.write_text(s)
    diff=subprocess.run(['git','-C',str(SOURCE),'diff','--','src/public/module_slepc_adapter.f90','src/lst/module_lst.f90'],capture_output=True,text=True,check=True).stdout
    (OUT/'solver_source.patch').write_text(diff)
    print(diff,flush=True)
    return baseline

def main():
    assert subprocess.run(['git','-C',str(SOURCE),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()=='87995d66f6618c095f0a8302d8a7e8544c620a76'
    baseline=patch();env=dict(os.environ,PETSC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',PETSC_ARCH='arch-linux-c-debug',SLEPC_DIR='/home/zaza/ProgramFiles/petsc-3.22.1',LD_LIBRARY_PATH='/usr/local/mpich-4.1.2/lib:'+os.environ.get('LD_LIBRARY_PATH',''))
    env['PATH']='/usr/local/mpich-4.1.2/bin:'+env.get('PATH','')
    command(['cmake','-S',SOURCE,'-B',SOURCE/'build','-DBUILD_TESTING=ON'],env)
    command(['cmake','--build',SOURCE/'build','--target','START_LST','test_slepc_adapter','test_neutral_numerics','--parallel','4'],env)
    command(['ctest','--test-dir',SOURCE/'build','--output-on-failure','-R','neutral_numerics|slepc_adapter_safety'],env)
    after={p:sha(ORIGINAL/p) for p in baseline};assert after==baseline
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),original=after,new_binary=str(SOURCE/'bin/START_LST'),new_binary_sha256=sha(SOURCE/'bin/START_LST'),branch='codex/full-beta-20261022',base='87995d66f6618c095f0a8302d8a7e8544c620a76',diff='outputs/solver_source.patch',git_status=subprocess.run(['git','-C',str(SOURCE),'status','--short'],capture_output=True,text=True,check=True).stdout)
    (OUT/'solver_build_integrity.json').write_text(json.dumps(result,indent=2)+'\n')
    print('BUILD INTEGRITY '+json.dumps(result),flush=True)

if __name__=='__main__':main()