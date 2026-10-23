"""Fix eigenfunction-export temporaries and expose an owned raw selected vector copy."""
import difflib, pathlib
ROOT=pathlib.Path(__file__).resolve().parent.parent
S=ROOT/'work/solver_source'
changes=[]

def save(rel,old,new):
    assert old!=new,rel
    (S/rel).write_text(new)
    changes.append(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before/'+rel,tofile='after/'+rel)))

p='src/lst/lst_solver_abstract.f90';old=(S/p).read_text();s=old
start=s.index('   subroutine normal_u(v)');end=s.index('   end subroutine normal_u',start)
block=s[start:end]
lines=[]
for line in block.splitlines():
    lines.append(line)
    if line.lstrip().startswith('call ') and line.rstrip().endswith(', ierr)'):
        lines.append('      call CheckNormalizationError(ierr)')
block='\n'.join(lines)+'\n'
block+='      call VecScatterDestroy(ctx, ierr)\n      call CheckNormalizationError(ierr)\n      call VecDestroy(v_all, ierr)\n      call CheckNormalizationError(ierr)\n\n'
s=s[:start]+block+s[end:]
helper='''
   !> PETSc failures must not be hidden by subsequent temporary-object calls.
   subroutine CheckNormalizationError(error)
      PetscErrorCode, intent(in) :: error
      if (error /= 0) then
         call PetscErrorF(error, __LINE__, "lst_solver_abstract.f90")
         call MPIU_Abort(PETSC_COMM_WORLD, error)
      end if
   end subroutine CheckNormalizationError

'''
s=s.replace('end module',helper+'end module');save(p,old,s)

p='src/public/module_slepc_adapter.f90';old=(S/p).read_text();s=old
s=s.replace('only: ieee_value, ieee_quiet_nan','only: ieee_value, ieee_quiet_nan, ieee_is_finite')
s=s.replace('      procedure(Finalize_interface), deferred :: Finalize\n','      procedure(Finalize_interface), deferred :: Finalize\n      procedure :: CopySelectedEigenVector\n')
marker='   !> Reuse the last selected eigenvector; fall back safely before any query.'
helper='''   !> Copy the eigenvector that GetEigenVal actually selected. Do not call a
   !> sorting getter here: the tracking policy may have changed since selection.
   !> The caller owns eigenVector and must destroy it when valid is true.
   subroutine CopySelectedEigenVector(this, eigenvalue, eigenVector, error, valid)
      class(slepc_type), intent(inout) :: this
      PetscScalar, intent(out) :: eigenvalue
      Vec, intent(out) :: eigenVector
      real(REAL64), intent(out) :: error
      logical, intent(out) :: valid
      integer :: selected, nconv
      valid = .false.
      eigenvalue = InvalidEigenvalue()
      eigenVector = PETSC_NULL_VEC
      error = huge(1.0_REAL64)
      if (.not. allocated(this%xr)) return
      if (.not. allocated(this%kr)) return
      if (.not. allocated(this%evIndex)) return
      nconv = size(this%kr)
      if (nconv == 0 .or. size(this%xr) /= nconv .or. size(this%evIndex) /= nconv) return
      selected = this%evIndex(nconv)
      if (selected < 1 .or. selected > nconv) return
      if (.not. ieee_is_finite(real(this%kr(selected), REAL64))) return
      if (.not. ieee_is_finite(aimag(this%kr(selected)))) return
      if (.not. ieee_is_finite(this%last_error)) return
      if (this%last_error < 0.0_REAL64 .or. this%last_error >= huge(1.0_REAL64)) return
      ! A Fortran null-input sentinel is not a writable output argument.
      eigenVector%v = 0
      call VecDuplicate(this%xr(selected), eigenVector, ierr)
      call CheckSolverConfiguration(ierr)
      call VecCopy(this%xr(selected), eigenVector, ierr)
      call CheckSolverConfiguration(ierr)
      eigenvalue = this%kr(selected)
      error = this%last_error
      valid = .true.
   end subroutine CopySelectedEigenVector

'''
assert marker in s;s=s.replace(marker,helper+marker);save(p,old,s)

p='src/lst/lst_eqn.f90';old=(S/p).read_text();s=old
s=s.replace('public :: type_lst_eqn, lst_eigen_error, lst_track_target, lst_is_tracking','public :: type_lst_eqn, lst_eigen_error, lst_track_target, lst_is_tracking, lst_copy_selected_eigenvector')
marker='   subroutine lst_setup(this)'
helper='''   !> Owned raw eigenvector copy, eigenvalue and residual of the already selected
   !> pair. This avoids export normalization/scatter and never reselects a mode.
   !> On valid=.true. the caller must call VecDestroy(eigenvector).
   subroutine lst_copy_selected_eigenvector(eigenvalue, eigenvector, error, valid)
      complex(REAL64), intent(out) :: eigenvalue
      Vec, intent(out) :: eigenvector
      real(REAL64), intent(out) :: error
      logical, intent(out) :: valid
      valid = .false.
      eigenvalue = cmplx(0.0_REAL64, 0.0_REAL64, REAL64)
      eigenvector = PETSC_NULL_VEC
      error = huge(1.0_REAL64)
      if (.not. allocated(Solver)) return
      call Solver%CopySelectedEigenVector(eigenvalue, eigenvector, error, valid)
   end subroutine lst_copy_selected_eigenvector

'''
assert marker in s;s=s.replace(marker,helper+marker);save(p,old,s)

out=ROOT/'outputs/raw_vector_resource_fix.patch'
assert not out.exists(),str(out)
out.write_text(''.join(changes))
print(''.join(changes),flush=True)
print('RAW API: call lst_copy_selected_eigenvector(eigenvalue, eigenvector, error, valid); valid means caller owns VecDestroy. No sorting/get_eigf/normal_u involved.',flush=True)