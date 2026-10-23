"""Add explicit one-shot PEP candidate retry; do not build shared source."""
import pathlib
ROOT=pathlib.Path(__file__).resolve().parent.parent
SRC=ROOT/'work/solver_source'
p=SRC/'src/public/module_slepc_adapter.f90'
s=p.read_text()
assert 'procedure :: PrepareCandidateRetry' not in s
s=s.replace('      procedure :: SelectEigenVectorByOverlap\n','      procedure :: SelectEigenVectorByOverlap\n      procedure :: PrepareCandidateRetry\n',1)
s=s.replace('      PEP :: solver !<多项式算子特征值问题求解器\n','''      PEP :: solver !<多项式算子特征值问题求解器
      integer :: retry_count = 0
      Vec :: retry_reference !< Owned copy, valid only with has_retry_reference
      logical :: has_retry_reference = .false.
      logical :: restore_dimensions = .false.
      PetscInt :: normal_nev, normal_ncv, normal_mpd
''',1)
start=s.index('      PetscInt :: nconv, nrow, guess_size',s.index('   subroutine PEPSolverSolve'))
end=s.index('      ! block\n',start)
s=s[:start]+'''      PetscInt :: nconv, nrow, guess_size, retry_ncv
      integer :: i, guess_index
      logical :: has_initialguess
      PetscScalar :: ki

      this%last_error = huge(1.0_REAL64)
      has_initialguess = .false.
      call PEPSetProblemType(this%solver, PEP_GENERAL, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPSetOperators(this%solver, nmat, mat, ierr)
      call CheckSolverConfiguration(ierr)
      call MatGetSize(mat(1), this%m, nrow, ierr)
      call CheckSolverConfiguration(ierr)
      target_value = target_input
      this%target = target_input
      call PEPSetTarget(this%solver, target_value, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPSetType(this%solver, PEPLINEAR, ierr)
      call CheckSolverConfiguration(ierr)
      ! A previous expansion must not remain sticky when options are absent.
      ! Production still supplies explicit -pep_nev 1 -pep_ncv 8 on normal solves.
      if (this%restore_dimensions) then
         call PEPSetDimensions(this%solver, this%normal_nev, this%normal_ncv, this%normal_mpd, ierr)
         call CheckSolverConfiguration(ierr)
         this%restore_dimensions = .false.
      end if
      call PEPSetFromOptions(this%solver, ierr)
      call CheckSolverConfiguration(ierr)
      if (this%retry_count > 0) then
         if (.not. this%has_retry_reference) then
            print *, 'Candidate retry rejected: no owned reference vector.'
            call CheckSolverConfiguration(PETSC_ERR_ARG_WRONGSTATE)
         end if
         call VecGetSize(this%retry_reference, guess_size, ierr)
         call CheckSolverConfiguration(ierr)
         if (guess_size /= this%m) then
            print *, 'Candidate retry rejected: reference/matrix size mismatch.', guess_size, this%m
            call ClearCandidateRetry(this)
            call CheckSolverConfiguration(PETSC_ERR_ARG_SIZ)
         end if
         call PEPGetDimensions(this%solver, this%normal_nev, this%normal_ncv, this%normal_mpd, ierr)
         call CheckSolverConfiguration(ierr)
         retry_ncv = max(24, 2*this%retry_count+8)
         call PEPSetDimensions(this%solver, this%retry_count, retry_ncv, PETSC_DECIDE, ierr)
         call CheckSolverConfiguration(ierr)
         this%restore_dimensions = .true.
         call PEPSetInitialSpace(this%solver, 1, this%retry_reference, ierr)
         call CheckSolverConfiguration(ierr)
         print *, 'Candidate retry consuming owned reference: nev=', this%retry_count, ' ncv=', retry_ncv
      else if (allocated(this%xr)) then
         if (size(this%xr) > 0) then
            guess_index = SelectedGuessIndex(this)
            call VecGetSize(this%xr(guess_index), guess_size, ierr)
            call CheckSolverConfiguration(ierr)
            if (guess_size == this%m) then
               call VecDuplicate(this%xr(guess_index), this%initialguess(1), ierr)
               call CheckSolverConfiguration(ierr)
               call VecCopy(this%xr(guess_index), this%initialguess(1), ierr)
               call CheckSolverConfiguration(ierr)
               call PEPSetInitialSpace(this%solver, 1, this%initialguess(1), ierr)
               call CheckSolverConfiguration(ierr)
               has_initialguess = .true.
            end if
         end if
      end if
      call PEPSolve(this%Solver, ierr)
      nconv = 0
      if (ierr == 0) then
         call PEPGetConverged(this%Solver, nconv, ierr)
         if (ierr /= 0) nconv = 0
      end if
      ! Both successful and failed PEPSolve returns consume this one-shot copy.
      call ClearCandidateRetry(this)
      if (has_initialguess) then
         call VecDestroy(this%initialguess(1), ierr)
         call CheckSolverConfiguration(ierr)
      end if
      call ReleaseEigenpairs(this)
      allocate (this%xr(nconv), this%kr(nconv))
      do i = 0, nconv - 1
         call MatCreateVecs(mat(1), this%xr(i + 1), PETSC_NULL_VEC, ierr)
         call CheckSolverConfiguration(ierr)
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, this%xr(i + 1), PETSC_NULL_VEC, ierr)
         call CheckSolverConfiguration(ierr)
      end do

'''+s[end:]
insert=s.index('   !> Reuse the last selected eigenvector; fall back safely before any query.')
s=s[:insert]+'''   !> Prepare exactly one expanded spatial PEP solve. The caller retains its
   !> reference; this object stores a true owned copy, replaced safely on repeat.
   !> Unsupported EPS/temporal use is an explicit configuration error.
   subroutine PrepareCandidateRetry(this, count, reference)
      class(slepc_type), intent(inout) :: this
      integer, intent(in) :: count
      Vec, intent(in) :: reference
      Vec :: created_vector
      PetscReal :: reference_norm
      select type (this)
      type is (slepcpep_type)
         if (this%lst_mode /= LST_SPATIAL_MODE) then
            print *, 'Candidate retry supports spatial PEP only.'
            call CheckSolverConfiguration(PETSC_ERR_SUP)
         end if
         if (count <= 0 .or. count > (huge(count)-8)/2) then
            print *, 'Candidate retry rejected: invalid requested count.', count
            call CheckSolverConfiguration(PETSC_ERR_ARG_OUTOFRANGE)
         end if
         if (PetscObjectIsNull(reference)) then
            print *, 'Candidate retry rejected: null reference.'
            call CheckSolverConfiguration(PETSC_ERR_ARG_NULL)
         end if
         call VecNorm(reference, NORM_2, reference_norm, ierr)
         call CheckSolverConfiguration(ierr)
         if (.not. ieee_is_finite(reference_norm) .or. reference_norm <= 0.0_REAL64) then
            print *, 'Candidate retry rejected: nonfinite or zero reference norm.'
            call CheckSolverConfiguration(PETSC_ERR_ARG_OUTOFRANGE)
         end if
         call VecDuplicate(reference, created_vector, ierr)
         call CheckSolverConfiguration(ierr)
         call VecCopy(reference, created_vector, ierr)
         call CheckSolverConfiguration(ierr)
         call ClearCandidateRetry(this)
         this%retry_reference = created_vector
         this%has_retry_reference = .true.
         this%retry_count = count
         print *, 'Candidate retry prepared: owned reference, count=', count
      class default
         print *, 'Candidate retry rejected: EPS is unsupported; spatial PEP required.'
         call CheckSolverConfiguration(PETSC_ERR_SUP)
      end select
   end subroutine PrepareCandidateRetry

   subroutine ClearCandidateRetry(this)
      class(slepcpep_type), intent(inout) :: this
      if (this%has_retry_reference) then
         call VecDestroy(this%retry_reference, ierr)
         call CheckSolverConfiguration(ierr)
      end if
      this%has_retry_reference = .false.
      this%retry_count = 0
   end subroutine ClearCandidateRetry

'''+s[insert:]
s=s.replace('            call VecDestroy(this%xr(i), ierr)\n','            call VecDestroy(this%xr(i), ierr)\n            call CheckSolverConfiguration(ierr)\n',1)
s=s.replace('      call EPSDestroy(this%solver, ierr)\n','      call EPSDestroy(this%solver, ierr)\n      call CheckSolverConfiguration(ierr)\n',1)
s=s.replace('      call ReleaseEigenpairs(this)\n      call PEPDestroy(this%solver, ierr)\n','      call ClearCandidateRetry(this)\n      this%restore_dimensions = .false.\n      call ReleaseEigenpairs(this)\n      call PEPDestroy(this%solver, ierr)\n      call CheckSolverConfiguration(ierr)\n',1)
p.write_text(s)
p=SRC/'src/lst/lst_eqn.f90';s=p.read_text()
s=s.replace('   public :: lst_select_eigenvector_by_overlap\n','   public :: lst_select_eigenvector_by_overlap, lst_prepare_candidate_retry\n',1)
marker='   subroutine lst_setup(this)\n'
assert marker in s
s=s.replace(marker,'''   !> One-shot spatial PEP candidate expansion. The solver copies reference;
   !> the caller keeps ownership and may destroy or reuse its own vector.
   subroutine lst_prepare_candidate_retry(count, reference)
      integer, intent(in) :: count
      Vec, intent(in) :: reference
      if (.not. allocated(Solver)) error stop 'Candidate retry requires an allocated LST solver.'
      call Solver%PrepareCandidateRetry(count, reference)
   end subroutine lst_prepare_candidate_retry

'''+marker,1)
p.write_text(s)
print('PATCHED one-shot retry adapter/lst_eqn; no build performed.')