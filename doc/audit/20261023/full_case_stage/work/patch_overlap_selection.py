"""Add explicit multi-eigenpair MAC selection; no default scalar selection changes."""
import difflib,pathlib
from full_driver import ROOT
S=ROOT/'work/solver_source';changes=[]

def save(path,old,new):
    assert old!=new
    (S/path).write_text(new)
    changes.append(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before/'+path,tofile='after/'+path)))

path='src/public/module_slepc_adapter.f90';old=(S/path).read_text();s=old
assert 'SelectEigenVectorByOverlap' not in s
s=s.replace('      logical :: track_target = .false. !< Select the eigenpair nearest the complex target for continuation','      logical :: track_target = .false. !< Select the eigenpair nearest the complex target for continuation\n      integer :: overlap_selected_index = 0 !< Explicit MAC choice for the next initial space only')
s=s.replace('      procedure :: CopySelectedEigenVector\n','      procedure :: CopySelectedEigenVector\n      procedure :: SelectEigenVectorByOverlap\n')
marker='   !> Reuse the last selected eigenvector; fall back safely before any query.'
helper='''   !> Explicitly select among this Solve's computed candidates by normalized
   !> Euclidean overlap. Residual/jump/overlap gates never relax. Scalar getters
   !> keep their original behavior unless this operation is explicitly called.
   subroutine SelectEigenVectorByOverlap(this, reference, reference_eigen, jump_tol, &
                                         error_tol, min_overlap, selected_eigen, owned_vector, error, overlap, valid)
      class(slepc_type), intent(inout) :: this
      Vec, intent(in) :: reference
      PetscScalar, intent(in) :: reference_eigen
      real(REAL64), intent(in) :: jump_tol, error_tol, min_overlap
      PetscScalar, intent(out) :: selected_eigen
      Vec, intent(out) :: owned_vector
      real(REAL64), intent(out) :: error, overlap
      logical, intent(out) :: valid
      PetscReal :: reference_norm, candidate_norm, candidate_error
      PetscScalar :: product
      real(REAL64) :: score, distance, best_score, second_score, best_distance, second_distance
      real(REAL64) :: best_error, second_error
      integer :: nconv, i, best_index, second_index, best_position, temporary, eligible

      selected_eigen = InvalidEigenvalue()
      owned_vector = PETSC_NULL_VEC
      error = huge(1.0_REAL64)
      overlap = 0.0_REAL64
      valid = .false.
      if (PetscObjectIsNull(reference)) return
      if (.not. ieee_is_finite(abs(reference_eigen))) return
      if (.not. ieee_is_finite(jump_tol) .or. jump_tol <= 0.0_REAL64) return
      if (.not. ieee_is_finite(error_tol) .or. error_tol <= 0.0_REAL64) return
      if (.not. ieee_is_finite(min_overlap) .or. min_overlap <= 0.0_REAL64 .or. min_overlap > 1.0_REAL64) return
      if (.not. allocated(this%xr) .or. .not. allocated(this%kr)) return
      nconv = size(this%kr)
      if (nconv == 0 .or. size(this%xr) /= nconv) return
      call VecNorm(reference, NORM_2, reference_norm, ierr)
      call CheckSolverConfiguration(ierr)
      if (.not. ieee_is_finite(reference_norm) .or. reference_norm <= 0.0_REAL64) return
      best_index = 0; second_index = 0; eligible = 0
      best_score = -1.0_REAL64; second_score = -1.0_REAL64
      best_distance = huge(1.0_REAL64); second_distance = huge(1.0_REAL64)
      best_error = huge(1.0_REAL64); second_error = huge(1.0_REAL64)
      do i = 1, nconv
         if (.not. ieee_is_finite(abs(this%kr(i)))) cycle
         distance = abs(this%kr(i) - reference_eigen)
         if (distance > jump_tol) cycle
         select type (this)
         type is (slepcpep_type)
            call PEPComputeError(this%solver, i-1, PEP_ERROR_ABSOLUTE, candidate_error, ierr)
         type is (slepceps_type)
            call EPSComputeError(this%solver, i-1, EPS_ERROR_ABSOLUTE, candidate_error, ierr)
         class default
            return
         end select
         call CheckSolverConfiguration(ierr)
         if (.not. ieee_is_finite(candidate_error)) cycle
         if (candidate_error < 0.0_REAL64 .or. candidate_error > error_tol) cycle
         call VecNorm(this%xr(i), NORM_2, candidate_norm, ierr)
         call CheckSolverConfiguration(ierr)
         if (.not. ieee_is_finite(candidate_norm) .or. candidate_norm <= 0.0_REAL64) cycle
         call VecDot(reference, this%xr(i), product, ierr)
         call CheckSolverConfiguration(ierr)
         score = (abs(product) / reference_norm / candidate_norm)**2
         if (.not. ieee_is_finite(score)) cycle
         score = min(1.0_REAL64, max(0.0_REAL64, score))
         eligible = eligible + 1
         if (score > best_score .or. (score == best_score .and. distance < best_distance)) then
            second_index = best_index; second_score = best_score
            second_distance = best_distance; second_error = best_error
            best_index = i; best_score = score; best_distance = distance; best_error = candidate_error
         else if (score > second_score .or. (score == second_score .and. distance < second_distance)) then
            second_index = i; second_score = score; second_distance = distance; second_error = candidate_error
         end if
      end do
      print *, 'MAC selection nconv=', nconv, ' eligible=', eligible, ' required=', min_overlap
      if (best_index == 0) then
         print *, 'MAC selection: no finite candidate passes residual, jump and norm gates.'
         return
      end if
      overlap = best_score
      print *, 'MAC best index=', best_index-1, ' eigen=', this%kr(best_index), &
               ' overlap=', best_score, ' residual=', best_error, ' jump=', best_distance
      if (second_index > 0) then
         print *, 'MAC second index=', second_index-1, ' eigen=', this%kr(second_index), &
                  ' overlap=', second_score, ' residual=', second_error, ' jump=', second_distance
      else
         print *, 'MAC second candidate: none within residual/jump/norm gates.'
      end if
      if (best_score < min_overlap) then
         print *, 'MAC accepted=', .false.
         return
      end if
      if (allocated(this%evIndex)) then
         if (size(this%evIndex) /= nconv) deallocate(this%evIndex)
      end if
      if (.not. allocated(this%evIndex)) then
         allocate(this%evIndex(nconv))
         this%evIndex = [(i, i=1,nconv)]
      end if
      best_position = 0
      do i = 1, nconv
         if (this%evIndex(i) == best_index) best_position = i
      end do
      if (best_position == 0) then
         this%evIndex = [(i, i=1,nconv)]
         best_position = best_index
      end if
      temporary = this%evIndex(nconv)
      this%evIndex(nconv) = best_index
      this%evIndex(best_position) = temporary
      this%last_error = best_error
      this%overlap_selected_index = best_index
      call this%CopySelectedEigenVector(selected_eigen, owned_vector, error, valid)
      print *, 'MAC accepted=', valid
   end subroutine SelectEigenVectorByOverlap

'''
assert marker in s;s=s.replace(marker,helper+marker)
oldguess='''      if (size(this%xr) == 0) return
      if (.not. allocated(this%evIndex)) return'''
newguess='''      if (size(this%xr) == 0) return
      candidate = this%overlap_selected_index
      if (candidate >= 1 .and. candidate <= size(this%xr)) then
         index = candidate
         return
      end if
      if (.not. allocated(this%evIndex)) return'''
assert oldguess in s;s=s.replace(oldguess,newguess)
marker='''      if (allocated(this%evIndex)) deallocate (this%evIndex)
      this%last_error = huge(1.0_REAL64)
   end subroutine ReleaseEigenpairs'''
assert marker in s;s=s.replace(marker,'''      if (allocated(this%evIndex)) deallocate (this%evIndex)
      this%last_error = huge(1.0_REAL64)
      this%overlap_selected_index = 0
   end subroutine ReleaseEigenpairs''')
save(path,old,s)

path='src/lst/lst_eqn.f90';old=(S/path).read_text();s=old
s=s.replace('   public :: type_lst_eqn, lst_eigen_error, lst_track_target, lst_is_tracking, lst_copy_selected_eigenvector','   public :: type_lst_eqn, lst_eigen_error, lst_track_target, lst_is_tracking, lst_copy_selected_eigenvector\n   public :: lst_select_eigenvector_by_overlap')
marker='   subroutine lst_setup(this)'
helper='''   !> Explicit multi-candidate overlap selection. owned_vector belongs to the
   !> caller only when valid=.true.; the caller must then VecDestroy it.
   subroutine lst_select_eigenvector_by_overlap(reference, reference_eigen, jump_tol, error_tol, min_overlap, &
                                              selected_eigen, owned_vector, error, overlap, valid)
      Vec, intent(in) :: reference
      complex(REAL64), intent(in) :: reference_eigen
      real(REAL64), intent(in) :: jump_tol, error_tol, min_overlap
      complex(REAL64), intent(out) :: selected_eigen
      Vec, intent(out) :: owned_vector
      real(REAL64), intent(out) :: error, overlap
      logical, intent(out) :: valid
      selected_eigen = cmplx(0.0_REAL64, 0.0_REAL64, REAL64)
      owned_vector = PETSC_NULL_VEC
      error = huge(1.0_REAL64); overlap = 0.0_REAL64; valid = .false.
      if (.not. allocated(Solver)) return
      call Solver%SelectEigenVectorByOverlap(reference, reference_eigen, jump_tol, error_tol, min_overlap, &
                                            selected_eigen, owned_vector, error, overlap, valid)
   end subroutine lst_select_eigenvector_by_overlap

'''
assert marker in s;s=s.replace(marker,helper+marker);save(path,old,s)
out=ROOT/'outputs/overlap_selection_api.patch';assert not out.exists();out.write_text(''.join(changes))
print(''.join(changes),flush=True)
print('OVERLAP SELECTION API WRITTEN; NOT BUILT pending root and profile changes.',flush=True)