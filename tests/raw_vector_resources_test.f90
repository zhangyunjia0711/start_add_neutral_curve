#include <slepc/finclude/slepc.h>
program raw_vector_copy_test
 use iso_fortran_env, only: REAL64
 use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
 use slepceps
 use slepcpep
 use petsc
 use mod_slepc_adapter
 use mod_lst_solver_abstract, only: normal_u
 implicit none
 type(slepceps_type) :: eps_adapter
 type(slepcpep_type) :: pep_adapter
 Mat :: matrices(3)
 Vec :: vector, raw, raw_saved, normalization_vector
 logical :: raw_valid
 PetscScalar :: raw_value
 real(REAL64) :: raw_error, raw_norm, saved_norm, selected_overlap
 PetscLogDouble :: memory_before, memory_after
 integer :: memory_pass
 PetscErrorCode :: ierr
 PetscScalar :: value, entry, default_value, tracked_value
 PetscScalar, allocatable :: values(:)
 Vec, allocatable :: vectors(:)
 PetscInt :: n, i, j, pass
 call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
 call eps_adapter%Create(1, 1)
 call pep_adapter%Create(1, 1)
 call eps_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (raw_valid) stop 101
 call pep_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (raw_valid) stop 102
 ! Before Solve and after an empty result, scalar getters must be safe.
 call eps_adapter%GetEigenVal(value)
 if (.not. ieee_is_nan(real(value))) stop 1
 allocate(eps_adapter%kr(0), eps_adapter%xr(0))
 call eps_adapter%GetEigenVector(value, vector)
 if (.not. ieee_is_nan(real(value))) stop 2
 call eps_adapter%GetEigenVals(values)
 if (size(values) /= 0) stop 3
 call eps_adapter%GetEigenVectors(values, vectors)
 if (size(values) /= 0 .or. size(vectors) /= 0) stop 4
 call pep_adapter%GetEigenVal(value)
 if (.not. ieee_is_nan(real(value))) stop 5
 allocate(pep_adapter%kr(0), pep_adapter%xr(0))
 call pep_adapter%GetEigenVector(value, vector)
 if (.not. ieee_is_nan(real(value))) stop 6
 call pep_adapter%GetEigenVals(values)
 if (size(values) /= 0) stop 7
 call pep_adapter%GetEigenVectors(values, vectors)
 if (size(values) /= 0 .or. size(vectors) /= 0) stop 8
 n = 14
 do j = 1, 3
   call MatCreate(PETSC_COMM_WORLD, matrices(j), ierr)
   call MatSetSizes(matrices(j), PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
   call MatSetType(matrices(j), MATAIJ, ierr)
   call MatSetUp(matrices(j), ierr)
 end do
 do i = 0, n-1
   entry = cmplx(real(i+1, REAL64), 0.001_REAL64*real(i+1, REAL64), REAL64)
   call MatSetValue(matrices(1), i, i, entry, INSERT_VALUES, ierr)
   entry = cmplx(1.0_REAL64, 0.0_REAL64, REAL64)
   call MatSetValue(matrices(2), i, i, entry, INSERT_VALUES, ierr)
   entry = cmplx(0.01_REAL64, 0.0_REAL64, REAL64)
   call MatSetValue(matrices(3), i, i, entry, INSERT_VALUES, ierr)
 end do
 do j = 1, 3
   call MatAssemblyBegin(matrices(j), MAT_FINAL_ASSEMBLY, ierr)
   call MatAssemblyEnd(matrices(j), MAT_FINAL_ASSEMBLY, ierr)
 end do
 ! GetEigenVal-only callers must leave a usable initial vector for later Solve.
 do pass = 1, 5
   call eps_adapter%Solve(2, matrices(1:2), cmplx(2.01_REAL64, 0.0_REAL64, REAL64))
   call eps_adapter%GetEigenVal(value)
   if (.not. ieee_is_finite(real(value))) stop 9
   if (eps_adapter%last_error > 1.0e-8_REAL64) stop 10
 end do
 ! Several eigenpairs converged: largest growth and closest target must differ.
 call eps_adapter%GetEigenVal(default_value)
 call eps_adapter%CopySelectedEigenVector(raw_value, raw_saved, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 103
 if (raw_error>1.e-8_REAL64 .or. raw_error/=eps_adapter%last_error) stop 104
 call VecNorm(raw_saved,NORM_2,saved_norm,ierr)
 if (ierr/=0) stop 105
 eps_adapter%track_target = .true.
 ! Copy must not reselect a different eigenpair after the tracking flag changes.
 call eps_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 106
 call VecScale(raw,cmplx(2.0_REAL64,0.0_REAL64,REAL64),ierr)
 if (ierr/=0) stop 107
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 108
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 109
 eps_adapter%track_target = .true.
 call eps_adapter%GetEigenVal(tracked_value)
 if (abs(real(tracked_value)-2.0_REAL64) > 1.0e-8_REAL64) stop 15
 if (abs(tracked_value-eps_adapter%target) >= abs(default_value-eps_adapter%target)) stop 16
 ! EXPLICIT MAC SELECT: choose the reference mode even when the scalar
 ! nearest-target getter selected another pair. Gates are not relaxed.
 call eps_adapter%SelectEigenVectorByOverlap(raw_saved,default_value, &
      2.0_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 121
 if (selected_overlap<.999999_REAL64 .or. raw_error/=eps_adapter%last_error) stop 122
 if (eps_adapter%overlap_selected_index<=0) stop 123
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 124
 ! The only candidate inside this tight lambda gate has a different vector.
 call eps_adapter%SelectEigenVectorByOverlap(raw_saved,tracked_value, &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 125
 call eps_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 126
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 127
 ! No candidate exists inside a remote tight gate; selection must stay intact.
 call eps_adapter%SelectEigenVectorByOverlap(raw_saved,cmplx(1000._REAL64,0._REAL64,REAL64), &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 128
 call eps_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 129
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 130
 ! Restore the original scalar query to check its unchanged default policy.
 ! End this explicit-MAC fixture before the pre-existing nearest-target solve.
 eps_adapter%overlap_selected_index = 0
 call eps_adapter%GetEigenVal(tracked_value)
 call eps_adapter%GetEigenVector(value, vector)
 if (abs(value-tracked_value) > 1.0e-12_REAL64) stop 17
 if (eps_adapter%last_error > 1.0e-8_REAL64) stop 11
 call eps_adapter%Solve(2, matrices(1:2), cmplx(2.02_REAL64, 0.0_REAL64, REAL64))
 call eps_adapter%GetEigenVal(value)
 if (abs(real(value)-2.0_REAL64) > 1.0e-8_REAL64) stop 18
 ! A copied vector must survive adapter replacement during the next Solve.
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 110
 call VecDestroy(raw_saved,ierr)
 if (ierr/=0) stop 111
 print *, 'EPS default growth selection:', default_value, 'target tracking:', tracked_value
 do pass = 1, 5
   call pep_adapter%Solve(3, matrices, cmplx(-2.01_REAL64, 0.0_REAL64, REAL64))
   call pep_adapter%GetEigenVal(value)
   if (.not. ieee_is_finite(real(value))) stop 12
   if (pep_adapter%last_error > 1.0e-8_REAL64) stop 13
 end do
 call pep_adapter%GetEigenVal(default_value)
 call pep_adapter%CopySelectedEigenVector(raw_value, raw_saved, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 103
 if (raw_error>1.e-8_REAL64 .or. raw_error/=pep_adapter%last_error) stop 104
 call VecNorm(raw_saved,NORM_2,saved_norm,ierr)
 if (ierr/=0) stop 105
 pep_adapter%track_target = .true.
 ! Copy must not reselect a different eigenpair after the tracking flag changes.
 call pep_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 106
 call VecScale(raw,cmplx(2.0_REAL64,0.0_REAL64,REAL64),ierr)
 if (ierr/=0) stop 107
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 108
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 109
 pep_adapter%track_target = .true.
 call pep_adapter%GetEigenVal(tracked_value)
 if (abs(tracked_value-pep_adapter%target) >= abs(default_value-pep_adapter%target)) stop 19
 if (abs(tracked_value-default_value) < 0.1_REAL64) stop 20
 ! EXPLICIT MAC SELECT: choose the reference mode even when the scalar
 ! nearest-target getter selected another pair. Gates are not relaxed.
 call pep_adapter%SelectEigenVectorByOverlap(raw_saved,default_value, &
      2.0_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 121
 if (selected_overlap<.999999_REAL64 .or. raw_error/=pep_adapter%last_error) stop 122
 if (pep_adapter%overlap_selected_index<=0) stop 123
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 124
 ! The only candidate inside this tight lambda gate has a different vector.
 call pep_adapter%SelectEigenVectorByOverlap(raw_saved,tracked_value, &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 125
 call pep_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 126
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 127
 ! No candidate exists inside a remote tight gate; selection must stay intact.
 call pep_adapter%SelectEigenVectorByOverlap(raw_saved,cmplx(1000._REAL64,0._REAL64,REAL64), &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 128
 call pep_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 129
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 130
 ! Restore the original scalar query to check its unchanged default policy.
 ! End this explicit-MAC fixture before the pre-existing nearest-target solve.
 pep_adapter%overlap_selected_index = 0
 call pep_adapter%GetEigenVal(tracked_value)
 call pep_adapter%GetEigenVector(value, vector)
 if (abs(value-tracked_value) > 1.0e-12_REAL64) stop 21
 if (pep_adapter%last_error > 1.0e-8_REAL64) stop 14
 call pep_adapter%Solve(3, matrices, cmplx(-2.02_REAL64, 0.0_REAL64, REAL64))
 call pep_adapter%GetEigenVal(value)
 if (abs(value-tracked_value) > 1.0e-10_REAL64) stop 22
 ! A copied vector must survive adapter replacement during the next Solve.
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 110
 call VecDestroy(raw_saved,ierr)
 if (ierr/=0) stop 111
 print *, 'PEP default growth selection:', default_value, 'target tracking:', tracked_value
 call eps_adapter%Finalize()
 call pep_adapter%Finalize()
 do j = 1, 3
   call MatDestroy(matrices(j), ierr)
 end do
 print *, 'adapter_safety_test PASS: empty getters, repeated solves, residuals, default and target branch selection'
 ! Repeated normalization must release temporary Vec and VecScatter objects.
 call VecCreateSeq(PETSC_COMM_SELF,2000,normalization_vector,ierr)
 if (ierr/=0) stop 112
 call VecSet(normalization_vector,cmplx(1._REAL64,0._REAL64,REAL64),ierr)
 if (ierr/=0) stop 113
 do memory_pass=1,10
   call normal_u(normalization_vector)
 end do
 call PetscMallocGetCurrentUsage(memory_before,ierr)
 if (ierr/=0) stop 114
 do memory_pass=1,300
   call normal_u(normalization_vector)
 end do
 call PetscMallocGetCurrentUsage(memory_after,ierr)
 if (ierr/=0) stop 115
 print *, 'NORMALIZATION MEMORY BEFORE/AFTER/DELTA:',memory_before,memory_after,memory_after-memory_before
 if (memory_after-memory_before>16384._REAL64) stop 116
 call VecDestroy(normalization_vector,ierr)
 if (ierr/=0) stop 117
 print *, 'RAW COPY AND NORMALIZATION MEMORY PASS'
 call SlepcFinalize(ierr)
end program raw_vector_copy_test
