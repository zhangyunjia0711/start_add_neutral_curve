#include <slepc/finclude/slepc.h>
program adapter_safety_test
 use iso_fortran_env, only: REAL64
 use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
 use slepceps
 use slepcpep
 use petsc
 use mod_slepc_adapter
 implicit none
 type(slepceps_type) :: eps_adapter
 type(slepcpep_type) :: pep_adapter
 Mat :: matrices(3)
 Vec :: vector
 PetscErrorCode :: ierr
 PetscScalar :: value, entry, default_value, tracked_value
 PetscScalar, allocatable :: values(:)
 Vec, allocatable :: vectors(:)
 PetscInt :: n, i, j, pass
 call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
 call eps_adapter%Create(1, 1)
 call pep_adapter%Create(1, 1)
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
 eps_adapter%track_target = .true.
 call eps_adapter%GetEigenVal(tracked_value)
 if (abs(real(tracked_value)-2.0_REAL64) > 1.0e-8_REAL64) stop 15
 if (abs(tracked_value-eps_adapter%target) >= abs(default_value-eps_adapter%target)) stop 16
 call eps_adapter%GetEigenVector(value, vector)
 if (abs(value-tracked_value) > 1.0e-12_REAL64) stop 17
 if (eps_adapter%last_error > 1.0e-8_REAL64) stop 11
 call eps_adapter%Solve(2, matrices(1:2), cmplx(2.02_REAL64, 0.0_REAL64, REAL64))
 call eps_adapter%GetEigenVal(value)
 if (abs(real(value)-2.0_REAL64) > 1.0e-8_REAL64) stop 18
 print *, 'EPS default growth selection:', default_value, 'target tracking:', tracked_value
 do pass = 1, 5
   call pep_adapter%Solve(3, matrices, cmplx(-2.01_REAL64, 0.0_REAL64, REAL64))
   call pep_adapter%GetEigenVal(value)
   if (.not. ieee_is_finite(real(value))) stop 12
   if (pep_adapter%last_error > 1.0e-8_REAL64) stop 13
 end do
 call pep_adapter%GetEigenVal(default_value)
 pep_adapter%track_target = .true.
 call pep_adapter%GetEigenVal(tracked_value)
 if (abs(tracked_value-pep_adapter%target) >= abs(default_value-pep_adapter%target)) stop 19
 if (abs(tracked_value-default_value) < 0.1_REAL64) stop 20
 call pep_adapter%GetEigenVector(value, vector)
 if (abs(value-tracked_value) > 1.0e-12_REAL64) stop 21
 if (pep_adapter%last_error > 1.0e-8_REAL64) stop 14
 call pep_adapter%Solve(3, matrices, cmplx(-2.02_REAL64, 0.0_REAL64, REAL64))
 call pep_adapter%GetEigenVal(value)
 if (abs(value-tracked_value) > 1.0e-10_REAL64) stop 22
 print *, 'PEP default growth selection:', default_value, 'target tracking:', tracked_value
 call eps_adapter%Finalize()
 call pep_adapter%Finalize()
 do j = 1, 3
   call MatDestroy(matrices(j), ierr)
 end do
 print *, 'adapter_safety_test PASS: empty getters, repeated solves, residuals, default and target branch selection'
 call SlepcFinalize(ierr)
end program adapter_safety_test
