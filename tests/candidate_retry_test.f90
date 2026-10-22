! Small diagonal QEP: no case files or full neutral scan are required.
! CTest also checks PETSc's final allocation dump and the explicit EPS rejection.
#include <slepc/finclude/slepc.h>
program candidate_retry_test
   use iso_fortran_env, only: REAL64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use slepcpep
   use slepceps
   use petsc
   use mod_slepc_adapter
   implicit none
   type(slepcpep_type) :: adapter
   type(slepceps_type) :: eps_adapter
   Mat :: matrices(3)
   Vec :: reference
   PetscScalar :: entry, eigenvalue, copied_eigenvalue
   PetscReal :: residual, reference_norm
   PetscLogDouble :: memory_before, memory_after
   PetscInt :: n, i, j, nev, ncv, mpd, memory_cycles
   PetscErrorCode :: error
   logical :: valid, reject_eps

   call SlepcInitialize(PETSC_NULL_CHARACTER, error)
   call checked(error)
   call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-reject_eps', reject_eps, valid, error)
   call checked(error)
   if (.not. valid) reject_eps = .false.
   call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-memory_cycles', memory_cycles, valid, error)
   call checked(error)
   if (.not. valid) memory_cycles = 0
   call adapter%Create(1, 2)
   call PEPSetDimensions(adapter%solver, 1, 8, PETSC_DECIDE, error)
   call checked(error)
   adapter%track_target = .true.
   n = 40
   do j = 1, 3
      call MatCreate(PETSC_COMM_WORLD, matrices(j), error)
      call checked(error)
      call MatSetSizes(matrices(j), PETSC_DECIDE, PETSC_DECIDE, n, n, error)
      call checked(error)
      call MatSetType(matrices(j), MATAIJ, error)
      call checked(error)
      call MatSetUp(matrices(j), error)
      call checked(error)
   end do
   do i = 0, n-1
      entry = cmplx(real(i+1, REAL64), .001_REAL64*real(i+1, REAL64), REAL64)
      call MatSetValue(matrices(1), i, i, entry, INSERT_VALUES, error)
      call checked(error)
      call MatSetValue(matrices(2), i, i, cmplx(1._REAL64, 0._REAL64, REAL64), INSERT_VALUES, error)
      call checked(error)
      call MatSetValue(matrices(3), i, i, cmplx(.01_REAL64, 0._REAL64, REAL64), INSERT_VALUES, error)
      call checked(error)
   end do
   do j = 1, 3
      call MatAssemblyBegin(matrices(j), MAT_FINAL_ASSEMBLY, error)
      call checked(error)
      call MatAssemblyEnd(matrices(j), MAT_FINAL_ASSEMBLY, error)
      call checked(error)
   end do
   call solve_once()
   call adapter%GetEigenVal(eigenvalue)
   if (.not. ieee_is_finite(abs(eigenvalue))) stop 201
   call adapter%CopySelectedEigenVector(copied_eigenvalue, reference, residual, valid)
   if (.not. valid .or. abs(copied_eigenvalue-eigenvalue)>1.e-12_REAL64) stop 202
   if (reject_eps) then
      call eps_adapter%Create(1, 1)
      call eps_adapter%PrepareCandidateRetry(8, reference)
      stop 203
   end if
   call adapter%PrepareCandidateRetry(8, reference)
   ! Destroy the caller's handle before Solve: the adapter must own a real copy.
   call VecDestroy(reference, error)
   call checked(error)
   call VecNorm(adapter%retry_reference, NORM_2, reference_norm, error)
   call checked(error)
   if (abs(reference_norm-1._REAL64)>1.e-10_REAL64) stop 204
   call solve_once()
   if (adapter%retry_count/=0 .or. adapter%has_retry_reference) stop 205
   if (.not. adapter%restore_dimensions) stop 206
   call PEPGetDimensions(adapter%solver, nev, ncv, mpd, error)
   call checked(error)
   if (nev/=8 .or. ncv/=24) stop 207
   call adapter%GetEigenVal(copied_eigenvalue)
   if (abs(copied_eigenvalue-eigenvalue)>1.e-9_REAL64) stop 208
   call solve_once()
   call PEPGetDimensions(adapter%solver, nev, ncv, mpd, error)
   call checked(error)
   if (nev/=1 .or. ncv/=8 .or. adapter%restore_dimensions) stop 209
   call adapter%GetEigenVal(copied_eigenvalue)
   call adapter%CopySelectedEigenVector(copied_eigenvalue, reference, residual, valid)
   if (.not. valid) stop 210
   call adapter%PrepareCandidateRetry(8, reference)
   call adapter%PrepareCandidateRetry(12, reference)
   call solve_once()
   call PEPGetDimensions(adapter%solver, nev, ncv, mpd, error)
   call checked(error)
   if (nev/=12 .or. ncv/=32) stop 211
   if (adapter%retry_count/=0 .or. adapter%has_retry_reference) stop 212
   call solve_once()
   call PEPGetDimensions(adapter%solver, nev, ncv, mpd, error)
   call checked(error)
   if (nev/=1 .or. ncv/=8) stop 213
   if (memory_cycles>0) then
      call PetscMallocGetCurrentUsage(memory_before, error)
      call checked(error)
      do i = 1, memory_cycles
         call solve_once()
      end do
      call PetscMallocGetCurrentUsage(memory_after, error)
      call checked(error)
      print *, 'NORMAL SOLVE MEMORY CYCLES/BEFORE/AFTER/DELTA:', memory_cycles, &
               memory_before, memory_after, memory_after-memory_before
      if (memory_after-memory_before>16384._REAL64) stop 217
      call PetscMallocGetCurrentUsage(memory_before, error)
      call checked(error)
      do i = 1, memory_cycles
         call adapter%PrepareCandidateRetry(8, reference)
         call solve_once()
         call solve_once()
      end do
      call PetscMallocGetCurrentUsage(memory_after, error)
      call checked(error)
      print *, 'RETRY PLUS RESTORE MEMORY CYCLES/BEFORE/AFTER/DELTA:', memory_cycles, &
               memory_before, memory_after, memory_after-memory_before
      if (memory_after-memory_before>16384._REAL64) stop 218
   end if
   ! Repeated replacement must release the old owned reference, even if no
   ! Solve consumes it. Finalize must also consume the final pending request.
   do i = 1, 10
      call adapter%PrepareCandidateRetry(8, reference)
   end do
   call PetscMallocGetCurrentUsage(memory_before, error)
   call checked(error)
   do i = 1, 2000
      call adapter%PrepareCandidateRetry(8, reference)
   end do
   call PetscMallocGetCurrentUsage(memory_after, error)
   call checked(error)
   print *, 'RETRY MEMORY BEFORE/AFTER/DELTA:', memory_before, memory_after, memory_after-memory_before
   if (memory_after-memory_before>1.e4_REAL64) stop 214
   call adapter%Finalize()
   if (adapter%retry_count/=0 .or. adapter%has_retry_reference .or. adapter%restore_dimensions) stop 215
   call VecNorm(reference, NORM_2, reference_norm, error)
   call checked(error)
   if (abs(reference_norm-1._REAL64)>1.e-10_REAL64) stop 216
   call VecDestroy(reference, error)
   call checked(error)
   do j = 1, 3
      call MatDestroy(matrices(j), error)
      call checked(error)
   end do
   print *, 'CANDIDATE RETRY PASS: owned copy, one-shot dimensions, replacement, restore and pending finalize'
   call SlepcFinalize(error)
   call checked(error)
contains
   subroutine solve_once()
      call adapter%Solve(3, matrices, cmplx(-2.01_REAL64, 0._REAL64, REAL64))
   end subroutine solve_once
   subroutine checked(code)
      PetscErrorCode, intent(in) :: code
      if (code/=0) stop 299
   end subroutine checked
end program candidate_retry_test
