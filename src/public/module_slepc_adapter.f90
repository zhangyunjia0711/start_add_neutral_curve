!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-23 11:16:39
!> @LastEditTime: 2022-12-12 16:18:53
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/module_slepc_adapter.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <slepc/finclude/slepc.h>
#include <petsc/finclude/petsc.h>

module mod_slepc_adapter

   use iso_fortran_env, only: REAL64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   ! use global_parameters, only: lst_mode
   use slepcpep
   use slepceps
   use petsc

   implicit none

   public slepcpep_type
   public slepceps_type
   public slepc_type
   private

#include "../include/start.f90"

   PetscErrorCode :: ierr

   !> SLEPC通用求解类
   type, abstract :: slepc_type

      PetscInt :: m !<问题矩阵维度
      Vec, allocatable :: xr(:)   !<特征向量
      PetscScalar, allocatable :: kr(:) !特征值
      Vec :: initialguess(1) !初始猜测向量
      integer, allocatable :: evIndex(:) !<特征值序号
      PetscScalar :: target !< 目标值
      ST :: st
      PC :: pc
      KSP :: ksp
      integer :: lst_mode
      real(REAL64) :: last_error = huge(1.0_REAL64) !< Absolute residual of the selected eigenpair
      logical :: track_target = .false. !< Select the eigenpair nearest the complex target for continuation
      integer :: overlap_selected_index = 0 !< Explicit MAC choice for the next initial space only

   contains

      

      procedure(create_interface), deferred :: Create
      procedure(solve_interface), deferred :: Solve
      procedure(GetEigenVal_interface), deferred :: GetEigenVal
      procedure(GetEigenVals_interface), deferred :: GetEigenVals
      procedure(GetEigenVector_interface), deferred :: GetEigenVector
      procedure(GetEigenVectors_interface), deferred :: GetEigenVectors
      procedure(Finalize_interface), deferred :: Finalize
      procedure :: CopySelectedEigenVector
      procedure :: SelectEigenVectorByOverlap
      procedure :: PrepareCandidateRetry

   end type slepc_type

   abstract interface

      subroutine create_interface(this, nev, lst_mode)

         import slepc_type
         class(slepc_type), intent(inout) :: this
         PetscInt, intent(in) :: nev
         integer, intent(in) :: lst_mode

      end subroutine create_interface

      ! subroutine solve_interface(this, nmat, mat, target_input)
      !   import slepc_type
      !   import abstract_matrix_object_type
      subroutine solve_interface(this, nmat, mat, target_input)

         import slepc_type
         import REAL64
         import tmat
         class(slepc_type), intent(inout) :: this
         integer, intent(in) :: nmat
         Mat, intent(in) :: mat(nmat)
         PetscScalar, intent(in) :: target_input

      end subroutine solve_interface

      !> 取特征值
      subroutine GetEigenVal_interface(this, eigenvalue)

         import slepc_type
         class(slepc_type), intent(inout) :: this
         PetscScalar, intent(out) :: eigenvalue

      end subroutine GetEigenVal_interface

      !> 取特征值
      subroutine GetEigenVals_interface(this, eigenvalue)

         import slepc_type
         class(slepc_type), intent(inout) :: this
         PetscScalar, intent(out), allocatable :: eigenvalue(:)

      end subroutine GetEigenVals_interface

      !> 取特征向量
      subroutine GetEigenVector_interface(this, eigenvalue, eigenVector)

         import slepc_type
         import REAL64
         import tVec
         class(slepc_type), intent(inout) :: this
         PetscScalar, intent(out) :: eigenvalue
         Vec, intent(out) :: eigenVector

      end subroutine GetEigenVector_interface

      !> 取特征向量
      subroutine GetEigenVectors_interface(this, eigenvalue, eigenVector)

         import slepc_type
         import REAL64
         import tVec
         class(slepc_type), intent(inout) :: this
         PetscScalar, intent(out), allocatable :: eigenvalue(:)
         Vec, intent(out), allocatable :: eigenVector(:)

      end subroutine GetEigenVectors_interface

      subroutine finalize_interface(this)

         import slepc_type
         class(slepc_type), intent(inout) :: this

      end subroutine finalize_interface

   end interface

   type, extends(slepc_type) :: slepceps_type
      EPS :: solver
   contains
      procedure :: Create => EPSSolverCreate !<创建求解器
      procedure :: Solve => EPSSolverSolve  !<求解问题
      procedure :: GetEigenVal => EPSSolverGetEigenVal !<取特征值
      procedure :: GetEigenVals => EPSSolverGetEigenVals !<取特征值
      procedure :: GetEigenVector => EPSSolverGetEigenVector !<取特征向量
      procedure :: GetEigenVectors => EPSSolverGetEigenVectors !<取特征向量
      procedure :: Finalize => EPSFinalize

   end type slepceps_type

   !> SLEPC多项式算子特征值问题求解器
   type, extends(slepc_type) :: slepcpep_type

      PEP :: solver !<多项式算子特征值问题求解器
      integer :: retry_count = 0
      Vec :: retry_reference !< Owned copy, valid only with has_retry_reference
      logical :: has_retry_reference = .false.
      logical :: restore_dimensions = .false.
      PetscInt :: normal_nev, normal_ncv, normal_mpd

   contains

      procedure :: Create => PEPSolverCreate !<创建求解器
      procedure :: Solve => PEPSolverSolve  !<求解问题
      procedure :: GetEigenVal => PEPSolverGetEigenVal !<取特征值
      procedure :: GetEigenVector => PEPSolverGetEigenVector !<取特征向量
      procedure :: GetEigenVals => PEPSolverGetEigenVals !<取特征值
      procedure :: GetEigenVectors => PEPSolverGetEigenVectors !<取特征向量
      procedure :: Finalize => PEPFinalize

   end type slepcpep_type

Contains

   !> 创建求解器
   subroutine PEPSolverCreate(this, nev, lst_mode)

      implicit none
      class(slepcpep_type), intent(inout) :: this
      integer, intent(in) :: nev !<求解特征值个数
      integer, intent(in) :: lst_mode

      call PEPCreate(PETSC_COMM_WORLD, this%solver, ierr)
      call CheckSolverConfiguration(ierr)
      this%lst_mode = lst_mode
      call PEPSetWhichEigenpairs(this%solver, PEP_TARGET_MAGNITUDE, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPSetType(this%solver, PEPLINEAR, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPSetConvergenceTest(this%solver, PEP_CONV_NORM, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPGetST(this%solver, this%st, ierr)
      call CheckSolverConfiguration(ierr)
      call SetDefaultDirectSolver(this)
      call PEPSetTolerances(this%solver, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call CheckSolverConfiguration(ierr)
      call PEPSetDimensions(this%solver, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
      call CheckSolverConfiguration(ierr)
      ! print*, nev, 'eigenvalues will be computed.'

   end subroutine PEPSolverCreate

   !> 创建求解器
   subroutine EPSSolverCreate(this, nev, lst_mode)

      implicit none
      class(slepceps_type), intent(inout) :: this
      integer, intent(in) :: nev !<求解特征值个数
      integer, intent(in) :: lst_mode

      call EPSCreate(PETSC_COMM_WORLD, this%solver, ierr)
      call CheckSolverConfiguration(ierr)
      this%lst_mode = lst_mode
      call EPSGetST(this%solver, this%st, ierr)
      call CheckSolverConfiguration(ierr)
      call EPSSetWhichEigenpairs(this%solver, EPS_TARGET_MAGNITUDE, ierr)
      call CheckSolverConfiguration(ierr)
      call SetDefaultDirectSolver(this)
      call EPSSetTolerances(this%solver, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call CheckSolverConfiguration(ierr)
      call EPSSetDimensions(this%solver, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
      call CheckSolverConfiguration(ierr)
      ! print*, nev, 'eigenvalues will be computed.'

   end subroutine EPSSolverCreate

   !> Set defaults only while the new ST/PC has not been set up. Runtime options
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
      PetscErrorCode :: report_error
      if (error /= 0) then
         ! PetscErrorF writes its first argument. Never give it an intent(in)
         ! alias, which may ultimately refer to a read-only PETSc constant.
         report_error = error
         call PetscErrorF(report_error, __LINE__, "module_slepc_adapter.f90")
         call MPIU_Abort(PETSC_COMM_WORLD, error)
      end if
   end subroutine CheckSolverConfiguration

   !> 求解特征值
   subroutine PEPSolverSolve(this, nmat, mat, target_input)

      class(slepcpep_type), intent(inout) :: this
      integer, intent(in) :: nmat
      Mat, intent(in) :: mat(nmat)
      complex(REAL64), intent(in) :: target_input
      PetscScalar :: target_value
      PetscInt :: nconv, nrow, guess_size, retry_ncv
      PetscErrorCode :: solve_error
      Vec :: retry_initialguess(1)
      integer :: i, guess_index
      logical :: has_initialguess
      PetscScalar :: ki

      this%last_error = huge(1.0_REAL64)
      has_initialguess = .false.
      ! PEPLINEAR's implicit setup allocates shell matrices and work vectors.
      ! Reset before rebuilding so repeated dimension/target changes release
      ! that setup. Our independent xr copies and one-shot reference survive.
      call PEPReset(this%solver, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      call PEPSetProblemType(this%solver, PEP_GENERAL, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      call PEPSetOperators(this%solver, nmat, mat, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      call MatGetSize(mat(1), this%m, nrow, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      target_value = target_input
      this%target = target_input
      call PEPSetTarget(this%solver, target_value, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      call PEPSetType(this%solver, PEPLINEAR, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      ! A previous expansion must not remain sticky when options are absent.
      ! Production still supplies explicit -pep_nev 1 -pep_ncv 8 on normal solves.
      if (this%restore_dimensions) then
         call PEPSetDimensions(this%solver, this%normal_nev, this%normal_ncv, this%normal_mpd, ierr)
         call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
         this%restore_dimensions = .false.
      end if
      call PEPSetFromOptions(this%solver, ierr)
      call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
      if (this%retry_count > 0) then
         if (.not. this%has_retry_reference) then
            print *, 'Candidate retry rejected: no owned reference vector.'
            call CheckPEPSolveConfiguration(this, PETSC_ERR_ARG_WRONGSTATE, has_initialguess)
         end if
         call VecGetSize(this%retry_reference, guess_size, ierr)
         call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
         if (guess_size /= this%m) then
            print *, 'Candidate retry rejected: reference/matrix size mismatch.', guess_size, this%m
            call CheckPEPSolveConfiguration(this, PETSC_ERR_ARG_SIZ, has_initialguess)
         end if
         call PEPGetDimensions(this%solver, this%normal_nev, this%normal_ncv, this%normal_mpd, ierr)
         call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
         retry_ncv = max(24, 2*this%retry_count+8)
         call PEPSetDimensions(this%solver, this%retry_count, retry_ncv, PETSC_DECIDE, ierr)
         call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
         this%restore_dimensions = .true.
         retry_initialguess(1) = this%retry_reference
         call PEPSetInitialSpace(this%solver, 1, retry_initialguess, ierr)
         call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
         print *, 'Candidate retry consuming owned reference: nev=', this%retry_count, ' ncv=', retry_ncv
      else if (allocated(this%xr)) then
         if (size(this%xr) > 0) then
            guess_index = SelectedGuessIndex(this)
            call VecGetSize(this%xr(guess_index), guess_size, ierr)
            call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
            if (guess_size == this%m) then
               call VecDuplicate(this%xr(guess_index), this%initialguess(1), ierr)
               call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
               has_initialguess = .true.
               call VecCopy(this%xr(guess_index), this%initialguess(1), ierr)
               call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
               call PEPSetInitialSpace(this%solver, 1, this%initialguess(1), ierr)
               call CheckPEPSolveConfiguration(this, ierr, has_initialguess)
            end if
         end if
      end if
      call PEPSolve(this%Solver, ierr)
      solve_error = ierr
      nconv = 0
      if (solve_error == 0) then
         call PEPGetConverged(this%Solver, nconv, ierr)
         solve_error = ierr
         if (solve_error /= 0) nconv = 0
      end if
      ! Both successful and failed PEPSolve returns consume this one-shot copy.
      call ClearCandidateRetry(this)
      if (has_initialguess) then
         call VecDestroy(this%initialguess(1), ierr)
         call CheckSolverConfiguration(ierr)
      end if
      ! Preserve the existing empty-result failure behavior, with an explicit
      ! diagnostic and with all one-shot resources already released.
      if (solve_error /= 0) print *, 'PEP solve returned an error; no candidates retained:', solve_error
      call ReleaseEigenpairs(this)
      allocate (this%xr(nconv), this%kr(nconv))
      do i = 0, nconv - 1
         call MatCreateVecs(mat(1), this%xr(i + 1), PETSC_NULL_VEC, ierr)
         call CheckSolverConfiguration(ierr)
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, this%xr(i + 1), PETSC_NULL_VEC, ierr)
         call CheckSolverConfiguration(ierr)
      end do

      ! block
      !  Integer :: i
      !  integer :: rank=0
      !  PetscScalar :: kr, ki
      !  Vec :: xr, xi
      !  PetscReal :: error
      !  PetscInt :: nev, maxit, its
      !  EPSType :: tname
      !  PetscReal :: tol
      !
      !  rank=0
      !
      !  call PEPGetIterationNumber(this%Solver, its, ierr)
      !  if(rank == 0) write(*, *)'iterations number:', its
      !  call PEPGetType(this%Solver, tname, ierr)
      !  if(rank==0)write(*, *)'solver type:', tname
      !  call PEPGetDimensions(this%Solver, nev, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      !  if(rank==0)write(*, *)'eigenvalues number:', nev
      !  call PEPGetTolerances(this%Solver, tol, maxit, ierr)
      !  if(rank==0)write(*, *)'tol=:', tol, 'maxit=:', maxit
      ! if(nconv>0)then
      !   print*, 'k', '      ||Ax-kx||/||kx||'
      ! endif
      ! do i=0, nconv-1
! comments by Jianxin 2018-11-05
      ! call PEPGetEigenPair(this%solver, i, kr, ki, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
      ! call PEPComputeError(this%Solver, i, EPS_ERROR_RELATIVE, error, ierr)
      ! write(*, '(I6, 4F20.7)') i, kr, error
! comments by Jianxin 2018-11-05
      ! enddo

      ! endblock

   end subroutine PEPSolverSolve

   !> 求解特征值
   subroutine EPSSolverSolve(this, nmat, mat, target_input)

      class(slepceps_type), intent(inout) :: this
      integer, intent(in) :: nmat
      Mat, intent(in) :: mat(nmat)
      complex(REAL64), intent(in) :: target_input
      PetscScalar :: target_value
      PetscInt :: nconv, nrow, guess_size
      integer :: i, guess_index
      logical :: has_initialguess
      PetscScalar :: ki

      this%last_error = huge(1.0_REAL64)
      has_initialguess = .false.
      call EPSSetProblemType(this%solver, EPS_GNHEP, ierr)
      call EPSSetType(this%solver, EPSARNOLDI, ierr)
      call EPSSetOperators(this%solver, mat(1), mat(2), ierr)
      call MatGetSize(mat(1), this%m, nrow, ierr)
      target_value = target_input
      this%target = target_input
      ! print*, 'The eps solver is selected!'
      ! print*, 'The target_value is', target_value
      call EPSSetTarget(this%solver, target_value, ierr)
      ! Direct-solver defaults were installed once in Create; preserve runtime overrides.
      call EPSSetFromOptions(this%solver, ierr)
      call CheckSolverConfiguration(ierr)
      if (allocated(this%xr)) then
         if (size(this%xr) > 0) then
            guess_index = SelectedGuessIndex(this)
            call VecGetSize(this%xr(guess_index), guess_size, ierr)
            if (guess_size == this%m) then
               call VecDuplicate(this%xr(guess_index), this%initialguess(1), ierr)
               call VecCopy(this%xr(guess_index), this%initialguess(1), ierr)
               call EPSSetInitialSpace(this%solver, 1, this%initialguess(1), ierr)
               has_initialguess = .true.
            end if
         end if
      end if
      ! print*, 'problem is setting....'
      ! print*, this%target

      ! block
      !
      !   PetscViewer::Viewer
      !
      !   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'LST_matrix1.mat', FILE_MODE_WRITE, &
      !   Viewer, ierr)
      !   call MatView(mat(1), viewer, ierr)
      !   call PetscViewerDestroy(Viewer, ierr)
      !   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'LST_matrix2.mat', FILE_MODE_WRITE, &
      !   Viewer, ierr)
      !   call MatView(mat(2), viewer, ierr)
      !   call PetscViewerDestroy(Viewer, ierr)
      !
      !   ! call SLEPCFinalize(ierr)
      !
      !   ! stop
      !
      !
      ! endblock

      call EPSSolve(this%Solver, ierr)
      nconv = 0
      if (ierr == 0) call EPSGetConverged(this%Solver, nconv, ierr)
      if (ierr /= 0) nconv = 0
      if (has_initialguess) call VecDestroy(this%initialguess(1), ierr)
      ! Destroy PETSc objects before replacing their Fortran handles.
      call ReleaseEigenpairs(this)
      allocate (this%xr(nconv), this%kr(nconv))
      do i = 0, nconv - 1
         call MatCreateVecs(mat(1), this%xr(i + 1), PETSC_NULL_VEC, ierr)
         ! GetEigenVal alone must also leave a valid vector for the next Solve.
         call EPSGetEigenPair(this%solver, i, this%kr(i + 1), ki, this%xr(i + 1), PETSC_NULL_VEC, ierr)
      end do

!     block
!      Integer :: i
!      integer :: rank=0
!      PetscScalar :: kr, ki
!      Vec :: xr, xi
!      PetscReal :: error
!      PetscInt :: nev, maxit, its
!      EPSType :: tname
!      PetscReal :: tol
!
!      rank=0
!
!      call EPSGetIterationNumber(this%Solver, its, ierr)
!      if(rank == 0) write(*, *)'iterations number:', its
!      call EPSGetType(this%Solver, tname, ierr)
!      if(rank==0)write(*, *)'solver type:', tname
!      call EPSGetDimensions(this%Solver, nev, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
!      if(rank==0)write(*, *)'eigenvalues number:', nev
!      call EPSGetTolerances(this%Solver, tol, maxit, ierr)
!      if(rank==0)write(*, *)'tol=:', tol, 'maxit=:', maxit
!     if(nconv>0)then
!       print*, 'k', '      ||Ax-kx||/||kx||'
!     endif
!     do i=0, nconv-1
! ! comments by JK 2018-11-05
!       call EPSGetEigenvalue(this%solver, i, kr, ki, ierr)
!       call EPSComputeError(this%Solver, i, EPS_ERROR_RELATIVE, error, ierr)
!       write(*, '(I6, 6F20.7)') i, kr, error
! ! comments by JK 2018-11-05
!     enddo
!
!     endblock

   end subroutine EPSSolverSolve

   !> 取特征值
   subroutine PEPSolverGetEigenVal(this, eigenvalue)

      use mod_sort
      implicit none
      class(slepcpep_type), intent(inout) :: this
      PetscScalar, intent(out) :: eigenvalue
      PetscScalar :: ki
      integer :: nconv, i
      integer, parameter :: SPATIAL = 1, TEMPORAL = 2, TEMPORAL_PERIODIC_EPS = 3

      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      if (nconv == 0) then
         eigenvalue = InvalidEigenvalue()
         return
      end if

      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

      print *, 'ncov=', nconv

      do i = 0, nconv - 1
         this%evIndex(i + 1) = i + 1
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
         call PEPComputeError(this%solver, i, EPS_ERROR_ABSOLUTE, error, ierr)
         print *, i, this%kr(i+1), 'The error is', error
         !  eigenvalue(i+1)=this%kr(i+1)
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      ! Keep the default growth-rate ordering unless continuation explicitly requests a target branch.
      if (this%track_target) this%evIndex = .sort.(-abs(this%kr - this%target))

      eigenvalue = this%kr(this%evIndex(nconv))
      call PEPComputeError(this%solver, this%evIndex(nconv) - 1, PEP_ERROR_ABSOLUTE, error, ierr)
      if (ierr == 0) this%last_error = error
      print *, 'pep_get_eigenvalue...'

    !!!Undo_set_initial_space

   end subroutine PEPSolverGetEigenVal

   !> 取特征值
   subroutine PEPSolverGetEigenVals(this, eigenvalue)

      use mod_sort
      implicit none
      class(slepcpep_type), intent(inout) :: this
      PetscScalar, intent(out), allocatable :: eigenvalue(:)
      PetscScalar :: ki
      integer :: nconv, i
      integer, parameter :: SPATIAL = 1, TEMPORAL = 2, TEMPORAL_PERIODIC_EPS = 3
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)

      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))
      if (nconv == 0) return

      do i = 0, nconv - 1
         this%evIndex(i + 1) = i + 1
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
         call PEPComputeError(this%solver, i, EPS_ERROR_ABSOLUTE, error, ierr)
         print *, i, this%kr(i+1), 'The error is', error
         !  eigenvalue(i+1)=this%kr(i+1)
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      do i = 1, nconv
         eigenvalue(i) = this%kr(this%evIndex(nconv + 1 - i))
      end do
      print *, 'pep_get_eigenvalues...'

    !!!Undo_set_initial_space

   end subroutine PEPSolverGetEigenVals

   !> 取特征值
   subroutine EPSSolverGetEigenVal(this, eigenvalue)

      use mod_sort
      implicit none
      class(slepceps_type), intent(inout) :: this
      PetscScalar, intent(out) :: eigenvalue
      PetscScalar :: ki
      integer :: nconv, i
      integer, parameter :: SPATIAL = 1, TEMPORAL = 2, TEMPORAL_PERIODIC_EPS = 3
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      if (nconv == 0) then
         eigenvalue = InvalidEigenvalue()
         return
      end if

      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

      do i = 0, nconv - 1
         this%evIndex(i + 1) = i + 1
         call EPSGetEigenPair(this%solver, i, this%kr(i + 1), ki, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
         call EPSComputeError(this%solver, i, EPS_ERROR_ABSOLUTE, error, ierr)
         print *, i, this%kr(i+1), 'The error is', error
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      ! Keep the default growth-rate ordering unless continuation explicitly requests a target branch.
      if (this%track_target) this%evIndex = .sort.(-abs(this%kr - this%target))

      eigenvalue = this%kr(this%evIndex(nconv))
      call EPSComputeError(this%solver, this%evIndex(nconv) - 1, EPS_ERROR_ABSOLUTE, error, ierr)
      if (ierr == 0) this%last_error = error

    !!!Undo_set_initial_space

   end subroutine EPSSolverGetEigenVal

   !> 取特征值
   subroutine EPSSolverGetEigenVals(this, eigenvalue)

      use mod_sort
      implicit none
      class(slepceps_type), intent(inout) :: this
      PetscScalar, intent(out), allocatable :: eigenvalue(:)
      PetscScalar :: ki
      integer :: nconv, i
      integer, parameter :: SPATIAL = 1, TEMPORAL = 2, TEMPORAL_PERIODIC_EPS = 3
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)

      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))
      if (nconv == 0) return

      do i = 0, nconv - 1
         this%evIndex(i + 1) = i + 1
         call EPSGetEigenPair(this%solver, i, this%kr(i + 1), ki, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
         call EPSComputeError(this%solver, i, EPS_ERROR_ABSOLUTE, error, ierr)
         print *, i, this%kr(i+1), 'The error is', error
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      do i = 1, nconv
         eigenvalue(i) = this%kr(this%evIndex(nconv + 1 - i))
      end do
    !!!Undo_set_initial_space

   end subroutine EPSSolverGetEigenVals

   !> 取特征向量
   subroutine PEPSolverGetEigenVector(this, eigenvalue, eigenVector)

      use mod_sort
      implicit none
      class(slepcpep_type), intent(inout) :: this
      PetscScalar, intent(out) :: eigenvalue
      Vec, intent(out) :: eigenVector
      integer :: i, nconv
      PetscScalar :: ki
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      if (nconv == 0) then
         eigenvalue = InvalidEigenvalue()
         eigenVector = PETSC_NULL_VEC
         return
      end if
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

      do i = 0, nconv - 1
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, this%xr(i + 1), PETSC_NULL_VEC, ierr)
         ! call PEPComputeError(this%solver, i, PEP_ERROR_ABSOLUTE, error, ierr)
         ! print *, i, this%kr(i+1),'The error is', error
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      ! Keep the default growth-rate ordering unless continuation explicitly requests a target branch.
      if (this%track_target) this%evIndex = .sort.(-abs(this%kr - this%target))

      eigenvalue = this%kr(this%evIndex(nconv))
      call PEPComputeError(this%solver, this%evIndex(nconv) - 1, PEP_ERROR_ABSOLUTE, error, ierr)
      if (ierr == 0) this%last_error = error
      eigenVector = this%xr(this%evIndex(nconv))

      print *, 'pep_get_eigenvector...'

      ! do i=1, nconv
      !   call VecGetArrayF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      !   eigenVector(i-1)=tmp
      !   call VecRestoreArrayReadF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      ! enddo

   end subroutine PEPSolverGetEigenVector

   !> 取特征向量
   subroutine PEPSolverGetEigenVectors(this, eigenvalue, eigenVector)

      use mod_sort
      implicit none
      class(slepcpep_type), intent(inout) :: this
      PetscScalar, intent(out), allocatable :: eigenvalue(:)
      Vec, intent(out), allocatable :: eigenVector(:)
      integer :: i, nconv
      PetscScalar :: ki
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      allocate (eigenVector(nconv))
      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))
      if (nconv == 0) return

      do i = 0, nconv - 1
         call PEPGetEigenPair(this%solver, i, this%kr(i + 1), ki, this%xr(i + 1), PETSC_NULL_VEC, ierr)
         ! call PEPComputeError(this%solver, i, PEP_ERROR_ABSOLUTE, error, ierr)
         ! print *, i, this%kr(i+1), 'The error is', error
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      do i = 1, nconv
         eigenvalue(i) = this%kr(this%evIndex(nconv + 1 - i))
         eigenVector(i) = this%xr(this%evIndex(nconv + 1 - i))
      end do

      print *, 'pep_get_eigenvectors...'

      ! do i=1, nconv
      !   call VecGetArrayF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      !   eigenVector(i-1)=tmp
      !   call VecRestoreArrayReadF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      ! enddo

   end subroutine PEPSolverGetEigenVectors

   !> 取特征向量
   subroutine EPSSolverGetEigenVector(this, eigenvalue, eigenVector)

      use mod_sort
      implicit none
      class(slepceps_type), intent(inout) :: this
      Vec, intent(out) :: eigenVector
      PetscScalar, intent(out) :: eigenvalue
      integer :: i, nconv
      PetscReal :: error

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      if (nconv == 0) then
         eigenvalue = InvalidEigenvalue()
         eigenVector = PETSC_NULL_VEC
         return
      end if
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

      do i = 0, nconv - 1
         call EPSGetEigenVector(this%solver, i, this%xr(i + 1), PETSC_NULL_VEC, ierr)
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      ! Keep the default growth-rate ordering unless continuation explicitly requests a target branch.
      if (this%track_target) this%evIndex = .sort.(-abs(this%kr - this%target))

      eigenvalue = this%kr(this%evIndex(nconv))
      call EPSComputeError(this%solver, this%evIndex(nconv) - 1, EPS_ERROR_ABSOLUTE, error, ierr)
      if (ierr == 0) this%last_error = error
      eigenVector = this%xr(this%evIndex(nconv))

   end subroutine EPSSolverGetEigenVector

   !> 取特征向量
   subroutine EPSSolverGetEigenVectors(this, eigenvalue, eigenVector)

      use mod_sort
      implicit none
      class(slepceps_type), intent(inout) :: this
      PetscScalar, intent(out), allocatable :: eigenvalue(:)
      Vec, intent(out), allocatable :: eigenVector(:)
      integer :: i, nconv

      this%last_error = huge(1.0_REAL64)
      nconv = ConvergedCount(this)
      allocate (eigenVector(nconv))
      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))
      if (nconv == 0) return

      do i = 0, nconv - 1
         call EPSGetEigenVector(this%solver, i, this%xr(i + 1), PETSC_NULL_VEC, ierr)
      end do

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         this%evIndex = .sort.aimag(this%kr)
      case (LST_SPATIAL_MODE)
         this%evIndex = .sort.aimag(-this%kr)
      end select

      do i = 1, nconv
         eigenvalue(i) = this%kr(this%evIndex(nconv + 1 - i))
         eigenVector(i) = this%xr(this%evIndex(nconv + 1 - i))
      end do

      ! do i=1, nconv
      !   call VecGetArrayF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      !   eigenVector(i-1)=tmp
      !   call VecRestoreArrayReadF90(this%xr(this%evIndex(nconv+1-i)), tmp, ierr)
      ! enddo

   end subroutine EPSSolverGetEigenVectors

   !> Copy the eigenvector that GetEigenVal actually selected. Do not call a
   !> sorting getter here: the tracking policy may have changed since selection.
   !> The caller owns eigenVector and must destroy it when valid is true.
   subroutine CopySelectedEigenVector(this, eigenvalue, eigenVector, error, valid)
      class(slepc_type), intent(inout) :: this
      PetscScalar, intent(out) :: eigenvalue
      Vec, intent(out) :: eigenVector
      real(REAL64), intent(out) :: error
      logical, intent(out) :: valid
      integer :: selected, nconv
      Vec :: created_vector
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
      ! Create into a normal local handle; the invalid-return null sentinel is
      ! an input sentinel and must not be passed as a PETSc output argument.
      call VecDuplicate(this%xr(selected), created_vector, ierr)
      call CheckSolverConfiguration(ierr)
      call VecCopy(this%xr(selected), created_vector, ierr)
      call CheckSolverConfiguration(ierr)
      eigenVector = created_vector
      eigenvalue = this%kr(selected)
      error = this%last_error
      valid = .true.
   end subroutine CopySelectedEigenVector

   !> Explicitly select among this Solve's computed candidates by normalized
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

   !> Release one-shot resources before a fatal configuration failure. The
   !> separate cleanup error variable preserves the original PETSc error.
   subroutine CheckPEPSolveConfiguration(this, error, has_initialguess)
      class(slepcpep_type), intent(inout) :: this
      PetscErrorCode, intent(in) :: error
      logical, intent(in) :: has_initialguess
      PetscErrorCode :: cleanup_error
      if (error == 0) return
      call ClearCandidateRetry(this)
      if (has_initialguess) then
         call VecDestroy(this%initialguess(1), cleanup_error)
         call CheckSolverConfiguration(cleanup_error)
      end if
      call CheckSolverConfiguration(error)
   end subroutine CheckPEPSolveConfiguration

   !> Prepare exactly one expanded spatial PEP solve. The caller retains its
   !> reference; this object stores a true owned copy, replaced safely on repeat.
   !> Unsupported EPS/temporal use is an explicit configuration error.
   subroutine PrepareCandidateRetry(this, count, reference_vector)
      class(slepc_type), intent(inout) :: this
      integer, intent(in) :: count
      Vec, intent(in) :: reference_vector
      Vec :: created_vector
      PetscReal :: reference_norm
      PetscErrorCode :: copy_error, cleanup_error
      select type (this)
      type is (slepcpep_type)
         if (this%lst_mode /= LST_SPATIAL_MODE) then
            print *, 'Candidate retry supports spatial PEP only.'
            call CheckPEPSolveConfiguration(this, PETSC_ERR_SUP, .false.)
         end if
         if (count <= 0 .or. count > (huge(count)-8)/2) then
            print *, 'Candidate retry rejected: invalid requested count.', count
            call CheckPEPSolveConfiguration(this, PETSC_ERR_ARG_OUTOFRANGE, .false.)
         end if
         if (PetscObjectIsNull(reference_vector)) then
            print *, 'Candidate retry rejected: null reference.'
            call CheckPEPSolveConfiguration(this, PETSC_ERR_ARG_NULL, .false.)
         end if
         call VecNorm(reference_vector, NORM_2, reference_norm, ierr)
         call CheckPEPSolveConfiguration(this, ierr, .false.)
         if (.not. ieee_is_finite(reference_norm) .or. reference_norm <= 0.0_REAL64) then
            print *, 'Candidate retry rejected: nonfinite or zero reference norm.'
            call CheckPEPSolveConfiguration(this, PETSC_ERR_ARG_OUTOFRANGE, .false.)
         end if
         call VecDuplicate(reference_vector, created_vector, ierr)
         call CheckPEPSolveConfiguration(this, ierr, .false.)
         call VecCopy(reference_vector, created_vector, ierr)
         copy_error = ierr
         if (copy_error /= 0) then
            call VecDestroy(created_vector, cleanup_error)
            call CheckSolverConfiguration(cleanup_error)
            call CheckPEPSolveConfiguration(this, copy_error, .false.)
         end if
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
      PetscErrorCode :: cleanup_error
      if (this%has_retry_reference) then
         call VecDestroy(this%retry_reference, cleanup_error)
         call CheckSolverConfiguration(cleanup_error)
      end if
      this%has_retry_reference = .false.
      this%retry_count = 0
   end subroutine ClearCandidateRetry

   !> Reuse the last selected eigenvector; fall back safely before any query.
   integer function SelectedGuessIndex(this) result(index)
      class(slepc_type), intent(in) :: this
      integer :: candidate
      index = 1
      if (.not. allocated(this%xr)) return
      if (size(this%xr) == 0) return
      candidate = this%overlap_selected_index
      if (candidate >= 1 .and. candidate <= size(this%xr)) then
         index = candidate
         return
      end if
      if (.not. allocated(this%evIndex)) return
      if (size(this%evIndex) /= size(this%xr)) return
      candidate = this%evIndex(size(this%evIndex))
      if (candidate >= 1 .and. candidate <= size(this%xr)) index = candidate
   end function SelectedGuessIndex

   !> Number of converged eigenpairs; safe before the first solve as well.
   integer function ConvergedCount(this) result(nconv)
      class(slepc_type), intent(in) :: this
      nconv = 0
      if (allocated(this%kr)) nconv = size(this%kr)
   end function ConvergedCount

   !> A failed scalar query must not index an empty eigenvalue array.
   function InvalidEigenvalue() result(value)
      PetscScalar :: value
      real(REAL64) :: nan_value
      nan_value = ieee_value(0.0_REAL64, ieee_quiet_nan)
      value = cmplx(nan_value, nan_value, kind=REAL64)
   end function InvalidEigenvalue

   !> Vec handles own PETSc objects and must be destroyed, not only deallocated.
   subroutine ReleaseEigenpairs(this)
      class(slepc_type), intent(inout) :: this
      integer :: i
      if (allocated(this%xr)) then
         do i = 1, size(this%xr)
            call VecDestroy(this%xr(i), ierr)
            call CheckSolverConfiguration(ierr)
         end do
         deallocate (this%xr)
      end if
      if (allocated(this%kr)) deallocate (this%kr)
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      this%last_error = huge(1.0_REAL64)
      this%overlap_selected_index = 0
   end subroutine ReleaseEigenpairs

   subroutine EPSFinalize(this)
      implicit none
      class(slepceps_type), intent(inout) :: this

      call ReleaseEigenpairs(this)
      call EPSDestroy(this%solver, ierr)
      call CheckSolverConfiguration(ierr)

   end subroutine EPSFinalize

   subroutine PEPFinalize(this)
      implicit none
      class(slepcpep_type), intent(inout) :: this

      call ClearCandidateRetry(this)
      this%restore_dimensions = .false.
      call ReleaseEigenpairs(this)
      call PEPDestroy(this%solver, ierr)
      call CheckSolverConfiguration(ierr)

   end subroutine PEPFinalize

end module mod_slepc_adapter
