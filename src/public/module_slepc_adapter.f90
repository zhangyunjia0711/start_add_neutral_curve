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
      Vec :: initialguess !初始猜测向量
      integer, allocatable :: evIndex(:) !<特征值序号
      PetscScalar :: target !< 目标值
      ST :: st
      PC :: pc
      KSP :: ksp
      integer :: lst_mode

   contains

      

      procedure(create_interface), deferred :: Create
      procedure(solve_interface), deferred :: Solve
      procedure(GetEigenVal_interface), deferred :: GetEigenVal
      procedure(GetEigenVals_interface), deferred :: GetEigenVals
      procedure(GetEigenVector_interface), deferred :: GetEigenVector
      procedure(GetEigenVectors_interface), deferred :: GetEigenVectors
      procedure(Finalize_interface), deferred :: Finalize

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
      ST :: st

      call PEPCreate(PETSC_COMM_WORLD, this%solver, ierr)
      this%lst_mode = lst_mode
      call PEPSetWhichEigenpairs(this%solver, PEP_TARGET_MAGNITUDE, ierr)
      call PEPSetType(this%solver, PEPLINEAR, ierr)
      call PEPSetConvergenceTest(this%solver, PEP_CONV_NORM, ierr)
      call PEPGetST(this%solver, st, ierr)
      call STSetType(st, STSINVERT, ierr)
      call PEPSetTolerances(this%solver, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call PEPSetDimensions(this%solver, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
      ! print*, nev, 'eigenvalues will be computed.'

   end subroutine PEPSolverCreate

   !> 创建求解器
   subroutine EPSSolverCreate(this, nev, lst_mode)

      implicit none
      class(slepceps_type), intent(inout) :: this
      integer, intent(in) :: nev !<求解特征值个数
      integer, intent(in) :: lst_mode
      ST :: st

      call EPSCreate(PETSC_COMM_WORLD, this%solver, ierr)
      this%lst_mode = lst_mode
      call EPSGetST(this%solver, st, ierr)
      call EPSSetWhichEigenpairs(this%solver, EPS_TARGET_MAGNITUDE, ierr)
      call STSetType(st, STSINVERT, ierr)
      call EPSSetTolerances(this%solver, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call EPSSetDimensions(this%solver, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
      ! print*, nev, 'eigenvalues will be computed.'

   end subroutine EPSSolverCreate

   !> 求解特征值
   subroutine PEPSolverSolve(this, nmat, mat, target_input)

      class(slepcpep_type), intent(inout) :: this
      integer, intent(in) :: nmat
      Mat, intent(in) :: mat(nmat)
      complex(REAL64), intent(in) :: target_input
      PetscScalar :: target_value
      PetscInt :: nconv, nrow
      integer :: i

      call PEPSetProblemType(this%solver, PEP_GENERAL, ierr)
      ! call PEPSetType(this%solver, PEPLINEAR, ierr)
      call PEPSetOperators(this%solver, nmat, mat, ierr)
      call MatGetSize(mat(1), this%m, nrow, ierr)
      target_value = target_input
      this%target = target_input
      ! print*, 'The pep solver is selected!'
      ! print*, 'The target_value is', target_value
      call PEPSetTarget(this%solver, target_value, ierr)
      call PEPSetType(this%solver, PEPLINEAR, ierr)
      if (allocated(this%xr)) then
         ! print*, 'Setting initial guess space....'
         call VecDuplicate(this%xr(size(this%xr)), this%initialguess, ierr)
         call VecCopy(this%xr(size(this%xr)), this%initialguess, ierr)
         call PEPSetInitialSpace(this%solver, 1, this%initialguess, ierr)
      end if
      !call PEPSetRefine(this%solver, PEP_REFINE_SIMPLE, 1, 1.0d-8, 2, PEP_REFINE_SCHEME_SCHUR, ierr)
      call PEPGetST(this%solver, this%st, ierr)
      call STSetType(this%st, STSINVERT, ierr)
      call STGetKSP(this%st, this%ksp, ierr)
      call KSPGetPC(this%ksp, this%pc, ierr); 
      call KSPSetType(this%ksp, KSPPREONLY, ierr)
      call PCSetType(this%pc, PCLU, ierr)
#ifdef __linux__
      call PCFactorSetMatSolverType(this%pc, MATSOLVERMKL_CPARDISO, ierr)
#elif __APPLE__
      call PCFactorSetMatSolverType(this%pc, MATSOLVERSUPERLU_DIST, ierr)
#endif
      call PEPSetFromOptions(this%solver, ierr)
      call PEPSolve(this%Solver, ierr)
      call PEPGetConverged(this%Solver, nconv, ierr)
      ! write(*, *)'ncov=:', nconv
      if (allocated(this%xr)) then
         if (size(this%xr) /= nconv) then
            deallocate (this%xr)
            allocate (this%xr(nconv))
         end if
      else
         allocate (this%xr(nconv))
      end if
      if (allocated(this%kr)) then
         if (size(this%kr) /= nconv) then
            deallocate (this%kr)
            allocate (this%kr(nconv))
         end if
      else
         allocate (this%kr(nconv))
      end if
      do i = 0, nconv - 1
         call MatCreateVecs(mat(1), this%xr(i + 1), PETSC_NULL_VEC, ierr)
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
      PetscInt :: nconv, nrow
      integer :: i

      call EPSSetProblemType(this%solver, EPS_GNHEP, ierr)
      call EPSSetType(this%solver, EPSARNOLDI, ierr)
      call EPSSetOperators(this%solver, mat(1), mat(2), ierr)
      ! call MatGetSize(mat(1), this%m, nrow, ierr)
      target_value = target_input
      this%target = target_input
      ! print*, 'The eps solver is selected!'
      ! print*, 'The target_value is', target_value
      call EPSSetTarget(this%solver, target_value, ierr)
      call EPSGetST(this%solver, this%st, ierr)
      call STSetType(this%st, STSINVERT, ierr)
      call STGetKSP(this%st, this%ksp, ierr)
      call KSPGetPC(this%ksp, this%pc, ierr); 
      call KSPSetType(this%ksp, KSPPREONLY, ierr)
      call PCSetType(this%pc, PCLU, ierr)
#ifdef __linux__
      call PCFactorSetMatSolverType(this%pc, MATSOLVERMKL_CPARDISO, ierr)
#elif __APPLE__
      call PCFactorSetMatSolverType(this%pc, MATSOLVERSUPERLU_DIST, ierr)
#endif
      call EPSSetFromOptions(this%solver, ierr)
      if (allocated(this%xr)) then
         ! print*, 'Setting initial guess space....'
         call VecDuplicate(this%xr(size(this%xr)), this%initialguess, ierr)
         call VecCopy(this%xr(size(this%xr)), this%initialguess, ierr)
         call EPSSetInitialSpace(this%solver, 1, this%initialguess, ierr)
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

      call EPSGetConverged(this%Solver, nconv, ierr)
      ! write(*, *)'ncov=:', nconv
      !print*, size(this%xr)
      !print*, size(this%kr)
      if (allocated(this%xr)) then
         if (size(this%xr) /= nconv) then
            deallocate (this%xr)
            allocate (this%xr(nconv))
         end if
      else
         allocate (this%xr(nconv))
      end if
      if (allocated(this%kr)) then
         if (size(this%kr) /= nconv) then
            deallocate (this%kr)
            allocate (this%kr(nconv))
         end if
      else
         allocate (this%kr(nconv))
      end if
      do i = 0, nconv - 1
         call MatCreateVecs(mat(1), this%xr(i + 1), PETSC_NULL_VEC, ierr)
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

      nconv = size(this%kr)

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

      eigenvalue = this%kr(this%evIndex(nconv))
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

      nconv = size(this%kr)

      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

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

      nconv = size(this%kr)

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

      eigenvalue = this%kr(this%evIndex(nconv))

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

      nconv = size(this%kr)

      allocate (eigenvalue(nconv))
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

      nconv = size(this%kr)
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

      eigenvalue = this%kr(this%evIndex(nconv))
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

      nconv = size(this%kr)
      allocate (eigenVector(nconv))
      allocate (eigenvalue(nconv))
      if (allocated(this%evIndex)) deallocate (this%evIndex)
      allocate (this%evIndex(nconv))

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

      nconv = size(this%kr)
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

      eigenvalue = this%kr(this%evIndex(nconv))
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

      nconv = size(this%kr)
      allocate (eigenVector(nconv))
      allocate (eigenvalue(nconv))
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

   subroutine EPSFinalize(this)
      implicit none
      class(slepceps_type), intent(inout) :: this

      call EPSDestroy(this%solver, ierr)

   end subroutine EPSFinalize

   subroutine PEPFinalize(this)
      implicit none
      class(slepcpep_type), intent(inout) :: this

      call PEPDestroy(this%solver, ierr)

   end subroutine PEPFinalize

end module mod_slepc_adapter
