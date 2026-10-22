!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-05 09:32:30
!> @LastEditTime: 2022-12-12 19:33:30
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lst/lst_solver_abstract.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lst_solver_abstract

   use petsc
   use iso_fortran_env, only: REAL64, REAL32
   use mod_lst_grid_profile, only: type_lst_grid_profile
   use mod_lst_baseflow_profile, only: type_lst_bf_profile

   use mod_solver_abstract

   implicit none
   private
   public :: type_lst_abstract
   public :: normal_u

#include "../include/start.f90"

   type, extends(type_solver_abstract), abstract :: type_lst_abstract
      integer :: lst_method = LST_LOCAL
      integer :: lst_mode = LST_TEMPORAL_MODE
      integer :: lst_nev = 1

   contains

      procedure(initial_interface), deferred :: initial
      procedure(setup_interface), deferred :: setup
      procedure(solve_interface), deferred :: solve
      procedure(get_eigf_interface), deferred :: get_eigf

      procedure :: set_grid => lst_set_grid
      procedure :: set_baseflow => lst_set_baseflow

      procedure :: set_parameters_real8 => lst_set_parameters_real8
      procedure :: set_parameters_real4 => lst_set_parameters_real4
      procedure :: set_parameters_int => lst_set_parameters_int

      procedure(finalize_interface), deferred :: finalize

   end type type_lst_abstract

   abstract interface

      subroutine initial_interface(this, jn_, Re, Ma, Te, yb_, ye_)

         import type_lst_abstract, real64
         implicit none
         class(type_lst_abstract), intent(inout) :: this
         integer, intent(in) :: jn_
         real(real64), intent(in) ::  Re, Ma, Te
         integer, intent(in), optional :: yb_, ye_
      end subroutine initial_interface

      subroutine setup_interface(this)
         import type_lst_abstract
         implicit none
         class(type_lst_abstract), intent(inout) :: this

      end subroutine setup_interface

      subroutine solve_interface(this, ome, alf, bat)

         import :: type_lst_abstract, real64
         implicit none
         class(type_lst_abstract), intent(in) :: this
         complex(REAL64), intent(inout) :: ome, alf, bat

      end subroutine solve_interface

      subroutine get_eigf_interface(this, eigf)

         import :: type_lst_abstract
         import :: tvec
         implicit none
         class(type_lst_abstract), intent(in) :: this
         Vec, intent(out) :: eigf

      end subroutine get_eigf_interface

      subroutine finalize_interface(this)

         import type_lst_abstract
         implicit none
         class(type_lst_abstract), intent(inout) :: this

      end subroutine finalize_interface

   end interface
contains

   subroutine lst_set_parameters_real8(this, keyword, value)

      implicit none
      class(type_lst_abstract), intent(inout) :: this
      character(*), intent(in) :: keyword
      real(REAL64), intent(in) ::  value

      select case (trim(keyword))
      case ("Ma", "MA", "ma")
         this%Ma = value
      case ("Re", "RE", "re")
         this%Re = value
      case ("Te", "TE", "te")
         this%Te = value
      case ("GAMMA", "Gamma", "gamma")
         this%GAMMA = value
      case ("Pr", "pr", "PR")
         this%Pr = value
      case default
         stop "the value is wrong in set_parameters in lst."
      end select

   end subroutine lst_set_parameters_real8

   subroutine lst_set_parameters_real4(this, keyword, value)

      implicit none
      class(type_lst_abstract), intent(inout) :: this
      character(*), intent(in) :: keyword
      real(REAL32), intent(in) ::  value

      select case (trim(keyword))
      case ("Ma", "MA", "ma")
         this%Ma = value
      case ("Re", "RE", "re")
         this%Re = value
      case ("Te", "TE", "te")
         this%Te = value
      case ("GAMMA", "Gamma", "gamma")
         this%GAMMA = value
      case ("Pr", "pr", "PR")
         this%Pr = value
      case default
         stop "the value is wrong in set_parameters in lst."
      end select

   end subroutine lst_set_parameters_real4

   subroutine lst_set_parameters_int(this, keyword, value)

      implicit none
      class(type_lst_abstract), intent(inout) :: this
      character(*), intent(in) :: keyword
      integer, intent(in) ::  value

      select case (trim(keyword))
      case ("bc_y0")
         this%bc_y0 = value
      case ("bc_yn")
         this%bc_yn = value
      case ("lst_mode")
         this%lst_mode = value
      case ("Coordinate_System")
         this%Coordinate_System = value
      case default
         stop "the value is wrong in set_parameters in lst."
      end select

   end subroutine lst_set_parameters_int

   subroutine lst_set_grid(this, grid)
      implicit none
      class(type_lst_abstract), intent(inout) :: this
      type(type_lst_grid_profile), intent(in) :: grid

      this%grid = grid

   end subroutine lst_set_grid

   subroutine lst_set_baseflow(this, baseflow)
      implicit none
      class(type_lst_abstract), intent(inout) :: this
      type(type_lst_bf_profile), intent(in) :: baseflow

      this%baseflow = baseflow

   end subroutine lst_set_baseflow

   subroutine normal_u(v)

      implicit none
      Vec :: v
      Vec :: v_all
      PetscScalar, pointer :: tmp_(:)
      PetscScalar :: maxu
      integer :: maxu_loc(1)
      integer :: jn
      PetscScalar, allocatable :: phi(:, :)
      VecScatter :: ctx
      PetscErrorCode :: ierr

      call VecScatterCreateToZero(v, ctx, v_all, ierr)
      call CheckNormalizationError(ierr)

      call VecScatterBegin(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call CheckNormalizationError(ierr)
      call VecScatterEnd(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call CheckNormalizationError(ierr)

      call VecGetArrayReadF90(v_all, tmp_, ierr)
      call CheckNormalizationError(ierr)
      jn = size(tmp_)/5
      allocate (phi(5, jn))
      phi(1:5, 1:jn) = reshape(tmp_, [5, jn])
      call VecRestoreArrayReadF90(v_all, tmp_, ierr)
      call CheckNormalizationError(ierr)

      maxu_loc = maxloc(abs(phi(2, :)))
      maxu = phi(2, maxu_loc(1))
      deallocate (phi)

      ! call MPI_bcast(maxu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call VecScale(v, 1.0d0/maxu, ierr)
      call CheckNormalizationError(ierr)

      call VecScatterDestroy(ctx, ierr)
      call CheckNormalizationError(ierr)
      call VecDestroy(v_all, ierr)
      call CheckNormalizationError(ierr)

   end subroutine normal_u


   !> PETSc failures must not be hidden by subsequent temporary-object calls.
   subroutine CheckNormalizationError(error)
      PetscErrorCode, intent(in) :: error
      PetscErrorCode :: report_error
      if (error /= 0) then
         report_error = error
         call PetscErrorF(report_error, __LINE__, "lst_solver_abstract.f90")
         call MPIU_Abort(PETSC_COMM_WORLD, error)
      end if
   end subroutine CheckNormalizationError

end module
