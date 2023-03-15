!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2023-03-15 09:49:31
!> @LastEditTime: 2023-03-15 09:50:24
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lpse/lpse_solver_abstract.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lpse_solver_abstract

   use petsc
   use iso_fortran_env, only: REAL64, REAL32

   use mod_solver_abstract

   implicit none
   private
   public :: type_lpse_abstract

#include "../include/start.f90"

   type, extends(type_solver_abstract), abstract :: type_lpse_abstract

      integer :: jb, je
      integer :: jn
      complex(real64) :: ome
      complex(real64) :: bat
      complex(real64) :: cycle_m = -1
      complex(real64), allocatable, dimension(:, :) :: phi0_array, phi1_array, phi2_array
      complex(REAL64), allocatable, dimension(:, :, :) :: B_0, D_0, V_0
      complex(REAL64), allocatable, dimension(:, :, :) :: D_2, G_0
      Vec :: phi
      DM :: dm
      Mat :: L_LPSE
      Vec :: RHS
      KSP :: solver
      PC :: SolverPC

      
      real(REAL64) :: lambda !< control d/dx order
      integer :: solve_alpha = SOLVE_ALPHA_BOTH
      real(REAL64) :: sigma_control = 0.0d0
      real(REAL64) :: sigma_p = 1.0d0
      real(REAL64) :: rd_factor = 1.3d0

   contains

      procedure(initial_interface), deferred :: initial
      procedure(setup_interface), deferred :: setup
      procedure(set_ic_interface), deferred :: set_ic
      procedure(solve_interface), deferred :: solve
      procedure(finalize_interface), deferred :: finalize

      !   procedure :: set_grid=> lst_set_grid
      !   procedure :: set_baseflow => lst_set_baseflow

      procedure :: set_parameters_real8 => lpse_set_parameters_real8
      procedure :: set_parameters_real4 => lpse_set_parameters_real4
      procedure :: set_parameters_int => lpse_set_parameters_int

   end type type_lpse_abstract

   abstract interface

      subroutine initial_interface(this, jn_, Re, Ma, Te, yb_, ye_)

         import type_lpse_abstract, real64
         implicit none
         class(type_lpse_abstract), intent(inout) :: this
         integer, intent(in) :: jn_
         real(real64), intent(in) ::  Re, Ma, Te
         integer, intent(in), optional :: yb_, ye_
      end subroutine initial_interface

      subroutine setup_interface(this)
         import type_lpse_abstract
         implicit none
         class(type_lpse_abstract), intent(inout) :: this

      end subroutine setup_interface

      subroutine solve_interface(this, alf, phi0, phi1, phi2)

         import :: type_lpse_abstract, real64, tvec
         implicit none
         class(type_lpse_abstract), intent(inout) :: this
         complex(REAL64), intent(inout) :: alf
         Vec, intent(inout) :: phi0
         Vec, intent(in), optional :: phi1, phi2

      end subroutine solve_interface

      subroutine set_ic_interface(this, phi)

         import type_lpse_abstract, tvec
         implicit none
         class(type_lpse_abstract), intent(inout) :: this
         Vec, intent(inout) :: phi

      end subroutine set_ic_interface

      subroutine finalize_interface(this)

         import type_lpse_abstract
         implicit none
         class(type_lpse_abstract), intent(inout) :: this

      end subroutine finalize_interface

   end interface
contains

   subroutine lpse_set_parameters_real8(this, keyword, value)

      implicit none
      class(type_lpse_abstract), intent(inout) :: this
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
      case ("ome")
         this%ome = value
      case ("bat")
         this%bat = value
         this%cycle_m = -1
      case ("cycle_m")
         this%cycle_m = value
      case ("sigma_control")
         this%sigma_control = value
      case ("sigma_p")
         this%sigma_p = value
      case ("rd_factor")
         this%rd_factor = value
      case default
         stop "the value is wrong in set_parameters in lpse."
      end select

   end subroutine lpse_set_parameters_real8

   subroutine lpse_set_parameters_real4(this, keyword, value)

      implicit none
      class(type_lpse_abstract), intent(inout) :: this
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
      case ("ome")
         this%ome = value
      case ("bat")
         this%bat = value
         this%cycle_m = -1
      case ("cycle_m")
         this%cycle_m = value
      case ("sigma_control")
         this%sigma_control = value
      case ("sigma_p")
         this%sigma_p = value
      case ("rd_factor")
         this%rd_factor = value
      case default
         stop "the value is wrong in set_parameters in lpse."
      end select

   end subroutine lpse_set_parameters_real4

   subroutine lpse_set_parameters_int(this, keyword, value)

      implicit none
      class(type_lpse_abstract), intent(inout) :: this
      character(*), intent(in) :: keyword
      integer, intent(in) ::  value

      select case (trim(keyword))
      case ("bc_y0")
         this%bc_y0 = value
      case ("bc_yn")
         this%bc_yn = value
      case ("Coordinate_System")
         this%Coordinate_System = value
      case ("solve_alpha")
         this%solve_alpha = value
      case default
         stop "the value is wrong in set_parameters in lpse."
      end select

   end subroutine lpse_set_parameters_int

   ! subroutine lpse_set_grid(this, grid)
   !         implicit none
   !         class(type_lpse_abstract),intent(inout) :: this
   !         type(type_lst_grid_profile),intent(in) :: grid

   !           this%grid=grid

   !   end subroutine lpse_set_grid

   !   subroutine lpse_set_baseflow(this, baseflow)
   !   implicit none
   !   class(type_lpse_abstract),intent(inout) :: this
   !   type(type_lst_bf_profile),intent(in) :: baseflow

   !         this%baseflow=baseflow

   !   end subroutine lpse_set_baseflow
end module
