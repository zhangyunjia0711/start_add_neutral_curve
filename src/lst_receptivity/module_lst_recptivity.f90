!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  module: mod_solve_lst
!> @breif module for solving lst receptivity.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2021-01-05 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2021-01-05
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lst_receptivity

   use global_parameters
   use mod_lst_baseflow_profile
   use mod_lst_grid_profile
   use mod_lst_solver_abstract
   use mod_lst_eqn
   use mod_alst_eqn
   use stdlib_string_type, only: char

   implicit none

   private
   public :: solve_lst_receptivity

#include "../include/start.f90"

contains

   subroutine solve_lst_receptivity()

      implicit none

      call solve_lst_receptivity_local()

   end subroutine solve_lst_receptivity

   subroutine solve_lst_receptivity_local()

      implicit none
      complex(REAL64) :: ome, alf, bat
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile
      type(type_lst_eqn) :: lst
      type(type_alst_eqn) :: alst
      Vec :: lst_eigf, alst_eigf, vec_local, source_local
      complex(REAL64) :: Cr

      integer :: k, ki
      real(REAL64) :: pk
      Mat :: mat_local
      complex(REAL64) :: Q, force_source, bc_force

      call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call lst%set_parameter("Pr", Pr)

      call alst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call alst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call alst%set_parameter("Pr", Pr)

      grid_profile = assemble_lst_grid_profile(iloc)
      baseflow_profile = assemble_lst_baseflow_profile(iloc)
      ome = omega(1) + ZI*omega(2)
      alf = alpha(1) + ZI*alpha(2)
      bat = beta(1) + ZI*beta(2)

      call lst%set_grid(grid_profile)
      call lst%set_baseflow(baseflow_profile)
      call lst%setup()

      call alst%set_grid(grid_profile)
      call alst%set_baseflow(baseflow_profile)
      call alst%setup()

      call lst%solve(ome, alf, bat)
      call alst%solve(ome, alf, bat)

      call MatDuplicate(lst%L_LST(2), MAT_COPY_VALUES, mat_local, ierr)

      call lst%get_eigf(lst_eigf)
      call alst%get_eigf(alst_eigf)

      call VecDuplicate(lst_eigf, vec_local, ierr)

        !! compute receptivity coefficient
        !! $C_r = i frac{<\hat phi_{\alpha_s}, \hat S> - B.C.}{Q}

      !Q
      call MatAXPY(mat_local, 2.0d0*alf, lst%L_LST(3), SUBSET_NONZERO_PATTERN, ierr)
      call MatMult(mat_local, lst_eigf, vec_local, ierr)
      call VecDot(alst_eigf, vec_local, Q, ierr)

      ! <\hat phi_{\alpha_s}, \hat S>
      call VecDot(alst_eigf, source_local, force_source, ierr)

      ! b.c.
      bc_force = integral_bc(baseflow_profile, grid_profile, alst_eigf, cmplx([0.0d0, 1.0d0, 0.0d0, 0.0d0]))

      Cr = (force_source - bc_force)/Q*ZI

      call lst%finalize()
      call alst%finalize()
      call grid_profile%finalize()
      call baseflow_profile%finalize()
      call MatDestroy(mat_local, ierr)
      call VecDestroy(lst_eigf, ierr)
      call VecDestroy(alst_eigf, ierr)
      call VecDestroy(vec_local, ierr)
      call VecDestroy(source_local, ierr)

   end subroutine solve_lst_receptivity_local

   function integral_bc(baseflow, grid, alst_eigf, bcwall) result(bc_force)
      implicit none
      type(type_lst_bf_profile), intent(in) :: baseflow
      type(type_lst_grid_profile), intent(in) :: grid
      Vec, intent(in) :: alst_eigf
      complex, intent(in) :: bcwall(4)
      complex(REAL64) :: bc_force

      complex(REAL64) :: u_w, v_w, w_w, T_w
      complex(REAL64) :: rho_star_w, du_star_w, dv_star_w, dw_star_w, dT_star_w
      real(REAL64) :: Rho0_w, mu0
      complex(REAL64), pointer :: tmp_(:)

      u_w = bcwall(1); v_w = bcwall(2); w_w = bcwall(3); T_w = bcwall(4)
      Rho0_w = baseflow%q0(1, 1)

      !! set mu0
      associate (T => baseflow%q0(5, 1), C1 => 110.4D0)
         mu0 = T*sqrt(T)*(1.0d0 + C1/Te)/(T + C1/Te)/Re
      end associate
      !!set rho_star_w, du_star_w, dv_star_w, dw_star_w, dT_star_w value
      call VecGetArrayReadF90(alst_eigf, tmp_, ierr)
      rho_star_w = tmp_(1)
      du_star_w = (-1.0d0*tmp_(12) + 4.0d0*tmp_(7) - 3.0d0*tmp_(2))/2.0d0*grid%eta_y(1)
      dv_star_w = (-1.0d0*tmp_(13) + 4.0d0*tmp_(8) - 3.0d0*tmp_(3))/2.0d0*grid%eta_y(1)
      dw_star_w = (-1.0d0*tmp_(14) + 4.0d0*tmp_(9) - 3.0d0*tmp_(4))/2.0d0*grid%eta_y(1)
      dT_star_w = (-1.0d0*tmp_(15) + 4.0d0*tmp_(10) - 3.0d0*tmp_(5))/2.0d0*grid%eta_y(1)
      call VecRestoreArrayReadF90(alst_eigf, tmp_, ierr)

      bc_force = Rho0_w*rho_star_w*v_w + &
                 mu0*du_star_w*u_w + &
                 3.0d0/4.0d0*mu0*dv_star_w*v_w + &
                 mu0*dw_star_w*w_w + &
                 mu0/Pr*dT_star_w*T_w

   end function integral_bc

end module mod_lst_receptivity
