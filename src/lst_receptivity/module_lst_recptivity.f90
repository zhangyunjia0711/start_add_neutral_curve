!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-11-01 11:44:33
!> @LastEditTime: 2022-12-12 19:32:44
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lst_receptivity/module_lst_recptivity.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
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

   Vec :: source_whole
   complex(REAL64), allocatable :: Cr(:)

#include "../include/start.f90"

contains

   subroutine solve_lst_receptivity()

      implicit none

      allocate (Cr(in))

      call load_source_term()
      call solve_lst_receptivity_local()
      call finalize_receptivity()

      deallocate (Cr)

   end subroutine solve_lst_receptivity

   subroutine load_source_term()

      implicit none
      PetscViewer :: viewer

      call VecDuplicate(Flowfield, source_whole, ierr)
      if (char(source_term_file) /= 'none') then
         call PetscViewerHDF5Open(PETSC_COMM_WORLD, char(source_term_file), FILE_MODE_READ, viewer, ierr)
         call PetscObjectSetName(source_whole, "source", ierr)
         call VecLoad(source_whole, viewer, ierr)
         call PetscViewerDestroy(viewer, ierr)
      else
         call VecZeroEntries(source_whole, ierr)
      end if

   end subroutine load_source_term

   subroutine finalize_receptivity()

      implicit none

      call DMRestoreGlobalVector(meshDA, source_whole, ierr)
      call DMRestoreGlobalVector(coordDA, Coord, ierr)
      call DMRestoreGlobalVector(meshDA, Flowfield, ierr)
      call DMRestoreLocalVector(coordDA, Coord_local, ierr)
      call DMRestoreLocalVector(meshDA, Flowfield_local, ierr)
      call DMDestroy(meshDA, ierr)
      call DMDestroy(coordDA, ierr)
      !   call DMDestroy(DA, ierr)

   end subroutine

   subroutine solve_lst_receptivity_local()

      implicit none
      complex(REAL64) :: ome, alf, bat
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile
      type(type_lst_eqn) :: lst
      type(type_alst_eqn) :: alst
      Vec :: lst_eigf, alst_eigf, vec_local, source_local

      integer :: i
      Mat :: mat_local
      complex(REAL64) :: Q, force_source, bc_force
      complex(REAL64) :: alf_iloc

      open (999, file='out/receptivity_coef.plt')

      call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call lst%set_parameter("Pr", Pr)
      lst%isothermal = isothermal

      call alst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call alst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call alst%set_parameter("Pr", Pr)
      alst%isothermal = isothermal

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

      alf_iloc = alf

      call MatDuplicate(lst%L_LST(2), MAT_COPY_VALUES, mat_local, ierr)

      call lst%get_eigf(lst_eigf)
      call alst%get_eigf(alst_eigf)

      call VecDuplicate(lst_eigf, vec_local, ierr)
      call VecDuplicate(lst_eigf, source_local, ierr)

      !! get_source_local_value
      call get_source_value_local(source_whole, source_local, iloc)

        !! compute receptivity coefficient
        !! $C_r = i frac{<\hat phi_{\alpha_s}, \hat S> - B.C.}{Q}

      !Q
      call MatAXPY(mat_local, 2.0d0*alf, lst%L_LST(3), SUBSET_NONZERO_PATTERN, ierr)
      call MatMult(mat_local, lst_eigf, vec_local, ierr)
      q = inner_dot(alst_eigf, vec_local, grid_profile%y)

      ! <\hat phi_{\alpha_s}, \hat S>
      force_source = inner_dot(alst_eigf, source_local, grid_profile%y)

      ! b.c.
      bc_force = integral_bc(baseflow_profile, grid_profile, alst_eigf, cmplx([0.0d0, 1.0d0, 0.0d0, 0.0d0]))

      Cr(iloc) = (force_source - bc_force)/Q*ZI

      do i = iloc, iend
         grid_profile = assemble_lst_grid_profile(iloc)
         baseflow_profile = assemble_lst_baseflow_profile(iloc)

         call lst%set_grid(grid_profile)
         call lst%set_baseflow(baseflow_profile)
         call lst%setup()

         call alst%set_grid(grid_profile)
         call alst%set_baseflow(baseflow_profile)
         call alst%setup()

         call lst%solve(ome, alf, bat)
         call alst%solve(ome, alf, bat)

         call lst%get_eigf(lst_eigf)
         call alst%get_eigf(alst_eigf)

         call get_source_value_local(source_whole, source_local, i)

        !! compute receptivity coefficient
        !! $C_r = i frac{<\hat phi_{\alpha_s}, \hat S> - B.C.}{Q}

         !Q
         call MatCopy(lst%L_LST(2), mat_local, DIFFERENT_NONZERO_PATTERN, ierr)
         call MatAXPY(mat_local, 2.0d0*alf, lst%L_LST(3), SUBSET_NONZERO_PATTERN, ierr)
         call MatMult(mat_local, lst_eigf, vec_local, ierr)
         q = inner_dot(alst_eigf, vec_local, grid_profile%y)

         ! <\hat phi_{\alpha_s}, \hat S>
         force_source = inner_dot(alst_eigf, source_local, grid_profile%y)

         ! b.c.
         bc_force = integral_bc(baseflow_profile, grid_profile, alst_eigf, cmplx(bc_disturb))

         Cr(i) = (force_source - bc_force)/Q*ZI
         write (999, *) i, Cr(i), abs(Cr(i))
      end do

      alf = alf_iloc

      do i = iloc, istart, -1
         grid_profile = assemble_lst_grid_profile(iloc)
         baseflow_profile = assemble_lst_baseflow_profile(iloc)

         call lst%set_grid(grid_profile)
         call lst%set_baseflow(baseflow_profile)
         call lst%setup()

         call alst%set_grid(grid_profile)
         call alst%set_baseflow(baseflow_profile)
         call alst%setup()

         call lst%solve(ome, alf, bat)
         call alst%solve(ome, alf, bat)

         call lst%get_eigf(lst_eigf)
         call alst%get_eigf(alst_eigf)

         call get_source_value_local(source_whole, source_local, i)

       !! compute receptivity coefficient
       !! $C_r = i frac{<\hat phi_{\alpha_s}, \hat S> - B.C.}{Q}

         !Q
         call MatCopy(lst%L_LST(2), mat_local, DIFFERENT_NONZERO_PATTERN, ierr)
         call MatAXPY(mat_local, 2.0d0*alf, lst%L_LST(3), SUBSET_NONZERO_PATTERN, ierr)
         call MatMult(mat_local, lst_eigf, vec_local, ierr)
         !call VecDot(alst_eigf, vec_local, Q, ierr)
         q = inner_dot(alst_eigf, vec_local, grid_profile%y)

         ! <\hat phi_{\alpha_s}, \hat S>
         force_source = inner_dot(alst_eigf, source_local, grid_profile%y)

         ! b.c.
         bc_force = integral_bc(baseflow_profile, grid_profile, alst_eigf, cmplx(bc_disturb))

         Cr(i) = (force_source - bc_force)/Q*ZI
         if(i<iloc) write (999, *) i, Cr(i), abs(Cr(i))
      end do

      call lst%finalize()
      call alst%finalize()
      call grid_profile%finalize()
      call baseflow_profile%finalize()
      call MatDestroy(mat_local, ierr)
      call VecDestroy(lst_eigf, ierr)
      call VecDestroy(alst_eigf, ierr)
      call VecDestroy(vec_local, ierr)
      call VecDestroy(source_local, ierr)

      close (999)

   end subroutine solve_lst_receptivity_local

   subroutine get_source_value_local(whole, local, iloc)

      Vec, intent(in) :: whole
      Vec, intent(inout) :: local
      integer, intent(in) :: iloc

      complex(REAL64), pointer :: array_whole(:, :, :), array_local(:)

      ! integer :: size_local, size_whole, in

      ! call VecGetSize(whole, size_whole, ierr)
      ! call VecGetSize(local, size_local, ierr)

      ! in = size_whole/size_local

      call DMDAVecGetArrayReadF90(meshDA, whole, array_whole, ierr)
      call VecGetArrayF90(local, array_local, ierr)

      array_local = pack(array_whole(:, iloc - 1, :), mask=.true.)

      call VecRestoreArrayF90(local, array_local, ierr)
      call DMDAVecRestoreArrayReadF90(meshDA, whole, array_whole, ierr)

   end subroutine get_source_value_local

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

      bc_force = Rho0_w*conjg(rho_star_w)*v_w + &
                 mu0*conjg(du_star_w)*u_w + &
                 4.0d0/3.0d0*mu0*conjg(dv_star_w)*v_w + &
                 mu0*conjg(dw_star_w)*w_w + &
                 mu0/Pr*conjg(dT_star_w)*T_w

   end function integral_bc

   function inner_dot(v1, v2, y) result(q)

      use stdlib_quadrature, only: trapz_weights
      implicit none
      
      Vec, intent(in) :: v1, v2
      real(REAL64), intent(in) :: y(:)
      complex(REAL64) :: q

      Vec :: tmp
      complex(REAL64), pointer:: array(:)
      integer :: j
      complex(REAL64) :: tmp_(0:jn-1)

      call VecDuplicate(v1, tmp, ierr)
      call VecCopy(v1, tmp, ierr)
      call VecConjugate(tmp, ierr)
      call VecPointwiseMult(tmp, tmp, v2, ierr)
      call VecGetArrayReadF90(tmp, array, ierr)

      do j=0, jn-1
         tmp_(j)=sum(array(j*5+1:j*5+5))
      enddo
      q = sum(trapz_weights(y)*tmp_)
      call VecRestoreArrayReadF90(tmp, array, ierr)

   end function inner_dot

end module mod_lst_receptivity
