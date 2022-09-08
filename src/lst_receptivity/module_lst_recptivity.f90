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

   subroutine write_lst_function(ome, alf, bat, ix, v)

      use stdlib_string_type
      use global_parameters
      implicit none
      complex(REAL64) :: alf, bat, ome
      integer :: ix
      Vec :: v

      type(string_type) :: fn
      character(255) :: tmp
      PetscScalar :: d
      Vec :: v_all
      VecScatter :: ctx
      PetscScalar, pointer :: tmp_(:)
      Complex(REAL64) :: phi(6, jn)
      character(40) :: name
      integer :: l, j, maxu_loc(1)
      real(REAL64) :: pe
      complex(REAL64):: maxu
      PetscViewer :: Viewer
      Vec :: tmpV

      ! select case (lst_mode)
      ! case (LST_TEMPORAL_MODE)
      !   ome=ome+(0.0d0, 1.0d-8)
      ! case (LST_SPATIAL_MODE)
      !   alf=alf+(0.0d0, 1.0d-8)
      ! end select

      call PetscPrintf(PETSC_COMM_WORLD, "Writing "//char(solname)//" Eigenfunction......\n", ierr)

      write (tmp, ('(I0)')) ix
      fn = solname//'_eigenfunction_i='//trim(adjustl(tmp))//'.plt'
      ! write(tmp, ('(F20.5)'))ome%re
      ! fn=fn//'_ome='//trim(adjustl(tmp))
      ! write(tmp, ('(F20.5)'))alf%re
      ! fn=fn//'_alf='//trim(adjustl(tmp))
      ! write(tmp, ('(F20.5)'))bat%re
      ! fn=fn//'_bat='//trim(adjustl(tmp))//'.plt'

      call VecScatterCreateToZero(v, ctx, v_all, ierr)
      ! call VecGetType(v, name, ierr)
      ! call VecGetType(v_all, name, ierr)
      ! call VecGetSize(v, l, ierr)
      ! call VecGetSize(v_all, l, ierr)

      call VecScatterBegin(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)

      if (rank == 0) then
         call VecGetArrayReadF90(v_all, tmp_, ierr)
         phi(1:5, 1:jn) = reshape(tmp_, [5, jn])
         call VecRestoreArrayReadF90(v_all, tmp_, ierr)

         maxu_loc = maxloc(abs(phi(2, :)))
         maxu = phi(2, maxu_loc(1))
         phi = phi/maxu
         pe = 1.0d0/(GAMMA*MA*MA)
         do j = 1, jn
            phi(6, j) = (qq(1, ix, j)*phi(5, j) + qq(5, ix, j)*phi(1, j))*pe
         end do

         ! block
         !   use mod_dx
         !   integer :: l1, l2
         !   complex(REAL64), dimension(5, jn) :: dphidj, dphidjj
         !   do l=1, 5
         !     call dx_fdm4nd_center(phi(l, :), ybg+1, yeg+1, dphidj(l, :), ybg+1, yeg+1)
         !     call dxx_fdm4nd_center(phi(l, :), ybg+1, yeg+1, dphidjj(l, :), ybg+1, yeg+1)
         !   enddo
         !
         !
         !   alf1=d
         !   print*, ome1, alf1
         !   ! print*, grid_profile%eta_y(149), grid_profile%eta_yy(149)
         !   ! print*, lbound(grid_profile%eta_y), maxu_loc
         !   print*, phi(:, 149)
         !   print*, dphidj(:, 149)
         !   print*, dphidjj(:, 149)
         !   do l1=1, 5
         !     do l2=1, 5
         !       print*, l1, l2, A0(l1, l2, 149)
         !     enddo
         !   enddo
         !   print*, matmul(-ZI*ome1*G0(:, :, 149)+ZI*alf1*A0(:, :, 149)+D0(:, :, 149)+ &
         !   alf1**2*V0xx(:, :, 149), phi(1:5, 149))+ &
         !   matmul(B0(:, :, 149)-ZI*alf1*V0xy(:, :, 149), dphidj(1:5, 149)*grid_profile%eta_y(149))+ &
         !   matmul(-V0yy(:,:,149), dphidj(:, 149)*grid_profile%eta_yy(149)+dphidjj(1:5, 149)*grid_profile%eta_y(149)**2)
         !
         !   print*, matmul(-ZI*ome1*G0(:, :, 149)+ZI*alf1*A0(:, :, 149)+D0(:, :, 149)+ &
         !   alf1**2*V0xx(:, :, 149), phi(1:5, 149))
         !   read(*, *)
         ! endblock

         call execute_command_line('mkdir -p out')
         open (1002, file='out/'//char(fn))
         write (1002, ('(3(A, 2E20.13, 3X))')) '#omega=', ome, 'alpha=', alf, 'beta=', bat
         write (1002, *) 'variables=y, |u|, u_r, u_i, |v|, v_r, v_i, &
   &                 |w|, w_r, w_i, |T|, T_r, T_i, |p|, p_r, p_i, &
   &                  |rho|, rho_r, rho_i, x, i'
         do j = 1, jn
            write (1002, '(20E20.7, 2X, I0)') yy(ix, j), &
               abs(phi(2, j)), phi(2, j)%re, phi(2, j)%im, &
               abs(phi(3, j)), phi(3, j)%re, phi(3, j)%im, &
               abs(phi(4, j)), phi(4, j)%re, phi(4, j)%im, &
               abs(phi(5, j)), phi(5, j)%re, phi(5, j)%im, &
               abs(phi(6, j)), phi(6, j)%re, phi(6, j)%im, &
               abs(phi(1, j)), phi(1, j)%re, phi(1, j)%im, &
               xx(ix, yb + 1), ix
         end do
         close (1002)
      end if

      call MPI_bcast(maxu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call VecScale(V, 1.0d0/maxu, ierr)

      call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'out/disturbance.petsc', FILE_MODE_WRITE, &
                                 Viewer, ierr)
      call VecView(V, Viewer, ierr)
      call PetscViewerDestroy(Viewer, ierr)

      call VecScatterDestroy(ctx, ierr)
      call VecDestroy(v_all, ierr)
      ! call VecDestroy(V, ierr)

   end subroutine write_lst_function

end module mod_lst_receptivity
