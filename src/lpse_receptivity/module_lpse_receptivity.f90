#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>
module mod_lpse_receptivity

   use global_parameters
   use mod_lpse_eqn
   use mod_alpse_eqn
   use mod_lst_grid_profile
   use mod_lst_baseflow_profile
   implicit none

   public :: lpse_receptivity
#include "../include/start.f90"

   complex(REAL64), allocatable :: lpse_alpha(:), lnA(:)
   complex(REAL64) :: ome, bat
   Vec, allocatable :: lpse_function(:), alpse_function(:)
   ! complex(REAL64), allocatable :: A(:, :, :, :)
   Mat, allocatable :: matA(:)

   DM :: DM_local
   type(type_lpse_eqn) :: lpse_eqn
   type(type_alpse_eqn) :: alpse_eqn
   complex(REAL64), allocatable :: lambda(:)

   Vec :: source_force
   complex(REAL64), allocatable :: wall_force(:, :)

contains

   subroutine lpse_receptivity()

      implicit none

      call initialize_lpse_receptivity()
      call solve_lpse_receptivity()
      call finalize_lpse_receptivity()

   end subroutine lpse_receptivity

   subroutine initialize_lpse_receptivity()

      implicit none
      integer :: i

      call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_local, ierr)
      call DMSetUp(DM_local, ierr)

      call DMGetGlobalVector(meshDA, source_force, ierr); call VecZeroEntries(source_force, ierr)
      allocate (wall_force(4, istart:iend)); wall_force = 0.0d0

      allocate (lpse_alpha(istart:iend))

      ome = omega(1)
      bat = beta(1)
      lpse_alpha(istart) = alpha(1) + ZI*alpha(2)

      allocate (lnA(istart:iend)); lnA = 0.0d0

      allocate (lpse_function(istart - 1:iend + 1))
      allocate (alpse_function(istart - 1:iend + 1))

      do i = istart - 1, iend + 1
         call DMGetGlobalVector(DM_local, lpse_function(i), ierr); call VecZeroEntries(lpse_function(i), ierr)
         call DMGetGlobalVector(DM_local, alpse_function(i), ierr); call VecZeroEntries(alpse_function(i), ierr)
      end do

      ! allocate (A(5, 5, 0:jn - 1, istart:iend))
      allocate (MatA(istart:iend))
      do i = istart, iend
         call DMCreateMatrix(DM_local, MatA(i), ierr); call MatZeroEntries(MatA(i), ierr)
      end do

      allocate (lambda(istart:iend)); lambda = 0.0d0

      print *, 'initialize_lpse_receptivity.'

   end subroutine initialize_lpse_receptivity

   subroutine solve_lpse_receptivity()

      implicit none

      call solve_lpse()
      call solve_alpse()

      call test_orth()
      call load_force_info()
      call compute_receptivity()

      print *, 'solve_lpse_receptivity.'

   end subroutine solve_lpse_receptivity

   subroutine load_force_info()
      implicit none

      print *, 'load_source_force'
      print *, 'load_bc_force'

   end subroutine load_force_info

   subroutine compute_receptivity()
      implicit none
      integer :: i

      type(type_lst_grid_profile) :: grid_profile
      complex(REAL64) :: source_force_term_part(istart:iend, 0:jn - 1)
      complex(REAL64) :: wall_term_part(istart:iend)
      complex(REAL64) :: freestream_term_part(istart:iend)

      call compute_source_term_part(source_force_term_part, alpse_function, source_force, lnA)
      call compute_wall_term_part(wall_term_part, alpse_function, wall_force)
      do i = istart, iend
         lambda(i) = 0.0
      end do

   end subroutine compute_receptivity

   subroutine compute_source_term_part(source_force_term_part, alpse_function, source_force, lnA)
      implicit none
      complex(REAL64), intent(inout) :: source_force_term_part(:, :)
      Vec, intent(in) :: alpse_function(:)
      Vec, intent(in) :: source_force
      complex(REAL64), intent(in) :: lnA(:)

      complex(REAL64), pointer :: array(:, :)
      complex(REAL64), pointer :: array2D(:, :, :)

      integer :: i, j

      call DMDAVecGetArrayF90(meshDA, source_force, array2D, ierr)
      do i = 0, iend - istart
         call DMDAVecGetArrayReadF90(DM_local, alpse_function(i + istart), array, ierr)
         do j = 0, jn - 1
            source_force_term_part(i + 1, j + 1) = dot_product((array(:, j)), array2D(:, i, j))/exp(lnA(i + 1))
         end do
         call DMDAVecRestoreArrayReadF90(DM_local, alpse_function(i + istart), array, ierr)
      end do
      call DMDAVecRestoreArrayF90(meshDA, source_force, array2D, ierr)

   end subroutine compute_source_term_part

   subroutine compute_wall_term_part(wall_term_part, alpse_function, wall_force)

      implicit none
      complex(REAL64), intent(inout) :: wall_term_part(:)
      Vec, intent(in) :: alpse_function(:)
      complex(REAL64), intent(in) :: wall_force(:, :)
      integer :: i
      complex(REAL64), pointer :: array(:, :)
      complex(REAL64) :: u_w, v_w, w_w, T_w
      complex(REAL64) :: rho_star_w, du_star_w, dv_star_w, dw_star_w, dT_star_w
      real(REAL64) :: mu0

      do i = 1, iend - istart + 1
         associate (T => qq(5, i + istart - 1, 1), C1 => 110.4D0)
            mu0 = T*sqrt(T)*(1.0d0 + C1/Te)/(T + C1/Te)/Re
         end associate

         u_w = wall_force(1, i); v_w = wall_force(2, i); w_w = wall_force(3, i); T_w = wall_force(4, i)

         call DMDAVecGetArrayReadF90(DM_local, alpse_function(i + istart), array, ierr)
         rho_star_w = array(0, 0)
         du_star_w = (-1.0d0*array(1, 2) + 4.0d0*array(1, 1) - 3.0d0*array(1, 0))/2.0d0*eta_y(i + istart - 1, 1)
         dv_star_w = (-1.0d0*array(2, 2) + 4.0d0*array(2, 1) - 3.0d0*array(2, 0))/2.0d0*eta_y(i + istart - 1, 1)
         dw_star_w = (-1.0d0*array(3, 2) + 4.0d0*array(3, 1) - 3.0d0*array(3, 0))/2.0d0*eta_y(i + istart - 1, 1)
         dT_star_w = (-1.0d0*array(4, 2) + 4.0d0*array(4, 1) - 3.0d0*array(4, 0))/2.0d0*eta_y(i + istart - 1, 1)
         call DMDAVecRestoreArrayReadF90(DM_local, alpse_function(i + istart), array, ierr)
         wall_term_part(i) = qq(1, i + istart - 1, 1)*rho_star_w*v_w + &
                             mu0*du_star_w*u_w + &
                             3.0d0/4.0d0*mu0*dv_star_w*v_w + &
                             mu0*dw_star_w*w_w + &
                             mu0/Pr*dT_star_w*T_w
      end do

   end subroutine compute_wall_term_part

   subroutine test_orth()

      implicit none
      Vec :: v_tmp
      type(type_lst_grid_profile) :: grid_profile
      integer :: i
      complex(REAL64) :: q
      REAL(REAL64) :: norm_inf1, norm_inf2

      call VecDuplicate(lpse_function(1), v_tmp, ierr)
      open (999, file='out/q.plt')
      do i = istart, iend
         call MatMult(MatA(i), lpse_function(i), v_tmp, ierr)
         grid_profile = assemble_lst_grid_profile(i)
         q = inner_dot_y(alpse_function(i), v_tmp, grid_profile%y)
         if (abs(q) > 1.0d2 .or. isnan(abs(q))) then
            call VecNorm(v_tmp, NORM_INFINITY, norm_inf1, ierr)
            call VecNorm(alpse_function(i), NORM_INFINITY, norm_inf2, ierr)
            print *, norm_inf1, norm_inf2
            error stop 'q is nan'
         end if
         write (999, '(2X, I4, 2X, 3(F20.8, 2X))') i, q, abs(q)
         ! write(999, *) i, q, abs(q)
      end do
      close (999)

   end subroutine test_orth

   function inner_dot_y(v1, v2, y) result(q)

      use stdlib_quadrature, only: trapz_weights
      implicit none

      Vec, intent(in) :: v1, v2
      real(REAL64), intent(in) :: y(:)
      complex(REAL64) :: q

      Vec :: tmp
      complex(REAL64), pointer:: array(:)
      integer :: j
      complex(REAL64) :: tmp_(0:jn - 1)

      call VecDuplicate(v1, tmp, ierr)
      call VecCopy(v1, tmp, ierr)
      call VecConjugate(tmp, ierr)
      call VecPointwiseMult(tmp, tmp, v2, ierr)
      call VecGetArrayReadF90(tmp, array, ierr)

      do j = 0, jn - 1
         tmp_(j) = sum(array(j*5 + 1:j*5 + 5))
      end do
      q = sum(trapz_weights(y)*tmp_)
      call VecRestoreArrayReadF90(tmp, array, ierr)

   end function inner_dot_y

   subroutine solve_lpse()

      implicit none

      call lpse_eqn%initial(jn, Re, Ma, Te, yb, ye)
      lpse_eqn%isothermal = isothermal
      call lpse_eqn%set_parameter("ome", omega(1))
      call lpse_eqn%set_parameter("bat", beta(1))
      call lpse_eqn%set_parameter("sigma_p", sigma_p)
      call lpse_eqn%set_parameter("sigma_control", sigma_control)
      call lpse_eqn%set_parameter("rd_factor", rd_factor)
      call lpse_eqn%set_parameter("Pr", Pr)

      call lpse_eqn%set_parameter("solve_alpha", SOLVE_ALPHA_BOTH)
      call solve_lpse_native()
      call lpse_eqn%finalize()

      print *, 'solve_lpse'

   end subroutine solve_lpse

   subroutine solve_alpse()
      implicit none

      call alpse_eqn%initial(jn, Re, Ma, Te, yb, ye)
      alpse_eqn%isothermal = isothermal
      call alpse_eqn%set_parameter("ome", omega(1))
      call alpse_eqn%set_parameter("bat", beta(1))
      call alpse_eqn%set_parameter("sigma_p", sigma_p)
      call alpse_eqn%set_parameter("sigma_control", sigma_control)
      call alpse_eqn%set_parameter("rd_factor", rd_factor)
      call alpse_eqn%set_parameter("Pr", Pr)

      call alpse_eqn%set_parameter("solve_alpha", SOLVE_ALPHA_NONE)
      call solve_lpse_adjoint()
      call alpse_eqn%finalize()
      print *, 'solve_alpse'
   end subroutine solve_alpse

   subroutine set_lpse_ic()

      use mod_lst_eqn
      implicit none
      integer :: i
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile
      integer :: iloc
      type(type_lst_eqn) :: lst

      !!! compute from lst directly
      iloc = istart
      call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call lst%set_parameter("Pr", Pr)
      lst%isothermal = isothermal

      grid_profile = assemble_lst_grid_profile(iloc)
      baseflow_profile = assemble_lst_baseflow_profile(iloc)

      call lst%set_grid(grid_profile)
      call lst%set_baseflow(baseflow_profile)
      call lst%setup()

      call lst%solve(ome, lpse_alpha(iloc), bat)
      call lst%get_eigf(lpse_function(iloc))
      call lst%finalize()

   end subroutine set_lpse_ic

   subroutine set_alpse_ic()

      use mod_alst_eqn
      implicit none
      integer :: i
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile
      integer :: iloc
      type(type_alst_eqn) :: alst

   !!! compute from lst directly
      iloc = iend
      call alst%set_parameter("lst_mode", LST_SPATIAL_MODE)
      call alst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
      call alst%set_parameter("Pr", Pr)
      alst%isothermal = isothermal

      grid_profile = assemble_lst_grid_profile(iloc)
      baseflow_profile = assemble_lst_baseflow_profile(iloc)

      call alst%set_grid(grid_profile)
      call alst%set_baseflow(baseflow_profile)
      call alst%setup()

      call alst%solve(ome, lpse_alpha(iloc), bat)
      call alst%get_eigf(alpse_function(iloc))
      call alst%finalize()

   end subroutine set_alpse_ic

   subroutine solve_lpse_native()

      implicit none
      integer :: i, j
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile
      complex(REAL64) :: A(5, 5, 0:jn - 1)
      integer :: idxm(4, 1), idxn(4, 1)

      grid_profile = assemble_lst_grid_profile(istart)
      baseflow_profile = assemble_lst_baseflow_profile(istart)

      call set_lpse_ic()
      print *, 'istart=', istart
      call lpse_eqn%set_ic(lpse_function(istart))
      do i = istart + 1, iend
         grid_profile = assemble_lst_grid_profile(i)
         baseflow_profile = assemble_lst_baseflow_profile(i)

         call lpse_eqn%set_grid(grid_profile)
         call lpse_eqn%set_baseflow(baseflow_profile)
         call lpse_eqn%setup()
         call lpse_eqn%get_matrix(lpse_alpha(i), "A_0", A)
         do j = 0, jn - 1
            idxm = -1; idxn = -1
            idxm(MatStencil_i, 1) = j; idxn(MatStencil_i, 1) = j
            call MatSetValuesBlockedStencil(MatA(i), 1, idxm, 1, idxn, transpose(A(:, :, j)), INSERT_VALUES, ierr)
         end do
         call MatAssemblyBegin(MatA(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(MatA(i), MAT_FINAL_ASSEMBLY, ierr)
         ! call solve_lpse_eqn(grid_profile, baseflow_profile, ome, alf, bat, i, phi(i), phi(i-1), phi(max(1,i-2)))
         lpse_alpha(i) = lpse_alpha(i - 1)
         call VecCopy(lpse_function(i - 1), lpse_function(i), ierr)
         call lpse_eqn%solve(lpse_alpha(i), lpse_function(i), lpse_function(i - 1), lpse_function(i - 2))

         lnA(i) = lnA(i - 1) + 0.5d0*(lpse_alpha(i) + lpse_alpha(i - 1))/grid_profile%xi_x(yb + 1)
         print *, 'ix=', i, 'is done.'
         print *, 'alpha=', lpse_alpha(i), "bat=", bat
      end do
      call PetscPrintf(PETSC_COMM_WORLD, 'Solve LPSE done!\n', ierr)

      call MatCopy(MatA(istart + 1), MatA(istart), SAME_NONZERO_PATTERN, ierr)

   end subroutine solve_lpse_native

   subroutine solve_lpse_adjoint()

      implicit none
      integer :: i
      type(type_lst_grid_profile) :: grid_profile
      type(type_lst_bf_profile) :: baseflow_profile

      call set_alpse_ic()

      grid_profile = assemble_lst_grid_profile(iend)
      baseflow_profile = assemble_lst_baseflow_profile(iend)

      call alpse_eqn%set_ic(alpse_function(iend))

      do i = iend - 1, istart, -1
         grid_profile = assemble_lst_grid_profile(i)
         baseflow_profile = assemble_lst_baseflow_profile(i)
         ! ome=wave_param(i)%ome
         ! alf=wave_param(i)%alf
         ! bat=wave_param(i)%bat

         call alpse_eqn%set_grid(grid_profile)
         call alpse_eqn%set_baseflow(baseflow_profile)
         call alpse_eqn%setup()

         call VecCopy(alpse_function(i + 1), alpse_function(i), ierr)
         call alpse_eqn%solve(lpse_alpha(i), alpse_function(i), alpse_function(i + 1), alpse_function(i + 2))

         print *, 'ix=', i, 'is done.'
         print *, 'alpha=', lpse_alpha(i), "bat=", bat
      end do
      call PetscPrintf(PETSC_COMM_WORLD, 'Solve APSE done!\n', ierr)

   end subroutine solve_lpse_adjoint

   subroutine finalize_lpse_receptivity()

      implicit none
      integer :: i

      deallocate (lpse_alpha, lnA)
      deallocate (lambda)
      do i = istart - 1, iend + 1
         call DMRestoreGlobalVector(DM_local, lpse_function(i), ierr)
         call DMRestoreGlobalVector(DM_local, alpse_function(i), ierr)
      end do

      do i = istart, iend
         call MatDestroy(MatA(i), ierr)
      end do

      deallocate (lpse_function, alpse_function)
      deallocate (MatA)

      call DMDestroy(DM_local, ierr)

      call DMRestoreGlobalVector(meshDA, source_force, ierr)
      call DMDestroy(meshDA, ierr)
      deallocate (wall_force)

      print *, 'finalize_lpse_receptivity.'

   end subroutine finalize_lpse_receptivity

end module mod_lpse_receptivity
