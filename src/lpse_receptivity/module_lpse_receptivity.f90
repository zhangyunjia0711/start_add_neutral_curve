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

      print *, 'initialize_lpse_receptivity.'

   end subroutine initialize_lpse_receptivity

   subroutine solve_lpse_receptivity()

      implicit none

      call solve_lpse()
      call solve_alpse()

      call test_orth()

      print *, 'solve_lpse_receptivity.'

   end subroutine solve_lpse_receptivity

   subroutine test_orth()

      implicit none

   end subroutine test_orth

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
      complex(REAL64) :: A(5, 5, jn)
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
         do j=0, jn-1
            idxm = 0; idxn = 0
            idxm(MatStencil_i, 1) = j; idxn(MatStencil_i, 1) = j
            call MatSetValuesBlockedStencil(MatA(i), 1, idxm, 1, idxn, transpose(A(:, :, j+1)), INSERT_VALUES, ierr)
         enddo
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

      call MatCopy(MatA(istart+1), MatA(istart), SAME_NONZERO_PATTERN, ierr)

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
      do i = istart - 1, iend + 1
         call DMRestoreGlobalVector(DM_local, lpse_function(i), ierr)
         call DMRestoreGlobalVector(DM_local, alpse_function(i), ierr)
      end do

      do i= istart, iend
         call MatDestroy(MatA(i), ierr)
      enddo

      deallocate (lpse_function, alpse_function)
      deallocate (MatA)

      call DMDestroy(DM_local, ierr)

      print *, 'finalize_lpse_receptivity.'

   end subroutine finalize_lpse_receptivity

end module mod_lpse_receptivity
