!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-29 15:17:05
!> @LastEditTime: 2023-03-07 09:30:48
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lpse/lpse_eqn.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ module lpse equations
!! \[
!! \phi=\phi(x, y, z, t) \\
!!   G_0 \frac{\partial \phi}{\partial t}+ A_0 \frac{\partial \phi}{h_1\partial x}
!! + B_0 \frac{\partial \phi}{\partial y}+ C_0 \frac{\partial \phi}{h_3\partial z}
!! + D_0 \phi= \\
!!  V_{0xx} \frac{\partial^2 \phi}{h_1^2 \partial x^2}
!! +V_{0xy} \frac{\partial^2 \phi}{h_1   \partial x\partial y}
!! +V_{0xz} \frac{\partial^2 \phi}{h_1h_3\partial x\partial z}
!! +V_{0yy} \frac{\partial^2 \phi}{      \partial y^2}
!! +V_{0yz} \frac{\partial^2 \phi}{h_3   \partial y\partial z}
!! +V_{0zz} \frac{\partial^2 \phi}{h_3^2 \partial z^2}
!! \]

!! \[
!! \phi=\hat\phi(y) {\rm e}^{i(\int_{x0}^ x \alpha d\xi +\beta z - \omega t)} \\
!!  L\hat\phi=0\\
!!  L=A\frac{\partial}{\partial x} + B\frac{\partial}{\partial y} + V\frac{\partial^2}{\partial y^2}-D
!! \]
module mod_lpse_eqn

   use petsc

   use mod_lst_baseflow_profile
   use mod_lst_grid_profile

   use iso_fortran_env, only: REAL64, REAL32

   use mod_lpse_solver_abstract, only: type_lpse_abstract

   implicit none

   private
   public :: type_lpse_eqn
   public :: A0, V0xx


#include "../include/start.f90"

   logical :: lpse_eqn_has_initialize = .FALSE.
   integer :: yb, ye

   PetscErrorCode :: ierr


   real(REAL64), allocatable, dimension(:, :, :) :: G0, A0, B0, C0, D0
   real(REAL64), allocatable, dimension(:, :, :) :: V0xx, V0xy, V0xz
   real(REAL64), allocatable, dimension(:, :, :) :: V0yy, V0yz, V0zz
   real(REAL64), allocatable, dimension(:, :, :) :: Aplus
   real(REAL64), allocatable, dimension(:, :) :: modify_Px_A

   type, extends(type_lpse_abstract) :: type_lpse_eqn

      integer :: iterator_times
      logical :: solved_lpse_matrix = .FALSE.
      complex(real64), allocatable :: A0_rec(:, :, :)

   contains

      procedure :: initial => lpse_intial
      procedure :: set_ic => lpse_set_ic
      procedure :: setup => lpse_setup
      procedure :: solve => lpse_solve
      procedure :: get_matrix => get_matrix_lpse
      procedure :: finalize => lpse_eqn_finalize

      procedure :: solve_baseflow_matrix
      procedure :: solve_bf_matrix_j
      procedure :: solve_lpse_matrix
      procedure :: solve_lpse_matrix_j
      procedure :: solve_lpse_matrix_j0
      procedure :: solve_lpse_matrix_jn
      procedure :: solve_lpse_matrix_bc
      procedure :: solve_lpse_phi
      procedure :: solve_lpse_bc
      procedure :: solve_alpha_local_new
      procedure :: solve_lpse_bc_j
      procedure :: solve_rhs_new
      procedure :: solve_left_matrix

   end type type_lpse_eqn

contains

   subroutine lpse_intial(this, jn_, Re, Ma, Te, yb_, ye_)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this
      integer, intent(in) :: jn_
      real(real64), intent(in) ::  Re, Ma, Te
      integer, intent(in), optional :: yb_, ye_

      associate (jb => this%jb, je => this%je, jn => this%jn)
         jn = jn_
         yb = yb_; ye = ye_
         jb = yb + 1; je = ye + 1

         this%Re = Re
         this%Ma = Ma
         this%Te = Te

         allocate (G0(5, 5, jb:je)); G0 = 0.0d0
         allocate (A0(5, 5, jb:je)); A0 = 0.0d0
         allocate (B0(5, 5, jb:je)); B0 = 0.0d0
         allocate (C0(5, 5, jb:je)); C0 = 0.0d0
         allocate (D0(5, 5, jb:je)); D0 = 0.0d0
         allocate (V0xx(5, 5, jb:je)); V0xx = 0.0d0
         allocate (V0xy(5, 5, jb:je)); V0xy = 0.0d0
         allocate (V0xz(5, 5, jb:je)); V0xz = 0.0d0
         allocate (V0yy(5, 5, jb:je)); V0yy = 0.0d0
         allocate (V0yz(5, 5, jb:je)); V0yz = 0.0d0
         allocate (V0zz(5, 5, jb:je)); V0zz = 0.0d0
         allocate (Aplus(5, 5, jb:je)); Aplus = 0.0d0

         ! associate(D_2 => this%D_2, B_0=>this%B_0, D_0=>this%D_0, V_0=>this%V_0)
         allocate (this%B_0(5, 5, jb:je)); this%B_0 = 0.0d0
         allocate (this%D_0(5, 5, jb:je)); this%D_0 = 0.0d0
         allocate (this%V_0(5, 5, jb:je)); this%V_0 = 0.0d0
         allocate (this%D_2(5, 5, jb:je)); this%D_2 = 0.0d0
         allocate (this%G_0(5, 5, jb:je)); this%G_0 = 0.0d0
         allocate (this%A0_rec(5, 5, jb:je)); this%A0_rec = 0.0d0
         ! end associate

         allocate (modify_Px_A(2, jb:je)); modify_Px_A = 0.0d0

         lpse_eqn_has_initialize = .TRUE.

         call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, this%dm, ierr)
         call DMSetUp(this%dm, ierr)
         call DMCreateMatrix(this%dm, this%L_LPSE, ierr)
         call DMGetGlobalVector(this%dm, this%RHS, ierr); call VecZeroEntries(this%RHS, ierr)
         call VecDuplicate(this%RHS, this%phi, ierr); call VecZeroEntries(this%phi, ierr)

         allocate (this%phi0_array(5, jn)); this%phi0_array = 0.0d0
         allocate (This%phi1_array(5, jn)); this%phi1_array = 0.0d0
         allocate (this%phi2_array(5, jn)); this%phi2_array = 0.0d0

         call KSPCreate(PETSC_COMM_WORLD, this%solver, ierr)
         call KSPGetPC(this%Solver, this%SolverPC, ierr)
         call KSPSetType(this%Solver, KSPPREONLY, ierr)
         ! call KSPGetPC(this%Solver, this%SolverPC, ierr)
         call PCSetType(this%SolverPC, PCLU, ierr)
! #ifdef __linux__
!          call PCFactorSetMatSolverType(this%SolverPC, MATSOLVERMKL_CPARDISO, ierr)
! #elif __APPLE__
         call PCFactorSetMatSolverType(this%SolverPC, MATSOLVERSUPERLU_DIST, ierr)
! #endif
         call KSPSetFromOptions(this%Solver, ierr)
         ! call KSPView(Solver, PETSC_VIEWER_STDOUT_SELF, ierr)
         call PetscPrintf(PETSC_COMM_WORLD, 'allocate LPSE!\n', ierr)

         this%lambda = 0.0d0 !< initial lambda is 1st order.

      end associate

   end subroutine lpse_intial

   subroutine lpse_setup(this)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this

      call this%solve_baseflow_matrix()
      this%solved_lpse_matrix = .false.
      if(.not. ((this%cycle_m%re+1)<=epsilon(1.0d0))) this%bat = this%cycle_m * this%grid%drz

   end subroutine lpse_setup

   subroutine lpse_set_ic(this, phi)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this
      Vec, intent(inout) ::  phi
      PetscScalar, pointer :: tmp(:, :)

      call VecCopy(phi, this%phi, ierr)
      call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
      this%phi0_array = tmp
      call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)
      this%phi1_array = 0.0d0
      this%phi2_array = 0.0d0
      this%lambda = 0.0d0

   end subroutine lpse_set_ic

   subroutine lpse_solve(this, alf, phi0, phi1, phi2)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this
      complex(REAL64), intent(inout) :: alf
      Vec, intent(inout) ::  phi0
      Vec, intent(in), optional :: phi1, phi2
      PetscScalar, pointer :: tmp(:, :)
      complex(REAL64) :: alf_old, delta

      call VecCopy(phi0, this%phi, ierr)
      call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
      this%phi0_array = tmp
      call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)

      if ((.not. present(phi1)) .and. (present(phi2))) then
         stop "Phi i-2 can't be given without Phi i-1."
      else if (((.not. present(phi1)) .and. (.not. present(phi2)))) then
         this%phi2_array = this%phi1_array
         this%phi1_array = this%phi0_array
      else if (present(phi1) .and. (.not. present(phi2))) then
         call DMDAVecGetArrayReadF90(this%dm, phi1, tmp, ierr)
         this%phi1_array = tmp
         call DMDAVecRestoreArrayReadF90(this%dm, phi1, tmp, ierr)
         this%lambda = 0.0d0
      else if (present(phi1) .and. present(phi2)) then
         call DMDAVecGetArrayReadF90(this%dm, phi1, tmp, ierr)
         this%phi1_array = tmp
         call DMDAVecRestoreArrayReadF90(this%dm, phi1, tmp, ierr)
         call DMDAVecGetArrayReadF90(this%dm, phi2, tmp, ierr)
         this%phi2_array = tmp
         call DMDAVecRestoreArrayReadF90(this%dm, phi2, tmp, ierr)
         this%lambda = 0.5d0
      end if

      if (sum(abs(this%phi2_array)) <= epsilon(1.0d0)) then  !phi2 is zero means 1 order is used
         this%lambda = 0.0d0
      else ! else 2 order is used
         this%lambda = 0.5d0
      end if

      ! this%bat = this%cycle_m*this%grid%drz

      alf_old = alf
      delta = 1.0d10
      ! read(*, *)
      this%iterator_times = -1

      ! print*, 'lpse'
      ! print*, alf
      ! print*, this%bat, this%cycle_m, this%grid%drz
      ! print*, this%ome
      
      do while (abs(delta) >= (epsilon(1.0e0)*100.0d0))

         call this%solve_lpse_matrix(alf)
         call this%solve_left_matrix()
         call this%solve_rhs_new()
         call this%solve_lpse_matrix_bc(alf)
         call this%solve_lpse_bc()
         call this%solve_lpse_phi()
         call this%solve_alpha_local_new(alf)
         ! print*, alf, alf_old, (alf-alf_old)
         delta = abs(alf - alf_old)
         alf_old = alf
         ! call VecWAXPY(dphi, (-1.0d0, 0.0d0), phi_old, phi_new, ierr)
         ! call VecNorm(dphi, NORM_2, delta, ierr)
         ! call VecEqual(phi_new, phi_old, flg, ierr)
         this%iterator_times = this%iterator_times + 1
      end do

      if (this%solve_alpha == SOLVE_ALPHA_BOTH) print *, "the iteration times is", this%iterator_times

      call VecCopy(this%phi, phi0, ierr)
      call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
      this%phi0_array = tmp
      call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)

   end subroutine lpse_solve

   subroutine solve_alpha_local_new(this, alf)

      use functional
      use mod_dx
      implicit none
      class(type_lpse_eqn) :: this
      complex(REAL64), intent(inout) :: alf
      integer :: l, j
      real(REAL64), dimension(this%jn) :: hx
      complex(REAL64), dimension(5, this%jn) :: dphidi, dphidj, dphidjj, dphidx
      complex(REAL64), dimension(this%jn) :: phi_ke_j, dphi_ke_j
      PetscScalar :: phi_ke, dphi_ke
      PetscScalar, pointer :: tmp(:, :)

      if (this%solve_alpha /= SOLVE_ALPHA_NONE) then
         phi_ke = 0.0d0; dphi_ke = 0.0d0

         call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
         this%phi0_array = tmp
         call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)
         dphidi = this%phi0_array - this%phi1_array
        !!dphi_dj, dphi_djj
         do l = 1, 5
            call dx_fdm4nd_center(this%phi0_array(l, :), this%jb, this%je, dphidj(l, :), this%jb, this%je)
            ! call dxx_fdm4nd_center(this%phi0_array(l, :), jb, je, dphidjj(l, :), jb, je)
         end do

         hx(:) = 1.0d0 + this%grid%drx*this%grid%y(:)

         do l = 1, 5
            dphidx(l, :) = dphidi(l, :)*this%grid%xi_x(:) + dphidj(l, :)*this%grid%eta_x(:)
         end do

         do j = 1, this%jn
            phi_ke_j(j) = (this%phi0_array(2, j)*conjg(this%phi0_array(2, j)) + &
                           this%phi0_array(3, j)*conjg(this%phi0_array(3, j)) + &
                           this%phi0_array(4, j)*conjg(this%phi0_array(4, j)))/hx(j)
            dphi_ke_j(j) = (dphidx(2, j)*conjg(this%phi0_array(2, j)) + &
                            dphidx(3, j)*conjg(this%phi0_array(3, j)) + &
                            dphidx(4, j)*conjg(this%phi0_array(4, j)))/hx(j)
         end do

         phi_ke = Trap(phi_ke_j, this%grid%y)
         dphi_ke = Trap(dphi_ke_j, this%grid%y)

         alf = alf - this%rd_factor*ZI*1.0d0/phi_ke*dphi_ke
         ! print*, phi_ke, dphi_ke
         if (this%solve_alpha == SOLVE_ALPHA_REAL) then
            alf = real(alf)
         end if
      else
         alf = alf
      end if

   end subroutine solve_alpha_local_new

   subroutine solve_rhs_new(this)

      use functional
      implicit none
      class(type_lpse_eqn) :: this

      integer :: bs
      PetscScalar :: rhs_j(5, this%jb:this%je), box(5*(this%je - this%jb + 1))
      integer :: j
      PetscScalar, dimension(5, this%jb:this%je) :: phi_iminus1, phi_iminus2

      !! RHS=lambda*xi_x*A(phi(i-2)-phi(i-1))+(1+lambda)*xi_x*A*phi(i-1)
      associate (RHS => this%RHS, lambda => this%lambda, A_0 => this%D_2, jb => this%jb, je => this%je, jn => this%jn)

         do j = jb, je
            rhs_j(:, j) = lambda*matmul(A_0(:, :, j), (this%phi1_array(:, j) - this%phi2_array(:, j))) &
                          + (1.0d0 + lambda)*matmul(A_0(:, :, j), this%phi1_array(:, j))
         end do

         ! print*, 'rhs=', sum(abs(rhs_j))
         ! print*, 'phi_m1=', sum(abs(this%phi1_array))
         ! print*, 'phi_m2=', sum(abs(this%phi2_array))

         box = reshape(rhs_j, [5*(je - jb + 1)])
         call VecSetValuesBlocked(RHS, je - jb + 1, arange(jb - 1, je - 1), box, INSERT_VALUES, ierr)
         call VecAssemblyBegin(RHS, ierr)
         call VecAssemblyEnd(RHS, ierr)
         call VecGetBlockSize(RHS, bs, ierr)
         ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)
         ! read(*, *)

      end associate

   end subroutine solve_rhs_new

   subroutine get_matrix_lpse(this, alf, mat_name, matrix)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this
      complex(REAL64), intent(in) :: alf
      character(len=*), intent(in) :: mat_name
      complex(REAL64), intent(out) :: matrix(5, 5, this%jn)

      select case (trim(mat_name))
      case ('A_0')
         if (.not. this%solved_lpse_matrix) call this%solve_lpse_matrix(alf)
         ! matrix = this%D_2
         matrix = this%A0_rec
         ! if(isnan(matrix(2,2,2)%re)) then
         !    print*, 'in lpse_eqn.f90'
         !    block
         !       integer :: l1, l2, j
         !       do j=2, 2
         !       do l1 = 1, 5
         !          do l2=1, 5
         !             print*, l1, l2, this%D_2(l1, l2, j)
         !          enddo
         !       enddo
         !    enddo
         !    endblock
         !    stop
         ! endif
      case default
         stop "A wrong matrix name is given in get_matrix@lpse_eqn."
      end select

   end subroutine get_matrix_lpse

   function Trap(f, y) result(integal)

      complex(REAL64), intent(in) :: f(:)
      Real(REAL64), intent(in) :: y(:)
      complex(REAL64) :: integal
      integer :: i

      integer :: num(1)

      num = ubound(f)
      integal = 0.0d0
      do i = 2, num(1)
         integal = integal + (f(i) + f(i - 1))*0.5d0*(y(i) - y(i - 1))
      end do

   end function Trap

   subroutine solve_lpse_phi(this)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this
      integer :: its
      character(20) :: ksptype

      call KSPSetOperators(this%Solver, this%L_LPSE, this%L_LPSE, ierr)
      call KSPSetFromOptions(this%solver, ierr)
      call KSPSetup(this%Solver, ierr)
      call KSPSolve(this%Solver, this%RHS, this%phi, ierr)

   end subroutine solve_lpse_phi

   subroutine solve_baseflow_matrix(this)

      implicit none
      class(type_lpse_eqn) :: this
      integer :: j, j_b, j_e

        !! internal point
      call this%solve_bf_matrix_j()

   end subroutine solve_baseflow_matrix

   subroutine solve_lpse_matrix_bc(this, alf)

      implicit none
      class(type_lpse_eqn) :: this
      complex(REAL64), intent(in) :: alf

        !! b.c. point
      if (this%jb == 1) call solve_lpse_matrix_j0(this, this%jb, alf)
      if (this%je == this%jn) call solve_lpse_matrix_jn(this, this%je, alf)

   end subroutine solve_lpse_matrix_bc
   !
   subroutine solve_bf_matrix_j(this)

      implicit none
      class(type_lpse_eqn) :: this
      real(REAL64), parameter :: s = 0.0d0
      real(REAL64) :: u, ux, uy, uz, uyy, uxx, uxy, uxz
      real(REAL64) :: v, vx, vy, vz, vyy, vxy, vyz
      real(REAL64) :: w, wx, wy, wz, wxx, wyy, wxz, wyz
      real(REAL64) :: rho, rhox, rhoy, rhoz
      real(REAL64) :: T, Tx, Ty, Tz, Tyy, Txx, Tzz, Tyz

      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      real(REAL64) :: divu, divuy, divux, divuz
      real(REAL64) :: p11, p12, p13, p22, p23, p33, p12y, p22y, p23y, phi
      real(REAL64) :: p11x, p12x, p13x, p13z, p23z, p33z
      real(REAL64) :: pe, cm, ck
      real(REAL64) :: rm, rmx, rmy, rmz, rmt, rmtt, rmtx, rmtz, rmty
      real(REAL64) :: rk, rkt, rktt, rky, rkz, rkty
      real(REAL64) :: cv, cp, rek
      real(REAL64) :: omega_p, Mx

      integer :: j
      real(REAL64), parameter :: Z2D3 = 2.d0/3.d0, D3 = 1.d0/3.d0

      integer :: jb_local, je_local

      associate (grid_profile => this%grid, bf_profile => this%baseflow, Re => this%Re, Ma => this%Ma, &
                 Te => this%Te, Pr => this%Pr, C1 => this%C1, C2 => this%C2, GAMMA => this%GAMMA)
         jb_local = this%jb !max(2, jb)
         je_local = this%je !min(je, jn-1)

         do j = jb_local, je_local

            G0(:, :, j) = 0.0d0
            A0(:, :, j) = 0.0d0
            B0(:, :, j) = 0.0d0
            C0(:, :, j) = 0.0d0
            D0(:, :, j) = 0.0d0
            V0xx(:, :, j) = 0.0d0
            V0xy(:, :, j) = 0.0d0
            V0xz(:, :, j) = 0.0d0
            V0yy(:, :, j) = 0.0d0
            V0yz(:, :, j) = 0.0d0
            V0zz(:, :, j) = 0.0d0
            modify_Px_A(2, j) = 0.0d0

            y = grid_profile%y(j)
            drx = grid_profile%drx
            drz = grid_profile%drz

            select case (this%Coordinate_System)
            case (COMOVING_CS)
               hx = 1.d0 + drx*y
               hx_y = drx
               d12 = hx_y/hx

               hz = 1.d0 + drz*y
               hz_y = drz
               d32 = hz_y/hz
            case (CYLINDRICAL_CS)
               hx = 1.d0
               hx_y = 0.0d0
               d12 = 0.0d0

               hz = y
               if (abs(hz) <= epsilon(1.0d0)) hz = epsilon(1.0d0)
               hz_y = 1.0d0
               d32 = hz_y/hz
            end select

            d132 = d12 + d32

            rho = bf_profile%q0(1, j)
            rhox = bf_profile%dq0x(1, j)/hx
            rhoy = bf_profile%dq0y(1, j)

            u = bf_profile%q0(2, j)
            ux = bf_profile%dq0x(2, j)/hx
            uy = bf_profile%dq0y(2, j)
            uxx = bf_profile%dq0xx(2, j)/hx/hx
            uyy = bf_profile%dq0yy(2, j)
            uxy = bf_profile%dq0xy(2, j)/hx

            v = bf_profile%q0(3, j)
            vx = bf_profile%dq0x(3, j)/hx
            vy = bf_profile%dq0y(3, j)
            vyy = bf_profile%dq0yy(3, j)
            vxy = bf_profile%dq0xy(3, j)/hx

            w = bf_profile%q0(4, j)
            wx = bf_profile%dq0x(4, j)/hx
            wy = bf_profile%dq0y(4, j)
            wxx = bf_profile%dq0xx(4, j)/hx/hx
            wyy = bf_profile%dq0yy(4, j)

            T = bf_profile%q0(5, j)
            Tx = bf_profile%dq0x(5, j)/hx
            Ty = bf_profile%dq0y(5, j)
            Txx = bf_profile%dq0xx(5, j)/hx/hx
            Tyy = bf_profile%dq0yy(5, j)

            rhoz = 0.0d0/hz
            uz = 0.0d0/hz
            vz = 0.0d0/hz
            wz = 0.0d0/hz
            Tz = 0.0d0/hz
            Tzz = 0.d00/hz/hz
            uxz = 0.0d0/hx/hz
            vyz = 0.0d0/hz
            wxz = 0.0d0/hx/hz
            wyz = 0.0d0/hz
            Tyz = 0.0d0/hz

            ! print*, j
            ! print*, Tx, ux, vx, wx, rhox, v
            ! read(*, *)

            ! if(j==1)then
            !     open(99999,file='baseflow_bin.dat', form='unformatted')
            ! endif
            ! if ( j>=2 .and. j<=501) then
            !     read(99999)rho, rhox, rhoy, u, ux, uy, uyy, uxy, v, vy, vyy, vxy, &
            !     w, wy, wyy, T, Tx, Ty, Tyy
            ! endif
            ! if(j==501)then
            !
            !     close(99999)
            ! endif
            !
            pe = 1.d0/(GAMMA*Ma*Ma)
            cm = C1/Te
            ck = C2/Te
            cv = pe/(GAMMA - 1.d0)
            cp = cv*GAMMA
            rek = re*pr/cp

            rm = T*sqrt(T)*(1.0d0 + cm)/(T + cm)/re
            rmt = rm*(1.5d0/T - 1.0d0/(T + cm))
            rmtt = rmt*(1.5d0/T - 1.0d0/(T + cm)) - rm*(1.5d0/T**2 - 1.0d0/(T + cm)**2)
            rmx = rmt*tx
            rmy = rmt*ty
            rmz = 0.0d0
            rmty = rmtt*ty
            rmtx = rmtt*tx
            rmtz = 0.0d0

            rk = T*sqrt(T)*(1.d0 + ck)/(T + ck)/rek
            rkt = rk*(1.5d0/T - 1.0d0/(T + ck))
            rktt = rkt*(1.5d0/T - 1.0d0/(T + ck)) - rk*(1.5d0/T**2 - 1.0d0/(T + ck)**2)
            rky = rkt*ty
            rkz = 0.0d0
            rkty = rktt*ty

            divu = ux + vy + v*d132
            divux = uxx + vxy + vx*d132
            divuy = uxy + vyy + vy*d132 - v*(d12**2 + d32**2)
            divuz = uxz + vyz + vz*d132  !! A fatal error occurred in divuz compuation..

            p11 = 2.0d0*(ux + d12*v) - z2d3*divu
            p12 = vx + uy - d12*u
            p13 = wx + uz
            p22 = 2.0d0*vy - z2d3*divu
            p23 = wy + vz - d32*w
            p33 = 2.0d0*(d32*v) - z2d3*divu

            p11x = 2.0d0*(uxx + d12*vx) - z2d3*divux
            p12x = uxy - d12*ux
            p13x = wxx + uxz
            p13z = wxz
            p23z = wyz - d32*wz
            p33z = 2.0d0*d32*vz - z2d3*divuz
            p12y = vxy + uyy - d12*uy + d12**2*u
            p22y = 2.0d0*vyy - z2d3*divuy
            p23y = wyy - d32*wy + d32**2*w

            phi = 2.0d0*(ux + d12*v)**2 + 2.0d0*vy**2 + 2.0d0*(d32*v)**2 + p12**2 + p13**2 + p23**2 - z2d3*divu**2

            Mx = u/sqrt(T)*Ma
            ! omega_p=this%sigma_p*min(1.0d0, GAMMA*Mx**2/(1.0d0+(GAMMA-1.0d0)*Mx**2))+this%sigma_control
            omega_p = this%sigma_p*min(1.0d0, Mx**2) + this%sigma_control
            ! print*, omega_p, j

            G0(1, 1, j) = 1.0d0
            G0(2, 2, j) = rho
            G0(3, 3, j) = rho
            G0(4, 4, j) = rho
            G0(5, 1, j) = -pe*t
            G0(5, 5, j) = rho*cp - pe*rho

            A0(1, 1, j) = u
            A0(1, 2, j) = rho
            A0(2, 1, j) = pe*T
            A0(2, 2, j) = rho*u - z2d3*rmx*(2.0d0 + s)
            A0(2, 3, j) = -rmy - rm*(2.0d0*d12 + z2d3*(s + 0.5d0)*d132)
            A0(2, 4, j) = -rmz
            A0(2, 5, j) = pe*rho - rmt*p11

            A0(3, 2, j) = rm*(-z2d3*(s - 1.0d0)*d12*2.0 + d12) - Z2D3*rmy*(s - 1.0d0)
            A0(3, 3, j) = rho*u - rmx
            A0(3, 5, j) = -rmt*p12
            A0(4, 2, j) = -z2d3*rmz*(s - 1.0d0)
            A0(4, 4, j) = rho*u - rmx
            A0(4, 5, j) = -rmt*p13
            A0(5, 1, j) = -pe*u*T
            A0(5, 2, j) = -2.0d0*rm*p11
            A0(5, 3, j) = -2.0d0*rm*p12
            A0(5, 4, j) = -2.0d0*rm*p13
            A0(5, 5, j) = rho*cp*u - pe*u*rho - 2.0d0*Tx*rkt

            modify_Px_A(1, j) = pe*T*(-1.0d0 + omega_p)
            modify_Px_A(2, j) = pe*rho*(-1.0d0 + omega_p)
            !
            ! block
            !           if(j==149)then
            !             print*, 'check op'
            !             print*, A0(5, 5, j)
            !             print*, omega_p, Mx**2
            !             print*, modify_Px_A(:, j)
            !             read(*, *)
            !           endif
            ! endblock
            !
            B0(1, 1, j) = v
            B0(1, 3, j) = rho
            B0(2, 2, j) = rho*v - rmy - d32*rm
            B0(2, 3, j) = z2d3*(-rmx)*(s - 1)
            B0(2, 5, j) = -rmt*p12
            B0(3, 1, j) = pe*T
            B0(3, 2, j) = -rmx
            B0(3, 3, j) = rho*v - Z2D3*(s + 2.0d0)*rm*d132 - (2.0d0 + s)*Z2d3*rmy
            B0(3, 4, j) = -rmz
            B0(3, 5, j) = pe*rho - rmt*p22
            B0(4, 3, j) = z2d3*(-rmz)*(s - 1)
            B0(4, 4, j) = rho*v - rm*d132 - rmy
            B0(4, 5, j) = -rmt*p23
            B0(5, 1, j) = -pe*v*t
            B0(5, 2, j) = -2.0d0*rm*p12
            B0(5, 3, j) = -2.0d0*rm*p22
            B0(5, 4, j) = -2.0d0*rm*p23
            B0(5, 5, j) = rho*cp*v - pe*v*rho - 2.0d0*rky - rk*d132

            C0(1, 1, j) = w
            C0(1, 4, j) = rho
            C0(2, 2, j) = rho*w - rmz
            C0(2, 4, j) = -z2d3*rmx*(s - 1.0d0)
            C0(2, 5, j) = -rmt*p13
            C0(3, 3, j) = rho*w - rmz
            !C0(3, 4, j)=7.0d0/3.0d0*rm*d32+Z2D3*rmy
            C0(3, 4, j) = (-z2d3*(s - 1.0d0)*d32*2.0d0 + d32)*rm - z2d3*rmy*(s - 1.0d0)
            C0(3, 5, j) = -rmt*p23
            C0(4, 1, j) = pe*t
            C0(4, 2, j) = -rmx
            C0(4, 3, j) = -rmy - D3*(2.0d0*s + 1.0d0)*rm*d132 - 2.0d0*rm*d32
            C0(4, 4, j) = rho*w - z2d3*rmz*(2.0d0 + s)
            C0(4, 5, j) = pe*rho - rmt*p33
            C0(5, 1, j) = -pe*w*T
            C0(5, 2, j) = -2.0d0*rm*p13
            C0(5, 3, j) = -2.0d0*rm*p23
            C0(5, 4, j) = -2.0d0*rm*p33
            C0(5, 5, j) = rho*cp*w - pe*w*rho - 2.0d0*rkz

            D0(1, 1, j) = divu
            D0(1, 2, j) = rhox
            D0(1, 3, j) = rhoy + rho*d132
            D0(1, 4, j) = rhoz
            D0(2, 1, j) = pe*tx + u*(ux + v*d12) + v*uy + w*uz
            D0(2, 2, j) = rho*(ux + v*d12) + d132*rm*d12 + rmy*d12
            D0(2, 3, j) = rho*uy + rho*u*d12 - z2d3*rmx*((s - 1.0d0)*d132 + 3.0d0*d12)
            D0(2, 4, j) = rho*uz
            D0(2, 5, j) = pe*rhox - (p11x + p12y + p13z + d132*p12)*rmt - rmtx*p11 - rmty*p12 - rmtz*p13
            D0(3, 1, j) = pe*Ty + u*(vx - u*d12) + v*vy + w*(vz - w*d32)
            D0(3, 2, j) = rho*(vx - 2.0d0*u*d12) + rmx*d12
            D0(3, 3, j) = rho*vy + 2.0d0*rm*d12*d12 + 2.0d0*rm*d32*d32 + Z2D3*rmy*d132
            D0(3, 4, j) = rho*(vz - 2.0d0*w*d32) + rmz*d32
            D0(3, 5, j) = pe*rhoy - rmtx*p12 - rmty*p22 - rmtz*p23 - rmt*(p12x + p22y + p23z + p22*d132)
            D0(4, 1, j) = Tz*pe + u*wx + v*wy + w*(wz + v*d32)
            D0(4, 2, j) = rho*wx
            D0(4, 3, j) = rho*(wy + w*d32) - z2d3*((s - 1)*d132 + 3*d32)*rmz
            D0(4, 4, j) = rho*(wz + d32*v) + rm*(d12 + 2.0d0*d32)*d32 + rmy*d32
            D0(4, 5, j) = pe*rhoz - (p23*d132 + p13x + p23y + p33z)*rmt - rmtx*p13 - rmty*p23 - rmtz*p33
            D0(5, 1, j) = cp*(u*tx + v*ty + w*tz) - pe*(u*tx + v*ty + w*tz)
            D0(5, 2, j) = rho*cp*tx - pe*(rho*tx + t*rhox) + 2.0d0*rm*p12*d12
            D0(5, 3, j) = rho*cp*Ty - pe*(rho*Ty + T*rhoy) - 2.0d0*rm*(p11*d12 + p33*d32)
            D0(5, 4, j) = Cp*rho*Tz - pe*rho*Tz - pe*T*rhoz + 2.0d0*rm*p23*d32
            D0(5, 5, j) = -pe*(u*rhox + v*rhoy + w*rhoz) - (Txx + Tyy + Tzz + Ty*d132)*rkt - rktt*(Tx*Tx + Ty*Ty + Tz*Tz) - phi*rmt

            V0xx(2, 2, j) = (2.0 + s)*Z2D3*rm
            V0xx(3, 3, j) = rm
            V0xx(4, 4, j) = rm
            V0xx(5, 5, j) = rk

            V0yy(2, 2, j) = rm
            V0yy(3, 3, j) = (2.0 + s)*Z2D3*rm
            V0yy(4, 4, j) = rm
            V0yy(5, 5, j) = rk

            V0zz(2, 2, j) = rm
            V0zz(3, 3, j) = rm
            V0zz(4, 4, j) = (2.0 + s)*Z2D3*rm
            V0zz(5, 5, j) = rk

            V0xy(2, 3, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0xy(3, 2, j) = D3*(2.0d0*s + 1.0d0)*rm

            V0xz(2, 4, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0xz(4, 2, j) = D3*(2.0d0*s + 1.0d0)*rm

            V0yz(3, 4, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0yz(4, 3, j) = D3*(2.0d0*s + 1.0d0)*rm

            ! call get_Aplus(A0(:, :, j), G0(:, :, j), Aplus(:, :, j), j)
            ! Aplus(:, :, j)=Aplus(:, :, j)/hx !Aplus(:, :, j)/hx
            A0(:, :, j) = A0(:, :, j)/hx
            C0(:, :, j) = C0(:, :, j)/hz
            V0xx(:, :, j) = V0xx(:, :, j)/hx**2
            V0xy(:, :, j) = V0xy(:, :, j)/hx
            V0xz(:, :, j) = V0xz(:, :, j)/hx/hz
            V0yz(:, :, j) = V0yz(:, :, j)/hz
            V0zz(:, :, j) = V0zz(:, :, j)/hz**2
            modify_Px_A(:, j) = modify_Px_A(:, j)/hx

         end do

      end associate

   end subroutine solve_bf_matrix_j

   ! subroutine get_Aplus(A, G, Aplus, j)

   !    implicit none
   !    real(REAL64), dimension(5, 5), intent(in) :: A, G
   !    real(REAL64), dimension(5, 5), intent(out) :: Aplus
   !    integer :: j
   !    integer :: ipiv(5), info
   !    real(REAL64), dimension(5, 5) :: At, Gt
   !    real(REAL64) :: work(100)
   !    real(REAL64) :: alfr(5), alfi(5), beta(5)
   !    real(REAL64) :: vl(5, 5), vr(5, 5)
   !    complex(REAL64) :: lambda_(5)
   !    real(REAL64) :: diag(5, 5)
   !    integer :: l1, l2
   !    integer :: i
   !    real(REAL64) :: maxA

   !    Aplus = 0.0d0; diag(5, 5) = 0.0d0
   !    work = 0.0d0
   !    alfr = 0.0d0; alfi = 0.0d0; beta = 0.0d0
   !    vl = 0.0d0; vr = 0.0d0; lambda_ = 0.0d0; diag = 0.0d0
   !    At = A; Gt = G

   !    call dggev('N', 'V', 5, At, 5, Gt, 5, alfr, alfi, &
   !               beta, vl, 5, vr, 5, work, 100, info)
   !    lambda_ = (alfr + ZI*alfi)/beta
   !    do i = 1, 5
   !       diag(i, i) = max(0.0d0, lambda_(i)%re)
   !    end do
   !    vl = vr
   !      !! vl is inv(vr)
   !    call dgetrf(5, 5, vl, 5, ipiv, info)
   !    call dgetri(5, vl, 5, ipiv, work, 5, info)
   !    At = matmul(G, vr)
   !    At = matmul(At, diag)
   !    At = matmul(At, vl)

   !    Aplus = At

   ! end subroutine get_Aplus

   subroutine solve_lpse_matrix_j0(this, jb, alf)

      implicit none
      class(type_lpse_eqn) :: this
      integer, intent(in) :: jb
      complex(REAL64), intent(in) :: alf

      real(REAL64) :: rho, rhox, rhoy
      real(REAL64) :: u, ux
      real(REAL64) :: v, vy
      real(REAL64) :: w
      real(REAL64) :: divu
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      associate (grid_profile => this%grid, bf_profile => this%baseflow, ome => this%ome, bat => this%bat)
         y = grid_profile%y(1)
         drx = grid_profile%drx
         drz = grid_profile%drz

         select case (this%Coordinate_System)
         case (COMOVING_CS)
            hx = 1.d0 + drx*y
            hx_y = drx
            d12 = hx_y/hx

            hz = 1.d0 + drz*y
            hz_y = drz
            d32 = hz_y/hz
         case (CYLINDRICAL_CS)
            hx = 1.d0
            hx_y = 0.0d0
            d12 = 0.0d0

            hz = y
            if (abs(hz) <= epsilon(1.0d0)) hz = epsilon(1.0d0)
            hz_y = 1.0d0
            d32 = hz_y/hz
         end select

         d132 = d12 + d32

         rho = bf_profile%q0(1, 1)
         u = bf_profile%q0(2, 1)
         v = bf_profile%q0(3, 1)
         w = bf_profile%q0(4, 1)
         rhox = bf_profile%dq0x(1, 1)/hx
         ux = bf_profile%dq0x(2, 1)/hx
         rhoy = bf_profile%dq0y(1, 1)
         vy = bf_profile%dq0y(3, 1)

         divu = ux + vy + v*d132

         associate (D_2 => this%D_2, B_0 => this%B_0, D_0 => this%D_0, V_0 => this%V_0)
            select case (this%bc_y0)
            case (BC_WALL)
               ! V0xx(2:5, :, 1)=0.0d0
               ! V0xy(2:5, :, 1)=0.0d0
               ! V0xz(2:5, :, 1)=0.0d0
               ! V0yy(2:5, :, 1)=0.0d0
               ! V0yz(2:5, :, 1)=0.0d0
               ! V0zz(2:5, :, 1)=0.0d0
               ! modify_Px_A(:, 1)=0.0d0
               ! modify_Px_D(:, 1)=0.d00
               ! Aplus(:, :, 1)=0.0d0
               ! Aplus(1, :, 1)=A0(1, :, 1)
               !
               ! G0(1, 1, 1)=1.0d0
               ! A0(1, 1, 1)=u
               ! B0(1, 1, 1)=v
               ! B0(1, 3, 1)=rho
               ! C0(1, 1, 1)=w
               ! C0(1, 4, 1)=rho
               ! D0(1, 1, 1)= divu
               ! D0(1, 2, 1)= rhox
               ! D0(1, 3, 1)= rhoy+rho*d132
               ! D0(1, 4, 1)=0.0d0
               !
               ! D0(2, 2, 1)=1.0d0
               ! D0(3, 3, 1)=1.0d0
               ! D0(4, 4, 1)=1.0d0
               ! D0(5, 5, 1)=1.0d0
               V0xx(:, :, 1) = 0.0d0
               V0xy(:, :, 1) = 0.0d0
               V0xz(:, :, 1) = 0.0d0
               V0yy(:, :, 1) = 0.0d0
               V0yz(:, :, 1) = 0.0d0
               V0zz(:, :, 1) = 0.0d0

               D_2(2:5, :, 1) = 0.0d0
               B_0(2:5, :, 1) = 0.0d0
               V_0(2:5, :, 1) = 0.0d0
               D_0(2:5, :, 1) = 0.0d0
               D_0(2, 2, 1) = 1.0d0
               D_0(3, 3, 1) = 1.0d0
               D_0(4, 4, 1) = 1.0d0
               if (this%isothermal) then
                  D_0(5, 5, 1) = 1.0d0
               else
                  B_0(5, 5, 1) = 1.0d0
               end if
            case (BC_POLE)
               V0xx(:, :, 1) = 0.0d0
               V0xy(:, :, 1) = 0.0d0
               V0xz(:, :, 1) = 0.0d0
               V0yy(:, :, 1) = 0.0d0
               V0yz(:, :, 1) = 0.0d0
               V0zz(:, :, 1) = 0.0d0
               D_2(:, :, 1) = 0.0d0
               B_0(:, :, 1) = 0.0d0
               V_0(:, :, 1) = 0.0d0
               D_0(:, :, 1) = 0.0d0
               if (abs(bat) <= epsilon(1.0d0)) then
                  B_0(1, 1, 1) = 1.0d0
                  B_0(2, 2, 1) = 1.0d0
                  B_0(5, 5, 1) = 1.0d0

                  D_0(3, 3, 1) = 1.0d0
                  D_0(4, 4, 1) = 1.0d0
               else
                  D_0(1, 1, 1) = 1.0d0
                  D_0(2, 2, 1) = 1.0d0
                  D_0(5, 5, 1) = 1.0d0

                  B_0(3, 3, 1) = 1.0d0
                  B_0(4, 4, 1) = 1.0d0
               end if
            case (BC_POROUS)
               block
                  use Complex_Bessel
                  use mod_porous_parameters

                  complex(REAL64) :: porous_A, porous_B

!!! d, r, d is parameters

                  associate (pe => 1.0d0/(this%GAMMA*this%MA**2), &
                     &       Rho0 => this%baseflow%q0(1, 1), &
                     &       T0 => this%baseflow%q0(5, 1))
                     ! porous_A = get_porous_A

                     call get_porous_AB(Rho0, T0, ome, porous_A, porous_B)
                     ! print*, porous_A, porous_B
                     ! print*, rho0, T0, pe
                     ! print*, 'bc_porous bc is selected'
                     V0xx(:, :, 1) = 0.0d0
                     V0xy(:, :, 1) = 0.0d0
                     V0xz(:, :, 1) = 0.0d0
                     V0yy(:, :, 1) = 0.0d0
                     V0yz(:, :, 1) = 0.0d0
                     V0zz(:, :, 1) = 0.0d0

                     D_2(2:5, :, 1) = 0.0d0
                     B_0(2:5, :, 1) = 0.0d0
                     V_0(2:5, :, 1) = 0.0d0
                     D_0(2:5, :, 1) = 0.0d0

                     D_0(2, 2, 1) = 1.0d0

                     D_0(3, 1, 1) = -porous_A*T0*pe
                     D_0(3, 3, 1) = 1.0d0
                     D_0(3, 5, 1) = -porous_A*Rho0*pe

                     D_0(4, 4, 1) = 1.0d0
                     D_0(5, 1, 1) = -porous_B*T0*pe
                     D_0(5, 5, 1) = 1.0d0 - porous_B*Rho0*pe

                  end associate
               end block
            end select

         end associate

      end associate

      ! A0(:, :, 1)=A0(:, :, 1)/hx
      ! C0(:, :, 1)=C0(:, :, 1)/hz

      ! call solve_lpse_matrix_j(bf_profile,grid_profile, 1, ome, alf, bat)

   end subroutine solve_lpse_matrix_j0

   subroutine solve_lpse_matrix_jn(this, je, alf)

      implicit none
      class(type_lpse_eqn) :: this
      integer, intent(in) :: je
      complex(REAL64), intent(in) :: alf

      real(REAL64) :: rho, u, w, rhoy
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      real(REAL64) :: alfr(5), alfi(5), beta(5)
      real(REAL64) :: vl(5, 5), vr(5, 5)
      real(REAL64) :: lambda_(5), abs_lambda_(5), diag_plus(5, 5), matrix_plus(5, 5)
      integer :: ipiv(5), info
      real(REAL64) :: work(100)
      integer :: l

      associate (grid_profile => this%grid, bf_profile => this%baseflow, ome => this%ome, bat => this%bat, &
         &       D_2 => this%D_2, B_0 => this%B_0, D_0 => this%D_0, V_0 => this%V_0, jn => this%jn)
         select case (this%bc_yn)
         case (BC_FARFIELD)

            V0xx(:, :, jn) = 0.0d0
            V0xy(:, :, jn) = 0.0d0
            V0xz(:, :, jn) = 0.0d0
            V0yy(:, :, jn) = 0.0d0
            V0yz(:, :, jn) = 0.0d0
            V0zz(:, :, jn) = 0.0d0
            D_2(:, :, jn) = 0.0d0
            B_0(:, :, jn) = 0.0d0
            V_0(:, :, jn) = 0.0d0
            D_0(:, :, jn) = 0.0d0
            D_0(1, 1, jn) = 1.0d0
            D_0(2, 2, jn) = 1.0d0
            if (abs(ome) <= epsilon(1.0d0) .and. abs(bat) <= epsilon(1.0d0)) then
               ! print*, 'base flow distoration with dv/dy=0 farfield.'
               ! print*, ome, bat
               B_0(3, 3, jn) = 1.0d0
            else
               D_0(3, 3, jn) = 1.0d0
            end if
            D_0(4, 4, jn) = 1.0d0
            D_0(5, 5, jn) = 1.0d0
         case (BC_NONREFLECTION)
         !!! G^(-1)B = S^(-1)\Lambda S
            call dggev('N', 'V', 5, real(B0(:, :, jn), kind=real64), 5, &
            &          real(G0(:, :, jn), kind=real64), 5, &
            &          alfr, alfi, beta, vl, 5, vr, 5, work, 100, info)

            lambda_ = (alfr)/beta
            ! print*, lambda_
            !!! G^(-1)B^(+) = S^(-1)\Lambda^(+) S
            diag_plus = 0.0d0
            do l = 1, 5
               abs_lambda_(l) = real(sqrt(lambda_(l)**2), kind=REAL64)
               diag_plus(l, l) = 0.5d0*(lambda_(l) + abs_lambda_(l))
               ! diag_plus(l, l) = lambda_(l)
            end do
            vl = vr !! vl is inv(vr)
            call dgetrf(5, 5, vl, 5, ipiv, info)
            call dgetri(5, vl, 5, ipiv, work, 5, info)
            !!!  B^(+) = G S^(-1)\Lambda^(+) S
            B0(:, :, jn) = matmul(G0(:, :, jn), vr)
            B0(:, :, jn) = matmul(B0(:, :, jn), diag_plus)
            B0(:, :, jn) = matmul(B0(:, :, jn), vl)
         end select

      end associate

   end subroutine solve_lpse_matrix_jn
   !
   subroutine solve_lpse_matrix(this, alf)
      implicit none
      class(type_lpse_eqn) :: this
      complex(REAL64), intent(in) :: alf

      integer :: j

      do j = this%jb, this%je
         call this%solve_lpse_matrix_j(j, alf)
      end do
      this%solved_lpse_matrix = .true.

   end subroutine solve_lpse_matrix

   subroutine solve_lpse_matrix_j(this, j, alf)

      implicit none
      class(type_lpse_eqn) :: this
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: alf
      complex(REAL64) :: dadx

      complex(REAL64), dimension(5, 5) :: A_, B_, V_

        !!d\alpha/dx
!        dadx=(alf-wave_param(2)%alf)*grid_profile%xi_x(jb)
      associate (grid_profile => this%grid, bf_profile => this%baseflow, ome => this%ome, bat => this%bat)
         associate (A_0 => this%D_2, B_0 => this%B_0, D_0 => this%D_0, V_0 => this%V_0, G_0 => this%G_0, &
                    xi_x => grid_profile%xi_x(j), &
                    xi_y => grid_profile%xi_y(j), &
                    eta_x => grid_profile%eta_x(j), &
                    eta_y => grid_profile%eta_y(j), &
                    xi_yy => grid_profile%xi_yy(j), &
                    eta_yy => grid_profile%eta_yy(j), &
                    A0_rec => this%A0_rec)

            D_0(:, :, j) = D0(:, :, j) - ZI*ome*G0(:, :, j) + ZI*alf*A0(:, :, j) &
                           + ZI*bat*C0(:, :, j) + alf**2*V0xx(:, :, j) &
                           + alf*bat*V0xz(:, :, j) + bat**2*V0zz(:, :, j)! &
            !-ZI*dadx*V0xx(:, :, j)
            G_0(:, :, j) = G0(:, :, j)
!          D_0(2, 1, j)=D_0(2, 1, j)+modify_Px_D(1, j)
!          D_0(2, 5, j)=D_0(2, 5, j)+modify_Px_D(2, j)
            A_0(:, :, j) = A0(:, :, j) - ZI*bat*V0xz(:, :, j) - 2.0d0*ZI*alf*V0xx(:, :, j)
            A0_rec(:, :, j) = A_0(:, :, j)
            ! A_0(:, :, j)=A0(:, :, j)-ZI*bat*V0xz(:, :, j)!-2.0d0*ZI*alf*V0xx(:, :, j)
            A_0(2, 1, j) = A_0(2, 1, j) + modify_Px_A(1, j)
            A_0(2, 5, j) = A_0(2, 5, j) + modify_Px_A(2, j)

            ! block
            !   integer :: l1, l2
            !   print*, j
            !   do l1=1,5
            !     do l2=1,5
            !       write(*, *)l1, l2, A_0(l1, l2, j)
            !     enddo
            !   enddo
            !   read(*, *)
            !
            ! endblock
            ! print*, j, A0(2,1,j), A0(2,1, j)+modify_px_A(1, j), A_0(2, 1, j)
            ! print*, j, A0(2,5,j), A0(2,5, j)+modify_px_A(2, j), A_0(2, 5, j)
            ! read(*, *)

            B_0(:, :, j) = B0(:, :, j) - ZI*alf*V0xy(:, :, j) - ZI*bat*V0yz(:, :, j)
            V_0(:, :, j) = -V0yy(:, :, j)

            ! block
            !   integer :: l1, l2
            !   if(j==170)then
            !     print*, 'matrixG'
            !     do l1=1, 5
            !       do l2=1, 5
            !         print*, l1, l2, G0(l1, l2, j)
            !       enddo
            !     enddo
            !
            !     print*, 'matrixA'
            !     do l1=1, 5
            !       do l2=1, 5
            !         print*, l1, l2, A0(l1, l2, j)
            !       enddo
            !     enddo
            !
            !     print*, 'matrixB'
            !     do l1=1, 5
            !       do l2=1, 5
            !         print*, l1, l2, B0(l1, l2, j)
            !       enddo
            !     enddo
            !
            !     print*, 'matrixC'
            !     do l1=1, 5
            !       do l2=1, 5
            !         print*, l1, l2, C0(l1, l2, j)
            !       enddo
            !     enddo
            !
            !     print*, 'matrixD'
            !     do l1=1, 5
            !       do l2=1, 5
            !         print*, l1, l2, D0(l1, l2, j)
            !       enddo
            !     enddo
            !
            !
            !   endif
            !
            ! endblock

            A_ = A_0(:, :, j)*xi_x + B_0(:, :, j)*xi_y + V_0(:, :, j)*xi_yy
            B_ = A_0(:, :, j)*eta_x + B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy
            V_ = V_0(:, :, j)*eta_y**2

            A_0(:, :, j) = A_
            B_0(:, :, j) = B_
            V_0(:, :, j) = V_

         end associate
      end associate

   end subroutine solve_lpse_matrix_j
   !
   ! subroutine solve_lpse_matrix_temporal_j(j, ome, alf, bat)
   !
   !   implicit none
   !   integer, intent(in) :: j
   !   complex(REAL64), intent(in) :: ome
   !   complex(REAL64), intent(in) :: alf
   !   complex(REAL64), intent(in) :: bat
   !
   !   D_0(:, :, j)=D_0(:, :, j)+alf*D_1(:, :, j)+alf**2*D_2(:, :, j)
   !   B_0(:, :, j)=B_0(:, :, j)+alf*B_1(:, :, j)
   !   D_1(:, :, j)=G_0(:, :, j)
   !
   ! end subroutine solve_lpse_matrix_temporal_j
   !
   ! subroutine solve_lpse_matrix_spatial_j(j, ome, alf, bat)
   !
   !   implicit none
   !   integer, intent(in) :: j
   !   complex(REAL64), intent(in) :: ome
   !   complex(REAL64), intent(in) :: alf
   !   complex(REAL64), intent(in) :: bat
   !
   !   D_0(:, :, j)= D_0(:, :, j)+ome*G_0(:, :, j)
   !
   ! end subroutine solve_lpse_matrix_spatial_j

   subroutine solve_lpse_bc(this)

      implicit none
      class(type_lpse_eqn), intent(in) :: this

      integer :: j

      do j = this%jb, this%je
         if (j == 1 .or. j == this%jn) then
            call this%solve_lpse_bc_j(j)
         end if
      end do

   end subroutine solve_lpse_bc
   !

   subroutine solve_left_matrix(this)

      implicit none
      class(type_lpse_eqn), intent(in) :: this
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      PetscScalar :: box1(5, -2:2, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, j, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      integer :: jb_local, je_local

      associate (A_0 => this%D_2, B_0 => this%B_0, D_0 => this%D_0, V_0 => this%V_0, jn => this%jn, &
            &    C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
         jb_local = this%jb !max(2, jb)
         je_local = this%je !min(je, jn-1)

         do j = jb_local, je_local
            idxm = 0; idxn = 0
            if (j == 1) then
               ljb = 0; lje = 2
               idxn(3, -2:2) = [0, 0, 1, 2, 3] - 1
               c_idx = -2
            elseif (j == 2) then
               ljb = -1; lje = 2
               idxn(3, -2:2) = [0, 1, 2, 3, 4] - 1
               c_idx = -1
            elseif (j == jn - 1) then
               ljb = -2; lje = 1
               idxn(3, -2:2) = [jn - 3, jn - 2, jn - 1, jn, 0] - 1
               c_idx = 1
            elseif (j == jn) then
               ljb = -2; lje = 0
               idxn(3, -2:2) = [jn - 2, jn - 1, jn, 0, 0] - 1
               c_idx = 2
            else
               ljb = -2; lje = 2
               idxn(3, -2:2) = [j - 2, j - 1, j, j + 1, j + 2] - 1
               c_idx = 0
            end if
            idxm(3, 1) = j - 1

            tmp_mat_j0 = 0.0d0
            box1 = 0.0d0

            do lj = ljb, lje
               tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 + this%lambda)*A_0(:, :, j) + D_0(:, :, j)) + &
               ! tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 - 1.0d0)*A_0(:, :, j) + D_0(:, :, j)) + &
                                      C1(lj, c_idx)*(B_0(:, :, j)) + &
                                      D1(lj, c_idx)*(V_0(:, :, j))

               box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
            end do

            call MatSetValuesBlockedStencil(this%L_LPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
         end do
      end associate

      call MatAssemblyBegin(this%L_LPSE, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(this%L_LPSE, MAT_FINAL_ASSEMBLY, ierr)

      ! if(screen_flg)then
      !   call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LPSE left matrix done!\n", ierr)
      ! endif

   end subroutine solve_left_matrix

   subroutine solve_lpse_bc_j(this, j)

      use functional
      implicit none
      class(type_lpse_eqn), intent(in) :: this
      integer, intent(in) :: j
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2)
      complex(REAL64) :: tmpj0(5), tmpjn(5)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      associate (A_0 => this%D_2, B_0 => this%B_0, D_0 => this%D_0, V_0 => this%V_0, jn => this%jn, &
            &    C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
         if (.not. (j == 1 .or. j == jn)) then
            print *, 'j=', j, ' is wrong in lpse bc setting...'
            stop
         end if

         idxm = 0; idxn = 0
         if (j == 1) then
            ljb = 0; lje = 2
            idxn(3, -2:2) = [0, 0, 1, 2, 3] - 1
            c_idx = -2
         elseif (j == jn) then
            ljb = -2; lje = 0
            idxn(3, -2:2) = [jn - 2, jn - 1, jn, 0, 0] - 1
            c_idx = 2
         end if
         idxm(3, 1) = j - 1

         tmp_mat_j0 = 0.0d0
         box1 = 0.0d0

         do lj = ljb, lje
            tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 + this%lambda)*A_0(:, :, j) + D_0(:, :, j)) + &
            ! tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 -1.0d0)*A_0(:, :, j) + D_0(:, :, j)) + &
                                   C1(lj, c_idx)*(B_0(:, :, j)) + &
                                   D1(lj, c_idx)*(V_0(:, :, j))
            box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
         end do
         ! block
         !   integer :: l1, l2
         
         !   if(j==jn)then
         !   print*, j
         !   print*, ljb, lje
         !   do l1=3, 3
         !     do l2=1, 5
         !       print*, l1, l2, tmp_mat_j0(l1, l2, 0)
         !     enddo
         !   enddo
         !   stop
         !   endif
         
         ! endblock

         call MatSetValuesBlockedStencil(this%L_LPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)

            !!j==1
         if (j == 1) then
            tmpj0 = 0.0d0
            select case (this%bc_y0)
            case (BC_WALL, BC_POROUS)
               call VecSetValues(this%RHS, 4, arange(1, 4), tmpj0(1:4), INSERT_VALUES, ierr)
            case (BC_POLE)
               call VecSetValues(this%RHS, 5, arange(0, 4), tmpj0(1:5), INSERT_VALUES, ierr)
            case default
               call VecSetValues(this%RHS, 5, arange(0, 4), tmpj0(1:5), INSERT_VALUES, ierr)
            end select
         elseif (j == this%jn) then
            !!j==jn
            tmpjn = 0.0d0
            call VecSetValues(this%RHS, 5, arange((jn - 1)*5 + 0, (jn - 1)*5 + 4), tmpjn, INSERT_VALUES, ierr)
         end if

      end associate

      call MatAssemblyBegin(this%L_LPSE, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(this%L_LPSE, MAT_FINAL_ASSEMBLY, ierr)
      call VecAssemblyBegin(this%RHS, ierr)
      call VecAssemblyEnd(this%RHS, ierr)

      ! if(screen_flg)then
      !   call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LPSE matrix and RHS BC done!\n", ierr)
      ! endif
      ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)

      ! block
      !    Vec :: tmp
      !    REAL(8) :: norm

      !    if (j == this%jn) then
      !       call VecDuplicate(this%phi, tmp, ierr)
      !       call MatMult(this%L_LPSE, this%phi, tmp, ierr)
      !       ! call VecView(this%phi, PETSC_VIEWER_STDOUT_SELF, ierr)
      !       ! call VecView(this%RHS, PETSC_VIEWER_STDOUT_SELF, ierr)
      !       ! call VecAXPY(tmp, (-1.0d0, 0.0d0), this%RHS, ierr)
      !       ! call VecSetValue(tmp , 2, (0.0d0, 0.0d0), INSERT_VALUES, ierr)
      !       ! call VecView(tmp, PETSC_VIEWER_STDOUT_SELF, ierr)
      !       call VecNorm(tmp, NORM_2, norm, ierr)
      !       print*, norm
      !       stop
      !    endif
      ! endblock

   end subroutine solve_lpse_bc_j

   subroutine lpse_eqn_finalize(this)

      implicit none
      class(type_lpse_eqn), intent(inout) :: this

      deallocate (G0)
      deallocate (A0)
      deallocate (B0)
      deallocate (C0)
      deallocate (D0)
      deallocate (V0xx)
      deallocate (V0xy)
      deallocate (V0xz)
      deallocate (V0yy)
      deallocate (V0yz)
      deallocate (V0zz)
      deallocate (Aplus)

      deallocate (this%B_0)
      deallocate (this%D_0)
      deallocate (this%V_0)
      deallocate (this%D_2)
      deallocate (this%G_0)
      deallocate (this%A0_rec)

      deallocate (modify_Px_A)

      ! deallocate(phi_m1)
      ! deallocate(phi_m2)
      lpse_eqn_has_initialize = .FALSE.

      deallocate (this%phi0_array)
      deallocate (this%phi1_array)
      deallocate (this%phi2_array)

      ! call VecDestroy(phi_new, ierr)
      ! call VecDestroy(phi_old, ierr)
      call DMRestoreGlobalVector(this%dm, this%RHS, ierr)
      call DMRestoreGlobalVector(this%dm, this%phi, ierr)
      call MatDestroy(this%L_LPSE, ierr)
      call KSPDestroy(this%solver, ierr)
      call DMDestroy(this%dm, ierr)

   end subroutine lpse_eqn_finalize

end module mod_lpse_eqn
