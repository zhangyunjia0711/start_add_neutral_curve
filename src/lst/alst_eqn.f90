!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  module: module_alst_eqn
!> @breif alst equations.
!  DESCRIPTION:
!>
!!
!  REVISION HISTORY:
!  2021-01-06 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2021-01-06
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_alst_eqn

   use petsc
   use iso_fortran_env, only: REAL64
   use mod_lst_solver_abstract

   use global_parameters, only: rank

   use mod_lst_grid_profile
   use mod_lst_baseflow_profile
   use mod_slepc_adapter, only: slepceps_type, slepcpep_type, slepc_type
   implicit none

#include "../include/start.f90"

   logical :: alst_eqn_has_initialize = .FALSE.
   integer :: iloc
   integer :: jb, je
   integer :: yb, ye
   integer :: jn

   DM :: DM_ALST
   Mat, allocatable :: L_ALST(:)
   class(slepc_type), allocatable :: Solver

   complex(REAL64), allocatable :: eig_d(:)

   complex(REAL64) :: ome1, alf1

   !> \f
    !! \phi=\phi(x, y, z, t)
    !!   G0 \frac{\partial \phi}{\partial t}+ A0 \frac{\partial \phi}{h_1\partial x}
    !! + B0 \frac{\partial \phi}{\partial y}+ C0 \frac{\partial \phi}{h_3\partial z}
    !! + D0 phi=
    !!  V_{0xx} \frac{\partial^2 \phi}{h_1^2 \partial x^2}
    !! +V_{0xy} \frac{\partial^2 \phi}{h_1   \partial x\partial y}
    !! +V_{0xz} \frac{\partial^2 \phi}{h_1h_3\partial x\partial z}
    !! +V_{0yy} \frac{\partial^2 \phi}{      \partial y^2}
    !! +V_{0yz} \frac{\partial^2 \phi}{h_3   \partial y\partial z}
    !! +V_{0zz} \frac{\partial^2 \phi}{h_3^2 \partial z^2}
    !! \f

   real(REAL64), allocatable, dimension(:, :, :) :: G0, A0, B0, C0, D0
   real(REAL64), allocatable, dimension(:, :, :) :: V0xx, V0xy, V0xz
   real(REAL64), allocatable, dimension(:, :, :) :: V0yy, V0yz, V0zz

   !> \f
    !! \phi=\hat\phi(y)
    !!  L\hat\phi=0
    !!  L=
    !!  (+D_0 + B_0\frac{\partial}{\parital y}+ V_0\frac{\partial^2}{\parital y^2})
    !!  +\omega*(G_0)
    !!  +alpha*(D_1+B_1\frac{\partial}{\parital y})
    !!  +alpha^2*(D_2)
    !! \f
   complex(REAL64), allocatable, dimension(:, :, :) :: G_0, B_0, D_0, V_0
   complex(REAL64), allocatable, dimension(:, :, :) ::      B_1, D_1
   complex(REAL64), allocatable, dimension(:, :, :) ::           D_2

   real(REAL64), allocatable, dimension(:, :, :) :: dB0, ddV0yy, dV0xy
   real(REAL64), allocatable, dimension(:, :, :) :: dV0yy, dV0yz
   complex(REAL64), allocatable, dimension(:, :, :) :: dB_0, dV_0, ddV_0
   complex(REAL64), allocatable, dimension(:, :, :) :: dB_1

   PetscErrorCode :: ierr

   type, extends(type_lst_abstract) :: type_alst_eqn

   contains
      procedure :: initial => alst_intial
      procedure :: setup => alst_setup
      procedure :: solve => alst_solve
      procedure :: get_eigf => alst_get_eigenfunction
      procedure :: finalize => alst_finalize

      procedure, private :: solve_baseflow_matrix
      procedure, private :: solve_bf_matrix_j

      procedure, private :: solve_alst_matrix
      procedure, private :: solve_alst_matrix_j
      procedure, private :: solve_alst_matrix_bc
      procedure, private :: solve_alst_matrix_j0
      procedure, private :: solve_alst_matrix_jn
      procedure, private :: solve_EVP_local
      procedure, private :: solve_EVP_matrix
      procedure, private :: solve_EVP_matrix_bc

   end type type_alst_eqn

   private
   ! public :: solve_alst_eqn
   ! public :: write_alst_function
   ! public :: initialize_alst_eqn
   ! public :: finalize_alst_eqn

   public :: type_alst_eqn

contains

   subroutine alst_intial(this, jn_, Re, Ma, Te, yb_, ye_)

      implicit none
      class(type_alst_eqn), intent(inout) :: this
      integer, intent(in) :: jn_
      real(real64), intent(in) ::  Re, Ma, Te
      integer, intent(in), optional :: yb_, ye_
      integer :: xb1, xl1

      if (.not. present(yb_)) then
         yb = 0
      else
         yb = yb_
      end if
      if (.not. present(ye_)) then
         ye = jn - 1
      else
         ye = ye_
      end if

      jn = jn_
      jb = yb + 1; je = ye + 1
      this%Re = Re
      this%Ma = Ma
      this%Te = Te

      allocate (G0(5, 5, jb:je))
      allocate (A0(5, 5, jb:je))
      allocate (B0(5, 5, jb:je))
      allocate (C0(5, 5, jb:je))
      allocate (D0(5, 5, jb:je))
      allocate (V0xx(5, 5, jb:je))
      allocate (V0xy(5, 5, jb:je))
      allocate (V0xz(5, 5, jb:je))
      allocate (V0yy(5, 5, jb:je))
      allocate (V0yz(5, 5, jb:je))
      allocate (V0zz(5, 5, jb:je))

      allocate (G_0(5, 5, jb:je))
      allocate (B_0(5, 5, jb:je))
      allocate (D_0(5, 5, jb:je))
      allocate (V_0(5, 5, jb:je))
      allocate (B_1(5, 5, jb:je))
      allocate (D_1(5, 5, jb:je))
      allocate (D_2(5, 5, jb:je))

      allocate (dB0(5, 5, jb:je))
      allocate (ddV0yy(5, 5, jb:je))
      allocate (dV0xy(5, 5, jb:je))
      allocate (dV0yy(5, 5, jb:je))
      allocate (dV0yz(5, 5, jb:je))

      allocate (dB_0(5, 5, jb:je))
      allocate (dV_0(5, 5, jb:je))
      allocate (ddV_0(5, 5, jb:je))
      allocate (dB_1(5, 5, jb:je))

      alst_eqn_has_initialize = .TRUE.

      call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_ALST, ierr)
      call DMSetUp(DM_ALST, ierr)
      if (this%lst_mode == LST_TEMPORAL_MODE) then
         allocate (L_ALST(2))
         call DMCreateMatrix(DM_ALST, L_ALST(1), ierr)
         call DMCreateMatrix(DM_ALST, L_ALST(2), ierr)
      else if (this%lst_mode == LST_SPATIAL_MODE) then
         allocate (L_ALST(3))
         call DMCreateMatrix(DM_ALST, L_ALST(1), ierr)
         call DMCreateMatrix(DM_ALST, L_ALST(2), ierr)
         call DMCreateMatrix(DM_ALST, L_ALST(3), ierr)
      end if

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         allocate (slepceps_type :: Solver)
      case (LST_SPATIAL_MODE)
         allocate (slepcpep_type::Solver)
      end select

      ! select case (lst_method)
      ! case (LST_LOCAL, LST_DOMAIN)
      call Solver%Create(1, this%lst_mode)
      ! case (LST_GLOBAL)
      !   call Solver%Create(lst_nev)
      ! end select

   end subroutine alst_intial

   subroutine alst_setup(this)
      implicit none
      class(type_alst_eqn), intent(inout) :: this

      call this%solve_baseflow_matrix()

   end subroutine alst_setup

   subroutine alst_solve(this, ome, alf, bat)

      implicit none
      class(type_alst_eqn), intent(in) :: this
      complex(REAL64), intent(inout) :: ome, alf, bat

      ! call this%solve_lst_matrix(ome, alf, bat)

      ! call this%solve_EVP_matrix()
      ! call this%solve_lst_matrix_bc(ome, alf, bat)
      ! call this%solve_EVP_matrix_bc()
      ! call this%solve_EVP_local(ome, alf, bat)
      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         ome = conjg(ome)
      case (LST_SPATIAL_MODE)
         alf = conjg(alf)
      end select

      call this%solve_alst_matrix(ome, alf, bat)
      call this%solve_EVP_matrix()
      call this%solve_alst_matrix_bc(ome, alf, bat)
      call this%solve_EVP_matrix_bc()

      call this%solve_EVP_local(ome, alf, bat)

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         ome = conjg(ome)
      case (LST_SPATIAL_MODE)
         alf = conjg(alf)
      end select

   end subroutine alst_solve

   subroutine alst_get_eigenfunction(this, eigf)

      use mod_lst_solver_abstract, only: normal_u
      implicit none
      class(type_alst_eqn), intent(in) :: this
      Vec, intent(out) :: eigf
      PetscScalar :: d

      call Solver%GetEigenVector(d, eigf)
      call normal_u(eigf)

   end subroutine alst_get_eigenfunction

   subroutine solve_baseflow_matrix(this)

      implicit none
      class(type_alst_eqn) :: this
      integer :: j, j_b, j_e

      ! if(jb==1) then
      !   j_b=jb+1
      ! else
      !   j_b=jb
      ! endif
      !
      ! if(je==jn) then
      !   j_e=je-1
      ! else
      !   j_e=je
      ! endif

        !! internal point
      call this%solve_bf_matrix_j(jb, je)
   end subroutine solve_baseflow_matrix

   subroutine solve_alst_matrix_bc(this, ome, alf, bat)

      implicit none
      class(type_alst_eqn) :: this
      complex(REAL64), intent(in) :: ome, alf, bat
        !! b.c. point
      if (jb == 1) call this%solve_alst_matrix_j0(jb, ome, alf, bat)
      if (je == jn) call this%solve_alst_matrix_jn(je, ome, alf, bat)

   end subroutine solve_alst_matrix_bc

   subroutine solve_bf_matrix_j(this, jb, je)

      implicit none
      class(type_alst_eqn) :: this
      integer, intent(in) :: jb, je

      real(REAL64) :: u, uy, uyy, w, wy, wyy
      real(REAL64) :: rho, rhoy, T, Ty, Tyy

      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d12y, d32y, d132, d132y

      real(REAL64) :: p12, p23, p12y, p23y, phi
      real(REAL64) :: pe, cm, ck
      real(REAL64) :: rm, rmy, rmt, rmtt, rmty, rmyy
      real(REAL64) :: rk, rkt, rktt, rky, rkty, rkyy
      real(REAL64) :: cv, cp, rek

      integer :: j
      real(REAL64), parameter :: Z2D3 = 2.d0/3.d0, D3 = 1.d0/3.d0

      associate (grid_profile => this%grid, bf_profile => this%baseflow, &
                 Ma => this%Ma, Re => this%Re, Te => this%Te, Pr => This%Pr, &
                 GAMMA => this%GAMMA, C1 => this%C1, C2 => this%C2)

         do j = jb, je

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

            dB0(:, :, j) = 0.0d0
            dV0yy(:, :, j) = 0.0d0
            dV0xy(:, :, j) = 0.0d0
            dV0yz(:, :, j) = 0.0d0
            ddV0yy(:, :, j) = 0.0d0

            u = bf_profile%q0(2, j)
            uy = bf_profile%dq0y(2, j)
            uyy = bf_profile%dq0yy(2, j)
            w = bf_profile%q0(4, j)
            wy = bf_profile%dq0y(4, j)
            wyy = bf_profile%dq0yy(4, j)
            T = bf_profile%q0(5, j)
            Ty = bf_profile%dq0y(5, j)
            Tyy = bf_profile%dq0yy(5, j)
            rho = bf_profile%q0(1, j)
            rhoy = bf_profile%dq0y(1, j)

            pe = 1.d0/(GAMMA*Ma*Ma)
            cm = C1/Te
            ck = C2/Te
            cv = pe/(GAMMA - 1.d0)
            cp = cv*GAMMA
            rek = re*pr/cp

            rm = T*sqrt(T)*(1.0d0 + cm)/(T + cm)/Re
            rmt = rm*(1.5d0/T - 1.0d0/(T + cm))
            rmtt = rmt*(1.5d0/T - 1.0d0/(T + cm)) - rm*(1.5d0/T**2 - 1.0d0/(T + cm)**2)
            rmy = rmt*ty
            rmty = rmtt*ty
            rmyy = rmtt*ty**2 + rmt*tyy

            rk = T*sqrt(T)*(1.d0 + ck)/(T + ck)/rek
            rkt = rk*(1.5d0/T - 1.0d0/(T + ck))
            rktt = rkt*(1.5d0/T - 1.0d0/(T + ck)) - rk*(1.5d0/T**2 - 1.0d0/(T + ck)**2)
            rky = rkt*ty
            rkty = rktt*ty
            rkyy = rktt*ty**2 + rkt*tyy

            y = grid_profile%y(j)
            drx = grid_profile%drx
            drz = grid_profile%drz

            select case (this%Coordinate_System)
            case (COMOVING_CS)
               hx = 1.d0 + drx*y
               hx_y = drx
               d12 = hx_y/hx
               d12y = -d12**2

               hz = 1.d0 + drz*y
               hz_y = drz
               d32 = hz_y/hz
               d32y = -d32**2
            case (CYLINDRICAL_CS)
               hx = 1.d0
               hx_y = 0.0d0
               d12 = hx_y/hx
               d12y = -d12**2

               hz = y
               hz_y = 1.0d0
               d32 = hz_y/hz
               d32y = -d32**2
            end select

            d132 = d12 + d32
            d132y = d12y + d32y

            p12 = uy - d12*u
            p23 = wy - d32*w

            p12y = uyy - d12*uy + d12**2*u
            p23y = wyy - d32*wy + d32**2*w

            phi = p12**2 + p23**2

            G0(1, 1, j) = 1.0d0
            G0(2, 2, j) = rho
            G0(3, 3, j) = rho
            G0(4, 4, j) = rho
            G0(5, 1, j) = -pe*t
            G0(5, 5, j) = rho*cp - pe*rho

            A0(1, 1, j) = u
            A0(1, 2, j) = rho
            A0(2, 1, j) = pe*T
            A0(2, 2, j) = rho*u
            A0(2, 3, j) = -2.0d0*d12*rm - d3*rm*d132 - rmy
            A0(2, 5, j) = pe*rho
            A0(3, 2, j) = 7.0d0/3.0d0*rm*d12 + Z2D3*rmy
            A0(3, 3, j) = rho*u
            A0(3, 5, j) = -rmt*p12
            A0(4, 4, j) = rho*u
            A0(5, 1, j) = -pe*u*T
            A0(5, 3, j) = -2.0d0*rm*p12
            A0(5, 5, j) = rho*cp*u - pe*u*rho

            B0(1, 3, j) = rho
            B0(2, 2, j) = -d132*rm - rmy
            B0(2, 5, j) = -rmt*p12
            B0(3, 1, j) = pe*T
            B0(3, 3, j) = -2.0d0*Z2D3*rm*d132 - 2.0d0*Z2d3*rmy
            B0(3, 5, j) = pe*rho
            B0(4, 4, j) = -rm*d132 - rmy
            B0(4, 5, j) = -rmt*p23
            B0(5, 2, j) = -2.0d0*rm*p12
            B0(5, 4, j) = -2.0d0*rm*p23
            B0(5, 5, j) = -2.0d0*rky - rk*d132

            dB0(1, 3, j) = rhoy
            dB0(2, 2, j) = -d132*rmy - d132y*rm - rmyy
            dB0(2, 5, j) = -rmty*p12 - rmt*p12y
            dB0(3, 1, j) = pe*Ty
            dB0(3, 3, j) = -2.0d0*Z2D3*rmy*d132 - 2.0d0*Z2D3*rm*d132y &
                           - 2.0d0*Z2d3*rmyy
            dB0(3, 5, j) = pe*rhoy
            dB0(4, 4, j) = -rmy*d132 - rm*d132y - rmyy
            dB0(4, 5, j) = -rmty*p23 - rmt*p23y
            dB0(5, 2, j) = -2.0d0*rmy*p12 - 2.0d0*rm*p12y
            dB0(5, 4, j) = -2.0d0*rmy*p23 - 2.0d0*rm*p23y
            dB0(5, 5, j) = -2.0d0*rkyy - rky*d132 - rk*d132y

            C0(1, 1, j) = w
            C0(1, 4, j) = rho
            C0(2, 2, j) = rho*w
            C0(3, 3, j) = rho*w
            C0(3, 4, j) = 7.0d0/3.0d0*rm*d32 + Z2D3*rmy
            C0(3, 5, j) = -rmt*p23
            C0(4, 1, j) = pe*t
            C0(4, 3, j) = -D3*rm*d12 - 7.0d0/3.0d0*rm*d32 - rmy
            C0(4, 4, j) = rho*w
            C0(4, 5, j) = pe*rho
            C0(5, 1, j) = -pe*w*T
            C0(5, 3, j) = -2.0d0*rm*p23
            C0(5, 5, j) = rho*cp*w - pe*w*rho

            D0(1, 3, j) = rhoy + rho*d132
            D0(2, 2, j) = (2.0d0*d12 + d32)*rm*d12 + rmy*d12
            D0(2, 3, j) = rho*uy + rho*u*d12
            D0(2, 5, j) = -(2.0d0*d12 + d32)*rmt*p12 - rmty*p12 - rmt*p12y
            D0(3, 1, j) = -u*u*d12 - w*w*d32 + pe*Ty
            D0(3, 2, j) = -2.0d0*rho*u*d12
            D0(3, 3, j) = 2.0d0*rm*d12*d12 + 2.0d0*rm*d32*d32 + Z2D3*rmy*d132
            D0(3, 4, j) = -2.0d0*rho*w*d32
            D0(3, 5, j) = pe*rhoy
            D0(4, 3, j) = rho*wy + rho*w*d32
            D0(4, 4, j) = rm*(d12 + 2.0d0*d32)*d32 + rmy*d32
            D0(4, 5, j) = -p23*(d12 + 2.0d0*d32)*rmt - p23y*rmt - rmty*p23
            D0(5, 2, j) = +2.0d0*rm*p12*d12
            D0(5, 3, j) = rho*cp*Ty - pe*(rho*Ty + T*rhoy)
            D0(5, 4, j) = 2.0d0*rm*p23*d32
            D0(5, 5, j) = -Ty*d132*rkt - phi*rmt - Tyy*rkt - rktt*Ty*Ty!-rkty*Ty

            V0xx(2, 2, j) = 2.0*Z2D3*rm
            V0xx(3, 3, j) = rm
            V0xx(4, 4, j) = rm
            V0xx(5, 5, j) = rk

            V0yy(2, 2, j) = rm
            V0yy(3, 3, j) = 2.0d0*Z2D3*rm
            V0yy(4, 4, j) = rm
            V0yy(5, 5, j) = rk

            dV0yy(2, 2, j) = rmy
            dV0yy(3, 3, j) = 2.0d0*Z2D3*rmy
            dV0yy(4, 4, j) = rmy
            dV0yy(5, 5, j) = rky

            ddV0yy(2, 2, j) = rmyy
            ddV0yy(3, 3, j) = 2.0d0*Z2D3*rmyy
            ddV0yy(4, 4, j) = rmyy
            ddV0yy(5, 5, j) = rkyy

            V0zz(2, 2, j) = rm
            V0zz(3, 3, j) = rm
            V0zz(4, 4, j) = 2.0d0*Z2D3*rm
            V0zz(5, 5, j) = rk

            V0xy(2, 3, j) = D3*rm
            V0xy(3, 2, j) = d3*rm

            dV0xy(2, 3, j) = D3*rmy
            dV0xy(3, 2, j) = D3*rmy

            V0xz(2, 4, j) = D3*rm
            V0xz(4, 2, j) = D3*rm

            V0yz(3, 4, j) = D3*rm
            V0yz(4, 3, j) = D3*rm

            dV0yz(3, 4, j) = D3*rmy
            dV0yz(4, 3, j) = D3*rmy

            A0(:, :, j) = A0(:, :, j)/hx
            C0(:, :, j) = C0(:, :, j)/hz
            V0xx(:, :, j) = V0xx(:, :, j)/hx**2
            V0xy(:, :, j) = V0xy(:, :, j)/hx
            dV0xy(:, :, j) = dV0xy(:, :, j)/hx
            V0xz(:, :, j) = V0xz(:, :, j)/hx/hz
            V0yz(:, :, j) = V0yz(:, :, j)/hz
            dV0yz(:, :, j) = dV0yz(:, :, j)/hz
            V0zz(:, :, j) = V0zz(:, :, j)/hz**2
            !
            ! G0(:, :, j)=conjg(transpose(G0(:, :, j)))
            ! A0(:, :, j)=conjg(transpose(A0(:, :, j)))
            ! B0(:, :, j)=conjg(transpose(B0(:, :, j)))
            ! C0(:, :, j)=conjg(transpose(C0(:, :, j)))
            ! V0xx(:, :, j)=conjg(transpose(V0xx(:, :, j)))
            ! V0xy(:, :, j)=conjg(transpose(V0xy(:, :, j)))
            ! V0xz(:, :, j)=conjg(transpose(V0xz(:, :, j)))
            ! V0yy(:, :, j)=conjg(transpose(V0yy(:, :, j)))
            ! V0yz(:, :, j)=conjg(transpose(V0yz(:, :, j)))
            ! V0zz(:, :, j)=conjg(transpose(V0zz(:, :, j)))
            ! dB0(:, :, j)=conjg(transpose(dB0(:, :, j)))
            ! dV0yy(:, :, j)=conjg(transpose(dV0yy(:, :, j)))
            ! ddV0yy(:, :, j)=conjg(transpose(ddV0yy(:, :, j)))
            ! dV0xy(:, :, j)=conjg(transpose(dV0xy(:, :, j)))
            ! dV0yz(:, :, j)=conjg(transpose(dV0yz(:, :, j)))

         end do
      end associate

   end subroutine solve_bf_matrix_j

   subroutine solve_alst_matrix_j0(this, jb, ome, alf, bat)

      implicit none
      class(type_alst_eqn) :: this
      integer, intent(in) :: jb
      complex(REAL64), intent(in) :: ome, alf, bat

      real(REAL64) :: rho, u, w, rhoy
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      associate (bf_profile => this%baseflow, grid_profile => this%grid)

         rho = bf_profile%q0(1, 1)
         u = bf_profile%q0(2, 1)
         w = bf_profile%q0(3, 1)
         rhoy = bf_profile%dq0y(1, 1)

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
            d12 = hx_y/hx

            hz = y
            hz_y = 1.0d0
            d32 = hz_y/hz
         end select

         d132 = d12 + d32

         select case (this%bc_y0)
         case (BC_WALL)

            G0(:, 2:5, 1) = 0.0d0
            A0(:, 2:5, 1) = 0.0d0
            B0(:, 2:5, 1) = 0.0d0
            C0(:, 2:5, 1) = 0.0d0
            D0(:, 2:5, 1) = 0.0d0
            V0xx(:, 2:5, 1) = 0.0d0
            V0xy(:, 2:5, 1) = 0.0d0
            V0xz(:, 2:5, 1) = 0.0d0
            V0yy(:, 2:5, 1) = 0.0d0
            V0yz(:, 2:5, 1) = 0.0d0
            V0zz(:, 2:5, 1) = 0.0d0
            dB0(:, 2:5, 1) = 0.0d0
            dV0yy(:, 2:5, 1) = 0.0d0
            dV0yz(:, 2:5, 1) = 0.0d0
            dV0xy(:, 2:5, 1) = 0.0d0
            ddV0yy(:, 2:5, 1) = 0.0d0

            ! G0(1, 1, 1)=1.0d0
            ! A0(1, 1, 1)=u
            ! A0(1, 2, 1)=rho
            ! B0(1, 3, 1)=rho
            ! C0(1, 1, 1)=w
            ! D0(1, 3, 1)= rhoy+rho*d132
            ! D0(1, 1, 1)=1.0d0
            D0(2, 2, 1) = 1.0d0
            D0(3, 3, 1) = 1.0d0
            D0(4, 4, 1) = 1.0d0
            if (this%isothermal) then
               D0(5, 5, 1) = 1.0d0
            else
               B0(5, 5, 1) = 1.0d0
            end if
         case (BC_POLE)
            G0(:, :, 1) = 0.0d0
            A0(:, :, 1) = 0.0d0
            B0(:, :, 1) = 0.0d0
            C0(:, :, 1) = 0.0d0
            D0(:, :, 1) = 0.0d0
            V0xx(:, :, 1) = 0.0d0
            V0xy(:, :, 1) = 0.0d0
            V0xz(:, :, 1) = 0.0d0
            V0yy(:, :, 1) = 0.0d0
            V0yz(:, :, 1) = 0.0d0
            V0zz(:, :, 1) = 0.0d0

            if (int(bat + 0.001d0) == 0) then
               B0(1, 1, 1) = 1.0d0
               B0(2, 2, 1) = 1.0d0
               B0(5, 5, 1) = 1.0d0

               D0(3, 3, 1) = 1.0d0
               D0(4, 4, 1) = 1.0d0
            else
               D0(1, 1, 1) = 1.0d0
               D0(2, 2, 1) = 1.0d0
               D0(5, 5, 1) = 1.0d0

               B0(3, 3, 1) = 1.0d0
               B0(4, 4, 1) = 1.0d0
            end if

         end select

         A0(:, :, 1) = A0(:, :, 1)/hx
         C0(:, :, 1) = C0(:, :, 1)/hz

         call this%solve_alst_matrix_j(1, ome, alf, bat)

      end associate

   end subroutine solve_alst_matrix_j0

   subroutine solve_alst_matrix_jn(this, je, ome, alf, bat)

      implicit none
      class(type_alst_eqn) :: this
      integer, intent(in) :: je
      complex(REAL64), intent(in) :: ome, alf, bat

      real(REAL64) :: rho, u, w, rhoy
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      select case (this%bc_yn)
      case (BC_FARFIELD)
         G0(:, :, jn) = 0.0d0
         A0(:, :, jn) = 0.0d0
         B0(:, :, jn) = 0.0d0
         C0(:, :, jn) = 0.0d0
         D0(:, :, jn) = 0.0d0
         V0xx(:, :, jn) = 0.0d0
         V0xy(:, :, jn) = 0.0d0
         V0xz(:, :, jn) = 0.0d0
         V0yy(:, :, jn) = 0.0d0
         V0yz(:, :, jn) = 0.0d0
         V0zz(:, :, jn) = 0.0d0
         dB0(:, :, jn) = 0.0d0
         dV0yy(:, :, jn) = 0.0d0
         dV0yz(:, :, jn) = 0.0d0
         dV0xy(:, :, jn) = 0.0d0
         ddV0yy(:, :, jn) = 0.0d0

         D0(1, 1, jn) = 1.0d0
         D0(2, 2, jn) = 1.0d0
         D0(3, 3, jn) = 1.0d0
         D0(4, 4, jn) = 1.0d0
         D0(5, 5, jn) = 1.0d0
      end select
      call this%solve_alst_matrix_j(jn, ome, alf, bat)

   end subroutine solve_alst_matrix_jn

   subroutine solve_alst_matrix(this, ome, alf, bat)
      implicit none
      class(type_alst_eqn) :: this
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      integer :: j

      do j = jb, je
         call this%solve_alst_matrix_j(j, ome, alf, bat)
      end do

   end subroutine solve_alst_matrix

   subroutine solve_alst_matrix_j(this, j, ome, alf, bat)

      implicit none
      class(type_alst_eqn) :: this
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      G_0(:, :, j) = -ZI*G0(:, :, j)

      D_0(:, :, j) = D0(:, :, j) + ZI*bat*C0(:, :, j) + bat**2*V0zz(:, :, j)
      B_0(:, :, j) = B0(:, :, j) - ZI*bat*V0yz(:, :, j)
      V_0(:, :, j) = -V0yy(:, :, j)

      dB_0(:, :, j) = dB0(:, :, j) - ZI*bat*dV0yz(:, :, j)
      dV_0(:, :, j) = -dV0yy(:, :, j)
      ddV_0(:, :, j) = -ddV0yy(:, :, j)
      D_0(:, :, j) = D_0(:, :, j) - dB_0(:, :, j) + ddV_0(:, :, j)
      B_0(:, :, j) = -B_0(:, :, j) + 2.0d0*dV_0(:, :, j)

      D_1(:, :, j) = ZI*A0(:, :, j) + bat*V0xz(:, :, j)
      B_1(:, :, j) = -ZI*V0xy(:, :, j)

      dB_1(:, :, j) = -ZI*dV0xy(:, :, j)
      D_1(:, :, j) = D_1(:, :, j) - dB_1(:, :, j)
      B_1(:, :, j) = -B_1(:, :, j)

      D_2(:, :, j) = +V0xx(:, :, j)

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call solve_alst_matrix_temporal_j(j, ome, alf, bat)
      case (LST_SPATIAL_MODE)
         call solve_alst_matrix_spatial_j(j, ome, alf, bat)
      end select

      G_0(:, :, j) = conjg(transpose(G_0(:, :, j)))
      D_0(:, :, j) = conjg(transpose(D_0(:, :, j)))
      B_0(:, :, j) = conjg(transpose(B_0(:, :, j)))
      V_0(:, :, j) = conjg(transpose(V_0(:, :, j)))
      D_1(:, :, j) = conjg(transpose(D_1(:, :, j)))
      B_1(:, :, j) = conjg(transpose(B_1(:, :, j)))
      D_2(:, :, j) = conjg(transpose(D_2(:, :, j)))

! block
!   integer :: l1, l2
!         if(jb==1)then
! do l1=1, 5
!   do l2=1, 5
!     print("(2I0, 22F12.7)"),l1, l2, G_0(l1, l2, 1), D_0(l1, l2, 1), &
!     B_0(l1, l2, 1), V_0(l1, l2, 1), D_1(l1, l2, 1), B_1(l1, l2, 1), D_2(l1, l2, 1)
!   enddo
! enddo
!
!         endif
! endblock
   end subroutine solve_alst_matrix_j

   subroutine solve_alst_matrix_temporal_j(j, ome, alf, bat)

      implicit none
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      D_0(:, :, j) = D_0(:, :, j) + alf*D_1(:, :, j) + alf**2*D_2(:, :, j)
      B_0(:, :, j) = B_0(:, :, j) + alf*B_1(:, :, j)
      D_1(:, :, j) = G_0(:, :, j)

      ! block
      !   integer :: l1, l2
      !   if(rank==0)then
      !     do l1=1, 5
      !       do l2=1, 5
      !         write(*, '(3I0, 6(F14.6, 2X))')rank, l1, l2, D_0(l1, l2, 1), B_0(l1, l2, 1), D_1(l1, l2, 1)
      !       enddo
      !     enddo
      !   endif
      !   if(rank==1)then
      !     do l1=1, 5
      !       do l2=1, 5
      !         write(*, '(3I0, 6(F14.6, 2X))')rank, l1, l2, D_0(l1, l2, jn), B_0(l1, l2, jn), D_1(l1, l2, jn)
      !       enddo
      !     enddo
      !   endif
      ! endblock

   end subroutine solve_alst_matrix_temporal_j

   subroutine solve_alst_matrix_spatial_j(j, ome, alf, bat)

      implicit none
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      D_0(:, :, j) = D_0(:, :, j) + ome*G_0(:, :, j)

   end subroutine solve_alst_matrix_spatial_j

   subroutine solve_EVP_matrix(this)

      implicit none
      class(type_alst_eqn) :: this

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call solve_EVP_matrix_temporal(this%grid)
      case (LST_SPATIAL_MODE)
         call solve_EVP_matrix_spatial(this%grid)
      end select

   end subroutine solve_EVP_matrix

   subroutine solve_EVP_matrix_bc(this)

      implicit none
      class(type_alst_eqn) :: this

      integer :: j

      do j = jb, je
         if (j == 1 .or. j == jn) then
            select case (this%lst_mode)
            case (LST_TEMPORAL_MODE)
               call solve_EVP_matrix_bc_temporal(j, this%grid)
            case (LST_SPATIAL_MODE)
               call solve_EVP_matrix_bc_spatial(j, this%grid)
            end select
         end if
      end do

   end subroutine solve_EVP_matrix_bc

   subroutine solve_EVP_matrix_temporal(grid_profile)

      implicit none
      type(type_lst_grid_profile), intent(in) :: grid_profile
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, j, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      do j = jb, je
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

         associate (eta_y => grid_profile%eta_y(j), eta_yy => grid_profile%eta_yy(j), &
                &  C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
            tmp_mat_j0 = 0.0d0
            tmp_mat_j1 = 0.0d0
            box1 = 0.0d0
            box2 = 0.0d0
            do lj = ljb, lje
               tmp_mat_j0(:, :, lj) = delta(lj)*(D_0(:, :, j)) + &
                                      C1(lj, c_idx)*(B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy) + &
                                      D1(lj, c_idx)*(V_0(:, :, j)*eta_y*eta_y)

               tmp_mat_j1(:, :, lj) = -delta(lj)*(D_1(:, :, j))
               box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
            end do
            box2(:, :) = transpose(tmp_mat_j1(:, :, 0))
         end associate

         call MatSetValuesBlockedStencil(L_ALST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
         call MatSetValuesBlocked(L_ALST(2), 1, j - 1, 1, j - 1, box2, INSERT_VALUES, ierr)
      end do

      do i = 1, 2
         call MatAssemblyBegin(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set temporal ALST EVP matrix done!\n", ierr)

      ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'ALST_matrix1.mat', FILE_MODE_WRITE, &
      ! Viewer, ierr)
      ! call MatView(L_ALST(1), viewer, ierr)
      ! call PetscViewerDestroy(Viewer, ierr)
      ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'ALST_matrix2.mat', FILE_MODE_WRITE, &
      ! Viewer, ierr)
      ! call MatView(L_ALST(2), viewer, ierr)
      ! call PetscViewerDestroy(Viewer, ierr)
      !
      ! call SLEPCFinalize(ierr)
      !
      ! stop

   end subroutine solve_EVP_matrix_temporal

   subroutine solve_EVP_matrix_spatial(grid_profile)

      implicit none
      type(type_lst_grid_profile), intent(in) :: grid_profile
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, -2:2, 5), box3(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2), tmp_mat_j2(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, j, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      do j = jb, je
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

         associate (eta_y => grid_profile%eta_y(j), eta_yy => grid_profile%eta_yy(j), &
                &  C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
            tmp_mat_j0 = 0.0d0
            tmp_mat_j1 = 0.0d0
            box1 = 0.0d0
            box2 = 0.0d0
            do lj = ljb, lje
               tmp_mat_j0(:, :, lj) = delta(lj)*(D_0(:, :, j)) + &
                                      C1(lj, c_idx)*(B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy) + &
                                      D1(lj, c_idx)*(V_0(:, :, j)*eta_y*eta_y)

               tmp_mat_j1(:, :, lj) = delta(lj)*(D_1(:, :, j)) + &
                                      C1(lj, c_idx)*(B_1(:, :, j)*eta_y)

               tmp_mat_j2(:, :, lj) = delta(lj)*(D_2(:, :, j))

               box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
               box2(:, lj, :) = transpose(tmp_mat_j1(:, :, lj))
            end do
            box3(:, :) = transpose(tmp_mat_j2(:, :, 0))
         end associate

         call MatSetValuesBlockedStencil(L_ALST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
         call MatSetValuesBlockedStencil(L_ALST(2), 1, idxm, 5, idxn, box2, INSERT_VALUES, ierr)
         call MatSetValuesBlocked(L_ALST(3), 1, j - 1, 1, j - 1, box3, INSERT_VALUES, ierr)
      end do

      do i = 1, 3
         call MatAssemblyBegin(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set spatial ALST EVP matrix done!\n", ierr)

   end subroutine solve_EVP_matrix_spatial

   subroutine solve_EVP_matrix_bc_temporal(j, grid_profile)

      implicit none
      integer, intent(in) :: j
      type(type_lst_grid_profile), intent(in) :: grid_profile
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2)

      integer :: i, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      if (.not. (j == 1 .or. j == jn)) then
         print *, 'j=', j, ' is wrong in EVP matrix bc setting...'
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

      associate (eta_y => grid_profile%eta_y(j), eta_yy => grid_profile%eta_yy(j), &
             &  C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
         tmp_mat_j0 = 0.0d0
         tmp_mat_j1 = 0.0d0
         box1 = 0.0d0
         box2 = 0.0d0
         do lj = ljb, lje
            tmp_mat_j0(:, :, lj) = delta(lj)*(D_0(:, :, j)) + &
                                   C1(lj, c_idx)*(B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy) + &
                                   D1(lj, c_idx)*(V_0(:, :, j)*eta_y*eta_y)

            tmp_mat_j1(:, :, lj) = -delta(lj)*(D_1(:, :, j))
            box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
         end do
         box2(:, :) = transpose(tmp_mat_j1(:, :, 0))
      end associate

      call MatSetValuesBlockedStencil(L_ALST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
      call MatSetValuesBlocked(L_ALST(2), 1, j - 1, 1, j - 1, box2, INSERT_VALUES, ierr)

      do i = 1, 2
         call MatAssemblyBegin(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set temporal ALST EVP matrix BC done!\n", ierr)

   end subroutine solve_EVP_matrix_bc_temporal

   subroutine solve_EVP_matrix_bc_spatial(j, grid_profile)

      implicit none
      integer, intent(in) :: j
      type(type_lst_grid_profile), intent(in) :: grid_profile
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, -2:2, 5), box3(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2), tmp_mat_j2(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, lj
      integer :: ljb, lje, c_idx
      MatStencil :: idxm(4, 1), idxn(4, -2:2)
      PetscViewer :: Viewer

      if (.not. (j == 1 .or. j == jn)) then
         print *, 'j=', j, ' is wrong in EVP matrix bc setting...'
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

      associate (eta_y => grid_profile%eta_y(j), eta_yy => grid_profile%eta_yy(j), &
             &  C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
         tmp_mat_j0 = 0.0d0
         tmp_mat_j1 = 0.0d0
         box1 = 0.0d0
         box2 = 0.0d0
         do lj = ljb, lje
            tmp_mat_j0(:, :, lj) = delta(lj)*(D_0(:, :, j)) + &
                                   C1(lj, c_idx)*(B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy) + &
                                   D1(lj, c_idx)*(V_0(:, :, j)*eta_y*eta_y)

            tmp_mat_j1(:, :, lj) = delta(lj)*(D_1(:, :, j)) + &
                                   C1(lj, c_idx)*(B_1(:, :, j)*eta_y)

            tmp_mat_j2(:, :, lj) = delta(lj)*(D_2(:, :, j))

            box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
            box2(:, lj, :) = transpose(tmp_mat_j1(:, :, lj))
         end do
         box3(:, :) = transpose(tmp_mat_j2(:, :, 0))
      end associate

      call MatSetValuesBlockedStencil(L_ALST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
      call MatSetValuesBlockedStencil(L_ALST(2), 1, idxm, 5, idxn, box2, INSERT_VALUES, ierr)
      call MatSetValuesBlocked(L_ALST(3), 1, j - 1, 1, j - 1, box3, INSERT_VALUES, ierr)

      do i = 1, 3
         call MatAssemblyBegin(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(L_ALST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set spatial ALST EVP matrix BC done!\n", ierr)

   end subroutine solve_EVP_matrix_bc_spatial

   subroutine solve_EVP_local(this, ome, alf, bat)

      implicit none
      class(type_alst_eqn) :: this
      complex(REAL64), intent(inout) :: ome
      complex(REAL64), intent(inout) :: alf
      complex(REAL64), intent(in) :: bat

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call Solver%Solve(2, L_ALST, ome)
         call Solver%GetEigenVal(ome)
         if (rank == 0) print('(A, 2F12.7)'), 'omega------------>', ome%re, -ome%im
      case (LST_SPATIAL_MODE)
         call Solver%Solve(3, L_ALST, alf)
         call Solver%GetEigenVal(alf)
         if (rank == 0) print('(A, 2F12.7)'), 'alpha------------>', alf%re, -alf%im
      end select

   end subroutine solve_EVP_local

   subroutine alst_finalize(this)

      implicit none
      class(type_alst_eqn), intent(inout):: this

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

      deallocate (G_0)
      deallocate (B_0)
      deallocate (D_0)
      deallocate (V_0)
      deallocate (B_1)
      deallocate (D_1)
      deallocate (D_2)

      ! alst_eqn_has_initialize=.FALSE.
      call Solver%Finalize()
      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call MatDestroy(L_ALST(1), ierr)
         call MatDestroy(L_ALST(2), ierr)
      case (LST_SPATIAL_MODE)
         call MatDestroy(L_ALST(1), ierr)
         call MatDestroy(L_ALST(2), ierr)
         call MatDestroy(L_ALST(3), ierr)
      end select

      deallocate (L_ALST)
      call DMDestroy(DM_ALST, ierr)

   end subroutine alst_finalize

end module mod_alst_eqn

! subroutine initialize_alst_eqn()

!   implicit none
!   integer :: xb1, xl1

!   jb=yb+1; je=ye+1

!   allocate(G0  (5, 5, jb:je))
!   allocate(A0  (5, 5, jb:je))
!   allocate(B0  (5, 5, jb:je))
!   allocate(C0  (5, 5, jb:je))
!   allocate(D0  (5, 5, jb:je))
!   allocate(V0xx(5, 5, jb:je))
!   allocate(V0xy(5, 5, jb:je))
!   allocate(V0xz(5, 5, jb:je))
!   allocate(V0yy(5, 5, jb:je))
!   allocate(V0yz(5, 5, jb:je))
!   allocate(V0zz(5, 5, jb:je))

!   allocate(G_0(5, 5, jb:je))
!   allocate(B_0(5, 5, jb:je))
!   allocate(D_0(5, 5, jb:je))
!   allocate(V_0(5, 5, jb:je))
!   allocate(B_1(5, 5, jb:je))
!   allocate(D_1(5, 5, jb:je))
!   allocate(D_2(5, 5, jb:je))

!   allocate(dB0   (5, 5, jb:je))
!   allocate(ddV0yy(5, 5, jb:je))
!   allocate(dV0xy (5, 5, jb:je))
!   allocate(dV0yy (5, 5, jb:je))
!   allocate(dV0yz (5, 5, jb:je))

!   allocate(dB_0 (5, 5, jb:je))
!   allocate(dV_0 (5, 5, jb:je))
!   allocate(ddV_0(5, 5, jb:je))
!   allocate(dB_1 (5, 5, jb:je))

!   alst_eqn_has_initialize=.TRUE.

!   call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_ALST, ierr)
!   call DMSetUp(DM_ALST, ierr)
!   if (lst_mode==LST_TEMPORAL_MODE)then
!     allocate(L_ALST(2))
!     call DMCreateMatrix(DM_ALST, L_ALST(1), ierr)
!     call DMCreateMatrix(DM_ALST, L_ALST(2), ierr)
!   else if(lst_mode==LST_SPATIAL_MODE)then
!     allocate(L_ALST(3))
!     call DMCreateMatrix(DM_ALST, L_ALST(1), ierr)
!     call DMCreateMatrix(DM_ALST, L_ALST(2), ierr)
!     call DMCreateMatrix(DM_ALST, L_ALST(3), ierr)
!   endif

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     allocate(slepceps_type :: Solver)
!   case (LST_SPATIAL_MODE)
!     allocate(slepcpep_type::Solver)
!   end select

!   select case (lst_method)
!   case (LST_LOCAL, LST_DOMAIN)
!     call Solver%Create(1)
!   case (LST_GLOBAL)
!     call Solver%Create(lst_nev)
!   end select

! end subroutine initialize_alst_eqn

! !> solve alst eqntions.
! subroutine solve_alst_eqn(grid_profile, baseflow_profile, ome, alf, bat, ix)

!   implicit none
!   type(type_lst_grid_profile), intent(in) :: grid_profile
!   type(type_lst_bf_profile), intent(in) :: baseflow_profile
!   complex(REAL64), intent(inout) :: ome, alf, bat
!   integer, intent(in) :: ix

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     ome=conjg(ome)
!   case (LST_SPATIAL_MODE)
!     alf=conjg(alf)
!   end select

!   ome1=ome; alf1=alf
!   iloc=ix
!   call solve_baseflow_matrix(baseflow_profile, grid_profile)
!   call solve_alst_matrix(ome, alf, bat)
!   call solve_EVP_matrix(grid_profile)
!   call solve_alst_matrix_bc(baseflow_profile, grid_profile, ome, alf, bat)
!   call solve_EVP_matrix_bc(grid_profile)

!   select case (lst_method)
!   case (LST_LOCAL, LST_DOMAIN)
!     call solve_EVP_local(ome, alf, bat)
!     select case (lst_mode)
!     case (LST_TEMPORAL_MODE)
!       ome=conjg(ome)
!     case (LST_SPATIAL_MODE)
!       alf=conjg(alf)
!     end select
!     if(rank==0)then
!       call write_alst_value_local(ome, alf, bat)
!     endif
!   case (LST_GLOBAL)
!     call solve_EVP_global(ome, alf, bat)
!     select case (lst_mode)
!     case (LST_TEMPORAL_MODE)
!       ome=conjg(ome)
!     case (LST_SPATIAL_MODE)
!       alf=conjg(alf)
!     end select
!     if(rank==0)then
!       call write_alst_value_global(ome, alf, bat)
!     endif
!   end select

! end subroutine solve_alst_eqn

! subroutine solve_EVP_global(ome, alf, bat)

!   implicit none
!   complex(REAL64), intent(inout) :: ome
!   complex(REAL64), intent(inout) :: alf
!   complex(REAL64), intent(in) :: bat
!   integer :: i

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     call Solver%Solve(2, L_ALST, ome)
!     call Solver%GetEigenVals(eig_d)
!     do i=1, size(eig_d)
!       if(rank==0) print('(A, 2F12.7)'), 'omega------------>', eig_d(i)%re, eig_d(i)%im
!     enddo
!   case (LST_SPATIAL_MODE)
!     call Solver%Solve(3, L_ALST, alf)
!     call Solver%GetEigenVals(eig_d)
!     do i=1, size(eig_d)
!       if(rank==0) print('(A, 2F12.7)'), 'alpha------------>', eig_d(i)%re, eig_d(i)%im
!     enddo
!   end select

! end subroutine solve_EVP_global

! subroutine write_alst_value_local(ome, alf, bat)

!   use global_parameters, only: xx
!   implicit none
!   complex(REAL64), intent(in) :: ome
!   complex(REAL64), intent(in) :: alf
!   complex(REAL64), intent(in) :: bat

!   write(1001, ('(I0, 12(2X, F20.8))'))iloc, xx(iloc, yb+1), Ma, Re, Te, ome, &
!           alf, bat, ome/alf
!   call flush(1001)

! end subroutine write_alst_value_local

! subroutine write_alst_value_global(ome, alf, bat)

!   use global_parameters, only: xx
!   implicit none
!   complex(REAL64), intent(in) :: ome
!   complex(REAL64), intent(in) :: alf
!   complex(REAL64), intent(in) :: bat
!   integer :: i

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     do i=1, size(eig_d)
!     write(1001, ('(I0, 12(2X, F20.8))'))iloc, xx(iloc, yb+1), Ma, Re, Te, eig_d(i), &
!             alf, bat, eig_d(i)/alf
!     enddo
!   case (LST_SPATIAL_MODE)
!     do i=1, size(eig_d)
!     write(1001, ('(I0, 12(2X, F20.8))'))iloc, xx(iloc, yb+1), Ma, Re, Te, ome, &
!             eig_d(i), bat, ome/eig_d(i)
!     enddo
!   end select

!   call flush(1001)

! end subroutine write_alst_value_global

! subroutine write_alst_function(grid_profile, baseflow_profile, ome, alf, bat, ix)

!   use stringifor
!   use global_parameters
!   implicit none
!   type(type_lst_grid_profile), intent(in) :: grid_profile
!   type(type_lst_bf_profile), intent(in) :: baseflow_profile
!   complex(REAL64), intent(inout) :: ome, alf, bat
!   integer, intent(in) :: ix

!   type(string) :: fn
!   character(255) :: tmp
!   PetscScalar :: d
!   Vec :: v, v_all
!   PetscViewer :: Viewer
!   VecScatter :: ctx
!   PetscScalar, pointer :: tmp_(:)
!   Complex(REAL64) :: phi(6, jn)
!   character(40) :: name
!   integer :: l, j, maxu_loc(1)
!   real(REAL64) :: pe
!   complex(REAL64):: maxu

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     ome=ome+(0.0d0, 1.0d-8)
!   case (LST_SPATIAL_MODE)
!     alf=alf+(0.0d0, 1.0d-8)
!   end select
!   call solve_alst_eqn(grid_profile, baseflow_profile, ome, alf, bat, ix)
!   call Solver%GetEigenVector(d, v)

!   call PetscPrintf(PETSC_COMM_WORLD, "Writing ALST Eigenfunction......\n", ierr)

!   write(tmp, ('(I0)'))ix
!   fn='ALST_eigenfunction_i='//trim(adjustl(tmp))
!   write(tmp, ('(F20.5)'))ome%re
!   fn=fn//'_ome='//trim(adjustl(tmp))
!   write(tmp, ('(F20.5)'))alf%re
!   fn=fn//'_alf='//trim(adjustl(tmp))
!   write(tmp, ('(F20.5)'))bat%re
!   fn=fn//'_bat='//trim(adjustl(tmp))//'.plt'

!   call VecScatterCreateToZero(v, ctx, v_all, ierr)
!   ! call VecGetType(v, name, ierr)
!   ! call VecGetType(v_all, name, ierr)
!   ! call VecGetSize(v, l, ierr)
!   ! call VecGetSize(v_all, l, ierr)

!   call VecScatterBegin(ctx, v, v_all,INSERT_VALUES, SCATTER_FORWARD, ierr)
!   call VecScatterEnd(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)

!   if(rank==0)then
!     call VecGetArrayReadF90(v_all, tmp_, ierr)
!     phi(1:5, 1:jn)=reshape(tmp_,[5, jn])
!     call VecRestoreArrayReadF90(v_all, tmp_, ierr)
!     maxu_loc=maxloc(abs(phi(2, :)))
!     maxu=phi(2, maxu_loc(1))
!     phi=phi/maxu
!     pe=1.0d0/(GAMMA*MA*MA)
!     do j=1, jn
!       phi(6, j)=(qq_0(ix, j, 1)*phi(5, j)+qq_0(ix, j, 5)*phi(1, j))*pe
!     enddo

!     ! block
!     !   use mod_dx
!     !   integer :: l1, l2
!     !   integer :: jloc
!     !   complex(REAL64), dimension(5, jn) :: dphidj
!     !   complex(REAL64), dimension(5, jn) :: dphidjj
!     !
!     !   jloc=149
!     !   do l=1, 5
!     !     call dx_fdm4nd_center(phi(l, :), 1, jn, dphidj(l, :), 1, jn)
!     !     call dxx_fdm4nd_center(phi(l, :), 1, jn, dphidjj(l, :), 1, jn)
!     !   enddo
!     !
!     !
!     !   alf1=d
!     !   print*, ome1, alf1
!     !   ! print*, grid_profile%eta_y(149), grid_profile%eta_yy(149)
!     !   ! print*, lbound(grid_profile%eta_y), maxu_loc
!     !   print*, phi(1:5, jloc)
!     !   print*, dphidj(1:5, jloc)
!     !   print*, dphidjj(1:5, jloc)
!     !   ! do l1=1, 5
!     !   !   do l2=1, 5
!     !   !     print*, l1, l2, C0(l1, l2, jloc), D0(l1, l2, jloc)
!     !   !   enddo
!     !   ! enddo
!     !   print*, matmul(D_0(:, :, jloc)+alf1*D_1(:, :, jloc)+alf1**2*D_2(:, :, jloc), phi(1:5, jloc))+ &
!     !   matmul(B_0(:, :, jloc)+alf1*B_1(:, :, jloc), dphidj(1:5, jloc)*grid_profile%eta_y(jloc)) +&
!     !   matmul(V_0(:, :, jloc), dphidj(:, jloc)*grid_profile%eta_yy(jloc)+dphidjj(1:5, jloc)*grid_profile%eta_y(jloc)**2)
!     !
!     !   ! matmul(-ZI*ome1*G0(:, :, 149)+ZI*alf1*A0(:, :, 149)+D0(:, :, 149)+ &
!     !   ! alf1**2*V0xx(:, :, 149), phi(1:5, 149))+ &
!     !   ! matmul(B0(:, :, 149)-ZI*alf1*V0xy(:, :, 149), dphidj(1:5, 149)*grid_profile%eta_y(149))+ &
!     !   ! matmul(-V0yy(:,:,149), dphidj(:, 149)*grid_profile%eta_yy(149)+dphidjj(1:5, 149)*grid_profile%eta_y(149)**2)
!     !
!     !   read(*, *)
!     ! endblock

!     open(1002, file='out/'//fn)
!     write(1002,('(3(A, 2E20.13, 3X))'))'#omega=', ome, 'alpha=', alf,'beta=', bat
!     write(1002, *)'variables=y, |u|, u_r, u_i, |v|, v_r, v_i, &
!                |w|, w_r, w_i, |T|, T_r, T_i, |p|, p_r, p_i, &
!                 |rho|, rho_r, rho_i, x, i'
!     do j=1, jn
!       write(1002, '(20E20.7, 2X, I0)')yy_0(ix, j), &
!       abs(phi(2, j)), phi(2, j)%re, phi(2, j)%im, &
!       abs(phi(3, j)), phi(3, j)%re, phi(3, j)%im, &
!       abs(phi(4, j)), phi(4, j)%re, phi(4, j)%im, &
!       abs(phi(5, j)), phi(5, j)%re, phi(5, j)%im, &
!       abs(phi(6, j)), phi(6, j)%re, phi(6, j)%im, &
!       abs(phi(1, j)), phi(1, j)%re, phi(1, j)%im, &
!       xx_0(ix, yb+1), ix
!     enddo
!     close(1002)
!   endif

!   call VecScale(V, 1.0d0/maxu, ierr)

!   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, &
!   'out/disturbance_adjoint.petsc', FILE_MODE_WRITE, Viewer, ierr)
!   call VecView(V, Viewer, ierr)
!   call PetscViewerDestroy(Viewer, ierr)

!   call VecScatterDestroy(ctx, ierr)
!   call VecDestroy(v_all, ierr)

! end subroutine write_alst_function

! subroutine finalize_alst_eqn()

!   implicit none

!   deallocate(G0  )
!   deallocate(A0  )
!   deallocate(B0  )
!   deallocate(C0  )
!   deallocate(D0  )
!   deallocate(V0xx)
!   deallocate(V0xy)
!   deallocate(V0xz)
!   deallocate(V0yy)
!   deallocate(V0yz)
!   deallocate(V0zz)

!   deallocate(G_0)
!   deallocate(B_0)
!   deallocate(D_0)
!   deallocate(V_0)
!   deallocate(B_1)
!   deallocate(D_1)
!   deallocate(D_2)

!   alst_eqn_has_initialize=.FALSE.
!   if(rank==0) close(1001)
!   call Solver%Finalize()
!   select case (this%lst_mode)
!   case (LST_TEMPORAL_MODE)
!     call MatDestroy(L_ALST(1), ierr)
!     call MatDestroy(L_ALST(2), ierr)
!   case (LST_SPATIAL_MODE)
!     call MatDestroy(L_ALST(1), ierr)
!     call MatDestroy(L_ALST(2), ierr)
!     call MatDestroy(L_ALST(3), ierr)
!   end select

!   deallocate(L_ALST)
!   call DMDestroy(DM_ALST, ierr)

! end subroutine finalize_alst_eqn