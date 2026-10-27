!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2023-03-23 13:01:53
!> @LastEditTime: 2023-05-18 15:57:14
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lst/lst_eqn.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lst_eqn

   use petsc

   !>常量 保留 use module
   use iso_fortran_env, only: REAL64, REAL32
   ! use global_parameters, only: yb, ye

   !>未来剔除的变量 以后逐渐废弃
   use global_parameters, only: rank
   ! use global_parameters, only: solname

   use mod_lst_grid_profile, only: type_lst_grid_profile
   use mod_lst_baseflow_profile, only: type_lst_bf_profile

   use mod_slepc_adapter, only: slepceps_type, slepcpep_type, slepc_type

   use mod_lst_solver_abstract

   implicit none

#include "../include/start.f90"

   !> LST 求解类型
   type, extends(type_lst_abstract) :: type_lst_eqn

      ! integer :: bc_y0=BC_WALL, bc_yn=BC_FARFIELD
      Mat, allocatable :: L_LST(:)
      ! integer :: lst_method=LST_LOCAL
      ! integer :: lst_mode=LST_TEMPORAL_MODE
      ! integer :: lst_nev=1
      ! integer :: Coordinate_System=COMOVING_CS

      ! real(REAL64) :: GAMMA=1.4d0, Pr=0.72d0
      ! real(REAL64) :: C1=110.4D0
      ! real(REAL64) :: C2=110.4D0
      ! real(REAL64) :: Re, Ma, Te

      ! logical :: has_initialized = .FALSE.

      ! type(type_lst_bf_profile) :: baseflow
      ! type(type_lst_grid_profile) :: grid

   contains

      ! generic :: set_parameter => lst_set_parameters_real8, &
      !                             lst_set_parameters_real4, &
      !                             lst_set_parameters_int
      procedure :: initial => lst_initial
      ! procedure :: set_grid => lst_set_grid
      ! procedure :: set_baseflow => lst_set_baseflow
      procedure :: setup => lst_setup
      procedure :: solve => lst_solve
      procedure :: get_eigf => lst_get_eigenfunction
      procedure :: get_matrix => get_matrix_lst

      ! procedure, private :: lst_set_parameters_real8
      ! procedure, private :: lst_set_parameters_real4
      ! procedure, private :: lst_set_parameters_int

      procedure, private :: solve_baseflow_matrix
      procedure, private :: solve_bf_matrix_j

      procedure, private :: solve_lst_matrix
      procedure, private :: solve_lst_matrix_j
      procedure, private :: solve_lst_matrix_j0
      procedure, private :: solve_lst_matrix_jn
      procedure, private :: solve_lst_matrix_bc

      procedure, private :: solve_evp_local
      procedure, private :: solve_evp_matrix
      procedure, private :: solve_evp_matrix_temporal
      procedure, private :: solve_evp_matrix_spatial
      procedure, private :: solve_EVP_matrix_bc
      procedure, private :: solve_evp_matrix_bc_temporal
      procedure, private :: solve_evp_matrix_bc_spatial

      procedure :: finalize => lst_finalize

   end type type_lst_eqn

   ! logical :: lst_eqn_has_initialize=.FALSE.

   ! integer :: iloc
   integer :: jb, je
   integer :: jn
   integer :: yb, ye

   DM :: DM_LST

   class(slepc_type), allocatable :: Solver

   complex(REAL64), allocatable :: eig_d(:)

   ! complex(REAL64) :: ome1, alf1

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

   real(REAL64), allocatable, dimension(:, :, :) :: G0, A0, B0, C0
   complex(REAL64), allocatable, dimension(:, :, :) :: D0
   real(REAL64), allocatable, dimension(:, :, :) :: V0xx, V0xy, V0xz
   real(REAL64), allocatable, dimension(:, :, :) :: V0yy, V0yz, V0zz
   logical :: solved_lst_matrix = .FALSE.

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

   PetscErrorCode :: ierr

   private

   public :: type_lst_eqn, lst_eigen_error, lst_track_target, lst_is_tracking, lst_copy_selected_eigenvector
   public :: lst_select_eigenvector_by_overlap, lst_prepare_candidate_retry

contains

   subroutine lst_track_target(enabled)
      logical, intent(in) :: enabled
      if (allocated(Solver)) Solver%track_target=enabled
   end subroutine lst_track_target

   logical function lst_is_tracking()
      lst_is_tracking=.false.
      if (allocated(Solver)) lst_is_tracking=Solver%track_target
   end function lst_is_tracking

   function lst_eigen_error() result(error)
      real(REAL64) :: error
      error=huge(1.0d0)
      if (allocated(Solver)) error=Solver%last_error
   end function lst_eigen_error

   !> Owned raw eigenvector copy, eigenvalue and residual of the already selected
   !> pair. This avoids export normalization/scatter and never reselects a mode.
   !> On valid=.true. the caller must call VecDestroy(eigenvector).
   subroutine lst_copy_selected_eigenvector(eigenvalue, eigenvector, error, valid)
      complex(REAL64), intent(out) :: eigenvalue
      Vec, intent(out) :: eigenvector
      real(REAL64), intent(out) :: error
      logical, intent(out) :: valid
      valid = .false.
      eigenvalue = cmplx(0.0_REAL64, 0.0_REAL64, REAL64)
      eigenvector = PETSC_NULL_VEC
      error = huge(1.0_REAL64)
      if (.not. allocated(Solver)) return
      call Solver%CopySelectedEigenVector(eigenvalue, eigenvector, error, valid)
   end subroutine lst_copy_selected_eigenvector

   !> Explicit multi-candidate overlap selection. owned_vector belongs to the
   !> caller only when valid=.true.; the caller must then VecDestroy it.
   subroutine lst_select_eigenvector_by_overlap(reference, reference_eigen, jump_tol, error_tol, min_overlap, &
                                              selected_eigen, owned_vector, error, overlap, valid)
      Vec, intent(in) :: reference
      complex(REAL64), intent(in) :: reference_eigen
      real(REAL64), intent(in) :: jump_tol, error_tol, min_overlap
      complex(REAL64), intent(out) :: selected_eigen
      Vec, intent(out) :: owned_vector
      real(REAL64), intent(out) :: error, overlap
      logical, intent(out) :: valid
      selected_eigen = cmplx(0.0_REAL64, 0.0_REAL64, REAL64)
      owned_vector = PETSC_NULL_VEC
      error = huge(1.0_REAL64); overlap = 0.0_REAL64; valid = .false.
      if (.not. allocated(Solver)) return
      call Solver%SelectEigenVectorByOverlap(reference, reference_eigen, jump_tol, error_tol, min_overlap, &
                                            selected_eigen, owned_vector, error, overlap, valid)
   end subroutine lst_select_eigenvector_by_overlap

   !> One-shot spatial PEP candidate expansion. The solver copies reference;
   !> the caller keeps ownership and may destroy or reuse its own vector.
   subroutine lst_prepare_candidate_retry(count, reference_vector)
      integer, intent(in) :: count
      Vec, intent(in) :: reference_vector
      if (.not. allocated(Solver)) error stop 'Candidate retry requires an allocated LST solver.'
      call Solver%PrepareCandidateRetry(count, reference_vector)
   end subroutine lst_prepare_candidate_retry

   subroutine lst_setup(this)
      implicit none
      class(type_lst_eqn), intent(inout) :: this

      call this%solve_baseflow_matrix()

   end subroutine lst_setup

   subroutine lst_initial(this, jn_, Re, Ma, Te, yb_, ye_)

      implicit none
      class(type_lst_eqn), intent(inout) :: this
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

      this%has_initialized = .TRUE.

#if PETSC_VERSION_LT(3,22,0)
      call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_LST, ierr)
#else
      call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER_ARRAY, DM_LST, ierr)
#endif 
      call DMSetUp(DM_LST, ierr)
      if (this%lst_mode == LST_TEMPORAL_MODE) then
         allocate (this%L_LST(2))
         call DMCreateMatrix(DM_LST, this%L_LST(1), ierr)
         call DMCreateMatrix(DM_LST, this%L_LST(2), ierr)
      else if (this%lst_mode == LST_SPATIAL_MODE) then
         allocate (this%L_LST(3))
         call DMCreateMatrix(DM_LST, this%L_LST(1), ierr)
         call DMCreateMatrix(DM_LST, this%L_LST(2), ierr)
         call DMCreateMatrix(DM_LST, this%L_LST(3), ierr)
      end if

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         allocate (slepceps_type :: Solver)
      case (LST_SPATIAL_MODE)
         allocate (slepcpep_type::Solver)
      end select

      call Solver%Create(1, this%lst_mode)

   end subroutine lst_initial

   subroutine lst_solve(this, ome, alf, bat)

      implicit none
      class(type_lst_eqn), intent(in) :: this
      complex(REAL64), intent(inout) :: ome, alf, bat

      call this%solve_lst_matrix(ome, alf, bat)

      call this%solve_EVP_matrix()
      call this%solve_lst_matrix_bc(ome, alf, bat)
      solved_lst_matrix =.TRUE.
      call this%solve_EVP_matrix_bc()
      call this%solve_EVP_local(ome, alf, bat)

   end subroutine lst_solve

   subroutine lst_get_eigenfunction(this, eigf)

      use mod_lst_solver_abstract, only: normal_u
      implicit none
      class(type_lst_eqn), intent(in) :: this
      Vec, intent(out) :: eigf
      PetscScalar :: d

      call Solver%GetEigenVector(d, eigf)
      call normal_u(eigf)

   end subroutine lst_get_eigenfunction

   subroutine get_matrix_lst(this, ome, alf, bat, mat_name, matrix)

      implicit none
      class(type_lst_eqn), intent(inout) :: this
      complex(REAL64), intent(in) :: ome, alf, bat
      character(len=*), intent(in) :: mat_name
      complex(REAL64), intent(out) :: matrix(5, 5, jn)

      select case (trim(mat_name))
      case ('A_0')
         if (.not. solved_lst_matrix) then
            call this%solve_lst_matrix(ome, alf, bat)
            call this%solve_lst_matrix_bc(ome, alf, bat)
         endif
            matrix = A0 - ZI*bat*V0xz- 2.0d0*ZI*alf*V0xx
      case default
         stop "A wrong matrix name is given in get_matrix@lst_eqn."
      end select

   end subroutine get_matrix_lst

   subroutine solve_baseflow_matrix(this)

      implicit none
      class(type_lst_eqn) :: this
      integer :: j, j_b, j_e

      !! internal point
      call this%solve_bf_matrix_j(jb, je)

   end subroutine solve_baseflow_matrix

   subroutine solve_lst_matrix_bc(this, ome, alf, bat)

      implicit none
      class(type_lst_eqn) :: this
      complex(REAL64), intent(in) :: ome, alf, bat

      !! b.c. point
      if (jb == 1) call this%solve_lst_matrix_j0(jb, ome, alf, bat)
      if (je == jn) call this%solve_lst_matrix_jn(je, ome, alf, bat)

   end subroutine solve_lst_matrix_bc

   subroutine solve_bf_matrix_j(this, jb, je)

      implicit none
      class(type_lst_eqn) :: this
      integer, intent(in) :: jb, je

      real(REAL64) :: u, uy, uyy, w, wy, wyy
      real(REAL64) :: rho, rhoy, T, Ty, Tyy

      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      real(REAL64) :: p12, p23, p12y, p23y, phi
      real(REAL64) :: pe, cm, ck
      real(REAL64) :: rm, rmy, rmt, rmtt, rmty
      real(REAL64) :: rk, rkt, rktt, rky, rkty
      real(REAL64) :: cv, cp, rek

      integer :: j
      real(REAL64), parameter :: Z2D3 = 2.d0/3.d0, D3 = 1.d0/3.d0

      associate (Re => this%Re, Ma => this%Ma, Te => this%Te, &
                 grid_profile => this%grid, bf_profile => this%baseflow, &
                 GAMMA => this%GAMMA, C1 => this%C1, C2 => this%C2, Pr => this%Pr)

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
               ! print*, 'comoving cs is used.'
            case (CYLINDRICAL_CS)
               hx = 1.d0
               hx_y = 0.0d0
               d12 = 0.0d0

               hz = y
               if(abs(hz)<= epsilon(1.0d0)) hz = epsilon(1.0d0)
               hz_y = 1.0d0
               d32 = hz_y/hz
               ! print*, 'cylinderical cs is used.', hz, d32
            end select

            d132 = d12 + d32

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

            rk = T*sqrt(T)*(1.d0 + ck)/(T + ck)/rek
            rkt = rk*(1.5d0/T - 1.0d0/(T + ck))
            rktt = rkt*(1.5d0/T - 1.0d0/(T + ck)) - rk*(1.5d0/T**2 - 1.0d0/(T + ck)**2)
            rky = rkt*ty
            rkty = rktt*ty

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
            D0(5, 5, j) = -Ty*d132*rkt - phi*rmt - Tyy*rkt - rktt*Ty*Ty !-rkty*Ty

            V0xx(2, 2, j) = 2.0*Z2D3*rm
            V0xx(3, 3, j) = rm
            V0xx(4, 4, j) = rm
            V0xx(5, 5, j) = rk

            V0yy(2, 2, j) = rm
            V0yy(3, 3, j) = 2.0d0*Z2D3*rm
            V0yy(4, 4, j) = rm
            V0yy(5, 5, j) = rk

            V0zz(2, 2, j) = rm
            V0zz(3, 3, j) = rm
            V0zz(4, 4, j) = 2.0d0*Z2D3*rm
            V0zz(5, 5, j) = rk

            V0xy(2, 3, j) = D3*rm
            V0xy(3, 2, j) = D3*rm

            V0xz(2, 4, j) = D3*rm
            V0xz(4, 2, j) = D3*rm

            V0yz(3, 4, j) = D3*rm
            V0yz(4, 3, j) = D3*rm

            A0(:, :, j) = A0(:, :, j)/hx
            C0(:, :, j) = C0(:, :, j)/hz
            V0xx(:, :, j) = V0xx(:, :, j)/hx**2
            V0xy(:, :, j) = V0xy(:, :, j)/hx
            V0xz(:, :, j) = V0xz(:, :, j)/hx/hz
            V0yz(:, :, j) = V0yz(:, :, j)/hz
            V0zz(:, :, j) = V0zz(:, :, j)/hz**2
         end do

      end associate

   end subroutine solve_bf_matrix_j

   subroutine solve_lst_matrix_j0(this, jb, ome, alf, bat)

      implicit none
      class(type_lst_eqn) :: this
      integer, intent(in) :: jb
      ! type(type_lst_grid_profile), intent(in) :: grid_profile
      ! type(type_lst_bf_profile), intent(in) ::  bf_profile
      complex(REAL64), intent(in) :: ome, alf, bat

      real(REAL64) :: rho, u, w, rhoy
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      real(REAL64) :: alfr(5), alfi(5), beta(5)
      real(REAL64) :: vl(5, 5), vr(5, 5)
      real(REAL64) :: lambda_(5), abs_lambda_(5), diag_plus(5, 5),diag_minus(5, 5), matrix_plus(5, 5)
      integer :: ipiv(5), info
      real(REAL64) :: work(100)
      integer :: l


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
            d12 = 0.0d0

            hz = y
            if(abs(hz)<=epsilon(1.0d0)) hz = epsilon(1.0d0)
            hz_y = 1.0d0
            d32 = hz_y/hz
         end select

         d132 = d12 + d32

         select case (this%bc_y0)
         case (BC_WALL)
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

            G0(1, 1, 1) = 1.0d0
            A0(1, 1, 1) = u
            A0(1, 2, 1) = rho
            B0(1, 3, 1) = rho
            C0(1, 1, 1) = w
            C0(1, 4, 1) = rho
            D0(1, 3, 1) = rhoy + rho*d132
            D0(2, 2, 1) = 1.0d0
            D0(3, 3, 1) = 1.0d0
            D0(4, 4, 1) = 1.0d0
            if (this%isothermal) then
               D0(5, 5, 1) = 1.0d0
            else
               B0(5, 5, 1) = 1.0d0
            end if
            print*, "wall bc is used."

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
               print*, "2D-POLE bc is used."
            else
               D0(1, 1, 1) = 1.0d0
               D0(2, 2, 1) = 1.0d0
               D0(5, 5, 1) = 1.0d0
               B0(3, 3, 1) = 1.0d0
               B0(4, 4, 1) = 1.0d0
               print*, "3D-POLE bc is used."
            end if
         
         case (BC_POROUS)
            block
               use Complex_Bessel
               use mod_porous_parameters

               complex(REAL64) :: porous_A, porous_B

!!! d, r, d is parameters
               
               associate (pe=> 1.0d0/(this%GAMMA*this%MA**2), &
                  &       Rho0 => this%baseflow%q0(1, 1), &
                  &       T0=>this%baseflow%q0(5, 1))
                  ! porous_A = get_porous_A

                  call get_porous_AB(Rho0, T0, ome, porous_A, porous_B)
                  ! print*, porous_A, porous_B
                  ! print*, rho0, T0, pe
                  ! print*, 'bc_porous bc is selected'
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
      
                  G0(1, 1, 1) = 1.0d0
                  A0(1, 1, 1) = u
                  A0(1, 2, 1) = rho
                  B0(1, 3, 1) = rho
                  C0(1, 1, 1) = w
                  C0(1, 4, 1) = rho
                  D0(1, 3, 1) = rhoy + rho*d132

                  D0(2, 2, 1) = 1.0d0
                  D0(3, 1, 1) = -porous_A*T0*pe
                  D0(3, 3, 1) = 1.0d0
                  D0(3, 5, 1) = -porous_A*Rho0*pe
                  D0(4, 4, 1) = 1.0d0
                  D0(5, 1, 1) = -porous_B*T0*pe
                  D0(5, 5, 1) = 1.0d0-porous_B*Rho0*pe
               end associate
            end block
         case (BC_NONREFLECTION)
   
            !!! G^(-1)B = S^(-1)\Lambda S
            call dggev('N', 'V', 5, real(B0(:, :, jn), kind=real64), 5, &
            &          real(G0(:, :, jn), kind=real64), 5, &
            &          alfr, alfi, beta, vl, 5, vr, 5, work, 100, info)
   
            lambda_ = (alfr)/beta
            print*, lambda_
            !!! G^(-1)B^(+) = S^(-1)\Lambda^(+) S
            diag_plus= 0.0d0
            diag_minus = 0.0d0
            do l = 1, 5
               abs_lambda_(l) = real(sqrt(lambda_(l)**2), kind=REAL64)
               diag_plus(l, l) = 0.5d0*(lambda_(l) + abs_lambda_(l))
               diag_minus(l, l) = 0.5d0*(lambda_(l) - abs_lambda_(l))
            end do
            vl = vr !! vl is inv(vr)
            call dgetrf(5, 5, vl, 5, ipiv, info)
            call dgetri(5, vl, 5, ipiv, work, 5, info)
            !!!  B^(+) = G S^(-1)\Lambda^(-) S
            B0(:, :, jn) = matmul(G0(:, :, jn), vr)
            B0(:, :, jn) = matmul(B0(:, :, jn), diag_minus)
            B0(:, :, jn) = matmul(B0(:, :, jn), vl)
         end select

         A0(:, :, 1) = A0(:, :, 1)/hx
         C0(:, :, 1) = C0(:, :, 1)/hz

      end associate

      call this%solve_lst_matrix_j(1, ome, alf, bat)

   end subroutine solve_lst_matrix_j0

   subroutine solve_lst_matrix_jn(this, je, ome, alf, bat)

      implicit none
      class(type_lst_eqn) :: this
      integer, intent(in) :: je
      ! type(type_lst_grid_profile), intent(in) :: grid_profile
      ! type(type_lst_bf_profile), intent(in) ::  bf_profile
      complex(REAL64), intent(in) :: ome, alf, bat

      real(REAL64) :: rho, u, w, rhoy
      real(REAL64) :: y, drx, drz
      real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

      real(REAL64) :: alfr(5), alfi(5), beta(5)
      real(REAL64) :: vl(5, 5), vr(5, 5)
      real(REAL64) :: lambda_(5), abs_lambda_(5), diag_plus(5, 5), matrix_plus(5, 5)
      integer :: ipiv(5), info
      real(REAL64) :: work(100)
      integer :: l

      real(REAL64), parameter :: epsilon = 1.0d-6

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
         D0(1, 1, jn) = 1.0d0
         D0(2, 2, jn) = 1.0d0
         D0(3, 3, jn) = 1.0d0
         D0(4, 4, jn) = 1.0d0
         D0(5, 5, jn) = 1.0d0
      case (BC_NONREFLECTION)

         ! associate(rho => this%baseflow%q0(1, jn), T=>this%baseflow%q0(5, jn), &
         !    & pe=> 1.0d0/(this%GAMMA*this%Ma**2))
            ! B0(1, 3, jn) = rho
            ! B0(2, 2, jn) = 0.0d0
            ! B0(2, 5, jn) = 0.0d0
            ! B0(3, 1, jn) = pe*T
            ! B0(3, 3, jn) = 0.0d0
            ! B0(3, 5, jn) = pe*rho
            ! B0(4, 4, jn) = 0.0d0
            ! B0(4, 5, jn) = 0.0d0
            ! B0(5, 2, jn) = 0.0d0
            ! B0(5, 4, jn) = 0.0d0
            ! B0(5, 5, jn) = 0.0d0
         ! end associate

         !!! G^(-1)B = S^(-1)\Lambda S
         call dggev('N', 'V', 5, real(B0(:, :, jn), kind=real64), 5, &
         &          real(G0(:, :, jn), kind=real64), 5, &
         &          alfr, alfi, beta, vl, 5, vr, 5, work, 100, info)

         lambda_ = (alfr)/beta
         print*, lambda_
         !!! G^(-1)B^(+) = S^(-1)\Lambda^(+) S
         diag_plus= 0.0d0
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
         ! block
         !    integer :: l1, l2
         !    do l1=1, 5
         !       do l2=1, 5
         !          print*, l1, l2, B0(l1, l2, jn), matrix_plus(l1, l2)
         !       enddo
         !    enddo 
         !    stop

         ! endblock
            !  matrix_minus = matmul(matrixG, vr)
            !  matrix_minus = matmul(matrix_minus, diag_minus)
            !  matrix_minus = matmul(matrix_minus, vl)
      end select
      call this%solve_lst_matrix_j(jn, ome, alf, bat)

   end subroutine solve_lst_matrix_jn

   subroutine solve_lst_matrix(this, ome, alf, bat)
      implicit none
      class(type_lst_eqn) :: this
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      integer :: j

      do j = jb, je
         call this%solve_lst_matrix_j(j, ome, alf, bat)
      end do

   end subroutine solve_lst_matrix

   subroutine solve_lst_matrix_j(this, j, ome, alf, bat)

      implicit none
      class(type_lst_eqn) :: this
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      G_0(:, :, j) = -ZI*G0(:, :, j)
      D_0(:, :, j) = D0(:, :, j) + ZI*bat*C0(:, :, j) + bat**2*V0zz(:, :, j)
      B_0(:, :, j) = B0(:, :, j) - ZI*bat*V0yz(:, :, j)
      V_0(:, :, j) = -V0yy(:, :, j)
      D_1(:, :, j) = ZI*A0(:, :, j) + bat*V0xz(:, :, j)
      B_1(:, :, j) = -ZI*V0xy(:, :, j)
      D_2(:, :, j) = +V0xx(:, :, j)

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call solve_lst_matrix_temporal_j(j, ome, alf, bat)
      case (LST_SPATIAL_MODE)
         call solve_lst_matrix_spatial_j(j, ome, alf, bat)
      end select

   end subroutine solve_lst_matrix_j

   subroutine solve_lst_matrix_temporal_j(j, ome, alf, bat)

      implicit none
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      D_0(:, :, j) = D_0(:, :, j) + alf*D_1(:, :, j) + alf**2*D_2(:, :, j)
      B_0(:, :, j) = B_0(:, :, j) + alf*B_1(:, :, j)
      D_1(:, :, j) = G_0(:, :, j)

   end subroutine solve_lst_matrix_temporal_j

   subroutine solve_lst_matrix_spatial_j(j, ome, alf, bat)

      implicit none
      integer, intent(in) :: j
      complex(REAL64), intent(in) :: ome
      complex(REAL64), intent(in) :: alf
      complex(REAL64), intent(in) :: bat

      D_0(:, :, j) = D_0(:, :, j) + ome*G_0(:, :, j)

   end subroutine solve_lst_matrix_spatial_j

   subroutine solve_EVP_matrix(this)

      implicit none
      class(type_lst_eqn) :: this

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call this%solve_EVP_matrix_temporal()
      case (LST_SPATIAL_MODE)
         call this%solve_EVP_matrix_spatial()
      end select

   end subroutine solve_EVP_matrix

   subroutine solve_EVP_matrix_bc(this)

      implicit none
      class(type_lst_eqn) :: this

      integer :: j

      do j = jb, je
         if (j == 1 .or. j == jn) then
            select case (this%lst_mode)
            case (LST_TEMPORAL_MODE)
               call this%solve_EVP_matrix_bc_temporal(j)
            case (LST_SPATIAL_MODE)
               call this%solve_EVP_matrix_bc_spatial(j)
            end select
         end if
      end do

   end subroutine solve_EVP_matrix_bc

   subroutine solve_EVP_matrix_temporal(this)

      implicit none
      class(type_lst_eqn) :: this
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2)

      integer :: i, j, lj
      integer :: ljb, lje, c_idx
      integer :: single_idx(1) 
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

         associate (eta_y => this%grid%eta_y(j), eta_yy => this%grid%eta_yy(j), &
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

         single_idx(1) = j - 1

         call MatSetValuesBlockedStencil(this%L_LST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
         call MatSetValuesBlocked(this%L_LST(2), 1, single_idx, 1, single_idx, box2, INSERT_VALUES, ierr)
      end do

      do i = 1, 2
         call MatAssemblyBegin(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set temporal LST EVP matrix done!\n", ierr)

   end subroutine solve_EVP_matrix_temporal

   subroutine solve_EVP_matrix_spatial(this)

      implicit none
      class(type_lst_eqn) :: this
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, -2:2, 5), box3(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2), tmp_mat_j2(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, j, lj
      integer :: ljb, lje, c_idx
      integer :: single_idx(1) 
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

         associate (eta_y => this%grid%eta_y(j), eta_yy => this%grid%eta_yy(j), &
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

         single_idx(1) = j - 1

         call MatSetValuesBlockedStencil(this%L_LST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
         call MatSetValuesBlockedStencil(this%L_LST(2), 1, idxm, 5, idxn, box2, INSERT_VALUES, ierr)
         call MatSetValuesBlocked(this%L_LST(3), 1, single_idx, 1, single_idx, box3, INSERT_VALUES, ierr)
      end do

      do i = 1, 3
         call MatAssemblyBegin(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LST EVP matrix done!\n", ierr)

   end subroutine solve_EVP_matrix_spatial

   subroutine solve_EVP_matrix_bc_temporal(this, j)

      implicit none
      class(type_lst_eqn) :: this
      integer, intent(in) :: j
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2)

      integer :: i, lj
      integer :: ljb, lje, c_idx
      integer :: single_idx(1) 
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

      associate (eta_y => this%grid%eta_y(j), eta_yy => this%grid%eta_yy(j), &
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

      single_idx(1) = j - 1 

      call MatSetValuesBlockedStencil(this%L_LST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
      call MatSetValuesBlocked(this%L_LST(2), 1, single_idx, 1, single_idx, box2, INSERT_VALUES, ierr)

      do i = 1, 2
         call MatAssemblyBegin(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set temporal LST EVP matrix BC done!\n", ierr)

   end subroutine solve_EVP_matrix_bc_temporal

   subroutine solve_EVP_matrix_bc_spatial(this, j)

      implicit none
      integer, intent(in) :: j
      class(type_lst_eqn) :: this
      real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

      complex(REAL64) :: box1(5, -2:2, 5), box2(5, -2:2, 5), box3(5, 5)
      complex(REAL64) :: tmp_mat_j0(5, 5, -2:2), tmp_mat_j1(5, 5, -2:2), tmp_mat_j2(5, 5, -2:2)
      ! real(REAL64) :: eta_y, eta_yy

      integer :: i, lj
      integer :: ljb, lje, c_idx
      integer :: single_idx(1)
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

      associate (eta_y => this%grid%eta_y(j), eta_yy => this%grid%eta_yy(j), &
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

      single_idx(1) = j - 1 

      call MatSetValuesBlockedStencil(this%L_LST(1), 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
      call MatSetValuesBlockedStencil(this%L_LST(2), 1, idxm, 5, idxn, box2, INSERT_VALUES, ierr)
      call MatSetValuesBlocked(this%L_LST(3), 1, single_idx, 1, single_idx, box3, INSERT_VALUES, ierr)

      do i = 1, 3
         call MatAssemblyBegin(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(this%L_LST(i), MAT_FINAL_ASSEMBLY, ierr)
      end do

      call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LST EVP matrix BC done!\n", ierr)

   end subroutine solve_EVP_matrix_bc_spatial

   subroutine solve_EVP_local(this, ome, alf, bat)

      implicit none
      class(type_lst_eqn) :: this
      complex(REAL64), intent(inout) :: ome
      complex(REAL64), intent(inout) :: alf
      complex(REAL64), intent(in) :: bat

      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call Solver%Solve(2, this%L_LST, ome)
         call Solver%GetEigenVal(ome)
         if (rank == 0) print('(A, 2F12.7)'), 'omega------------>', ome%re, ome%im
      case (LST_SPATIAL_MODE)
         ! block
         !    Vec:: eigfunc, tmp
         !    PetscViewer :: viewer

         !    Mat :: L_tmp
         !    Real(8) :: norm

         !    call VecCreate(PETSC_COMM_WORLD, eigfunc, ierr)
         !    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "out/disturbance.petsc",  FILE_MODE_READ, viewer, ierr)
         !    call VecLoad(eigfunc, viewer, ierr)
         !    call PetscViewerDestroy(viewer, ierr)

         !    call MatDuplicate(this%L_LST(1), MAT_COPY_VALUES, L_tmp, ierr)
         !    call MatAXPY(L_tmp, alf, this%L_LST(2), SAME_NONZERO_PATTERN, ierr)
         !    call MatAXPY(L_tmp, alf**2, this%L_LST(3), SAME_NONZERO_PATTERN, ierr)
         !    call VecDuplicate(eigfunc, tmp, ierr)
         !    call MatMult(L_tmp, eigfunc, tmp, ierr)
         !    call VecNorm(tmp, NORM_2, norm, ierr)
         !    print*, norm
         ! endblock
         call Solver%Solve(3, this%L_LST, alf)
         call Solver%GetEigenVal(alf)
         if (rank == 0) print('(A, 2F12.7)'), 'alpha------------>', alf%re, alf%im
      end select

   end subroutine solve_EVP_local

   subroutine lst_finalize(this)

      implicit none
      class(type_lst_eqn), intent(inout) :: this

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

      ! this%has_initialized=.FALSE.
      ! if(rank==0) close(1001)
      call Solver%Finalize()
      select case (this%lst_mode)
      case (LST_TEMPORAL_MODE)
         call MatDestroy(this%L_LST(1), ierr)
         call MatDestroy(this%L_LST(2), ierr)
      case (LST_SPATIAL_MODE)
         call MatDestroy(this%L_LST(1), ierr)
         call MatDestroy(this%L_LST(2), ierr)
         call MatDestroy(this%L_LST(3), ierr)
      end select

      deallocate (this%L_LST)
      call DMDestroy(DM_LST, ierr)
      solved_lst_matrix=.FALSE.

   end subroutine lst_finalize

end module mod_lst_eqn

! subroutine initialize_lst_eqn()

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

!   lst_eqn_has_initialize=.TRUE.

!   call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_LST, ierr)
!   call DMSetUp(DM_LST, ierr)
!   if (lst_mode==LST_TEMPORAL_MODE)then
!     allocate(this%L_LST(2))
!     call DMCreateMatrix(DM_LST, this%L_LST(1), ierr)
!     call DMCreateMatrix(DM_LST, this%L_LST(2), ierr)
!   else if(lst_mode==LST_SPATIAL_MODE)then
!     allocate(this%L_LST(3))
!     call DMCreateMatrix(DM_LST, this%L_LST(1), ierr)
!     call DMCreateMatrix(DM_LST, this%L_LST(2), ierr)
!     call DMCreateMatrix(DM_LST, this%L_LST(3), ierr)
!   endif

!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     allocate(slepceps_type :: Solver)
!   case (LST_SPATIAL_MODE)
!     allocate(slepcpep_type::Solver)
!   end select

!   ! select case (lst_method)
!   ! case (LST_LOCAL, LST_DOMAIN)
!     call Solver%Create(1)
!   ! case (LST_GLOBAL)
!   !   call Solver%Create(lst_nev)
!   ! end select

! end subroutine initialize_lst_eqn

! !> solve lst eqntions.
! subroutine solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, ix)

!   implicit none
!   type(type_lst_grid_profile), intent(in) :: grid_profile
!   type(type_lst_bf_profile), intent(in) :: baseflow_profile
!   complex(REAL64), intent(inout) :: ome, alf, bat
!   integer, intent(in) :: ix

!   ! ome1=ome; alf1=alf
!   iloc=ix
!   call solve_baseflow_matrix(baseflow_profile, grid_profile)
!   call solve_lst_matrix(ome, alf, bat)
!   call solve_EVP_matrix(grid_profile)
!   call solve_lst_matrix_bc(baseflow_profile, grid_profile, ome, alf, bat)
!   call solve_EVP_matrix_bc(grid_profile)

!   ! select case (lst_method)
!   ! case (LST_LOCAL, LST_DOMAIN)
!   call solve_EVP_local(ome, alf, bat)
!   ! if(rank==0)then
!   !   call write_lst_value_local(ome, alf, bat)
!   ! endif
!   ! case (LST_GLOBAL)
!   !   call solve_EVP_global(ome, alf, bat)
!   !   if(rank==0)then
!   !     call write_lst_value_global(ome, alf, bat)
!   !   endif
!   ! end select

! end subroutine solve_lst_eqn

! subroutine solve_EVP_global(tome, alf, bat)

!   implicit none
!   complex(REAL64), intent(inout) :: ome
!   complex(REAL64), intent(inout) :: alf
!   complex(REAL64), intent(in) :: bat
!   integer :: i

!   select case (this%lst_mode)
!   case (LST_TEMPORAL_MODE)
!     call Solver%Solve(2, this%L_LST, ome)
!     call Solver%GetEigenVals(eig_d)
!     do i=1, size(eig_d)
!       if(rank==0) print('(A, 2F12.7)'), 'omega------------>', eig_d(i)%re, eig_d(i)%im
!     enddo
!   case (LST_SPATIAL_MODE)
!     call Solver%Solve(3, this%L_LST, alf)
!     call Solver%GetEigenVals(eig_d)
!     do i=1, size(eig_d)
!       if(rank==0) print('(A, 2F12.7)'), 'alpha------------>', eig_d(i)%re, eig_d(i)%im
!     enddo
!   end select

! end subroutine solve_EVP_global

! subroutine write_lst_value_local(ome, alf, bat)

!   use global_parameters, only: xx
!   implicit none
!   complex(REAL64), intent(in) :: ome
!   complex(REAL64), intent(in) :: alf
!   complex(REAL64), intent(in) :: bat

! end subroutine write_lst_value_local

! subroutine write_lst_value_global(ome, alf, bat)

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

! end subroutine write_lst_value_global

! subroutine finalize_lst_eqn()

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

!   lst_eqn_has_initialize=.FALSE.
!   if(rank==0) close(1001)
!   call Solver%Finalize()
!   select case (lst_mode)
!   case (LST_TEMPORAL_MODE)
!     call MatDestroy(this%L_LST(1), ierr)
!     call MatDestroy(this%L_LST(2), ierr)
!   case (LST_SPATIAL_MODE)
!     call MatDestroy(this%L_LST(1), ierr)
!     call MatDestroy(this%L_LST(2), ierr)
!     call MatDestroy(this%L_LST(3), ierr)
!   end select

!   deallocate(this%L_LST)
!   call DMDestroy(DM_LST, ierr)

! end subroutine finalize_lst_eqn
