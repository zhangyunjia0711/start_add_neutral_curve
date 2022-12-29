!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-09 10:45:42
!> @LastEditTime: 2022-12-12 19:34:26
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: adjoint linearized parabolized stability equations solver.
!> @FilePath: /START_all/src/lpse/alpse_eqn.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>
!@ apse equations solver
module mod_alpse_eqn

    use petsc
    ! use global_parameters, only: iother=>iloc
    ! use global_parameters, only: phi_t=>phi
    ! use global_parameters, only: wave_param_t=> wave_param
    ! use global_parameters

    use iso_fortran_env, only: REAL64

    ! use global_parameters, only: Re, Ma, Te, Pr, C1, C2, GAMMA
    ! use global_parameters, only: solve_alpha
    ! use global_parameters, only: sigma_control, sigma_p
    ! use mod_lpse_grid_profile
    ! use mod_lpse_baseflow_profile
    use mod_lst_baseflow_profile, only: type_lst_bf_profile => type_lst_bf_profile
    use mod_lst_grid_profile, only: type_lst_grid_profile => type_lst_grid_profile

    use mod_lpse_solver_abstract

    implicit none

    private
    public :: type_alpse_eqn


#include "../include/start.f90"
    logical :: alpse_eqn_has_initialize = .FALSE.
    integer :: iloc
    integer :: jb, je
    integer :: jn
    integer :: yb, ye

    Mat :: L_ALPSE
    Vec :: RHS
    KSP :: solver
    PC :: SolverPC
    PetscErrorCode :: ierr

    complex(REAL64) :: alf_old

    real(REAL64), allocatable, dimension(:, :, :) :: G0, A0, B0, C0, D0
    real(REAL64), allocatable, dimension(:, :, :) :: V0xx, V0xy, V0xz
    real(REAL64), allocatable, dimension(:, :, :) :: V0yy, V0yz, V0zz
    real(REAL64), allocatable, dimension(:, :, :) :: Aplus
    real(REAL64), allocatable, dimension(:, :, :) :: dA0, dB0, dV0yy, ddV0yy, dV0xy, dV0yz

    real(REAL64), allocatable, dimension(:, :) :: modify_Px_A
  !> \f
  !! \phi=\hat\phi(y)
  !!  L\hat\phi=0
  !!  L=
  !!  (+D_0 + B_0\frac{\partial}{\parital y}+ V_0\frac{\partial^2}{\parital y^2})
  !!  +\omega*(G_0)
  !!  +alpha*(D_1+B_1\frac{\partial}{\parital y})
  !!  +alpha^2*(D_2)
  !! \f
    complex(REAL64), allocatable, dimension(:, :, :) :: B_0, D_0, V_0
    complex(REAL64), allocatable, dimension(:, :, :) :: D_2
    complex(REAL64), allocatable, dimension(:, :, :) :: dA_0, DB_0, dV_0, ddV_0

    complex(REAL64), allocatable, dimension(:, :) :: phi_m1, phi_m2

    real(REAL64) :: lambda !< control d/dx order

    type, extends(type_lpse_abstract) :: type_alpse_eqn

    contains

        procedure :: initial => alpse_intial
        procedure :: set_ic => alpse_set_ic
        procedure :: setup => alpse_setup
        procedure :: solve => alpse_solve
        procedure :: finalize => alpse_eqn_finalize

    end type type_alpse_eqn

contains

    subroutine alpse_intial(this, jn_, Re, Ma, Te, yb_, ye_)

        implicit none
        class(type_alpse_eqn), intent(inout) :: this
        integer, intent(in) :: jn_
        real(real64), intent(in) ::  Re, Ma, Te
        integer, intent(in), optional :: yb_, ye_

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
        allocate (dA0(5, 5, jb:je)); dA0 = 0.0d0
        allocate (dB0(5, 5, jb:je)); dB0 = 0.0d0
        allocate (ddV0yy(5, 5, jb:je)); ddV0yy = 0.0d0
        allocate (dV0yy(5, 5, jb:je)); dV0yy = 0.0d0
        allocate (dV0xy(5, 5, jb:je)); dV0xy = 0.0d0
        allocate (dV0yz(5, 5, jb:je)); dV0yz = 0.0d0

        allocate (B_0(5, 5, jb:je)); B_0 = 0.0d0
        allocate (D_0(5, 5, jb:je)); D_0 = 0.0d0
        allocate (V_0(5, 5, jb:je)); V_0 = 0.0d0
        allocate (D_2(5, 5, jb:je)); D_2 = 0.0d0
        allocate (dA_0(5, 5, jb:je)); dA_0 = 0.0d0
        allocate (dB_0(5, 5, jb:je)); dB_0 = 0.0d0
        allocate (dV_0(5, 5, jb:je)); dV_0 = 0.0d0
        allocate (ddV_0(5, 5, jb:je)); ddV_0 = 0.0d0

        allocate (modify_Px_A(2, jb:je)); modify_Px_A = 0.0d0

        alpse_eqn_has_initialize = .TRUE.

        call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, this%dm, ierr)
        call DMSetUp(this%dm, ierr)
        call DMCreateMatrix(this%dm, L_ALPSE, ierr)
        call DMGetGlobalVector(this%dm, RHS, ierr); call VecZeroEntries(RHS, ierr)
        call VecDuplicate(RHS, this%phi, ierr); call VecZeroEntries(this%phi, ierr)

        allocate (this%phi0_array(5, jn)); this%phi0_array = 0.0d0
        allocate (This%phi1_array(5, jn)); this%phi1_array = 0.0d0
        allocate (this%phi2_array(5, jn)); this%phi2_array = 0.0d0

        call KSPCreate(PETSC_COMM_WORLD, solver, ierr)
        call KSPGetPC(Solver, SolverPC, ierr)
        call KSPSetType(Solver, KSPPREONLY, ierr)
        call KSPGetPC(Solver, SolverPC, ierr)
        call PCSetType(SolverPC, PCLU, ierr)
#ifdef __linux__
        call PCFactorSetMatSolverType(SolverPC, MATSOLVERMKL_CPARDISO, ierr)
#elif __APPLE__
        call PCFactorSetMatSolverType(SolverPC, MATSOLVERSUPERLU_DIST, ierr)
#endif

        call KSPSetFromOptions(Solver, ierr)
        ! call KSPView(Solver, PETSC_VIEWER_STDOUT_SELF, ierr)
        call PetscPrintf(PETSC_COMM_WORLD, 'allocate ALPSE!\n', ierr)

        lambda = 0.0d0 !< initial lambda is 1st order.

    end subroutine alpse_intial

    subroutine alpse_setup(this)

        implicit none
        class(type_alpse_eqn), intent(inout) :: this

        call solve_baseflow_matrix(this)

    end subroutine alpse_setup

    subroutine alpse_set_ic(this, phi)

        implicit none
        class(type_alpse_eqn), intent(inout) :: this
        Vec, intent(inout) ::  phi
        PetscScalar, pointer :: tmp(:, :)

        call VecCopy(phi, this%phi, ierr)
        call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
        this%phi0_array = tmp
        call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)
        this%phi1_array = 0.0d0
        this%phi2_array = 0.0d0
        lambda = 0.0d0

    end subroutine alpse_set_ic

    subroutine alpse_solve(this, alf, phi0, phi1, phi2)

        implicit none
        class(type_alpse_eqn), intent(inout) :: this
        complex(REAL64), intent(inout) :: alf
        Vec, intent(inout) ::  phi0
        Vec, intent(in), optional :: phi1, phi2
        PetscScalar, pointer :: tmp(:, :)
        complex(REAL64) :: alf_old, delta

        if (this%solve_alpha /= SOLVE_ALPHA_NONE) then
            print *, 'The initial alf is', alf
        end if

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
            lambda = 0.0d0
        else if (present(phi1) .and. present(phi2)) then
            call DMDAVecGetArrayReadF90(this%dm, phi1, tmp, ierr)
            this%phi1_array = tmp
            call DMDAVecRestoreArrayReadF90(this%dm, phi1, tmp, ierr)
            call DMDAVecGetArrayReadF90(this%dm, phi2, tmp, ierr)
            this%phi2_array = tmp
            call DMDAVecRestoreArrayReadF90(this%dm, phi2, tmp, ierr)
            lambda = 0.5d0
        end if

        if (sum(abs(this%phi2_array)) <= epsilon(1.0d0)) then  !phi2 is zero means 1 order is used
            lambda = 0.0d0
        else ! else 2 order is used
            lambda = 0.5d0
        end if

        alf_old = alf
        delta = 1.0d10
        ! read(*, *)

        block
            integer :: sum

            sum = 0
            do while (abs(delta) >= 1.0d-8)

                call solve_alpse_matrix(this, alf)
                call solve_left_matrix(this)
                call solve_rhs_new(this)
                call solve_alpse_matrix_bc(this, alf)
                call solve_alpse_bc(this)
                call solve_alpse_phi(this)
                call solve_alpha_local_new(this, alf)
                ! print*, alf, alf_old, (alf-alf_old)
                delta = abs(alf - alf_old)
                alf_old = alf
                sum = sum + 1
                ! call VecWAXPY(dphi, (-1.0d0, 0.0d0), phi_old, phi_new, ierr)
                ! call VecNorm(dphi, NORM_2, delta, ierr)
                ! call VecEqual(phi_new, phi_old, flg, ierr)
            end do
            if (this%solve_alpha /= SOLVE_ALPHA_NONE) print *, 'the iteration times is', sum
        end block
        call VecCopy(this%phi, phi0, ierr)
        call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
        this%phi0_array = tmp
        call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)

    end subroutine alpse_solve

    subroutine solve_alpha_local_new(this, alf)

        use functional
        use mod_dx
        implicit none
        class(type_alpse_eqn) :: this
        complex(REAL64), intent(inout) :: alf
        integer :: l, j
        real(REAL64), dimension(jn) :: hx
        complex(REAL64), dimension(5, jn) :: dphidi, dphidj, dphidjj, dphidx
        complex(REAL64), dimension(jn) :: phi_ke_j, dphi_ke_j
        PetscScalar :: phi_ke, dphi_ke
        PetscScalar, pointer :: tmp(:, :)

        if (this%solve_alpha /= SOLVE_ALPHA_NONE) then
            phi_ke = 0.0d0; dphi_ke = 0.0d0

            call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
            this%phi0_array = tmp
            call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)
            dphidi = this%phi1_array - this%phi0_array
    !!dphi_dj, dphi_djj
            do l = 1, 5
                call dx_fdm4nd_center(this%phi0_array(l, :), jb, je, dphidj(l, :), jb, je)
                !call dxx_fdm4nd_center(this%phi0_array(l, :), jb, je, dphidjj(l, :), jb, je)
            end do

            hx(:) = 1 + this%grid%drx*this%grid%y(:)
            do l = 1, 5
                dphidx(l, :) = dphidi(l, :)*this%grid%xi_x(:) + dphidj(l, :)*this%grid%eta_x(:)
            end do

            do j = 1, jn
                phi_ke_j(j) = (this%phi0_array(2, j)*conjg(this%phi0_array(2, j)) + &
                               this%phi0_array(3, j)*conjg(this%phi0_array(3, j)) + &
                               this%phi0_array(4, j)*conjg(this%phi0_array(4, j)))*hx(j)
                dphi_ke_j(j) = (dphidx(2, j)*conjg(this%phi0_array(2, j)) + &
                                dphidx(3, j)*conjg(this%phi0_array(3, j)) + &
                                dphidx(4, j)*conjg(this%phi0_array(4, j)))*hx(j)
            end do

            phi_ke = Trap(phi_ke_j, this%grid%y)
            dphi_ke = Trap(dphi_ke_j, this%grid%y)

            alf = alf - this%rd_factor*conjg(ZI*1.0d0/phi_ke*dphi_ke)
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
        class(type_alpse_eqn) :: this

        integer :: bs
        PetscScalar :: rhs_j(5, jb:je), box(5*(je - jb + 1))
        integer :: j
        PetscScalar, dimension(5, jb:je) :: phi_iminus1, phi_iminus2

  !! RHS=lambda*xi_x*A(phi(i-2)-phi(i-1))+(1+lambda)*xi_x*A*phi(i-1)

        do j = jb, je
            associate (A_0 => D_2)
                rhs_j(:, j) = lambda*matmul(A_0(:, :, j), (this%phi2_array(:, j) - this%phi1_array(:, j))) &
                              - (1.0d0 + lambda)*matmul(A_0(:, :, j), this%phi1_array(:, j))
            end associate
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

    end subroutine solve_rhs_new

    subroutine alpse_eqn_finalize(this)

        implicit none
        class(type_alpse_eqn), intent(inout) :: this

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
        deallocate (dA0)
        deallocate (dB0)
        deallocate (dV0yy)
        deallocate (ddV0yy)
        deallocate (dV0xy)
        deallocate (dV0yz)

        deallocate (B_0)
        deallocate (D_0)
        deallocate (V_0)
        deallocate (D_2)
        deallocate (dB_0)
        deallocate (dA_0)
        deallocate (dV_0)
        deallocate (ddV_0)
        deallocate (modify_Px_A)

        ! deallocate(phi_m1)
        ! deallocate(phi_m2)
        alpse_eqn_has_initialize = .FALSE.

        deallocate (this%phi0_array)
        deallocate (this%phi1_array)
        deallocate (this%phi2_array)

        ! call VecDestroy(phi_new, ierr)
        ! call VecDestroy(phi_old, ierr)
        call DMRestoreGlobalVector(this%dm, RHS, ierr)
        call DMRestoreGlobalVector(this%dm, this%phi, ierr)
        call MatDestroy(L_ALPSE, ierr)
        call KSPDestroy(solver, ierr)
        call DMDestroy(this%dm, ierr)

    end subroutine alpse_eqn_finalize

    subroutine solve_baseflow_matrix(this)

        implicit none
        class(type_alpse_eqn) :: this
        integer :: j, j_b, j_e

    !! internal point
        call solve_bf_matrix_j(this, jb, je)

    end subroutine solve_baseflow_matrix

    subroutine solve_bf_matrix_j(this, jb, je)

        implicit none
        class(type_alpse_eqn) :: this
        integer, intent(in) :: jb, je

        real(REAL64), parameter :: s = 0.0d0
        real(REAL64) :: u, ux, uy, uz, uyy, uxx, uxy, uxz
        real(REAL64) :: v, vx, vy, vz, vyy, vxy, vxz
        real(REAL64) :: w, wx, wy, wz, wxx, wyy, wxz, wyz
        real(REAL64) :: rho, rhox, rhoy, rhoz
        real(REAL64) :: T, Tx, Ty, Tz, Tyy, Txx, Tzz, Tyz

        real(REAL64) :: y, drx, drz
        real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132, d12y, d32y, d132y

        real(REAL64) :: divu, divuy, divux, divuz
        real(REAL64) :: p11, p12, p13, p22, p23, p33, p12y, p22y, p23y, phi
        real(REAL64) :: p11x, p12x, p13x, p13z, p23z, p33z
        real(REAL64) :: pe, cm, ck
        real(REAL64) :: rm, rmx, rmy, rmz, rmt, rmtt, rmtx, rmtz, rmty, rmyy, rmyz
        real(REAL64) :: rk, rkt, rktt, rky, rkz, rkty, rkyy
        real(REAL64) :: cv, cp, rek
        real(REAL64) :: omega_p, Mx

        integer :: j
        real(REAL64), parameter :: Z2D3 = 2.d0/3.d0, D3 = 1.d0/3.d0

        integer :: jb_local, je_local

        jb_local = jb !max(2, jb)
        je_local = je !min(je, jn-1)

        associate (grid_profile => this%grid, bf_profile => this%baseflow, &
                   Re => this%Re, Ma => this%Ma, Pr => this%Pr, Te => this%Te, &
                   C1 => this%C1, C2 => this%C2, GAMMA => this%GAMMA)
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
                d12y = -d12**2

                hz = 1.d0 + drz*y
                hz_y = drz
                d32 = hz_y/hz
                d32y = -d32**2

            case (CYLINDRICAL_CS)
                hx = 1.d0
                hx_y = 0.0d0
                d12 = 0.0d0
                d12y = -d12**2

                hz = y
                if(abs(hz)<= epsilon(1.0d0)) hz=epsilon(1.0d0)
                hz_y = 1.0d0
                d32 = hz_y/hz
                d32y = -d32**2
            end select

            d132 = d12 + d32
            d132y = d12y + d32y

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
            vxz = 0.0d0/hx/hz
            wxz = 0.0d0/hx/hz
            wyz = 0.0d0/hz
            Tyz = 0.0d0/hz

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
            rmyy = rmtt*ty**2 + rmt*tyy
            rmyz = rmtt*ty*tz + rmt*tyz

            rk = T*sqrt(T)*(1.d0 + ck)/(T + ck)/rek
            rkt = rk*(1.5d0/T - 1.0d0/(T + ck))
            rktt = rkt*(1.5d0/T - 1.0d0/(T + ck)) - rk*(1.5d0/T**2 - 1.0d0/(T + ck)**2)
            rky = rkt*ty
            rkz = 0.0d0
            rkty = rktt*ty
            rkyy = rktt*ty**2 + rkt*tyy

            divu = ux + vy + v*d132
            divux = uxx + vxy + vx*d132
            divuy = uxy + vyy + vy*d132 - v*(d12**2 + d32**2)
            divuz = uxz + vxz + vz*d132

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

            A0(3, 2, j) = rm*(-z2d3*(s - 1.0d0)*d132 + d12) - Z2D3*rmy*(s - 1.0d0)
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

            dA0(1, 1, j) = ux
            dA0(1, 2, j) = rhox
            dA0(2, 1, j) = pe*Tx
            dA0(2, 2, j) = rhox*u + rho*ux
            dA0(2, 3, j) = 0.0d0
            dA0(2, 4, j) = 0.0d0
            dA0(2, 5, j) = pe*rhox
            dA0(3, 2, j) = 0.0d0
            dA0(3, 3, j) = rhox*u + rho*ux
            dA0(3, 5, j) = 0.0d0
            dA0(4, 2, j) = 0.0d0
            dA0(4, 4, j) = rhox*u + rho*ux
            dA0(4, 5, j) = 0.0d0
            dA0(5, 1, j) = -pe*ux*T - pe*u*Tx
            dA0(5, 2, j) = 0.0d0
            dA0(5, 3, j) = 0.0d0
            dA0(5, 4, j) = 0.0d0
            dA0(5, 5, j) = rhox*cp*u - pe*ux*rho + rho*cp*ux - pe*u*rhox

            modify_Px_A(1, j) = pe*T*(-1.0d0 + omega_p)
            modify_Px_A(2, j) = pe*rho*(-1.0d0 + omega_p)

            B0(1, 1, j) = v
            B0(1, 3, j) = rho
            B0(2, 2, j) = rho*v - rmy - d32*rm
            B0(2, 3, j) = z2d3*(-rmx)*(s - 1)
            B0(2, 5, j) = -rmt*p12
            B0(3, 1, j) = pe*T
            B0(3, 2, j) = -rmx
            B0(3, 3, j) = rho*v - Z2D3*(2*s + 1)*rm*d132 - (2.0d0 + s)*Z2d3*rmy
            B0(3, 4, j) = -rmz
            B0(3, 5, j) = pe*rho - rmt*p22
            B0(4, 3, j) = z2d3*(-rmz)*(s - 1)
            B0(4, 4, j) = rho*v - rm*d12 - rmy
            B0(4, 5, j) = -rmt*p23
            B0(5, 1, j) = -pe*v*t
            B0(5, 2, j) = -2.0d0*rm*p12
            B0(5, 3, j) = -2.0d0*rm*p22
            B0(5, 4, j) = -2.0d0*rm*p23
            B0(5, 5, j) = rho*cp*v - pe*v*rho - 2.0d0*rky - rk*d132

            dB0(1, 1, j) = vy
            dB0(1, 3, j) = rhoy
            dB0(2, 2, j) = rhoy*v + rho*vy - rmyy - d32y*rm - d32*rmy
            dB0(2, 3, j) = 0.0d0
            dB0(2, 5, j) = -rmty*p12 - rmt*p12y
            dB0(3, 1, j) = pe*Ty
            dB0(3, 2, j) = 0.0d0
            dB0(3, 3, j) = rhoy*v + rho*vy - Z2D3*(2*s + 1)*(rmy*d132 + rm*d132y) - (2.0d0 + s)*Z2d3*rmyy
            dB0(3, 4, j) = -rmyz
            dB0(3, 5, j) = pe*rhoy - rmty*p22 - rmt*p22y
            dB0(4, 3, j) = z2d3*(-rmyz)*(s - 1)
            dB0(4, 4, j) = rhoy*v + rho*vy - rmy*d12 - rm*d12y - rmyy
            dB0(4, 5, j) = -rmty*p23 - rmt*p23y
            dB0(5, 1, j) = -pe*(vy*t + v*ty)
            dB0(5, 2, j) = -2.0d0*(rmy*p12 + rm*p12y)
            dB0(5, 3, j) = -2.0d0*(rmy*p22 + rm*p22y)
            dB0(5, 4, j) = -2.0d0*(rmy*p23 + rm*p23y)
            dB0(5, 5, j) = rhoy*cp*v + rho*cp*vy - pe*(vy*rho + v*rhoy) - 2.0d0*rkyy - (rky*d132 + rk*d132y)

            C0(1, 1, j) = w
            C0(1, 4, j) = rho
            C0(2, 2, j) = rho*w - rmz
            C0(2, 4, j) = -z2d3*rmx*(s - 1.0d0)
            C0(2, 5, j) = -rmt*p13
            C0(3, 3, j) = rho*w - rmz
            !C0(3, 4, j)=7.0d0/3.0d0*rm*d32+Z2D3*rmy
            C0(3, 4, j) = (-z2d3*(s - 1.0d0)*d132 + d32)*rm - z2d3*rmy*(s - 1.0d0)
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
            D0(3, 3, j) = rho*vy - z2d3*(s - 1)*d132*(d132*rm + rmy)
            D0(3, 4, j) = rho*(vz - 2.0d0*w*d32) + rmz*d32
            D0(3, 5, j) = pe*rhoy - rmtx*p12 - rmty*p22 - rmtz*p23 - rmt*(p12x + p22y + p23z + p22*d132)
            D0(4, 1, j) = Tz*pe + u*wx + v*wy + w*(wz + v*d32)
            D0(4, 2, j) = rho*wx
            D0(4, 3, j) = rho*(wy + w*d32) - z2d3*((s - 1)*d132 + 3*d32)*rmz
            D0(4, 4, j) = rho*(wz + d32*v) + rm*d132*d32 + rmy*d32
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

            dV0yy(2, 2, j) = rmy
            dV0yy(3, 3, j) = (2.0 + s)*Z2D3*rmy
            dV0yy(4, 4, j) = rmy
            dV0yy(5, 5, j) = rky

            ddV0yy(2, 2, j) = rmyy
            ddV0yy(3, 3, j) = (2.0 + s)*Z2D3*rmyy
            ddV0yy(4, 4, j) = rmyy
            ddV0yy(5, 5, j) = rkyy

            V0zz(2, 2, j) = rm
            V0zz(3, 3, j) = rm
            V0zz(4, 4, j) = (2.0 + s)*Z2D3*rm
            V0zz(5, 5, j) = rk

            V0xy(2, 3, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0xy(3, 2, j) = D3*(2.0d0*s + 1.0d0)*rm
            dV0xy(2, 3, j) = D3*(2.0d0*s + 1.0d0)*rmy
            dV0xy(3, 2, j) = D3*(2.0d0*s + 1.0d0)*rmy

            V0xz(2, 4, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0xz(4, 2, j) = D3*(2.0d0*s + 1.0d0)*rm

            V0yz(3, 4, j) = D3*(2.0d0*s + 1.0d0)*rm
            V0yz(4, 3, j) = D3*(2.0d0*s + 1.0d0)*rm
            dV0yz(3, 4, j) = D3*(2.0d0*s + 1.0d0)*rmy
            dV0yz(4, 3, j) = D3*(2.0d0*s + 1.0d0)*rmy

            ! call get_Aplus(A0(:, :, j), G0(:, :, j), Aplus(:, :, j), j)
            ! Aplus(:, :, j)=Aplus(:, :, j)/hx
            A0(:, :, j) = A0(:, :, j)/hx
            dA0(:, :, j) = dA0(:, :, j)/hx
            C0(:, :, j) = C0(:, :, j)/hz
            V0xx(:, :, j) = V0xx(:, :, j)/hx**2
            V0xy(:, :, j) = V0xy(:, :, j)/hx
            dV0xy(:, :, j) = dV0xy(:, :, j)/hx
            V0xz(:, :, j) = V0xz(:, :, j)/hx/hz
            V0yz(:, :, j) = V0yz(:, :, j)/hz
            dV0yz(:, :, j) = dV0yz(:, :, j)/hz
            V0zz(:, :, j) = V0zz(:, :, j)/hz**2
            modify_Px_A(:, j) = modify_Px_A(:, j)/hx

        end do
        end associate

    end subroutine solve_bf_matrix_j

    subroutine solve_alpse_matrix(this, alf)
        implicit none
        class(type_alpse_eqn) :: this
        complex(REAL64), intent(in) :: alf

        integer :: j

        do j = jb, je
            call solve_alpse_matrix_j(this, j, alf)
        end do

        ! do  j=jb, je
        !   call solve_alpse_matrix_j(j)
        ! enddo

    end subroutine solve_alpse_matrix

    subroutine solve_alpse_matrix_j(this, j, alf)

        implicit none
        class(type_alpse_eqn) :: this
        integer, intent(in) :: j
        complex(REAL64), intent(in) :: alf

        complex(REAL64), dimension(5, 5) :: A_, B_, V_, dAx_

    !!d\alpha/dx
        ! dadx=(alf-wave_param(2)%alf)*grid_profile%xi_x(jb)

        ! if(j==121)then
        !   print*, alf, wave_param(2)%alf, dadx
        !   read(*, *)
        ! endif

        associate (grid_profile => this%grid, ome => this%ome, bat => this%bat)
            associate (A_0 => D_2, xi_x => grid_profile%xi_x(j), &
                       xi_y => grid_profile%xi_y(j), &
                       eta_x => grid_profile%eta_x(j), &
                       eta_y => grid_profile%eta_y(j), &
                       xi_yy => grid_profile%xi_yy(j), &
                       eta_yy => grid_profile%eta_yy(j))

                D_0(:, :, j)=D0(:, :, j)-ZI*ome*G0(:, :, j)+ZI*alf*A0(:, :, j) &
                             +ZI*bat*C0(:, :, j)+alf**2*V0xx(:, :, j) &
                             +alf*bat*V0xz(:, :, j)+bat**2*V0zz(:, :, j) !&
                !  -ZI*dadx*V0xx(:, :, j)

!          D_0(2, 1, j)=D_0(2, 1, j)+modify_Px_D(1, j)
!          D_0(2, 5, j)=D_0(2, 5, j)+modify_Px_D(2, j)
                ! A_0(:, :, j)=Aplus(:, :, j)-2.0d0*ZI*alf*V0xx(:, :, j)-ZI*bat*V0xz(:, :, j)
                A_0(:, :, j) = A0(:, :, j) - 2.0d0*ZI*alf*V0xx(:, :, j) - ZI*bat*V0xz(:, :, j)
                A_0(2, 1, j) = A_0(2, 1, j) + modify_Px_A(1, j)
                A_0(2, 5, j) = A_0(2, 5, j) + modify_Px_A(2, j)
                dA_0(:, :, j) = dA0(:, :, j)

                B_0(:, :, j) = B0(:, :, j) - ZI*alf*V0xy(:, :, j) - ZI*bat*V0yz(:, :, j)
                dB_0(:, :, j) = dB0(:, :, j) - ZI*alf*dV0xy(:, :, j) - ZI*bat*dV0yz(:, :, j)
                V_0(:, :, j) = -V0yy(:, :, j)
                dV_0(:, :, j) = -dV0yy(:, :, j)
                ddV_0(:, :, j) = -ddV0yy(:, :, j)

                A_0(:, :, j) = -A_0(:, :, j)
                B_0(:, :, j) = -B_0(:, :, j) + 2.0d0*dV_0(:, :, j)

                D_0(:, :, j) = D_0(:, :, j) - dA_0(:, :, j) - dB_0(:, :, j) + ddV_0(:, :, j)

                D_0(:, :, j) = conjg(transpose(D_0(:, :, j)))
                B_0(:, :, j) = conjg(transpose(B_0(:, :, j)))
                V_0(:, :, j) = conjg(transpose(V_0(:, :, j)))
                A_0(:, :, j) = conjg(transpose(A_0(:, :, j)))

                A_ = A_0(:, :, j)*xi_x + B_0(:, :, j)*xi_y + V_0(:, :, j)*xi_yy
                B_ = A_0(:, :, j)*eta_x + B_0(:, :, j)*eta_y + V_0(:, :, j)*eta_yy
                V_ = V_0(:, :, j)*eta_y**2

                A_0(:, :, j) = A_
                B_0(:, :, j) = B_
                V_0(:, :, j) = V_

            end associate
        end associate

    end subroutine solve_alpse_matrix_j

    subroutine solve_left_matrix(this)

        implicit none
        type(type_alpse_eqn), intent(in) :: this
        real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

        PetscScalar :: box1(5, -2:2, 5)
        complex(REAL64) :: tmp_mat_j0(5, 5, -2:2)
        ! real(REAL64) :: eta_y, eta_yy

        integer :: i, j, lj
        integer :: ljb, lje, c_idx
        MatStencil :: idxm(4, 1), idxn(4, -2:2)
        PetscViewer :: Viewer

        integer :: jb_local, je_local

        associate (grid_profile => this%grid)
            jb_local = jb !max(2, jb)
            je_local = je !min(je, jn-1)

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

                associate (A_0 => D_2, C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
                    tmp_mat_j0 = 0.0d0
                    box1 = 0.0d0
                    do lj = ljb, lje
                        tmp_mat_j0(:, :, lj) = delta(lj)*(-(1.0d0 + lambda)*A_0(:, :, j) + D_0(:, :, j)) + &
                                               C1(lj, c_idx)*(B_0(:, :, j)) + &
                                               D1(lj, c_idx)*(V_0(:, :, j))

                        box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
                    end do
                end associate

                call MatSetValuesBlockedStencil(L_ALPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
            end do

            call MatAssemblyBegin(L_ALPSE, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(L_ALPSE, MAT_FINAL_ASSEMBLY, ierr)

            ! if(screen_flg)then
            !   call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LPSE left matrix done!\n", ierr)
            ! endif
        end associate

    end subroutine solve_left_matrix

    subroutine solve_alpse_matrix_bc(this, alf)

        implicit none
        class(type_alpse_eqn) :: this
        complex(REAL64), intent(in) :: alf

    !! b.c. point
        if (jb == 1) call solve_alpse_matrix_j0(this, jb, alf)
        if (je == jn) call solve_alpse_matrix_jn(this, je, alf)

    end subroutine solve_alpse_matrix_bc

    subroutine solve_alpse_matrix_j0(this, jb, alf)

        implicit none
        integer, intent(in) :: jb
        class(type_alpse_eqn) :: this
        complex(REAL64), intent(in) :: alf

        real(REAL64) :: rho, rhox, rhoy
        real(REAL64) :: u, ux
        real(REAL64) :: v, vy
        real(REAL64) :: w
        real(REAL64) :: divu
        real(REAL64) :: y, drx, drz
        real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

        associate (grid_profile => this%grid, bf_profile => this%baseflow)
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
                if(abs(hz)<= epsilon(1.0d0)) hz = epsilon(1.0d0)
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

            select case (this%bc_y0)
            case (BC_WALL)
                V0xx(2:5, :, 1) = 0.0d0
                V0xy(2:5, :, 1) = 0.0d0
                V0xz(2:5, :, 1) = 0.0d0
                V0yy(2:5, :, 1) = 0.0d0
                V0yz(2:5, :, 1) = 0.0d0
                V0zz(2:5, :, 1) = 0.0d0
                D_2(2:5, :, 1) = 0.0d0
                B_0(2:5, :, 1) = 0.0d0
                V_0(2:5, :, 1) = 0.0d0
                D_0(2:5, :, 1) = 0.0d0
                D_0(2, 2, 1) = 1.0d0
                D_0(3, 3, 1) = 1.0d0
                D_0(4, 4, 1) = 1.0d0
                if(this%isothermal)then
                    D_0(5, 5, 1) = 1.0d0
                else
                    B_0(5, 5, 1) = 1.0d0
                endif
            case (BC_POLE)
                D_2(:, :, 1) = 0.0d0
                B_0(:, :, 1) = 0.0d0
                V_0(:, :, 1) = 0.0d0
                D_0(:, :, 1) = 0.0d0
                V0xx(:, :, 1) = 0.0d0
                V0xy(:, :, 1) = 0.0d0
                V0xz(:, :, 1) = 0.0d0
                V0yy(:, :, 1) = 0.0d0
                V0yz(:, :, 1) = 0.0d0
                V0zz(:, :, 1) = 0.0d0
                if (int(this%bat + 0.001d0) == 0) then
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
            end select

        end associate

        ! A0(:, :, 1)=A0(:, :, 1)/hx
        ! C0(:, :, 1)=C0(:, :, 1)/hz

        ! call solve_lpse_matrix_j(bf_profile,grid_profile, 1, ome, alf, bat)

    end subroutine solve_alpse_matrix_j0

    subroutine solve_alpse_matrix_jn(this, je, alf)

        implicit none
        integer, intent(in) :: je
        class(type_alpse_eqn) :: this
        complex(REAL64), intent(in) :: alf

        real(REAL64) :: rho, u, w, rhoy
        real(REAL64) :: y, drx, drz
        real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132

        select case (this%bc_yn)
        case (BC_FARFIELD, BC_NONREFLECTION)
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
            D_0(3, 3, jn) = 1.0d0
            D_0(4, 4, jn) = 1.0d0
            D_0(5, 5, jn) = 1.0d0
        end select

    end subroutine solve_alpse_matrix_jn

    subroutine solve_alpse_bc(this)

        implicit none
        class(type_alpse_eqn) :: this

        integer :: j

        do j = jb, je
            if (j == 1 .or. j == jn) then
                call solve_alpse_bc_j(this, j)
            end if
        end do

    end subroutine solve_alpse_bc

    subroutine solve_alpse_bc_j(this, j)

        use functional
        implicit none
        class(type_alpse_eqn) :: this
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

        associate (A_0 => D_2, C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
            tmp_mat_j0 = 0.0d0
            box1 = 0.0d0
            do lj = ljb, lje
                tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 + lambda)*A_0(:, :, j) + D_0(:, :, j)) + &
                                       C1(lj, c_idx)*(B_0(:, :, j)) + &
                                       D1(lj, c_idx)*(V_0(:, :, j))
                box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
            end do
            ! block
            !   integer :: l1, l2
            !
            !   if(j==jn)then
            !   print*, j
            !   print*, ljb, lje
            !   do l1=1, 5
            !     do l2=1, 5
            !       print*, l1, l2, tmp_mat_j0(l1, l2, 0)
            !     enddo
            !   enddo
            !   stop
            !   endif
            !
            ! endblock

        end associate

        call MatSetValuesBlockedStencil(L_ALPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)

    !!j==1
        if (j == 1) then
            tmpj0 = 0.0d0
            select case (this%bc_y0)
            case (BC_WALL)
                call VecSetValues(RHS, 4, arange(1, 4), tmpj0(1:4), INSERT_VALUES, ierr)            
            case (BC_POLE)
                call VecSetValues(RHS, 5, arange(0, 4), tmpj0(1:5), INSERT_VALUES, ierr)                        
            case default
                call VecSetValues(RHS, 5, arange(0, 4), tmpj0(1:5), INSERT_VALUES, ierr)                                    
            end select
        elseif (j == jn) then
    !!j==jn
            tmpjn = 0.0d0
            call VecSetValues(RHS, 5, arange((jn - 1)*5 + 0, (jn - 1)*5 + 4), tmpjn, INSERT_VALUES, ierr)
        end if

        call MatAssemblyBegin(L_ALPSE, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(L_ALPSE, MAT_FINAL_ASSEMBLY, ierr)
        call VecAssemblyBegin(RHS, ierr)
        call VecAssemblyEnd(RHS, ierr)

        ! if(screen_flg)then
        !   call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LPSE matrix and RHS BC done!\n", ierr)
        ! endif
        ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)

    end subroutine solve_alpse_bc_j

    subroutine solve_alpse_phi(this)

        implicit none
        class(type_alpse_eqn) :: this
        integer :: its
        character(20) :: ksptype

        call KSPSetOperators(Solver, L_ALPSE, L_ALPSE, ierr)
        call KSPSetFromOptions(solver, ierr)
        call KSPSetup(Solver, ierr)
        call KSPSolve(Solver, RHS, this%phi, ierr)

    end subroutine solve_alpse_phi

    subroutine get_Aplus(A, G, Aplus, j)

        implicit none
        real(REAL64), dimension(5, 5), intent(in) :: A, G
        real(REAL64), dimension(5, 5), intent(out) :: Aplus
        integer :: j
        integer :: ipiv(5), info
        real(REAL64), dimension(5, 5) :: At, Gt
        real(REAL64) :: work(100)
        real(REAL64) :: alfr(5), alfi(5), beta(5)
        real(REAL64) :: vl(5, 5), vr(5, 5)
        complex(REAL64) :: lambda_(5)
        real(REAL64) :: diag(5, 5)
        integer :: l1, l2
        integer :: i
        real(REAL64) :: maxA

        Aplus = 0.0d0; diag(5, 5) = 0.0d0
        work = 0.0d0
        alfr = 0.0d0; alfi = 0.0d0; beta = 0.0d0
        vl = 0.0d0; vr = 0.0d0; lambda_ = 0.0d0; diag = 0.0d0
        At = A; Gt = G

        call dggev('N', 'V', 5, At, 5, Gt, 5, alfr, alfi, &
                   beta, vl, 5, vr, 5, work, 100, info)
        lambda_ = (alfr + ZI*alfi)/beta
        do i = 1, 5
            diag(i, i) = max(0.0d0, lambda_(i)%re)
        end do
        vl = vr
    !! vl is inv(vr)
        call dgetrf(5, 5, vl, 5, ipiv, info)
        call dgetri(5, vl, 5, ipiv, work, 5, info)
        At = matmul(G, vr)
        At = matmul(At, diag)
        At = matmul(At, vl)

        Aplus = At

    end subroutine get_Aplus

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

end module mod_alpse_eqn

! subroutine initialize_alpse_eqn(grid_profile, baseflow_profile, ome, alf, bat)

!     implicit none

!     integer :: i, j
!     integer :: iounit, ios
!     PetscViewer :: Viewer
!     type(type_lst_grid_profile), intent(in) :: grid_profile
!     type(type_lst_bf_profile), intent(in) :: baseflow_profile
!     complex(REAL64), intent(in) :: ome, alf, bat

!     jb=yb+1; je=ye+1

!     allocate(G0  (5, 5, jb:je));G0  =0.0d0
!     allocate(A0  (5, 5, jb:je));A0  =0.0d0
!     allocate(B0  (5, 5, jb:je));B0  =0.0d0
!     allocate(C0  (5, 5, jb:je));C0  =0.0d0
!     allocate(D0  (5, 5, jb:je));D0  =0.0d0
!     allocate(V0xx(5, 5, jb:je));V0xx=0.0d0
!     allocate(V0xy(5, 5, jb:je));V0xy=0.0d0
!     allocate(V0xz(5, 5, jb:je));V0xz=0.0d0
!     allocate(V0yy(5, 5, jb:je));V0yy=0.0d0
!     allocate(V0yz(5, 5, jb:je));V0yz=0.0d0
!     allocate(V0zz(5, 5, jb:je));V0zz=0.0d0
!     allocate(Aplus(5, 5, jb:je)); Aplus=0.0d0
!     allocate(dA0(5, 5, jb:je)); dA0=0.0d0
!     allocate(dB0(5, 5, jb:je)); dB0=0.0d0
!     allocate(ddV0yy(5, 5, jb:je)); ddV0yy=0.0d0
!     allocate(dV0yy(5, 5, jb:je)); dV0yy=0.0d0
!     allocate(dV0xy(5, 5, jb:je)); dV0xy=0.0d0
!     allocate(dV0yz(5, 5, jb:je)); dV0yz=0.0d0

!     allocate(B_0(5, 5, jb:je));B_0=0.0d0
!     allocate(D_0(5, 5, jb:je));D_0=0.0d0
!     allocate(V_0(5, 5, jb:je));V_0=0.0d0
!     allocate(D_2(5, 5, jb:je));D_2=0.0d0
!     allocate(dA_0(5, 5, jb:je)) ; dA_0=0.0d0
!     allocate(dB_0(5, 5, jb:je)) ; dB_0=0.0d0
!     allocate(dV_0(5, 5, jb:je)) ; dV_0=0.0d0
!     allocate(ddV_0(5, 5, jb:je)) ; ddV_0=0.0d0

!     allocate(modify_Px_A(2, jb:je)); modify_Px_A=0.0d0

!     alpse_eqn_has_initialize=.TRUE.

!     call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_ALPSE, ierr)

!     call DMSetUp(DM_ALPSE, ierr)
!     call DMCreateMatrix(DM_ALPSE, L_ALPSE, ierr)
!     call DMGetGlobalVector(DM_ALPSE, RHS, ierr)

!     call VecDuplicate(RHS, phi_new, ierr)
!     call VecDuplicate(RHS, phi_old, ierr)
!     do i=1, 3
!       call VecDuplicate(RHS, phi(i), ierr)
!     enddo

!     allocate(phi_m1(5, jb:je))
!     allocate(phi_m2(5, jb:je))
!     phi_m1=0.0d0
!     phi_m2=0.0d0

!     call VecZeroEntries(RHS, ierr)
!     call VecZeroEntries(phi_old, ierr)
!     call VecZeroEntries(phi_new, ierr)

!     call KSPCreate(PETSC_COMM_WORLD, solver, ierr)
!     ! call KSPSetInitialGuessNonzero(Solver, PETSC_TRUE, ierr)
!     ! call KSPSetTolerances(Solver, 1.0d-16, 1.0d-16, 1.0d8, 20000, ierr)
!     ! call KSPSetType(Solver, KSPFGMRES, ierr)
!     ! call KSPGetPC(Solver, SolverPC, ierr)
!     ! call PCSetType(SolverPC, PCKSP, ierr)
!     call KSPSetType(Solver, KSPPREONLY, ierr)
!     call KSPGetPC(Solver, SolverPC, ierr)
!     call PCSetType(SolverPC, PCLU, ierr)
! #ifdef __linux__
!     call PCFactorSetMatSolverType(SolverPC, MATSOLVERMKL_CPARDISO, ierr)
! #elif __APPLE__
!     call PCFactorSetMatSolverType(SolverPC, MATSOLVERSUPERLU_DIST, ierr)
! #endif

!     call PetscPrintf(PETSC_COMM_WORLD, 'allocate ALPSE!\n', ierr)

!     lambda=0.0d0 !< initial lambda is 1st order.

!     call execute_command_line("mkdir -p out")
!     ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "out/disturbance_lpse.petsc", FILE_MODE_WRITE, Viewer, ierr)
!     ! call VecView(phi(2), Viewer, ierr)
!     ! call PetscViewerDestroy(Viewer, ierr)
!     ! print*, wave_param(1)%alf%re, wave_param(1)%alf%im
!     ! print*, wave_param(2)%alf%re, wave_param(2)%alf%im
!     wave_param(2)%alf=conjg(alf)
!     wave_param(2)%ome=ome
!     wave_param(2)%bat=bat
!     wave_param(3)=wave_param(2)

!     ! call solve_baseflow_matrix(baseflow_profile, grid_profile)
!     ! call solve_alpse_matrix(baseflow_profile, grid_profile, ome, alf, bat)

!   end subroutine initialize_alpse_eqn

!   subroutine solve_alpse_eqn(grid_profile, baseflow_profile, ome, alf, bat, ix, phi0, phi1, phi2)

!     implicit none
!     type(type_lst_grid_profile), intent(in) :: grid_profile
!     type(type_lst_bf_profile), intent(in) :: baseflow_profile
!     complex(REAL64), intent(inout) :: ome, alf, bat
!     Vec, intent(inout) :: phi0
!     Vec, intent(in) :: phi1, phi2
!     integer, intent(in) :: ix
!     PetscViewer :: Viewer
!     REAL(REAL64) :: delta
!     Vec :: dphi
!     logical :: flg

!     if(solve_alpha/=SOLVE_ALPHA_NONE)then
!       print*, 'start solve alpse at i=', ix
!       print*, 'The initial alf is', alf
!     endif
!     iloc=ix

!     flg=.TRUE.

!     ome1=ome; alf1=alf

!     wave_param(3)=wave_param(2) !<initial wave parameters
!     wave_param(3)%alf=alf
!     wave_param(3)%bat=bat

!     ! print*, wave_param(3)%alf, wave_param(3)%ome, wave_param(3)%bat

!     call VecCopy(phi0, phi(3), ierr)
!     call VecCopy(phi1, phi(2), ierr)
!     call VecCopy(phi2, phi(1), ierr)

!     call solve_baseflow_matrix(baseflow_profile, grid_profile)
!     call get_phi_before_iloc(ix)

!     delta=1.0d10; alf_old=0.0d0

!     do while(abs(delta)>=1.0d-8)
!       call VecCopy(phi_new, phi_old, ierr)

!       call solve_alpse_matrix(baseflow_profile, grid_profile, ome, alf, bat)
!       call solve_left_matrix(grid_profile)
!       call solve_rhs(grid_profile)
!       call solve_alpse_matrix_bc(baseflow_profile, grid_profile, ome, alf, bat)
!       call solve_alpse_bc(grid_profile)
!       call solve_alpse_phi(phi_new)
!       call solve_alpha_local(alf, grid_profile)

!       delta=abs(alf-alf_old)
!       if(screen_flg) then
!         print*, alf, alf_old, delta
!         read(*, *)
!       endif
!     enddo

!     call VecCopy(phi_new, phi(3), ierr)
!     call VecCopy(phi(3), phi0, ierr)
!     wave_param(3)%alf=alf
!     lambda=0.5d0

!     ! associate(A_0=>D_2)
!     !   if(ix==iend)then
!     !     open(999,file='in/matrix_alpse_bc0.mat', form='unformatted')
!     !     write(999)A_0(:, 1, 1), B_0(:, 1, 1), D_0(:, 1, 1), V_0(:, 1, 1)
!     !     close(999)
!     !   endif
!     ! endassociate

!   end subroutine solve_alpse_eqn

!   subroutine finalize_alpse_eqn()

!     implicit none
!     integer :: i

!     deallocate(G0  )
!     deallocate(A0  )
!     deallocate(B0  )
!     deallocate(C0  )
!     deallocate(D0  )
!     deallocate(V0xx)
!     deallocate(V0xy)
!     deallocate(V0xz)
!     deallocate(V0yy)
!     deallocate(V0yz)
!     deallocate(V0zz)
!     deallocate(Aplus)
!     deallocate(dA0)
!     deallocate(dB0)
!     deallocate(dV0yy)
!     deallocate(ddV0yy)
!     deallocate(dV0xy)
!     deallocate(dV0yz)

!     deallocate(B_0)
!     deallocate(D_0)
!     deallocate(V_0)
!     deallocate(D_2)
!     deallocate(dB_0)
!     deallocate(dA_0)
!     deallocate(dV_0)
!     deallocate(ddV_0)
!     ! deallocate(A_0_backward)

!     deallocate(modify_Px_A)

!     deallocate(phi_m1)
!     deallocate(phi_m2)
!     alpse_eqn_has_initialize=.FALSE.

!     call VecDestroy(phi_new, ierr)
!     call VecDestroy(phi_old, ierr)
!     call VecDestroy(RHS, ierr)
!     call KSPDestroy(solver, ierr)
!     call DMDestroy(DM_ALPSE, ierr)

!   end subroutine finalize_alpse_eqn

!   subroutine get_phi_before_iloc(ix)

!       implicit none
!       integer, intent(in) :: ix
!       PetscScalar, pointer :: tmp(:, :)

!       call DMDAVecGetArrayReadF90(DM_ALPSE, phi(2), tmp, ierr)
!       phi_m1=tmp
!       call DMDAVecRestoreArrayReadF90(DM_ALPSE, phi(2), tmp, ierr)

!       if(ix==iend-1)then
!         phi_m2=0.0d0
!       else
!         call DMDAVecGetArrayReadF90(DM_ALPSE, phi(1), tmp, ierr)
!         phi_m2=tmp
!         call DMDAVecRestoreArrayReadF90(DM_ALPSE, phi(1), tmp, ierr)
!       endif

!   end subroutine get_phi_before_iloc

! subroutine solve_rhs(grid_profile)

!   use functional
!   implicit none
!   type(type_lst_grid_profile), intent(in) :: grid_profile

!   integer :: bs
!   PetscScalar :: rhs_j(5,jb:je),  box(5*(je-jb+1))
!   integer :: j
!   PetscScalar, dimension(5, jb:je) :: phi_iminus1, phi_iminus2

!   !! RHS=lambda*xi_x*A(phi(i-2)-phi(i-1))+(1+lambda)*xi_x*A*phi(i-1)

!   do j=jb, je
!     associate(A_0=>D_2)
!       rhs_j(:, j)=lambda*matmul(A_0(:, :, j), (phi_m2(:, j)-phi_m1(:, j))) &
!            -(1.0d0+lambda)*matmul(A_0(:, :, j), phi_m1(:, j))
!     end associate
!   enddo
!   box=reshape(rhs_j, [5*(je-jb+1)])
!   call VecSetValuesBlocked(RHS, je-jb+1, arange(jb-1, je-1), box, INSERT_VALUES, ierr)
!   call VecAssemblyBegin(RHS, ierr)
!   call VecAssemblyEnd(RHS, ierr)
!   call VecGetBlockSize(RHS, bs, ierr)
!   ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)

! end subroutine solve_rhs

! subroutine solve_alpha_local(alf, grid_profile)

!   use functional
!   use global_parameters, only: solve_alpha
!   use global_parameters, only: SOLVE_ALPHA_BOTH, SOLVE_ALPHA_NONE, SOLVE_ALPHA_REAL
!   use mod_dx
!   implicit none
!   complex(REAL64), intent(inout) :: alf
!   type(type_lst_grid_profile), intent(in) :: grid_profile

!   Vec ::  phi_new_local, phi_i_1_local
!   integer :: l, j

!   PetscScalar, pointer :: tmp(:, :)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: phi_i0 !< phi(i)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: phi_i1 !< phi(i-1)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: dphidi !< phi(i-1)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: dphidj !< phi(i-1)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: dphidjj !< phi(i-1)
!   PetscScalar, dimension(5, ybg+1:yeg+1) :: dphidx !< phi(i-1)
!   Real(REAL64), dimension(ybg+1:yeg+1) :: xi_x_local, eta_x_local
!   PetscScalar, dimension(ybg+1:yeg+1) :: phi_ke_j, dphi_ke_j
!   real(REAL64), dimension(ybg+1:yeg+1) :: hx
!   PetscScalar :: phi_ke, dphi_ke, phi_ke_total, dphi_ke_total
!   PetscScalar, pointer :: grid(:, :, :)

!   if(solve_alpha/=SOLVE_ALPHA_NONE)then
!   ! xi_x_local=0.0d0; eta_x_local=0.0d0
!     phi_ke_total=0.0d0; dphi_ke_total=0.0d0

!     call DMGetLocalVector(DM_ALPSE, phi_new_local, ierr)
!     call DMGetLocalVector(DM_ALPSE, phi_i_1_local, ierr)

!     call DMGlobalToLocalBegin(DM_ALPSE, phi_new, INSERT_VALUES, phi_new_local, ierr)
!     call DMGlobalToLocalEnd(DM_ALPSE, phi_new, INSERT_VALUES, phi_new_local, ierr)
!     ! call DMGlobalToLocalBegin(DM_ALPSE, phi(3), INSERT_VALUES, phi_new_local, ierr)
!     ! call DMGlobalToLocalEnd(DM_ALPSE, phi(3), INSERT_VALUES, phi_new_local, ierr)

!     call DMGlobalToLocalBegin(DM_ALPSE, phi(2), INSERT_VALUES, phi_i_1_local, ierr)
!     call DMGlobalToLocalEnd(DM_ALPSE, phi(2), INSERT_VALUES, phi_i_1_local, ierr)
!     call DMDAVecGetArrayReadF90(DM_ALPSE, phi_new_local, tmp, ierr)
!     phi_i0=tmp
!     call DMDAVecRestoreArrayReadF90(DM_ALPSE, phi_new_local, tmp, ierr)
!     call DMDAVecGetArrayReadF90(DM_ALPSE, phi_i_1_local, tmp, ierr)
!     phi_i1=tmp
!     call DMDAVecRestoreArrayReadF90(DM_ALPSE, phi_i_1_local, tmp, ierr)
!     call DMRestoreLocalVector(DM_ALPSE, phi_new_local, ierr)
!     call DMRestoreLocalVector(DM_ALPSE, phi_i_1_local, ierr)

!     !!dphi_di
!     dphidi=phi_i1-phi_i0
!     !!dphi_dj
!     do l=1, 5
!       call dx_fdm4nd_center(phi_i0(l, :), ybg+1, yeg+1, dphidj(l, :), ybg+1, yeg+1)
!       call dxx_fdm4nd_center(phi_i0(l, :), ybg+1, yeg+1, dphidjj(l, :), ybg+1, yeg+1)

!     enddo

!     xi_x_local(:)=grid_profile%xi_x
!     eta_x_local(:)=grid_profile%eta_x
!     hx(:)=1+grid_profile%drx*grid_profile%y(:)
!     ! do j=ybg+1, yb+1
!     !   xi_x_local(ybg+1:yb+1)=xi_x_local(yb+1)
!     !   eta_x_local(ybg+1:yb+1)=eta_x_local(yb+1)
!     !   hz(ybg+1:yb+1)=hz(yb+1)
!     ! enddo
!     ! do j=ye+1, yeg+1
!     !   xi_x_local(ye+1:yeg+1)=xi_x_local(ye+1)
!     !   eta_x_local(ye+1:yeg+1)=eta_x_local(ye+1)
!     !   hz(ye+1:yeg+1)=hz(ye+1)
!     ! enddo

!     do l=1, 5
!       dphidx(l, :)=dphidi(l, :)*xi_x_local+dphidj(l, :)*eta_x_local
!     enddo

!     ! block
!     !   integer :: jloc
!     !   integer :: l1, l2
!     !
!     !   ! ome1=2.43d0; alf1=(2.6628045336877375,-6.06710590373065392E-002)
!     !   jloc=149
!     !
!     !   print*, lambda
!     !       print*, phi_i0(:, jloc)
!     !       print*, phi_i1(:, jloc)
!     !       print*, dphidx(:, jloc)
!     !       print*, dphidj(:, jloc)
!     !       print*, dphidjj(:, jloc)
!     !       do l1=1, 5
!     !         do l2=1, 5
!     !           print*, l1, l2, C0(l1, l2, jloc), D0(l1, l2, jloc)
!     !         enddo
!     !       enddo
!     !
!     !       ! print*, matmul(D_0(:, :, jloc), phi_i0(1:5, jloc))+ &
!     !       ! matmul(B_0(:, :, jloc), dphidj(:, jloc)*grid_profile%eta_y(jloc)) +&
!     !       ! matmul(V_0(:, :, jloc), dphidj(:, jloc)*grid_profile%eta_yy(jloc)+dphidjj(1:5, jloc)*grid_profile%eta_y(jloc)**2)
!     !
!     !       print*, matmul(D_0(:, :, jloc), phi_i0(1:5, jloc))+ &
!     !       matmul(B_0(:, :, jloc), dphidj(:, jloc)) +&
!     !       matmul(V_0(:, :, jloc), dphidjj(:, jloc))+ &
!     !       matmul(D_2(:, :, jloc), dphidi(:, jloc))
!     !
!     !       read(*, *)
!     ! endblock

!     do j=ybg+1, yeg+1
!       phi_ke_j(j)=(phi_i0(2, j)*conjg(phi_i0(2, j))+ &
!                phi_i0(3, j)*conjg(phi_i0(3, j))+ &
!                phi_i0(4, j)*conjg(phi_i0(4, j)))*hx(j)
!       dphi_ke_j(j)=(dphidx(2, j)*conjg(phi_i0(2, j))+ &
!                dphidx(3, j)*conjg(phi_i0(3, j))+ &
!                dphidx(4, j)*conjg(phi_i0(4, j)))*hx(j)
!     enddo

!     call DMDAVecGetArrayReadF90(CoordDA, Coord_local, grid, ierr)
!       phi_ke=Trap(phi_ke_j, real(grid(1, iloc-1, ybg:yeg)))
!       dphi_ke=Trap(dphi_ke_j, real(grid(1, iloc-1, ybg:yeg)))
!     call DMDAVecRestoreArrayReadF90(CoordDA, Coord_local, grid, ierr)

!     alf_old=alf

!     if(rank /= 0)then
!       call MPI_Send(phi_ke, 1, MPI_DOUBLE_COMPLEX , 0, 0, PETSC_COMM_WORLD, ierr)
!       call MPI_Send(dphi_ke, 1, MPI_DOUBLE_COMPLEX , 0, 1, PETSC_COMM_WORLD, ierr)
!     else
!       phi_ke_total=phi_ke
!       dphi_ke_total=dphi_ke
!       do j=1, size-1
!         call MPI_Recv(phi_ke, 1, MPI_DOUBLE_COMPLEX, j, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!         phi_ke_total=phi_ke_total+phi_ke
!         call MPI_Recv(dphi_ke, 1, MPI_DOUBLE_COMPLEX, j, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!         dphi_ke_total=dphi_ke_total+dphi_ke
!       enddo
!       if(screen_flg) print*, 'o:', phi_ke_total, 'd:', dphi_ke_total
!       alf=alf-conjg(ZI*1.0d0/phi_ke_total*dphi_ke_total)
!       if(screen_flg)then
!         print*, dphi_ke_total, phi_ke_total, dphi_ke_total/phi_ke_total
!         print*, 'The new alpha=', alf
!       endif
!     endif
!     call MPI_bcast(alf, 1, MPI_DOUBLE_COMPLEX, 0, PETSC_COMM_WORLD, ierr)
!     if(solve_alpha==SOLVE_ALPHA_REAL)then
!       alf=real(alf)
!     endif
!   else
!     alf=alf
!     alf_old=alf
!     if(screen_flg)then
!       print*, 'The new alpha=', alf
!     endif
!   endif

! end subroutine solve_alpha_local
