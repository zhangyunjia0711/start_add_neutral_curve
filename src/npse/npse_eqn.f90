!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  module: module_lpse_eqn
!> @breif lpse equations.
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

module mod_npse_eqn

    use petsc

    use mod_lst_baseflow_profile
    use mod_lst_grid_profile

    use iso_fortran_env, only: REAL64, REAL32

    use mod_lpse_eqn, only: type_lpse_eqn

    implicit none

#include "../include/start.f90"

    PetscErrorCode :: ierr

    private
    ! public :: solve_lpse_eqn
    ! public :: write_lpse_function
    ! public :: initialize_lpse_eqn
    ! public :: finalize_lpse_eqn
    public :: type_npse_eqn

    real(REAL64), allocatable, dimension(:, :, :) :: G0, A0, B0, C0, D0
    real(REAL64), allocatable, dimension(:, :, :) :: V0xx, V0xy, V0xz
    real(REAL64), allocatable, dimension(:, :, :) :: V0yy, V0yz, V0zz
    real(REAL64), allocatable, dimension(:, :, :) :: Aplus
    real(REAL64), allocatable, dimension(:, :) :: modify_Px_A

    ! public :: DM_LPSE

    type, extends(type_lpse_eqn) :: type_npse_eqn

        PetscScalar, allocatable :: force(:, :)
        PetscReal, allocatable :: dt(:)
        real(REAL64) :: normrhs
        Mat :: L_NPSE
        type(type_lst_bf_profile) :: baseflow_with_distortion
        complex(REAL64), allocatable, dimension(:, :, :) :: B_0_npse, D_0_npse, V_0_npse
        complex(REAL64), allocatable, dimension(:, :, :) :: D_2_npse, G_0_npse

    contains

        procedure :: initial_force
        procedure :: set_force
        procedure :: set_baseflow_npse
        procedure :: npse_setup
        procedure :: solve_baseflow_matrix_npse
        procedure :: solve_bf_matrix_j_npse
        procedure :: solve_npse_matrix
        procedure :: solve_npse_matrix_j
        procedure :: solve_dt => solve_dt_npse
        procedure :: solve_left_matrix_npse
        procedure :: solve_rhs_new_npse
        procedure :: solve_npse => npse_solve
        procedure :: solve_npse_bc
        procedure :: solve_npse_bc_j
        procedure :: solve_npse_phi

        procedure :: finalize_force

        procedure :: solve_rhs_new

        !   procedure :: get=>get_array
        !   procedure :: set=>set_array

    end type type_npse_eqn

contains

    subroutine initial_force(this)

        implicit none
        class(type_npse_eqn), intent(inout) :: this

        ! jn=jn_
        allocate (this%force(5, this%jn)); this%force = 0.0d0
        allocate (this%dt(this%jn)); this%dt = MAXNUM
        call DMCreateMatrix(this%dm, this%L_NPSE, ierr)
        allocate (G0(5, 5, this%jn)); G0 = 0.0d0
        allocate (A0(5, 5, this%jn)); A0 = 0.0d0
        allocate (B0(5, 5, this%jn)); B0 = 0.0d0
        allocate (C0(5, 5, this%jn)); C0 = 0.0d0
        allocate (D0(5, 5, this%jn)); D0 = 0.0d0
        allocate (V0xx(5, 5, this%jn)); V0xx = 0.0d0
        allocate (V0xy(5, 5, this%jn)); V0xy = 0.0d0
        allocate (V0xz(5, 5, this%jn)); V0xz = 0.0d0
        allocate (V0yy(5, 5, this%jn)); V0yy = 0.0d0
        allocate (V0yz(5, 5, this%jn)); V0yz = 0.0d0
        allocate (V0zz(5, 5, this%jn)); V0zz = 0.0d0
        allocate (Aplus(5, 5, this%jn)); Aplus = 0.0d0
        allocate (modify_Px_A(2, this%jn)); modify_Px_A = 0.0d0
        allocate (this%B_0_npse(5, 5, this%jn)); this%B_0_npse = 0.0d0
        allocate (this%D_0_npse(5, 5, this%jn)); this%D_0_npse = 0.0d0
        allocate (this%V_0_npse(5, 5, this%jn)); this%V_0_npse = 0.0d0
        allocate (this%D_2_npse(5, 5, this%jn)); this%D_2_npse = 0.0d0
        allocate (this%G_0_npse(5, 5, this%jn)); this%G_0_npse = 0.0d0

    end subroutine initial_force

    subroutine set_force(this, force)

        implicit none
        class(type_npse_eqn), intent(inout) :: this
        complex(REAL64), intent(in) :: force(:, :)

        this%force = force

    end subroutine set_force

    subroutine set_baseflow_npse(this, baseflow)
        implicit none
        class(type_npse_eqn), intent(inout) :: this
        type(type_lst_bf_profile), intent(in) :: baseflow

        this%baseflow_with_distortion = baseflow

    end subroutine set_baseflow_npse

    subroutine npse_setup(this)

        implicit none
        class(type_npse_eqn), intent(inout) :: this

        call this%solve_baseflow_matrix_npse()

    end subroutine npse_setup

    subroutine solve_baseflow_matrix_npse(this)

        implicit none
        class(type_npse_eqn) :: this
        integer :: j, j_b, j_e

        !! internal point
        call this%solve_bf_matrix_j_npse()

    end subroutine solve_baseflow_matrix_npse

    subroutine solve_bf_matrix_j_npse(this)

        implicit none
        class(type_npse_eqn) :: this
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

        associate (grid_profile => this%grid, bf_profile => this%baseflow_with_distortion, Re => this%Re, Ma => this%Ma, &
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

    end subroutine solve_bf_matrix_j_npse

    subroutine solve_npse_matrix(this, alf)
        implicit none
        class(type_npse_eqn) :: this
        complex(REAL64), intent(in) :: alf

        integer :: j

        do j = this%jb, this%je
            call this%solve_npse_matrix_j(j, alf)
        end do

    end subroutine solve_npse_matrix

    subroutine solve_npse_matrix_j(this, j, alf)

        implicit none
        class(type_npse_eqn) :: this
        integer, intent(in) :: j
        complex(REAL64), intent(in) :: alf
        complex(REAL64) :: dadx

        complex(REAL64), dimension(5, 5) :: A_, B_, V_

        !!d\alpha/dx
!        dadx=(alf-wave_param(2)%alf)*grid_profile%xi_x(jb)
        associate (grid_profile => this%grid, bf_profile => this%baseflow_with_distortion, ome => this%ome, bat => this%bat)
            associate (A_0 => this%D_2_npse, B_0 => this%B_0_npse, D_0 => this%D_0_npse, V_0 => this%V_0_npse, &
                       G_0 => this%G_0_npse, &
                       xi_x => grid_profile%xi_x(j), &
                       xi_y => grid_profile%xi_y(j), &
                       eta_x => grid_profile%eta_x(j), &
                       eta_y => grid_profile%eta_y(j), &
                       xi_yy => grid_profile%xi_yy(j), &
                       eta_yy => grid_profile%eta_yy(j))

                D_0(:, :, j) = D0(:, :, j) - ZI*ome*G0(:, :, j) + ZI*alf*A0(:, :, j) &
                               + ZI*bat*C0(:, :, j) + alf**2*V0xx(:, :, j) &
                               + alf*bat*V0xz(:, :, j) + bat**2*V0zz(:, :, j)! &
                !-ZI*dadx*V0xx(:, :, j)
                G_0(:, :, j) = G0(:, :, j)
!          D_0(2, 1, j)=D_0(2, 1, j)+modify_Px_D(1, j)
!          D_0(2, 5, j)=D_0(2, 5, j)+modify_Px_D(2, j)
                A_0(:, :, j) = A0(:, :, j) - ZI*bat*V0xz(:, :, j) - 2.0d0*ZI*alf*V0xx(:, :, j)
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

    end subroutine solve_npse_matrix_j

    subroutine npse_solve(this, alf, phi0, phi1, phi2, lambda)

        use stdlib_optval, only: optval
        implicit none
        class(type_npse_eqn), intent(inout) :: this
        complex(REAL64), intent(inout) :: alf
        Vec, intent(inout) ::  phi0
        Vec, intent(in), optional :: phi1, phi2
        real(REAL64), optional :: lambda
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

        this%lambda = optval(lambda, this%lambda)

        alf_old = alf
        delta = 1.0d10
        ! read(*, *)
        this%iterator_times = -1

        ! call this%npse_setup()
        do while (abs(delta) >= (epsilon(1.0e0)*100.0d0))

            call this%solve_lpse_matrix(alf)
            call this%solve_npse_matrix(alf)

            call this%solve_left_matrix_npse()

            call this%solve_left_matrix()

            call this%solve_lpse_matrix_bc(alf)

            call this%solve_rhs_new_npse()

            call this%solve_lpse_bc()
            call this%solve_npse_bc()

            call this%solve_npse_phi()
            call this%solve_alpha_local_new(alf)
            ! print*, alf, alf_old, (alf-alf_old)
            delta = abs(alf - alf_old)
            alf_old = alf
            ! call VecWAXPY(dphi, (-1.0d0, 0.0d0), phi_old, phi_new, ierr)
            ! call VecNorm(dphi, NORM_2, delta, ierr)
            ! call VecEqual(phi_new, phi_old, flg, ierr)
            this%iterator_times = this%iterator_times + 1
        end do

        if (this%solve_alpha /= SOLVE_ALPHA_NONE) print *, "the iteration times is", this%iterator_times

        call VecCopy(this%phi, phi0, ierr)
        call DMDAVecGetArrayReadF90(this%dm, this%phi, tmp, ierr)
        this%phi0_array = tmp
        call DMDAVecRestoreArrayReadF90(this%dm, this%phi, tmp, ierr)

    end subroutine npse_solve

    subroutine solve_rhs_new(this)

        use functional
        implicit none
        class(type_npse_eqn) :: this

        integer :: bs
        PetscScalar :: rhs_j(5, this%jb:this%je), box(5*(this%je - this%jb + 1))
        integer :: j

      !! RHS=lambda*xi_x*A(phi(i-2)-phi(i-1))+(1+lambda)*xi_x*A*phi(i-1)

        associate (lambda => this%lambda, RHS => this%RHS, A_0 => this%D_2, jb => this%jb, je => this%je)
            do j = jb, je
                ! associate (A_0 => D_2)
                rhs_j(:, j) = lambda*matmul(A_0(:, :, j), (this%phi1_array(:, j) - this%phi2_array(:, j))) &
                              + (1.0d0 + lambda)*matmul(A_0(:, :, j), this%phi1_array(:, j)) + this%force(:, j)
                ! end associate
            end do

            this%normrhs = sqrt(sum(rhs_j(:, jb + 1:je - 1)*conjg(rhs_j(:, jb + 1:je - 1))))

            box = reshape(rhs_j, [5*(je - jb + 1)])
            call VecSetValuesBlocked(RHS, je - jb + 1, arange(jb - 1, je - 1), box, INSERT_VALUES, ierr)
            call VecAssemblyBegin(RHS, ierr)
            call VecAssemblyEnd(RHS, ierr)
            call VecGetBlockSize(RHS, bs, ierr)
            ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)
            !   read(*, *)
        end associate

    end subroutine solve_rhs_new

    subroutine solve_left_matrix_npse(this)

        implicit none
        class(type_npse_eqn), intent(in) :: this
        real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

        PetscScalar :: box1(5, -2:2, 5)
        complex(REAL64) :: tmp_mat_j0(5, 5, -2:2)
        ! real(REAL64) :: eta_y, eta_yy

        integer :: i, j, lj
        integer :: ljb, lje, c_idx
        MatStencil :: idxm(4, 1), idxn(4, -2:2)
        PetscViewer :: Viewer

        integer :: jb_local, je_local

        associate (A_0 => this%D_2_npse, B_0 => this%B_0_npse, D_0 => this%D_0_npse, V_0 => this%V_0_npse, jn => this%jn, &
                   G_0 => this%G_0_npse, C1 => FDM_C1_4ORD_CENTER, D1 => FDM_D1_4ORD_CENTER)
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
                    tmp_mat_j0(:, :, lj) = delta(lj)*((1.0d0 + this%lambda)*A_0(:, :, j) + D_0(:, :, j) &
                                                    &  + G_0(:, :, j)/this%dt(j)) + &
                                           C1(lj, c_idx)*(B_0(:, :, j)) + &
                                           D1(lj, c_idx)*(V_0(:, :, j))
                    box1(:, lj, :) = transpose(tmp_mat_j0(:, :, lj))
                end do

                call MatSetValuesBlockedStencil(this%L_NPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)
            end do
        end associate

        call MatAssemblyBegin(this%L_NPSE, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(this%L_NPSE, MAT_FINAL_ASSEMBLY, ierr)

    end subroutine solve_left_matrix_npse

    subroutine solve_npse_bc(this)

        implicit none
        class(type_npse_eqn), intent(in) :: this

        integer :: j

        do j = this%jb, this%je
            if (j == 1 .or. j == this%jn) then
                call this%solve_npse_bc_j(j)
            end if
        end do

    end subroutine solve_npse_bc

    subroutine solve_npse_bc_j(this, j)

        use functional
        implicit none
        class(type_npse_eqn), intent(in) :: this
        integer, intent(in) :: j
        real(REAL64), parameter :: delta(-2:2) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]

        complex(REAL64) :: box1(5, -2:2, 5)
        complex(REAL64) :: tmp_mat_j0(5, 5, -2:2)
        complex(REAL64) :: tmpj0(4), tmpjn(5)
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

            call MatSetValuesBlockedStencil(this%L_NPSE, 1, idxm, 5, idxn, box1, INSERT_VALUES, ierr)

        end associate

        call MatAssemblyBegin(this%L_NPSE, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(this%L_NPSE, MAT_FINAL_ASSEMBLY, ierr)

        ! if(screen_flg)then
        !   call PetscPrintf(PETSC_COMM_WORLD, "Set spatial LPSE matrix and RHS BC done!\n", ierr)
        ! endif
        ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)

    end subroutine solve_npse_bc_j

    subroutine solve_rhs_new_npse(this)

        use functional
        implicit none
        class(type_npse_eqn) :: this

        ! integer :: bs
        PetscScalar :: rhs_j(5, this%jb:this%je), box(5*(this%je - this%jb + 1))
        PetscScalar :: rhs1_j(5, this%jb:this%je)
        integer :: j

        Vec :: tmp

        call VecDuplicate(this%phi, tmp, ierr)
      !! RHS=lambda*xi_x*A(phi(i-2)-phi(i-1))+(1+lambda)*xi_x*A*phi(i-1)

        associate (lambda => this%lambda, RHS => this%RHS, A_0 => this%D_2, jb => this%jb, je => this%je)
            do j = jb, je
                ! associate (A_0 => D_2)
                rhs1_j(:, j) = lambda*matmul(A_0(:, :, j), (this%phi1_array(:, j) - this%phi2_array(:, j))) &
                            & + (1.0d0 + lambda)*matmul(A_0(:, :, j), this%phi1_array(:, j)) + this%force(:, j)
                rhs_j(:, j) = rhs1_j(:, j) !+  matmul(this%G_0(:, :, j)/this%dt(j), this%phi0_array(:, j))
                ! end associate
            end do

            this%normrhs = sqrt(sum(rhs1_j(:, jb + 1:je - 1)*conjg(rhs1_j(:, jb + 1:je - 1))))

            box = reshape(rhs_j, [5*(je - jb + 1)])
            call VecSetValuesBlocked(RHS, je - jb + 1, arange(jb - 1, je - 1), box, INSERT_VALUES, ierr)
            call VecAssemblyBegin(RHS, ierr)
            call VecAssemblyEnd(RHS, ierr)
            call MatMult(this%L_LPSE, this%phi, tmp, ierr)
            call VecAXPY(RHS, (-1.0d0, 0.0d0), tmp, ierr)
            ! call VecGetBlockSize(RHS, bs, ierr)
            ! call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, ierr)
            !   read(*, *)
        end associate

        call VecDestroy(tmp, ierr)

    end subroutine solve_rhs_new_npse

    subroutine solve_npse_phi(this)

        implicit none
        class(type_npse_eqn), intent(inout) :: this
        integer :: its
        character(20) :: ksptype
        Vec :: delta_phi

        call VecDuplicate(this%phi, delta_phi, ierr)
        call KSPSetOperators(this%Solver, this%L_NPSE, this%L_NPSE, ierr)
        call KSPSetFromOptions(this%solver, ierr)
        call KSPSetup(this%Solver, ierr)
        call KSPSolve(this%Solver, this%RHS, delta_phi, ierr)
        call VecAXPY(this%phi, (1.0d0, 0.0d0), delta_phi, ierr)
        call VecDestroy(delta_phi, ierr)

    end subroutine solve_npse_phi

    subroutine solve_dt_npse(this, cfl_number)
        implicit none
        class(type_npse_eqn), intent(inout) :: this
        real(REAL64), intent(in) :: cfl_number
        integer :: j
        real(REAL64) :: acoustic, t1, t2, t3
        real(REAL64) :: cm, miu, mu
        real(REAL64) :: dxi, deta
        real(REAL64) :: ubar, vbar, wbar

        do j = 2, this%jn - 1
            associate (                          &
                &  rho => this%baseflow%q0(1, j), &
                &  u => this%baseflow%q0(2, j), &
                &  v => this%baseflow%q0(3, j), &
                &  w => this%baseflow%q0(4, j), &
                &  T => this%baseflow%q0(5, j), &
                & xi_x => this%grid%xi_x(j), &
                & xi_y => this%grid%xi_y(j), &
                & eta_x => this%grid%eta_x(j), &
                & eta_y => this%grid%eta_y(j) &
                )

                cm = this%C1/this%Te
                miu = T*sqrt(T)*(1.0d0 + cm)/(T + cm)
                mu = miu*max(4.0d0/3.0d0, this%GAMMA/this%Pr)/this%Re/rho

                acoustic = this%Ma/sqrt(t)

                dxi = sqrt(xi_x**2 + xi_y**2)
                deta = sqrt(eta_x**2 + eta_y**2)

                ubar = (xi_x*u + xi_y*v)/dxi
                vbar = (eta_x*u + eta_y*v)/deta
                wbar = w

                t1 = abs(ubar) + acoustic + 2.0d0*dxi*mu
                t2 = abs(vbar) + acoustic + 2.0d0*deta*mu
                t3 = abs(wbar) + acoustic + 2.0d0*1.0d0*mu

                this%dt(j) = cfl_number/(dxi*t1 + deta*t2 + t3)
            end associate
        end do

    end subroutine solve_dt_npse

    subroutine finalize_force(this)

        implicit none
        class(type_npse_eqn), intent(inout) :: this

        ! jn=jn_
        deallocate (this%force)
        deallocate (this%dt)

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
        deallocate (modify_Px_A)

        deallocate (this%B_0_npse)
        deallocate (this%D_0_npse)
        deallocate (this%V_0_npse)
        deallocate (this%D_2_npse)
        deallocate (this%G_0_npse)

        call MatDestroy(this%L_NPSE, ierr)

    end subroutine finalize_force

    ! subroutine set_cfl(this, cfl)
    !     implicit none
    !     class(type_npse_eqn), intent(inout) :: this
    !     real(REAL64), intent(in) :: cfl

    !     this%cfl_factor=cfl
    ! end subroutine set_cfl

end module mod_npse_eqn
