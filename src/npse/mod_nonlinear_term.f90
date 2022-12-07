
!#include <mkl_dfti.f90>

#include <slepc/finclude/slepc.h>
#include <petsc/finclude/petsc.h>

! module MKL_fft

!     use MKL_DFTI

!     implicit none

!     integer :: MNSpectral(2)
!     integer :: status

! contains

!     subroutine SetSpectralScheme(handle)

!         implicit none
!         type(DFTI_DESCRIPTOR), pointer :: handle
!         real*8  :: ScaleFactor

!         ScaleFactor = 1.0d0/dble(MNSpectral(1)*MNSpectral(2))
!         status = DftiCreateDescriptor(handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, MNSpectral(1:2))

!         ! Status = DftiSetValue (handle, DFTI_PLACEMENT, DFTI_INPLACE)
!         Status = DftiSetValue(handle, DFTI_FORWARD_SCALE, ScaleFactor)
!         status = DftiSetValue(handle, DFTI_NUMBER_OF_USER_THREADS, 4)
!         Status = DftiCommitDescriptor(handle)

!     end subroutine SetSpectralScheme

!     !---inv=0,for forward translation,1 for backward translation
!     subroutine InitialSpectralScheme(handle, Minput, Ninput)

!         implicit none
!         type(DFTI_DESCRIPTOR), pointer :: handle
!         integer :: Ninput, Minput
!         integer ::inv
!         real*8::ScaleFactor

!         MNSpectral(1) = Minput
!         MNSpectral(2) = Ninput
!         call SetSpectralScheme(handle)

!     end subroutine InitialSpectralScheme

!     subroutine SpectralForward(handle, X)
!         implicit none
!         type(DFTI_DESCRIPTOR), pointer :: handle
!         complex*16::X(MNSpectral(1), MNSpectral(2))
!         !complex*16,allocatable::XX(:) !(MNSpectral(1)*MNSpectral(2))
!         !Equivalence(X,XX)
!         status = DftiComputeForward(handle, X(1:MNSpectral(1), 1))
!     end subroutine

!     subroutine SpectralBackward(handle, X)
!         implicit none
!         type(DFTI_DESCRIPTOR), pointer :: handle
!         complex*16::X(MNSpectral(1), MNSpectral(2))
!         !complex*16,allocatable::XX(:) !(MNSpectral(1)*MNSpectral(2))
!         !Equivalence(X,XX)
!         status = DftiComputeBackward(handle, X(1:MNSpectral(1), 1))

!     end subroutine

!     subroutine FinalizeSpectralScheme(handle)
!         type(DFTI_DESCRIPTOR), pointer :: handle

!         Status = DftiFreeDescriptor(handle)

!     end subroutine FinalizeSpectralScheme

! end module MKL_fft

module mod_nonlinear_term

    use iso_fortran_env, only: REAL64, REAL32
    use stdlib_error, only: error_stop

    use mod_lst_grid_profile
    use mod_lst_baseflow_profile
    use mod_solver_abstract

    use petsc
    !use MKL_fft
! #ifdef MKL
    use mod_fft
! #else
    use mod_petsc_fft
! #endif

    implicit none

    private
    public :: type_nonlinear_term

#include "../include/start.f90"

    type, extends(type_solver_abstract) :: type_nonlinear_term

        integer :: m
        integer :: n
        integer :: msp, nsp
        integer :: jn

        DM :: dm

        real(REAL64) :: ome1, ome2
        real(REAL64) :: bat1, bat2
        complex(REAL64), allocatable :: alf(:, :)
        complex(REAL64), allocatable :: exp_integral_alf(:, :)

        complex(REAL64), allocatable :: f_spectrum(:, :, :, :)
        real(REAL64), allocatable :: f_physics(:, :, :, :)

        type(type_point_spectrum), allocatable :: dis_profile_spectrum(:)

        PetscScalar, allocatable :: phi0_array(:, :, :, :), phi1_array(:, :, :, :), phi2_array(:, :, :, :)

        real(REAL64) :: lambda = 0.0d0
        ! type(DFTI_DESCRIPTOR), pointer :: handle
        type(fft_type) :: handle_new

    contains

        procedure :: initial => initial_nonliner_term
        procedure :: set_dis => set_dis_nonliner_term
        procedure :: setup => setup_nonliner_term
        procedure :: get_force => get_force_nonlinear_term
        procedure :: finalize => finalize_nonliner_term

        procedure, private :: compute_dis_profile
        procedure, private :: compute_nonlinear_term
        procedure, private :: compute_force
        procedure, private :: ifft_spt_phy
        procedure, private :: ifft_var
        procedure, private :: fft_phy_spt
        procedure, private :: fft_var

        procedure :: set_parameters_real4 => nonliner_term_set_parameters_real4
        procedure :: set_parameters_real8 => nonliner_term_set_parameters_real8
        procedure :: set_parameters_int => nonliner_term_set_parameters_int

    end type type_nonlinear_term

    type :: type_point_spectrum

        integer :: m
        integer :: n
        complex(REAL64), dimension(:, :), allocatable :: rho_a0, u_a0, v_a0, w_a0, T_a0
        complex(REAL64), dimension(:, :), allocatable :: rho_t, u_t, v_t, w_t, T_t
        complex(REAL64), dimension(:, :), allocatable :: rho_x, u_x, v_x, w_x, T_x
        complex(REAL64), dimension(:, :), allocatable :: rho_y, u_y, v_y, w_y, T_y
        complex(REAL64), dimension(:, :), allocatable :: rho_z, u_z, v_z, w_z, T_z
        complex(REAL64), dimension(:, :), allocatable :: u_xx, v_xx, w_xx, T_xx
        complex(REAL64), dimension(:, :), allocatable :: u_yy, v_yy, w_yy, T_yy
        complex(REAL64), dimension(:, :), allocatable :: u_zz, v_zz, w_zz, T_zz

        complex(REAL64), dimension(:, :), allocatable :: u_xy, v_xy
        complex(REAL64), dimension(:, :), allocatable :: u_xz, w_xz
        complex(REAL64), dimension(:, :), allocatable :: v_yz, w_yz

    end type type_point_spectrum

    type :: type_point_physics

        real(REAL64), dimension(:, :), allocatable :: rho_a0, u_a0, v_a0, w_a0, T_a0
        real(REAL64), dimension(:, :), allocatable :: rho_t, u_t, v_t, w_t, T_t
        real(REAL64), dimension(:, :), allocatable :: rho_x, u_x, v_x, w_x, T_x
        real(REAL64), dimension(:, :), allocatable :: rho_y, u_y, v_y, w_y, T_y
        real(REAL64), dimension(:, :), allocatable :: rho_z, u_z, v_z, w_z, T_z
        real(REAL64), dimension(:, :), allocatable :: u_xx, v_xx, w_xx, T_xx
        real(REAL64), dimension(:, :), allocatable :: u_yy, v_yy, w_yy, T_yy
        real(REAL64), dimension(:, :), allocatable :: u_zz, v_zz, w_zz, T_zz

        real(REAL64), dimension(:, :), allocatable :: u_xy, v_xy
        real(REAL64), dimension(:, :), allocatable :: u_xz, w_xz
        real(REAL64), dimension(:, :), allocatable :: v_yz, w_yz

    end type type_point_physics

    type(type_point_physics) :: point_physics

    complex(REAL64), allocatable :: dphidi (:, :, :, :)
    complex(REAL64), allocatable :: dphidj (:, :, :, :)
    complex(REAL64), allocatable :: dphidij(:, :, :, :)
    complex(REAL64), allocatable :: dphidjj(:, :, :, :)
    complex(REAL64), allocatable :: dphidx (:, :, :, :)
    complex(REAL64), allocatable :: dphidy (:, :, :, :)
    complex(REAL64), allocatable :: dphidxx(:, :, :, :)
    complex(REAL64), allocatable :: dphidxy(:, :, :, :)
    complex(REAL64), allocatable :: dphidyy(:, :, :, :)


contains

    subroutine initial_nonliner_term(this, jn_, Re, Ma, Te, yb_, ye_)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        integer, intent(in) :: jn_
        real(real64), intent(in) ::  Re, Ma, Te
        integer, intent(in), optional :: yb_, ye_
        integer :: j
        PetscErrorCode :: ierr

        this%jn = jn_
        this%Re = Re
        this%Ma = Ma
        this%Te = Te

        call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn_, 5, 2, PETSC_NULL_INTEGER, this%dm, ierr)
        call DMSetUp(this%dm, ierr)

        associate (m => this%m, n => this%n, msp => this%msp, nsp => this%nsp)
            allocate (this%alf(0:m, -n:n))
            allocate (this%exp_integral_alf(0:m, -n:n))
            allocate (this%f_spectrum(5, this%jn, 0:m, -n:n))
            allocate (this%f_physics(5, this%jn, 0:this%msp - 1, 0:this%nsp - 1))
            allocate (this%dis_profile_spectrum(jn_))
            allocate (this%phi0_array(5, this%jn, 0:m, -n:n)); this%phi0_array = 0.0d0
            allocate (this%phi1_array(5, this%jn, 0:m, -n:n)); this%phi1_array = 0.0d0
            allocate (this%phi2_array(5, this%jn, 0:m, -n:n)); this%phi2_array = 0.0d0
            do j = 1, jn_
                this%dis_profile_spectrum(j)%m = m
                this%dis_profile_spectrum(j)%n = n
                allocate (this%dis_profile_spectrum(j)%rho_a0(0:m, -n:n)); this%dis_profile_spectrum(j)%rho_a0 = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_a0(0:m, -n:n)); this%dis_profile_spectrum(j)%u_a0 = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_a0(0:m, -n:n)); this%dis_profile_spectrum(j)%v_a0 = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_a0(0:m, -n:n)); this%dis_profile_spectrum(j)%w_a0 = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_a0(0:m, -n:n)); this%dis_profile_spectrum(j)%T_a0 = 0.0d0
                allocate (this%dis_profile_spectrum(j)%rho_t(0:m, -n:n)); this%dis_profile_spectrum(j)%rho_t = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_t(0:m, -n:n)); this%dis_profile_spectrum(j)%u_t = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_t(0:m, -n:n)); this%dis_profile_spectrum(j)%v_t = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_t(0:m, -n:n)); this%dis_profile_spectrum(j)%w_t = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_t(0:m, -n:n)); this%dis_profile_spectrum(j)%T_t = 0.0d0
                allocate (this%dis_profile_spectrum(j)%rho_x(0:m, -n:n)); this%dis_profile_spectrum(j)%rho_x = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_x(0:m, -n:n)); this%dis_profile_spectrum(j)%u_x = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_x(0:m, -n:n)); this%dis_profile_spectrum(j)%v_x = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_x(0:m, -n:n)); this%dis_profile_spectrum(j)%w_x = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_x(0:m, -n:n)); this%dis_profile_spectrum(j)%T_x = 0.0d0
                allocate (this%dis_profile_spectrum(j)%rho_y(0:m, -n:n)); this%dis_profile_spectrum(j)%rho_y = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_y(0:m, -n:n)); this%dis_profile_spectrum(j)%u_y = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_y(0:m, -n:n)); this%dis_profile_spectrum(j)%v_y = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_y(0:m, -n:n)); this%dis_profile_spectrum(j)%w_y = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_y(0:m, -n:n)); this%dis_profile_spectrum(j)%T_y = 0.0d0
                allocate (this%dis_profile_spectrum(j)%rho_z(0:m, -n:n)); this%dis_profile_spectrum(j)%rho_z = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_z(0:m, -n:n)); this%dis_profile_spectrum(j)%u_z = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_z(0:m, -n:n)); this%dis_profile_spectrum(j)%v_z = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_z(0:m, -n:n)); this%dis_profile_spectrum(j)%w_z = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_z(0:m, -n:n)); this%dis_profile_spectrum(j)%T_z = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_xx(0:m, -n:n)); this%dis_profile_spectrum(j)%u_xx = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_xx(0:m, -n:n)); this%dis_profile_spectrum(j)%v_xx = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_xx(0:m, -n:n)); this%dis_profile_spectrum(j)%w_xx = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_xx(0:m, -n:n)); this%dis_profile_spectrum(j)%T_xx = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_yy(0:m, -n:n)); this%dis_profile_spectrum(j)%u_yy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_yy(0:m, -n:n)); this%dis_profile_spectrum(j)%v_yy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_yy(0:m, -n:n)); this%dis_profile_spectrum(j)%w_yy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_yy(0:m, -n:n)); this%dis_profile_spectrum(j)%T_yy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_zz(0:m, -n:n)); this%dis_profile_spectrum(j)%u_zz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_zz(0:m, -n:n)); this%dis_profile_spectrum(j)%v_zz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_zz(0:m, -n:n)); this%dis_profile_spectrum(j)%w_zz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%T_zz(0:m, -n:n)); this%dis_profile_spectrum(j)%T_zz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_xy(0:m, -n:n)); this%dis_profile_spectrum(j)%u_xy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_xy(0:m, -n:n)); this%dis_profile_spectrum(j)%v_xy = 0.0d0
                allocate (this%dis_profile_spectrum(j)%u_xz(0:m, -n:n)); this%dis_profile_spectrum(j)%u_xz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_xz(0:m, -n:n)); this%dis_profile_spectrum(j)%w_xz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%v_yz(0:m, -n:n)); this%dis_profile_spectrum(j)%v_yz = 0.0d0
                allocate (this%dis_profile_spectrum(j)%w_yz(0:m, -n:n)); this%dis_profile_spectrum(j)%w_yz = 0.0d0
            end do

            allocate (point_physics%rho_a0(0:msp - 1, 0:nsp - 1)); point_physics%rho_a0 = 0.0d0
            allocate (point_physics%u_a0(0:msp - 1, 0:nsp - 1)); point_physics%u_a0 = 0.0d0
            allocate (point_physics%v_a0(0:msp - 1, 0:nsp - 1)); point_physics%v_a0 = 0.0d0
            allocate (point_physics%w_a0(0:msp - 1, 0:nsp - 1)); point_physics%w_a0 = 0.0d0
            allocate (point_physics%T_a0(0:msp - 1, 0:nsp - 1)); point_physics%T_a0 = 0.0d0
            allocate (point_physics%rho_t(0:msp - 1, 0:nsp - 1)); point_physics%rho_t = 0.0d0
            allocate (point_physics%u_t(0:msp - 1, 0:nsp - 1)); point_physics%u_t = 0.0d0
            allocate (point_physics%v_t(0:msp - 1, 0:nsp - 1)); point_physics%v_t = 0.0d0
            allocate (point_physics%w_t(0:msp - 1, 0:nsp - 1)); point_physics%w_t = 0.0d0
            allocate (point_physics%T_t(0:msp - 1, 0:nsp - 1)); point_physics%T_t = 0.0d0
            allocate (point_physics%rho_x(0:msp - 1, 0:nsp - 1)); point_physics%rho_x = 0.0d0
            allocate (point_physics%u_x(0:msp - 1, 0:nsp - 1)); point_physics%u_x = 0.0d0
            allocate (point_physics%v_x(0:msp - 1, 0:nsp - 1)); point_physics%v_x = 0.0d0
            allocate (point_physics%w_x(0:msp - 1, 0:nsp - 1)); point_physics%w_x = 0.0d0
            allocate (point_physics%T_x(0:msp - 1, 0:nsp - 1)); point_physics%T_x = 0.0d0
            allocate (point_physics%rho_y(0:msp - 1, 0:nsp - 1)); point_physics%rho_y = 0.0d0
            allocate (point_physics%u_y(0:msp - 1, 0:nsp - 1)); point_physics%u_y = 0.0d0
            allocate (point_physics%v_y(0:msp - 1, 0:nsp - 1)); point_physics%v_y = 0.0d0
            allocate (point_physics%w_y(0:msp - 1, 0:nsp - 1)); point_physics%w_y = 0.0d0
            allocate (point_physics%T_y(0:msp - 1, 0:nsp - 1)); point_physics%T_y = 0.0d0
            allocate (point_physics%rho_z(0:msp - 1, 0:nsp - 1)); point_physics%rho_z = 0.0d0
            allocate (point_physics%u_z(0:msp - 1, 0:nsp - 1)); point_physics%u_z = 0.0d0
            allocate (point_physics%v_z(0:msp - 1, 0:nsp - 1)); point_physics%v_z = 0.0d0
            allocate (point_physics%w_z(0:msp - 1, 0:nsp - 1)); point_physics%w_z = 0.0d0
            allocate (point_physics%T_z(0:msp - 1, 0:nsp - 1)); point_physics%T_z = 0.0d0
            allocate (point_physics%u_xx(0:msp - 1, 0:nsp - 1)); point_physics%u_xx = 0.0d0
            allocate (point_physics%v_xx(0:msp - 1, 0:nsp - 1)); point_physics%v_xx = 0.0d0
            allocate (point_physics%w_xx(0:msp - 1, 0:nsp - 1)); point_physics%w_xx = 0.0d0
            allocate (point_physics%T_xx(0:msp - 1, 0:nsp - 1)); point_physics%T_xx = 0.0d0
            allocate (point_physics%u_yy(0:msp - 1, 0:nsp - 1)); point_physics%u_yy = 0.0d0
            allocate (point_physics%v_yy(0:msp - 1, 0:nsp - 1)); point_physics%v_yy = 0.0d0
            allocate (point_physics%w_yy(0:msp - 1, 0:nsp - 1)); point_physics%w_yy = 0.0d0
            allocate (point_physics%T_yy(0:msp - 1, 0:nsp - 1)); point_physics%T_yy = 0.0d0
            allocate (point_physics%u_zz(0:msp - 1, 0:nsp - 1)); point_physics%u_zz = 0.0d0
            allocate (point_physics%v_zz(0:msp - 1, 0:nsp - 1)); point_physics%v_zz = 0.0d0
            allocate (point_physics%w_zz(0:msp - 1, 0:nsp - 1)); point_physics%w_zz = 0.0d0
            allocate (point_physics%T_zz(0:msp - 1, 0:nsp - 1)); point_physics%T_zz = 0.0d0
            allocate (point_physics%u_xy(0:msp - 1, 0:nsp - 1)); point_physics%u_xy = 0.0d0
            allocate (point_physics%v_xy(0:msp - 1, 0:nsp - 1)); point_physics%v_xy = 0.0d0
            allocate (point_physics%u_xz(0:msp - 1, 0:nsp - 1)); point_physics%u_xz = 0.0d0
            allocate (point_physics%w_xz(0:msp - 1, 0:nsp - 1)); point_physics%w_xz = 0.0d0
            allocate (point_physics%v_yz(0:msp - 1, 0:nsp - 1)); point_physics%v_yz = 0.0d0
            allocate (point_physics%w_yz(0:msp - 1, 0:nsp - 1)); point_physics%w_yz = 0.0d0

        end associate

        ! call InitialSpectralScheme(this%handle, 4*this%m, 4*this%n)
        call this%handle_new%initialize(this%msp, this%nsp)

        allocate(dphidi(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidj(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidij(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidjj(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidx(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidy(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidxx(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidxy(5, this%jn, 0:this%m, -this%n:this%n))
        allocate(dphidyy(5, this%jn, 0:this%m, -this%n:this%n))

    end subroutine initial_nonliner_term

    subroutine set_dis_nonliner_term(this, alf, lnA0_, phi0, phi1, phi2, lambda)

        use stdlib_optval, only: optval
        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        complex(REAL64), intent(in) :: alf(0:this%m, -this%n:this%n)
        Complex(REAL64), intent(in) :: lnA0_(0:this%m, -this%n:this%n)
        Vec, intent(inout) :: phi0(0:this%m, -this%n:this%n)
        Vec, intent(inout) :: phi1(0:this%m, -this%n:this%n)
        Vec, intent(inout) :: phi2(0:this%m, -this%n:this%n)
        real(REAL64), optional :: lambda

        integer :: m, n

        do n = -this%n, this%n
            do m = 0, this%m
                this%phi0_array(:, :, m, n) = get_array(phi0(m, n), this%dm)
                this%phi1_array(:, :, m, n) = get_array(phi1(m, n), this%dm)
                this%phi2_array(:, :, m, n) = get_array(phi2(m, n), this%dm)
            end do
        end do

        this%exp_integral_alf = exp(ZI*lnA0_)
        this%alf = alf

        if (sum(abs(this%phi2_array)) <= epsilon(1.0d0)) then  !phi2 is zero means 1 order is used
            this%lambda = 0.0d0
        else ! else 2 order is used
            this%lambda = 0.5d0
        end if

        this%lambda = optval(lambda, this%lambda)
        
    end subroutine set_dis_nonliner_term

    subroutine setup_nonliner_term(this)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this

        call this%compute_dis_profile()
        call this%compute_nonlinear_term()

    end subroutine setup_nonliner_term

    subroutine get_force_nonlinear_term(this, force)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        complex(real64), intent(inout) :: force(5, this%jn, 0:this%m, -this%n:this%n)

        force = this%f_spectrum

    end subroutine get_force_nonlinear_term

    subroutine compute_nonlinear_term(this)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this

        integer :: j, l

!!        ! $omp parallel
!!        ! $omp do
        do j = 1, this%jn
            ! print *, 'jloc=', j
            ! print *, this%dis_profile_spectrum(j)%v_a0(1, 0)
            call this%ifft_spt_phy(this%dis_profile_spectrum(j), point_physics)
            ! print *, point_physics%v_a0(0, 0)
            call this%compute_force(j, point_physics, this%f_physics(:, j, :, :))
            call this%fft_phy_spt(this%f_physics(:, j, :, :), this%f_spectrum(:, j, :, :))
            do l = 1, 5
                this%f_spectrum(l, j, :, :) = this%f_spectrum(l, j, :, :)/this%exp_integral_alf(:, :)
            end do
        end do
!!        !$omp enddo
!!        !$omp end parallel

    end subroutine compute_nonlinear_term

    subroutine fft_phy_spt(this, f_phy, f_spt)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        real(REAL64), intent(inout) :: f_phy(5, 0:this%msp - 1, 0:this%nsp - 1)
        complex(REAL64), intent(inout) :: f_spt(5, 0:this%m, -this%n:this%n)
        integer :: l

        do l = 1, 5
            ! print *, l
            ! print*, 'f_phy:'
            ! print*, f_phy(l, :, :)
            ! read(*, *)
            call this%fft_var(f_phy(l, :, :), f_spt(l, :, :))
            ! print*, 'f_spt:'
            ! print*, f_spt(l, :, :)
            ! read(*, *)
        end do

    end subroutine fft_phy_spt

    subroutine fft_var(this, cc, aa)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        real(real64), intent(inout) :: cc(0:this%msp - 1, 0:this%nsp - 1)
        complex(real64), intent(inout) :: aa(0:this%m, -this%n:this%n)

        integer :: msp, nsp, nma, nna
        complex(REAL64) :: bb(0:this%msp - 1, 0:this%nsp - 1)
        integer :: i, j, ii, jj

        msp = this%msp
        nsp = this%nsp

        nma = msp
        nna = nsp

        bb = cc

        call this%handle_new%fft(bb(0:msp - 1, 0:nsp - 1), bb(0:msp - 1, 0:nsp - 1))

        aa = 0.0d0

        ! print*, 'bb:'
        ! print*, bb
        ! read(*, *)
        ! aa(0, 0) = bb(0, 0)              ! m=0,n=0

        ! do i = 1, this%m
        !     ii = nma - i                 ! m=\0,n=0
        !     aa(i, 0) = bb(ii, 0)        !(-m,0)
        ! end do

        ! do j = 1, this%n                       ! m=0,n=\0
        !     jj = nna - j
        !     aa(0, j) = bb(0, j)                !(0, n)
        !     aa(0, -j) = conjg(bb(0, j))        !(0,-n) !!may be not used in -this%n:-1
        ! end do

        ! do i = 1, this%m                       !m<0
        !     ii = nma - i
        !     do j = 1, this%n
        !         jj = nna - j
        !         aa(i, j) = bb(ii, j)          !(-m,n)
        !         aa(i, -j) = bb(ii, jj)          !(-m,-n)
        !     end do
        ! end do
        i = 0
        aa(0, 0) = bb(0, 0)
        do j = 1, this%n
            jj = nna -j
            aa(0,  j) = bb(0, j )
            aa(0, -j) = bb(0, jj)
        enddo
        do i = 1, this%m
            ii = nma - i
            aa(i, 0) = bb(ii, 0)
            do j= 1, this%n
                jj = nna - j
                aa(i, j) = bb(ii, j)
                aa(i, -j) = bb(ii, jj)
            enddo
        enddo

    end subroutine fft_var

    subroutine compute_force(this, jloc, point_physics, force)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        integer, intent(in) :: jloc
        type(type_point_physics), intent(in) :: point_physics
        real(REAL64), intent(inout) :: force(:, :, :)

        real(REAL64) :: gmf, gmf1, gm1m2, pe, paraR
        real(REAL64) :: hx, hz, hx_y, hz_y, d12, d32, d132
        real(real64), parameter :: para4b3 = 4.0d0/3.0d0
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: divu, divux, divuy, divuz
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: s11, s22, s33, s12, s13, s23
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: s11x, s12y, s13z
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: s12x, s22y, s23z
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: s13x, s23y, s33z
        real(REAL64), dimension(0:this%msp - 1, 0:this%nsp - 1) :: miu, miux, miuy, miuz
        real(REAL64) :: cm, miu0, miu0T, miu0TT, miu0Tx, miu0Ty, miu0Tz
        real(REAL64) :: s110, s220, s330, s120, s130, s230

        real(REAL64) :: omv, Mx

        gmf = -1.d0/this%Gamma
        gmf1 = gmf + 1.0d0
        gm1M2 = ((this%Gamma - 1.d0)*this%Ma*this%Ma)
        pe = 1.d0/(this%Gamma*this%Ma*this%Ma)
        ParaR = gm1M2*Para4b3

        associate (drx => this%grid%drx, drz => this%grid%drz, y => this%grid%y, &
              &rho => point_physics%rho_a0, rhox => point_physics%rho_x, rhoy => point_physics%rho_y, rhoz => point_physics%rho_z, &
                 & u => point_physics%u_a0, v => point_physics%v_a0, w => point_physics%w_a0, &
                 & ut => point_physics%u_t, vt => point_physics%v_t, wt => point_physics%w_t, &
                 & rhot => point_physics%rho_t, Tt => point_physics%T_t, &
                 & ux => point_physics%u_x, uy => point_physics%u_y, uz => point_physics%u_z, &
                 & vx => point_physics%v_x, vy => point_physics%v_y, vz => point_physics%v_z, &
                 & wx => point_physics%w_x, wy => point_physics%w_y, wz => point_physics%w_z, &
                 & uxx => point_physics%u_xx, vyy => point_physics%v_yy, wzz => point_physics%w_zz, &
                 & uyy => point_physics%u_yy, uzz => point_physics%u_zz, &
                 & vxx => point_physics%v_xx, vzz => point_physics%v_zz, &
                 & wxx => point_physics%w_xx, wyy => point_physics%w_yy, &
                 & uxy => point_physics%u_xy, vxy => point_physics%v_xy, &
                 & uxz => point_physics%u_xz, wxz => point_physics%w_xz, &
                 & vyz => point_physics%v_yz, wyz => point_physics%w_yz, &
                 & T => point_physics%T_a0, Tx => point_physics%T_x, Ty => point_physics%T_y, Tz => point_physics%T_z, &
                 & Txx => point_physics%T_xx, Tyy => point_physics%T_yy, Tzz => point_physics%T_zz, &
                 & u0 => this%baseflow%q0(2, jloc), u0x => this%baseflow%dq0x(2, jloc), u0y => this%baseflow%dq0y(2, jloc), &
                 & v0 => this%baseflow%q0(3, jloc), v0x => this%baseflow%dq0x(3, jloc), v0y => this%baseflow%dq0y(3, jloc), &
                 & w0 => this%baseflow%q0(4, jloc), w0x => this%baseflow%dq0x(4, jloc), w0y => this%baseflow%dq0y(4, jloc), &
                 & T0 => this%baseflow%q0(5, jloc), T0x => this%baseflow%dq0x(5, jloc), T0y => this%baseflow%dq0y(5, jloc), &
                 & rho0 => this%baseflow%q0(1, jloc), rho0x => this%baseflow%dq0x(1, jloc), rho0y => this%baseflow%dq0y(1, jloc))

            select case (this%Coordinate_System)
            case (COMOVING_CS)
                hx = 1.d0 + drx*y(jloc)
                hx_y = drx
                d12 = hx_y/hx

                hz = 1.d0 + drz*y(jloc)
                hz_y = drz
                d32 = hz_y/hz
            case (CYLINDRICAL_CS)
                hx = 1.d0
                hx_y = 0.0d0
                d12 = 0.0d0

                hz = y(jloc)
                hz_y = 1.0d0
                d32 = hz_y/hz
            end select

            d132 = d12 + d32

            divu = ux/hx + vy + wz/hz + v*d132
            divux = uxx/hx/hx + vxy/hx + wxz/hx/hz + vx/hx*d132
            divuy = uxy/hx - d12*ux/hx + vyy + wyz/hz - d32*wz/hz + vy*d132 - v*(d12**2 + d32**2)
            divuz = uxz/hx/hz + vyz + vz/hz*d132 + wzz/hz/hz

            s11 = ux/hx + v*d12
            s22 = vy
            s33 = wz/hz + v*d32
            s12 = 0.5d0*(vx/hx - u*d12 + uy)
            s13 = 0.5d0*(wx/hx + uz/hz)
            s23 = 0.5d0*(wy + vz/hz - w*d32)

            s11x = uxx/hx/hx + vx/hx*d12
            s12y = 0.5d0*(vxy/hx - d12*vx/hx + uyy - uy*d12 - u*d12**2)
            s13z = 0.5d0*(wxz/hx/hz + uzz/hz/hz)

            s12x = 0.5d0*(vxx/hx/hx + uxy/hx - ux/hx*d12)
            s22y = vyy
            s23z = 0.5d0*(wyz/hz + vzz/hz/hz - wz/hz*d32)

            s13x = 0.5d0*(wxx/hx/hx + uxz/hx/hz)
            s23y = 0.5d0*(wyy + vyz/hz - d32*vz/hz - wy*d32 - w*d32**2)
            s33z = wzz/hz/hz + vz/hz*d32

            cm = this%C1/this%Te
            miu0 = T0*sqrt(T0)*(1.0d0 + cm)/(T0 + cm)/this%Re
            miu0T = miu0*(1.5d0/T0 - 1.0d0/(T0 + cm))
            miu0TT = miu0T*(1.5d0/T0 - 1.0d0/(T0 + cm)) - miu0*(1.5d0/T0**2 - 1.0d0/(T0 + cm)**2)
            miu0Tx = miu0TT*T0x
            miu0Ty = miu0TT*T0y
            miu0Tz = miu0TT*0.0d0
            miu = miu0T*T
            miux = miu0Tx/hx*T + miu0T*Tx/hx
            miuy = miu0Ty*T + miu0T*Ty
            miuz = miu0Tz/hz*T + miu0T*Tz/hz

            s110 = u0x/hx + v0*d12
            s220 = v0y
            s330 = v0*d32
            s120 = 0.5d0*(v0x/hx - d12*u0 + u0y)
            s130 = 0.5d0*(w0x/hx)
            s230 = 0.5d0*(w0y - w0*d32)

            Mx = u0/sqrt(T0)*this%Ma
            omv = min(1.0d0, Mx**2)

            ! print*, 1,ux(0, 0), vy(0, 0), wz(0, 0)
            ! print*, 2,divu(0, 0)
            ! print*, 3,divux(0, 0), divuy(0, 0), divuz(0, 0)
            ! print*, 4,s11(0, 0), s22(0, 0), s33(0, 0)
            ! print*, 5,s12(0, 0), s13(0, 0), s23(0, 0)
            ! print *,6, s11x(0, 0), s12y(0, 0), s13z(0, 0)
            ! print*, uy(0, 0)
            ! print*, uxx(0, 0), vxx(0, 0), wxx(0, 0)
            ! print*, Txx(0, 0)
            ! print*, uyy(0, 0), vyy(0, 0), wyy(0, 0)
            ! print*, Tyy(0, 0)
            ! print*, uzz(0, 0), vzz(0, 0), wzz(0, 0)
            ! print*, Tzz(0, 0)
            ! print*, vxx(0, 0), uxy(0, 0), s12x(0, 0)
            ! print*, 0.5d0*(vxx(0, 0)/hx/hx - uxy(0, 0)/hx)
            ! print*, ux(0, 0)/hx*d12

            ! print *,s12x(0, 0)!, s22y(0, 0), s23z(0, 0)
            ! print *,8, s13x(0, 0), s23y(0, 0), s33z(0, 0)
            ! print *,9, miu(0, 0)
            ! print *,10, miux(0, 0), miuy(0, 0), miuz(0, 0)
            ! print *,11, s110!, s220, s330
            ! print *,12, s120, s130, s230
            ! print *,13, miu0, miu0T, miu0TT
            ! print *,14, miu0Tx!, miu0Ty

            ! block !! check variables
            !     print *, 'check variables bf'
            !     print *, miu0, miu0T, miu0Tx, miu0Ty, miu0Tz
            !     print *, 'check medium variables depend'
            !     print *, hx, hz, d12, d32
            !     print *, v(0, 0)
            !     print *, ux(0, 0), vy(0, 0), wz(0, 0)
            !     print *, uxy(0, 0), vyy(0, 0), wyz(0, 0)
            !     print *, 'check medium variables'
            !     print *, divu(0, 0), divux(0, 0), divuy(0, 0), divuz(0, 0)
            !     print *, s11(0, 0), s22(0, 0), s33(0, 0)
            !     print *, s12(0, 0), s13(0, 0), s23(0, 0)
            !     print *, s11x(0, 0), s12y(0, 0), s13z(0, 0)
            !     print *, s12x(0, 0), s22y(0, 0), s23z(0, 0)
            !     print *, s13x(0, 0), s23y(0, 0), s33z(0, 0)
            !     print *, miu(0, 0), miux(0, 0), miuy(0, 0), miuz(0, 0)

            !     print *, 'check variables dis'
            !     print *, rho(0, 0), u(0, 0), v(0, 0), w(0, 0), T(0, 0)
            !     print *, rhot(0, 0), ut(0, 0), vt(0, 0), wt(0, 0), Tt(0, 0)
            !     print *, rhox(0, 0), ux(0, 0), vx(0, 0), wx(0, 0), Tx(0, 0)
            !     print *, rhoy(0, 0), uy(0, 0), vy(0, 0), wy(0, 0), Ty(0, 0)
            !     print *, rhoz(0, 0), uz(0, 0), vz(0, 0), wz(0, 0), Tz(0, 0)
            !     print *, uxx(0, 0), uyy(0, 0), uzz(0, 0)
            !     print *, vxx(0, 0), vyy(0, 0), vzz(0, 0)
            !     print *, wxx(0, 0), wyy(0, 0), wzz(0, 0)
            !     print *, Txx(0, 0), Tyy(0, 0), Tzz(0, 0)
            !     ! read (*, *)
            ! end block

            ! print *, 'omv=', omv
            ! read (*, *)

            force(1, :, :) = -(rho*ux/hx + u*rhox/hx + rho*vy + v*rhoy + rho*wz/hz + w*rhoz/hz + &
                           &  rho*v*d132)

            force(2, :, :) = -rho*(ut + u0*ux/hx + v0*uy + w0*uz/hz + u0*v*d12 + &
                           & +u0x/hx*u + u0y*v + d12*v0*u + u*ux/hx + v*uy + w*uz/hz + u*v*d12) &
                           & - rho0*(u*ux/hx + v*uy + w*uz/hz + u*v*d12) &
                           & - pe*(T*rhox/hx + rho*Tx/hx)*omv &
                           & + 2.d0*(miux*s11 + miu*s11x + miuy*s12 + miu*s12y + miuz*s13 + miu*s13z &
                           & + miu*(2.d0*d12 + d32)*s12) &
                           & - 2.d0/3.d0*(miux*divu + miu*divux)
            ! print *, "force2="
            ! print *, force(2, 1, 1)
            ! read (*, *)
            force(3, :, :) = -rho*(vt + u0*vx/hx + v0*vy + w0*vz/hz - u0*u*d12 - w0*w*d32 &
                           & + u*v0x/hx + v*v0y - d12*u0*u - d32*w0*w &
                           & + u*vx/hx + v*vy + w*vz/hz - u*u*d12 - w*w*d32) &
                           & - rho0*(u*vx/hx + v*vy + w*vz/hz - u*u*d12 - w*w*d32) &
                           & - pe*(rho*Ty + T*rhoy) +&
                           & +2.d0*(miux*s12 + miu*s12x + miuy*s22 + miu*s22y + miuz*s23 + miu*s23z &
                           & + miu*s22*d12 + miu*s22*d32 - miu*s11*d12 - miu*s33*d32) &
                           & - 2.d0/3.d0*(miuy*divu + miu*divuy)

            force(4, :, :) = -rho*(wt + u0*wx/hx + v0*wy + w0*wz/hz + w0*v*d32 &
                           & + u*w0x/hx + v*w0y + d32*v0*w &
                           & + u*wx/hx + v*wy + w*wz/hz + w*v*d32) &
                           & - rho0*(u*wx/hx + v*wy + w*wz/hz + w*v*d32) &
                           & - pe*(rho*Tz/hz + T*rhoz/hz) &
                           & + 2.d0*(miux*s13 + miu*s13x + miuy*s23 + miu*s23y + miuz*s33 + miu*s33z &
                           & + miu*s23*d12 + +miu*s23*d32 + miu*s23*d32) &
                           & - 2.d0/3.d0*(miuz*divu + miu*divuz)

            force(5, :, :) = rho*gmf*(Tt + (u*T0x/hx + v*T0y) + (u0*Tx/hx + v0*Ty + w0*Tz/hz) &
                           & + (u*Tx/hx + v*Ty + w*Tz/hz)) + rho0*gmf*(u*Tx/hx + v*Ty + w*Tz/hz) &
                           & + gmf1*T*(rhot + (u*rho0x/hx + v*rho0y) + (u0*rhox/hx + v0*rhoy + w0*rhoz/hz) &
                           & + (u*rhox/hx + v*rhoy + w*rhoz/hz)) + gmf1*T0*(u*rhox/hx + v*rhoy + w*rhoz/hz) &
                           & + (miux*Tx/hx + miuy*Ty + miuz*Tz/hz + miu*(Txx/hx/hx + Tyy + Tzz/hz/hz) &
                           & + miu*(d12*Ty + d32*Ty))/this%pr &
                           & + gm1M2*(4.d0*miu*(s110*s11 + s220*s22 + s330*s33 + 2.d0*(s120*s12 + s130*s13 &
                           & + s230*s23)) - 4.d0/3.d0*miu*(u0x/hx + v0y + +v0*(d12 + d32))*divu) &
                           & + gm1M2*(miu + miu0)*(2.d0*(s11*s11 + s22*s22 + s33*s33 &
                           &  + 2.d0*(s12*s12 + s13*s13 + s23*s23)) - 2.d0/3.d0*divu*divu)

            force(5, :, :) = force(5, :, :)/((this%GAMMA - 1.0d0)*this%Ma**2)
        end associate

    end subroutine

    subroutine ifft_spt_phy(this, point_spt, point_phy)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        type(type_point_spectrum), intent(in) :: point_spt
        type(type_point_physics), intent(inout) :: point_phy

        call this%ifft_var(point_spt%rho_a0, point_phy%rho_a0)
        call this%ifft_var(point_spt%u_a0, point_phy%u_a0)
        ! print*, 'in:', point_spt%v_a0(1, 0)
        ! print*, 'in done============'
        call this%ifft_var(point_spt%v_a0, point_phy%v_a0)
        ! print*, 'out:', point_phy%v_a0(0, 0)
        ! print*, 'out done============'
        call this%ifft_var(point_spt%w_a0, point_phy%w_a0)
        call this%ifft_var(point_spt%T_a0, point_phy%T_a0)
        call this%ifft_var(point_spt%rho_t, point_phy%rho_t)
        call this%ifft_var(point_spt%u_t, point_phy%u_t)
        call this%ifft_var(point_spt%v_t, point_phy%v_t)
        call this%ifft_var(point_spt%w_t, point_phy%w_t)
        call this%ifft_var(point_spt%T_t, point_phy%T_t)
        call this%ifft_var(point_spt%rho_x, point_phy%rho_x)
        call this%ifft_var(point_spt%u_x, point_phy%u_x)
        call this%ifft_var(point_spt%v_x, point_phy%v_x)
        call this%ifft_var(point_spt%w_x, point_phy%w_x)
        call this%ifft_var(point_spt%T_x, point_phy%T_x)
        call this%ifft_var(point_spt%rho_y, point_phy%rho_y)
        call this%ifft_var(point_spt%u_y, point_phy%u_y)
        call this%ifft_var(point_spt%v_y, point_phy%v_y)
        call this%ifft_var(point_spt%w_y, point_phy%w_y)
        call this%ifft_var(point_spt%T_y, point_phy%T_y)
        call this%ifft_var(point_spt%rho_z, point_phy%rho_z)
        call this%ifft_var(point_spt%u_z, point_phy%u_z)
        call this%ifft_var(point_spt%v_z, point_phy%v_z)
        call this%ifft_var(point_spt%w_z, point_phy%w_z)
        call this%ifft_var(point_spt%T_z, point_phy%T_z)
        call this%ifft_var(point_spt%u_xx, point_phy%u_xx)
        call this%ifft_var(point_spt%v_xx, point_phy%v_xx)
        call this%ifft_var(point_spt%w_xx, point_phy%w_xx)
        call this%ifft_var(point_spt%T_xx, point_phy%T_xx)
        call this%ifft_var(point_spt%u_yy, point_phy%u_yy)
        call this%ifft_var(point_spt%v_yy, point_phy%v_yy)
        call this%ifft_var(point_spt%w_yy, point_phy%w_yy)
        call this%ifft_var(point_spt%T_yy, point_phy%T_yy)
        call this%ifft_var(point_spt%u_zz, point_phy%u_zz)
        call this%ifft_var(point_spt%v_zz, point_phy%v_zz)
        call this%ifft_var(point_spt%w_zz, point_phy%w_zz)
        call this%ifft_var(point_spt%T_zz, point_phy%T_zz)
        call this%ifft_var(point_spt%u_xy, point_phy%u_xy)
        call this%ifft_var(point_spt%v_xy, point_phy%v_xy)
        call this%ifft_var(point_spt%u_xz, point_phy%u_xz)
        call this%ifft_var(point_spt%w_xz, point_phy%w_xz)
        call this%ifft_var(point_spt%v_yz, point_phy%v_yz)
        call this%ifft_var(point_spt%w_yz, point_phy%w_yz)

    end subroutine ifft_spt_phy

    subroutine ifft_var(this, aa, cc)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        complex(REAL64), intent(in) :: aa(0:this%m, -this%n:this%n)
        real(REAL64), intent(inout) :: cc(0:this%msp - 1, 0:this%nsp - 1)
        complex(REAL64) :: bb(0:this%msp - 1, 0:this%nsp - 1)
        integer :: nma, nna, msp, nsp

        integer :: shapeaa(2)

        integer :: i, ii, j, jj

        msp = this%msp
        nsp = this%nsp
        nma = msp
        nna = nsp

        bb = 0.0d0

        bb(0, 0) = aa(0, 0) !m=0, n=0
        do i = 1, this%m
            ii = nma - i
            bb(ii, 0) = aa(i, 0) ! (-m, 0)
            bb(i, 0) = conjg(aa(i, 0)) !(-m, 0)
        end do

        do j = 1, this%n
            jj = nna - j
            bb(0, j) = aa(0, j)  !(0, n)
            bb(0, jj) = conjg(aa(0, j)) !(n, -n)  !!! @TODO this line may be not used for -n
        end do

        !! m<0, n>0
        do i = 1, this%m
            ii = nma - i
            do j = 1, this%n
                jj = nna - j
                bb(ii, j) = aa(i, j) !(-m, n)
                bb(ii, jj) = aa(i, -j) !(m, -n)
            end do
        end do

        !! m<0, n<0
        do i = 1, this%m
            ii = nma - i
            do j = 1, this%n
                jj = nna - j
                bb(i, jj) = conjg(bb(ii, j))  !(-m, -n)
                bb(i, j) = conjg(bb(ii, jj)) !(m, n)
            end do
        end do

        !call SpectralBackward(this%handle, bb(0:msp - 1, 0:nsp - 1))
        ! print*, bb(1, 0)
        call this%handle_new%ifft(bb(0:msp - 1, 0:nsp - 1), bb(0:msp - 1, 0:nsp - 1))
        cc = real(bb, REAL64)
        ! print*, cc(0, 0)

    end subroutine ifft_var

    subroutine compute_dis_profile(this)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this

        integer :: m, n, j, l

        associate (profile => this%dis_profile_spectrum, exp_ialf => this%exp_integral_alf, &
                   phi0 => this%phi0_array, phi1 => this%phi1_array, phi2 => this%phi2_array, &
                   xi_x => this%grid%xi_x, eta_x => this%grid%eta_x, &
                   xi_y => this%grid%xi_y, eta_y => this%grid%eta_y, &
                   xi_xx => this%grid%xi_xx, eta_xx => this%grid%eta_xx, &
                   xi_xy => this%grid%xi_xy, eta_xy => this%grid%eta_xy, &
                   xi_yy => this%grid%xi_yy, eta_yy => this%grid%eta_yy, &
                   alf => this%alf, &
                   ome1 => this%ome1, bat1 => this%bat1, &
                   ome2 => this%ome2, bat2 => this%bat2)
            !   do n = -this%n, this%n
            !       do m = 0, this%m
            !           do j = 1, this%jn

            !! set profile phi_a0
            do concurrent(j=1:this%jn, m=0:this%m, n=-this%n:this%n)
                profile(j)%rho_a0(m, n) = phi0(1, j, m, n)*exp_ialf(m, n)
                profile(j)%u_a0(m, n) = phi0(2, j, m, n)*exp_ialf(m, n)
                profile(j)%v_a0(m, n) = phi0(3, j, m, n)*exp_ialf(m, n)
                profile(j)%w_a0(m, n) = phi0(4, j, m, n)*exp_ialf(m, n)
                profile(j)%T_a0(m, n) = phi0(5, j, m, n)*exp_ialf(m, n)
            end do
            !       end do
            !   end do

            block
                use mod_dx

            !!D phi /D xi
                do concurrent(l=1:5, m=0:this%m, n=-this%n:this%n, j=1:this%jn)
                    dphidi(l, j, m, n) = (phi0(l, j, m, n)*(1.0d0 + this%lambda) + &
                                        & phi1(l, j, m, n)*(-1.0d0 - 2.0d0*this%lambda) + &
                                        & phi2(l, j, m, n)*(this%lambda))
                end do
            !!D phi /D eta, D^2phi /D eta^2
                !  do n = -this%n, this%n
                !      do m = 0, this%m
                !          do l = 1, 5
                do concurrent(l=1:5, m=0:this%m, n=-this%n:this%n)
                    call dx_fdm4nd_center(phi0(l, :, m, n), 1, this%jn, dphidj(l, :, m, n), 1, this%jn)
                    call dxx_fdm4nd_center(phi0(l, :, m, n), 1, this%jn, dphidjj(l, :, m, n), 1, this%jn)
                end do
                !          end do
                !      end do
                !  end do
            !! D^2 phi /D eta D xi
                !  do n = -this%n, this%n
                !      do m = 0, this%m
                !          do l = 1, 5
                do concurrent(l=1:5, m=0:this%m, n=-this%n:this%n)
                    call dx_fdm4nd_center(dphidi(l, :, m, n), 1, this%jn, dphidij(l, :, m, n), 1, this%jn)
                end do
                !          end do
                !      end do
                !  end do

                !  do n = -this%n, this%n
                !      do m = 0, this%m
                !          do j = 1, this%jn
                do concurrent(j=1:this%jn, m=0:this%m, n=-this%n:this%n)
                    dphidx(:, j, m, n) = dphidi(:, j, m, n)*this%grid%xi_x(j) + dphidj(:, j, m, n)*eta_x(j)
                    dphidy(:, j, m, n) = dphidi(:, j, m, n)*this%grid%xi_y(j) + dphidj(:, j, m, n)*eta_y(j)

                    dphidxy(:, j, m, n) = dphidi(:, j, m, n)*xi_xy(j) + dphidj(:, j, m, n)*eta_xy(j) + &
                                        & dphidjj(:, j, m, n)*(eta_y(j)*eta_x(j)) + &
                                        & dphidij(:, j, m, n)*(xi_x(j)*eta_y(j) + xi_y(j)*eta_x(j))
                    dphidyy(:, j, m, n) = dphidi(:, j, m, n)*xi_yy(j) + dphidj(:, j, m, n)*eta_yy(j) + &
                                        & dphidjj(:, j, m, n)*(eta_y(j)*eta_y(j)) + &
                                        & dphidij(:, j, m, n)*(2.0d0*xi_y(j)*eta_y(j))
                    dphidxx(:, j, m, n) = dphidi(:, j, m, n)*xi_xx(j) + dphidj(:, j, m, n)*eta_xx(j) + &
                                        & dphidjj(:, j, m, n)*(eta_x(j)*eta_x(j)) + &
                                        & dphidij(:, j, m, n)*(2.0d0*xi_x(j)*eta_x(j))
                end do
                !          end do
                !      end do
                !  end do
                ! block
                !     print *, 'check dis profile pre compute'
                !     print *, dphidjj(2, 1, 1, 0), dphidj(2, 1, 1, 0), dphidyy(2, 1, 1, 0)
                !     print *, dphidjj(3, 1, 1, 0), dphidj(3, 1, 1, 0), dphidyy(3, 1, 1, 0)
                !     print *, dphidjj(4, 1, 1, 0), dphidj(4, 1, 1, 0), dphidyy(4, 1, 1, 0)
                !     print *, dphidjj(5, 1, 1, 0), dphidj(5, 1, 1, 0), dphidyy(5, 1, 1, 0)
                !     ! read(*, *)
                ! end block

                do concurrent(j=1:this%jn, m=0:this%m, n=-this%n:this%n)
                !! dphi'dx_spt, dphi'dxx_spt
                    profile(j)%rho_x(m, n) = dphidx(1, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%rho_a0(m, n)
                    profile(j)%u_x(m, n) = dphidx(2, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%u_a0(m, n)
                    profile(j)%v_x(m, n) = dphidx(3, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%v_a0(m, n)
                    profile(j)%w_x(m, n) = dphidx(4, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%w_a0(m, n)
                    profile(j)%T_x(m, n) = dphidx(5, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%T_a0(m, n)
                    ! auxx(j, m, n) = aua0(j, m, n)*aia*aia + 2.d0*aia*aux(j, m, n)*a0 + ai*alfx(m, n)*aua0(j, m, n)
                    profile(j)%u_xx(m, n) = profile(j)%u_a0(m, n)*(ZI*alf(m, n))**2 + &
                                          & 2.0d0*(ZI*alf(m, n))*dphidx(2, j, m, n)*exp_ialf(m, n)
                    profile(j)%v_xx(m, n) = profile(j)%v_a0(m, n)*(ZI*alf(m, n))**2 + &
                    & 2.0d0*(ZI*alf(m, n))*dphidx(3, j, m, n)*exp_ialf(m, n)
                    profile(j)%w_xx(m, n) = profile(j)%w_a0(m, n)*(ZI*alf(m, n))**2 + &
                    & 2.0d0*(ZI*alf(m, n))*dphidx(4, j, m, n)*exp_ialf(m, n)
                    profile(j)%T_xx(m, n) = profile(j)%T_a0(m, n)*(ZI*alf(m, n))**2 + &
                    & 2.0d0*(ZI*alf(m, n))*dphidx(5, j, m, n)*exp_ialf(m, n)
                end do

                !! dphi'dt_spt, dphi'dz_spt
                do concurrent(j=1:this%jn, m=0:this%m, n=-this%n:this%n)
                    profile(j)%rho_t(m, n) = -ZI*(ome1*real(m, REAL64)+ome2*real(n, REAL64))*profile(j)%rho_a0(m, n)
                    profile(j)%u_t(m, n) = -ZI*(ome1*real(m, REAL64)+ome2*real(n, REAL64))*profile(j)%u_a0(m, n)
                    profile(j)%v_t(m, n) = -ZI*(ome1*real(m, REAL64)+ome2*real(n, REAL64))*profile(j)%v_a0(m, n)
                    profile(j)%w_t(m, n) = -ZI*(ome1*real(m, REAL64)+ome2*real(n, REAL64))*profile(j)%w_a0(m, n)
                    profile(j)%T_t(m, n) = -ZI*(ome1*real(m, REAL64)+ome2*real(n, REAL64))*profile(j)%T_a0(m, n)
                    profile(j)%rho_z(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%rho_a0(m, n)
                    profile(j)%u_z(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%u_a0(m, n)
                    profile(j)%v_z(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%v_a0(m, n)
                    profile(j)%w_z(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%w_a0(m, n)
                    profile(j)%T_z(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%T_a0(m, n)

                    profile(j)%u_zz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%u_z(m, n)
                    profile(j)%v_zz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%v_z(m, n)
                    profile(j)%w_zz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%w_z(m, n)
                    profile(j)%T_zz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%T_z(m, n)
                    !! dphi'dxz
                    profile(j)%u_xz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%u_x(m, n)
                    profile(j)%w_xz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%w_x(m, n)
                end do

                do concurrent(j=1:this%jn, m=0:this%m, n=-this%n:this%n)
                   !! dphi'dy_spt, d^2phi'dyy_spt
                    profile(j)%rho_y(m, n) = dphidy(1, j, m, n)*exp_ialf(m, n)
                    profile(j)%u_y(m, n) = dphidy(2, j, m, n)*exp_ialf(m, n)
                    profile(j)%v_y(m, n) = dphidy(3, j, m, n)*exp_ialf(m, n)
                    profile(j)%w_y(m, n) = dphidy(4, j, m, n)*exp_ialf(m, n)
                    profile(j)%T_y(m, n) = dphidy(5, j, m, n)*exp_ialf(m, n)
                    profile(j)%u_yy(m, n) = dphidyy(2, j, m, n)*exp_ialf(m, n)
                    profile(j)%v_yy(m, n) = dphidyy(3, j, m, n)*exp_ialf(m, n)
                    profile(j)%w_yy(m, n) = dphidyy(4, j, m, n)*exp_ialf(m, n)
                    profile(j)%T_yy(m, n) = dphidyy(5, j, m, n)*exp_ialf(m, n)
                    !! dphi'dxy_spt, dphi'dyz_spt
                    profile(j)%u_xy(m, n) = dphidxy(2, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%u_y(m, n)
                    profile(j)%v_xy(m, n) = dphidxy(3, j, m, n)*exp_ialf(m, n) + ZI*alf(m, n)*profile(j)%v_y(m, n)
                    profile(j)%v_yz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%v_y(m, n)
                    profile(j)%w_yz(m, n) = ZI*(bat1*real(m, REAL64)+bat2*real(n, REAL64))*profile(j)%w_y(m, n)
                end do
                ! block
                !     print *, 'check dis profile'
                !     print *, profile(2)%u_yy(1, 0), dphidyy(2, 1, 1, 0)
                !     print *, profile(2)%v_yy(1, 0), dphidyy(3, 1, 1, 0)
                !     print *, profile(2)%w_yy(1, 0), dphidyy(4, 1, 1, 0)
                !     print *, profile(2)%T_yy(1, 0), dphidyy(5, 1, 1, 0)
                !     ! read(*, *)
                ! end block
            end block

        end associate

    end subroutine compute_dis_profile

    subroutine finalize_nonliner_term(this)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        integer :: j

        do j = 1, this%jn
            this%dis_profile_spectrum(j)%m = 0
            this%dis_profile_spectrum(j)%n = 0
            deallocate (this%dis_profile_spectrum(j)%rho_a0)
            deallocate (this%dis_profile_spectrum(j)%u_a0)
            deallocate (this%dis_profile_spectrum(j)%v_a0)
            deallocate (this%dis_profile_spectrum(j)%w_a0)
            deallocate (this%dis_profile_spectrum(j)%T_a0)
            deallocate (this%dis_profile_spectrum(j)%rho_t)
            deallocate (this%dis_profile_spectrum(j)%u_t)
            deallocate (this%dis_profile_spectrum(j)%v_t)
            deallocate (this%dis_profile_spectrum(j)%w_t)
            deallocate (this%dis_profile_spectrum(j)%T_t)
            deallocate (this%dis_profile_spectrum(j)%rho_x)
            deallocate (this%dis_profile_spectrum(j)%u_x)
            deallocate (this%dis_profile_spectrum(j)%v_x)
            deallocate (this%dis_profile_spectrum(j)%w_x)
            deallocate (this%dis_profile_spectrum(j)%T_x)
            deallocate (this%dis_profile_spectrum(j)%rho_y)
            deallocate (this%dis_profile_spectrum(j)%u_y)
            deallocate (this%dis_profile_spectrum(j)%v_y)
            deallocate (this%dis_profile_spectrum(j)%w_y)
            deallocate (this%dis_profile_spectrum(j)%T_y)
            deallocate (this%dis_profile_spectrum(j)%rho_z)
            deallocate (this%dis_profile_spectrum(j)%u_z)
            deallocate (this%dis_profile_spectrum(j)%v_z)
            deallocate (this%dis_profile_spectrum(j)%w_z)
            deallocate (this%dis_profile_spectrum(j)%T_z)
            deallocate (this%dis_profile_spectrum(j)%u_xx)
            deallocate (this%dis_profile_spectrum(j)%v_xx)
            deallocate (this%dis_profile_spectrum(j)%w_xx)
            deallocate (this%dis_profile_spectrum(j)%T_xx)
            deallocate (this%dis_profile_spectrum(j)%u_yy)
            deallocate (this%dis_profile_spectrum(j)%v_yy)
            deallocate (this%dis_profile_spectrum(j)%w_yy)
            deallocate (this%dis_profile_spectrum(j)%T_yy)
            deallocate (this%dis_profile_spectrum(j)%u_zz)
            deallocate (this%dis_profile_spectrum(j)%v_zz)
            deallocate (this%dis_profile_spectrum(j)%w_zz)
            deallocate (this%dis_profile_spectrum(j)%T_zz)
            deallocate (this%dis_profile_spectrum(j)%u_xy)
            deallocate (this%dis_profile_spectrum(j)%v_xy)
            deallocate (this%dis_profile_spectrum(j)%u_xz)
            deallocate (this%dis_profile_spectrum(j)%w_xz)
            deallocate (this%dis_profile_spectrum(j)%v_yz)
            deallocate (this%dis_profile_spectrum(j)%w_yz)
        end do

        deallocate (point_physics%rho_a0)
        deallocate (point_physics%u_a0)
        deallocate (point_physics%v_a0)
        deallocate (point_physics%w_a0)
        deallocate (point_physics%T_a0)
        deallocate (point_physics%rho_t)
        deallocate (point_physics%u_t)
        deallocate (point_physics%v_t)
        deallocate (point_physics%w_t)
        deallocate (point_physics%T_t)
        deallocate (point_physics%rho_x)
        deallocate (point_physics%u_x)
        deallocate (point_physics%v_x)
        deallocate (point_physics%w_x)
        deallocate (point_physics%T_x)
        deallocate (point_physics%rho_y)
        deallocate (point_physics%u_y)
        deallocate (point_physics%v_y)
        deallocate (point_physics%w_y)
        deallocate (point_physics%T_y)
        deallocate (point_physics%rho_z)
        deallocate (point_physics%u_z)
        deallocate (point_physics%v_z)
        deallocate (point_physics%w_z)
        deallocate (point_physics%T_z)
        deallocate (point_physics%u_xx)
        deallocate (point_physics%v_xx)
        deallocate (point_physics%w_xx)
        deallocate (point_physics%T_xx)
        deallocate (point_physics%u_yy)
        deallocate (point_physics%v_yy)
        deallocate (point_physics%w_yy)
        deallocate (point_physics%T_yy)
        deallocate (point_physics%u_zz)
        deallocate (point_physics%v_zz)
        deallocate (point_physics%w_zz)
        deallocate (point_physics%T_zz)
        deallocate (point_physics%u_xy)
        deallocate (point_physics%v_xy)
        deallocate (point_physics%u_xz)
        deallocate (point_physics%w_xz)
        deallocate (point_physics%v_yz)
        deallocate (point_physics%w_yz)

        deallocate (this%alf)
        deallocate (this%exp_integral_alf)
        deallocate (this%f_spectrum)
        deallocate (this%f_physics)
        deallocate (this%dis_profile_spectrum)
        deallocate (this%phi0_array)
        deallocate (this%phi1_array)
        deallocate (this%phi2_array)
        call this%handle_new%Finalize()
        this%m = 0
        this%n = 0
        this%ome1 = 0
        this%bat1 = 0
        this%ome2 = 0
        this%bat2 = 0
        this%jn = 0

        deallocate(dphidi )
        deallocate(dphidj )
        deallocate(dphidij)
        deallocate(dphidjj)
        deallocate(dphidx )
        deallocate(dphidy )
        deallocate(dphidxx)
        deallocate(dphidxy)
        deallocate(dphidyy)

    end subroutine finalize_nonliner_term

    subroutine nonliner_term_set_parameters_real8(this, keyword, value)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        character(*), intent(in) :: keyword
        real(REAL64), intent(in) ::  value

        select case (trim(keyword))
        case ("Ma", "MA", "ma")
            this%Ma = value
        case ("Re", "RE", "re")
            this%Re = value
        case ("Te", "TE", "te")
            this%Te = value
        case ("GAMMA", "Gamma", "gamma")
            this%GAMMA = value
        case ("Pr", "pr", "PR")
            this%Pr = value
        case ("ome1")
            this%ome1 = value
        case ("bat1")
            this%bat1 = value
        case ("ome2")
            this%ome2 = value
        case ("bat2")
            this%bat2 = value            
        case default
            call error_stop("the value is wrong in set_parameters in nonlinear term.")
        end select

    end subroutine nonliner_term_set_parameters_real8

    subroutine nonliner_term_set_parameters_real4(this, keyword, value)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        character(*), intent(in) :: keyword
        real(REAL32), intent(in) ::  value

        select case (trim(keyword))
        case ("Ma", "MA", "ma")
            this%Ma = value
        case ("Re", "RE", "re")
            this%Re = value
        case ("Te", "TE", "te")
            this%Te = value
        case ("GAMMA", "Gamma", "gamma")
            this%GAMMA = value
        case ("Pr", "pr", "PR")
            this%Pr = value
        case ("ome1")
            this%ome1 = value
        case ("bat1")
            this%bat1 = value
        case ("ome2")
            this%ome2 = value
        case ("bat2")
            this%bat2 = value
        case default
            call error_stop("the value is wrong in set_parameters in nonlinear term.")
        end select

    end subroutine nonliner_term_set_parameters_real4

    subroutine nonliner_term_set_parameters_int(this, keyword, value)

        implicit none
        class(type_nonlinear_term), intent(inout) :: this
        character(*), intent(in) :: keyword
        integer, intent(in) ::  value

        select case (trim(keyword))
        case ("Coordinate_System")
            this%Coordinate_System = value
        case ("m")
            this%m = value
            if (this%m == 0) then
                this%msp = 1
            else
                this%msp = 4*(this%m + 1)
            end if
        case ("n")
            this%n = value
            if (this%n == 0) then
                this%nsp = 1
            else
                this%nsp = 4*(this%n + 1)
            end if
        case ("jn")
            this%jn = value
        case default
            call error_stop("the value is wrong in set_parameters in nonlinear term.")
        end select

    end subroutine nonliner_term_set_parameters_int

    function get_array(V, dm) result(array)

        Vec, intent(in) :: V
        DM, intent(in) :: dm
        PetscScalar, pointer :: tmp(:, :)
        PetscScalar, allocatable :: array(:, :)

        PetscErrorCode :: ierr

        call DMDAVecGetArrayReadF90(dm, V, tmp, ierr)
        allocate (array, source=tmp)
        call DMDAVecRestoreArrayReadF90(dm, V, tmp, ierr)

    end function get_array

end module mod_nonlinear_term

