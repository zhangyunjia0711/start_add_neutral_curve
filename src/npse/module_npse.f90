!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-12 08:16:46
!> @LastEditTime: 2022-12-12 19:30:37
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/npse/module_npse.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_npse

    ! use global_parameters
    use iso_fortran_env, only: REAL64
    use mod_lst_baseflow_profile
    use mod_lst_grid_profile

    use global_parameters, only: sigma_control, sigma_p
    use global_parameters, only: rd_factor
    use global_parameters, only: omega_force_factor, omega_phi_factor, ix_force, ix_phi
    use global_parameters, only: atol_nonlinear, rtol_nonlinear
    use global_parameters, only: maxit_alpha, maxit_nonlinear
    use global_parameters, only: cfl_number, ix_cfl

    use global_parameters, only: in, jn
    use global_parameters, only: xx
    use global_parameters, only: Re, Ma, Te, Pr
    use global_parameters, only: GAMMA
    use global_parameters, only: yb, ye!, xb, xe

    use global_parameters, only: istart, iend
    use global_parameters, only: m_global, n_global
    use global_parameters, only: theta1, theta2
    use global_parameters, only: inlet_dis
    use global_parameters, only: ndis
    use global_parameters, only: isothermal

    ! use global_parameters, only: phi
    ! use global_parameters, only: wave_param
    ! use global_parameters, only: lnA

    use global_parameters, only: meshDA
    use global_parameters, only: rank

    use petsc

    ! use mod_lpse_grid_profile, only: type_lpse_grid_domain=> type_lpse_grid_profile
    use mod_npse_eqn
    use mod_nonlinear_term

    implicit none

#include "../include/start.f90"

    ! complex(REAL64) :: ome, alf, bat
    DM :: DM_local
    type(type_npse_eqn) :: npse_eqn
    PetscErrorCode :: ierr
    ! type(string_type) :: solname

    complex(REAL64), allocatable :: alf(:, :, :)
    ! complex(REAL64) :: ome
    ! complex(REAL64) :: bat
    integer(kind(type_solve_alpha_method)), allocatable :: solve_alpha(:, :)

    ! type(type_waveparam), allocatable :: wave_param(:)
    complex(REAL64), allocatable :: lnA0(:, :, :)  !! 当前站位的 \int_0^x \alpha d\xi
    complex(REAL64), allocatable :: lnA1(:, :) !! 上一站位的 \int_0^x \alpha d\xi
    Vec, dimension(:), allocatable :: phi(:, :, :), phi_old(:, :)
    Vec :: delta_baseflow(3)
    complex(REAL64), allocatable :: nonlinear_term(:, :, :, :)
    complex(REAL64), allocatable :: nonlinear_term_old(:, :, :, :)
    real(REAL64), allocatable :: VecDistance(:, :)
    type(type_nonlinear_term) :: nonlinear_solver

    private
    public :: npse

    ! DM :: DM_LPSE

contains

    subroutine npse()

        implicit none

        ! solname = 'NPSE'
        !DM_local=>lpse_eqn%dm
        call PetscPrintf(PETSC_COMM_WORLD, 'Nonlinear PSE is selected!\n', ierr)

        call initialize_npse()

        call solve_npse_native()

        call finalize_npse()

    end subroutine npse

    subroutine initialize_npse()

        implicit none
        integer :: i, m, n

        call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_local, ierr)
        call DMSetUp(DM_local, ierr)

        allocate (phi(istart - 1:iend, 0:m_global, -n_global:n_global))
        allocate (phi_old(0:m_global, -n_global:n_global))
        allocate (alf(istart:iend, 0:m_global, -n_global:n_global)); alf= 0.0d0
        allocate (VecDistance(0:m_global, -n_global:n_global))
        allocate (nonlinear_term(5, jn, 0:m_global, -n_global:n_global)); nonlinear_term = 0.0d0
        allocate (nonlinear_term_old(5, jn, 0:m_global, -n_global:n_global)); nonlinear_term_old = 0.0d0

        call nonlinear_solver%set_parameter("m", m_global)
        call nonlinear_solver%set_parameter("n", n_global)
        call nonlinear_solver%set_parameter("ome1", theta1%re)
        call nonlinear_solver%set_parameter("ome2", theta2%re)
        call nonlinear_solver%set_parameter("bat1", theta1%im)
        call nonlinear_solver%set_parameter("bat2", theta2%im)
        call nonlinear_solver%initial(jn, Re, Ma, Te)

        allocate (solve_alpha(0:m_global, -n_global:n_global))
        do n = -n_global, n_global
            do m = 0, m_global
                do i = istart - 1, iend
                    call DMGetGlobalVector(DM_local, phi(i, m, n), ierr)
                    call VecZeroEntries(phi(i, m, n), ierr)
                end do
                call DMGetGlobalVector(DM_local, phi_old(m, n), ierr)
                call VecZeroEntries(phi_old(m, n), ierr)
            end do
        end do
        do i = 1, 3
            call DMGetGlobalVector(DM_local, delta_baseflow(i), ierr)
            call VecZeroEntries(delta_baseflow(i), ierr)
        end do

        allocate (lnA0(istart:iend, 0:m_global, -n_global:n_global)); lnA0 = 0.0d0
        allocate (lnA1(0:m_global, -n_global:n_global)); lnA1 = 0.0

        call npse_eqn%initial(jn, Re, Ma, Te, yb, ye)
        call npse_eqn%set_parameter("sigma_p", sigma_p)
        call npse_eqn%set_parameter("sigma_control", sigma_control)
        call npse_eqn%set_parameter("rd_factor", rd_factor)
        call npse_eqn%set_parameter("Pr", Pr)
        npse_eqn%isothermal = isothermal
        call npse_eqn%initial_force()

        solve_alpha = SOLVE_ALPHA_NONE
        do i = 1, ndis
            alf(IStart, inlet_dis(i)%m, inlet_dis(i)%n) = inlet_dis(i)%alf
            solve_alpha(inlet_dis(i)%m, inlet_dis(i)%n) = inlet_dis(i)%solve_alpha
        end do

        call update_alpha(istart, assemble_lst_grid_profile(istart))

        call set_npse_ic()

    end subroutine initialize_npse

    subroutine set_npse_ic()

        use mod_lst_solver_abstract
        use mod_lst_eqn

        implicit none
        integer :: i
        PetscViewer :: Viewer
        type(type_lst_eqn) :: lst
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        integer :: iloc
        complex(REAL64) :: ome, bat

        iloc = istart
    !!! compute from lst directly
        call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
        call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
        call lst%set_parameter("Pr", Pr)
        lst%isothermal = isothermal
        
        grid_profile = assemble_lst_grid_profile(iloc)
        baseflow_profile = assemble_lst_baseflow_profile(iloc)

        call lst%set_grid(grid_profile)
        call lst%set_baseflow(baseflow_profile)
        call lst%setup()

        do i = 1, ndis
            ome = theta1%re*inlet_dis(i)%m+theta2%re*inlet_dis(i)%n
            bat = theta1%im*inlet_dis(i)%m+theta2%im*inlet_dis(i)%n
            print*, i, ome%re, bat%re, alf(iloc, inlet_dis(i)%m, inlet_dis(i)%n)
            call lst%solve(ome, alf(iloc, inlet_dis(i)%m, inlet_dis(i)%n), bat)
            call lst%get_eigf(phi(iloc, inlet_dis(i)%m, inlet_dis(i)%n))
            call VecScale(phi(iloc, inlet_dis(i)%m, inlet_dis(i)%n), inlet_dis(i)%amp, ierr)
            print *, 'disturbance at inlet', iloc, inlet_dis(i)%m, inlet_dis(i)%n
            print *, 'drx=', grid_profile%drx, 'drz=', grid_profile%drz
            print *, 'alpha=', alf(iloc, inlet_dis(i)%m, inlet_dis(i)%n)%re, alf(iloc, inlet_dis(i)%m, inlet_dis(i)%n)%im
        end do

    end subroutine set_npse_ic

    subroutine solve_npse_native()

        use omp_lib
        implicit none
        integer :: i, m, n, l
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        PetscViewer :: Viewer
        integer, allocatable :: iterator_times(:, :)
        real(REAL64) :: omega_force, omega_phi
        integer :: itimes_nonlinear, itimes_alpha
        real(REAL64)::Norm_rhs, norm_rhs_mn(0:m_global, -n_global:n_global)
        real(REAL64) :: lambda
        type(type_lst_bf_profile) :: delta_baseflow_profile
        logical :: flg_mfd = .true.

        allocate (iterator_times(0:m_global, -n_global:n_global))
        do n = -n_global, n_global
            do m = 0, m_global
                call npse_eqn%set_ic(phi(istart, m, n))
            end do
        end do

        do i = istart + 1, iend

            grid_profile = assemble_lst_grid_profile(i)
            baseflow_profile = assemble_lst_baseflow_profile(i)

            !! generate baseflow distoration and use  basflow and its distoration generate the new baseflow
            !! meanwhile phi(0,0) is set zero because it is a part of base flow.

            if (i >= ix_force) then
                omega_force = omega_force_factor
            else
                omega_force = 1.0d0
            end if
            if (i >= ix_phi) then
                omega_phi = omega_phi_factor
            else
                omega_phi = 1.0d0
            end if

            if (i == istart + 1) then
                lambda = 0.0d0
            else
                lambda = 0.5d0
            end if

            ! call solve_lpse_eqn(grid_profile, baseflow_profile, ome, alf, bat, i, phi(i), phi(i-1), phi(max(1,i-2)))

            !! initial_ilocation_value
            do n = -n_global, n_global
                do m = 0, m_global
                    ! if (i == istart + 1) then
                    call VecCopy(phi(i - 1, m, n), phi(i, m, n), ierr)
                    alf(i, m, n) = alf(i - 1, m, n)
                    ! elseif (i == istart + 2) then
                    !     call VecWAXPY(phi(i, m, n), (-1.0d0, 0.0d0), phi(i - 2, m, n), phi(i - 1, m, n), ierr)
                    !     call VecAXPY(phi(i, m, n), (1.0d0, 0.0d0), phi(i - 1, m, n), ierr)
                    !     alf(i, m, n) = 2.0d0*alf(i - 1, m, n) - alf(i - 2, m, n)
                    ! elseif (i >= istart + 3) then
                    !     call VecWAXPY(phi(i, m, n), (-1.0d0, 0.0d0), phi(i - 3, m, n), phi(i - 2, m, n), ierr)
                    !     call VecScale(phi(i, m, n), (0.5d0, 0.0d0), ierr)
                    !     call VecAXPY(phi(i, m, n), (-1.5d0, 0.0d0), phi(i - 2, m, n), ierr)
                    !     call VecAXPY(phi(i, m, n), (2.5d0, 0.0d0), phi(i - 1, m, n), ierr)
                    !     alf(i, m, n) = 2.5d0*alf(i - 1, m, n) - 2.0d0*alf(i - 2, m, n) + 0.5d0*(alf(i - 3, m, n))
                    ! end if
                    ! alf(i, m, n) = alf(i - 1, m, n)
                    lnA0(i, m, n) = lnA1(m, n) + 0.5d0*(alf(i, m, n) + alf(i - 1, m, n))/grid_profile%xi_x(yb + 1)
                end do
            end do

            call npse_eqn%set_grid(grid_profile)
            call nonlinear_solver%set_grid(grid_profile)
            call npse_eqn%set_baseflow(baseflow_profile)
            call nonlinear_solver%set_baseflow(baseflow_profile)
            call npse_eqn%setup()

            ! if (flg_mfd) then
            !     call VecCopy(phi(i, 0, 0), delta_baseflow(1), ierr)
            !     call VecCopy(phi(i - 1, 0, 0), delta_baseflow(2), ierr)
            !     call VecCopy(phi(i - 2, 0, 0), delta_baseflow(3), ierr)
            !     delta_baseflow_profile = get_baseflow_profile(DM_Local, delta_baseflow, grid_profile, lambda)
            !     call npse_eqn%set_baseflow_npse(baseflow_profile + delta_baseflow_profile)
            ! else
            !     call npse_eqn%set_baseflow_npse(baseflow_profile)
            ! end if
            ! call npse_eqn%npse_setup()

            if (i < ix_cfl) then
                call npse_eqn%solve_dt(MAXNUM)
            else
                call npse_eqn%solve_dt(cfl_number)
            end if

            call update_alpha(i, grid_profile)
            lnA0(i, :, :) = lnA1(:, :) + 0.5d0*(alf(i, :, :) + alf(i - 1, :, :))/grid_profile%xi_x(yb + 1)
            call nonlinear_solver%set_dis(alf(i, :, :), lnA0(i, :, :), phi(i, :, :), phi(i - 1, :, :), phi(i - 2, :, :))
            call nonlinear_solver%setup()
            call nonlinear_solver%get_force(nonlinear_term)
            nonlinear_term_old = 0.0d0
            ! call nonlinear_solver%set_dis(alf(i, :, :), lnA0(i, :, :), phi(i, :, :), phi(i - 1, :, :), phi(i - 2, :, :))
            ! call nonlinear_solver%setup()
            ! call nonlinear_solver%get_force(nonlinear_term)
            ! nonlinear_term_old = nonlinear_term

            iterator_times = 1
            itimes_alpha = 0
            do while ((sum(iterator_times) > 0) .and. itimes_alpha < maxit_alpha)

                !!! Get wave number with nonlinear force
                lnA0(i, :, :) = lnA1(:, :) + 0.5d0*(alf(i, :, :) + alf(i - 1, :, :))/grid_profile%xi_x(yb + 1)
                nonlinear_term_old = nonlinear_term
                call nonlinear_solver%set_dis(alf(i, :, :), lnA0(i, :, :), phi(i, :, :), phi(i - 1, :, :), phi(i - 2, :, :))
                call nonlinear_solver%setup()
                call nonlinear_solver%get_force(nonlinear_term)

                nonlinear_term = omega_force*nonlinear_term + (1.0d0 - omega_force)*nonlinear_term_old

                iterator_times = 0
                do l = 1, ndis
                    m = inlet_dis(l)%m; n = inlet_dis(l)%n
                    call npse_eqn%set_parameter("ome", theta1%re*m+theta2%re*n)
                    call npse_eqn%set_parameter("bat", theta1%im*m+theta2%im*n)
                    call npse_eqn%set_parameter("solve_alpha", solve_alpha(m, n))

                    call npse_eqn%set_force(nonlinear_term(:, :, m, n))
                    call npse_eqn%solve(alf(i, m, n), phi(i, m, n), phi(i - 1, m, n), phi(i - 2, m, n))
                    iterator_times(m, n) = iterator_times(m, n) + npse_eqn%iterator_times
                end do

                !! update streamwise wave number alpha
                call update_alpha(i, grid_profile)
                lnA0(i, :, :) = lnA1(:, :) + 0.5d0*(alf(i, :, :) + alf(i - 1, :, :))/grid_profile%xi_x(yb + 1)

                !! solve the shape function with the computed wave number
                VecDistance = 1.0d0
                itimes_nonlinear = 0
                Norm_rhs = 0.0d10

                do while (((sqrt(sum(VecDistance**2)) >=  &
                            max(atol_nonlinear, rtol_nonlinear*Norm_rhs))) .and. &
                            ! .and. &
                            (itimes_nonlinear < maxit_nonlinear))

                    nonlinear_term_old = nonlinear_term
                    if (flg_mfd) then
                        call VecCopy(phi(i, 0, 0), delta_baseflow(1), ierr)
                        call VecCopy(phi(i - 1, 0, 0), delta_baseflow(2), ierr)
                        call VecCopy(phi(i - 2, 0, 0), delta_baseflow(3), ierr)
                        delta_baseflow_profile = get_baseflow_profile(DM_Local, delta_baseflow, grid_profile, lambda)
                        call npse_eqn%set_baseflow_npse(baseflow_profile + delta_baseflow_profile)
                    else
                        call npse_eqn%set_baseflow_npse(baseflow_profile)
                    end if
                    call npse_eqn%npse_setup()

                    do n = -n_global, n_global
                        do m = 0, m_global
                            call npse_eqn%set_parameter("ome", theta1%re*m+theta2%re*n)
                            call npse_eqn%set_parameter("bat", theta1%im*m+theta2%im*n)
                            call npse_eqn%set_parameter("solve_alpha", SOLVE_ALPHA_NONE)

                            !   call VecCopy(phi(i - 1, m, n), phi(i, m, n), ierr)
                            call npse_eqn%set_force(nonlinear_term(:, :, m, n))

                            call VecCopy(phi(i, m, n), phi_old(m, n), ierr)

                            call npse_eqn%solve_npse(alf(i, m, n), phi(i, m, n), phi(i - 1, m, n), phi(i - 2, m, n), lambda)
                            VecDistance(m, n) = VecsNorm2(phi(i, m, n), phi_old(m, n))
                            !! use phi_new to store V_new-V_old
                            call VecAXPY(phi(i, m, n), (-1.0d0, 0.0d0), phi_old(m, n), ierr)
                            call VecAYPX(phi(i, m, n), cmplx(omega_phi, 0.0d0, REAL64), phi_old(m, n), ierr)

                            iterator_times(m, n) = iterator_times(m, n) + npse_eqn%iterator_times

                            lnA0(i, m, n) = lnA1(m, n) + 0.5d0*(alf(i, m, n) + alf(i - 1, m, n))/grid_profile%xi_x(yb + 1)
                            ! norm_rhs_mn(m, n) = npse_eqn%normrhs
                            call VecNorm(phi(i, m, n), NORM_2, norm_rhs_mn(m, n), ierr)
                            !   print*, 'the process is', omp_get_thread_num()

                        end do
                    end do

                    Norm_rhs = sqrt((sum(norm_rhs_mn**2)))
                    print "(I5, 1X, A, E13.5, 1X, A, E13.5)", itimes_nonlinear, 'The 2-Norm is', &
                        sqrt(sum(VecDistance**2)), &
                        'The relative error is', sqrt(sum(VecDistance**2))/Norm_rhs

                    call nonlinear_solver%set_dis(alf(i, :, :), lnA0(i, :, :), phi(i, :, :), phi(i - 1, :, :), phi(i - 2, :, :))
                    call nonlinear_solver%setup()
                    call nonlinear_solver%get_force(nonlinear_term)
                    nonlinear_term = omega_force*nonlinear_term + (1.0d0 - omega_force)*nonlinear_term_old
                    ! nonlinear_term_old = nonlinear_term

                    if (sqrt(sum(VecDistance**2)) >= 1.0d20) then
                        print *, '[WARNING:] The NPSE is devergence at ix=', i
                        print *, 'You can adjust omega_force_factor or/and omega_phi_factor for a more stable impliment.'
                        goto 1000
                    end if
                    ! if(sum(iterator_times)==0) VecDistance = VecDistance*10d2
                    itimes_nonlinear = itimes_nonlinear + 1
                    if (itimes_nonlinear > maxit_nonlinear) then
                        print *, '[WARNING] the iteration for nonlinear larger than max times.'
                        goto 1000
                    endif

                end do


                ! print *, 'ix=', i, 'is done.'
                print *, 'The whole iteration times is', sum(iterator_times)
                print "(A,I2,A,I2,A, F13.7,A,F14.7, A, F10.5)", 'alpha(', inlet_dis(1)%m, ', ', inlet_dis(1)%n, ')=', & !,,A, F10.5
                    alf(i, inlet_dis(1)%m, inlet_dis(1)%n)%re, ' ', alf(i, inlet_dis(1)%m, inlet_dis(1)%n)%im, &
                    "im,   bat=", theta1%im*inlet_dis(1)%m+theta2%im*inlet_dis(1)%n

                itimes_alpha = itimes_alpha + 1

            end do

            if (itimes_alpha == maxit_alpha) print *, '[WARNING] the iteration for alpha larger than max times.'
            print *, 'ix=', i, 'solve npse is done.'
            lnA1(:, :) = lnA0(i, :, :)

            if (mod(i, 10) == 0) call write_npse_result()
        end do

1000    continue
        call PetscPrintf(PETSC_COMM_WORLD, 'Solve NPSE done!\n', ierr)

        call write_npse_result()

        deallocate (iterator_times)

    end subroutine solve_npse_native

    logical function isConvergence()

        implicit none
        real(REAL64) :: maxnonlinear, maxnonlinear_old
        logical :: first_step = .true.

        maxnonlinear = maxval(abs(nonlinear_term))
        maxnonlinear_old = maxval(abs(nonlinear_term_old))

        if (abs(maxnonlinear - maxnonlinear_old)/maxnonlinear <= 1.0d-3 .and. (.not. first_step)) then
            isConvergence = .True.
            first_step = .true.
        else
            isConvergence = .False.
            first_step = .false.
        end if
        print *, 'nonlinear variation:', abs(maxnonlinear - maxnonlinear_old)/maxnonlinear

    end function

    subroutine update_alpha(i, grid)

        implicit none
        integer, intent(in) :: i
        type(type_lst_grid_profile), intent(in) :: grid
        complex(REAL64):: alf_tmp
        integer :: m, n
!! 对atnpse需要进行修改 修改更新波数的方法和方式以及判断标准
        if (abs(inlet_dis(1)%m*theta1%re+inlet_dis(1)%n*theta2%re) >= epsilon(1.0d0)) then  !! 非横流
            alf_tmp = alf(i, inlet_dis(1)%m, inlet_dis(1)%n)/(inlet_dis(1)%m*theta1%re+inlet_dis(1)%n*theta2%re)
            do n = -n_global, n_global
                do m = 0, m_global
                    if (solve_alpha(m, n) == SOLVE_ALPHA_REAL) then
                        alf(i, m, n) = alf(i, m, n) !+ ZI*aimag(alf_tmp)*m
                    elseif (solve_alpha(m, n) == SOLVE_ALPHA_NONE) then
                        alf(i, m, n) = alf_tmp*(m*theta1%re+n*theta2%re)
                    end if
                    if (i > istart) lnA0(i, m, n) = lnA1(m, n) + 0.5d0*(alf(i, m, n) + alf(i - 1, m, n))/grid%xi_x(yb + 1)
                end do
            end do
        else  !! 横流
            alf_tmp = alf(i, inlet_dis(1)%m, inlet_dis(1)%n)/(inlet_dis(1)%m*theta1%im+inlet_dis(1)%n*theta2%im)
            do n = -n_global, n_global
                do m = 0, m_global
                    if (solve_alpha(m, n) == SOLVE_ALPHA_REAL) then
                        alf(i, m, n) = alf(i, m, n) !+ ZI*aimag(alf_tmp)*m
                    elseif (solve_alpha(m, n) == SOLVE_ALPHA_NONE) then
                        alf(i, m, n) = alf_tmp*(m*theta1%im+n*theta2%im)
                    end if
                    if (i > istart) lnA0(i, m, n) = lnA1(m, n) + 0.5d0*(alf(i, m, n) + alf(i - 1, m, n))/grid%xi_x(yb + 1)
                end do
            end do
        end if

    end subroutine update_alpha

    real(REAL64) function VecsNorm2(Vec1, Vec2)

        implicit none
        Vec, intent(in) :: Vec1
        Vec, intent(in) :: Vec2
        Vec :: v
        PetscErrorCode :: ierr
        Real(REAL64) :: Norm1, Norm2
        Complex(REAL64), dimension(10) :: v1, v2

        call VecDuplicate(Vec1, v, ierr)
        ! call VecGetValues(Vec1, 10, [1,2,3,4,5,6,7,8,9,10], v1, ierr)
        ! call VecGetValues(Vec2, 1, [1,2,3,4,5,6,7,8,9,10], v2, ierr)
        call VecWAXPY(v, (-1.0d0, 0.0d0), Vec2, Vec1, ierr)
        call VecNorm(Vec2, NORM_2, Norm2, ierr)
        call VecNorm(v, NORM_2, VecsNorm2, ierr)
        ! print*, vecsNorm2, vecNorm1
        ! VecsNorm2=abs(VecsNorm2/VecNorm1)
        call VecDestroy(v, ierr)

    end function VecsNorm2

    subroutine finalize_npse()
        implicit none

        integer :: i, m, n

        call npse_eqn%finalize_force()
        call npse_eqn%finalize()
        call nonlinear_solver%finalize()
        deallocate (nonlinear_term)
        deallocate (nonlinear_term_old)
        deallocate (solve_alpha)
        deallocate (lnA0)
        deallocate (lnA1)
        deallocate (alf)
        do n = -n_global, n_global
            do m = 0, m_global
                do i = istart - 1, iend
                    call DMRestoreGlobalVector(DM_local, phi(i, m, n), ierr)
                end do
                call DMRestoreGlobalVector(DM_local, phi_old(m, n), ierr)
            end do
        end do
        do i = 1, 3
            call DMRestoreGlobalVector(DM_local, delta_baseflow(i), ierr)
        end do
        deallocate (phi)
        deallocate (phi_old)
        deallocate (VecDistance)
        call DMDestroy(DM_local, ierr)

    end subroutine finalize_npse

    subroutine write_npse_result()

        use global_parameters, only: prefix, gridfile, flowfile, title
        use stdlib_string_type
        implicit none

        PetscViewer :: Viewer
        integer :: iounit, ios
        type(string_type) :: name, strings
        type(string_type) :: astring(6)
        integer :: i, m, n
        Vec :: phi_domain
        PetscScalar, pointer :: tmp_domain(:, :, :), tmp_local(:, :)
        integer(8) :: hdf_file_id

        REAL(REAL64), parameter :: epsilon = 1.0d-8
    !! Write Number
        ! name=solname//'_i='//chars(istart)//'_i='//chars(iend)//'_alpha.plt'
        astring(1) = "NPSE"
        astring(2) = string_type(istart)!; astring(2)=astring(2)%slice(2, max(2, ceiling(log10(istart + 0.1d0) + 1)))
        astring(2) = 'istart='//astring(2)
        astring(3) = string_type(iend)!; astring(3)=astring(3)%slice(2, max(2, ceiling(log10(iend + 0.1d0) + 1)))
        astring(3) = 'iend='//astring(3)
        call DMGetGlobalVector(meshDA, phi_domain, ierr)

        astring(6) = 'phi.h5'
        name = astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(6) !name%join(astring(1:3), sep='_')//'_'//astring(6)
        call PetscViewerHDF5Open(PETSC_COMM_WORLD, char(prefix//name), FILE_MODE_WRITE, Viewer, ierr)
        strings = "Writing the shape functions phi to "//name//'.\n'
        call PetscPrintf(PETSC_COMM_WORLD, char(strings), ierr)

        call PetscViewerHDF5WriteAttribute(viewer, "/", "Problem", PETSC_STRING, "NPSE", ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "Grid", PETSC_STRING, char(gridfile), ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "Flow", PETSC_STRING, char(flowfile), ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "Title", PETSC_STRING, char(title), ierr)

        call PetscViewerHDF5WriteAttribute(viewer, "/", "Ma", PETSC_DOUBLE, Ma, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "Re", PETSC_DOUBLE, Re, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "Te", PETSC_DOUBLE, Te, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "omega1", PETSC_DOUBLE, theta1%re, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "beta1", PETSC_DOUBLE, theta1%im, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "omega2", PETSC_DOUBLE, theta2%re, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "beta2", PETSC_DOUBLE, theta2%im, ierr)

        call PetscViewerHDF5WriteAttribute(viewer, "/", "jn", PETSC_INT, jn, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "istart", PETSC_INT, istart, ierr)
        call PetscViewerHDF5WriteAttribute(viewer, "/", "iend", PETSC_INT, iend, ierr)

        do n = -n_global, n_global
            do m = 0, m_global

                astring(4) = string_type(m) !astring(4)=astring(4)%slice(2, max(2, ceiling(log10(m + 0.1d0) + 1)))
                astring(4) = 'm='//astring(4)
                astring(5) = string_type(n)!; astring(5)=astring(5)%slice(2, max(2, ceiling(log10(n + 0.1d0) + 1)))
                astring(5) = 'n='//astring(5)
                !  astring(6) = 'alpha.plt'
                name = astring(4)//"_"//astring(5)
                !  strings = "Writing the wave numbers alpha to "//name//'\n'
                !  call PetscPrintf(PETSC_COMM_WORLD, strings//'', ierr)

               !! Write phi
                call DMDAVecGetArrayF90(meshDA, phi_domain, tmp_domain, ierr)
                do i = istart, iend
                    call DMDAVecGetArrayReadF90(DM_local, phi(i, m, n), tmp_local, ierr)
                    tmp_domain(:, i - 1, :) = tmp_local(:, :)*exp(ZI*lnA0(i, m, n))
                    call DMDAVecRestoreArrayReadF90(DM_local, phi(i, m, n), tmp_local, ierr)
                end do
                call DMDAVecRestoreArrayF90(meshDA, phi_domain, tmp_domain, ierr)

                call PetscObjectSetName(phi_domain, char(name), ierr)
                call VecView(phi_domain, Viewer, ierr)

                !! Write Alpha
                !  call PetscViewerHDF5GetFileId(Viewer, hdf_file_id, ierr)
                !  call PetscViewerHDF5WriteAttribute(viewer, name//'', "alpha", PETSC_COMPLEX, alf(istart:iend), ierr)
                !     open (newunit=iounit, file=prefix//name, iostat=ios, action="write")
                !     if (ios /= 0) stop "Error opening file "

                !     write (iounit, '(A)') '#variables=i, x, Ma, Re, Te, omega_r, omega_i, &
                ! &                alpha_r, alpha_i, beta_r, beta_i, c_r, c_i'
                !     do i = istart, iend
                !         write (iounit, ('(I0, 12(2X, F20.8))')) &
                !             i, xx(i, yb + 1), Ma, Re, Te, ome0*m, &
                !             alf(i, m, n), bat0*n, ome0*m/(alf(i, m, n) + epsilon)
                !     end do
                !     close (iounit)
                !     call PetscPrintf(PETSC_COMM_WORLD, 'Write '//name//' done.\n', ierr)

            end do
        end do
        call PetscViewerDestroy(Viewer, ierr)

        call PetscPrintf(PETSC_COMM_WORLD, 'Write NPSE result file '// &
                         char(astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(6))// &
                         ' done.\n', ierr)
        call VecDestroy(phi_domain, ierr)

    end subroutine write_npse_result

end module mod_npse
