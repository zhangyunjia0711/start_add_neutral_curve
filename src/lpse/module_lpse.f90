!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  module: mod_solve_lpse
!> @breif module for solving lpse.
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

module mod_lpse

    ! use global_parameters
    use iso_fortran_env, only: REAL64
    use mod_lst_baseflow_profile
    use mod_lst_grid_profile
    use global_parameters, only: sigma_control, sigma_p
    use global_parameters, only: rd_factor

    use global_parameters, only: alpha, beta, omega
    use global_parameters, only: solve_alpha
    use global_parameters, only: in, jn
    use global_parameters, only: xx
    use global_parameters, only: Re, Ma, Te, Pr
    use global_parameters, only: lpse_IsAdjoint
    use global_parameters, only: yb, ye!, xb, xe

    use global_parameters, only: istart, iend
    use global_parameters, only: isothermal
    use global_parameters, only: bc_y0, bc_yn
    use global_parameters, only: Coordinate_System

    ! use global_parameters, only: phi
    ! use global_parameters, only: wave_param
    ! use global_parameters, only: lnA

    use global_parameters, only: meshDA
    use global_parameters, only: rank

    use stdlib_string_type
    use petsc

    ! use mod_lpse_grid_profile, only: type_lpse_grid_domain=> type_lpse_grid_profile
    use mod_lpse_eqn
    use mod_alpse_eqn
    use mod_lpse_solver_abstract

    implicit none

#include "../include/start.f90"

    ! complex(REAL64) :: ome, alf, bat
    DM :: DM_local
    class(type_lpse_abstract), allocatable :: lpse_eqn
    PetscErrorCode :: ierr
    type(string_type) :: solname

    ! type :: type_waveparam
    complex(REAL64), allocatable :: alf(:)
    complex(REAL64) :: ome
    complex(REAL64) :: bat
    ! end type

    ! type(type_waveparam), allocatable :: wave_param(:)
    complex(REAL64), allocatable :: lnA(:)
    Vec, dimension(:), allocatable :: phi

    private
    public :: lpse

    ! DM :: DM_LPSE

contains

    subroutine lpse()

        implicit none

        if (lpse_IsAdjoint) then
            solname = 'ALPSE'
            allocate (type_alpse_eqn :: lpse_eqn)
            call PetscPrintf(PETSC_COMM_WORLD, 'Adjoint PSE is selected!\n', ierr)
        else
            solname = 'LPSE'
            allocate (type_lpse_eqn :: lpse_eqn)
            !DM_local=>lpse_eqn%dm
            call PetscPrintf(PETSC_COMM_WORLD, 'Linear PSE is selected!\n', ierr)
        end if
        call lpse_eqn%initial(jn, Re, Ma, Te, yb, ye)
        lpse_eqn%isothermal = isothermal
        call lpse_eqn%set_parameter("ome", omega(1))
        call lpse_eqn%set_parameter("bat", beta(1))
        call lpse_eqn%set_parameter("sigma_p", sigma_p)
        call lpse_eqn%set_parameter("sigma_control", sigma_control)
        call lpse_eqn%set_parameter("rd_factor", rd_factor)
        call lpse_eqn%set_parameter("Pr", Pr)
        call lpse_eqn%set_parameter("Coordinate_System", Coordinate_System)
        call lpse_eqn%set_parameter("bc_y0", bc_y0)
        call lpse_eqn%set_parameter("bc_yn", bc_yn)

        call initialize_lpse()

        if (lpse_IsAdjoint) then
            call lpse_eqn%set_parameter("solve_alpha", solve_alpha)
            call solve_lpse_adjoint()
        else
            call lpse_eqn%set_parameter("solve_alpha", solve_alpha)
            call solve_lpse_native()
        end if

        call finalize_lpse()

    end subroutine lpse

    subroutine initialize_lpse()

        implicit none
        integer :: i

        call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 2, PETSC_NULL_INTEGER, DM_local, ierr)
        call DMSetUp(DM_local, ierr)

        allocate (phi(istart-1:iend+1))
        allocate (alf(istart:iend))
        do i = istart-1, iend+1
            call DMGetGlobalVector(DM_local, phi(i), ierr)
            call VecZeroEntries(phi(i), ierr)
        end do
        ome = omega(1)
        bat = beta(1)
        alf(istart) = alpha(1) + ZI*alpha(2)
        alf(iend) = alpha(1) + ZI*alpha(2)
        allocate (lnA(istart:iend)); lnA = 0.0d0

        call set_lpse_ic()
        ! do i=xb+2, xe+1
        !   call VecCopy(phi(xb+1), phi(i), ierr)
        ! enddo

    end subroutine initialize_lpse

    subroutine set_lpse_ic()

        use mod_lst_solver_abstract
        use mod_lst_eqn
        use mod_alst_eqn
        implicit none
        integer :: i
        PetscViewer :: Viewer
        integer :: iounit, ios
        type(string_type) :: name, strings
        type(string_type) :: astring(4)
        integer :: idump
        real(REAL64) :: rdump
        class(type_lst_abstract), allocatable :: lst
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        integer :: iloc

    !! read wave numbers form files when alpse is computed without wave number
        if (lpse_IsAdjoint .and. solve_alpha == SOLVE_ALPHA_NONE) then

            astring(1) = 'LPSE'
            astring(2) = string_type(istart); astring(2) = 'istart='//astring(2)
            astring(3) = string_type(iend); astring(3) = 'iend='//astring(3)
            astring(4) = 'alpha.plt'
            name = astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(4)
            strings = "Reading the wave numbers alpha from "//name//'.\n'
            call PetscPrintf(PETSC_COMM_WORLD, char(strings), ierr)

            open (newunit=iounit, file='./out/'//char(name), iostat=ios, action="read")
            if (ios /= 0) stop "Error opening wavenumber file which is required in alpse."
            read (iounit, *)
            do i = istart, iend

                read (iounit, *) &
                    idump, rdump, rdump, rdump, rdump, omega(1), omega(2), &
                    alpha(1), alpha(2), beta(1), beta(2), rdump, rdump

                alpha(2) = alpha(2) !! conjg
                alf(i) = alpha(1) + alpha(2)*ZI

                ! lnA(i)=0.0d0
            end do
            close (iounit)
        end if

        ! !!! Form Files
        ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "in/"//solname//"_disturbance+1.petsc", FILE_MODE_READ, Viewer, ierr)
        ! if( lpse_IsAdjoint) then
        !   call VecLoad(phi(iend), Viewer, ierr)
        ! else
        !   call VecLoad(phi(istart), Viewer, ierr)
        ! endif
        ! call PetscViewerDestroy(Viewer, ierr)

        !!! compute from lst directly
        if (.not. lpse_IsAdjoint) then
            allocate (type_lst_eqn :: lst)
            iloc = istart
        else
            allocate (type_alst_eqn :: lst)
            iloc = iend
        end if
        call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
        call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
        call lst%set_parameter("Pr", Pr)
        call lst%set_parameter("Coordinate_System", Coordinate_System)
        call lst%set_parameter("bc_y0", bc_y0)
        call lst%set_parameter("bc_yn", bc_yn)

        lst%isothermal = isothermal

        grid_profile = assemble_lst_grid_profile(iloc)
        baseflow_profile = assemble_lst_baseflow_profile(iloc)

        call lst%set_grid(grid_profile)
        call lst%set_baseflow(baseflow_profile)
        call lst%setup()

        call lst%solve(ome, alf(iloc), bat)
        call lst%get_eigf(phi(iloc))
        call lst%finalize()

    end subroutine set_lpse_ic

    subroutine solve_lpse_native()

        implicit none
        integer :: i
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        PetscViewer :: Viewer

        call lpse_eqn%set_ic(phi(istart))
        do i = istart + 1, iend
            grid_profile = assemble_lst_grid_profile(i)
            baseflow_profile = assemble_lst_baseflow_profile(i)

            call lpse_eqn%set_grid(grid_profile)
            call lpse_eqn%set_baseflow(baseflow_profile)
            call lpse_eqn%setup()
            ! call solve_lpse_eqn(grid_profile, baseflow_profile, ome, alf, bat, i, phi(i), phi(i-1), phi(max(1,i-2)))
            alf(i) = alf(i - 1)
            call VecCopy(phi(i - 1), phi(i), ierr)
            call lpse_eqn%solve(alf(i), phi(i), phi(i-1), phi(i-2))

            lnA(i) = lnA(i - 1) + 0.5d0*(alf(i) + alf(i - 1))/grid_profile%xi_x(yb + 1)
            if (rank == 0) then
                print *, 'ix=', i, 'is done.', 'x=', xx(i, 1)
                print *, 'alpha=', alf(i), "bat=", bat
            end if
        end do
        call PetscPrintf(PETSC_COMM_WORLD, 'Solve LPSE done!\n', ierr)

        call write_lpse_result()

        ! call finalize_lpse_eqn()
        call lpse_eqn%finalize()

    end subroutine solve_lpse_native

    subroutine solve_lpse_adjoint()

        implicit none
        integer :: i
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        integer :: iounit, ios
        type(string_type) :: name, strings
        type(string_type) :: astring(4)
        integer :: idump
        real(REAL64) :: rdump

        print *, 'solving adjoint LPSE!'

        grid_profile = assemble_lst_grid_profile(iend)
        baseflow_profile = assemble_lst_baseflow_profile(iend)

        call lpse_eqn%set_ic(phi(iend))

        do i = iend - 1, istart, -1
            grid_profile = assemble_lst_grid_profile(i)

            baseflow_profile = assemble_lst_baseflow_profile(i)
            ! ome=wave_param(i)%ome
            ! alf=wave_param(i)%alf
            ! bat=wave_param(i)%bat

            call lpse_eqn%set_grid(grid_profile)
            call lpse_eqn%set_baseflow(baseflow_profile)
            call lpse_eqn%setup()

            if (solve_alpha /= SOLVE_ALPHA_NONE) then
                print *, 'start solve alpse at i=', i
                alf(i) = alf(i + 1)
            end if
            call VecCopy(phi(i + 1), phi(i), ierr)
            call lpse_eqn%solve(alf(i), phi(i), phi(i+1), phi(i+2))

            if (rank == 0) then
                print *, 'ix=', i, 'is done.'
                print *, 'alpha=', alf(i), "bat=", bat
            end if
        end do
        call PetscPrintf(PETSC_COMM_WORLD, 'Solve APSE done!\n', ierr)

        call write_lpse_result()

        call lpse_eqn%finalize()

    end subroutine solve_lpse_adjoint

    subroutine finalize_lpse()
        implicit none

        integer :: i
        deallocate (lnA)
        deallocate (alf)
        do i = istart-1, iend+1
            call DMRestoreGlobalVector(DM_local, phi(i), ierr)
        end do
        deallocate (phi)
        call DMDestroy(DM_local, ierr)

    end subroutine finalize_lpse

    subroutine write_lpse_result()

        use global_parameters, only: GAMMA, Ma, qq
        implicit none
        PetscViewer :: Viewer
        integer :: iounit, ios
        type(string_type) :: name, strings
        type(string_type) :: astring(4)
        integer :: i
        Vec :: phi_domain
        PetscScalar, pointer :: tmp_domain(:, :, :), tmp_local(:, :)
        PetscScalar, allocatable :: A_domain(:, :, :)
        real(REAL64) :: pe

        REAL(REAL64), parameter :: epsilon = 1.0d-8
    !! Write Number
        ! name=solname//'_i='//chars(istart)//'_i='//chars(iend)//'_alpha.plt'
        astring(1) = solname
        astring(2) = string_type(istart); astring(2) = 'istart='//astring(2)
        astring(3) = string_type(iend); astring(3) = 'iend='//astring(3)
        astring(4) = 'alpha.plt'
        name = astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(4)
        strings = "Writing the wave numbers alpha to "//name//'\n'
        call PetscPrintf(PETSC_COMM_WORLD, char(strings), ierr)

        open (newunit=iounit, file='./out/'//char(name), iostat=ios, action="write")
        if (ios /= 0) stop "Error opening file "

        write (iounit, '(A)') '#variables=i, x, Ma, Re, Te, omega_r, omega_i, &
    &                alpha_r, alpha_i, beta_r, beta_i, c_r, c_i'
        do i = istart, iend
            write (iounit, ('(I0, 12(2X, F20.8))')) &
                i, xx(i, yb + 1), Ma, Re, Te, ome, &
                alf(i), bat, ome/(alf(i) + epsilon)
        end do
        close (iounit)
        call PetscPrintf(PETSC_COMM_WORLD, 'Write wavenum done.\n', ierr)

    !! Write phi
        call DMGetGlobalVector(meshDA, phi_domain, ierr)
        call DMDAVecGetArrayF90(meshDA, phi_domain, tmp_domain, ierr)
        allocate(A_domain, source = tmp_domain)
        do i = istart, iend
            call DMDAVecGetArrayReadF90(DM_local, phi(i), tmp_local, ierr)
            tmp_domain(:, i - 1, :) = tmp_local(:, :) ! *exp(ZI*lnA(i))
            A_domain(:, i - 1, :) = tmp_local(:, :) *exp(ZI*lnA(i))
            call DMDAVecRestoreArrayReadF90(DM_local, phi(i), tmp_local, ierr)
        end do
        call DMDAVecRestoreArrayF90(meshDA, phi_domain, tmp_domain, ierr)

        astring(4) = 'phi.h5'
        name = astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(4)
        strings = "Writing the shape phi phi to "//name//'\n'
        call PetscPrintf(PETSC_COMM_WORLD, char(strings), ierr)
        call PetscViewerHDF5Open(PETSC_COMM_WORLD, "./out/"//char(name), FILE_MODE_WRITE, Viewer, ierr)
        call PetscObjectSetName(phi_domain, "phi", ierr)
        call VecView(phi_domain, Viewer, ierr)
        call VecDestroy(phi_domain, ierr)
        call PetscViewerDestroy(Viewer, ierr)

        !!! Write Amp
        astring(4) = 'Amp.plt'
        name = astring(1)//"_"//astring(2)//"_"//astring(3)//"_"//astring(4)
        strings = "Writing the amplitude to "//name//'\n'
        call PetscPrintf(PETSC_COMM_WORLD, char(strings), ierr)

        open (newunit=iounit, file='./out/'//char(name), iostat=ios, action="write")
        if (ios /= 0) stop "Error opening file "

        write (iounit, '(A)') '#variables=x, A_p, A_u, A_T'
        pe = GAMMA*Ma**2
        do i = istart, iend
            write (iounit, ('(4(2X, F20.8))')) &
                xx(i, yb + 1), &
                (abs((A_domain(0, i-1, 0)*qq(5, i, 1)+A_domain(4, i-1, 0)*qq(1, i, 1))*pe)), &
                maxval(abs(A_domain(1, i-1, :))), maxval(abs(A_domain(4, i-1, :)))
        end do
        close (iounit)
        call PetscPrintf(PETSC_COMM_WORLD, 'Write amplitude done.\n', ierr)


        call PetscPrintf(PETSC_COMM_WORLD, 'Write phi done.\n', ierr)

        ! print*, solve_alpha
        ! if(solve_alpha /= 0)then
        !   call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "./out/"//astring(1)//"_outlet.petsc", FILE_MODE_WRITE, Viewer, ierr)
        !   if(lpse_IsAdjoint)then
        !     call VecView(phi(istart), Viewer, ierr)
        !   else
        !     call VecView(phi(iend), Viewer, ierr)
        !   endif
        !   call PetscViewerDestroy(Viewer, ierr)
        ! endif

    end subroutine write_lpse_result

end module mod_lpse
