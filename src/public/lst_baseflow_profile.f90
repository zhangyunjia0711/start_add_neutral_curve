!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-05 09:32:33
!> @LastEditTime: 2022-12-12 16:18:01
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/lst_baseflow_profile.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lst_baseflow_profile

    use global_parameters
    use mod_wall_local_coordinates, only: wall_local_enabled, wall_local_flow
    implicit none

    public :: type_lst_bf_profile
    public :: assemble_lst_baseflow_profile
    public :: get_baseflow_profile
    public :: operator(+)
    private

    interface operator(+)
      module procedure :: add_baseflow_profile
    end interface

    type :: type_lst_bf_profile

        logical :: has_assemble = .FALSE.
        integer :: jb, je
        integer :: ix
        real(REAL64), allocatable, dimension(:, :) :: q0
        real(REAL64), allocatable, dimension(:, :) :: dq0y
        real(REAL64), allocatable, dimension(:, :) :: dq0x
        real(REAL64), allocatable, dimension(:, :) :: dq0xx
        real(REAL64), allocatable, dimension(:, :) :: dq0yy
        real(REAL64), allocatable, dimension(:, :) :: dq0xy

    contains

        procedure :: print => print_lst_baseflow_profile
        procedure :: finalize => finalize_lst_baseflow_profile

    end type type_lst_bf_profile

contains

    function assemble_lst_baseflow_profile(ix) result(profile)

        implicit none
        integer, intent(in) :: ix
        type(type_lst_bf_profile) :: profile

        integer :: jb, je

        jb = yb + 1; je = ye + 1
        if (.not. profile%has_assemble) then
            allocate (profile%q0(5, jb:je))
            allocate (profile%dq0x(5, jb:je))
            allocate (profile%dq0y(5, jb:je))
            allocate (profile%dq0xx(5, jb:je))
            allocate (profile%dq0yy(5, jb:je))
            allocate (profile%dq0xy(5, jb:je))
        end if

        profile%ix = ix

        profile%has_assemble = .TRUE.
        profile%jb = yb + 1; profile%je = ye + 1

        profile%q0 = qq(:, ix, jb:je)
        profile%dq0x = dq0dx(:, ix, jb:je)
        profile%dq0y = dq0dy(:, ix, jb:je)
        profile%dq0xx = dq0dxx(:, ix, jb:je)
        profile%dq0yy = dq0dyy(:, ix, jb:je)
        profile%dq0xy = dq0dxy(:, ix, jb:je)

        if (wall_local_enabled()) then
            call wall_local_flow(ix,profile%q0,profile%dq0y,profile%dq0yy)
            ! The direct local operator does not use these derivatives. Do not
            ! expose old global derivatives as if they were local/nonparallel.
            profile%dq0x=0.0_real64
            profile%dq0xx=0.0_real64
            profile%dq0xy=0.0_real64
        endif

    end function assemble_lst_baseflow_profile

    function get_baseflow_profile(dm, vecs, grid, lambda) result(bf_profile)

        use mod_lst_grid_profile, only: type_lst_grid_profile
        implicit none
        DM, intent(in) :: dm
        Vec, intent(in) :: vecs(3)
        type(type_lst_grid_profile), intent(in) :: grid
        real(REAL64), intent(in) :: lambda
        type(type_lst_bf_profile) :: bf_profile
        integer :: jb, je
        integer :: l, j

        PetscScalar, pointer :: phi0_array(:, :), phi1_array(:, :), phi2_array(:, :)

        call DMDAVecGetArrayReadF90(dm, Vecs(1), phi0_array, ierr)
        call DMDAVecGetArrayReadF90(dm, Vecs(2), phi1_array, ierr)
        call DMDAVecGetArrayReadF90(dm, Vecs(3), phi2_array, ierr)

        jb = yb + 1; je = ye + 1

        if (.not. bf_profile%has_assemble) then
            allocate (bf_profile%q0(5, jb:je))
            allocate (bf_profile%dq0x(5, jb:je))
            allocate (bf_profile%dq0y(5, jb:je))
            allocate (bf_profile%dq0xx(5, jb:je))
            allocate (bf_profile%dq0yy(5, jb:je))
            allocate (bf_profile%dq0xy(5, jb:je))
        end if

        block
            use mod_dx
            complex(REAL64) :: dphidi(5, jb:je)
            complex(REAL64) :: dphidj(5, jb:je)
            complex(REAL64) :: dphidjj(5, jb:je)
            complex(REAL64) :: dphidij(5, jb:je)

            bf_profile%q0 = phi0_array
            dphidi(:, :) = (phi0_array(:, :)*(1.0d0 + lambda) + &
                          & phi1_array(:, :)*(-1.0d0 - 2.0d0*lambda) + &
                          & phi2_array(:, :)*lambda)

            do l = 1, 5
                call dx_fdm4nd_center(phi0_array(l - 1, :), jb, je, dphidj(l, :), jb, je)
                call dxx_fdm4nd_center(phi0_array(l - 1, :), jb, je, dphidjj(l, :), jb, je)
                call dx_fdm4nd_center(dphidi(l, :), jb, je, dphidij(l, :), jb, je)
            end do

            do j = jb, je
                bf_profile%dq0x(:, j) = dphidi(:, j)*grid%xi_x(j) + dphidj(:, j)*grid%eta_x(j)
                bf_profile%dq0y(:, j) = dphidi(:, j)*grid%xi_y(j) + dphidj(:, j)*grid%eta_y(j)

                bf_profile%dq0xy(:, j) = dphidi(:, j)*grid%xi_xy(j) + dphidj(:, j)*grid%eta_xy(j) + &
                                       & dphidjj(:, j)*(grid%eta_y(j)*grid%eta_x(j)) + &
                                       & dphidij(:, j)*(grid%xi_x(j)*grid%eta_y(j) + grid%xi_y(j)*grid%eta_x(j))
                bf_profile%dq0yy(:, j) = dphidi(:, j)*grid%xi_yy(j) + dphidj(:, j)*grid%eta_yy(j) + &
                                       & dphidjj(:, j)*(grid%eta_y(j)*grid%eta_y(j)) + &
                                       & dphidij(:, j)*(2.0d0*grid%xi_y(j)*grid%eta_y(j))
                bf_profile%dq0xx(:, j) = dphidi(:, j)*grid%xi_xx(j) + dphidj(:, j)*grid%eta_xx(j) + &
                                    & dphidjj(:, j)*(grid%eta_x(j)*grid%eta_x(j)) + &
                                    & dphidij(:, j)*(2.0d0*grid%xi_x(j)*grid%eta_x(j))
            end do

        end block
        call DMDAVecRestoreArrayReadF90(dm, Vecs(1), phi0_array, ierr)
        call DMDAVecRestoreArrayReadF90(dm, Vecs(2), phi1_array, ierr)
        call DMDAVecRestoreArrayReadF90(dm, Vecs(3), phi2_array, ierr)

    end function get_baseflow_profile

    function add_baseflow_profile(profile1, profile2) result(profile)

        implicit none
        type(type_lst_bf_profile), intent(in) :: profile1, profile2
        type(type_lst_bf_profile) :: profile

        if (.not. profile%has_assemble) then
            allocate (profile%q0, source=profile1%q0)
            allocate (profile%dq0x, source=profile1%q0)
            allocate (profile%dq0y, source=profile1%q0)
            allocate (profile%dq0xx, source=profile1%q0)
            allocate (profile%dq0yy, source=profile1%q0)
            allocate (profile%dq0xy, source=profile1%q0)
        end if

        profile%q0   =profile1%q0   +profile2%q0   
        profile%dq0x =profile1%dq0x +profile2%dq0x 
        profile%dq0y =profile1%dq0y +profile2%dq0y 
        profile%dq0xx=profile1%dq0xx+profile2%dq0xx
        profile%dq0yy=profile1%dq0yy+profile2%dq0yy
        profile%dq0xy=profile1%dq0xy+profile2%dq0xy

    end function add_baseflow_profile

    subroutine print_lst_baseflow_profile(profile, filename)

        implicit none
        class(type_lst_bf_profile), intent(in) :: profile
        character(*), intent(in) :: filename
        character(3) :: fn
        integer :: j, l

        write (fn, '(I3.3)') profile%ix
        open (99, file=filename//fn//'.plt')
        do j = profile%jb, profile%je
            write (99, '(2X, 25(2X, E16.6))') yy(1, j), (profile%q0(l, j), l=1, 5), &
                (profile%dq0y(l, j), l=1, 5), &
                (profile%dq0yy(l, j), l=1, 5), &
                (profile%dq0x(l, j), l=1, 5), &
                (profile%dq0xy(l, j), l=1, 5), &
                (profile%dq0xx(l, j), l=1, 5)

        end do
        close (99)

    end subroutine print_lst_baseflow_profile

    subroutine finalize_lst_baseflow_profile(profile)

        implicit none
        class(type_lst_bf_profile) :: profile

        profile%has_assemble = .FALSE.
        profile%jb = 0; profile%je = 0; profile%ix = -1
        deallocate (profile%q0)
        deallocate (profile%dq0x)
        deallocate (profile%dq0y)
        deallocate (profile%dq0xx)
        deallocate (profile%dq0yy)
        deallocate (profile%dq0xy)

    end subroutine finalize_lst_baseflow_profile

end module mod_lst_baseflow_profile
