!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-09 10:45:42
!> @LastEditTime: 2022-12-12 16:19:20
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/solver_abstract.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>
module mod_solver_abstract

    use petsc
    use iso_fortran_env, only: REAL64, REAL32
    use mod_lst_grid_profile, only: type_lst_grid_profile
    use mod_lst_baseflow_profile, only: type_lst_bf_profile
    
    implicit none
    private
    public :: type_solver_abstract

#include "../include/start.f90"

    type, abstract :: type_solver_abstract
        integer :: bc_y0 = BC_WALL, bc_yn = BC_FARFIELD
        integer :: Coordinate_System = COMOVING_CS

        real(REAL64) :: GAMMA = 1.4d0, Pr = 0.72d0
        real(REAL64) :: C1 = 110.4D0
        real(REAL64) :: C2 = 110.4D0
        real(REAL64) :: Re, Ma, Te

        logical :: has_initialized = .FALSE.
        logical :: isothermal = .TRUE.

        type(type_lst_bf_profile) :: baseflow
        type(type_lst_grid_profile) :: grid
    contains

        generic :: set_parameter => set_parameters_real8, &
                                  & set_parameters_real4, &
                                  & set_parameters_int

        procedure(initial_interface), deferred :: initial
        procedure(setup_interface), deferred :: setup

        procedure :: set_grid
        procedure :: set_baseflow

        procedure(finalize_interface), deferred :: finalize

        procedure(interface_set_parameters_real8), deferred :: set_parameters_real8
        procedure(interface_set_parameters_real4), deferred :: set_parameters_real4
        procedure(interface_set_parameters_int), deferred :: set_parameters_int

    end type type_solver_abstract

    abstract interface

        subroutine interface_set_parameters_real8(this, keyword, value)

            import type_solver_abstract, real64
            implicit none
            class(type_solver_abstract), intent(inout) :: this
            character(*), intent(in) :: keyword
            real(REAL64), intent(in) :: value

        end subroutine interface_set_parameters_real8

        subroutine interface_set_parameters_real4(this, keyword, value)

            import type_solver_abstract, real32
            implicit none
            class(type_solver_abstract), intent(inout) :: this
            character(*), intent(in) :: keyword
            real(REAL32), intent(in) :: value

        end subroutine interface_set_parameters_real4

        subroutine interface_set_parameters_int(this, keyword, value)

            import type_solver_abstract
            implicit none
            class(type_solver_abstract), intent(inout) :: this
            character(*), intent(in) :: keyword
            integer, intent(in) :: value

        end subroutine interface_set_parameters_int

        subroutine initial_interface(this, jn_, Re, Ma, Te, yb_, ye_)

            import type_solver_abstract, real64
            implicit none
            class(type_solver_abstract), intent(inout) :: this
            integer, intent(in) :: jn_
            real(real64), intent(in) ::  Re, Ma, Te
            integer, intent(in), optional :: yb_, ye_

        end subroutine initial_interface

        subroutine setup_interface(this)
            import type_solver_abstract
            implicit none
            class(type_solver_abstract), intent(inout) :: this

        end subroutine setup_interface

        subroutine finalize_interface(this)

            import type_solver_abstract
            implicit none
            class(type_solver_abstract), intent(inout) :: this

        end subroutine finalize_interface

    end interface
contains

    subroutine set_grid(this, grid)
        implicit none
        class(type_solver_abstract), intent(inout) :: this
        type(type_lst_grid_profile), intent(in) :: grid

        this%grid = grid

    end subroutine set_grid

    subroutine set_baseflow(this, baseflow)
        implicit none
        class(type_solver_abstract), intent(inout) :: this
        type(type_lst_bf_profile), intent(in) :: baseflow

        this%baseflow = baseflow

        ! block !check baseflow

        !         call baseflow%print('baseflow')
        !         ! read(*, *)

        ! endblock

    end subroutine set_baseflow

end module mod_solver_abstract
