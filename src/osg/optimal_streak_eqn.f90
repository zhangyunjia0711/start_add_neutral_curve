!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-29 15:17:05
!> @LastEditTime: 2023-03-15 10:34:48
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/osg/optimal_streak_eqn.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------


#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ module optimal streak equations solver.
module mod_optimal_streak_eqn

    use iso_fortran_env, only: REAL64
    use petsc
    use mod_lst_grid_profile
    use mod_lst_baseflow_profile
    use mod_lpse_eqn
    use mod_alpse_eqn
    use mod_solver_abstract
    use mod_lpse_solver_abstract
    use mod_lst_eqn
    use global_parameters, only: GAMMA, drz

    implicit none

    private
    public :: type_optimal_streak_eqn

#include "../include/start.f90"

    PetscErrorCode :: ierr

    type(type_lpse_eqn)  :: lpse
    type(type_alpse_eqn) :: alpse
    type(type_lst_eqn) :: lst

    type :: type_optimal_streak_eqn !! optimal streak equations solver container.

        Vec :: q0, q1
        Vec :: q0_star, q1_star
        complex(REAL64), allocatable :: alf(:)
        integer :: jn, i_x0, i_x1
        DM :: dm

        REAL(REAL64) :: Ma

        type(type_lst_grid_profile) :: grid_x0, grid_x1
        type(type_lst_bf_profile) :: baseflow_x0, baseflow_x1

        complex(REAL64), allocatable, dimension(:, :, :) :: A_x0
        complex(REAL64), allocatable, dimension(:, :, :) :: A_x1
        Real(REAL64), allocatable, dimension(:, :) :: M_x0, M_x1

    contains

        procedure :: initial => initial_optimal_streak_eqn !! initialize.
        procedure :: setup => setup_optimal_streak_eqn     !! setup.
        procedure :: set_ic => set_ic_optimal_streak       !! set initial condition.
        procedure :: solve => solve_optimal_streak         !! solve the equations.
        procedure :: get_eigf => get_eigf_optimal_streak   !! get eigenfunctions.
        procedure :: finalize => finalize_optimal_streak_eqn  !! finalize.

        procedure, private :: get_lpse_inlet
        procedure, private :: get_alpse_inlet

    end type type_optimal_streak_eqn

contains

    subroutine initial_optimal_streak_eqn(this, jn, Re, Ma, Te, Pr, yb, ye, isothermal)
        !! initial optimal streak equations.
        class(type_optimal_streak_eqn), intent(inout) :: this 
        integer, intent(in) :: jn, yb, ye  !! jn, lbound and ubound.
        real(real64), intent(in) :: Re, Ma, Te, Pr !! flow parameters.
        logical, intent(in) :: isothermal  !! temperature wall bc.

        this%jn = jn
        this%Ma = Ma
        allocate (this%A_x0(5, 5, jn))
        allocate (this%A_x1(5, 5, jn))
        allocate (this%M_x0(5, jn))
        allocate (this%M_x1(5, jn))

        call DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, jn, 5, 1, PETSC_NULL_INTEGER, this%dm, ierr)
        call DMSetUp(this%dm, ierr)
        call DMGetGlobalVector(this%dm, this%q0, ierr); call VecZeroEntries(this%q0, ierr)
        call VecDuplicate(this%q0, this%q1, ierr); call VecZeroEntries(this%q1, ierr)
        call VecDuplicate(this%q0, this%q0_star, ierr); call VecZeroEntries(this%q0_star, ierr)
        call VecDuplicate(this%q0, this%q1_star, ierr); call VecZeroEntries(this%q1_star, ierr)

        call lst%set_parameter("lst_mode", LST_SPATIAL_MODE)
        call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
        call lst%set_parameter("Pr", Pr)
        lst%isothermal = isothermal
        ! allocate(type_lpse_eqn :: lpse)
        call lpse%initial(jn, Re, Ma, Te, yb, ye)
        ! call lpse%set_parameter("solve_alpha", SOLVE_ALPHA_NONE)
        call lpse%set_parameter("solve_alpha", SOLVE_ALPHA_REAL)
        call lpse%set_parameter("Pr", Pr)
        lpse%isothermal = isothermal
        ! allocate(type_alpse_eqn :: alpse)
        call alpse%initial(jn, Re, Ma, Te, yb, ye)
        call alpse%set_parameter("solve_alpha", SOLVE_ALPHA_NONE)
        call alpse%set_parameter("Pr", Pr)
        alpse%isothermal = isothermal

    end subroutine initial_optimal_streak_eqn

    subroutine setup_optimal_streak_eqn(this, i_x0_, grid_x0, baseflow_x0, i_x1_, grid_x1, baseflow_x1, alf_, bat_, ome_)
    !! setup optimal streak equations.
        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        integer, intent(in) :: i_x0_  !! x0 index.
        integer, intent(in) :: i_x1_  !! x1 index.
        type(type_lst_grid_profile), intent(in) :: grid_x0, grid_x1 !! grid profile at x0 and x1.
        type(type_lst_bf_profile), intent(in) :: baseflow_x0, baseflow_x1 !! baseflow profile at x0 and x1.
        complex(real64), intent(in) :: alf_ !! \(\alpha\) in equations.
        real(REAL64), intent(in) :: bat_ !! \(\beta\) in equations.
        real(REAL64), intent(in) :: ome_ !! \(\omega\) in equations.

        this%i_x0 = i_x0_; this%i_x1 = i_x1_
        this%grid_x0 = grid_x0; this%grid_x1 = grid_x1
        this%baseflow_x0 = baseflow_x0; this%baseflow_x1 = baseflow_x1
        if (allocated(this%alf)) deallocate (this%alf)
        allocate (this%alf(i_x0_:i_x1_))
        this%alf(i_x0_) = alf_

        if(abs(drz(i_x0_))<=epsilon(1.0d0))then
            call lpse%set_parameter("bat", bat_)
        else
            call lpse%set_parameter("cycle_m", bat_/drz(i_x0_))
        endif
        if(abs(drz(i_x1_))<=epsilon(1.0d0))then
            call alpse%set_parameter("bat", bat_)
        else
            call alpse%set_parameter("cycle_m", bat_/drz(i_x0_))
        endif

        call lpse%set_parameter("ome", ome_)
        call alpse%set_parameter("ome", ome_)
        call lpse%set_parameter("rd_factor", 1.0d0)
        call alpse%set_parameter("rd_factor", 1.0d0)

        call setM(this, this%baseflow_x0, this%M_x0)
        call setM(this, this%baseflow_x1, this%M_x1)

        call lpse%set_grid(this%grid_x0)
        call lpse%set_baseflow(this%baseflow_x0)
        call lpse%setup()
        call lpse%get_matrix(this%alf(i_x0_), "A_0", this%A_x0)

        call lpse%set_grid(this%grid_x1)
        call lpse%set_baseflow(this%baseflow_x1)
        call lpse%setup()
        call lpse%get_matrix(this%alf(i_x1_), "A_0", this%A_x1)

    end subroutine setup_optimal_streak_eqn

    subroutine setM(this, baseflow, M)

        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        type(type_lst_bf_profile), intent(in) :: baseflow
        Real(REAL64), dimension(5, this%jn), intent(inout) :: M
        integer :: j
        REAL(REAL64) :: pe, cv

        pe = 1.0d0/(GAMMA*this%Ma**2)
        cv = pe/(GAMMA - 1.d0)

        do j = 1, this%jn
            associate (rho => baseflow%q0(1, j), t => baseflow%q0(5, j))
                M(1, j) = t/rho*pe
                M(2, j) = rho
                M(3, j) = rho
                M(4, j) = rho
                M(5, j) = cv*rho/T
            end associate
        end do

    end subroutine setM

    subroutine set_ic_optimal_streak(this, vec)

        use mod_lst_baseflow_profile
        use mod_lst_grid_profile
        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        Vec, optional :: vec
        PetscRandom :: rctx
        PetscInt :: locs(9)
        PetscScalar :: zeros(9)
        REAL(real64) :: norm
        complex(real64) :: alf, ome, bat
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile

        if (present(vec)) then
            call VecCopy(vec, this%q0, ierr)
        else
            call PetscRandomCreate(PETSC_COMM_WORLD, rctx, ierr)
            call PetscRandomSetFromOptions(rctx, ierr)
            call VecSetRandom(this%q0, rctx, ierr)
            call PetscRandomDestroy(rctx, ierr)
            locs=[2, 3, 4, 5, this%jn*5-4, this%jn*5-3, this%jn*5-2, this%jn*5-1, this%jn*5]-1
            zeros=[0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
            call VecSetValues(this%q0, 9, locs, zeros, INSERT_VALUES, ierr)
            call VecNorm(this%q0, NORM_2, norm, ierr)
            call VecScale(this%q0, dcmplx(1.0d0/norm), ierr)
            ! grid_profile = assemble_lst_grid_profile(this%i_x0)
            ! baseflow_profile = assemble_lst_baseflow_profile(this%i_x0)

            ! call lst%set_grid(grid_profile)
            ! call lst%set_baseflow(baseflow_profile)
            ! call lst%setup()
            ! alf = this%alf(this%i_x0) + (0.0d0, 1.0d-7)
            ! ome = lpse%ome
            ! bat = lpse%bat
            ! call lst%solve(ome, alf, bat)
            ! call lst%get_eigf(this%q0)
        end if

    end subroutine set_ic_optimal_streak

    subroutine solve_optimal_streak(this, G_factor)

        use CMD_Progress
        use stdlib_string_type
        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        REAL(real64), intent(out) :: G_factor
        ! Vec, inteng(inout) :: q0
        complex(REAL64) :: alf

        Vec :: q
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        integer :: i
        Real(REAL64) :: Norm
        type(CLS_CMD_Progress) :: progress
        REAL(real64) :: Norm_x0, Norm_x1
        PetscScalar, pointer :: tmp(:, :)
        integer :: sum

        call Progress%Set(N=(this%i_x1 - this%i_x0 + 1)*2, L=40)!// 1700次，显示长度25
        progress%M = "#" !// 已完成部分的字符，不是必须
        progress%O = "." !// 未完成部分的字符，不是必须

        G_factor = epsilon(1.0d0)
        call VecDuplicate(this%q0, q, ierr)
        sum = 0
        do
            sum = sum + 1
            progress%Prefix = string_type(sum)
            progress%Prefix = "optimal streak loop "//char(progress%Prefix)//':'

            call VecCopy(this%q0, q, ierr)
            ! call VecNorm(q, NORM_2, Norm, ierr)
            ! print*, 'initial LPSE at', this%i_x0, ' the 2-Norm is', Norm

            call lpse%set_ic(this%q0)
            ! call VecNorm(q, NORM_2, Norm, ierr)
            ! print*, 'Initial 2-Norm is ', Norm
            do i = this%i_x0 + 1, this%i_x1
                if (mod((i - this%i_x0), 30) == 0) call progress%Put((i - this%i_x0), CMD_PROGRESS_ABSOLUTE)
                grid_profile = assemble_lst_grid_profile(i)
                baseflow_profile = assemble_lst_baseflow_profile(i)
                call lpse%set_grid(grid_profile)
                call lpse%set_baseflow(baseflow_profile)
                call lpse%setup()
                this%alf(i) = this%alf(i - 1)
                ! this%alf(i)=0.0d0
                call lpse%solve(this%alf(i), q)
                ! call VecNorm(q, NORM_2, Norm, ierr)
                ! print*, 'Solving LPSE at', i, ' the 2-Norm is', Norm
            end do
            call VecCopy(q, this%q1, ierr)

            block !! check Norm_E
                call DMDAVecGetArrayReadF90(this%dm, this%q0, tmp, ierr)
                call ArrayNorm(this%jn, this%grid_x0%y, tmp, this%M_x0, Norm_x0)
                call DMDAVecRestoreArrayReadF90(this%dm, this%q0, tmp, ierr)
                call DMDAVecGetArrayReadF90(this%dm, this%q1, tmp, ierr)
                call ArrayNorm(this%jn, this%grid_x1%y, tmp, this%M_x1, Norm_x1)
                call DMDAVecRestoreArrayReadF90(this%dm, this%q1, tmp, ierr)
                ! if (abs(Norm_x1/Norm_x0 - G_factor)/G_factor <= (epsilon(1.0e0)*100.0e0)) exit
                if (abs(Norm_x1/Norm_x0 - G_factor)/G_factor <= 1d-3) then
                    print *, 'Growth_factor:', Norm_x1/Norm_x0, 'G_rtol is:', abs(Norm_x1/Norm_x0 - G_factor)/G_factor
                    exit
                endif
                ! read(*, *)
            end block

            call this%get_alpse_inlet()
            call VecCopy(this%q1_star, q, ierr)
            ! call VecNorm(q, NORM_2, Norm, ierr)
            ! print*, 'initial ALPSE at', this%i_x1, ' the 2-Norm is', Norm

            call alpse%set_ic(this%q1_star)

            do i = this%i_x1 - 1, this%i_x0, -1
                if (mod((this%i_x1 - i), 30) == 0) call progress%Put(this%i_x1*2 - i - this%i_x0, CMD_PROGRESS_ABSOLUTE)
                grid_profile = assemble_lst_grid_profile(i)
                baseflow_profile = assemble_lst_baseflow_profile(i)
                call alpse%set_grid(grid_profile)
                call alpse%set_baseflow(baseflow_profile)
                call alpse%setup()
                this%alf(i) = this%alf(i + 1)
                ! call VecNorm(q, NORM_2, Norm, ierr)
                ! print*, 'Solving ALPSE at', i, 'the 2-Norm first is', Norm
                ! if(i==this%i_x1-1) then
                !    print*, maxval(abs(alpse%phi0_array))
                !    print*, maxval(abs(alpse%phi1_array))
                !    print*, maxval(abs(alpse%phi2_array))
                !    read(*, *)
                ! endif
                call alpse%solve(this%alf(i), q)
                ! call VecNorm(q, NORM_2, Norm, ierr)
                ! print*, 'Solving ALPSE at', i, 'the 2-Norm then is', Norm
            end do
            call VecCopy(q, this%q0_star, ierr)

            print *, 'Growth_factor:', Norm_x1/Norm_x0, 'G_rtol is:', abs(Norm_x1/Norm_x0 - G_factor)/G_factor
            G_factor = Norm_x1/Norm_x0

            call this%get_lpse_inlet()
            ! call VecView(this%q0, PETSC_VIEWER_STDOUT_SELF, ierr)
            ! read(*, *)
        end do

        G_factor = Norm_x1/Norm_x0
        ! VecCopy

    end subroutine

    subroutine get_eigf_optimal_streak(this, eigf_x0, eigf_x1)

        implicit none
        class(type_optimal_streak_eqn), intent(in) :: this
        Vec, intent(inout) :: eigf_x0, eigf_x1

        call VecDuplicate(this%q0, eigf_x0, ierr)
        call VecDuplicate(this%q1, eigf_x1, ierr)
        call VecCopy(this%q0, eigf_x0, ierr)
        call VecCopy(this%q0, eigf_x1, ierr)

    end subroutine get_eigf_optimal_streak

    subroutine get_lpse_inlet(this)

        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        PetscScalar, pointer :: tmp(:, :)
        PetscScalar :: q0_star(5, this%jn)
        PetscScalar :: A_0_h(5, 5)
        real(REAL64) :: E
        integer :: j

        call DMDAVecGetArrayReadF90(this%dm, this%q0_star, tmp, ierr)
        q0_star = tmp
        call DMDAVecRestoreArrayReadF90(this%dm, this%q0_star, tmp, ierr)

        q0_star(1, 1) = 0.0d0
        do j = 1, this%jn
            A_0_h = conjg(transpose(this%A_x0(:, :, j)))
            q0_star(:, j) = matmul(A_0_h, q0_star(:, j))
            q0_star(:, j) = q0_star(:, j)/this%M_x0(:, j)
        end do

       !! set BCs lpse
        q0_star(2:5, 1) = 0.0d0
        q0_star(1:5, this%jn) = 0.0d0

       !! normal_with disturbance
        call ArrayNorm(this%jn, this%grid_x0%y, q0_star, this%M_x0, E)
        q0_star = q0_star/sqrt(E)

        call DMDAVecGetArrayF90(this%dm, this%q0, tmp, ierr)
        tmp = q0_star
        call DMDAVecRestoreArrayF90(this%dm, this%q0, tmp, ierr)

    end subroutine get_lpse_inlet

    subroutine get_alpse_inlet(this)

        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this
        PetscScalar, pointer :: tmp(:, :)

        PetscScalar :: q1(5, this%jn)
        PetscScalar :: A_0_h(5, 5)
        integer :: info, ipiv(5)
        integer :: j

        call DMDAVecGetArrayReadF90(this%dm, this%q1, tmp, ierr)
        q1 = tmp
        call DMDAVecRestoreArrayReadF90(this%dm, this%q1, tmp, ierr)

        do j = 1, this%jn
            A_0_h = conjg(transpose(this%A_x1(:, :, j)))
            call zgesv(5, 1, A_0_h, 5, ipiv, q1(:, j), 5, info)
            q1(:, j) = this%M_x1(:, j)*q1(:, j)
        end do

        ! set_BCs alpse
        q1(1:5, 1) = 0.0d0
        q1(1:5, this%jn) = 0.0d0
      !!
        call DMDAVecGetArrayF90(this%dm, this%q1_star, tmp, ierr)
        tmp = q1/maxval(abs(q1(2, :)))
        ! print*, maxval(abs(Tmp))
        ! read(*, *)
        call DMDAVecRestoreArrayF90(this%dm, this%q1_star, tmp, ierr)
        ! call VecView(this%q1_star, PETSC_VIEWER_STDOUT_SELF, ierr)
        ! read(*, *)

    end subroutine get_alpse_inlet

    subroutine finalize_optimal_streak_eqn(this)

        implicit none
        class(type_optimal_streak_eqn), intent(inout) :: this

        deallocate (this%A_x0)
        deallocate (this%A_x1)
        deallocate (this%M_x0)
        deallocate (this%M_x1)

        call DMRestoreGlobalVector(this%dm, this%q0, ierr)
        call DMRestoreGlobalVector(this%dm, this%q1, ierr)
        call DMRestoreGlobalVector(this%dm, this%q0_star, ierr)
        call DMRestoreGlobalVector(this%dm, this%q1_star, ierr)

        call DMDestroy(this%dm, ierr)

    end subroutine finalize_optimal_streak_eqn

    subroutine ArrayNorm(jn, y, q, M, Norm)

        implicit none
        integer, intent(in) :: jn
        real(REAL64), intent(in) :: y(jn)
        complex(REAL64), intent(in) :: q(5, jn)
        real(REAL64), intent(in) :: M(5, jn)
        real(real64), intent(out) :: Norm
        integer :: j
        real(real64) :: energy_local(jn)

        do j = 1, jn
            energy_local(j) = dot_product(q(:, j), M(:, j)*q(:, j))
        end do
        Norm = dIntergal(jn, energy_local, y)

    end subroutine ArrayNorm

    !> 实数积分函数
    function dIntergal(n, Var, Y) Result(IntValue)

        implicit none

        integer, intent(in) :: n !< 离散点数量
        real(REAL64), intent(in) :: Var(n) !< 被积函数的函数值
        real(REAL64), intent(in) :: Y(n) !< 被积函数函数值对应的自变量
        real(REAL64) :: IntValue !< 积分值
        integer :: i

        IntValue = 0.0d0
        do i = 1, n - 1
            IntValue = IntValue + 0.5d0*(Var(i) + Var(i + 1))*(Y(i + 1) - Y(i))
        end do

    end function dIntergal

end module mod_optimal_streak_eqn
