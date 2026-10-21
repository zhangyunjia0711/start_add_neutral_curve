!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-09 10:45:42
!> @LastEditTime: 2022-12-12 19:33:48
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lst/module_lst.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lst

    use global_parameters
    use mod_lst_baseflow_profile
    use mod_lst_grid_profile
    use mod_lst_solver_abstract
    use mod_lst_eqn
    use mod_neutral_curve, only: run_neutral_curve, neutral_batch_enabled, neutral_exit_code
    use mod_neutral_numerics, only: NC_SUCCESS, NC_CLOSED, NC_BOUNDARY
    use mod_alst_eqn
    use stdlib_string_type, only: char

    implicit none

    private
    public :: solve_lst

#include "../include/start.f90"

contains

    subroutine solve_lst()

        implicit none

        if (lst_IsAdjoint) then
            solname = string_type('ALST')
        else
            solname = string_type('LST')
        end if

        ! select case (lst_method)
        ! case (LST_GLOBAL)
        ! call solve_lst_global()
        ! case (LST_LOCAL)
        call solve_lst_local()
        ! case (LST_DOMAIN)
        ! call solve_lst_domain()
        ! end select

    end subroutine solve_lst

    subroutine solve_lst_local()

        implicit none
        complex(REAL64) :: ome, alf, bat
        type(type_lst_grid_profile) :: grid_profile
        type(type_lst_bf_profile) :: baseflow_profile
        class(type_lst_abstract), allocatable :: lst
        Vec :: eigf

        integer :: k, ki, input_status, curve_status
        real(REAL64) :: pk
        !integer, parameter :: comm=PETSC_COMM_WORLD
        if (lst_IsAdjoint) then
            ! call initialize_alst_eqn()
            allocate (type_alst_eqn :: lst)
        else
            allocate (type_lst_eqn :: lst)
        end if
        call lst%set_parameter("lst_mode", lst_mode)
        call lst%set_parameter("Coordinate_System", Coordinate_System)
        call lst%initial(jn, Re=Re, Ma=Ma, Te=Te, yb_=yb, ye_=ye)
        call lst%set_parameter("Pr", Pr)
        call lst%set_parameter("bc_y0", bc_y0)
        call lst%set_parameter("bc_yn", bc_yn)
        lst%isothermal = isothermal

        k = 0; ki = 0; pk = 0.0d0
        if (rank == 0) then
            call execute_command_line('mkdir -p out')
            open (1001, file='./out/'//char(solname)//'_LOCAL_eigenvalues.plt', access='append')
            write (1001, *) '#variables=i, x, Ma, Re, Te, omega_r, omega_i, &
    &                alpha_r, alpha_i, beta_r, beta_i, c_r, c_i'
        end if

        call PetscPrintf(PETSC_COMM_WORLD, "==========Compute "//char(solname)//" local!==========\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=-1: end "//char(solname)//".\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 0: output eigen function.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 1: Modify local Re.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 2: Modify alpha_r.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 3: Modify beta_r.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 4: Modify omega_r.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 5: Modify alpha_i.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 6: Modify beta_i.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 7: Modify omega_i.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k= 8: Modify iloc.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=31: neutral curves across stations (safeguarded secant).\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=32: neutral curves across stations (bisection).\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=33: growth contours and neutral curves.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=36: arclength neutral curve at a station.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  k=37: neutral branches in a parameter plane.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  neutral settings: [neutral] in the input config file.\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "===  delta_k=ki*pk.(eg. k, ki, pk=4, 20, 0.001 means delta_omega_r=20*0.001)\n", ierr)
        call PetscPrintf(PETSC_COMM_WORLD, "======================================\n", ierr)

        grid_profile = assemble_lst_grid_profile(iloc)
        baseflow_profile = assemble_lst_baseflow_profile(iloc)
        ome = omega(1) + ZI*omega(2)
        alf = alpha(1) + ZI*alpha(2)
        bat = beta(1) + ZI*beta(2)

        call lst%set_grid(grid_profile)
        call lst%set_baseflow(baseflow_profile)
        call lst%setup()

        ! call grid_profile%print("grid.plt")
        ! call baseflow_profile%print("baseflow.plt")
        if (rank == 0) then
            write (*, "(2X,'iloc=',I4, 6(2X, A, F17.5))") iloc, 'x=', xx(iloc, ye + 1), &
                'drx=', drx(iloc), 'drz=', drz(iloc)
            write (*, "(6(2X, A, F17.5))") 'Re=', Re, 'Ma=', Ma, 'Te=', Te
            write (*, "((2X, A, F17.5,'+', F17.5, 'i'))") 'omega=', ome
            write (*, "((2X, A, F17.5,'+', F17.5, 'i'))") 'alpha=', alf
            write (*, "((2X, A, F17.5,'+', F17.5, 'i'))") 'beta =', bat
        end if
        ! if(lst_IsAdjoint)then
        !   call solve_alst_eqn(grid_profile, baseflow_profile, ome, alf, bat, iloc)
        ! else
        call lst%solve(ome, alf, bat)
        ! endif
        if (rank == 0) then
            write (1001, ('(I0, 12(2X, F20.8))')) iloc, xx(iloc, yb + 1), Ma, Re, Te, ome, &
                alf, bat, ome/alf
            call flush (1001)
        end if

        if (neutral_batch_enabled()) then
            call run_neutral_curve(lst, ome, alf, bat, result_status=curve_status)
            if (.not.any(curve_status==[NC_SUCCESS,NC_CLOSED,NC_BOUNDARY])) neutral_exit_code=1
            k=-1
        end if
        do while (k /= -1)
            if (rank == 0) then
                write (*, '("k, ki, pk=", $)')
                read (*, *, iostat=input_status) k, ki, pk
                if (input_status/=0) then
                    k=-1; ki=0; pk=0.0d0
                end if
            end if
            call MPI_bcast(k, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(ki, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(pk, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            if (k==-1) exit
            if (any(k==[31,32,33,36,37])) then
                call run_neutral_curve(lst, ome, alf, bat, task_override=k)
                cycle
            end if
            do while (ki >= 1 .or. k == 0)
                select case (k)
                case (0)
                    ! if(lst_IsAdjoint)then
                    !   call write_alst_function(grid_profile, baseflow_profile, ome, alf, bat, iloc)
                    ! else
                    ! iloc = iloc + int(pk)
                    ! grid_profile = assemble_lst_grid_profile(iloc)
                    ! baseflow_profile = assemble_lst_baseflow_profile(iloc)
                    ! call lst%set_grid(grid_profile)
                    ! call lst%set_baseflow(baseflow_profile)
                    ! call lst%setup()
                    ! endif
                    ! call baseflow_profile%print('baseflow_profile.plt')
                    k = 999; ki = 0
                case (1)
                    Re = Re + pk
                    call lst%set_parameter("Re", Re)
                    call lst%setup()
                case (2)
                    alf = alf + pk
                case (3)
                    bat = bat + pk
                case (4)
                    ome = ome + pk
                case (5)
                    alf = alf + pk*ZI
                case (6)
                    bat = bat + pk*ZI
                case (7)
                    ome = ome + pk*ZI
                case (8)
                    iloc = iloc + int(pk)
                    grid_profile = assemble_lst_grid_profile(iloc)
                    baseflow_profile = assemble_lst_baseflow_profile(iloc)
                    call lst%set_grid(grid_profile)
                    call lst%set_baseflow(baseflow_profile)
                    call lst%setup()
                end select

                if (rank == 0) then
                    write (*, "(2X,'iloc=',I4, 6(2X, A, F17.5))") iloc, 'x=', xx(iloc, ye + 1), &
                        'drx=', drx(iloc), 'drz=', drz(iloc)
                    write (*, "(6(2X, A, F17.5))") 'Re=', Re, 'Ma=', Ma, 'Te=', Te
                    write (*, "((2X, A, F17.8,'+', F17.8, 'i'))") 'omega=', ome
                    write (*, "((2X, A, F17.8,'+', F17.8, 'i'))") 'alpha=', alf
                    write (*, "((2X, A, F17.8,'+', F17.8, 'i'))") 'beta =', bat
                end if

                if (k < 500) then
                    ! if(lst_IsAdjoint)then
                    !   call solve_alst_eqn(grid_profile, baseflow_profile, ome, alf, bat, iloc)
                    ! else
                    call lst%solve(ome, alf, bat)
                    ! endif
                    if (rank == 0) then
                        write (1001, ('(I0, 12(2X, F20.8))')) iloc, xx(iloc, yb + 1), Ma, Re, Te, ome, &
                            alf, bat, ome/alf
                        call flush (1001)
                    end if
                else
                  !!!!!! 不知道为什么连续调用两次求解特征值就会报错，说主对角线值为零 
                  !!!!!! 所以考虑到求解问题的特点干脆我就直接关掉再复求特征值的部分~
                    call lst%get_eigf(eigf)
                    call write_lst_function(ome, alf, bat, iloc, eigf)
                end if
                ki = ki - 1
            end do
            !call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, iloc)
        end do

        call PetscPrintf(PETSC_COMM_WORLD, "LST_LOCAL is done.\n", ierr)
        ! if(lst_IsAdjoint)then
        ! call finalize_alst_eqn()!'
        ! endif

        if (rank==0) close(1001)
        call lst%finalize()
        call grid_profile%finalize()
        call baseflow_profile%finalize()

    end subroutine solve_lst_local

    ! subroutine solve_lst_global()

    !   implicit none
    !   type(type_lst_grid_profile) :: grid_profile
    !   type(type_lst_bf_profile) :: baseflow_profile
    !   complex(REAL64) :: ome, alf, bat

    !   call initialize_lst_eqn()

    !   if(rank==0)then
    !       call execute_command_line('mkdir -p out')
    !       open(1001, file='./out/'//solname//'_GLOBAL_eigenvalues.plt', access='append')
    !       write(1001, *)'#variables=i, x, Ma, Re, Te, omega_r, omega_i, &
    !               alpha_r, alpha_i, beta_r, beta_i, c_r, c_i'
    !   endif

    !   grid_profile = assemble_lst_grid_profile(iloc)
    !   baseflow_profile = assemble_lst_baseflow_profile(iloc)

    !   ome=omega(1)+ZI*omega(2)
    !   alf=alpha(1)+ZI*alpha(2)
    !   bat=beta (1)+ZI*beta (2)

    !   call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, iloc)

    !   call grid_profile%finalize()
    !   call baseflow_profile%finalize()

    ! end subroutine solve_lst_global

    !@TODO should be rewrite for the class for additional features.
! block

!   subroutine solve_lst_domain()

!     implicit none

!     if(lst_IsAdjoint)then
!       call initialize_alst_eqn()
!     else
!       call initialize_lst_eqn()
!     endif

!     ! call solve_neutral_curve()
!     call solve_integral_N()

!     print*, 'lst_domain done!'

!   end subroutine solve_lst_domain

!   subroutine solve_integral_N()

!     implicit none
!     integer, parameter :: iome_b=1229, iome_e=1675
!     integer, parameter :: nline=100

!     type(type_lst_grid_profile) :: grid_profile
!     type(type_lst_bf_profile) :: baseflow_profile
!     complex(REAL64) :: alf_contours((iome_b+1)/2:(iome_e+1)/2, istart:iend)
!     real(REAL64) :: N_contours((iome_b+1)/2:(iome_e+1)/2, istart:iend)
!     integer :: iounit
!     integer :: line, iome, i, iome_
!     real(REAL64) :: omer, omei, alfr, alfi, batr, bati
!     complex(REAL64) :: ome, alf, bat
!     iome_=0

!     alf_contours=0.0d0; N_contours=0.0d0

!     open(unit=iounit, file='result3.txt')
!     do line=1, nline
!       read(iounit, *)iloc, omer, omei, alfr, alfi, batr, bati
!       ome=cmplx(omer, omei); alf=cmplx(alfr, alfi); bat=cmplx(batr, bati)
!       iome=(int(omer*1000)+1)/2
!       if(iome==iome_) then
!         cycle
!       else
!         iome_=iome
!       endif
!       print*, line, iome, iloc, istart
!       do i=iloc, iend, 2
!         grid_profile = assemble_lst_grid_profile(i)
!         baseflow_profile = assemble_lst_baseflow_profile(i)
!         alf%im=alf%im*1.05d0
!         call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)
!         alf_contours(iome, i)=alf
!         if((i/=istart) .and. (abs(aimag(alf_contours(iome, max(i-1, istart))))>epsilon(1.0d0)))then
!           N_contours(iome, i)=(-aimag(alf_contours(iome, i-1))-aimag(alf_contours(iome, i)))*0.5d0+N_contours(iome, i-1)
!         endif
!         call grid_profile%finalize()
!         call baseflow_profile%finalize()
!         if(i>10 .and. &
!         (maxval(aimag(alf_contours(iome, i:max(i-10, 1):-1)))-minval(aimag(alf_contours(iome, i:max(i-10, 1):-1))))<=1d-8)exit
!         if(N_contours(iome, i)<=-1.0d-8) exit
!       enddo
!     enddo
!     close(iounit)

!     if(rank==0)then
!       open(unit=iounit, file='contours3.plt')
!       write(iounit, *)'variables="<i>x</i>", "<greek>w</greek>", "-<greek>a</greek><sub>i</sub>", "<i>N</i>"'
!       write(iounit, *)'zone i=', (iend-istart)/2+1, 'j=', (iome_e-iome_b)/2+1
!       do iome=iome_b, iome_e, 2
!         do i=istart, iend, 2
!           write(iounit, *)real(i-1, REAL64), iome/1000.0d0, -alf_contours((iome+1)/2, i)%im, N_contours((iome+1)/2, i)
!         enddo
!       enddo
!       close(iounit)
!     endif

!   end subroutine solve_integral_N

!   subroutine solve_neutral_curve()
!     implicit none
!     integer :: i
!     type(type_lst_grid_profile) :: grid_profile
!     type(type_lst_bf_profile) :: baseflow_profile
!     complex(REAL64) :: ome_, alf_, bat_
!     complex(REAL64) :: ome, alf, bat
!     complex(REAL64) :: ome_max, alf_max, bat_max
!     real(REAL64) :: gr_max

!     real(REAL64), parameter :: delta_ome = 0.002d0

!     ome=omega(1)+ZI*omega(2)
!     alf=alpha(1)+ZI*alpha(2)
!     bat=beta (1)+ZI*beta (2)

!     gr_max=0.0d0

!     do i=iloc, istart, -1

!       if(rank==0) print*, 'solving lst at location i=', i
!       grid_profile = assemble_lst_grid_profile(i)
!       baseflow_profile = assemble_lst_baseflow_profile(i)

!       ome_=ome; alf_=alf; bat_=bat

!       call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)

!       gr_max=-aimag(alf)

!       do while (.not. (aimag(alf)>0))
!         ome=ome+delta_ome
!         call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)
!         if(-aimag(alf)>gr_max)then
!           gr_max=-aimag(alf)
!           alf_max=alf; bat_max=bat; ome_max=ome
!         endif
!       end do

!       ome=ome_; alf=alf_; bat=bat_

!       do while (.not. (aimag(alf)>0))

!         ome=ome-delta_ome

!         call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)

!         if(-aimag(alf)>gr_max)then
!           gr_max=-aimag(alf)
!           alf_max=alf; bat_max=bat; ome_max=ome
!         endif
!       end do

!       if(abs(gr_max)<epsilon(1.0d0)) then
!         print*, 'stop at', abs(gr_max)
!         exit
!       endif

!       ome=ome_max; alf=alf_max; bat=bat_max; gr_max=0.0d0

!     enddo

!     ome=omega(1)+ZI*omega(2)
!     alf=alpha(1)+ZI*alpha(2)
!     bat=beta (1)+ZI*beta (2)

!     gr_max=0.0d0

!     do i=iloc+1, iend, 1

!       if(rank==0) print*, 'solving lst at location i=', i
!       grid_profile = assemble_lst_grid_profile(i)
!       baseflow_profile = assemble_lst_baseflow_profile(i)

!       ome_=ome; alf_=alf; bat_=bat

!       call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)

!       gr_max=-aimag(alf)

!       do while (.not. (aimag(alf)>0))
!         ome=ome+delta_ome
!         call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)
!         if(-aimag(alf)>gr_max)then
!           gr_max=-aimag(alf)
!           alf_max=alf; bat_max=bat; ome_max=ome
!         endif
!       end do

!       ome=ome_; alf=alf_; bat=bat_

!       do while (.not. (aimag(alf)>0))

!         ome=ome-delta_ome

!         call solve_lst_eqn(grid_profile, baseflow_profile, ome, alf, bat, i)

!         if(-aimag(alf)>gr_max)then
!           gr_max=-aimag(alf)
!           alf_max=alf; bat_max=bat; ome_max=ome
!         endif
!       end do

!       if(abs(gr_max)<epsilon(1.0d0)) then
!         print*, 'stop at', abs(gr_max)
!         exit
!       endif

!       ome=ome_max; alf=alf_max; bat=bat_max; gr_max=0.0d0

!     enddo

!     call grid_profile%finalize()
!     call baseflow_profile%finalize()

!   end subroutine solve_neutral_curve
! endblock

    subroutine write_lst_function(ome, alf, bat, ix, v)

        use stdlib_string_type
        use global_parameters
        implicit none
        complex(REAL64) :: alf, bat, ome
        integer :: ix
        Vec :: v

        type(string_type) :: fn
        character(255) :: tmp
        PetscScalar :: d
        Vec :: v_all
        VecScatter :: ctx
        PetscScalar, pointer :: tmp_(:)
        Complex(REAL64) :: phi(6, jn)
        character(40) :: name
        integer :: l, j, maxu_loc(1)
        real(REAL64) :: pe
        complex(REAL64):: maxu
        PetscViewer :: Viewer
        Vec :: tmpV

        ! select case (lst_mode)
        ! case (LST_TEMPORAL_MODE)
        !   ome=ome+(0.0d0, 1.0d-8)
        ! case (LST_SPATIAL_MODE)
        !   alf=alf+(0.0d0, 1.0d-8)
        ! end select

        call PetscPrintf(PETSC_COMM_WORLD, "Writing "//char(solname)//" Eigenfunction......\n", ierr)

        write (tmp, ('(I0)')) ix
        fn = solname//'_eigenfunction_i='//trim(adjustl(tmp))//'.plt'
        ! write(tmp, ('(F20.5)'))ome%re
        ! fn=fn//'_ome='//trim(adjustl(tmp))
        ! write(tmp, ('(F20.5)'))alf%re
        ! fn=fn//'_alf='//trim(adjustl(tmp))
        ! write(tmp, ('(F20.5)'))bat%re
        ! fn=fn//'_bat='//trim(adjustl(tmp))//'.plt'

        call VecScatterCreateToZero(v, ctx, v_all, ierr)
        ! call VecGetType(v, name, ierr)
        ! call VecGetType(v_all, name, ierr)
        ! call VecGetSize(v, l, ierr)
        ! call VecGetSize(v_all, l, ierr)

        call VecScatterBegin(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, v, v_all, INSERT_VALUES, SCATTER_FORWARD, ierr)

        if (rank == 0) then
            call VecGetArrayReadF90(v_all, tmp_, ierr)
            phi(1:5, 1:jn) = reshape(tmp_, [5, jn])
            call VecRestoreArrayReadF90(v_all, tmp_, ierr)

            maxu_loc = maxloc(abs(phi(2, :)))
            maxu = phi(2, maxu_loc(1))
            phi = phi/maxu
            pe = 1.0d0/(GAMMA*MA*MA)
            do j = 1, jn
                phi(6, j) = (qq(1, ix, j)*phi(5, j) + qq(5, ix, j)*phi(1, j))*pe
            end do

            ! block
            !   use mod_dx
            !   integer :: l1, l2
            !   complex(REAL64), dimension(5, jn) :: dphidj, dphidjj
            !   do l=1, 5
            !     call dx_fdm4nd_center(phi(l, :), ybg+1, yeg+1, dphidj(l, :), ybg+1, yeg+1)
            !     call dxx_fdm4nd_center(phi(l, :), ybg+1, yeg+1, dphidjj(l, :), ybg+1, yeg+1)
            !   enddo
            !
            !
            !   alf1=d
            !   print*, ome1, alf1
            !   ! print*, grid_profile%eta_y(149), grid_profile%eta_yy(149)
            !   ! print*, lbound(grid_profile%eta_y), maxu_loc
            !   print*, phi(:, 149)
            !   print*, dphidj(:, 149)
            !   print*, dphidjj(:, 149)
            !   do l1=1, 5
            !     do l2=1, 5
            !       print*, l1, l2, A0(l1, l2, 149)
            !     enddo
            !   enddo
            !   print*, matmul(-ZI*ome1*G0(:, :, 149)+ZI*alf1*A0(:, :, 149)+D0(:, :, 149)+ &
            !   alf1**2*V0xx(:, :, 149), phi(1:5, 149))+ &
            !   matmul(B0(:, :, 149)-ZI*alf1*V0xy(:, :, 149), dphidj(1:5, 149)*grid_profile%eta_y(149))+ &
            !   matmul(-V0yy(:,:,149), dphidj(:, 149)*grid_profile%eta_yy(149)+dphidjj(1:5, 149)*grid_profile%eta_y(149)**2)
            !
            !   print*, matmul(-ZI*ome1*G0(:, :, 149)+ZI*alf1*A0(:, :, 149)+D0(:, :, 149)+ &
            !   alf1**2*V0xx(:, :, 149), phi(1:5, 149))
            !   read(*, *)
            ! endblock

            call execute_command_line('mkdir -p out')
            open (1002, file='out/'//char(fn))
            write (1002, ('(3(A, 2E20.13, 3X))')) '#omega=', ome, 'alpha=', alf, 'beta=', bat
            write (1002, *) 'variables=y, |u|, u_r, u_i, |v|, v_r, v_i, &
      &                 |w|, w_r, w_i, |T|, T_r, T_i, |p|, p_r, p_i, &
      &                  |rho|, rho_r, rho_i, x, i'
            do j = 1, jn
                write (1002, '(20E20.7, 2X, I0)') yy(ix, j), &
                    abs(phi(2, j)), phi(2, j)%re, phi(2, j)%im, &
                    abs(phi(3, j)), phi(3, j)%re, phi(3, j)%im, &
                    abs(phi(4, j)), phi(4, j)%re, phi(4, j)%im, &
                    abs(phi(5, j)), phi(5, j)%re, phi(5, j)%im, &
                    abs(phi(6, j)), phi(6, j)%re, phi(6, j)%im, &
                    abs(phi(1, j)), phi(1, j)%re, phi(1, j)%im, &
                    xx(ix, yb + 1), ix
            end do
            close (1002)
        end if

        call MPI_bcast(maxu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        call VecScale(V, 1.0d0/maxu, ierr)

        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'out/disturbance.petsc', FILE_MODE_WRITE, &
                                   Viewer, ierr)
        call VecView(V, Viewer, ierr)
        call PetscViewerDestroy(Viewer, ierr)

        call VecScatterDestroy(ctx, ierr)
        call VecDestroy(v_all, ierr)
        ! call VecDestroy(V, ierr)

    end subroutine write_lst_function

end module mod_lst
