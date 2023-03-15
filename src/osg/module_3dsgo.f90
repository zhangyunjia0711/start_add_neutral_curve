!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-29 15:17:05
!> @LastEditTime: 2023-03-15 11:41:36
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/osg/module_3dsgo.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_3dsog

   use global_parameters
   use mod_lst_baseflow_profile
   use mod_lst_grid_profile
   ! use mod_lpse_eqn
   use mod_optimal_streak_eqn

   implicit none
   integer :: jb, je

   type(type_optimal_streak_eqn) :: solver
   type(type_lst_grid_profile) :: grid_x0, grid_x1
   type(type_lst_bf_profile) :: baseflow_x0, baseflow_x1

   private
   public :: solve_3dogs_new

contains

   subroutine solve_3dogs_new()
      implicit none
      real(REAL64) :: G_factor

      call initial_3dsgo()

      ! grid_x0=assemble_lst_grid_profile(i_x0)
      ! grid_x1=assemble_lst_grid_profile(i_x1)
      ! baseflow_x0=assemble_lst_baseflow_profile(i_x0)
      ! baseflow_x1=assemble_lst_baseflow_profile(i_x1)
      ! call solver%setup(i_x0, grid_x0, baseflow_x0, i_x1, grid_x1, baseflow_x1, dcmplx(initial_alf), initial_bat, initial_ome)
      ! call solver%set_ic()
      ! call solver%solve(G_factor)
      call solve_3dsog()

      call finalize_3dsgo()

   end subroutine solve_3dogs_new

   subroutine initial_3dsgo()

      implicit none

      jb = yb + 1; je = ye + 1

      call solver%initial(jn, Re, Ma, Te, Pr, yb, ye, isothermal)

   end subroutine initial_3dsgo

   subroutine finalize_3dsgo()

      implicit none

      jb = yb + 1; je = ye + 1

      call solver%finalize()
      call grid_x0%finalize()
      call grid_x1%finalize()
      call baseflow_x0%finalize()
      call baseflow_x1%finalize()

   end subroutine finalize_3dsgo

   subroutine solve_3dsog()

      implicit none
      real(REAL64) :: G_factor

      integer :: k, ki
      real(REAL64) :: pk
      Vec :: streak_eigf_x0, streak_eigf_x1

      grid_x0 = assemble_lst_grid_profile(i_x0)
      grid_x1 = assemble_lst_grid_profile(i_x1)
      baseflow_x0 = assemble_lst_baseflow_profile(i_x0)
      baseflow_x1 = assemble_lst_baseflow_profile(i_x1)
      if(flg_cycle)  initial_bat = initial_cycle_m*drz(i_x0)
      print *, 'i_x0=', i_x0, 'i_x1=', i_x1
      if(flg_cycle)  print *, 'cycle_m=', initial_cycle_m
      print *, 'beta=', initial_bat
      print *, 'omega=', initial_ome
      print *, 'alpha=', initial_alf
      call solver%setup(i_x0, grid_x0, baseflow_x0, i_x1, grid_x1, baseflow_x1, dcmplx(initial_alf), initial_bat, initial_ome)
      call solver%set_ic()
      call solver%solve(G_factor)

      open (999, file='out/optimal_growth_streak_results.plt', access='append')
      if(.not. flg_cycle)then
         write (999, *) '#variables= i_x0, i_x1, alf, bat, G'
         write (999, *) i_x0, i_x1, initial_alf, initial_bat, G_factor
      else
         write (999, *) '#variables= i_x0, i_x1, alf, cycle_m, G'
         write (999, *) i_x0, i_x1, initial_alf, initial_cycle_m, G_factor
      endif
      flush (999)

      call PetscPrintf(PETSC_COMM_WORLD, "==========Compute 3D optimal growth streak local!==========\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  k=-1: end 3D optimal growth streak.\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  k= 0: Output the optimal streak funciton at i_x0.\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  k= 1: Modify i_x0.\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  k= 2: Modify i_x1.\n", ierr)
      if(.not. flg_cycle)then
         call PetscPrintf(PETSC_COMM_WORLD, "===  k= 3: Modify beta.\n", ierr)
      else
         call PetscPrintf(PETSC_COMM_WORLD, "===  k= 3: Modify cycle_m.\n", ierr)
      endif
      call PetscPrintf(PETSC_COMM_WORLD, "===  k= 4: Modify omega.\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  k= 5: Modify alpha.\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "===  delta_k=ki*pk.(eg. k, ki, pk=4, 20, 0.001 means delta_omega_r=20*0.001)\n", ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "======================================\n", ierr)

      k = 0; ki = 0; pk = 0.0d0
      do while (k /= -1)
         print *, 'i_x0=', i_x0, 'i_x1=', i_x1
         print *, 'beta=', initial_cycle_m*drz(i_x0)
         print *, 'omega=', initial_ome
         print *, 'alpha=', initial_alf
         write (*, '("k, ki, pk=", $)')
         read (*, *) k, ki, pk
         if (k == -1) then
            print *, 'quit!'
            exit
         elseif (k < 0 .or. k > 5) then
            print *, 'k value is wrong, and input other value!'
            cycle
         end if
         do while (ki >= 1 .or. k == 0)
            select case (k)
            case (0)
               print *, 'getting eigf form solver...'
               call solver%get_eigf(streak_eigf_x0, streak_eigf_x1)
               print *, 'writing the eigf'
               call write_optimal_streak_function(streak_eigf_x0, 'x0')
               call write_optimal_streak_function(streak_eigf_x1, 'x1')
               k = 999; ki = 0
            case (1) !! adjust i_x0
               i_x0 = i_x0 + int(pk)
            case (2) !! adjust i_x1
               i_x1 = i_x1 + int(pk)
            case (3) !! adjust beta
               if(flg_cycle)then
                  initial_cycle_m = initial_cycle_m + pk
                  initial_bat = initial_cycle_m * drz(i_x0)
               else
                  initial_bat = initial_bat + pk
               endif
            case (4) !! adjust omega
               initial_ome = initial_ome + pk
               ! initial_alf = initial_ome/c_phase
            case (5) !! adjust alpha
               initial_alf = initial_alf + pk
            end select

            print *, 'i_x0=', i_x0, 'i_x1=', i_x1
            if(flg_cycle)  print *, 'cycle_m=', initial_cycle_m
            print *, 'beta=', initial_bat
            print *, 'omega=', initial_ome
            print *, 'alpha=', initial_alf

            if (k < 500) then
               grid_x0 = assemble_lst_grid_profile(i_x0)
               grid_x1 = assemble_lst_grid_profile(i_x1)
               baseflow_x0 = assemble_lst_baseflow_profile(i_x0)
               baseflow_x1 = assemble_lst_baseflow_profile(i_x1)
               call solver%setup(i_x0, grid_x0, baseflow_x0, i_x1, grid_x1, baseflow_x1, &
               &    dcmplx(initial_alf), initial_bat, initial_ome)
               call solver%set_ic(solver%q0)
               ! call solver%set_ic()
               call solver%solve(G_factor)
               if(.not. flg_cycle)then
                  write (999, *) i_x0, i_x1, initial_alf, initial_bat, G_factor
               else
                  write (999, *) i_x0, i_x1, initial_alf, initial_cycle_m, G_factor
               endif
               flush (999)

               if(isNAN(G_factor)) then
                  error stop "G_factor is NaN, please recompute the problem."
               endif
               flush (999)
            end if
            ki = ki - 1
         end do
      end do

      close (999)

   end subroutine solve_3dsog

   subroutine write_optimal_streak_function(v, name)

      use stdlib_string_type
      use global_parameters
      implicit none
      Vec :: v
      character(len=2) :: name

      type(string_type) :: filename
      Vec :: v_all
      VecScatter :: ctx
      PetscScalar, pointer :: tmp_(:)
      Complex(REAL64) :: phi(6, jn)
      integer :: l, j
      real(REAL64) :: pe
      PetscViewer :: Viewer
      ! Vec :: tmpV

      ! select case (lst_mode)
      ! case (LST_TEMPORAL_MODE)
      !   ome=ome+(0.0d0, 1.0d-8)
      ! case (LST_SPATIAL_MODE)
      !   alf=alf+(0.0d0, 1.0d-8)
      ! end select

      call PetscPrintf(PETSC_COMM_WORLD, "Writing Optimal streak......\n", ierr)

      filename = 'optimal_streak_'//string_type(i_x0)//'_'//string_type(i_x1)//'_'//name
      call PetscPrintf(PETSC_COMM_WORLD, "the file's prefix is "//char(filename)//"......\n", ierr)

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

         ! maxu_loc=maxloc(abs(phi(2, :)))
         ! maxu=phi(2, maxu_loc(1))
         ! phi=phi/maxu
         pe = 1.0d0/(GAMMA*MA*MA)
         do j = 1, jn
            phi(6, j) = (solver%baseflow_x0%q0(1, j)*phi(5, j) + solver%baseflow_x0%q0(5, j)*phi(1, j))*pe
         end do

         call execute_command_line('mkdir -p out')
         open (1002, file='out/'//char(filename)//'.plt', access='append')
         write (1002, ('(3(A, E20.13, 3X))')) '#omega=', real(initial_ome), 'alpha=', real(initial_alf), 'beta=', real(initial_bat)
         write (1002, *) 'variables=y, |u|, u_r, u_i, |v|, v_r, v_i, &
   &                 |w|, w_r, w_i, |T|, T_r, T_i, |p|, p_r, p_i, &
   &                  |rho|, rho_r, rho_i'
         do j = 1, jn
            write (1002, '(19E20.7)') solver%grid_x0%y(j), &
               abs(phi(2, j)), phi(2, j)%re, phi(2, j)%im, &
               abs(phi(3, j)), phi(3, j)%re, phi(3, j)%im, &
               abs(phi(4, j)), phi(4, j)%re, phi(4, j)%im, &
               abs(phi(5, j)), phi(5, j)%re, phi(5, j)%im, &
               abs(phi(6, j)), phi(6, j)%re, phi(6, j)%im, &
               abs(phi(1, j)), phi(1, j)%re, phi(1, j)%im
         end do
         close (1002)
      end if

      call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'out/'//char(filename)//'.petsc', FILE_MODE_WRITE, &
                                 Viewer, ierr)
      call VecView(V, Viewer, ierr)
      call PetscViewerDestroy(Viewer, ierr)

      call VecScatterDestroy(ctx, ierr)
      call VecDestroy(v_all, ierr)

   end subroutine write_optimal_streak_function

end module mod_3dsog

