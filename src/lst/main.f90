!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-05 09:32:30
!> @LastEditTime: 2022-12-12 19:33:38
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: main program for lst.
!> @FilePath: /START_all/src/lst/main.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ LST solver.
program LST_prog

  use petsc !, only: comm=>PETSC_COMM_WORLD
  use slepcsys
  use global_parameters, only:rank
  use mod_cfgio_adapter
  use mod_petsc_loader
  use mod_metrics
  use mod_flow_grad
  use mod_lst
  use mod_neutral_curve, only: neutral_exit_code, neutral_batch_enabled

  implicit none

  character(len=256) :: cfg_file

  logical :: neutral_requested
  integer :: comm_size
  PetscBool :: set

  call SLEPCInitialize(PETSC_NULL_CHARACTER, ierr)
  ! CHKERRA(ierr)

  call mpi_comm_rank(PETSC_COMM_WORLD, rank, ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-f',cfg_file,set,ierr)
  ! CHKERRA(ierr)

  if(.not. set) then
    print*, 'should use -f option to determin the config file.'
    stop
  endif

  call cfg_loader(trim(cfg_file))
  neutral_requested=.false.
  if (rank==0) neutral_requested=neutral_batch_enabled()
  call MPI_Bcast(neutral_requested,1,MPI_LOGICAL,0,PETSC_COMM_WORLD,ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,comm_size,ierr)
  if (neutral_requested .and. comm_size/=1) then
    call PetscPrintf(PETSC_COMM_WORLD,"Neutral curves require a single MPI process.\n",ierr)
    call SLEPCFinalize(ierr)
    stop 1
  end if

  call load_info(PETSC_COMM_WORLD)

  call metrics(PETSC_COMM_WORLD)

  call flow_grad(PETSC_COMM_WORLD)

  call solve_lst()

  call SLEPCFinalize(ierr)
  if (neutral_exit_code/=0) stop 1

  ! CHKERRA(ierr)

end program LST_prog
