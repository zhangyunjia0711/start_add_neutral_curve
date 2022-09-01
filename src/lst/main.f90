!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  File: main.f90
!> @file
!> @breif main procedure of PreComp.
!  DESCRIPTION:
!>
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

program LST

  use petsc !, only: comm=>PETSC_COMM_WORLD
  use slepcsys
  use global_parameters, only:rank
  use mod_cfgio_adapter
  use mod_petsc_loader
  use mod_metrics
  use mod_flow_grad
  use mod_lst

  implicit none

  character(len=256) :: cfg_file

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

  call load_info(PETSC_COMM_WORLD)

  call metrics(PETSC_COMM_WORLD)

  call flow_grad(PETSC_COMM_WORLD)

  call solve_lst()

  call SLEPCFinalize(ierr)

  ! CHKERRA(ierr)

end program LST
