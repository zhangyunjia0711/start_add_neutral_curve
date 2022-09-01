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

program NPSE_PROG

   use petsc
   use slepcsys
   use global_parameters, only: rank
   use mod_cfgio_adapter
   use mod_petsc_loader
   use mod_metrics
   use mod_flow_grad
   use mod_npse

   implicit none

   character(len=256) :: cfg_file

   PetscBool :: set

   call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
   CHKERRA(ierr)

   call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-f', cfg_file, set, ierr)
   CHKERRA(ierr)

   if (.not. set) then
      print *, 'should use -f option to determin the config file.'
      stop
   end if

   call cfg_loader(trim(cfg_file))

   call load_info(PETSC_COMM_WORLD)

   call metrics(PETSC_COMM_WORLD)

   call flow_grad(PETSC_COMM_WORLD)

   call npse()

   call SLEPCFinalize(ierr)

   CHKERRA(ierr)

end program NPSE_PROG
