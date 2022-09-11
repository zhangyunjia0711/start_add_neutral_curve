#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_lpse_receptivity

    use global_parameters
    implicit none

    public :: lpse_receptivity

    complex(REAL64), allocatable :: lpse_alpha(:)
    Vec, allocatable :: lpse_function(:), alpse_function(:)

    
contains

    subroutine lpse_receptivity()

        implicit none

        call initialize_lpse_receptivity()
        call solve_lpse_receptivity()
        call finalize_lpse_receptivity()

    end subroutine lpse_receptivity

    subroutine initialize_lpse_receptivity()

        implicit none
        integer :: i

        allocate(lpse_alpha(istart:iend))
        allocate(lpse_function(istart:iend))
        allocate(alpse_function(istart:iend))

        do i=istart, iend
            call DMGetGlobalVector(meshDA, lpse_function(i), ierr)
            call DMGetGlobalVector(meshDA, alpse_function(i), ierr)
        enddo
        print*, 'initialize_lpse_receptivity.'
        

    end subroutine initialize_lpse_receptivity

    subroutine solve_lpse_receptivity()
        implicit none

        print*, 'solve_lpse_receptivity.'
        
    end subroutine solve_lpse_receptivity

    subroutine finalize_lpse_receptivity()

        implicit none
        integer :: i

        deallocate(lpse_alpha)
        do i=istart, iend
            call DMRestoreGlobalVector(meshDA, lpse_function(i), ierr)
            call DMRestoreGlobalVector(meshDA, alpse_function(i), ierr)
        enddo

        deallocate(lpse_function, alpse_function)

        print*, 'finalize_lpse_receptivity.'
        

    end subroutine finalize_lpse_receptivity
    
end module mod_lpse_receptivity