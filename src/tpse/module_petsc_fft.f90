module mod_petsc_fft
#include <petsc/finclude/petsc.h>

#ifndef MKL
   use petsc
   use iso_fortran_env, only: REAL64
   implicit none
   private

   public :: fft_type

   type :: fft_type

      Mat :: handle
      Vec :: forward_in
      Vec :: forward_out
      Vec :: backward_out
      logical :: initialized = .FALSE.
      integer :: shape(3) = 1

   contains

      generic :: initialize => initialize_1d, initialize_2d
      generic :: fft => fft_1d, fft_2d
      generic :: ifft => ifft_1d, ifft_2d
      procedure :: finalize

      procedure, private :: initialize_1d
      procedure, private :: initialize_2d

      procedure, private :: fft_1d
      procedure, private :: fft_2d
      procedure, private :: ifft_1d
      procedure, private :: ifft_2d

   end type fft_type

   PetscErrorCode :: ierr

contains

   subroutine initialize_1d(this, n1)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1

      call MatCreateFFT(PETSC_COMM_SELF, 1, [n1], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1

   end subroutine initialize_1d

   subroutine initialize_2d(this, n1, n2)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1, n2

      call MatCreateFFT(PETSC_COMM_SELF, 2, [n1, n2], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1; this%shape(2) = n2

   end subroutine initialize_2d

   subroutine fft_1d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:)
      complex(REAL64), intent(inout) :: vec_out(:)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      call VecGetArrayF90(this%forward_in, array_in, ierr)
      array_in(:) = vec_in(:)
      call VecRestoreArrayF90(this%forward_in, array_in, ierr)

      call MatMult(this%handle, this%forward_in, this%forward_out, ierr)

      call VecGetArrayReadF90(this%forward_out, array_out, ierr)
      vec_out(:) = array_out(:)/real(product(this%shape), kind=REAL64)
      call VecRestoreArrayReadF90(this%forward_out, array_out, ierr)

   end subroutine fft_1d

   subroutine fft_2d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :)
      complex(REAL64), intent(inout) :: vec_out(:, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      call VecGetArrayF90(this%forward_in, array_in, ierr)
      array_in(:) = pack(vec_in(:, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_in, array_in, ierr)

      call MatMult(this%handle, this%forward_in, this%forward_out, ierr)

      call VecGetArrayReadF90(this%forward_out, array_out, ierr)
      vec_out(:, :) = reshape(array_out(:)/real(product(this%shape), kind=REAL64), this%shape(1:2))
      call VecRestoreArrayReadF90(this%forward_out, array_out, ierr)

   end subroutine fft_2d

   subroutine ifft_1d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:)
      complex(REAL64), intent(inout) :: vec_out(:)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      call VecGetArrayF90(this%forward_out, array_in, ierr)
      array_in(:) = vec_in(:)
      call VecRestoreArrayF90(this%forward_out, array_in, ierr)

      call MatMultTranspose(this%handle, this%forward_out, this%backward_out, ierr)

      call VecGetArrayReadF90(this%backward_out, array_out, ierr)
      vec_out(:) = array_out(:)
      call VecRestoreArrayReadF90(this%backward_out, array_out, ierr)

   end subroutine ifft_1d

   subroutine ifft_2d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :)
      complex(REAL64), intent(inout) :: vec_out(:, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      call VecGetArrayF90(this%forward_out, array_in, ierr)
      array_in(:) = pack(vec_in(:, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_out, array_in, ierr)

      call MatMultTranspose(this%handle, this%forward_out, this%backward_out, ierr)

      call VecGetArrayReadF90(this%backward_out, array_out, ierr)
      vec_out(:, :) = reshape(array_out(:), this%shape(1:2))
      call VecRestoreArrayReadF90(this%backward_out, array_out, ierr)

   end subroutine ifft_2d

   subroutine finalize(this)

      class(fft_type), intent(inout) :: this

      call MatDestroy(this%handle, ierr)
      call VecDestroy(this%forward_in, ierr)
      call VecDestroy(this%forward_out, ierr)
      call VecDestroy(this%backward_out, ierr)
      this%initialized = .false.
      this%shape = 1

   end subroutine finalize
#endif
end module mod_petsc_fft
