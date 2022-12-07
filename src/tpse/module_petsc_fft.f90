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
      integer :: shape(4) = 1

   contains

      generic :: initialize => initialize_1d, initialize_2d, initialize_3d, initialize_4d
      generic :: fft => fft_1d, fft_2d, fft_3d, fft_4d
      generic :: ifft => ifft_1d, ifft_2d, ifft_3d, ifft_4d
      procedure :: finalize

      procedure, private :: initialize_1d
      procedure, private :: initialize_2d
      procedure, private :: initialize_3d
      procedure, private :: initialize_4d

      procedure, private :: fft_1d
      procedure, private :: fft_2d
      procedure, private :: fft_3d
      procedure, private :: fft_4d
      procedure, private :: ifft_1d
      procedure, private :: ifft_2d
      procedure, private :: ifft_3d
      procedure, private :: ifft_4d

   end type fft_type

   PetscErrorCode :: ierr

contains

   subroutine initialize_1d(this, n1)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1

      print*, 'initial 1d fft'

      call MatCreateFFT(PETSC_COMM_SELF, 1, [n1], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1

   end subroutine initialize_1d

   subroutine initialize_2d(this, n1, n2)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1, n2

      print*, 'initial 2d fft'
      call MatCreateFFT(PETSC_COMM_SELF, 2, [n1, n2], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1; this%shape(2) = n2

   end subroutine initialize_2d

   subroutine initialize_3d(this, n1, n2, n3)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1, n2, n3

      print*, 'initial 3d fft'
      call MatCreateFFT(PETSC_COMM_SELF, 3, [n1, n2, n3], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1; this%shape(2) = n2; this%shape(3) = n3

   end subroutine initialize_3d

   subroutine initialize_4d(this, n1, n2, n3, n4)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) ::  n1, n2, n3, n4

      print*, 'initial 4d fft'

      call MatCreateFFT(PETSC_COMM_SELF, 4, [n1, n2, n3, n4], MATFFTW, this%handle, ierr)
      call MatCreateVecsFFTW(this%handle, this%forward_in, this%forward_out, this%backward_out, ierr)
      this%initialized = .TRUE.

      this%shape(1) = n1; this%shape(2) = n2; this%shape(3) = n3; this%shape(4) = n4

   end subroutine initialize_4d

   subroutine fft_1d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:)
      complex(REAL64), intent(inout) :: vec_out(:)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'fft 1d'
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

      print*, 'fft 2d'

      call VecGetArrayF90(this%forward_in, array_in, ierr)
      array_in(:) = pack(vec_in(:, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_in, array_in, ierr)

      call MatMult(this%handle, this%forward_in, this%forward_out, ierr)

      call VecGetArrayReadF90(this%forward_out, array_out, ierr)
      vec_out(:, :) = reshape(array_out(:)/real(product(this%shape), kind=REAL64), this%shape(1:2))
      call VecRestoreArrayReadF90(this%forward_out, array_out, ierr)

   end subroutine fft_2d

   subroutine fft_3d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :, :)
      complex(REAL64), intent(inout) :: vec_out(:, :, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'fft 3d'

      call VecGetArrayF90(this%forward_in, array_in, ierr)
      array_in(:) = pack(vec_in(:, :, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_in, array_in, ierr)

      call MatMult(this%handle, this%forward_in, this%forward_out, ierr)

      call VecGetArrayReadF90(this%forward_out, array_out, ierr)
      vec_out(:, :, :) = reshape(array_out(:)/real(product(this%shape), kind=REAL64), this%shape(1:3))
      call VecRestoreArrayReadF90(this%forward_out, array_out, ierr)

   end subroutine fft_3d

   subroutine fft_4d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :, :, :)
      complex(REAL64), intent(inout) :: vec_out(:, :, :, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'fft 4d'

      call VecGetArrayF90(this%forward_in, array_in, ierr)
      array_in(:) = pack(vec_in(:, :, :, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_in, array_in, ierr)

      call MatMult(this%handle, this%forward_in, this%forward_out, ierr)

      call VecGetArrayReadF90(this%forward_out, array_out, ierr)
      vec_out(:, :, :, :) = reshape(array_out(:)/real(product(this%shape), kind=REAL64), this%shape(1:4))
      call VecRestoreArrayReadF90(this%forward_out, array_out, ierr)

   end subroutine fft_4d

   subroutine ifft_1d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:)
      complex(REAL64), intent(inout) :: vec_out(:)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'ifft 1d'

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

      print*, 'ifft 2d'

      call VecGetArrayF90(this%forward_out, array_in, ierr)
      array_in(:) = pack(vec_in(:, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_out, array_in, ierr)

      call MatMultTranspose(this%handle, this%forward_out, this%backward_out, ierr)

      call VecGetArrayReadF90(this%backward_out, array_out, ierr)
      vec_out(:, :) = reshape(array_out(:), this%shape(1:2))
      call VecRestoreArrayReadF90(this%backward_out, array_out, ierr)

   end subroutine ifft_2d

   subroutine ifft_3d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :, :)
      complex(REAL64), intent(inout) :: vec_out(:, :, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'ifft 3d'

      call VecGetArrayF90(this%forward_out, array_in, ierr)
      array_in(:) = pack(vec_in(:, :, :), mask=.true.)
      call VecRestoreArrayF90(this%forward_out, array_in, ierr)

      call MatMultTranspose(this%handle, this%forward_out, this%backward_out, ierr)

      call VecGetArrayReadF90(this%backward_out, array_out, ierr)
      vec_out(:, :, :) = reshape(array_out(:), this%shape(1:3))
      call VecRestoreArrayReadF90(this%backward_out, array_out, ierr)

   end subroutine ifft_3d

   subroutine ifft_4d(this, vec_in, vec_out)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: vec_in(:, :, :, :)
      complex(REAL64), intent(inout) :: vec_out(:, :, :, :)
      complex(REAL64), pointer :: array_in(:), array_out(:)

      print*, 'ifft 4d'

      call VecGetArrayF90(this%forward_out, array_in, ierr)
      array_in(:) = pack(vec_in(:, :, :, :), mask=.true.)
      print*, shape(array_in)
      call VecRestoreArrayF90(this%forward_out, array_in, ierr)

      print*, shape(vec_in)

      call MatMultTranspose(this%handle, this%forward_out, this%backward_out, ierr)

      stop
      call VecGetArrayReadF90(this%backward_out, array_out, ierr)
      vec_out(:, :, :, :) = reshape(array_out(:), this%shape(1:4))
      call VecRestoreArrayReadF90(this%backward_out, array_out, ierr)

   end subroutine ifft_4d

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
