
!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  Module: fft
!
!  DESCRIPTION:
!> @breif FFT(傅里叶变换).
!>
!!
!  REVISION HISTORY:
!  yyyy-mm-dd - Initial Version
!  TODO_yyyy-mm-dd - TODO_describe_appropriate_changes - TODO_name
!> @author
!> Liu Jianxin
!> @date 2017-11-10
!------------------------------------------------------------------------------
#ifdef MKL
#include <mkl_dfti.f90>

module mod_fft

   use iso_fortran_env, only: REAL64
   use MKL_DFTI
   ! use global_parameters
   implicit none

   private

   integer :: status
   ! integer, parameter :: REAL64=8

   type, public :: fft_type

      private
      type(DFTI_DESCRIPTOR), pointer :: handle
      integer :: n1, n2
      integer, public :: rstrides(3)
      integer, public :: cstrides(3)
      logical, public :: isInitial
      logical :: isRZ

   Contains

      !Procedure :: Initial
      procedure :: Finalize
      procedure :: HasInitial

      generic :: fft => fft_zz, fft_zz_inplace, fft_zz_2d
      generic :: ifft => ifft_zz, ifft_zz_inplace, ifft_zz_2d
      generic :: Initialize => Initial_1d, initial_2d

      Procedure, private :: Initial_1d
      procedure, private :: initial_2d
      ! procedure, private :: fft_rz
      procedure, private :: fft_zz
      procedure, private :: fft_zz_inplace
      procedure, private :: fft_zz_2d
      ! procedure, private :: ifft_zr
      procedure, private :: ifft_zz
      procedure, private :: ifft_zz_2d
      procedure, private :: ifft_zz_inplace

   end type fft_type

contains

   function HasInitial(this) result(Has)

      implicit none
      class(fft_type), intent(inout) :: this
      logical :: Has

      has = this%isInitial

   end function HasInitial

   subroutine initial_1d(this, n)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) :: n

      call this%finalize()

      status = DftiCreateDescriptor(this%handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n)
      status = DftiSetValue(this%handle, DFTI_BACKWARD_SCALE, 1.0_REAL64/real(n, REAL64))

      status = DftiCommitDescriptor(this%handle)

      this%isInitial = .True.
      this%n1=n

   end subroutine initial_1d

   subroutine initial_2d(this, n1, n2)

      implicit none
      class(fft_type), intent(inout) :: this
      integer, intent(in) :: n1, n2

      call this%finalize()

      status = DftiCreateDescriptor(this%handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, [n1, n2])
      status = DftiSetValue(this%handle, DFTI_FORWARD_SCALE, 1.0_REAL64/real(n1*n2, REAL64))
      ! status = DftiSetValue(this%handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      ! status = DftiSetValue(this%handle, DFTI_CONJUGATE_EVEN_STORAGE,                   &
      !   &                    DFTI_COMPLEX_COMPLEX)

      status = DftiCommitDescriptor(this%handle)

      ! this%rstrides = [0, 1, n1]
      ! this%cstrides = [0, 1, INT(n1/2) + 1]
      this%n1=n1
      this%n2=n2

      ! status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%rstrides)
      ! status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%cstrides)
      ! status = DftiCommitDescriptor(this%handle)
      ! this%isRZ = .True.
      this%isInitial = .True.

   end subroutine initial_2d

   subroutine fft_zz_2d(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: Ain(1:this%n1, 1:this%n2)
      complex(REAL64), intent(out) :: Aout(1:this%n1, 1:this%n2)
      integer :: dim_Ain(2), dim_Aout(2)

      if (.not. this%isRZ) then
         ! status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%rstrides)
         ! status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%cstrides)
         ! status = DftiCommitDescriptor(this%handle)
         this%isRZ = .True.
      end if

      dim_Ain = shape(Ain); dim_Aout = shape(Aout)

   !    if (dim_Ain(1) /= this%rstrides(3)) then
   !       write (*, *) 'the first rank in real Ain is wrong: n1=', dim_Ain(1), &
   !   &   'not', this%rstrides(3)
   !       stop
   !    end if
   !    if (dim_Aout(1) /= this%cstrides(3)) then
   !       write (*, *) 'the first rank in complex Aout is wrong: n1=', dim_Aout(1), &
   !   &   'not', this%cstrides(3)
   !       stop
   !    end if
      Aout=Ain
      status = DftiComputeForward(this%handle, Aout(:, 1))

   end subroutine fft_zz_2d

   subroutine ifft_zz_2d(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: Ain(1:this%n1, 1:this%n2)
      complex(REAL64), intent(out) :: Aout(1:this%n1, 1:this%n2)
      integer :: dim_Ain(2), dim_Aout(2)
 
      ! print*, 'n1=', this%n1, 'n2=', this%n2

      if (this%isRZ) then
         ! status = DftiSetValue(this%handle, DFTI_INPUT_STRIDES, this%cstrides)
         ! status = DftiSetValue(this%handle, DFTI_OUTPUT_STRIDES, this%rstrides)
         ! status = DftiCommitDescriptor(this%handle)
         this%isRZ = .False.
      end if

      dim_Ain = shape(Ain); dim_Aout = shape(Aout)

      ! block
      !    integer :: l1, l2
      ! print*, 'A_in:'
      ! do l1 = 1, ubound(Aout, dim=1)
      !    do l2=1, ubound(Aout, dim=2)
      !       print*, l1, l2, Aout(l1, l2)
      !    enddo
      ! end do

      Aout=Ain

      status = DftiComputeBackward(this%handle, Aout(:, 1))
      ! print*, 'A_out:'
      ! do l1 = 1, ubound(Aout, dim=1)
      !    do l2=1, ubound(Aout, dim=2)
      !       print*, l1, l2, Aout(l1, l2)
      !    enddo
      ! end do
      ! endblock

   end subroutine ifft_zz_2d

   subroutine fft_zz(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: Ain(:)
      complex(REAL64), intent(out) :: Aout(:)

      Aout = Ain
      status = DftiComputeBackward(this%handle, Aout)

   end subroutine fft_zz

   subroutine ifft_zz(this, Ain, Aout)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(in) :: Ain(:)
      complex(REAL64), intent(out) :: Aout(:)

      Aout = Ain
      status = DftiComputeForward(this%handle, Aout)

   end subroutine ifft_zz

   subroutine fft_zz_inplace(this, Ain)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(inout) :: Ain(:)

      status = DftiComputeBackward(this%handle, Ain)

   end subroutine fft_zz_inplace

   subroutine ifft_zz_inplace(this, Ain)

      implicit none
      class(fft_type), intent(inout) :: this
      complex(REAL64), intent(inout) :: Ain(:)

      status = DftiComputeForward(this%handle, Ain)

   end subroutine ifft_zz_inplace

   ! subroutine fft_rz(this, Ain, Aout)

   !    implicit none
   !    class(fft_type), intent(inout) :: this
   !    real(REAL64), intent(in) :: Ain(:)
   !    complex(REAL64), intent(out) :: Aout(:)

   !    Aout = Ain
   !    status = DftiComputeBackward(this%handle, Aout)

   ! end subroutine fft_rz

   ! subroutine ifft_zr(this, Ain, Aout)

   !    implicit none
   !    class(fft_type), intent(inout) :: this
   !    complex(REAL64), intent(in) :: Ain(:)
   !    real(REAL64), intent(out) :: Aout(:)
   !    complex(REAL64), allocatable :: tmp(:)

   !    if (.not. allocated(tmp)) then
   !       allocate (tmp, source=Ain)
   !    else
   !       tmp = Ain
   !    end if
   !    status = DftiComputeForward(this%handle, tmp)
   !    Aout = real(tmp, REAL64)

   ! end subroutine ifft_zr

   subroutine finalize(this)
      implicit none
      class(fft_type), intent(inout) :: this

      status = DftiFreeDescriptor(this%handle)
      this%isInitial = .False.

   end subroutine finalize

end module mod_fft
#endif