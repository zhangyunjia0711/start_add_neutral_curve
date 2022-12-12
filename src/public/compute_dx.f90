!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-05 09:32:33
!> @LastEditTime: 2022-12-12 16:16:41
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/compute_dx.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
module mod_dx

   use iso_fortran_env, only: REAL64
   implicit none

   public :: dx_fdm4nd_center, dxx_fdm4nd_center
   private

   interface dx_fdm4nd_center

      module procedure dx_fdm4nd_center_r8
      module procedure dx_fdm4nd_center_z8

   end interface

   interface dxx_fdm4nd_center

      module procedure dxx_fdm4nd_center_r8
      module procedure dxx_fdm4nd_center_z8

   end interface

   interface dx_fdm4nd_center_normal

      module procedure dx_fdm4nd_center_normal_r8
      module procedure dx_fdm4nd_center_normal_z8

   end interface

   interface dx_fdm4nd_center_head

      module procedure dx_fdm4nd_center_head_r8
      module procedure dx_fdm4nd_center_head_z8

   end interface

   interface dx_fdm4nd_center_tail

      module procedure dx_fdm4nd_center_tail_r8
      module procedure dx_fdm4nd_center_tail_z8

   end interface

   interface dx_fdm4nd_center_body

      module procedure dx_fdm4nd_center_body_r8
      module procedure dx_fdm4nd_center_body_z8

   end interface

   interface dxx_fdm4nd_center_normal

      module procedure dxx_fdm4nd_center_normal_r8
      module procedure dxx_fdm4nd_center_normal_z8

   end interface

   interface dxx_fdm4nd_center_head

      module procedure dxx_fdm4nd_center_head_r8
      module procedure dxx_fdm4nd_center_head_z8

   end interface

   interface dxx_fdm4nd_center_tail

      module procedure dxx_fdm4nd_center_tail_r8
      module procedure dxx_fdm4nd_center_tail_z8

   end interface

   interface dxx_fdm4nd_center_body

      module procedure dxx_fdm4nd_center_body_r8
      module procedure dxx_fdm4nd_center_body_z8

   end interface

contains

   pure subroutine dx_fdm4nd_center_r8(f, ib1, ie1, df, ib2, ie2)

      implicit none

      real(REAL64), intent(in) :: f(:)
      real(REAL64), intent(out) :: df(:)
      integer, intent(in) :: ib1, ie1, ib2, ie2

      integer :: n1, n2
      integer :: n

      n1 = ie1 - ib1 + 1; n2 = ie2 - ib2 + 1

      if (n1 == n2) then
         call dx_fdm4nd_center_normal(f, df, n1)
      elseif (n1 /= n2) then
         if ((ib1 == ib2) .and. (ie1 /= ie2)) then
            call dx_fdm4nd_center_head(f, df, n1)
         elseif ((ie1 == ie2) .and. (ib1 /= ib2)) then
            call dx_fdm4nd_center_tail(f, df, n1)
         elseif ((ib1 /= ib2) .and. (ie1 /= ie2) .and. ((ib2 - ib1) == (ie1 - ie2))) then
            call dx_fdm4nd_center_body(f, df, n1)
         ! else
         !    print *, 'there is some mistake with dx_fdm_4nd'
         end if
      end if

   end subroutine dx_fdm4nd_center_r8

   pure subroutine dx_fdm4nd_center_normal_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(n)
      real(REAL64), intent(out) :: df(n)

      integer :: i

      df(1) = (-1.0d0*f(3) + 4.0d0*f(2) - 3.0d0*f(1))/2.0d0
      df(2) = (-f(4) + 6.0d0*f(3) - 3.0d0*f(2) - 2.0d0*f(1))/6.0d0
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do
      df(n - 1) = (2.0d0*f(n) + 3.0d0*f(n - 1) - 6.0d0*f(n - 2) + f(n - 3))/6.0d0
      df(n) = (3.0d0*f(n) - 4.0d0*f(n - 1) + 1.0d0*f(n - 2))/2.0d0

   end subroutine dx_fdm4nd_center_normal_r8

   pure subroutine dx_fdm4nd_center_head_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(1:n - 2)

      integer :: i

      df(1) = (-1.0d0*f(3) + 4.0d0*f(2) - 3.0d0*f(1))/2.0d0
      df(2) = (-f(4) + 6.0d0*f(3) - 3.0d0*f(2) - 2.0d0*f(1))/6.0d0
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do

   end subroutine dx_fdm4nd_center_head_r8

   pure subroutine dx_fdm4nd_center_tail_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(3:n)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do
      df(n - 1) = (2.0d0*f(n) + 3.0d0*f(n - 1) - 6.0d0*f(n - 2) + f(n - 3))/6.0d0
      df(n) = (3.0d0*f(n) - 4.0d0*f(n - 1) + 1.0d0*f(n - 2))/2.0d0

   end subroutine dx_fdm4nd_center_tail_r8

   pure subroutine dx_fdm4nd_center_body_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(3:n - 2)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do

   end subroutine dx_fdm4nd_center_body_r8

   pure subroutine dx_fdm4nd_center_z8(f, ib1, ie1, df, ib2, ie2)

      implicit none

      complex(REAL64), intent(in) :: f(:)
      complex(REAL64), intent(out) :: df(:)
      integer, intent(in) :: ib1, ie1, ib2, ie2

      integer :: n1, n2
      integer :: n

      n1 = ie1 - ib1 + 1; n2 = ie2 - ib2 + 1

      if (n1 == n2) then
         call dx_fdm4nd_center_normal(f, df, n1)
      elseif (n1 /= n2) then
         if ((ib1 == ib2) .and. (ie1 /= ie2)) then
            call dx_fdm4nd_center_head(f, df, n1)
         elseif ((ie1 == ie2) .and. (ib1 /= ib2)) then
            call dx_fdm4nd_center_tail(f, df, n1)
         elseif ((ib1 /= ib2) .and. (ie1 /= ie2) .and. ((ib2 - ib1) == (ie1 - ie2))) then
            call dx_fdm4nd_center_body(f, df, n1)
         ! else
         !    print *, 'there is some mistake with dx_fdm_4nd'
         end if
      end if

   end subroutine dx_fdm4nd_center_z8

   pure subroutine dx_fdm4nd_center_normal_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(n)
      complex(REAL64), intent(out) :: df(n)

      integer :: i

      df(1) = (-1.0d0*f(3) + 4.0d0*f(2) - 3.0d0*f(1))/2.0d0
      df(2) = (-f(4) + 6.0d0*f(3) - 3.0d0*f(2) - 2.0d0*f(1))/6.0d0
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do
      df(n - 1) = (2.0d0*f(n) + 3.0d0*f(n - 1) - 6.0d0*f(n - 2) + f(n - 3))/6.0d0
      df(n) = (3.0d0*f(n) - 4.0d0*f(n - 1) + 1.0d0*f(n - 2))/2.0d0

   end subroutine dx_fdm4nd_center_normal_z8

   pure subroutine dx_fdm4nd_center_head_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(1:n - 2)

      integer :: i

      df(1) = (-1.0d0*f(3) + 4.0d0*f(2) - 3.0d0*f(1))/2.0d0
      df(2) = (-f(4) + 6.0d0*f(3) - 3.0d0*f(2) - 2.0d0*f(1))/6.0d0
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do

   end subroutine dx_fdm4nd_center_head_z8

   pure subroutine dx_fdm4nd_center_tail_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(3:n)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do
      df(n - 1) = (2.0d0*f(n) + 3.0d0*f(n - 1) - 6.0d0*f(n - 2) + f(n - 3))/6.0d0
      df(n) = (3.0d0*f(n) - 4.0d0*f(n - 1) + 1.0d0*f(n - 2))/2.0d0

   end subroutine dx_fdm4nd_center_tail_z8

   pure subroutine dx_fdm4nd_center_body_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(3:n - 2)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 8.0d0*f(i + 1) - 8.0d0*f(i - 1) + f(i - 2))/12.0d0
      end do

   end subroutine dx_fdm4nd_center_body_z8

   pure subroutine dxx_fdm4nd_center_r8(f, ib1, ie1, df, ib2, ie2)

      implicit none

      real(REAL64), intent(in) :: f(ib1:ie1)
      real(REAL64), intent(out) :: df(ib2:ie2)
      integer, intent(in) :: ib1, ie1, ib2, ie2

      integer :: n1, n2
      integer :: n

      n1 = ie1 - ib1 + 1; n2 = ie2 - ib2 + 1

      if (n1 == n2) then
         call dxx_fdm4nd_center_normal(f, df, n1)
      elseif (n1 /= n2) then
         if ((ib1 == ib2) .and. (ie1 /= ie2)) then
            call dxx_fdm4nd_center_head(f, df, n1)
         elseif ((ie1 == ie2) .and. (ib1 /= ib2)) then
            call dxx_fdm4nd_center_tail(f, df, n1)
         elseif ((ib1 /= ib2) .and. (ie1 /= ie2) .and. ((ib2 - ib1) == (ie1 - ie2))) then
            call dxx_fdm4nd_center_body(f, df, n1)
         ! else
         !    print *, 'There is some mistake with dx_fdm_4nd'
         end if
      end if

   end subroutine dxx_fdm4nd_center_r8

   pure subroutine dxx_fdm4nd_center_z8(f, ib1, ie1, df, ib2, ie2)

      implicit none

      complex(REAL64), intent(in) :: f(ib1:ie1)
      complex(REAL64), intent(out) :: df(ib2:ie2)
      integer, intent(in) :: ib1, ie1, ib2, ie2

      integer :: n1, n2
      integer :: n

      n1 = ie1 - ib1 + 1; n2 = ie2 - ib2 + 1

      if (n1 == n2) then
         call dxx_fdm4nd_center_normal(f, df, n1)
      elseif (n1 /= n2) then
         if ((ib1 == ib2) .and. (ie1 /= ie2)) then
            call dxx_fdm4nd_center_head(f, df, n1)
         elseif ((ie1 == ie2) .and. (ib1 /= ib2)) then
            call dxx_fdm4nd_center_tail(f, df, n1)
         elseif ((ib1 /= ib2) .and. (ie1 /= ie2) .and. ((ib2 - ib1) == (ie1 - ie2))) then
            call dxx_fdm4nd_center_body(f, df, n1)
         ! else
            ! print *, 'There is some mistake with dx_fdm_4nd'
         end if
      end if

   end subroutine dxx_fdm4nd_center_z8

   pure subroutine dxx_fdm4nd_center_normal_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(n)
      real(REAL64), intent(out) :: df(n)

      ! real(REAL64), parameter :: d1(7)=[1.0d0/90.0d0,-3.0d0/20.0d0,3.0d0/2.0d0,-49.0d0/18.0d0,3.0d0/2.0d0,-3.0d0/20.0d0,1.0d0/90.0d0]

      integer :: i

      df(1) = f(1) - 2.0d0*f(2) + f(3)
      df(2) = f(1) - 2.0d0*f(2) + f(3)
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
         ! df(i)=sum(d1*f(i-3:i+3))
      end do
      df(n - 1) = f(n) - 2.0d0*f(n - 1) + f(n - 2)
      df(n) =  f(n) - 2.0d0*f(n - 1) + f(n - 2)

   end subroutine dxx_fdm4nd_center_normal_r8

   pure subroutine dxx_fdm4nd_center_head_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(1:n - 2)

      integer :: i

      df(1) = f(1) - 2.0d0*f(2) + f(3)
      df(2) = f(1) - 2.0d0*f(2) + f(3)
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
      end do

   end subroutine dxx_fdm4nd_center_head_r8

   pure subroutine dxx_fdm4nd_center_tail_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(3:n)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
      end do
      df(n - 1) = f(n) - 2.0d0*f(n - 1) + f(n - 2)
      df(n) = f(n) - 2.0d0*f(n - 1) + f(n - 2)

   end subroutine dxx_fdm4nd_center_tail_r8

   pure subroutine dxx_fdm4nd_center_body_r8(f, df, n)

      implicit none
      integer, intent(in) :: n
      real(REAL64), intent(in) :: f(1:n)
      real(REAL64), intent(out) :: df(3:n - 2)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
         ! df(i)=sum(dd1(:, 0)*f(i-2:i+2))
      end do

   end subroutine dxx_fdm4nd_center_body_r8

   pure subroutine dxx_fdm4nd_center_normal_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(n)
      complex(REAL64), intent(out) :: df(n)

      integer :: i

      df(1) = f(1) - 2.0d0*f(2) + f(3)
      df(2) = f(1) - 2.0d0*f(2) + f(3)
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
         ! df(i)=sum(d1*f(i-3:i+3))
      end do
      df(n - 1) = f(n) - 2.0d0*f(n - 1) + f(n - 2)
      df(n) = f(n) - 2.0d0*f(n - 1) + f(n - 2)

   end subroutine dxx_fdm4nd_center_normal_z8

   pure subroutine dxx_fdm4nd_center_head_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(1:n - 2)

      integer :: i

      df(1) = f(1) - 2.0d0*f(2) + f(3)
      df(2) = f(1) - 2.0d0*f(2) + f(3)
      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
      end do

   end subroutine dxx_fdm4nd_center_head_z8

   pure subroutine dxx_fdm4nd_center_tail_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(3:n)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
      end do
      df(n - 1) = f(n) - 2.0d0*f(n - 1) + f(n - 2)
      df(n) = f(n) - 2.0d0*f(n - 1) + f(n - 2)

   end subroutine dxx_fdm4nd_center_tail_z8

   pure subroutine dxx_fdm4nd_center_body_z8(f, df, n)

      implicit none
      integer, intent(in) :: n
      complex(REAL64), intent(in) :: f(1:n)
      complex(REAL64), intent(out) :: df(3:n - 2)

      integer :: i

      do i = 3, n - 2
         df(i) = (-f(i + 2) + 16.0d0*f(i + 1) - 30.0d0*f(i) + 16.0d0*f(i - 1) - f(i - 2))/12.0d0
         ! df(i)=sum(dd1(:, 0)*f(i-2:i+2))
      end do

   end subroutine dxx_fdm4nd_center_body_z8

end module mod_dx
