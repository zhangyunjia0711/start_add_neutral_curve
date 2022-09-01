!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code FooPSE
!------------------------------------------------------------------------------
!
!  module: mod_grid_profile
!> @breif modules for set the grid profile.
!  DESCRIPTION:
!> In LST, the y_grid in normal direction must be set and used in lst_eqn.
!! It contains the y distribution, the metrics, Lame coef and etc.
!!
!  REVISION HISTORY:
!  2021-01-05 - Initial Version
!> @author
!> Liu Jianxin
!> \date 2021-01-05
!------------------------------------------------------------------------------

module mod_lst_grid_profile

   use global_parameters
   implicit none

   type :: type_lst_grid_profile

      logical :: has_assemble = .FALSE.
      integer :: jb, je

      real(REAL64) :: drx = 0.0d0, drz = 0.0d0

      real(REAL64), allocatable, dimension(:) :: y
      real(REAL64), allocatable, dimension(:) :: xi_x
      real(REAL64), allocatable, dimension(:) :: xi_y
      real(REAL64), allocatable, dimension(:) :: eta_x
      real(REAL64), allocatable, dimension(:) :: eta_y
      real(REAL64), allocatable, dimension(:) :: xi_xx
      real(REAL64), allocatable, dimension(:) :: xi_xy
      real(REAL64), allocatable, dimension(:) :: xi_yy
      real(REAL64), allocatable, dimension(:) :: eta_xx
      real(REAL64), allocatable, dimension(:) :: eta_xy
      real(REAL64), allocatable, dimension(:) :: eta_yy

   contains

      procedure :: print => print_lst_grid_profile
      procedure :: finalize => finalize_lst_grid_profile

   end type type_lst_grid_profile

   private
   public :: type_lst_grid_profile
   public :: assemble_lst_grid_profile

contains

   function assemble_lst_grid_profile(ix) result(profile)

      implicit none
      integer, intent(in) :: ix
      type(type_lst_grid_profile) :: profile
      integer :: jb, je

      jb = yb + 1; je = ye + 1

      if (.not. profile%has_assemble) then
         allocate (profile%y(jb:je))
         allocate (profile%xi_x(jb:je))
         allocate (profile%xi_y(jb:je))
         allocate (profile%eta_x(jb:je))
         allocate (profile%eta_y(jb:je))
         allocate (profile%xi_xx(jb:je))
         allocate (profile%xi_xy(jb:je))
         allocate (profile%xi_yy(jb:je))
         allocate (profile%eta_xx(jb:je))
         allocate (profile%eta_xy(jb:je))
         allocate (profile%eta_yy(jb:je))
      end if

      profile%has_assemble = .TRUE.
      profile%jb = yb + 1; profile%je = ye + 1

      profile%drx = drx(ix); profile%drz = drz(ix)

      profile%y = yy(ix, jb:je)
      profile%xi_x = xi_x(ix, :)
      profile%xi_y = xi_y(ix, :)
      profile%eta_x = eta_x(ix, :)
      profile%eta_y = eta_y(ix, :)
      profile%xi_xx = xi_xx(ix, :)
      profile%xi_xy = xi_xy(ix, :)
      profile%xi_yy = xi_yy(ix, :)
      profile%eta_xx = eta_xx(ix, :)
      profile%eta_xy = eta_xy(ix, :)
      profile%eta_yy = eta_yy(ix, :)

   end function assemble_lst_grid_profile

   subroutine print_lst_grid_profile(profile, filename)

      implicit none
      class(type_lst_grid_profile), intent(in) :: profile
      character(*), intent(in) :: filename
      character(2) :: fn
      integer :: j

      write (fn, '(I2.2)') rank
      open (99, file=filename//fn)
      do j = profile%jb, profile%je
         write (99, '(2X, I5, 10(2X, E16.6))') j, profile%y(j), profile%eta_y(j), profile%eta_yy(j)
      end do
      close (99)

   end subroutine print_lst_grid_profile

   subroutine finalize_lst_grid_profile(profile)
      implicit none
      class(type_lst_grid_profile) :: profile

      profile%has_assemble = .FALSE.
      profile%jb = 0; profile%je = 0
      profile%drx = 0.0d0; profile%drz = 0.0d0

      deallocate (profile%y)
      deallocate (profile%xi_x)
      deallocate (profile%xi_y)
      deallocate (profile%eta_x)
      deallocate (profile%eta_y)
      deallocate (profile%xi_xx)
      deallocate (profile%xi_xy)
      deallocate (profile%xi_yy)
      deallocate (profile%eta_xx)
      deallocate (profile%eta_xy)
      deallocate (profile%eta_yy)

   end subroutine finalize_lst_grid_profile

end module mod_lst_grid_profile
