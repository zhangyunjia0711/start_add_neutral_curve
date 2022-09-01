#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_metrics

   use global_parameters

   implicit none

   private
   public :: metrics

   Vec :: XIX, XIY, ETAX, ETAY
   Vec :: XIX_local, XIY_local, ETAX_local, ETAY_local

   real(REAL64), dimension(:, :), allocatable :: x_xi, y_xi, x_eta, y_eta
   real(REAL64), dimension(:, :), allocatable :: jacobi

   integer :: ib, ie, jb, je
   integer :: ibg, ieg, jbg, jeg

contains

   subroutine metrics(comm)

      implicit none
      PetscInt, INTENT(in) :: comm

      ib = xb + 1; ie = xe + 1; jb = yb + 1; je = ye + 1
      ibg = xbg + 1; ieg = xeg + 1; jbg = ybg + 1; jeg = yeg + 1

      call allocate_metrics()

      call compute_contravariant_metrics()

      call compute_convariant_metrics()

      ! call print_metrics()

      call deallocate_local_metrics()

      call PetscPrintf(comm, "Compute metrics!\n", ierr)

   end subroutine metrics

   subroutine allocate_metrics()

      implicit none

      call DMGetGlobalVector(DA, XIX, ierr)
      call DMGetLocalVector(DA, XIX_local, ierr)

      call VecDuplicate(XIX, XIY, ierr)
      call VecDuplicate(XIX, ETAX, ierr)
      call VecDuplicate(XIX, ETAY, ierr)

      call VecDuplicate(XIX_local, XIY_local, ierr)
      call VecDuplicate(XIX_local, ETAX_local, ierr)
      call VecDuplicate(XIX_local, ETAY_local, ierr)

      allocate (xi_x(ib:ie, jb:je))
      allocate (xi_y(ib:ie, jb:je))
      allocate (eta_x(ib:ie, jb:je))
      allocate (eta_y(ib:ie, jb:je))
      allocate (xi_xx(ib:ie, jb:je))
      allocate (xi_yy(ib:ie, jb:je))
      allocate (xi_xy(ib:ie, jb:je))
      allocate (eta_xx(ib:ie, jb:je))
      allocate (eta_yy(ib:ie, jb:je))
      allocate (eta_xy(ib:ie, jb:je))

      allocate (x_xi(ib:ie, jb:je))
      allocate (y_xi(ib:ie, jb:je))
      allocate (x_eta(ib:ie, jb:je))
      allocate (y_eta(ib:ie, jb:je))

      allocate (jacobi(ib:ie, jb:je))

   end subroutine allocate_metrics

   subroutine compute_contravariant_metrics()

      use mod_dx, only: dx_fdm4nd_center
      implicit none

      integer :: i, j

      do j = jb, je
         call dx_fdm4nd_center(xx(:, j), ibg, ieg, x_xi(:, j), ib, ie)
         call dx_fdm4nd_center(yy(:, j), ibg, ieg, y_xi(:, j), ib, ie)
      end do
      do i = ib, ie
         call dx_fdm4nd_center(xx(i, :), jbg, jeg, x_eta(i, :), jb, je)
         call dx_fdm4nd_center(yy(i, :), jbg, jeg, y_eta(i, :), jb, je)
      end do

   end subroutine compute_contravariant_metrics

   subroutine compute_convariant_metrics()

      use mod_dx, only: dx_fdm4nd_center
      implicit none
      integer :: i, j
      PetscScalar, pointer :: tmp(:, :)
      real(REAL64), allocatable, dimension(:, :) :: xi_x_local
      real(REAL64), allocatable, dimension(:, :) :: xi_y_local
      real(REAL64), allocatable, dimension(:, :) :: eta_x_local
      real(REAL64), allocatable, dimension(:, :) :: eta_y_local
      real(REAL64), allocatable, dimension(:, :) :: xi_x_xi
      real(REAL64), allocatable, dimension(:, :) :: xi_x_eta
      real(REAL64), allocatable, dimension(:, :) :: xi_y_xi
      real(REAL64), allocatable, dimension(:, :) :: xi_y_eta
      real(REAL64), allocatable, dimension(:, :) :: eta_x_xi
      real(REAL64), allocatable, dimension(:, :) :: eta_x_eta
      real(REAL64), allocatable, dimension(:, :) :: eta_y_xi
      real(REAL64), allocatable, dimension(:, :) :: eta_y_eta

      do j = jb, je
         do i = ib, ie
            jacobi(i, j) = 1.0d0/(x_xi(i, j)*y_eta(i, j) - y_xi(i, j)*x_eta(i, j))
         end do
      end do

    !! 1nd convariant metrics
      do j = jb, je
         do i = ib, ie
            xi_x(i, j) = y_eta(i, j)*jacobi(i, j)
            xi_y(i, j) = -x_eta(i, j)*jacobi(i, j)
            eta_x(i, j) = -y_xi(i, j)*jacobi(i, j)
            eta_y(i, j) = x_xi(i, j)*jacobi(i, j)
         end do
      end do

    !!!!!!!!!!!!!!!!!!!!!!!!!
    !! 2nd convariant metrics
    !!!!!!!!!!!!!!!!!!!!!!!!!

    !! send and receive
      call DMDAVecGetArrayF90(DA, XIX, tmp, ierr)
      tmp(:, :) = xi_x(:, :)
      call DMDAVecRestoreArrayF90(DA, XIX, tmp, ierr)
      call DMGlobalToLocalBegin(DA, XIX, INSERT_VALUES, XIX_local, ierr)
      call DMGlobalToLocalEnd(DA, XIX, INSERT_VALUES, XIX_local, ierr)

      call DMDAVecGetArrayF90(DA, XIY, tmp, ierr)
      tmp(:, :) = xi_y(:, :)
      call DMDAVecRestoreArrayF90(DA, XIY, tmp, ierr)
      call DMGlobalToLocalBegin(DA, XIY, INSERT_VALUES, XIY_local, ierr)
      call DMGlobalToLocalEnd(DA, XIY, INSERT_VALUES, XIY_local, ierr)

      call DMDAVecGetArrayF90(DA, ETAX, tmp, ierr)
      tmp(:, :) = eta_x(:, :)
      call DMDAVecRestoreArrayF90(DA, ETAX, tmp, ierr)
      call DMGlobalToLocalBegin(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)
      call DMGlobalToLocalEnd(DA, ETAX, INSERT_VALUES, ETAX_local, ierr)

      call DMDAVecGetArrayF90(DA, ETAY, tmp, ierr)
      tmp(:, :) = eta_y(:, :)
      call DMDAVecRestoreArrayF90(DA, ETAY, tmp, ierr)
      call DMGlobalToLocalBegin(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)
      call DMGlobalToLocalEND(DA, ETAY, INSERT_VALUES, ETAY_local, ierr)

      allocate (xi_x_local(ibg:ieg, jbg:jeg))
      allocate (xi_y_local(ibg:ieg, jbg:jeg))
      allocate (eta_x_local(ibg:ieg, jbg:jeg))
      allocate (eta_y_local(ibg:ieg, jbg:jeg))

      call DMDAVecGetArrayReadF90(DA, XIX_local, tmp, ierr)
      xi_x_local = tmp
      call DMDAVecRestoreArrayReadF90(DA, XIX_local, tmp, ierr)

      call DMDAVecGetArrayReadF90(DA, XIY_local, tmp, ierr)
      xi_y_local = tmp
      call DMDAVecRestoreArrayReadF90(DA, XIY_local, tmp, ierr)

      call DMDAVecGetArrayReadF90(DA, ETAX_local, tmp, ierr)
      eta_x_local = tmp
      call DMDAVecRestoreArrayReadF90(DA, ETAX_local, tmp, ierr)
      call DMDAVecGetArrayReadF90(DA, ETAY_local, tmp, ierr)
      eta_y_local = tmp
      call DMDAVecRestoreArrayReadF90(DA, ETAY_local, tmp, ierr)

      allocate (xi_x_xi(ib:ie, jb:je))
      allocate (xi_x_eta(ib:ie, jb:je))
      allocate (xi_y_xi(ib:ie, jb:je))
      allocate (xi_y_eta(ib:ie, jb:je))
      allocate (eta_x_xi(ib:ie, jb:je))
      allocate (eta_x_eta(ib:ie, jb:je))
      allocate (eta_y_xi(ib:ie, jb:je))
      allocate (eta_y_eta(ib:ie, jb:je))

      do j = jb, je
         call dx_fdm4nd_center(xi_x_local(:, j), ibg, ieg, xi_x_xi(:, j), ib, ie)
         call dx_fdm4nd_center(xi_y_local(:, j), ibg, ieg, xi_y_xi(:, j), ib, ie)
         call dx_fdm4nd_center(eta_x_local(:, j), ibg, ieg, eta_x_xi(:, j), ib, ie)
         call dx_fdm4nd_center(eta_y_local(:, j), ibg, ieg, eta_y_xi(:, j), ib, ie)
      end do
      do i = ib, ie
         call dx_fdm4nd_center(xi_x_local(i, :), jbg, jeg, xi_x_eta(i, :), jb, je)
         call dx_fdm4nd_center(xi_y_local(i, :), jbg, jeg, xi_y_eta(i, :), jb, je)
         call dx_fdm4nd_center(eta_x_local(i, :), jbg, jeg, eta_x_eta(i, :), jb, je)
         call dx_fdm4nd_center(eta_y_local(i, :), jbg, jeg, eta_y_eta(i, :), jb, je)
      end do

      do j = jb, je
         do i = ib, ie
            xi_xx(i, j) = xi_x(i, j)*xi_x_xi(i, j) + eta_x(i, j)*xi_x_eta(i, j)
            xi_yy(i, j) = xi_y(i, j)*xi_y_xi(i, j) + eta_y(i, j)*xi_y_eta(i, j)
            xi_xy(i, j) = xi_y(i, j)*xi_x_xi(i, j) + eta_y(i, j)*xi_x_eta(i, j)
            eta_xx(i, j) = xi_x(i, j)*eta_x_xi(i, j) + eta_x(i, j)*eta_x_eta(i, j)
            eta_yy(i, j) = xi_y(i, j)*eta_y_xi(i, j) + eta_y(i, j)*eta_y_eta(i, j)
            eta_xy(i, j) = xi_y(i, j)*eta_x_xi(i, j) + eta_y(i, j)*eta_x_eta(i, j)
         end do
      end do

      deallocate (xi_x_xi)
      deallocate (xi_x_eta)
      deallocate (xi_y_xi)
      deallocate (xi_y_eta)
      deallocate (eta_x_xi)
      deallocate (eta_x_eta)
      deallocate (eta_y_xi)
      deallocate (eta_y_eta)
      deallocate (xi_x_local)
      deallocate (xi_y_local)
      deallocate (eta_x_local)
      deallocate (eta_y_local)

   end subroutine compute_convariant_metrics

   subroutine deallocate_local_metrics()

      implicit none

      deallocate (x_xi)
      deallocate (y_xi)
      deallocate (x_eta)
      deallocate (y_eta)
      deallocate (jacobi)

      call VecDestroy(XIX, ierr)
      call VecDestroy(XIY, ierr)
      call VecDestroy(ETAX, ierr)
      call VecDestroy(ETAY, ierr)
      call VecDestroy(XIX_local, ierr)
      call VecDestroy(XIY_local, ierr)
      call VecDestroy(ETAX_local, ierr)
      call VecDestroy(ETAY_local, ierr)

   end subroutine deallocate_local_metrics

   subroutine print_metrics()

      implicit none

      if (rank == 1) then
         print *, rank
         print *, xi_x(xb + 250, yb + 20), xi_y(xb + 250, yb + 20)
         print *, eta_x(xb + 250, yb + 20), eta_y(xb + 250, yb + 20)
      end if

   end subroutine print_metrics

end module mod_metrics
