#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>
 
module mod_flow_grad

    use global_parameters
    implicit none

    Vec :: dq0dI
    Vec :: dq0dI_local
    Vec :: dq0dJ
    Vec :: dq0dJ_local

    integer :: ib, ie, jb, je
    integer :: ibg, ieg, jbg, jeg

    private
    public :: flow_grad

contains

    subroutine flow_grad(comm)

        implicit none
        PetscInt, INTENT(in) :: comm

        ib = xb + 1; ie = xe + 1; jb = yb + 1; je = ye + 1
        ibg = xbg + 1; ieg = xeg + 1; jbg = ybg + 1; jeg = yeg + 1

        call allocate_flow_grad()

        call compute_convariant_flow_grad()

        call deallocate_local_flow_grad()

        call PetscPrintf(comm, "Compute flow grads done!\n", ierr)

    end subroutine flow_grad

    subroutine allocate_flow_grad()

        implicit none

        call DMGetGlobalVector(meshDA, dq0dI, ierr)
        call VecDuplicate(dq0dI, dq0dJ, ierr)

        call DMGetLocalVector(meshDA, dq0dI_local, ierr)
        call VecDuplicate(dq0dI_local, dq0dJ_local, ierr)

        allocate (dq0dx(5, ibg:ieg, jbg:jeg))
        allocate (dq0dy(5, ibg:ieg, jbg:jeg))
        allocate (dq0dxx(5, ibg:ieg, jbg:jeg))
        allocate (dq0dyy(5, ibg:ieg, jbg:jeg))
        allocate (dq0dxy(5, ibg:ieg, jbg:jeg))

    end subroutine allocate_flow_grad

    subroutine compute_convariant_flow_grad()
!!@TODO_the second ddx should be rewrite for that value is wrong.
!! This is because that the method of ddx is different between CFD and START.
!! It cause a wrong value of ddy in START, but, it seems that the eigenvalues is
!! still right.

        use mod_dx, only: dx_fdm4nd_center, dxx_fdm4nd_center
        implicit none
        real(REAL64), allocatable, dimension(:, :, :) :: dq0_i, dq0_ii, dq0_j, dq0_jj, dq0_ij
        real(REAL64), allocatable, dimension(:, :, :) :: dq0_x_local, dq0_y_local
        ! real(REAL64), allocatable, dimension(:, :, :) :: dq0_ij_debug
        PetscScalar, pointer :: tmp(:, :, :)

        integer :: i, j, l

        allocate (dq0_i(5, ib:ie, jb:je))
        allocate (dq0_ii(5, ib:ie, jb:je))
        allocate (dq0_j(5, ib:ie, jb:je))
        allocate (dq0_jj(5, ib:ie, jb:je))
        allocate (dq0_ij(5, ib:ie, jb:je))
        ! allocate(dq0_ij_debug(5, ib:ie, jb:je))

        ! !> dq0di, dq0dii
        ! do j=jb, je
        !   do l=1, 5
        !     call dx_fdm4nd_center (qq(l, :, j), ibg, ieg, dq0_i (l, :, j), ib, ie)
        !     ! call dxx_fdm4nd_center(qq(l, :, j), ibg, ieg, dq0_ii(l, :, j), ib, ie)
        !   enddo
        ! enddo
        !
        ! !> dq0dj, dq0djj
        ! do i=ib, ie
        !   do l=1, 5
        !     call dx_fdm4nd_center (qq(l, i, :), jbg, jeg, dq0_j(l, i, :), jb, je)
        !     ! call dxx_fdm4nd_center(qq(l, i, :), jbg, jeg, dq0_jj(l, i, :), jb, je)
        !   enddo
        ! enddo
        !
        ! !! 1nd grad
        ! do l=1, 5
        !   call convariant_grad_1nd(dq0dx(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), xi_x, eta_x)
        !   call convariant_grad_1nd(dq0dy(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), xi_y, eta_y)
        ! enddo

        !> 1nd grad

        do l = 1, 5
            call dfdx(qq(l, :, :), dq0dx(l, ib:ie, jb:je))
            call dfdy(qq(l, :, :), dq0dy(l, ib:ie, jb:je))
        end do

        !> send and receive dq0dx, dq0dy

        call DMDAVecGetArrayF90(meshDA, dq0dI, tmp, ierr)
        tmp(:, :, :) = dq0dx(:, ib:ie, jb:je)
        call DMDAVecRestoreArrayF90(meshDA, dq0dI, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMGlobalToLocalEnd(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)

        call DMDAVecGetArrayF90(meshDA, dq0dJ, tmp, ierr)
        tmp(:, :, :) = dq0dy(:, ib:ie, jb:je)
        call DMDAVecRestoreArrayF90(meshDA, dq0dJ, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, dq0dJ, INSERT_VALUES, dq0dJ_local, ierr)
        call DMGlobalToLocalEnd(meshDA, dq0dJ, INSERT_VALUES, dq0dJ_local, ierr)

        allocate (dq0_x_local(5, ibg:ieg, jbg:jeg))
        allocate (dq0_y_local(5, ibg:ieg, jbg:jeg))

        call DMDAVecGetArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        dq0_x_local(:, :, :) = tmp(:, :, :)
        dq0dx = tmp
        call DMDAVecRestoreArrayReadF90(meshDA, dq0dI_local, tmp, ierr)

        call DMDAVecGetArrayReadF90(meshDA, dq0dJ_local, tmp, ierr)
        dq0_y_local(:, :, :) = tmp(:, :, :)
        dq0dy = tmp
        call DMDAVecRestoreArrayReadF90(meshDA, dq0dJ_local, tmp, ierr)

    !!2nd grad
        do l = 1, 5
            call dfdx(dq0_x_local(l, :, :), dq0dxx(l, :, :))
            call dfdy(dq0_y_local(l, :, :), dq0dyy(l, :, :))
            call dfdx(dq0_y_local(l, :, :), dq0dxy(l, :, :))
        end do

        call DMDAVecGetArrayF90(meshDA, dq0dI, tmp, ierr)
        tmp(:, :, :) = dq0dxx(:, ib:ie, jb:je)
        call DMDAVecRestoreArrayF90(meshDA, dq0dI, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMGlobalToLocalEnd(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMDAVecGetArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        dq0dxx = tmp
        call DMDAVecRestoreArrayReadF90(meshDA, dq0dI_local, tmp, ierr)

        call DMDAVecGetArrayF90(meshDA, dq0dI, tmp, ierr)
        tmp(:, :, :) = dq0dyy(:, ib:ie, jb:je)
        call DMDAVecRestoreArrayF90(meshDA, dq0dI, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMGlobalToLocalEnd(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMDAVecGetArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        dq0dyy = tmp
        call DMDAVecRestoreArrayReadF90(meshDA, dq0dI_local, tmp, ierr)

        call DMDAVecGetArrayF90(meshDA, dq0dI, tmp, ierr)
        tmp(:, :, :) = dq0dxy(:, ib:ie, jb:je)
        call DMDAVecRestoreArrayF90(meshDA, dq0dI, tmp, ierr)
        call DMGlobalToLocalBegin(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMGlobalToLocalEnd(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        call DMDAVecGetArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        dq0dxy = tmp
        call DMDAVecRestoreArrayReadF90(meshDA, dq0dI_local, tmp, ierr)

        ! !> send and receive dq0di, dq0dj
        !
        ! ! call DMDAVecGetArrayF90(meshDA, dq0dI, tmp, ierr)
        ! ! tmp(:, :, :)=dq0_i(:, :, :)
        ! ! call DMDAVecRestoreArrayF90(meshDA, dq0dI, tmp, ierr)
        ! ! call DMGlobalToLocalBegin(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        !
        ! call DMDAVecGetArrayF90(meshDA, dq0dJ, tmp, ierr)
        ! tmp(:, :, :)=dq0_j(:, :, :)
        ! call DMDAVecRestoreArrayF90(meshDA, dq0dJ, tmp, ierr)
        ! ! call DMGlobalToLocalEnd(meshDA, dq0dI, INSERT_VALUES, dq0dI_local, ierr)
        ! call DMGlobalToLocalBegin(meshDA, dq0dJ, INSERT_VALUES, dq0dJ_local, ierr)
        ! call DMGlobalToLocalEnd(meshDA, dq0dJ, INSERT_VALUES, dq0dJ_local, ierr)
        !
        ! allocate(dq0_j_local(5, ibg:ieg, jbg:jeg))
        ! ! allocate(dq0_j_local(5, ibg:ieg, jbg:jeg))
        !
        ! ! call DMDAVecGetArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        ! ! dq0_i_local(:, :, :)=tmp(:, :, :)
        ! ! call DMDAVecRestoreArrayReadF90(meshDA, dq0dI_local, tmp, ierr)
        ! !
        ! call DMDAVecGetArrayReadF90(meshDA, dq0dJ_local, tmp, ierr)
        ! dq0_j_local(:, :, :)=tmp(:, :, :)
        ! call DMDAVecRestoreArrayReadF90(meshDA, dq0dJ_local, tmp, ierr)
        !
        ! do j=jb, je
        !   do l=1, 5
        !      call dx_fdm4nd_center(dq0_j_local(l, :, j), ibg, ieg, dq0_ij(l, :, j), ib, ie)
        !   enddo
        ! enddo
        !
        ! ! do i=ib, ie
        ! !   do l=1, 5
        ! !      call dx_fdm4nd_center(dq0_i_local(l, i, :), jbg, jeg, dq0_ij_debug(l, i, :), jb, je)
        ! !   enddo
        ! ! enddo
        !
        ! !> get_convariant_grad
        ! !! 1nd grad
        ! do l=1, 5
        !   call convariant_grad_1nd(dq0dx(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), xi_x, eta_x)
        !   call convariant_grad_1nd(dq0dy(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), xi_y, eta_y)
        !   call convariant_grad_2nd(dq0dxx(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), &
        !                            dq0_ii(l, :, :), dq0_jj(l, :, :), dq0_ij(l, :, :), &
        !                            xi_x, xi_x, eta_x, eta_x, xi_xx, eta_xx)
        !   call convariant_grad_2nd(dq0dyy(l, :, :), dq0_i(l, :, :), dq0_j(l, :, :), &
        !                            dq0_ii(l, :, :), dq0_jj(l, :, :), dq0_ij(l, :, :), &
        !                            xi_y, xi_y, eta_y, eta_y, xi_yy, eta_yy)
        ! enddo
        !
        ! call convariant_grad_2nd(du0xy(:, :), dq0_i(2, :, :), dq0_j(2, :, :), &
        !                          dq0_ii(2, :, :), dq0_jj(2, :, :), dq0_ij(2, :, :), &
        !                          xi_x, xi_y, eta_x, eta_y, xi_xy, eta_xy)
        !
        !
        ! call convariant_grad_2nd(dv0xy(:, :), dq0_i(3, :, :), dq0_j(3, :, :), &
        !                          dq0_ii(3, :, :), dq0_jj(3, :, :), dq0_ij(3, :, :), &
        !                          xi_x, xi_y, eta_x, eta_y, xi_xy, eta_xy)

        ! block
        !   use global_parameters, only: cc1=>FDM_C1_4ORD_CENTER, dd1=>FDM_D1_4ORD_CENTER
        !   real(REAL64), parameter :: ee1(3)=[1.0d0, -2.0d0,+1.0d0]
        !     print*, dq0dy(2, 2, 37), dq0dyy(2, 2, 37)
        !     print('(5E20.7)'), sum(cc1(:,0)*dq0_j(2, 2, 34:38)) ,sum(cc1(:,0)*dq0_j(2, 2, 35:39)), sum(cc1(:,0)*dq0_j(2, 2, 36:40))
        !     print('(5E20.7)'), sum(dd1(:,0)*qq(2, 2, 34:38)), sum(dd1(:,0)*qq(2, 2, 35:39)), sum(dd1(:,0)*qq(2, 2, 36:40))
        !     print('(5E20.7)'), sum(ee1(:)*qq(2, 2, 35:37)), sum(ee1(:)*qq(2, 2, 36:38)), sum(ee1(:)*qq(2, 2, 37:39))
        !
        !     print*, 'coffs'
        !     print('(5E20.7)'), dq0_jj(2, 2, 35:39)
        !     print('(5E20.7)'), dq0_j(2, 2, 35:39)
        !     print('(5E20.7)'), qq(2, 2, 35:39)
        !     print('(5E20.7)'), eta_y(2, 35:39)
        !     print('(5E20.7)'), eta_yy(2, 35:39)
        !     print('(5E20.7)'),sum(cc1(:,0)*dq0_j(2, 2, 33:37)),sum(cc1(:,0)*dq0_j(2, 2, 34:38)), &
        !     sum(cc1(:,0)*dq0_j(2, 2, 35:39)), sum(cc1(:,0)*dq0_j(2, 2, 36:40)), sum(cc1(:,0)*dq0_j(2, 2, 37:41))
        !     print('(5E20.7)'), dq0_jj(2, 2, 35:39)
        !
        !     ! print('(4E20.7)'), xi_y(2, 37), eta_y(2, 37), xi_yy(2, 37), eta_yy(2, 37)
        !     ! print*, eta_yy(2, 37)*dq0_j(2, 2, 37)+eta_y(2, 37)**2*dq0_jj(2, 2, 37)
        !     ! print*, dq0dxx(2, 2, 37), dq0_ii(2, 2, 37), dq0_jj(2, 2, 37)
        !     ! stop
        !
        ! endblock

        deallocate (dq0_i)
        deallocate (dq0_ii)
        deallocate (dq0_j)
        deallocate (dq0_jj)
        deallocate (dq0_ij)

    end subroutine compute_convariant_flow_grad

    subroutine deallocate_local_flow_grad()

        implicit none

        call VecDestroy(dq0dJ, ierr)
        call VecDestroy(dq0dJ_local, ierr)

    end subroutine deallocate_local_flow_grad

    elemental subroutine convariant_grad_1nd(dfdx, dfdi, dfdj, didx, djdx)

        implicit none

        real(REAL64), intent(out) :: dfdx
        real(REAL64), intent(in)  :: dfdi, dfdj, didx, djdx

        dfdx = dfdi*didx + dfdj*djdx

    end subroutine convariant_grad_1nd

    elemental subroutine convariant_grad_2nd(dfdxy, dfdi, dfdj, dfdii, dfdjj, dfdij, &
                                             didx, didy, djdx, djdy, didxy, djdxy)

        implicit none

        real(REAL64), intent(out) :: dfdxy
        real(REAL64), intent(in)  :: dfdi, dfdj, dfdii, dfdjj, dfdij
        real(REAL64), intent(in)  :: didx, didy, djdx, djdy, didxy, djdxy

        dfdxy = dfdi*didxy + dfdj*djdxy + dfdii*didx*didy + dfdjj*djdx*djdy + &
                dfdij*(didx*djdy + didy*djdx)

    end subroutine convariant_grad_2nd

    subroutine dfdx(f, df)

        use mod_dx, only: dx_fdm4nd_center
        implicit none
        real(REAL64), intent(inout) :: f(:, :)
        real(REAL64), intent(inout) :: df(:, :)
        real(REAL64) :: df_i(ib:ie, jb:je), df_j(ib:ie, jb:je)
        integer :: i, j

        !dfdi
        do j = jb, je
            call dx_fdm4nd_center(f(:, j), ibg, ieg, df_i(:, j), ib, ie)
        end do

        !> dfdj
        do i = ib, ie
            call dx_fdm4nd_center(f(i, :), jbg, jeg, df_j(i, :), jb, je)
        end do

        do j = jb, je
            do i = ib, ie
                call convariant_grad_1nd(df(i, j), df_i(i, j), df_j(i, j), xi_x(i, j), eta_x(i, j))
            end do
        end do

    end subroutine dfdx

    subroutine dfdy(f, df)

        use mod_dx, only: dx_fdm4nd_center
        implicit none

        real(REAL64), intent(inout) :: f(:, :)
        real(REAL64), intent(inout) :: df(:, :)
        real(REAL64) :: df_i(ib:ie, jb:je), df_j(ib:ie, jb:je)
        integer :: i, j

        !dfdi
        do j = jb, je
            call dx_fdm4nd_center(f(:, j), ibg, ieg, df_i(:, j), ib, ie)
        end do

        !> dfdj
        do i = ib, ie
            call dx_fdm4nd_center(f(i, :), jbg, jeg, df_j(i, :), jb, je)
        end do

        do j = jb, je
            do i = ib, ie
                call convariant_grad_1nd(df(i, j), df_i(i, j), df_j(i, j), xi_y(i, j), eta_y(i, j))
            end do
        end do

    end subroutine dfdy

end module mod_flow_grad
