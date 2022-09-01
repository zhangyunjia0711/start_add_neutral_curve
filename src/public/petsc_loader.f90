#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

module mod_petsc_loader

   use petsc
   use global_parameters
   implicit none

   public :: load_info
   ! private

   PetscViewer :: Viewer
   PetscScalar, pointer :: flow(:, :, :), grid(:, :, :)

contains

   subroutine load_info(comm)

      implicit none

      PetscInt, INTENT(in) :: comm


      call read_curvature()
      call read_mesh_2d(comm)
      call load_mesh_info(comm)
      call load_flow_info(comm)
      call PetscPrintf(comm, "Load done!\n", ierr)

   end subroutine load_info

   subroutine plot3d_load()

      implicit none

      integer :: i, j, l

      print *, "Plot3D format is selected and the file format is translating..."

      open (11, file=char(gridfile), form='unformatted')
      read (11)
      read (11) in, jn
      print *, 'in=', in, 'jn=', jn
      allocate (xx_0(in, jn), yy_0(in, jn))
      read (11) xx_0, yy_0
      close (11)

      print *, 'read grid is done!'

      open (12, file=char(flowfile), form='unformatted')
      read (12)
      read (12) in, jn, ln
      print *, 'in=', in, 'jn=', jn, 'ln=', ln
      allocate (qq_0(in, jn, 5))

      read (12) (((qq_0(i, j, l), i=1, in), j=1, jn), l=1, 5)
      close (12)

      print *, 'read base flow is done!'

   end subroutine plot3d_load

   subroutine read_curvature()
      implicit none
      integer :: ilength, i
      logical :: alive
  
      allocate(drx(in), drz(in))
      inquire(file=char(curfile), exist=alive)
      print*, 'reading curvature file!', char(curfile)
      if(alive)then
        open(13, file=char(curfile))
        read(13, *)
        read(13, *)ilength
        do i=1, ilength
          read(13, *)drx(i), drz(i)
        enddo
        ! do i=1, in
        !     read(13)drx(i), drz(i)
        !     print*, i, drx(i), drz(i)
        ! enddo
        close(13)
      else
        drx=drxz(1); drz=drxz(2)
      endif
      print*, 'reading curvature file done!'
  
   end subroutine read_curvature

   subroutine read_mesh_2d(comm)

      implicit none
      PetscInt, INTENT(in) :: comm
      Vec :: Coord_0
      VecScatter :: ctx0, ctx
      PetscScalar, pointer :: tmp_(:)
      real(REAL64) :: loader(5, in, jn)
      integer :: i, j, l

      call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
      &                 DMDA_STENCIL_BOX, in, jn, 1, PETSC_DECIDE, &
      &                 1, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, DA, ierr)
      call DMSetUp(DA, ierr)

      ! call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
      ! &                 DMDA_STENCIL_BOX, in, jn,  1, PETSC_DECIDE, &
      ! &                 2, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, coordDA, ierr)
      ! call DMSetUp(coordDA, ierr)

      ! call DMGetCoordinates(DA, Coord, ierr)
      ! call VecZeroEntries(Coord, ierr)
      ! call DMGetCoordinatesLocal(DA, Coord_local, ierr)
      ! call DMGetGlobalVector(coordDA, Coord, ierr)
      ! call DMGetLocalVector(coordDA, Coord_local, ierr)

      ! call VecScatterCreateToZero(Coord, ctx0, Coord_0, ierr)
      ! call VecScatterCreate(Coord_0, PETSC_NULL_IS, Coord, PETSC_NULL_IS,ctx, ierr)

      !   block
      !     integer :: isize
      !   call VecGetSize(Coord, isize, ierr)
      !   print*, rank,'isize_coord=', isize
      !   call VecGetSize(Coord, isize, ierr)
      !   print*, rank,'isize_coord_0=', isize
      !
      !   if(rank==0)then
      !     call VecGetArrayF90(Coord_0, tmp_, ierr)
      !     ! print*, shape(tmp_)
      !     ! do i=1, in
      !     ! do j=1, jn
      !     ! do l=1, 2
      !     do i=1, in
      !       do j=1, jn
      !         loader(1, i, j)=xx_0(i, j)
      !         loader(2, i, j)=yy_0(i, j)
      !       enddo
      !     enddo
      !     print*, loader(1, 1:5, 2)
      !     tmp_=reshape(loader(1:2, 1:in, 1:jn), [2*in*jn])
      !     ! enddo
      !     ! enddo
      !     ! enddo
      !     call VecRestoreArrayF90(Coord_0, tmp_, ierr)
      !     ! call VecView(Coord_0, PETSC_VIEWER_STDOUT_SELF, ierr)
      !
      !   endif
      ! endblock
      !   call VecScatterBegin(ctx, coord_0, Coord, INSERT_VALUES, SCATTER_FORWARD, ierr)
      !   call VecScatterEnd(ctx, coord_0, Coord, INSERT_VALUES, SCATTER_FORWARD, ierr)
      !
      !   if(rank==0)then
      !
      !   endif

      call DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &
      &                 DMDA_STENCIL_BOX, in, jn, 1, PETSC_DECIDE, &
      &                 5, 2, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, meshDA, ierr)
      call DMSetUp(meshDA, ierr)
      call DMGetCoordinateDM(MeshDA, CoordDA, ierr)

      call DMGetGlobalVector(CoordDA, Coord, ierr)
      call DMSetCoordinates(DA, Coord, ierr)
      call DMGetGlobalVector(meshDA, Flowfield, ierr)
      call DMGetLocalVector(meshDA, Flowfield_local, ierr)

      ! call PetscViewerBinaryOpen(comm, "in/grid.petsc",FILE_MODE_READ, Viewer, ierr)

      ! call VecLoad(Coord, Viewer, ierr)
      ! ! stop
      ! ! call DMGlobalToLocalBegin(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
      ! ! call DMGlobalToLocalEnd(coordDA, Coord, INSERT_VALUES, Coord_local, ierr)
      ! call PetscViewerDestroy(Viewer, ierr)

      call PetscViewerHDF5Open(comm, char(gridfile), FILE_MODE_READ, viewer, ierr)
      call PetscObjectSetName(coord, "grid", ierr)
      call VecLoad(coord, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)

      ! call PetscViewerBinaryOpen(comm, "in/flow.petsc",FILE_MODE_READ, Viewer, ierr)
      call PetscViewerHDF5Open(comm, char(flowfile), FILE_MODE_READ, viewer, ierr)
      call PetscObjectSetName(Flowfield, "baseflow", ierr)
      call VecLoad(Flowfield, Viewer, ierr)
      call PetscViewerDestroy(Viewer, ierr)

      call DMGlobalToLocalBegin(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
      call DMGlobalToLocalEnd(meshDA, Flowfield, INSERT_VALUES, Flowfield_local, ierr)
      call DMGetLocalVector(coordDA, coord_local, ierr)
      call DMGlobalToLocalBegin(CoordDA, Coord, INSERT_VALUES, Coord_local, ierr)
      call DMGlobalToLocalEnd(CoordDA, Coord, INSERT_VALUES, Coord_local, ierr)

      call get_layout()

   end subroutine read_mesh_2d

   subroutine get_layout()

      implicit none
      call DMDAGetGhostCorners(coordDA, xbg, ybg, PETSC_NULL_INTEGER,                       &
      &                  xlg, ylg, PETSC_NULL_INTEGER, ierr)
      xeg = xbg + xlg - 1; yeg = ybg + ylg - 1
      call DMDAGetCorners(coordDA, xb, yb, PETSC_NULL_INTEGER,                       &
      &                  xl, yl, PETSC_NULL_INTEGER, ierr)
      xe = xb + xl - 1; ye = yb + yl - 1

   end subroutine get_layout

   subroutine load_mesh_info(comm)

      implicit none
      PetscInt, INTENT(in) :: comm
      integer :: i, j, l

      allocate (xx(xbg + 1:xeg + 1, ybg + 1:yeg + 1))
      allocate (yy(xbg + 1:xeg + 1, ybg + 1:yeg + 1))

      call DMDAVecGetArrayReadF90(coordDA, coord_local, grid, ierr)

      do j = ybg, yeg
         do i = xbg, xeg
            xx(i + 1, j + 1) = grid(0, i, j)
            yy(i + 1, j + 1) = grid(1, i, j)
         end do
      end do

      call DMDAVecRestoreArrayReadF90(coordDA, coord_local, grid, ierr)

   end subroutine load_mesh_info

   subroutine load_flow_info(comm)

      implicit none
      PetscInt, INTENT(in) :: comm
      integer :: i, j, l

      allocate (qq(5, xbg + 1:xeg + 1, ybg + 1:yeg + 1))
      call DMDAVecGetArrayReadF90(meshDA, Flowfield_local, flow, ierr)

      do l = 0, 4
         do j = ybg, yeg
            do i = xbg, xeg
               qq(l + 1, i + 1, j + 1) = flow(l, i, j)
            end do
         end do
      end do

      call DMDAVecRestoreArrayReadF90(meshDA, Flowfield_local, flow, ierr)

   end subroutine load_flow_info

end module mod_petsc_loader
