!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-08 13:45:30
!> @LastEditTime: 2022-12-12 19:33:01
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/lst/cfgio_adapter.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

!@ config file adapter for lst solver.
module mod_cfgio_adapter

   use cfgio_mod
   use global_parameters
   use stdlib_string_type

   implicit none

   public :: cfg_loader

   private
   character(len=256) :: Dir

#include "../include/start.f90"

contains

   !> read infomation from cfg_file
    !! @param[in] cfgfn filename of cfgio
   subroutine cfg_loader(cfgfn)

      implicit none
      character(len=*), intent(in) :: cfgfn
      integer :: ierr

      if (rank == 0) then
         cfg = parse_cfg(trim(cfgfn))
         call cfg_general(cfg)
         call cfg_filename(cfg)
         call cfg_domain(cfg)
         call cfg_freestream(cfg)
         call cfg_lst(cfg)
         call cfg_bc(cfg)
         call cfg%write()
         call cfg%write(trim(cfgfn))
      end if
      ! call bcast()
      ! if(mpi_ierr) then
      !     print*, 'There is some fatal errors in ini file.'
      !     call MPI_Finalize(ierr)
      !     stop
      ! endif

   end subroutine cfg_loader

   subroutine cfg_general(cfg)
      implicit none
      type(cfg_t) :: cfg
      character(len=256) :: cs

      if (cfg%has_key("general", "coordinate_system")) then
         call cfg%get("general", "coordinate_system", cs)
         select case (trim(cs))
         case ('comving')
            Coordinate_System = COMOVING_CS
            call PetscPrintf(PETSC_COMM_WORLD, "comoving cs is used.\n", ierr)
         case ('cylindrical')
            Coordinate_System = CYLINDRICAL_CS
            call PetscPrintf(PETSC_COMM_WORLD, "cylinderical cs is used.\n", ierr)
         case default
            call PetscPrintf(PETSC_COMM_WORLD, "the coordinate system is wrong, and comoving cs is used!\n", ierr)
            Coordinate_System = COMOVING_CS
         end select
      else
         Coordinate_System = COMOVING_CS
      end if

   end subroutine cfg_general

   !> read filenames
   !> @param[in] cfg cfg_type
   subroutine cfg_filename(cfg)

      implicit none
      type(cfg_t) :: cfg
      character(len=256) :: grid
      character(len=256) :: flow

      if (cfg%has_key("filenames", "dir")) then
         call cfg%get("filenames", "dir", dir)
      else
         dir = "./"
      end if
      call cfg%get("filenames", "grid", grid)
      call cfg%get("filenames", "flow", flow)

      gridfile = trim(dir)//trim(grid)
      flowfile = trim(dir)//trim(flow)

   end subroutine cfg_filename

   subroutine cfg_domain(cfg)

      implicit none
      type(cfg_t) :: cfg
      real(REAL64), allocatable :: tmp(:)
      integer :: npar
      character(256)::chardump

      call cfg%get("domain", "in", in)
      call cfg%get("domain", "jn", jn)

      if (cfg%has_key("domain", "curvature")) then
         npar = 2
         call cfg%get("domain", "curvature", tmp, npar)
         drxz = tmp(1:2)
      else if (cfg%has_key("filenames", "curvature")) then
         call cfg%get("filenames", "curvature", chardump)
         if (trim(chardump) /= 'null') then
            curfile = trim(dir)//trim(chardump)
         else
            drxz = 0.0d0
         end if
      else
         drxz = 0.0d0
      end if

      if (cfg%has_key("domain", "istart")) then
         call cfg%get("domain", "istart", istart)
      else
         print *, '[domain]istart is not given, use istart=1.'
         istart = 1
         call cfg%set("domain", "istart", istart)
      end if
      if (cfg%has_key("domain", "iend")) then
         call cfg%get("domain", "iend", iend)
      else
         print *, '[domain]iend is not given, use iend=in.'
         iend = in
         call cfg%set("domain", "iend", iend)
      end if
      if (cfg%has_key("domain", "iloc")) then
         call cfg%get("domain", "iloc", iloc)
      else
         print *, '[domain]iloc is not given, use iloc=in-2.'
         iloc = in - 2
         call cfg%set("domain", "iloc", iloc)
      end if

   end subroutine cfg_domain

   !> 来流信息提取函数
   !> @param[in] cfg cfg类
   !> @return 返回配置文件freestream字段相关信息
   subroutine cfg_freestream(cfg)

      implicit none
      type(cfg_t) :: cfg

      if (cfg%has_key("freestream", "ma")) then
         call cfg%get("freestream", "ma", Ma)
      else
         print *, 'No Mach number is input,'
         Ma = 0.001d0
         print *, 'The Ma is set to 0.001.'
      end if
      if (cfg%has_key("freestream", "re")) then
         call cfg%get("freestream", "re", Re)
      else
         print *, 'No Reynolds number is input,'
         Re = huge(1.0d0)
         print *, 'The Re is set to 1.0d20.'
      end if
      if (cfg%has_key("freestream", "te")) then
         call cfg%get("freestream", "te", Te)
      else
         print *, 'No freestream temperature is input,'
         Te = 300.0d0
         print *, 'The Te is set to 300K.'
      end if
      if (cfg%has_key("freestream", "pr")) then
         call cfg%get("freestream", "pr", Pr)
      else
         print *, 'No freestream Prandtl number is input,'
         Pr = 0.72d0
         print *, 'The Pr is set to 0.72.'
      end if

   end subroutine cfg_freestream

   subroutine cfg_lst(cfg)

      implicit none
      type(cfg_t) :: cfg
      character(len=256):: method_local, mode_local
      type(string_type) :: method_str, mode_str
      integer :: npar
      real(REAL64), allocatable, dimension(:) :: alf_in, ome_in, bat_in

      if (cfg%has_key("lst", "method")) then
         call cfg%get("lst", "method", method_local)
         method_str = trim(method_local)
         ! method_local=upper(method_local)
         select case (char(to_upper(method_str)))
         case ('LST_GLOBAL')
            lst_method = LST_GLOBAL
            print *, 'LST_GLOBAL is used.'
         case ('LST_LOCAL')
            lst_method = LST_LOCAL
            print *, 'LST_LOCAL is used.'
         case ('LST_DOMAIN')
            lst_method = LST_DOMAIN
            print *, 'LST_LOCAL is used.'
         case default
            print *, '[lst]method keyword is wrong.'
            print *, 'The keryword should be lst_global, lst_local, lst_domain.'
            mpi_ierr = .TRUE.
         end select
      else
         print *, '[lst]method is not given, use method=LST_GLOBAL'
         lst_method = LST_GLOBAL
         method_local = 'LST_GLOBAL'
         call cfg%set("lst", "method", trim(method_local))
      end if

      if (lst_method == LST_GLOBAL) then
         if (cfg%has_key("lst", "nev")) then
            call cfg%get("lst", "nev", lst_nev)
         else
            print *, '[Warning:]When select LST_GLOBAL, nev should be given.'
            print *, '[lst]nev is not given, use nev=<num>'
            mpi_ierr = .TRUE.
         end if
      end if

      if (cfg%has_key("lst", "mode")) then
         call cfg%get("lst", "mode", mode_local)
         mode_str = trim(mode_local)
         select case (char(to_upper(mode_str)))
         case ('LST_TEMPORAL_MODE')
            lst_mode = LST_TEMPORAL_MODE
            print *, 'LST_TEMPORAL_MODE is used.'
         case ('LST_SPATIAL_MODE')
            lst_mode = LST_SPATIAL_MODE
            print *, 'LST_SPATIAL_MODE is used.'
         case default
            print *, char(to_upper(mode_str))
            print *, '[LST]mode keyword is wrong.'
            print *, 'The keryword may be lst_spatial_mode, lst_temporal_mode.'
            mpi_ierr = .TRUE.
         end select
      else
         print *, '[lst]method is not given, use mode=LST_TEMPORAL_MODE'
         lst_mode = LST_TEMPORAL_MODE
         mode_local = 'LST_TEMPORAL_MODE'
         call cfg%set("lst", "mode", trim(mode_local))
      end if

      npar = 2

      if (cfg%has_key('lst', "alpha")) then
         call cfg%get('lst', "alpha", alf_in, npar)
         alpha = alf_in
      else
         print *, '[lst]keyword alpha should be given as alpha=Re(alpha), Im(alpha).'
         mpi_ierr = .TRUE.
      end if
      if (cfg%has_key('lst', "omega")) then
         call cfg%get('lst', "omega", ome_in, npar)
         omega = ome_in
      else
         print *, '[lst]keyword omega should be given as omega=Re(omega), Im(omega).'
         mpi_ierr = .TRUE.
      end if
      if (cfg%has_key('lst', "beta")) then
         call cfg%get('lst', "beta", bat_in, npar)
         beta = bat_in
      else
         print *, '[lst]keyword beta should be given as beta=Re(beta), Im(beta).'
         mpi_ierr = .TRUE.
      end if

      if (cfg%has_key('lst', 'adjoint')) then
         call cfg%get('lst', 'adjoint', lst_IsAdjoint)
         if (lst_IsAdjoint) then
            print *, 'ALST is solved.'
         else
            print *, 'LST is solved.'
         end if
      else
         print *, '[lst]keyword adjoint is not given, use adjoint=.False.'
         call cfg%set('lst', 'adjoint', .FALSE.)
         lst_IsAdjoint = .False.
      end if

   end subroutine cfg_lst

   subroutine cfg_bc(cfg)

      implicit none
      type(cfg_t) :: cfg
      character(256) :: bc_local
      type(string_type) :: bc_str

      if (cfg%has_key("bc", "y_0")) then
         call cfg%get("bc", 'y_0', bc_local)
         bc_str = trim(bc_local)
         select case (char(to_upper(bc_str)))
         case ('BC_WALL')
            print *, 'bc wall is used at y0.'
            bc_y0 = BC_WALL
            if (cfg%has_key("bc", "temperature_type")) then
               call cfg%get("bc", "temperature_type", isothermal)
            else
               isothermal = .TRUE.
            end if
         case ('BC_POLE')
            print *, 'bc pole is used at y0.'
            bc_y0 = BC_POLE
         case ('BC_POROUS')
            print *, 'bc porous is used at y0.'
            bc_y0 = BC_POROUS
            block
               use mod_porous_parameters
               call cfg%get("bc", "porous_porosity", porosity)
               call cfg%get("bc", "porous_radius", radius)
               call cfg%get("bc", "porous_depth", depth)
            end block
         case default
            print *, '[bc]y_0 keyword is wrong.'
            print *, 'The keryword may be BC_WALL.'
            mpi_ierr = .TRUE.
         end select
      else
         print *, '[bc]y_0 is not given, use y_0=BC_WALL'
         call cfg%set("bc", 'y_0', 'BC_WALL')
         bc_y0 = BC_WALL
      end if

      if (cfg%has_key("bc", "y_inf")) then
         call cfg%get("bc", 'y_inf', bc_local)
         bc_str = trim(bc_local)
         select case (char(to_upper(bc_str)))
         case ('BC_FARFIELD')
            print *, 'bc farfield is used at yn.'
            bc_yn = BC_FARFIELD
         case ("BC_NONREFLECTION")
            print *, 'bc nonreflection is used at yn.'
            bc_yn = BC_NONREFLECTION
         case default
            print *, '[bc]y_inf keyword is wrong.'
            print *, 'The keryword may be BC_FARFIELD.'
            mpi_ierr = .TRUE.
         end select
      else
         bc_yn = BC_FARFIELD
         print *, '[bc]y_inf is not given, use y_inf=BC_FARFIELD'
         call cfg%set("BC", 'y_inf', 'BC_FARFIELD')
      end if

   end subroutine cfg_bc

   ! subroutine bcast()

   !     implicit none
   !     integer :: ntmp(100)
   !     real(REAL64) :: rtmp(100)
   !     character(256) :: ctmp(100)

   !     integer :: ierr

   !     ntmp=0

   !     ntmp(1)=in
   !     ntmp(2)=jn

   !     ntmp(4)=istart
   !     ntmp(5)=iend
   !     ntmp(6)=iloc
   !     ntmp(7)=lst_method
   !     ntmp(8)=lst_mode
   !     ntmp(9)=lst_nev
   !     ntmp(10) =bc_y0
   !     ntmp(11)=bc_yn

   !     call MPI_Bcast(ntmp(1), 100, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   !     in=ntmp(1)
   !     jn=ntmp(2)

   !     istart=ntmp(4)
   !     iend=ntmp(5)
   !     iloc=ntmp(6)
   !     lst_method=ntmp(7)
   !     lst_mode=ntmp(8)
   !     lst_nev=ntmp(9)
   !     bc_y0=ntmp(10)
   !     bc_yn=ntmp(11)

   !     rtmp=0

   !     rtmp(1)=Ma
   !     rtmp(2)=Re
   !     rtmp(3)=Te
   !     rtmp(4)=alpha(1)
   !     rtmp(5)=alpha(2)
   !     rtmp(6)=omega(1)
   !     rtmp(7)=omega(2)
   !     rtmp(8)=beta (1)
   !     rtmp(9)=beta (2)

   !     call MPI_bcast(rtmp(1), 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   !     Ma=rtmp(1)
   !     Re=rtmp(2)
   !     Te=rtmp(3)
   !     alpha(1)=rtmp(4)
   !     alpha(2)=rtmp(5)
   !     omega(1)=rtmp(6)
   !     omega(2)=rtmp(7)
   !     beta (1)=rtmp(8)
   !     beta (2)=rtmp(9)

   !     ctmp(1)=char(gridfile)
   !     ctmp(2)=char(flowfile)

   !     call MPI_bcast(rtmp(1), 100*256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

   !     gridfile=trim(ctmp(1))
   !     flowfile=trim(ctmp(2))

   !     call MPI_Bcast(lst_IsAdjoint, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   !     call MPI_bcast(mpi_ierr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

   ! end subroutine bcast

end module mod_cfgio_adapter
