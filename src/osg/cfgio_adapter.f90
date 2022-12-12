!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-09 10:45:42
!> @LastEditTime: 2022-12-12 19:29:04
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/osg/cfgio_adapter.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

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

        if(rank==0)then
            cfg=parse_cfg(trim(cfgfn))
            if(cfg%has_key("general", "screen")) &
            call cfg%get("general", "screen", screen_flg)
            call cfg_filename(cfg)
            call cfg_domain(cfg)
            call cfg_freestream(cfg)
            call cfg_3dsog(cfg)
            call cfg_lpse(cfg)
            !call cfg_bc(cfg)
            call cfg%write()
            call cfg%write(trim(cfgfn))
        endif
        ! call bcast()
        ! if(mpi_ierr) then
        !     print*, 'There is some fatal errors in ini file.'
        !     call MPI_Finalize(ierr)
        !     stop
        ! endif

    end subroutine cfg_loader

    !> read filenames
    !> @param[in] cfg cfg_type
    subroutine cfg_filename(cfg)

        implicit none
        type(cfg_t) :: cfg
        character(len=256) :: grid
        character(len=256) :: flow

        if(cfg%has_key("filenames", "dir")) then
          call cfg%get("filenames", "dir", dir)
        else
          Dir="./"
        endif
        call cfg%get("filenames", "grid", grid)
        call cfg%get("filenames", "flow", flow)

        gridfile=string_type(trim(dir)//trim(grid))
        flowfile=string_type(trim(dir)//trim(flow))

    end subroutine cfg_filename

    subroutine cfg_domain(cfg)

      implicit none
      type(cfg_t) :: cfg
      real(REAL64), allocatable :: tmp(:)
      integer :: npar
      character(256)::chardump

      call cfg%get("domain", "in", in)
      call cfg%get("domain", "jn", jn)

      if(cfg%has_key("domain", "curvature"))then
        npar=2
        call cfg%get("domain", "curvature", tmp, npar)
        drxz=tmp(1:2)
      else if(cfg%has_key("filenames", "curvature"))then
        call cfg%get("filenames", "curvature", chardump)
        curfile=trim(dir)//trim(chardump)
      else
        drxz=0.0d0
      endif

      if(cfg%has_key("domain", "istart"))then
        call cfg%get("domain", "istart", istart)
      else
        print*, '[domain]istart is not given, use istart=1.'
        istart=1
        call cfg%set("domain", "istart", istart)
      endif
      if(cfg%has_key("domain", "iend"))then
        call cfg%get("domain", "iend", iend)
      else
        print*, '[domain]iend is not given, use iend=in.'
        iend=in
        call cfg%set("domain", "iend", iend)
      endif
      if(cfg%has_key("domain", "iloc"))then
        call cfg%get("domain", "iloc", iloc)
      else
        print*, '[domain]iloc is not given, use iloc=in-2.'
        iloc=in-2
        call cfg%set("domain", "iloc", iloc)
      endif

    end subroutine cfg_domain

    !> 来流信息提取函数
    !> @param[in] cfg cfg类
    !> @return 返回配置文件freestream字段相关信息
    subroutine cfg_freestream(cfg)

        implicit none
        type(cfg_t) :: cfg

        if(cfg%has_key("freestream", "ma"))then
          call cfg%get("freestream", "ma", Ma)
        else
          print*, 'No Mach number is input,'
          Ma=0.001d0
          print*, 'The Ma is set to 0.001.'
        endif
        if(cfg%has_key("freestream", "re"))then
          call cfg%get("freestream", "re", Re)
        else
          print*, 'No Reynolds number is input,'
          Re=huge(1.0d0)
          print*, 'The Re is set to 1.0d20.'
        endif
        if(cfg%has_key("freestream", "te"))then
          call cfg%get("freestream", "te", Te)
        else
          print*, 'No freestream temperature is input,'
          Te=300.0d0
          print*, 'The Te is set to 300K.'
        endif
        if(cfg%has_key("freestream", "pr"))then
          call cfg%get("freestream", "pr", Pr)
        else
          print*, 'No freestream Prandtl number is input,'
          Pr=0.72d0
          print*, 'The Pr is set to 0.72.'
        endif
        
    end subroutine cfg_freestream

    subroutine cfg_lpse(cfg)

      implicit none
      type(cfg_t) :: cfg

      !! this parameter control the dp/dx in LPSE, generally it can be set 1.
      if(cfg%has_key('lpse', 'sigma_p'))then
        call cfg%get('lpse', 'sigma_p', sigma_p)
      else
        sigma_p=1.0d0
      endif

      if(sigma_p>99.0d0) then
        sigma_p=0.0d0
        sigma_control=1.0d0
      else
        sigma_control=0.0d0
      endif

    end subroutine cfg_lpse

    subroutine cfg_3dsog(cfg)

      implicit none
      type(cfg_t) :: cfg
      character(len=256):: method_local, mode_local
      integer :: npar
      real(REAL64), allocatable, dimension(:) :: alf_in, ome_in, bat_in

      if(cfg%has_key('3dsog', 'i_x0'))then
        call cfg%get('3dsog', 'i_x0', i_x0)
      else
        i_x0=istart
        call cfg%set('3dsog', 'i_x0', i_x0)
      endif

      if(cfg%has_key('3dsog', 'i_x1'))then
        call cfg%get('3dsog', 'i_x1', i_x1)
      else
        i_x1=iend
        call cfg%set('3dsog', 'i_x1', i_x1)
      endif

      if(cfg%has_key('3dsog', 'initial_alpha'))then
        call cfg%get('3dsog', 'initial_alpha', initial_alf)
      else
        stop '[3dsog]initial_alpha should be given.'
      endif

      if(cfg%has_key('3dsog', 'initial_beta'))then
        call cfg%get('3dsog', 'initial_beta', initial_bat)
      else
        stop '[3dsog]initial_beta should be given.'
      endif

      if(cfg%has_key('3dsog', 'initial_omega'))then
        call cfg%get('3dsog', 'initial_omega', initial_ome)
      else
        initial_ome=0.0d0
        print*, '[3dsog]The initial_omega is not given, and the default value zero is adoped.'
      endif

    end subroutine cfg_3dsog

    subroutine cfg_bc(cfg)

        implicit none
        type(cfg_t) :: cfg
        character(256) :: bc_local
        type(string_type) :: bc_str

        if(cfg%has_key("bc", "y_0"))then
            call cfg%get("bc", 'y_0', bc_local)
            bc_str=trim(bc_local)
            select case (char(to_upper(bc_str)))
            case ('BC_WALL')
                print*, 'bc wall is used at y0.'
                bc_y0=BC_WALL
                if(cfg%has_key("bc", "temperature_type"))then
                  call cfg%get("bc", "temperature_type", isothermal)
                else
                  isothermal = .TRUE.
                endif
            case default
                print*, '[bc]y_0 keyword is wrong.'
                print*, 'The keryword may be BC_WALL.'
                mpi_ierr=.TRUE.
            end select
        else
            print*, '[bc]y_0 is not given, use y_0=BC_WALL'
            call cfg%set("BC", 'y_0', 'BC_WALL')
            bc_y0=BC_WALL
        endif

        if(cfg%has_key("bc", "y_inf"))then
            call cfg%get("bc", 'y_inf', bc_local)
            bc_str=trim(bc_local)
            select case (char(to_upper(bc_str)))
            case ('BC_FARFIELD')
                print*, 'bc farfield is used at yn.'
                bc_yn=BC_FARFIELD
            case default
                print*, '[bc]y_inf keyword is wrong.'
                print*, 'The keryword may be BC_FARFIELD.'
                mpi_ierr=.TRUE.
            end select
        else
            bc_yn=BC_FARFIELD
            print*, '[bc]y_inf is not given, use y_inf=BC_FARFIELD'
            call cfg%set("bc", 'y_inf', 'BC_FARFIELD')
        endif

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
    !     ntmp(7)=0
    !     ntmp(8)=0
    !     ntmp(9)=0
    !     ntmp(10) =bc_y0
    !     ntmp(11)=bc_yn

    !     call MPI_Bcast(ntmp(1), 100, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    !     in=ntmp(1)
    !     jn=ntmp(2)

    !     istart=ntmp(4)
    !     iend=ntmp(5)
    !     iloc=ntmp(6)
    !     ! lpse_method=ntmp(7)
    !     ! lpse_mode=ntmp(8)
    !     ! lpse_nev=ntmp(9)
    !     bc_y0=ntmp(10)
    !     bc_yn=ntmp(11)

    !     rtmp=0

    !     rtmp(1) =Ma
    !     rtmp(2) =Re
    !     rtmp(3) =Te
    !     rtmp(4) =alpha(1)
    !     rtmp(5) =alpha(2)
    !     rtmp(6) =omega(1)
    !     rtmp(7) =omega(2)
    !     rtmp(8) =beta (1)
    !     rtmp(9) =beta (2)
    !     rtmp(10)=sigma_p

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
    !     sigma_p =rtmp(10)

    !     ctmp(1)=char(gridfile)
    !     ctmp(2)=char(flowfile)

    !     call MPI_bcast(ctmp(1), 100*256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    !     gridfile=trim(ctmp(1))
    !     flowfile=trim(ctmp(2))

    !     call MPI_Bcast(lpse_IsAdjoint, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    !     call MPI_Bcast(screen_flg, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    !     call MPI_bcast(mpi_ierr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    ! end subroutine bcast

end module mod_cfgio_adapter
