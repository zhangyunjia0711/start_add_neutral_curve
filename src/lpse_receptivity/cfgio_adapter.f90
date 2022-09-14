!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!
!  File: cfgio_adapter.f90
!> @file
!> @breif cfgio代码适配器文件.
!  DESCRIPTION:
!>
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
      logical :: set

      if (rank == 0) then
         cfg = parse_cfg(trim(cfgfn))
         call cfg_general(cfg)
         call cfg_filename(cfg)
         call cfg_domain(cfg)
         call cfg_freestream(cfg)
         call cfg_lpse(cfg)
         call cfg_bc(cfg)
         call cfg_receptivity(cfg)
         call cfg%write()
         call cfg%write(trim(cfgfn))
      end if
      ! call bcast()
      ! if(mpi_ierr) then
      !     print*, 'There is some fatal errors in ini file.'
      !     call MPI_Finalize(ierr)
      !     stop
      ! endif

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-rd_factor', rd_factor, set, ierr)

   end subroutine cfg_loader

   subroutine cfg_general(cfg)
      implicit none
      type(cfg_t) :: cfg
      character(len=256) :: cs

      if (cfg%has_key("general", "screen")) &
         call cfg%get("general", "screen", screen_flg)

      if (cfg%has_key("general", "coordinate_system")) then
         call cfg%get("general", "coordinate_system", cs)
         select case (trim(cs))
         case ('comving')
            Coordinate_System = COMOVING_CS
         case ('cylindrical')
            Coordinate_System = CYLINDRICAL_CS
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
         Dir = "./"
      end if
      call cfg%get("filenames", "grid", grid)
      call cfg%get("filenames", "flow", flow)

      gridfile = string_type(trim(dir)//trim(grid))
      flowfile = string_type(trim(dir)//trim(flow))

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
         curfile = trim(dir)//trim(chardump)
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

   subroutine cfg_lpse(cfg)

      implicit none
      type(cfg_t) :: cfg
      character(len=256):: method_local, mode_local, cdump
      type(string_type) :: method_str, mode_str
      integer :: npar
      real(REAL64), allocatable, dimension(:) :: alf_in, ome_in, bat_in

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

      ! if(cfg%has_key('lpse', 'adjoint'))then
      !   call cfg%get('lpse', 'adjoint', lpse_IsAdjoint)
      !   if(lpse_IsAdjoint)then
      !     print*, 'ALPSE is solved.'
      !   else
      !     print*, 'LPSE is solved.'
      !   endif
      ! else
      !   print*, '[lpse]keyword adjoint is not given, use adjoint=.False.'
      !   call cfg%set('lpse', 'adjoint', .FALSE.)
      !   lpse_IsAdjoint=.False.
      ! endif

      ! !! this parameter control the dp/dx in LPSE, generally it can be set 1.
      ! if(cfg%has_key('lpse', 'sigma_p'))then
      !   call cfg%get('lpse', 'sigma_p', sigma_p)
      ! else
      !   sigma_p=1.0d0
      ! endif

      ! if(sigma_p>99.0d0) then
      !   sigma_p=0.0d0
      !   sigma_control=1.0d0
      ! else
      !   sigma_control=0.0d0
      ! endif
      sigma_p = 1.0d0; sigma_control = 0.0d0

      if (cfg%has_key('lpse', 'solve_alpha')) then
         call cfg%get('lpse', 'solve_alpha', cdump)
         select case (trim(cdump))
         case ('real')
            solve_alpha = SOLVE_ALPHA_REAL
         case ('both')
            solve_alpha = SOLVE_ALPHA_BOTH
         case ('none')
            solve_alpha = SOLVE_ALPHA_NONE
         end select
      else
         solve_alpha = SOLVE_ALPHA_BOTH
      end if

   end subroutine cfg_lpse

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
         case default
            print *, '[bc]y_0 keyword is wrong.'
            print *, 'The keryword may be BC_WALL.'
            mpi_ierr = .TRUE.
         end select
      else
         print *, '[bc]y_0 is not given, use y_0=BC_WALL'
         call cfg%set("BC", 'y_0', 'BC_WALL')
         bc_y0 = BC_WALL
      end if

      if (cfg%has_key("bc", "y_inf")) then
         call cfg%get("bc", 'y_inf', bc_local)
         bc_str = trim(bc_local)
         select case (char(to_upper(bc_str)))
         case ('BC_FARFIELD')
            print *, 'bc farfield is used at yn.'
            bc_yn = BC_FARFIELD
         case default
            print *, '[bc]y_inf keyword is wrong.'
            print *, 'The keryword may be BC_FARFIELD.'
            mpi_ierr = .TRUE.
         end select
      else
         bc_yn = BC_FARFIELD
         print *, '[bc]y_inf is not given, use y_inf=BC_FARFIELD'
         call cfg%set("bc", 'y_inf', 'BC_FARFIELD')
      end if

   end subroutine cfg_bc

   subroutine cfg_receptivity(cfg)

      implicit none
      type(cfg_t), intent(inout) :: cfg
      character(len=256) :: chartmp
      real(REAL64), allocatable :: rtmp(:)
      integer :: npar

      if (cfg%has_key("receptivity", "source")) then
         call cfg%get("receptivity", "source", chartmp)

         select case (trim(chartmp))
         case ("none", "zero")
            source_term_file = 'none'
         case default
            source_term_file = trim(dir)//trim(chartmp)
         end select
      else
         call PetscPrintf(PETSC_COMM_WORLD, "No source term in receptivity.\n", ierr)
         source_term_file = 'none'
      end if

      if (cfg%has_key("receptivity", "wall")) then
         call cfg%get("receptivity", "wall", chartmp)
         select case (trim(chartmp))
         case ("none", "zero")
            wall_force_file = 'none'
         case default
            wall_force_file = trim(dir)//trim(chartmp)
         end select
      else
         call PetscPrintf(PETSC_COMM_WORLD, "No wall bc in receptivity.\n", ierr)
         wall_force_file = 'none'
      end if

   end subroutine cfg_receptivity

end module mod_cfgio_adapter
