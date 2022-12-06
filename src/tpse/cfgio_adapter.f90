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
            ! call cfg_lpse(cfg)
            call cfg_npse(cfg)
            call cfg_bc(cfg)
            call cfg%write()
            call cfg%write("config.ini.see")
            ! call cfg%write(trim(cfgfn))
        end if
        ! call bcast()
        ! if (mpi_ierr) then
        !     print *, 'There is some fatal errors in ini file.'
        !     call MPI_Finalize(ierr)
        !     stop
        ! end if

        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-rd_factor', rd_factor, set, ierr)
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-rtol_nonlinear', rtol_nonlinear, set, ierr)
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atol_nonlinear', atol_nonlinear, set, ierr)
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-maxit_alf', maxit_alpha, set, ierr)
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-maxit_nonlinear', maxit_nonlinear, set, ierr)
        
    end subroutine cfg_loader

    subroutine cfg_general(cfg)
        implicit none
        type(cfg_t) :: cfg
        character(len=256) :: cs
        character(len=256) :: chardump

        if (cfg%has_key("general", "screen")) &
            call cfg%get("general", "screen", screen_flg)
         
         if (cfg%has_key("general", "title")) then
            call cfg%get("general", "title", chardump)
            title=trim(chardump)
         else
            title='untitled case'
         endif

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
        character(len=256) :: flow, pre

        if (cfg%has_key("filenames", "dir")) then
            call cfg%get("filenames", "dir", dir)
        else
            Dir = "./"
        end if
        call cfg%get("filenames", "grid", grid)
        call cfg%get("filenames", "flow", flow)

        gridfile = trim(dir)//trim(grid)
        flowfile = trim(dir)//trim(flow)

        if (cfg%has_key("filenames", "output prefix"))then
         call cfg%get("filenames", "output prefix", pre)
         prefix=trim(pre)
        else
         prefix='./out/'
        endif

    end subroutine cfg_filename

    subroutine cfg_domain(cfg)

        implicit none
        type(cfg_t) :: cfg
        real(REAL64), allocatable :: tmp(:)
        integer :: npar
        character(256)::char_dump

        call cfg%get("domain", "in", in)
        call cfg%get("domain", "jn", jn)

        if (cfg%has_key("domain", "curvature")) then
            npar = 2
            call cfg%get("domain", "curvature", tmp, npar)
            drxz = tmp(1:2)
        else if (cfg%has_key("filenames", "curvature")) then
            call cfg%get("filenames", "curvature", char_dump)
            curfile = trim(dir)//trim(char_dump)
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
            Re = 1.0d20
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
            call cfg%get("freestream", "pr", Te)
        else
            print *, 'No freestream Prandtl number is input,'
            Pr = 0.72d0
            print *, 'The Pr is set to 0.72.'
        end if

    end subroutine cfg_freestream

    subroutine cfg_npse(cfg)

        use stdlib_string_type
        implicit none
        type(cfg_t) :: cfg
        integer :: npar
        integer, allocatable :: integer_tmp(:)
        real(REAL64), allocatable :: real_tmp(:)
        real(REAL64) :: amp, angle
        character(len=256) :: cdump
        type(string_type) :: string_char
        character(len=256) :: format0
        real(REAL64) :: ome0, bat0
        integer :: i

        call cfg%get('tpse', 'm_global', m_global); m_global=m_global-1
        call cfg%get('tpse', 'n0_global', n0_global); n0_global=n0_global-1

        call cfg%get('tpse', 'theta1', theta1)
        call cfg%get('tpse', 'theta2', theta2)
        format0 = '(A, I4, 1X, A, I4, 1X, A, 2F10.5)'
        print trim(format0), 'The computational range is', 0, 'to', m_global, 'with theta1=', theta1
        print trim(format0), 'The computational range is', -n0_global, 'to', n0_global, 'with theta2=', theta2


        npar = 2

        call cfg%get('tpse', 'ndis', ndis)
        print *, 'The inlet disturbance number is', ndis
        print '(3(A4), 4(A10), (A14))', 'no.', 'm', 'n', 'alpha_r', 'alpha_i', 'amp', 'angle', 'solve_alpha'
        allocate (inlet_dis(ndis))
        do i = 1, ndis
            string_char = string_type(i)
            call cfg%get('tpse', 'dis_'//char(string_char), integer_tmp, npar)
            inlet_dis(i)%m = integer_tmp(1)
            inlet_dis(i)%n = integer_tmp(2)
            call cfg%get('tpse', 'dis_'//char(string_char)//'_alpha', real_tmp, npar)
            inlet_dis(i)%alf = real_tmp(1) + real_tmp(2)*ZI
            call cfg%get('tpse', 'dis_'//char(string_char)//'_amp', amp)
            if (cfg%has_key('tpse', 'dis_'//char(string_char)//'_angle')) then
                call cfg%get('tpse', 'dis_'//char(string_char)//'_angle', angle)
            else
                angle = 0.0d0
            end if
            inlet_dis(i)%amp = amp*cos(angle) + amp*sin(angle)*ZI
            if (cfg%has_key('tpse', 'dis_'//char(string_char)//'_solve_alpha')) then
                call cfg%get('tpse', 'dis_'//char(string_char)//'_solve_alpha', cdump)
                select case (trim(cdump))
                case ('real')
                    inlet_dis(i)%solve_alpha = SOLVE_ALPHA_REAL
                case ('both')
                    inlet_dis(i)%solve_alpha = SOLVE_ALPHA_BOTH
                case ('none')
                    inlet_dis(i)%solve_alpha = SOLVE_ALPHA_NONE
                case default
                    stop 'disturbance at inlet is wrong with a undefined solve alpha type.'
                end select
            else
                if (i == 1) then
                    inlet_dis(i)%solve_alpha = SOLVE_ALPHA_BOTH
                    cdump = 'both'
                else
                    inlet_dis(i)%solve_alpha = SOLVE_ALPHA_NONE
                    cdump = 'none'
                end if
            end if
            print '(3I4, 4(F10.4), A10)', i, inlet_dis(i)%m, inlet_dis(i)%n, &
                inlet_dis(i)%alf%re, inlet_dis(i)%alf%im, &
                amp, angle, trim(cdump)
        end do

        if (cfg%has_key('lpse', 'sigma_p')) then
            call cfg%get('lpse', 'sigma_p', sigma_p)
        else
            sigma_p = 1.0d0
        end if

        if (sigma_p > 99.0d0) then
            sigma_p = 0.0d0
            sigma_control = 1.0d0
        else
            sigma_control = 0.0d0
        end if

        if(cfg%has_key("tpse", "omega_force_factor")) then
            call cfg%get("tpse", "omega_force_factor", omega_force_factor)
        else
            omega_force_factor = 1.0d0
        endif

        if(cfg%has_key("tpse", "ix_force")) then
            call cfg%get("tpse", "ix_force", ix_force)
        else
            ix_force = istart
        endif

        if(cfg%has_key("tpse", "omega_phi_factor")) then
            call cfg%get("tpse", "omega_phi_factor", omega_phi_factor)
        else
            omega_phi_factor = 1.0d0
        endif

        if(cfg%has_key("tpse", "ix_phi")) then
            call cfg%get("tpse", "ix_phi", ix_phi)
        else
            ix_phi = istart
        endif

        if(cfg%has_key("tpse", "cfl_number")) then
            call cfg%get("tpse", "cfl_number", cfl_number)
        else
            cfl_number = MAXNUM 
        endif

        if(cfg%has_key("tpse", "ix_cfl")) then
            call cfg%get("tpse", "ix_cfl", ix_cfl)
        else
            ix_cfl = istart
        endif

    end subroutine cfg_npse

    subroutine cfg_lpse(cfg)

        use stdlib_string_type
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

        if (cfg%has_key('lpse', 'adjoint')) then
            call cfg%get('lpse', 'adjoint', lpse_IsAdjoint)
            if (lpse_IsAdjoint) then
                print *, 'ALPSE is solved.'
            else
                print *, 'LPSE is solved.'
            end if
        else
            print *, '[lpse]keyword adjoint is not given, use adjoint=.False.'
            call cfg%set('lpse', 'adjoint', .FALSE.)
            lpse_IsAdjoint = .False.
        end if

      !! this parameter control the dp/dx in LPSE, generally it can be set 1.
        if (cfg%has_key('lpse', 'sigma_p')) then
            call cfg%get('lpse', 'sigma_p', sigma_p)
        else
            sigma_p = 1.0d0
        end if

        if (sigma_p > 99.0d0) then
            sigma_p = 0.0d0
            sigma_control = 1.0d0
        else
            sigma_control = 0.0d0
        end if

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
                if(cfg%has_key("bc", "temperature_type"))then
                    call cfg%get("bc", "temperature_type", isothermal)
                  else
                    isothermal = .TRUE.
                  endif  
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

    ! subroutine bcast()

    !     implicit none
    !     integer :: ntmp(100)
    !     real(REAL64) :: rtmp(100)
    !     character(256) :: ctmp(100)

    !     integer :: ierr

    !     ntmp = 0

    !     ntmp(1) = in
    !     ntmp(2) = jn

    !     ntmp(4) = istart
    !     ntmp(5) = iend
    !     ntmp(6) = iloc
    !     ntmp(7) = 0
    !     ntmp(8) = 0
    !     ntmp(9) = 0
    !     ntmp(10) = bc_y0
    !     ntmp(11) = bc_yn
    !     ntmp(12) = Coordinate_System

    !     call MPI_Bcast(ntmp(1), 100, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    !     in = ntmp(1)
    !     jn = ntmp(2)

    !     istart = ntmp(4)
    !     iend = ntmp(5)
    !     iloc = ntmp(6)
    !     ! lpse_method=ntmp(7)
    !     ! lpse_mode=ntmp(8)
    !     ! lpse_nev=ntmp(9)
    !     bc_y0 = ntmp(10)
    !     bc_yn = ntmp(11)
    !     Coordinate_System = ntmp(12)

    !     rtmp = 0

    !     rtmp(1) = Ma
    !     rtmp(2) = Re
    !     rtmp(3) = Te
    !     rtmp(4) = alpha(1)
    !     rtmp(5) = alpha(2)
    !     rtmp(6) = omega(1)
    !     rtmp(7) = omega(2)
    !     rtmp(8) = beta(1)
    !     rtmp(9) = beta(2)
    !     rtmp(10) = sigma_p

    !     call MPI_bcast(rtmp(1), 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !     Ma = rtmp(1)
    !     Re = rtmp(2)
    !     Te = rtmp(3)
    !     alpha(1) = rtmp(4)
    !     alpha(2) = rtmp(5)
    !     omega(1) = rtmp(6)
    !     omega(2) = rtmp(7)
    !     beta(1) = rtmp(8)
    !     beta(2) = rtmp(9)
    !     sigma_p = rtmp(10)

    !     ctmp(1) = char(gridfile)
    !     ctmp(2) = char(flowfile)

    !     call MPI_bcast(ctmp(1), 100*256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    !     gridfile = trim(ctmp(1))
    !     flowfile = trim(ctmp(2))

    !     ! call MPI_Bcast(lpse_IsAdjoint, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    !     call MPI_Bcast(screen_flg, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    !     call MPI_bcast(mpi_ierr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    ! end subroutine bcast

end module mod_cfgio_adapter
