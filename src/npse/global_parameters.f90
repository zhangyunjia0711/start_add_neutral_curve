!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-29 15:17:05
!> @LastEditTime: 2023-04-04 09:52:11
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/npse/global_parameters.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ global parameters for npse solver.
module global_parameters

   use iso_fortran_env, only: REAL64
   use petsc
   use cfgio_mod
   use stdlib_string_type, only: string_type
   implicit none

   public

   type(cfg_t) :: cfg
   logical :: mpi_ierr = .FALSE.                            !<cfg类
   logical :: screen_flg = .TRUE.
   ! integer(kind(type_lst_mode)) :: lst_mode = LST_SPATIAL_MODE
   real(REAL64) :: GAMMA = 1.4D0
   real(REAL64) :: Pr = 0.72D0

   ! type(string_type) :: solname
   integer :: Coordinate_System

   type :: type_inlet_dis
      integer :: m
      integer :: n
      complex(REAL64) :: alf
      complex(REAL64) :: amp
      integer :: solve_alpha
   end type


   DM :: meshDA, coordDA, DA
   Vec :: Coord, Flowfield
   Vec :: Coord_local, Flowfield_local
   PetscErrorCode :: ierr

   type(string_type) :: gridfile, flowfile, curfile, title
   type(string_type) :: prefix
   PetscInt :: in, jn, ln
   PetscInt :: istart, iend, iloc
   Real(REAL64) :: drxz(2)
   Real(REAL64) :: Ma, Re, Te

   type(type_inlet_dis), allocatable :: inlet_dis(:)
   integer :: ndis
   real(REAL64) :: alpha(2), beta(2), omega(2)
   complex(REAL64) :: theta1, theta2
   integer :: m_global, n_global

   integer :: bc_y0, bc_yn

   real(REAL64), dimension(:, :), allocatable :: xx_0, yy_0
   real(REAL64), dimension(:, :, :), allocatable :: qq_0
   real(REAL64), dimension(:, :), allocatable :: xx, yy
   real(REAL64), dimension(:, :, :), allocatable :: qq
   real(REAL64), dimension(:), allocatable :: drx, drz

   real(REAL64), dimension(:, :), allocatable :: xi_x, xi_y, eta_x, eta_y
   real(REAL64), dimension(:, :), allocatable :: xi_xx, xi_yy, xi_xy
   real(REAL64), dimension(:, :), allocatable :: eta_xx, eta_yy, eta_xy

   real(REAL64), dimension(:, :, :), allocatable :: dq0dx, dq0dy
   real(REAL64), dimension(:, :, :), allocatable :: dq0dxx, dq0dyy
   real(REAL64), dimension(:, :, :), allocatable :: dq0dxy

   real(REAL64) :: sigma_p, sigma_control
   real(REAL64) :: rd_factor = 1.3d0

   ! integer(kind(type_solve_alpha_method)) :: solve_alpha

   logical :: lpse_IsAdjoint = .FALSE.

   integer :: solve_alpha
   integer :: size = 1, rank = 0
   integer :: xbg, ybg, xlg, ylg, xeg, yeg
   integer :: xb, yb, xl, yl, xe, ye

   logical :: isothermal = .TRUE.

   real(REAL64) :: omega_force_factor=1.0d0
   real(REAL64) :: omega_phi_factor=1.0d0
   real(REAL64) :: atol_nonlinear=1.0d-6
   real(REAL64) :: rtol_nonlinear=1.0d-3
   integer :: ix_force, ix_phi
   integer :: maxit_alpha=200
   integer :: maxit_nonlinear=5000

   real(REAL64) :: cfl_number=huge(1.0d0)
   integer :: ix_cfl

end module global_parameters
