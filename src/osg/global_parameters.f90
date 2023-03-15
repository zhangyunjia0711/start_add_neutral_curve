!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-29 15:17:05
!> @LastEditTime: 2023-03-15 11:35:15
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/osg/global_parameters.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ global parameters for optimal streak growth solver.
module global_parameters

  use iso_fortran_env, only: REAL64
  use petsc
  use cfgio_mod
  use stdlib_string_type, only: string_type
  implicit none

  public

  type(cfg_t) :: cfg
  logical :: mpi_ierr=.FALSE.                            !<cfg类
  logical :: screen_flg=.TRUE.
  ! integer :: lst_mode=LST_SPATIAL_MODE

  ! type(string) :: solname
  integer :: Coordinate_System
  logical :: flg_cycle= .FALSE.

  real(REAL64) :: GAMMA=1.4D0
  real(REAL64) :: Pr=0.72D0

  DM :: meshDA, coordDA, DA
  Vec :: Coord, Flowfield
  Vec :: Coord_local, Flowfield_local
  PetscErrorCode :: ierr

  type(string_type) :: gridfile, flowfile, curfile
  PetscInt :: in, jn, ln
  PetscInt :: istart, iend, iloc
  Real(REAL64) :: drxz(2)
  Real(REAL64) :: Ma, Re, Te

  real(REAL64) :: alpha(2), beta(2), omega(2)
  integer :: cycle_m
  logical :: lpse_IsAdjoint

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
  real(REAL64) :: rd_factor=1.3d0

  integer :: solve_alpha

  integer :: size=1, rank=0
  integer :: xbg, ybg, xlg, ylg, xeg, yeg
  integer :: xb, yb, xl, yl, xe, ye
  logical :: isothermal = .TRUE.

  integer :: i_x0, i_x1
  real(REAL64) :: initial_alf, initial_bat, initial_ome
  real(REAL64) :: initial_cycle_m
  real(REAL64) :: c_phase=1.0d0
  logical :: initial_flg

end module global_parameters
