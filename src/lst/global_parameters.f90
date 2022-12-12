!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-09 10:45:42
!> @LastEditTime: 2022-12-12 19:33:08
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: global parameters for linear stability theory equations.
!> @FilePath: /START_all/src/lst/global_parameters.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------


#include <petsc/finclude/petsc.h>
#include <slepc/finclude/slepc.h>

!@ global parameters for lst solver.
module global_parameters

  use iso_fortran_env, only: REAL64
  use petsc
  use cfgio_mod
  use stdlib_string_type, only: string_type
  implicit none

  public

  type(cfg_t) :: cfg
  logical :: mpi_ierr=.FALSE.                            !<cfg类

  real(REAL64) :: GAMMA=1.4D0
  real(REAL64) :: Pr=0.72D0

  type(string_type) :: solname
  integer :: Coordinate_System

  DM :: meshDA, coordDA, DA
  Vec :: Coord, Flowfield
  Vec :: Coord_local, Flowfield_local
  PetscErrorCode :: ierr

  type(string_type) :: gridfile, flowfile, curfile
  PetscInt :: in, jn, ln
  PetscInt :: istart, iend, iloc
  Real(REAL64) :: drxz(2)
  Real(REAL64) :: Ma, Re, Te

  integer :: lst_method, lst_mode, lst_nev
  real(REAL64) :: alpha(2), beta(2), omega(2)
  logical :: lst_IsAdjoint

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

  integer :: size=1, rank=0
  integer :: xbg, ybg, xlg, ylg, xeg, yeg
  integer :: xb, yb, xl, yl, xe, ye
  logical :: isothermal = .TRUE.

end module global_parameters
