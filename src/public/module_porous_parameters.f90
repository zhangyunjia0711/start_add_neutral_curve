!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-12-08 13:45:30
!> @LastEditTime: 2022-12-12 16:18:37
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/module_porous_parameters.f90
!> @Copyright (c) 2022 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
module mod_porous_parameters

    use iso_fortran_env
    use global_parameters, only: Ma, Re, Te, Pr, GAMMA
    use Complex_Bessel
    implicit none

    real(REAL64) :: C1 = 110.4D0
    real(REAL64) :: C2 = 110.4D0
    complex(REAL64), parameter :: ZI=(0.0d0, 1.0d0)


    real(REAL64) :: porosity, radius, depth
    
contains

    subroutine get_porous_AB(Rho0, T0, ome, porous_A, porous_B)

        implicit none
        real(REAL64), intent(in) :: Rho0, T0
        complex(REAL64), intent(in) :: ome
        complex(REAL64), intent(out) :: porous_A, porous_B
        complex(REAL64) :: Lambda
        complex(REAL64) :: Z_1, Y_1, Z_0
        complex(REAL64) :: kv, kt
        real(REAL64) :: cm, mu0

        cm = C1/Te
        mu0 = T0*sqrt(T0)*(1.0d0 + cm)/(T0 + cm)/Re

        kv = radius * sqrt(ZI*ome*Rho0/mu0)
        kt = kv*sqrt(Pr)
        Z_1=ZI*ome/T0 * bessel_jn(0, kv) / bessel_jn(2, kv)
        Y_1 = -ZI*ome*Ma**2*(GAMMA+(GAMMA-1.0d0)*bessel_jn(2, kt)/bessel_jn(0, kt))
        Z_0 = sqrt(Z_1/Y_1)

        Lambda = sqrt(Z_1*Y_1)

        porous_A = porosity/Z_0*tanh(Lambda*depth)
        porous_B = 0.0d0!-porosity*(GAMMA-1)*Ma**2*T0*bessel_jn(2, kt)/bessel_jn(0, kt)

    end

end module mod_porous_parameters