real(REAL64), parameter :: MAXNUM=huge(1.0d0)
real(REAL64), parameter :: C1=110.4D0
real(REAL64), parameter :: C2=110.4D0

complex(REAL64), parameter :: ZI=(0.0d0, 1.0d0)

enum, bind(C)
   enumerator :: type_bc
   enumerator :: BC_WALL
   enumerator :: BC_FARFIELD
   enumerator :: BC_POLE
end enum

enum, bind(C)
   enumerator :: type_solve_alpha_method
   enumerator :: SOLVE_ALPHA_BOTH
   enumerator :: SOLVE_ALPHA_REAL
   enumerator :: SOLVE_ALPHA_NONE
end enum

enum, bind(C)
   enumerator :: type_coordinate_system
   enumerator :: COMOVING_CS
   enumerator :: CYLINDRICAL_CS
end enum

enum, bind(C)
   enumerator :: type_lst_mode
   enumerator :: LST_TEMPORAL_MODE
   enumerator :: LST_SPATIAL_MODE
end enum

enum, bind(C)
   enumerator :: type_lst_method
   enumerator :: LST_LOCAL
   enumerator :: LST_GLOBAL
   enumerator :: LST_DOMAIN
end enum

real(REAL64), parameter, dimension(-2:2, -2:2) :: FDM_C1_4ORD_CENTER=reshape( [&
0.0d0       ,        0.0d0, -3.0d0/2.0d0, 2.0d0      ,  -1.0d0/2.0d0, &
0.0d0       , -1.0d0/3.0d0, -1.0d0/2.0d0, 1.0d0      ,  -1.0d0/6.0d0, &
1.0d0/12.0d0, -2.0d0/3.0d0, 0.0d0       , 2.0d0/3.0d0, -1.0d0/12.0d0, &
1.0d0/6.0d0 ,       -1.0d0,  1.0d0/2.0d0, 1.0d0/3.0d0,         0.0d0, &
1.0d0/2.0d0 ,       -2.0d0,  3.0d0/2.0d0, 0.0d0      ,         0.0d0  &
], [5, 5])
real(REAL64), parameter, dimension(-2:2, -2:2) :: FDM_D1_4ORD_CENTER=reshape( [&
0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0, &
0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
-1.0d0/12.0d0, 4.0d0/3.0d0, -15.0d0/6.0d0, 4.0d0/3.0d0, -1.0d0/12.0d0, &
0.0d0        ,       1.0d0,        -2.0d0,       1.0d0,         0.0d0, &
0.0d0        ,       0.0d0,         0.0d0,       0.0d0,         0.0d0  &
], [5, 5])

