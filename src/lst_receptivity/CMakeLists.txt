# LST
file(GLOB  SRC_FILES_LST_RECEPTIVITY
    "${PROJECT_SOURCE_DIR}/src/lst_receptivity/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lst_receptivity/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/alst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LST_Receptivity ${SRC_FILES_LST_RECEPTIVITY})
target_link_libraries(START_LST_Receptivity PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK) 
install(TARGETS START_LST_Receptivity DESTINATION bin)
