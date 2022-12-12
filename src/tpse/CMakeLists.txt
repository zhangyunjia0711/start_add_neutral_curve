# tPSE
file(GLOB  SRC_FILES_NPSE
    "${PROJECT_SOURCE_DIR}/src/tpse/*.f90"
    "${PROJECT_SOURCE_DIR}/src/tpse/*.F90"
    "${PROJECT_SOURCE_DIR}/src/npse/module_petsc_fft.f90"
    "${PROJECT_SOURCE_DIR}/src/npse/module_fft.f90"
    "${PROJECT_SOURCE_DIR}/src/npse/npse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_tPSE ${SRC_FILES_NPSE})
target_link_libraries(START_tPSE PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK)  
install(TARGETS START_tPSE DESTINATION bin)
