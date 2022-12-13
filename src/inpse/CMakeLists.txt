# inpse
file(GLOB  SRC_FILES_INPSE
    "${PROJECT_SOURCE_DIR}/src/inpse/*.f90"
    "${PROJECT_SOURCE_DIR}/src/inpse/*.F90"
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
add_executable(START_iNPSE ${SRC_FILES_INPSE})
target_link_libraries(START_iNPSE PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK)  
install(TARGETS START_iNPSE DESTINATION bin)
