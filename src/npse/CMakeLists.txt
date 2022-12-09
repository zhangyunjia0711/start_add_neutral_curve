# NPSE
file(GLOB  SRC_FILES_NPSE
    "${PROJECT_SOURCE_DIR}/src/npse/*.f90"
    "${PROJECT_SOURCE_DIR}/src/npse/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_NPSE ${SRC_FILES_NPSE})
target_link_libraries(START_NPSE PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK)  
install(TARGETS START_NPSE DESTINATION bin)
