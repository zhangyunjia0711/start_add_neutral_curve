# LPSE
file(GLOB  SRC_FILES_LPSE
    "${PROJECT_SOURCE_DIR}/src/lpse/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lst/alst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LPSE ${SRC_FILES_LPSE})
target_link_libraries(START_LPSE PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK)  
install(TARGETS START_LPSE DESTINATION bin)
