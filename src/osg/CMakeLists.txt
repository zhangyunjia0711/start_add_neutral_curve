# osg
file(GLOB  SRC_FILES_OSG
    "${PROJECT_SOURCE_DIR}/src/osg/*.f90"
    "${PROJECT_SOURCE_DIR}/src/osg/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/alpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_OSG ${SRC_FILES_OSG})
target_link_libraries(START_OSG PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK)  
install(TARGETS START_OSG DESTINATION bin)
