# LST
file(GLOB  SRC_FILES_LST
    "${PROJECT_SOURCE_DIR}/src/lst/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/*.F90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LST ${SRC_FILES_LST})
target_link_libraries(START_LST PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib) 
install(TARGETS START_LST DESTINATION bin)
