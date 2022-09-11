# LPSE
file(GLOB  SRC_FILES_LPSE_RECEPTIVITY
    "${PROJECT_SOURCE_DIR}/src/lpse_receptivity/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse_receptivity/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/alpse_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lpse/lpse_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/alst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LPSE_Receptivity ${SRC_FILES_LPSE_RECEPTIVITY})
target_link_libraries(START_LPSE_Receptivity PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib)  
install(TARGETS START_LPSE_Receptivity DESTINATION bin)
