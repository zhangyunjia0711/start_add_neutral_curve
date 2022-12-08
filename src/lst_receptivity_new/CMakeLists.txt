# LPSE
file(GLOB  SRC_FILES_LPSE_RECEPTIVITY_NEW
    "${PROJECT_SOURCE_DIR}/src/lst_receptivity_new/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lst_receptivity_new/*.F90"
    "${PROJECT_SOURCE_DIR}/src/lst/alst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_eqn.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LST_Receptivity_new ${SRC_FILES_LPSE_RECEPTIVITY_NEW})
target_link_libraries(START_LST_Receptivity_new PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib)  
install(TARGETS START_LST_Receptivity_new DESTINATION bin)
