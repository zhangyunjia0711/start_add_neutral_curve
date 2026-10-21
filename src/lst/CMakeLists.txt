# LST
file(GLOB  SRC_FILES_LST
    "${PROJECT_SOURCE_DIR}/src/lst/*.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/*.F90"
    "${PROJECT_SOURCE_DIR}/src/public/*.f90"
    "${PROJECT_SOURCE_DIR}/src/public/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    )
add_executable(START_LST ${SRC_FILES_LST})
target_link_libraries(START_LST PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib LAPACK::LAPACK) 
install(TARGETS START_LST DESTINATION bin)

# Run via ctest --test-dir <build> -R neutral_numerics.
if(BUILD_TESTING)
  add_executable(test_neutral_numerics
    "${PROJECT_SOURCE_DIR}/src/lst/module_neutral_numerics.f90"
    "${PROJECT_SOURCE_DIR}/tests/test_neutral_numerics.f90")
  set_target_properties(test_neutral_numerics PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/neutral_test_modules")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/neutral_test_modules")
  add_test(NAME neutral_numerics COMMAND test_neutral_numerics)
  add_executable(test_slepc_adapter
    "${PROJECT_SOURCE_DIR}/src/public/module_sort.f90"
    "${PROJECT_SOURCE_DIR}/src/public/module_slepc_adapter.f90"
    "${PROJECT_SOURCE_DIR}/tests/adapter_safety_test.f90")
  set_target_properties(test_slepc_adapter PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/adapter_test_modules")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/adapter_test_modules")
  target_link_libraries(test_slepc_adapter PkgConfig::PETSC PkgConfig::SLEPC)
  add_test(NAME slepc_adapter_safety COMMAND test_slepc_adapter)
endif()
