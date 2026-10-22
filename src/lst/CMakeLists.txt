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
  add_executable(test_wall_local_profiles
    "${PROJECT_SOURCE_DIR}/src/public/compute_dx.f90"
    "${PROJECT_SOURCE_DIR}/src/public/wall_local_coordinates.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/global_parameters.f90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/string_conv_mod.f90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/cfgio_mod.f90"
    "${PROJECT_SOURCE_DIR}/tests/wall_local_profiles_test.f90")
  set_target_properties(test_wall_local_profiles PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/wall_local_test_modules")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/wall_local_test_modules")
  target_include_directories(test_wall_local_profiles PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/wall_local_test_modules")
  target_link_libraries(test_wall_local_profiles PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib)
  add_test(NAME wall_local_profiles COMMAND test_wall_local_profiles)
  add_test(NAME wall_local_reject_adjoint COMMAND test_wall_local_profiles reject_adjoint)
  add_test(NAME wall_local_reject_nonlocal COMMAND test_wall_local_profiles reject_nonlocal)
  set_tests_properties(wall_local_reject_adjoint wall_local_reject_nonlocal PROPERTIES WILL_FAIL TRUE)

  # These deterministic fixtures use small diagonal problems, not case scans.
  add_executable(test_slepc_candidate_retry
    "${PROJECT_SOURCE_DIR}/src/public/module_sort.f90"
    "${PROJECT_SOURCE_DIR}/src/public/module_slepc_adapter.f90"
    "${PROJECT_SOURCE_DIR}/tests/candidate_retry_test.f90")
  set_target_properties(test_slepc_candidate_retry PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/candidate_retry_test_modules")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/candidate_retry_test_modules")
  target_include_directories(test_slepc_candidate_retry PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/candidate_retry_test_modules")
  target_link_libraries(test_slepc_candidate_retry PkgConfig::PETSC PkgConfig::SLEPC)
  add_executable(test_raw_vector_resources
    "${PROJECT_SOURCE_DIR}/src/public/module_sort.f90"
    "${PROJECT_SOURCE_DIR}/src/public/module_slepc_adapter.f90"
    "${PROJECT_SOURCE_DIR}/src/public/compute_dx.f90"
    "${PROJECT_SOURCE_DIR}/src/public/wall_local_coordinates.f90"
    "${PROJECT_SOURCE_DIR}/src/public/lst_grid_profile.f90"
    "${PROJECT_SOURCE_DIR}/src/public/lst_baseflow_profile.f90"
    "${PROJECT_SOURCE_DIR}/src/public/solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/global_parameters.f90"
    "${PROJECT_SOURCE_DIR}/src/lst/lst_solver_abstract.f90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/string_conv_mod.f90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/cfgio_mod.f90"
    "${PROJECT_SOURCE_DIR}/tests/raw_vector_resources_test.f90")
  set_target_properties(test_raw_vector_resources PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/raw_vector_test_modules")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/raw_vector_test_modules")
  target_include_directories(test_raw_vector_resources PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/raw_vector_test_modules")
  target_link_libraries(test_raw_vector_resources PkgConfig::PETSC PkgConfig::SLEPC PkgConfig::fortran_stdlib)
  add_test(NAME slepc_candidate_retry COMMAND "${CMAKE_COMMAND}"
    "-DPROGRAM=$<TARGET_FILE:test_slepc_candidate_retry>" -DMODE=retry
    -P "${PROJECT_SOURCE_DIR}/tests/run_vector_regression.cmake")
  add_test(NAME slepc_candidate_retry_reject_eps COMMAND "${CMAKE_COMMAND}"
    "-DPROGRAM=$<TARGET_FILE:test_slepc_candidate_retry>" -DMODE=reject_eps
    -P "${PROJECT_SOURCE_DIR}/tests/run_vector_regression.cmake")
  add_test(NAME raw_vector_resources COMMAND "${CMAKE_COMMAND}"
    "-DPROGRAM=$<TARGET_FILE:test_raw_vector_resources>" -DMODE=raw
    -P "${PROJECT_SOURCE_DIR}/tests/run_vector_regression.cmake")
  set_tests_properties(slepc_candidate_retry slepc_candidate_retry_reject_eps raw_vector_resources
    PROPERTIES TIMEOUT 40 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endif()
