#
#  This is a sample CMakeLists.txt suitable for direct use with a new PETSc application or
#  you can add fragments of material below to an existing applications CMakeLists.txt
#
#  Usage:  To build ex1 from ex1.c
#    rm -fr build
#    mkdir build
#    cd build
#    cmake ..
#    make
#    ./ex1
#
#  By default it gets both the compiler information and the library information from PETSc
#
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_TYPE RELEASE)
if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL DEBUG))
   message("Debug mode: ${CMAKE_Fortran_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL RELEASE))
   message("Release mode: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig:$ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)

# Remove the lines below if you do not wish to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-cpp -fopenmp -fbacktrace -fbounds-check -std=legacy)
  add_link_options(-fopenmp -llapack)
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  add_compile_options(-cpp -DCOMPLEX_PETSC -traceback -CB -qmkl)
  add_link_options(-qmkl)
endif()

#  tells CMake to build the application ex1 from the source file ex1.c
#  this must appear AFTER the compilers are set
project(START)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_PREFIX $ENV{HOME})
find_package(PkgConfig REQUIRED)
# find_package(fortran_stdlib REQUIRED)
pkg_search_module(fortran_stdlib REQUIRED IMPORTED_TARGET fortran_stdlib)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET slepc)

add_subdirectory(src/lst)
add_subdirectory(src/lpse)
add_subdirectory(src/npse)
add_subdirectory(src/tpse)
add_subdirectory(src/osg)
add_subdirectory(src/lst_receptivity)
add_subdirectory(src/lpse_receptivity)
add_subdirectory(src/lst_receptivity_new)

