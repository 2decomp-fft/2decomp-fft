cmake_minimum_required(VERSION 3.12)                        
cmake_policy(SET CMP0074 NEW)

project(2decomp LANGUAGES Fortran)
enable_testing()

option(ENABLE_OPENACC "Allow user to activate/deactivate OpenACC support" OFF)
option(ENABLE_CUDA "Allow user to activate/deactivate CUDA support" OFF)

if (ENABLE_CUDA)
  enable_language(CUDA)
endif()

set(AUTHOR "Stefano Rolfo;Charles Moulinec;Paul Bartholomew")
set(AUTHOR_DETAILS "stefano.rolfo@stfc.ac.uk;charles.moulinec@stfc.ac.uk;p.bartholomew@epcc.ed.ac.uk")
set(DESCRIPTION "Building 2decomp&fft using cmake")

message(STATUS "building ${PROJECT_NAME}")

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/mod)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/opt" CACHE PATH "..." FORCE)
endif()

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# Find the modules included with Xcompact
#SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Dev Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

include(${CMAKE_SOURCE_DIR}/cmake/gpu.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/mpi.cmake)

# FFT options
set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw mkl)
include(${CMAKE_SOURCE_DIR}/cmake/fft/fft.cmake)

# Create a static library for the fft
add_subdirectory(src)

# Add tests
add_subdirectory(examples)

# # Create an example dir with all input.i3d example files
# option(BUILD_TESTING "Build with tests" ON)
# set(test_dir "${PROJECT_BINARY_DIR}/Test")
# message(STATUS "Before test main ${test_dir}")
# if (${BUILD_TESTING})
#   file(MAKE_DIRECTORY ${test_dir})
#   include(CTest)
#   message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
#   message(STATUS "MPI EXEC: ${MPIEXEC_EXECUTABLE}")
#   file(MAKE_DIRECTORY ${test_dir}/App)
#   add_test(NAME TestAPP COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d WORKING_DIRECTORY ${test_dir}/App)
#   #if (PFUNIT_FOUND)
#   #  file(MAKE_DIRECTORY ${test_dir}/Verif_x_00)
#   #  add_test(NAME Test_dx_00 COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/verif_x_00 WORKING_DIRECTORY ${test_dir}/Verif_x_00)
#   #endif (PFUNIT_FOUND)
# endif()

